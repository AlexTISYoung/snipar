#!/usr/bin/env python
"""
This script performs imputation of missing parental genotypes from observed genotypes in a family. 
It can impute missing parents from families with no genotyped parents but at least two genotyped siblings, or one genotyped parent and one or more genotyped offspring. 
To specify the siblings, one can either provide a pedigree file (``--pedigree`` option) or
the relatedness inference output from KING with the ``--related`` ``--degree 1`` options along with age and sex information.

The pedigree file is a plain text file with header and columns: FID (family ID), IID (individual ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).
Note that individuals are assumed to have unique individual IDS (IID). Siblings are identified through individuals that have the same FID and the same FATHER_ID and MOTHER_ID.

Use the ``--king`` option to provide the KING relatedness inference output (usually has suffix .kin0) and the ``--agesex`` option to provide the age & sex information. The script
constructs a pedigree from this information and outputs it in the HDF5 output.

For phased BGEN or VCF input, a duoHMM-style model infers which offspring haplotype came from
the observed parent across the chromosome before imputing the other parent. Uncertain assignments
use the original imputation rules. Reuse calls from ``ibd.py`` with ``--parent-origin FILE`` or an
``.origin.hdf5`` sidecar next to the ``--ibd`` segments. Calls are checked against the input phased
genotypes and matched by sample, SNP ID, and position. The sidecar stores phased genotype evidence
for the called duos. Use ``--no_parent_origin`` to run the original phased method.

Use ``--ibd-posterior`` to supply marker-level IBD probabilities. A sibling posterior file
beside the ``--ibd`` segments is picked up automatically. For two-sibling families, imputation
averages over genotype-compatible IBD states and uses the segment estimate when probabilities
are unavailable. Larger families use hard states. The output preserves hard-state dosages
in ``robust_imputed_par_gts`` and reports the posterior-weighted allele count separately in
``num_expected_observed_parental_alleles``. ``posterior_num_sibships`` counts successful
two-sibling posterior imputations at each marker.

**Usage**

    ``impute.py --ibd chr_@.ibd --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4``

Args:
@parser@

Results:
    HDF5 files
        For each chromosome i, an HDF5 file is created at outprefix{i}. This file contains imputed genotypes, the position of SNPs, columns of resulting bim file, contents of resulting bim file, pedigree table and, family ids
        of the imputed parents, under the keys 'imputed_par_gts', 'pos', 'bim_columns', 'bim_values', 'pedigree' and, 'families', 'parental_status' respectively. There are also other details of the imputation in the resulting file. For more explanation see the documentation of snipar.imputation.impute_from_sibs.impute
        
"""
import logging
from snipar.imputation.preprocess_data import (
    ImputationGenotypes, frequency_model_aic, prepare_data, prepare_gts, prepare_imputation_pedigree,
)
from snipar.imputation.impute_from_sibs import ImputationFamilies, impute
from snipar.imputation.ibd_posterior import SiblingIBDPosteriors
from snipar.pedigree import create_pedigree
from snipar.utilities import parse_obsfiles
from snipar.genotype_cli import add_vcf_options, with_genotypes
import argparse
import h5py
import numpy as np
import pandas as pd
import os
import psutil
from contextlib import ExitStack, nullcontext
from multiprocessing import get_context
from time import time
from pathlib import Path
from tempfile import TemporaryDirectory, TemporaryFile
from types import SimpleNamespace
from snipar.utilities import get_parser_doc
from snipar.map import get_map_positions
from snipar.numrange import parseNumRange, NumRangeAction
from snipar.resources import (
    MIB, available_cpus, get_memory_budget, limit_native_threads, plan_batches, worker_environment,
)

_MARKER_ROWS = {"bim_values"}
_FAMILY_MARKER_COLUMNS = {"imputed_par_gts", "num_observed_parental_alleles",
                          "robust_imputed_par_gts", "num_expected_observed_parental_alleles"}
_MARKER_COLUMNS = _FAMILY_MARKER_COLUMNS | {"maf_coefs"}
_MARKER_VECTORS = {
    "pos", "non_duplicates", "ratio_ibd0", "mendelian_error_ratio",
    "estimated_genotyping_error", "parent_ratio_backup", "sib_ratio_backup",
    "posterior_num_sibships",
    "standard_f", "maf_TSS", "maf_RSS1", "maf_RSS2", "maf_R2_1", "maf_R2_2",
    "maf_larger1", "maf_less0",
}


def _write_block(file, values, target, total, bim_dtype, compression, compression_opts, order=None):
    """Write retained SNPs in bounded row tiles to the final HDF5 layout."""
    for key, source in values.items():
        if key in _MARKER_COLUMNS:
            axis = 1
        elif key in _MARKER_ROWS or key in _MARKER_VECTORS:
            axis = 0
        else:
            if key not in file:
                file[key] = source
            continue
        source = np.asarray(source)
        if key not in file:
            shape = list(source.shape)
            shape[axis] = total
            if key in _FAMILY_MARKER_COLUMNS:
                columns = min(128, source.shape[1])
                max_rows = max(1, MIB // (2 * columns))
                row_chunks = (shape[0] + max_rows - 1) // max_rows
                rows = (shape[0] + row_chunks - 1) // row_chunks
                file.create_dataset(key, shape, dtype="float16", chunks=(rows, columns),
                                    compression=compression, compression_opts=compression_opts)
            else:
                dtype = bim_dtype if key == "bim_values" else source.dtype
                file.create_dataset(key, shape, dtype=dtype)
        dest = file[key]
        if axis == 0:
            dest[target] = source if order is None else source[order]
        else:
            # Ask NumPy to round results to float16 before passing them to
            # HDF5, the same conversion used by the original direct writer.
            rows = max(1, 8 * MIB // max(1, source.shape[1] * dest.dtype.itemsize))
            for start in range(0, len(source), rows):
                end = min(start + rows, len(source))
                tile = source[start:end]
                if order is not None:
                    tile = tile[:, order]
                dest[start:end, target] = np.ascontiguousarray(tile, dtype=dest.dtype)


def _plan_imputation(reader, marker_count, family_count, max_sibs, threads, chunks, budget, n_duos=0,
                     posterior_pairs=None, *, local_duos=None):
    """Plan core and phase-context memory, with independent BGEN decoder tiles."""
    samples = reader.n_samples
    if reader.bgen is not None:
        # The probability decoder can need two full-sample arrays per
        # thread, independently of the retained sample subset.
        all_samples = len(reader.bgen.samples)
        decoder_bytes = max(1, 40 * all_samples)
        reader.read_threads = min(threads, max(1, budget // (decoder_bytes * 8)))
        reader.probability_bytes = min(128 * MIB, max(16 * samples, budget // 16))
        fixed_decode = decoder_bytes * reader.read_threads + reader.probability_bytes
    else:
        fixed_decode = 0
    if reader.with_phases:
        # Two haplotypes, a dosage, and at most 100 neighbouring SNPs
        # are retained. The core output and PC fit exclude the halo.
        phased_work_per_marker = 4 * samples
        pairs = max(1, max_sibs * (max_sibs - 1) // 2)
        phased_work_per_marker += 16 * threads * (pairs + max_sibs)
        halo_bytes = min(marker_count - 1, 100) * phased_work_per_marker
    else:
        phased_work_per_marker = 0
        halo_bytes = 0
    # Float64 imputed dosages and float32 observed-allele counts are
    # rounded into HDF5 with a small, independent staging buffer.
    marker_bytes = max(1, 12 * family_count + 4 * samples + 64 * threads + phased_work_per_marker)
    if n_duos:
        # Cached calls may coexist with phase validation and locally inferred duos.
        marker_bytes += 2 * n_duos
        if local_duos is None:
            local_duos = n_duos
        else:
            marker_bytes += 8 * n_duos
        if local_duos:
            marker_bytes += 32 * min(threads, local_duos)
            fixed_decode += 128 * local_duos + 32 * marker_count
    if posterior_pairs is not None:
        marker_bytes += 12 * family_count + 72 * posterior_pairs + 4 * threads
    if reader.matched_pcs is not None:
        marker_bytes += 72 * samples
    _, width = plan_batches(marker_count, marker_bytes, batch_size=(marker_count + chunks - 1) // chunks,
                            fixed_bytes=16 * MIB + fixed_decode + halo_bytes, memory_limit=budget)
    return width


def _parent_origin_pairs(sibships, genotype_rows, *, include_identifiers=False):
    pair_index = {}
    identifiers = [] if include_identifiers else None
    for father, mother, children, has_father, has_mother in sibships[
        ["FATHER_ID", "MOTHER_ID", "IID", "has_father", "has_mother"]
    ].itertuples(index=False, name=None):
        if has_father != has_mother:
            parent_id = father if has_father else mother
            parent = genotype_rows[parent_id]
            for child in children:
                pair = parent, genotype_rows[child]
                if pair not in pair_index:
                    pair_index[pair] = len(pair_index)
                    if identifiers is not None:
                        identifiers.append((parent_id, child))
    return np.array(list(pair_index), dtype=np.intp).reshape(-1, 2), pair_index, identifiers


def _imputation_alleles(bim, markers):
    if "allele_ids" in bim:
        return [value.split(",") for value in bim["allele_ids"].to_numpy()[markers]]
    return bim[["allele1", "allele2"]].to_numpy()[markers]


def _parent_origin_recombination(bim, markers, positions, mapfile):
    from snipar.imputation.duohmm import recombination_probabilities

    genetic_positions = None
    if mapfile is not None:
        genetic_positions = get_map_positions(
            mapfile, SimpleNamespace(pos=positions), min_map_prop=None, extrapolate=True,
        )
        genetic_positions -= min(0., genetic_positions[0])
    elif "morgans" in bim:
        # The third PLINK BIM field is in cM, as in SNIPAR's IBD input.
        cm = bim["morgans"].to_numpy(dtype=np.float64)[markers]
        if cm.size > 1 and np.any(cm != cm[0]):
            genetic_positions = cm
    return recombination_probabilities(positions, genetic_positions)


def _write_origin_checkpoints(genotypes, markers, pairs, recombination, boundaries, threads, checkpoint):
    from snipar.imputation.duohmm import advance_parental_origin

    alpha = None
    for first, stop in zip(boundaries[:-2], boundaries[1:-1]):
        phased, dosages = genotypes.read_calls(markers[first:stop])
        del dosages
        alpha = advance_parental_origin(phased, pairs, recombination[first:stop], alpha, threads=threads)
        checkpoint.write(alpha.tobytes())
        del phased


def _select_pc_model(genotypes, markers, width):
    """Choose one PC model for all blocks without holding the AIC sample in RAM."""
    sample = markers[np.linspace(0, len(markers) - 1, min(1000, len(markers)), dtype=int)]
    aic = np.zeros(genotypes.matched_pcs.shape[1] + 1)
    count = 0
    for start in range(0, len(sample), width):
        _, dosages = genotypes.read_calls(sample[start:start + width], include_phased=False)
        if not genotypes.pc_present.all():
            dosages = dosages[genotypes.pc_present]
        scores, terms = frequency_model_aic(dosages, genotypes.matched_pcs)
        aic += scores
        count += terms
    chosen = int(np.argmin(aic + count * 2 * np.arange(1, len(aic) + 1)))
    genotypes.matched_pcs = genotypes.matched_pcs[:, :chosen]


def _run_chromosome(pedigree, prepared_pedigree, pcs, pc_ids, args):
    with limit_native_threads(1):
        return run_imputation(dict(args, pedigree=pedigree, prepared_pedigree=prepared_pedigree,
                                   pcs=pcs, pc_ids=pc_ids))


def _init_chromosomes(pedigree, prepared_pedigree, pcs, pc_ids):
    global _chromosome_inputs
    _chromosome_inputs = pedigree, prepared_pedigree, pcs, pc_ids


def _run_chromosome_worker(args):
    return _run_chromosome(*_chromosome_inputs, args)


def run_imputation(data):
    """Runs the imputation and returns the consumed time
    Args:
        data : dict
            a dictionary with these keys and values:
            Keys:
                pedigree: pd.Dataframe
                    The standard pedigree table
                
                control : bool
                    Duplicates offsprings of families with more than one offspring and both parents and add '_' to the start of their FIDs.
                    These can be used for testing the imputation. The tests.test_imputation.imputation_test uses these.

                use_backup : bool, optional
                    Whether it should use backup imputation where there is no ibd information available.

                phased_address: str, optional
                    Address of the phased or unphased BGEN file (does not include '.bgen').
                
                unphased_address: str, optional
                    Address of the bed file (does not include '.bed'). Only one of unphased_address and phased_address is necessary.

                ibd_address : str
                    address of the ibd file. The king segments file should be accompanied with an allsegs file.

                ibd_is_king : boolean
                    Whether the ibd segments are in king format or snipar format.
                
                pcs : np.array[float], optional
                    A two-dimensional array containing pc scores for all individuals and SNPs respectively.

                pc_ids : set, optional
                    Set of ids of individuals in the pcs.

                find_optimal_pc : bool, optional
                    It will use Akaike information criterion to find the optimal number of PCs to use for MAF estimation.

                output_address: str
                    The address to write the result of imputation on. The default value for output_address is 'parent_imputed_chr'.

                start: int, optional
                    This function can do the imputation on a slice of each chromosome. If specified, his is the start of that slice(it is inclusive).

                end: int, optional
                    This function can do the imputation on a slice of each chromosome. If specified, this is the exclusive end of that slice.
                    
                bim: str, optional
                    Address of a bim file containing positions of SNPs if the address is different from Bim file of genotypes.

                fam: str, optional
                    Address of a fam file containing positions of SNPs if the address is different from Fam file of genotypes.

                threads: int, optional
                    Maximum OpenMP threads for this chromosome. The default uses the available CPU allowance.

                chunks: int
                    Minimum number of variant blocks. The memory budget can require more blocks.

                parent_origin: bool, optional
                    Infer parent of origin for phased data. Defaults to True.

                parent_origin_source: str, optional
                    Origin sidecar to reuse for matching phased duos and markers.

                parent_origin_map: str, optional
                    Genetic map for local parental-origin inference, with pposition and gposition columns in bp and cM.

                output_compression: str, optional
                    Optional compression algorithm used in writing the output as an hdf5 file. It can be either gzip or lzf. None means no compression.

                output_compression_opts: int, optional
                    Additional settings for the optional compression algorithm. Take a look at the create_dataset function of h5py library for more information. None means no compression setting.

                chromosome: str
                    name of the chromosome
                
                pedigree_nan: str
                    The value representing NaN in the pedigreee.

                silent_progress: bool
                    Hides the percentage of progress from logging
    Returns:
        float
            time consumed by the imputation.
    """
    chunks = data["chunks"]
    pedigree = data["pedigree"]
    control = data["control"]
    use_backup = data["use_backup"]
    phased_address = data.get("phased_address")
    unphased_address = data.get("unphased_address")
    pcs = data["pcs"]
    pc_ids = data["pc_ids"]
    find_optimal_pc = data["find_optimal_pc"]
    ibd_address = data.get("ibd_address")
    ibd_is_king = data.get("ibd_is_king")
    origin_path = data.get("parent_origin_source")
    posterior_path = data.get("ibd_posterior")
    if posterior_path is None and ibd_address is not None and not ibd_is_king:
        candidate = Path(f"{ibd_address}.posterior.hdf5")
        if candidate.is_file():
            posterior_path = str(candidate)
    output_address = data["output_address"]
    start = data.get("start")
    end = data.get("end")
    bim = data.get("bim")
    fam = data.get("fam")
    threads = data.get("threads")
    output_compression = data.get("output_compression")
    output_compression_opts = data.get("output_compression_opts")
    chromosome = data.get("chromosome")
    pedigree_nan = data.get("pedigree_nan")
    silent_progress = data.get("silent_progress")
    logging.info("Processing %s, %s", phased_address, unphased_address)
    if chunks < 1 or (threads is not None and threads < 1):
        raise ValueError("chunks and threads must be positive")
    threads = min(threads or available_cpus(), available_cpus())
    sibships, ibd, bim, chromosomes, ped_ids, pedigree_output = prepare_data(
        pedigree, phased_address, unphased_address, ibd_address, ibd_is_king, bim, fam,
        control, chromosome=chromosome, pedigree_nan=pedigree_nan,
        prepared_pedigree=data.get("prepared_pedigree"),
        has_ibd_posterior=posterior_path is not None,
    )
    if sibships.empty:
        raise ValueError("No families have the required sibling or parental genotypes for imputation")
    threads = min(threads, len(sibships))
    start_time = time()
    start = 0 if start is None else start
    end = len(bim) if end is None else end
    if not 0 <= start < end <= len(bim):
        raise ValueError("Imputation requires 0 <= start < end <= number of SNPs")
    markers = np.flatnonzero(~bim["id"].duplicated().to_numpy())
    markers = markers[(markers >= start) & (markers < end)]
    if not len(markers):
        raise ValueError("No unique SNPs in the requested slice")
    bim_dtype = bim.iloc[markers].to_numpy().astype("S").dtype
    output = Path(f"{output_address}.hdf5")
    posterior_source = nullcontext(None)
    if posterior_path is not None:
        posterior_path = str(posterior_path).replace("@", str(chromosomes[0]))
        posterior_source = SiblingIBDPosteriors(posterior_path, sibships, chromosome=chromosomes[0],
                                               include_larger_families=ibd_address is None)

    with ImputationGenotypes(phased_address, unphased_address, bim, ped_ids, pcs, pc_ids,
                             read_threads=threads) as genotypes, posterior_source as posteriors, ExitStack() as stack:
        pair_posterior_index = None if posteriors is None else posteriors.pair_index
        prepared = ImputationFamilies(sibships, genotypes.id_rows, ibd, genotypes.n_samples,
                                      posterior_index=pair_posterior_index)
        parent_origin = genotypes.with_phases and data.get("parent_origin", True) and sibships["single_parent"].any()
        if origin_path is not None and (not genotypes.with_phases or not data.get("parent_origin", True)):
            raise ValueError("--parent-origin requires phased genotypes and parental-origin inference")
        marker_order = None
        n_duos = n_local_duos = 0
        cached_origin = None
        if parent_origin:
            from snipar.imputation.duohmm import call_parental_origin

            source_path = origin_path
            if source_path is None and ibd_address is not None and not ibd_is_king:
                candidate = Path(f"{ibd_address}.origin.hdf5")
                if candidate.is_file():
                    source_path = str(candidate)
            pairs, pair_index, duo_ids = _parent_origin_pairs(
                sibships, genotypes.id_rows, include_identifiers=source_path is not None,
            )
            n_duos = len(pairs)
            positions = bim["coordinate"].to_numpy(dtype=np.int64)[markers]
            if np.any(positions[1:] < positions[:-1]):
                marker_order = np.argsort(positions, kind="stable")
                markers, positions = markers[marker_order], positions[marker_order]
            local_pairs = pairs
            if source_path is not None:
                from snipar.imputation.origin_io import OriginReader

                source_path = str(source_path).replace("@", str(chromosomes[0]))
                saved_origin = stack.enter_context(OriginReader(source_path, chromosomes[0], duo_ids))
                del duo_ids
                matched_markers = saved_origin.align(positions, bim["id"].to_numpy()[markers],
                                                     _imputation_alleles(bim, markers))
                if len(saved_origin.target_rows) and matched_markers:
                    cached_origin = saved_origin
                    local_pairs = pairs[saved_origin.missing_rows]
                    logging.info("Reusing parental origin for %d duos at %d of %d markers from %s",
                                 len(saved_origin.target_rows), matched_markers, len(markers), source_path)
                elif origin_path is not None:
                    raise ValueError(f"No matching parent-offspring pairs and markers in {source_path}")
                else:
                    logging.info("No matching parental origin calls in %s", source_path)
            n_local_duos = len(local_pairs)
            if n_local_duos:
                recombination = _parent_origin_recombination(bim, markers, positions, data.get("parent_origin_map"))
                logging.info("Inferring parent of origin for %d parent-offspring pairs", n_local_duos)
        if posteriors is not None:
            positions = bim["coordinate"].to_numpy(dtype=np.int64)[markers]
            identifiers = bim["id"].to_numpy(dtype="S")[markers]
            posteriors.align(positions, identifiers)
            logging.info("Matched %d sibling pairs to IBD posteriors", len(pair_posterior_index))
        genotypes.id_rows = None
        del ibd, ped_ids
        budget = get_memory_budget()
        if data.get("worker_budget") is not None:
            budget = min(budget, data["worker_budget"] - psutil.Process().memory_info().rss)
        budget = max(0, budget)
        width = _plan_imputation(genotypes, len(markers), len(sibships),
                                 int(sibships["sib_count"].max()), threads, chunks, budget, n_duos,
                                 None if posteriors is None else len(pair_posterior_index),
                                 local_duos=None if cached_origin is None else n_local_duos)
        count = max(min(chunks, len(markers)), (len(markers) + width - 1) // width)
        boundaries = np.linspace(0, len(markers), count + 1, dtype=int)
        logging.info("Using %d variant blocks (%d SNPs maximum, %.0f MiB budget, %d threads)",
                     count, width, budget / MIB, threads)
        if find_optimal_pc and genotypes.matched_pcs is not None:
            _select_pc_model(genotypes, markers, width)

        # The phased rolling IBD comparison uses 50 unique markers on
        # either side. Every internal block includes this context.
        half_window = 50 if genotypes.with_phases else 0
        with TemporaryDirectory(prefix=".snipar-impute-", dir=output.parent) as directory:
            temporary_output = Path(directory) / "output.hdf5"
            # Carry the HMM across block boundaries. Store one alpha per duo and block
            # on disk, then smooth in reverse order while writing the final output.
            scratch = TemporaryFile(dir=directory) if n_local_duos and count > 1 else nullcontext(None)
            with scratch as checkpoint, h5py.File(temporary_output, "w") as file:
                blocks = range(count)
                beta = None
                if checkpoint is not None:
                    logging.info("Computing parent-of-origin forward checkpoints for %d blocks", count)
                    _write_origin_checkpoints(genotypes, markers, local_pairs, recombination, boundaries,
                                              threads, checkpoint)
                    checkpoint_bytes = n_local_duos * 4 * np.dtype(np.float64).itemsize
                    blocks = reversed(blocks)
                for completed, block in enumerate(blocks, 1):
                    first, stop = boundaries[block:block + 2]
                    first_read = max(0, first - half_window)
                    stop_read = min(len(markers), stop + half_window)
                    retained = slice(first - first_read, stop - first_read)
                    if count > 1:
                        logging.info("Imputing variant block %d (%d/%d)", block + 1, completed, count)
                    phased, dosages, _, pos, freqs, metadata = prepare_gts(
                        phased_address, unphased_address, bim, pedigree_output, None,
                        chromosomes, start, end, marker_indices=markers[first_read:stop_read],
                        prepared=genotypes, retained=retained,
                    )
                    origin_args = {}
                    if parent_origin:
                        haplotypes = phased[:, retained, :]
                        if n_local_duos:
                            alpha = None
                            if block:
                                checkpoint.seek((block - 1) * checkpoint_bytes)
                                alpha = np.frombuffer(checkpoint.read(checkpoint_bytes), dtype=np.float64)
                                alpha = alpha.reshape(n_local_duos, 4)
                            origin, parent_haplotype, beta = call_parental_origin(
                                haplotypes, local_pairs, recombination[first:stop], alpha, beta, threads=threads,
                            )
                        if cached_origin is not None:
                            saved = np.full((n_duos, stop - first), -1, dtype=np.int8)
                            copies = np.full_like(saved, -1)
                            if n_local_duos:
                                saved[cached_origin.missing_rows] = origin
                                copies[cached_origin.missing_rows] = parent_haplotype
                            cached_origin.read_into(first, stop, haplotypes, pairs, saved, copies)
                            origin, parent_haplotype = saved, copies
                            del saved, copies
                        del haplotypes
                        origin_args = dict(parent_origin=origin, parent_haplotype=parent_haplotype,
                                           parent_origin_index=pair_index)
                    posterior_block = None if posteriors is None else posteriors.read(first, stop)
                    impute(sibships, None, phased, dosages, None, pos, metadata, str(chromosomes),
                           freqs, threads=threads, silent_progress=silent_progress or count > 1,
                           use_backup=use_backup, prepared=prepared, return_genotypes=False,
                           output_start=retained.start, output_stop=retained.stop, ibd_posteriors=posterior_block,
                           posterior_hard_fallback=ibd_address is None, **origin_args)
                    if parent_origin:
                        del origin, parent_haplotype, origin_args
                    del phased, dosages, freqs, posterior_block
                    target, order = slice(first, stop), None
                    if marker_order is not None:
                        targets = marker_order[first:stop]
                        order = np.argsort(targets)
                        target = targets[order]
                    _write_block(file, metadata, target, len(markers), bim_dtype,
                                 output_compression, output_compression_opts, order=order)
                    # Do not carry the completed genotype planes into the next read.
                    del metadata
            os.replace(temporary_output, output)
    return time() - start_time

def _load_pedigree(args):
    if args.pedigree:
        logging.info("Reading pedigree")
        pedigree = pd.read_csv(args.pedigree, sep=r"\s+", dtype=str, keep_default_na=False)
    else:
        logging.info("Creating pedigree")
        pedigree = create_pedigree(args.king, args.agesex)
    return pedigree[["FID", "IID", "FATHER_ID", "MOTHER_ID"]]


def _vcf_samples(args, available_ids):
    pedigree = _load_pedigree(args)
    pedigree_ids = set(pedigree["IID"])
    pedigree_ids.update(pedigree["FATHER_ID"], pedigree["MOTHER_ID"])
    pedigree_ids.intersection_update(available_ids)
    prepared = prepare_imputation_pedigree(pedigree, pedigree_ids, control=args.c,
                                           pedigree_nan=args.pedigree_nan)
    if prepared.sibships.empty:
        raise ValueError("No families have the required sibling or parental genotypes for imputation")
    # The output pedigree describes the original observations, including
    # complete trios whose genotypes are not needed for parental imputation.
    args._vcf_pedigree = prepared
    return {iid.decode() for iid in prepared.sample_ids}


@with_genotypes(output='auto', vcf_samples=_vcf_samples)
def main(args):
    """"Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    if args.bgen is None and args.bed is None:
        raise Exception("You should supplement the code with at least one genotype address") 
    
    from importlib.metadata import PackageNotFoundError, version
    try:
        logging.info("snipar version: %s", version("snipar"))
    except PackageNotFoundError:
        logging.info("Running snipar from an uninstalled source tree")

    prepared_pedigree = getattr(args, "_vcf_pedigree", None)
    pedigree = _load_pedigree(args) if prepared_pedigree is None else None
    pedigree_rows = len(pedigree) if prepared_pedigree is None else len(prepared_pedigree.table) - 1
    logging.info("Pedigree loaded.")

    pcs = None
    pc_ids = None
    if args.pcs:
        logging.info("loading pcs ...")
        #ordering should be the same
        pcs = pd.read_csv(args.pcs, sep=r"\s+", dtype=str, keep_default_na=False)
        pc_ids = pcs.iloc[:, 1].to_numpy(dtype="S")
        pcs = pcs.iloc[:, 2:].to_numpy(dtype=float)
        if not args.pc_num is None:
            if args.pc_num > pcs.shape[1]:
                raise Exception(f"pc_num={args.pc_num} more than {pcs.shape} PCs available in pcs")
            else:
                pcs = pcs[:,:args.pc_num]
        logging.info("pcs loaded")
    
    chromosomes = [None]
    if args.chr_range:
        chromosomes = args.chr_range
        if args.bed:
            if not ("@" in args.bed):
                raise Exception("with chr_range bed address requires a @ wildcard")
        elif args.bgen:
            if not ("@" in args.bgen):
                raise Exception("with chr_range bgen address requires a @ wildcard")
    else:
        if (args.bed and "@" in args.bed):
            files, chromosomes = parse_obsfiles(args.bed, "bed", False)
            chromosomes = chromosomes.astype('str').tolist()
        elif(args.bgen and "@" in args.bgen):
            files, chromosomes = parse_obsfiles(args.bgen, "bgen", False)
            chromosomes = chromosomes.astype('str').tolist()

    if args.bed and args.bgen:
        if ("@" in args.bed) ^ ("@" in args.bgen):
            raise Exception("Can not have @ for just one of the bed and bgen")

    if args.ibd is None:
        if args.ibd_is_king:
            raise Exception("can not use 'ibd_is_king' without any ibd file")

    def none_tansform(a, b, c):
        if a is not None:
            return a.replace(b, c)
        return None
    def output_for_chromosome(chrom):
        if len(chromosomes) > 1 and "@" not in args.out:
            return f"{args.out}{chrom}"
        return none_tansform(args.out, "@", str(chrom))
    if args.processes < 1 or (args.threads is not None and args.threads < 1):
        raise ValueError("processes and threads must be positive")
    cpu_count = available_cpus()
    budget = get_memory_budget()
    # Subprocesses share neither Python metadata nor native scratch under
    # spawn. Give each a bounded share, including the cost of its metadata.
    fixed_per_worker = 256 * MIB + pedigree_rows * 1024
    if pcs is not None:
        fixed_per_worker += 2 * pcs.nbytes
    if args.bgen:
        fixed_per_worker += 550 * pedigree_rows
    processes = min(args.processes, len(chromosomes), cpu_count, max(1, budget // fixed_per_worker))
    threads = min(args.threads or cpu_count, max(1, cpu_count // processes))
    inputs = [{"control":args.c,
            "use_backup":args.use_backup,
            "phased_address": none_tansform(args.bgen, "@", str(chromosome)),
            "unphased_address": none_tansform(args.bed, "@", str(chromosome)),
            "ibd_address": none_tansform(args.ibd, "@", str(chromosome)),            
            "ibd_is_king": args.ibd_is_king,
            "ibd_posterior": getattr(args, "ibd_posterior", None),
            "find_optimal_pc": args.find_optimal_pc,
            "output_address":output_for_chromosome(chromosome),
            "start": args.start,
            "end": args.end,
            "bim": none_tansform(args.bim, "@", str(chromosome)),
            "fam": none_tansform(args.fam, "@", str(chromosome)),
            "threads": threads,
            "worker_budget": budget // processes if processes > 1 else None,
            "chunks": args.chunks,
            "parent_origin": not getattr(args, "no_parent_origin", False),
            "parent_origin_source": getattr(args, "parent_origin", None),
            "parent_origin_map": none_tansform(getattr(args, "parent_origin_map", None), "@", str(chromosome)),
            "output_compression":args.output_compression,
            "output_compression_opts":args.output_compression_opts,
            "chromosome":chromosome,
            "pedigree_nan":args.pedigree_nan,
            'silent_progress':args.silent_progress
            }
            for chromosome in chromosomes]
    start_time = time()
    logging.info("Allocating %d chromosome processes and %d threads per process", processes, threads)
    chromosome_inputs = pedigree, prepared_pedigree, pcs, pc_ids
    if processes > 1:
        with worker_environment(1), get_context("spawn").Pool(
            processes, initializer=_init_chromosomes, initargs=chromosome_inputs,
        ) as pool:
            pool.map(_run_chromosome_worker, inputs, chunksize=1)
    else:
        for job in inputs:
            _run_chromosome(*chromosome_inputs, job)
    logging.info("Imputation time: %.2f seconds", time() - start_time)

parser = argparse.ArgumentParser()
parser.add_argument('-c',
                    action='store_true',
                    help = "Duplicates offsprings of families with more than one offspring and both parents and add '_' to the start of their FIDs. These can be used for testing the imputation. The tests.test_imputation.imputation_test uses these.")
parser.add_argument('-silent_progress',
                    action='store_true',
                    help = "Hides the percentage of progress from logging")
parser.add_argument('-use_backup',
                    action='store_true',
                    help = "Whether it should use backup imputation where there is no ibd information available")
parser.add_argument('--ibd',
                    type=str,
                    help='Address of the IBD file without suffix. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome(chr_range is an optional parameters for this script).')
parser.add_argument('--ibd_is_king',
                    action='store_true',
                    help='If not provided the ibd input is assumed to be in snipar. Otherwise its in king format with an allsegs file')
parser.add_argument('--ibd-posterior', help='Sibling IBD posterior HDF5 file. Use @ for chromosome-specific files.')
parser.add_argument('--bgen',
                    type=str, help='Prefix of phased or unphased .bgen genotypes. Use @ to substitute chromosome numbers.')
parser.add_argument('--bed',
                    type=str,help='Address of the unphased genotypes in .bed format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome(chr_range is an optional parameters for this script).')
parser.add_argument('--chr_range',
                    type=parseNumRange,
                    nargs='*',
                    action=NumRangeAction,
                    help='number of the chromosomes to be imputed. Should be a series of ranges with x-y format or integers.')
parser.add_argument('--bim',
                    type=str,
                    default = None,
                    help='Address of a bim file containing positions of SNPs if the address is different from Bim file of genotypes')
parser.add_argument('--fam',
                    type=str,
                    default = None,
                    help='Address of a fam file containing positions of SNPs if the address is different from fam file of genotypes')
parser.add_argument('--out',
                    type=str,
                    default = "parent_imputed",
                    help="Writes the result of imputation for chromosome i to outprefix{i}")
parser.add_argument('--start',
                    type=int,
                    help='The script can do the imputation on a slice of each chromosome. This is the start of that slice(it is inclusive)',
                    default=None)
parser.add_argument('--end',
                    type=int,
                    help='The script can do the imputation on a slice of each chromosome. Tthis is the exclusive end of that slice.',
                    default=None)
parser.add_argument('--pedigree',
                    type=str,
                    default = None,
                    help="Address of the pedigree file. Pedigree file is ' ' separated with columns 'FID', 'IID', 'FATHER_ID', 'MOTHER_ID'. Default NaN value of Pedigree file is '0'. If your NaN value is something else be sure to specify it with --pedigree_nan option.")
parser.add_argument('--king',
                    type=str,
                    default = None,
                    help="""Address of a kinship file in KING format. The tab-delimited file has columns "FID1", "ID1", "FID2", "ID2", and "InfType".
            Each row represents a relationship between two individuals. InfType column states the relationship between two individuals.
            The only relationships that matter for this script are full sibling and parent-offspring which are shown by 'FS' and 'PO' respectively.
            This file is used in creating a pedigree file and can be generated using KING.""")
parser.add_argument('--agesex',
                    type=str,
                    default = None,
                    help="""Address of the agesex file. This is a " " separated CSV with columns "FID", "IID", "FATHER_ID", "MOTHER_ID", "sex", "age".
            Each row contains the age and sex of one individual. Male and Female sex should be represented with 'M' and 'F'.
            Age column is used for distinguishing between parent and child in a parent-offsring relationship inferred from the kinship file.
            ID1 is a parent of ID2 if there is a 'PO' relationship between them and 'ID1' is at least 12 years older than ID2.""")
parser.add_argument('--pcs',
                    type=str,
                    default = None,
                    help='Address of the PCs file with header "FID IID PC1 PC2 ...". If provided MAFs will be estimated from PCs')
parser.add_argument('--pc_num',
                    type=int,
                    default = None,
                    help='Number of PCs to consider')
parser.add_argument('-find_optimal_pc',
                    action='store_true',
                    help='It will use Akaike information criterion to find the optimal number of PCs to use for MAF estimation.')
parser.add_argument('--threads',
                    type=int,
                    default=None,
                    help='Maximum OpenMP threads per chromosome. By default use the available CPU allowance.')
parser.add_argument('--processes',
                    type=int,
                    default=1,
                    help='Maximum concurrent chromosomes, subject to available CPU and memory.')
parser.add_argument('--chunks',
                    type=int,
                    default=1,
                    help='Minimum number of variant blocks per chromosome. Available memory may require more blocks.')
origin_options = parser.add_mutually_exclusive_group()
origin_options.add_argument('--no_parent_origin', action='store_true',
                            help='Use the original phased imputation method without parental-origin inference.')
origin_options.add_argument('--parent-origin', metavar='FILE',
                            help='Reuse duoHMM labels from an IBD origin HDF5 sidecar. Use @ for chromosome numbers.')
parser.add_argument('--parent_origin_map', type=str,
                    help='Genetic map for uncached parental-origin inference with pposition (bp) and gposition (cM) columns. '
                         'Use @ for chromosome numbers. The default uses BIM cM if available or 1 cM per Mb.')
parser.add_argument('--output_compression',
                    type=str,
                    default=None,
                    help='Optional compression algorithm used in writing the output as an hdf5 file. It can be either gzip or lzf')
parser.add_argument('--output_compression_opts',
                    type=int,
                    default=None,
                    help='Additional settings for the optional compression algorithm. Take a look at the create_dataset function of h5py library for more information.')
parser.add_argument('--pedigree_nan',
                    type=str,
                    default='0',
                    help='The value representing NaN in the pedigreee.')
add_vcf_options(parser, bgen_quality=True)
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
    args=parser.parse_args()
    main(args)
