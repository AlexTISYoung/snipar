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
    ImputationGenotypes, frequency_model_aic, prepare_data, prepare_gts,
)
from snipar.imputation.impute_from_sibs import ImputationFamilies, impute
from snipar.pedigree import create_pedigree
from snipar.utilities import parse_obsfiles
from snipar.genotype_cli import add_vcf_options, with_genotypes
import argparse
import h5py
import numpy as np
import pandas as pd
import os
import psutil
from multiprocessing import get_context
from time import time
from pathlib import Path
from tempfile import TemporaryDirectory
from snipar.utilities import get_parser_doc
from snipar.numrange import parseNumRange, NumRangeAction
from snipar.resources import (
    MIB, available_cpus, get_memory_budget, limit_native_threads, plan_batches, worker_environment,
)

_MARKER_ROWS = {"bim_values"}
_MARKER_COLUMNS = {"imputed_par_gts", "num_observed_parental_alleles", "maf_coefs"}
_MARKER_VECTORS = {
    "pos", "non_duplicates", "ratio_ibd0", "mendelian_error_ratio",
    "estimated_genotyping_error", "parent_ratio_backup", "sib_ratio_backup",
    "standard_f", "maf_TSS", "maf_RSS1", "maf_RSS2", "maf_R2_1", "maf_R2_2",
    "maf_larger1", "maf_less0",
}


def _write_block(file, values, first, stop, total, bim_dtype, compression, compression_opts):
    """Append retained SNPs in bounded row tiles to the final HDF5 layout."""
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
            if key in {"imputed_par_gts", "num_observed_parental_alleles"}:
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
            dest[first:stop] = source
        else:
            # Ask NumPy to round results to float16 before passing them to
            # HDF5, the same conversion used by the original direct writer.
            rows = max(1, 8 * MIB // max(1, source.shape[1] * dest.dtype.itemsize))
            for start in range(0, len(source), rows):
                end = min(start + rows, len(source))
                dest[start:end, first:stop] = np.ascontiguousarray(source[start:end], dtype=dest.dtype)


def _plan_imputation(reader, marker_count, family_count, max_sibs, threads, chunks, budget):
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
    if reader.matched_pcs is not None:
        marker_bytes += 72 * samples
    _, width = plan_batches(marker_count, marker_bytes, batch_size=(marker_count + chunks - 1) // chunks,
                            fixed_bytes=16 * MIB + fixed_decode + halo_bytes, memory_limit=budget)
    return width


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


def _run_chromosome(pedigree, pcs, pc_ids, args):
    with limit_native_threads(1):
        return run_imputation(dict(args, pedigree=pedigree, pcs=pcs, pc_ids=pc_ids))


def _init_chromosomes(pedigree, pcs, pc_ids):
    global _chromosome_inputs
    _chromosome_inputs = pedigree, pcs, pc_ids


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

    with ImputationGenotypes(phased_address, unphased_address, bim, ped_ids, pcs, pc_ids,
                             read_threads=threads) as genotypes:
        prepared = ImputationFamilies(sibships, genotypes.id_rows, ibd, genotypes.n_samples)
        genotypes.id_rows = None
        del ibd, ped_ids
        budget = get_memory_budget()
        if data.get("worker_budget") is not None:
            budget = min(budget, data["worker_budget"] - psutil.Process().memory_info().rss)
        budget = max(0, budget)
        width = _plan_imputation(genotypes, len(markers), len(sibships),
                                 int(sibships["sib_count"].max()), threads, chunks, budget)
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
            with h5py.File(temporary_output, "w") as file:
                for i, (first, stop) in enumerate(zip(boundaries[:-1], boundaries[1:]), 1):
                    first_read = max(0, first - half_window)
                    stop_read = min(len(markers), stop + half_window)
                    retained = slice(first - first_read, stop - first_read)
                    if count > 1:
                        logging.info("Imputing variant block %d/%d", i, count)
                    phased, dosages, _, pos, freqs, metadata = prepare_gts(
                        phased_address, unphased_address, bim, pedigree_output, None,
                        chromosomes, start, end, marker_indices=markers[first_read:stop_read],
                        prepared=genotypes, retained=retained,
                    )
                    impute(sibships, None, phased, dosages, None, pos, metadata, str(chromosomes),
                           freqs, threads=threads, silent_progress=silent_progress or count > 1,
                           use_backup=use_backup, prepared=prepared, return_genotypes=False,
                           output_start=retained.start, output_stop=retained.stop)
                    del phased, dosages, freqs
                    _write_block(file, metadata, first, stop, len(markers), bim_dtype,
                                 output_compression, output_compression_opts)
                    # Do not carry the completed genotype planes into the next read.
                    del metadata
            os.replace(temporary_output, output)
    return time() - start_time

#does the imputation and writes the results
@with_genotypes(output='auto')
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

    #fids starting with _ are reserved for control
    #Families should not have grandparents
    if not args.pedigree:
        logging.info("creating pedigree ...")
        pedigree = create_pedigree(args.king, args.agesex)
    else:
        logging.info("reading pedigree ...")
        pedigree = pd.read_csv(args.pedigree, sep=r"\s+", dtype=str, keep_default_na=False)
    pedigree = pedigree[['FID', 'IID', 'FATHER_ID', 'MOTHER_ID']]
    logging.info("pedigree loaded.")
    
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
    fixed_per_worker = 256 * MIB + len(pedigree) * 1024
    if pcs is not None:
        fixed_per_worker += 2 * pcs.nbytes
    if args.bgen:
        fixed_per_worker += 550 * len(pedigree)
    processes = min(args.processes, len(chromosomes), cpu_count, max(1, budget // fixed_per_worker))
    threads = min(args.threads or cpu_count, max(1, cpu_count // processes))
    inputs = [{"control":args.c,
            "use_backup":args.use_backup,
            "phased_address": none_tansform(args.bgen, "@", str(chromosome)),
            "unphased_address": none_tansform(args.bed, "@", str(chromosome)),
            "ibd_address": none_tansform(args.ibd, "@", str(chromosome)),            
            "ibd_is_king": args.ibd_is_king,
            "find_optimal_pc": args.find_optimal_pc,
            "output_address":output_for_chromosome(chromosome),
            "start": args.start,
            "end": args.end,
            "bim": none_tansform(args.bim, "@", str(chromosome)),
            "fam": none_tansform(args.fam, "@", str(chromosome)),
            "threads": threads,
            "worker_budget": budget // processes if processes > 1 else None,
            "chunks": args.chunks,
            "output_compression":args.output_compression,
            "output_compression_opts":args.output_compression_opts,
            "chromosome":chromosome,
            "pedigree_nan":args.pedigree_nan,
            'silent_progress':args.silent_progress
            }
            for chromosome in chromosomes]
    start_time = time()
    logging.info("Allocating %d chromosome processes and %d threads per process", processes, threads)
    if processes > 1:
        with worker_environment(1), get_context("spawn").Pool(
            processes, initializer=_init_chromosomes, initargs=(pedigree, pcs, pc_ids),
        ) as pool:
            pool.map(_run_chromosome_worker, inputs, chunksize=1)
    else:
        for job in inputs:
            _run_chromosome(pedigree, pcs, pc_ids, job)
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
