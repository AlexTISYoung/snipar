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
from snipar.imputation.preprocess_data import *
from snipar.imputation.impute_from_sibs import *
from snipar.pedigree import create_pedigree
from snipar.utilities import parse_obsfiles
import argparse
import h5py
import random
import pandas as pd
import os
from multiprocessing import Pool
from time import time
from pathlib import Path
from tempfile import TemporaryDirectory
from snipar.utilities import get_parser_doc
from snipar.numrange import parseNumRange, NumRangeAction
random.seed(1567924)

_MARKER_ROWS = {"bim_values"}
_MARKER_COLUMNS = {"imputed_par_gts", "num_observed_parental_alleles", "maf_coefs"}
_MARKER_VECTORS = {
    "pos", "non_duplicates", "ratio_ibd0", "mendelian_error_ratio",
    "estimated_genotyping_error", "parent_ratio_backup", "sib_ratio_backup",
    "standard_f", "maf_TSS", "maf_RSS1", "maf_RSS2", "maf_R2_1", "maf_R2_2",
    "maf_larger1", "maf_less0",
}


def _merge_chunks(chunks, filename, variant_count, compression, compression_opts):
    """Copy each retained marker span directly into the destination datasets."""
    bim_width = 0
    for name, _, _ in chunks:
        with h5py.File(name, "r") as data:
            bim_width = max(bim_width, data["bim_values"].dtype.itemsize)

    with h5py.File(filename, "w") as out:
        cursor = 0
        for name, selected, offset in chunks:
            with h5py.File(name, "r") as data:
                stop = cursor + selected.stop - selected.start
                for key, src in data.items():
                    if key in _MARKER_COLUMNS:
                        axis = 1
                    elif key in _MARKER_ROWS or key in _MARKER_VECTORS:
                        axis = 0
                    else:
                        if key not in out:
                            out[key] = src[()]
                        continue
                    if key not in out:
                        shape = list(src.shape)
                        shape[axis] = variant_count
                        if key == "maf_coefs":
                            out.create_dataset(key, shape, dtype=src.dtype,
                                               maxshape=(None, variant_count), fillvalue=np.nan)
                        elif key in {"imputed_par_gts", "num_observed_parental_alleles"}:
                            out.create_dataset(key, shape, dtype=src.dtype, chunks=True,
                                               compression=compression, compression_opts=compression_opts)
                        else:
                            dtype = f"S{bim_width}" if key == "bim_values" else src.dtype
                            out.create_dataset(key, shape, dtype=dtype)
                    dest = out[key]
                    if axis == 1:
                        if key == "maf_coefs" and src.shape[0] > dest.shape[0]:
                            dest.resize(src.shape[0], axis=0)
                        dest[:src.shape[0], cursor:stop] = src[:, selected]
                    elif key == "non_duplicates":
                        dest[cursor:stop] = src[selected] + offset
                    else:
                        dest[cursor:stop] = src[selected]
                cursor = stop


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
                    Address of the bed file (does not include '.bgen'). Only one of unphased_address and phased_address is necessary.
                
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
                    Number of the threads to be used. This should not exceed number of the available cores. The default number of the threads is one.

                chunks: int
                    Number of chunks load data in(each process).

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
    logging.info("processing " + str(phased_address) + "," + str(unphased_address))
    sibships, ibd, bim, chromosomes, ped_ids, pedigree_output = prepare_data(pedigree, phased_address, unphased_address, ibd_address, ibd_is_king, bim, fam, control, chromosome = chromosome, pedigree_nan=pedigree_nan)
    if chunks < 1:
        raise ValueError("chunks must be positive")
    if sibships.empty:
        raise ValueError("No families have the required sibling or parental genotypes for imputation")
    start_time = time()
    start = 0 if start is None else start
    end = len(bim) if end is None else end
    if not 0 <= start < end <= len(bim):
        raise ValueError("Imputation requires 0 <= start < end <= number of SNPs")

    def run_slice(first, last, out, compression=None, options=None, marker_indices=None):
        phased, dosages, ids, pos, freqs, metadata = prepare_gts(
            phased_address, unphased_address, bim, pedigree_output, ped_ids,
            chromosomes, first, last, pcs, pc_ids, find_optimal_pc,
            marker_indices=marker_indices,
        )
        impute(sibships, ids, phased, dosages, ibd, pos, metadata, str(chromosomes),
               freqs, out, threads=threads, output_compression=compression,
               output_compression_opts=options, silent_progress=silent_progress,
               use_backup=use_backup)

    if chunks == 1:
        run_slice(start, end, output_address, output_compression, output_compression_opts)
    else:
        _, markers = np.unique(bim["id"], return_index=True)
        markers.sort()
        markers = markers[(markers >= start) & (markers < end)]
        if not len(markers):
            raise ValueError("No unique SNPs in the requested slice")
        if find_optimal_pc and pcs is not None:
            sample = markers[np.linspace(0, len(markers) - 1, min(1000, len(markers)), dtype=int)]
            *_, model = prepare_gts(phased_address, unphased_address, bim, pedigree_output,
                                    ped_ids, chromosomes, start, end, pcs, pc_ids,
                                    find_optimal_pc=True, marker_indices=sample)
            pcs = pcs[:, :model["maf_x"].shape[1] - 1]
            find_optimal_pc = False
        # The phased IBD calculation uses a window of 50 unique markers on
        # either side. Include that context at each internal chunk boundary.
        half_window = 50 if phased_address else 0
        boundaries = np.linspace(0, len(markers), min(chunks, len(markers)) + 1, dtype=int)
        output = Path(f"{output_address}.hdf5")
        with TemporaryDirectory(prefix=".snipar-impute-", dir=output.parent) as directory:
            files = []
            for i, (first, last) in enumerate(zip(boundaries[:-1], boundaries[1:])):
                logging.info("imputing chunk %d/%d ...", i + 1, len(boundaries) - 1)
                first_read = max(0, first - half_window)
                last_read = min(len(markers), last + half_window)
                chunk_start = int(markers[first_read])
                chunk_end = int(markers[last_read - 1]) + 1
                address = str(Path(directory) / f"chunk{i}")
                run_slice(chunk_start, chunk_end, address, marker_indices=markers[first_read:last_read])
                selected = slice(first - first_read, last - first_read)
                files.append((f"{address}.hdf5", selected, chunk_start - start))
            merged = Path(directory) / "merged.hdf5"
            _merge_chunks(files, merged, len(markers), output_compression, output_compression_opts)
            os.replace(merged, output)
    return time() - start_time

#does the imputation and writes the results
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
    inputs = [{"pedigree": pedigree,
            "control":args.c,
            "use_backup":args.use_backup,
            "phased_address": none_tansform(args.bgen, "@", str(chromosome)),
            "unphased_address": none_tansform(args.bed, "@", str(chromosome)),
            "ibd_address": none_tansform(args.ibd, "@", str(chromosome)),            
            "ibd_is_king": args.ibd_is_king,
            "pcs": pcs,
            "pc_ids": pc_ids,
            "find_optimal_pc": args.find_optimal_pc,
            "output_address":output_for_chromosome(chromosome),
            "start": args.start,
            "end": args.end,
            "bim": none_tansform(args.bim, "@", str(chromosome)),
            "fam": none_tansform(args.fam, "@", str(chromosome)),
            "threads": args.threads,
            "chunks": args.chunks,
            "output_compression":args.output_compression,
            "output_compression_opts":args.output_compression_opts,
            "chromosome":chromosome,
            "pedigree_nan":args.pedigree_nan,
            'silent_progress':args.silent_progress
            }
            for chromosome in chromosomes]
    #TODO output more information about the imputation inside the hdf5 filehf
    if args.processes > 1:
        with Pool(args.processes) as pool:
            logging.info("staring process pool")        
            consumed_time = pool.map(run_imputation, inputs)
            logging.info("imputation time: "+str(np.sum(consumed_time)))
    else:
        start_time = time()
        for args in inputs:
            run_imputation(args)
        end_time = time()
        logging.info(f"imputation time: {end_time-start_time}")

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
                    type=str,help='Address of the phased genotypes in .bgen format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome(chr_range is an optional parameters for this script).')
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
                    default=1, 
                    help='Number of the threads to be used. This should not exceed number of the available cores. The default number of the threads is one.')
parser.add_argument('--processes',
                    type=int,
                    default=1,
                    help='Number of processes for imputation chromosomes. Each chromosome is done on one process.')
parser.add_argument('--chunks',
                    type=int,
                    default=1,
                    help='Number of chunks load data in(each process).')
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
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
    args=parser.parse_args()
    main(args)
