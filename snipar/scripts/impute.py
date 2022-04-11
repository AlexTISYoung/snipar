#!/usr/bin/env python
"""Runs the sib-regression for the specified source and writes the result in an HDF5 file.

This script performs imputation of missing parental genotypes from observed genotypes in a family. It can impute missing parents from families with no genotyped parents but at least two genotyped siblings, or one genotyped parent and one or more genotyped offspring. To specify the siblings, one can either provide a pedigree file (--pedigree option) or
 the relatedness inference output from KING with the --related --degree 1 options along with age and sex information.

The pedigree file is a plain text file with header and columns: FID (family ID), IID (individual ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).
Note that individuals are assumed to have unique individual IDS (IID). Siblings are identified through individuals that have the same FID and the same FATHER_ID and MOTHER_ID.

Use the --king option to provide the KING relatedness inference output (usually has suffix .kin0) and the --agesex option to provide the age & sex information. The script
constructs a pedigree from this information and outputs it in the HDF5 output.

Args:
    -c : optional
        Duplicates offsprings of families with more than one offspring and both parents and add '_' to the start of their FIDs.
        These can be used for testing the imputation. The tests.test_imputation.imputation_test uses these.

    -silent_progress : bool, optional
        Hides the percentage of progress from logging

    -use_backup : bool, optional
        Whether it should use backup imputation where there is no ibd infomation available.

    --ibd : str
        Address of the IBD file without suffix. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome(chr_range is an optional parameters for this script).

    --ibd_is_king
        If not provided the ibd input is assumed to be in snipar. Otherwise its in king format with an allsegs file

    --bgen : str
        Address of the phased genotypes in .bgen format(should not include '.bgen'). If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome(chr_range is an optional parameters for this script).

    --bed : str
        Address of the unphased genotypes in .bed format(should not include '.bed'). If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome(chr_range is an optional parameters for this script).

    --chr_range : int, optional
        number of the chromosomes to be imputed. Should be a series of ranges with x-y format or integers.'
    
    --bim : str
        Address of a bim file containing positions of SNPs if the address is different from Bim file of genotypes

    --fam : str
        Address of a fam file containing positions of SNPs if the address is different from fam file of genotypes

    --out: str, optional
        The script writes the result of imputation to this path. If it contains '@', result of imputation for chromosome i will be written to a similar path where @ has been replaced with i. The default value for output_address is 'parent_imputed_chr'.    

    --start: int, optional
        The script can do the imputation on a slice of each chromosome. This is the start of that slice(it is inclusive).

    --end: int, optional
        The script can do the imputation on a slice of each chromosome. This is the end of that slice(it is inclusive).

    --pedigree : string, optional
        Address of the pedigree file. Pedigree file is a ' ' seperated csv with columns 'FID', 'IID', 'FATHER_ID', 'MOTHER_ID'.

    --king : string, optional
        Address of a kinship file in KING format. kinship file is a '\t' seperated csv with columns "FID1", "ID1", "FID2", "ID2, "InfType".
        Each row represents a relationship between two individuals. InfType column states the relationship between two individuals.
        The only relationships that matter for this script are full sibling and parent-offspring which are shown by 'FS' and 'PO' respectively.
        This file is used in creating a pedigree file and can be generated using KING.

    --agesex : string, optional
        Address of the agesex file. This is a " " seperated CSV with columns "FID", "IID", "FATHER_ID", "MOTHER_ID", "sex", "age".
        Each row contains the age and sex of one individual. Male and Female sex should be represented with 'M' and 'F'.
        Age column is used for distinguishing between parent and child in a parent-offsring relationship inferred from the kinship file.
        ID1 is a parent of ID2 if there is a 'PO' relationship between them and 'ID1' is at least 12 years older than ID2.

    --pcs : str, optional
        Address of the PCs file with header "FID IID PC1 PC2 ...". If provided MAFs will be estimated from PCs

    --pc_num : int, optional
        Number of PCs to consider.

    -find_optimal_pc: optional
        It will use Akaike information criterion to find the optimal number of PCs to use for MAF estimation.
    
    --threads : int, optional
        Number of the threads to be used. This should not exceed number of the available cores. The default number of the threads is one.

    --processes: int, optional
        Number of processes for imputation chromosomes. Each chromosome is done on one process.

    --chunks: int, optional
        Number of chunks load data in(each process).

    --output_compression: str, optional
        Optional compression algorithm used in writing the output as an hdf5 file. It can be either gzip or lzf.

    --output_compression_opts': int, optional
        Additional settings for the optional compression algorithm. Take a look at the create_dataset function of h5py library for more information.

    --pedigree_nan: str, optional
        The value representing NaN in the pedigreee. Default is '0'

Results:
    HDF5 files
        For each chromosome i, an HDF5 file is created at outprefix{i}. This file contains imputed genotypes, the position of SNPs, columns of resulting bim file, contents of resulting bim file, pedigree table and, family ids
        of the imputed parents, under the keys 'imputed_par_gts', 'pos', 'bim_columns', 'bim_values', 'pedigree' and, 'families', 'parental_status' respectively. There are also other details of the imputation in the resulting file. For more explanation see the documentation of snipar.imputation.impute_from_sibs.impute
        
"""
import logging
from snipar.imputation.preprocess_data import *
from snipar.imputation.impute_from_sibs import *
from snipar.utilities import parse_obsfiles
import argparse
import h5py
import random
import pandas as pd
import os
from multiprocessing import Pool
from time import time
from snipar.utilities import NumRangeAction, parseNumRange
random.seed(1567924)

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
                    Whether it should use backup imputation where there is no ibd infomation available.

                phased_address: str, optional
                    Address of the bed file (does not inlude '.bgen'). Only one of unphased_address and phased_address is neccessary.
                
                unphased_address: str, optional
                    Address of the bed file (does not inlude '.bed'). Only one of unphased_address and phased_address is neccessary.                

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
                    This function can do the imputation on a slice of each chromosome. If specified, his is the end of that slice(it is inclusive).
                    
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
    number_of_snps = len(bim)
    start_time = time()
    #Doing imputation chunk by chunk
    if start is None:
        start = 0
    if end is None:
        end = number_of_snps

    if chunks > 1:
        for i in range(chunks):
            logging.info(f"imputing chunk {i+1}/{chunks}...")
            chunk_output_address = f"{output_address}_chunk{i}"
            interval = ((end-start+chunks-1)//chunks)
            chunk_start = start+i*interval
            chunk_end = min(start+(i+1)*interval, end)
            phased_gts, unphased_gts, iid_to_bed_index, pos, freqs, hdf5_output_dict = prepare_gts(phased_address, unphased_address, bim, pedigree_output, ped_ids, chromosomes, chunk_start, chunk_end, pcs, pc_ids, find_optimal_pc)
            imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, phased_gts, unphased_gts, ibd, pos, hdf5_output_dict, str(chromosomes), freqs, chunk_output_address, threads = threads, output_compression=output_compression, output_compression_opts=output_compression_opts, silent_progress=silent_progress, use_backup=use_backup)
            logging.info(f"imputing chunk {i}/{chunks} done")
        
        #Merging chunks one by one and then removing the outputs
        logging.info("merging chunks...")
        logging.info(f"merging chunks 1/{chunks}...")
        chunk_output_address = f"{output_address}_chunk{0}.hdf5"
        hdf5_results = {}
        with h5py.File(chunk_output_address, "r") as hf:
            for key, val in hf.items():
                hdf5_results[key] = np.array(val)
        os.remove(chunk_output_address)
        for i in range(1, chunks):
            logging.info(f"merging chunks {i+1}/{chunks}...")
            chunk_output_address = f"{output_address}_chunk{i}.hdf5"
            with h5py.File(chunk_output_address, "r") as hf:
                for key in ["imputed_par_gts",
                            "pos",
                            "ratio_ibd0",
                            "mendelian_error_ratio",
                            "estimated_genotyping_error",]:
                    #TODO fix max estimator loggigng
                    new_val = np.array(hf[key])
                    hdf5_results[key] = np.hstack((hdf5_results[key], new_val))
                for key in ["bim_values"]:
                    #TODO fix max estimator loggigng
                    new_val = np.array(hf[key])
                    hdf5_results[key] = np.vstack((hdf5_results[key], new_val))
                non_duplicates = np.array(hf["non_duplicates"])
                non_duplicates = hdf5_results["non_duplicates"][-1] + non_duplicates + 1
                hdf5_results["non_duplicates"] = np.hstack((hdf5_results["non_duplicates"], non_duplicates))
            os.remove(chunk_output_address)
        logging.info(f"writing results of the merge")
        #Writing the merged output
        with h5py.File(f"{output_address}.hdf5", "w") as hf:
            for key, val in hdf5_results.items():
                if key=='imputed_par_gts':
                    hf.create_dataset(key, val.shape, dtype = 'float16', chunks = True, compression = output_compression, compression_opts=output_compression_opts, data = val)
                else:
                    hf[key] = val
        logging.info(f"merging chunks done")
    elif chunks == 1:
        phased_gts, unphased_gts, iid_to_bed_index, pos, freqs, hdf5_output_dict = prepare_gts(phased_address, unphased_address, bim, pedigree_output, ped_ids, chromosomes, start, end, pcs, pc_ids, find_optimal_pc)
        imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, phased_gts, unphased_gts, ibd, pos, hdf5_output_dict, str(chromosomes), freqs, output_address, threads = threads, output_compression=output_compression, output_compression_opts=output_compression_opts, silent_progress=silent_progress, use_backup=use_backup)
    else:
        raise Exception("invalid chunks, chunks should be a positive integer")  
    end_time = time()
    return (end_time-start_time)

#does the imputation and writes the results
def main(args):
    """"Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    if args.bgen is None and args.bed is None:
        raise Exception("You should supplement the code with at least one genotype address") 
    
    try:
        dir_name = os.path.dirname(os.path.realpath(__file__))
        import git
        repo = git.Repo(dir_name)
        logging.info(f"Last commit is: {repo.head.commit}")
        logging.info(f"summary is: {repo.head.commit.summary}")
        logging.info(f"Active branch is: {repo.active_branch.name}")        
    except :
        pass

    #fids starting with _ are reserved for control
    #Families should not have grandparents
    if not args.pedigree:
        logging.info("creating pedigree ...")
        pedigree = create_pedigree(args.king, args.agesex)
    else:
        logging.info("reading pedigree ...")
        pedigree = pd.read_csv(args.pedigree, delim_whitespace=True)
    pedigree = pedigree[['FID', 'IID', 'FATHER_ID', 'MOTHER_ID']]
    logging.info("pedigree loaded.")
    
    pcs = None
    pc_ids = None
    if args.pcs:
        logging.info("loading pcs ...")
        #ordering should be the same
        pcs = pd.read_csv(f"{args.pcs}", delim_whitespace=True)
        pc_ids = pcs.iloc[:, 1].values.astype("S").tolist()
        pcs = pcs.iloc[:, 2:].values
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

    if args.bgen and chromosomes==[None]:
        raise Exception("bgen files should be used with chr_range argument")

    if args.ibd is None:
        if args.ibd_is_king:
            raise Exception("can not use 'ibd_is_king' without any ibd file")

    def none_tansform(a, b, c):
        if a is not None:
            return a.replace(b, c)
        return None
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
            "output_address":none_tansform(args.out, "@", str(chromosome)),
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
                    help = "Whether it should use backup imputation where there is no ibd infomation available")                    
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
                    help='Start index of SNPs to perform imputation for in genotype file (starting at zero). Should be used with end parameter.',
                    default=None)
parser.add_argument('--end',
                    type=int,
                    help='End index of SNPs to perform imputation for in genotype file (goes from 0 to (end-1). Should be used with start parameter.',
                    default=None)
parser.add_argument('--pedigree',
                    type=str,
                    default = None,
                    help='Pedigree file with siblings sharing family ID')
parser.add_argument('--king',
                    type=str,
                    default = None,
                    help='Address of the king file')
parser.add_argument('--agesex',
                    type=str,
                    default = None,
                    help='Address of the agesex file with header "FID IID age sex"')
parser.add_argument('--pcs',
                    type=str,
                    default = None,
                    help='Address of the PCs file with header "FID IID PC1 PC2 ...". The IID ordering should be the same as fam file. If provided MAFs will be estimated from PCs')
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
                    help='Number of the cores to be used')
parser.add_argument('--processes',
                    type=int,
                    default=1,
                    help='Number of processes for imputation chromosomes. Each chromosome is done on one process.')
parser.add_argument('--chunks',
                    type=int,
                    default=1,
                    help='Number of chunks in each process')
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
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
    args=parser.parse_args()
    main(args)
