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

    IBD : str
        Address of a file containing IBD statuses for all SNPs.
        This is a '\t seperated CSV with these columns: "chr", "ID1", "ID2", "IBDType", "StartSNP", "StopSNP".
        Each line states an IBD segment between a pair on individuals. This can be generated using King software

    --bgen : str
        Address of the phased genotypes in .bgen format(should not include '.bgen'). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of [from_chr, to_chr) for each chromosome(from_chr and to_chr are two optional parameters for this script).

    --bed : str
        Address of the unphased genotypes in .bed format(should not include '.bed'). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of [from_chr, to_chr) for each chromosome(from_chr and to_chr are two optional parameters for this script).

    bim : str
        Address of a bim file containing positions of SNPs if the address is different from Bim file of genotypes

    --from_chr : int, optional
        Which chromosome (<=). Should be used with to_chr parameter.
        

    --to_chr : int, optional
        Which chromosome (<). Should be used with from_chr parameter.

    --output_address: str, optional
        The script writes the result of imputation to this path. If it contains '~', result of imputation for chromosome i will be written to a similar path where ~ has been replaced with i. The default value for output_address is 'parent_imputed_chr'.

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

Results:
    HDF5 files
        For each chromosome i, an HDF5 file is created at outprefix{i}. This file contains imputed genotypes, the position of SNPs, columns of resulting bim file, contents of resulting bim file, pedigree table and, family ids
        of the imputed parents, under the keys 'imputed_par_gts', 'pos', 'bim_columns', 'bim_values', 'pedigree' and, 'families', 'parental_status' respectively.
        
"""
import logging
from sibreg.bin.preprocess_data import *
from sibreg.bin.impute_from_sibs import *
import time
import argparse
import h5py
import random
import pandas as pd
import os
from multiprocessing import Pool
random.seed(1567924)

def run_imputation(data):
    """Runs the imputation and returns the consumed time
    Args:
        data : dict
            a dictionary with these keys and values:
            Keys:
                pedigree: pd.Dataframe
                    The standard pedigree table

                unphased_address: str, optional
                    Address of the bed file (does not inlude '.bed'). Only one of unphased_address and phased_address is neccessary.
                
                phased_address: str, optional
                    Address of the bed file (does not inlude '.bgen'). Only one of unphased_address and phased_address is neccessary.


                ibd_pd: pd.Dataframe
                    IBD segments. Only needs to contain information about this chromosome. Should contain these columns: ID1, ID2, IBDType, Chr, start_coordinate, stop_coordinate

                output_address: str
                    The address to write the result of imputation on. The default value for output_address is 'parent_imputed_chr'.

                chunks: int
                    Number of chunks load data in(each process).

                start: int, optional
                    This function can do the imputation on a slice of each chromosome. If specified, his is the start of that slice(it is inclusive).

                end: int, optional
                    This function can do the imputation on a slice of each chromosome. If specified, his is the end of that slice(it is inclusive).
                    
                bim: str, optional
                    Address of a bim file containing positions of SNPs if the address is different from Bim file of genotypes.

                threads: int, optional
                    Number of the threads to be used. This should not exceed number of the available cores. The default number of the threads is one.

                output_compression: str, optional
                    Optional compression algorithm used in writing the output as an hdf5 file. It can be either gzip or lzf. None means no compression.

                output_compression_opts': int, optional
                    Additional settings for the optional compression algorithm. Take a look at the create_dataset function of h5py library for more information. None means no compression setting.
    Returns:
        float
            time consumed byt the imputation.
    """
    chunks = data["chunks"]
    pedigree = data["pedigree"]
    phased_address = data.get("phased_address")
    unphased_address = data.get("unphased_address")
    ibd_pd = data["ibd_pd"]
    output_address = data["output_address"]
    start = data.get("start")
    end = data.get("end")
    bim = data.get("bim")
    threads = data.get("threads")
    output_compression = data.get("output_compression")
    output_compression_opts = data.get("output_compression_opts")
    chromosome = data.get("chromosome")
    pedigree_nan = data.get("pedigree_nan")
    logging.info("processing " + str(phased_address) + "," + str(unphased_address))
    sibships, ibd, bim, chromosomes, ped_ids, pedigree_output = prepare_data(pedigree, phased_address, unphased_address, ibd_pd, bim, chromosome = chromosome, pedigree_nan=pedigree_nan)
    logging.info(f"after prepare_data ibd is {list(ibd_pd.items())[:2]}")
    number_of_snps = len(bim)
    start_time = time.time()
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
            phased_gts, unphased_gts, iid_to_bed_index, pos, freqs, hdf5_output_dict = prepare_gts(phased_address, unphased_address, bim, pedigree_output, ped_ids, chromosomes, chunk_start, chunk_end)
            imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, phased_gts, unphased_gts, ibd, pos, hdf5_output_dict, str(chromosomes), freqs, chunk_output_address, threads = threads, output_compression=output_compression, output_compression_opts=output_compression_opts)
            logging.info(f"imputing chunk {i}/{chunks} done")
        
        #Merging chunks one by one and then removing the outputs
        logging.info("merging chunks...")
        logging.info(f"merging chunks 1/{chunks}...")
        chunk_output_address = f"{output_address}_chunk{0}.hdf5"
        with h5py.File(chunk_output_address, "r") as hf:
            imputed_par_gts = np.array(hf['imputed_par_gts'])
            families = np.array(hf['families'])
            parental_status = np.array(hf['parental_status'])
            pos = np.array(hf['pos'])
            bim_columns = np.array(hf["bim_columns"])
            bim_values = np.array(hf["bim_values"])
            pedigree = np.array(hf["pedigree"])
            ratio_ibd0 = np.array(hf["ratio_ibd0"])
            non_duplicates = np.array(hf["non_duplicates"])
        os.remove(chunk_output_address)
        for i in range(1, chunks):
            logging.info(f"merging chunks {i+1}/{chunks}...")
            chunk_output_address = f"{output_address}_chunk{i}.hdf5"
            with h5py.File(chunk_output_address, "r") as hf:
                new_imputed_par_gts = np.array(hf['imputed_par_gts'])
                new_pos = np.array(hf['pos'])
                new_bim_values = np.array(hf["bim_values"])
                new_ratio_ibd0 = np.array(hf["ratio_ibd0"])
                new_non_duplicates = np.array(hf["non_duplicates"])
                new_non_duplicates = non_duplicates[-1] + new_non_duplicates + 1
                imputed_par_gts = np.hstack((imputed_par_gts, new_imputed_par_gts))
                pos = np.hstack((pos, new_pos))
                non_duplicates = np.hstack((non_duplicates, new_non_duplicates))
                bim_values = np.vstack((bim_values, new_bim_values))
                ratio_ibd0 = np.hstack((ratio_ibd0, new_ratio_ibd0))
            os.remove(chunk_output_address)
        logging.info(f"writing results of the merge")
        #Writing the merged output
        with h5py.File(f"{output_address}.hdf5", "w") as hf:
            hf.create_dataset('imputed_par_gts', imputed_par_gts.shape, dtype = 'float16', chunks = True, compression = output_compression, compression_opts=output_compression_opts, data = imputed_par_gts)
            hf['families'] = np.array(families, dtype='S')
            hf['parental_status'] = parental_status
            hf['pos'] = pos
            hf["bim_columns"] = bim_columns
            hf["bim_values"] = bim_values
            hf["pedigree"] =  pedigree
            hf["ratio_ibd0"] = ratio_ibd0
            hf["non_duplicates"] = non_duplicates
        logging.info(f"merging chunks done")
    elif chunks == 1:
        phased_gts, unphased_gts, iid_to_bed_index, pos, freqs, hdf5_output_dict = prepare_gts(phased_address, unphased_address, bim, pedigree_output, ped_ids, chromosomes, start, end)
        imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, phased_gts, unphased_gts, ibd, pos, hdf5_output_dict, str(chromosomes), freqs, output_address, threads = threads, output_compression=output_compression, output_compression_opts=output_compression_opts)
    else:
        raise Exception("invalid chunks, chunks should be a positive integer")  
    end_time = time.time()
    return (end_time-start_time)

#does the imputation and writes the results
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
    parser = argparse.ArgumentParser()
    parser.add_argument('-c',
                        action='store_true')
    parser.add_argument('ibd',
                        type=str,
                        help='IBD file, should contain these columns: ID1, ID2, IBDType, Chr, start_coordinate, stop_coordinate')                        
    parser.add_argument('--bgen',
                        type=str,help='Address of the phased genotypes in .bgen format. If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of [from_chr, to_chr) for each chromosome(from_chr and to_chr are two optional parameters for this script).')
    parser.add_argument('--bed',
                        type=str,help='Address of the unphased genotypes in .bed format. If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of [from_chr, to_chr) for each chromosome(from_chr and to_chr are two optional parameters for this script). Should be supplemented with from_chr and to_chr.')
    parser.add_argument('--from_chr',
                        type=int,                        
                        help='Which chromosome (<=). Should be used with to_chr parameter.')
    parser.add_argument('--to_chr',
                        type=int,
                        help='Which chromosome (<). Should be used with from_chr parameter.')
    parser.add_argument('--bim',
                        type=str,
                        default = None,
                        help='Address of a bim file containing positions of SNPs if the address is different from Bim file of genotypes')
    parser.add_argument('--output_address',
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

    args=parser.parse_args()
    if args.bgen is None and args.bed is None:
        raise Exception("You should supplement the code with at least one genotype address") 
        
    #fids starting with _ are reserved for control
    #Families should not have grandparents
    if not args.pedigree:
        logging.info("creating pedigree ...")
        pedigree = create_pedigree(args.king, args.agesex)
    else:
        pedigree = pd.read_csv(args.pedigree, delim_whitespace=True)
    logging.info("pedigree loaded.")

    if args.c:
        logging.info("Adding control to the pedigree ...")
        pedigree = add_control(pedigree)
        logging.info("Control Added.")
    
    logging.info("Loading ibd ...")
    ibd_pd = pd.read_csv(f"{args.ibd}.segments.gz", delim_whitespace=True).astype(str)
    print(set(ibd_pd.columns.values.tolist()))
    if not {"ID1", "ID2", "IBDType", "Chr", "start_coordinate", "stop_coordinate",}.issubset(set(ibd_pd.columns.values.tolist())):
        raise Exception("Invalid ibd columns. Columns must be: ID1, ID2, IBDType, Chr, start_coordinate, stop_coordinate")
    logging.info("ibd loaded.")
    logging.info(f"ibd is {ibd_pd}")
    logging.info(f"ibd0 is {ibd_pd[ibd_pd[IBDType]==0]}")
    if (args.from_chr is not None) and (args.to_chr is not None):
        chromosomes = [str(chromosome) for chromosome in range(args.from_chr, args.to_chr)]
    else:
        chromosomes = [None]

    if (args.bed and "~" in args.bed) or (args.bgen and "~" in args.bgen) or (args.output_address and "~" in args.output_address):
        if args.to_chr is None or args.from_chr is None:
            raise Exception("no chromosome range specified for the wildcard ~ in the address")

    if args.bgen:
        if args.to_chr is None or args.from_chr is None:
            raise Exception("Chromosome range should be specified with unphased genotype")

    def none_tansform(a, b, c):
        if a is not None:
            return a.replace(b, c)
        return None
    inputs = [{"pedigree": pedigree,
            "phased_address": none_tansform(args.bgen, "~", str(chromosome)),
            "unphased_address": none_tansform(args.bed, "~", str(chromosome)),
            "ibd_pd": ibd_pd,
            "output_address":none_tansform(args.output_address, "~", str(chromosome)),
            "start": args.start,
            "end": args.end,
            "bim": none_tansform(args.bim, "~", str(chromosome)),
            "threads": args.threads,
            "chunks": args.chunks,
            "output_compression":args.output_compression,
            "output_compression_opts":args.output_compression_opts,
            "chromosome":chromosome,
            "pedigree_nan":args.pedigree_nan,
            }
            for chromosome in chromosomes]
    #TODO output more information about the imputation inside the hdf5 filehf
    with Pool(args.processes) as pool:
        logging.info("staring process pool")
        consumed_time = pool.map(run_imputation, inputs)
        logging.info("imputation time: "+str(np.sum(consumed_time)))
