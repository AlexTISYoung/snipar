"""Runs the sib-regression for the specified source and writes the result in an HDF5 file.

This script performs imputation of missing parental genotypes for families with at least two siblings and no genotyped parents. It takes as input the genotypes
 of the siblings in a .bed file and the IBD segments output by KING (with --ibdsegs option). To specify the siblings, one can either provide a pedigree file (--pedigree option) or
 the relatedness inference output from KING with the --related --degree 1 options along with age and sex information.

The pedigree file is a plain text file with header and columns: FID (family ID), IID (individual ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).
Note that individuals are assumed to have unique individual IDS (IID). Siblings are identified through individuals that have the same FID and the same FATHER_ID and MOTHER_ID.

Use the --king option to provide the KING relatedness inference output (usually has suffix .kin0) and the --agesex option to provide the age & sex information. The script
constructs a pedigree from this information and outputs it in the HDF5 output.

Args:
    -c : optional
        Duplicates offsprings of families with more than one offspring and both parents and add '_' to the start of their FIDs.
        These can be used for testing the imputation. The tests.test_imputation.imputation_test uses these.

    from_chr : int
        Does the imputation for chromosomes with chromosome number bigger than or equal to this.

    to_chr : int
        Does the imputation for chromosomes with chromosome number less than this.

    IBD : str
            Address of a file containing IBD statuses for all SNPs.
        This is a '\t seperated CSV with these columns: "chr", "ID1", "ID2", "IBDType", "StartSNP", "StopSNP".
        Each line states an IBD segment between a pair on individuals. This can be generated using King software

    genotypes_prefix : str
        Prefix for the address of genotype bed files. Address of the bed file for chromosome i would be genotypes_prefix{i}.

    bim : str
        Address of a bim file containing positions of SNPs if the address is different from Bim file of genotypes

    --out_prefix: str, optional
        The script writes the result of imputation for chromosome i to outprefix{i}. the default value for out_prefix is 'parent_imputed_chr'.

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

Results:
    HDF5 files
        For each chromosome i, an HDF5 file is created at outprefix{i}. This file contains imputed genotypes, the position of SNPs, SNP ids, pedigree table and, family ids
        of the imputed parents, under the keys 'imputed_par_gts', 'pos', 'sid', 'pedigree' and, 'families', 'parental_status' respectively.
        
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
random.seed(1567924)

#does the imputation and writes the results
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', action='store_true')
    parser.add_argument('from_chr',type=int,help='Which chromosome (<=)')
    parser.add_argument('to_chr',type=int,help='Which chromosome (<)')
    parser.add_argument('ibd',type=str,help='IBD file')
    parser.add_argument('genotypes_prefix',type=str,help='prefix of genotypes in .bed format')
    parser.add_argument('--bim',type=str,default = None, help='Address of a bim file containing positions of SNPs if the address is different from Bim file of genotypes')
    parser.add_argument('--out_prefix',type=str,default = "parent_imputed_chr", help="Writes the result of imputation for chromosome i to outprefix{i}")
    parser.add_argument('--start', type=int,
                        help='Start index of SNPs to perform imputation for in genotype file (starting at zero)',
                        default=None)
    parser.add_argument('--end', type=int,
                        help='End index of SNPs to perform imputation for in genotype file (goes from 0 to (end-1)',
                        default=None)
    parser.add_argument('--pedigree',type=str,default = None, help='Pedigree file with siblings sharing family ID')
    parser.add_argument('--king',type=str,default = None, help='Address of the king file')
    parser.add_argument('--agesex',type=str,default = None, help='Address of the agesex file with header "FID IID age sex"')
    parser.add_argument('--threads', type=int, default=1, help='Number of the cores to be used')

    args=parser.parse_args()
    #fids starting with _ are reserved for control
    #Families should not have grandparents
    if not args.pedigree:
        logging.info("creating pedigree ...")
        pedigree = create_pedigree(args.king, args.agesex)
    else:
        pedigree = pd.read_csv(args.pedigree, sep = " ")
    logging.info("pedigree loaded.")

    if args.c:
        logging.info("Adding control to the pedigree ...")
        pedigree = add_control(pedigree)
        logging.info("Control Added.")
    

    consumed_time = 0
    logging.info("Loading ibd ...")
    ibd_pd = pd.read_csv(args.ibd, sep = "\t")
    logging.info("ibd loaded.")
    for chromosome in range(args.from_chr, args.to_chr):
        logging.info(str(chromosome) + " is chromosome")
        bed_address = args.genotypes_prefix
        if "~" in bed_address:
            bed_address = bed_address.replace("~", str(chromosome))
        sibships, iid_to_bed_index, gts, ibd, pos, hdf5_output_dict = prepare_data(pedigree, bed_address, ibd_pd, chromosome, args.start, args.end, args.bim)
        gts = gts.astype(float)
        pos = pos.astype(int)
        start_time = time.time()
        address = args.out_prefix
        if address is None:
            file_address = "outputs/parent_imputed_chr"+str(chromosome)
        else:
            file_address = address+str(chromosome)
        imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos, hdf5_output_dict, file_address, threads=args.threads)
        end_time = time.time()
        consumed_time += (end_time-start_time)
    logging.info("imputation time: "+str(consumed_time))
