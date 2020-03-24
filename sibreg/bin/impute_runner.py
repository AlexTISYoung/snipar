"""Runs the sib-regression for the specified source and writes the result in an HDF5 file.

This scripts uses genotype bed files and KING IBD output with a form of pedigree file and imputes the sum of parents for families with more than one sibling and no parents.
From options --pedigree and --agesex,--king one should be given to this script because the script needs the pedigree file
and if it's not present, it has to construct one.

Parameters
----------
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

Results
-------
HDF5 files
    For each chromosome i, an HDF5 file is created at outprefix{i}. This file contains imputed genotypes, the position of SNPs, SNP ids, pedigree table and, family ids
    of the imputed parents, under the keys 'imputed_par_gts', 'pos', 'sid', 'pedigree' and, 'families' respectively.
"""

import logging
from impute_from_sibs import *
import time
import argparse
import h5py
import random
import pandas as pd
import os
from preprocess_data import *
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

    args=parser.parse_args()
    #fids starting with _ are reserved for control
    #Families should not have grandparents
    if not args.pedigree:
        logging.info("creating pedigree ...")
        pedigree = create_pedigree(args.king, args.agesex)
    else:
        pedigree = pd.read_csv(args.pedigree, sep = " ")

    if args.c:
        logging.info("Adding control to the pedigree ...")
        pedigree = add_control(pedigree)

    consumed_time = 0
    ibd = pd.read_csv(args.ibd, sep = "\t")
    for chromosome in range(args.from_chr, args.to_chr):
        print(chromosome, " is chromosome")
        sibships, iid_to_bed_index, gts, ibd, pos, hdf5_output_dict = prepare_data(pedigree, args.genotypes_prefix+str(chromosome), ibd, chromosome, args.start, args.end, args.bim)
        gts = gts.astype(float)
        pos = pos.astype(int)
        start_time = time.time()
        if args.out_prefix is None:
            imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos, hdf5_output_dict, "test_data/parent_imputed_chr"+str(chromosome))
        else:
            imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos, hdf5_output_dict, args.out_prefix+str(chromosome))
        end_time = time.time()
        consumed_time += (end_time-start_time)

    logging.info("imputation time: "+str(consumed_time))
