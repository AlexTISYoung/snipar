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
    parser.add_argument('--bim',type=str,default = None, help='Bim file giving positions of SNPs in KING IBD file if different from Bim file of genotypes')
    parser.add_argument('--out_prefix',type=str,default = None, help="Writes the result of imputation for chromosome i to outprefix{i}")
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
