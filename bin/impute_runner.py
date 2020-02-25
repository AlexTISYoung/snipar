import logging
# from impute_from_sibs import *
from cython_impute_from_sibs import *
import time
import argparse
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
parser = argparse.ArgumentParser()
parser.add_argument('chr',type=int,help='Which chromosome (integer)')
parser.add_argument('ibd',type=str,help='IBD file')
parser.add_argument('genotypes',type=str,help='Genotypes in .bed format')
parser.add_argument('ped',type=str,help='Pedigree file with siblings sharing family ID')
parser.add_argument('out',type=str,help='Prefix of hdf5 output of imputed parental genotypes')
parser.add_argument('--king',action='store_true',default=False,help='IBD segments file in KING format (default 23andme)')
parser.add_argument('--bim',type=str,default = None, help='Bim file giving positions of SNPs in KING IBD file if different from Bim file of genotypes')
parser.add_argument('--start', type=int,
                    help='Start index of SNPs to perform imputation for in genotype file (starting at zero)',
                    default=0)
parser.add_argument('--end', type=int,
                    help='End index of SNPs to perform imputation for in genotype file (goes from 0 to (end-1)',
                    default=None)
args=parser.parse_args()
sibships, iid_to_bed_index, gts, ibd, pos = prepare_data(args)
start_time = time.time()
imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos)
end_time = time.time()
logging.info("imputation time: "+str(end_time-start_time))
# with open(args.out, "w") as f:
#     json.dump(imputed_par_gts , f)

# [16057417, 37309959, 1]
# ('4834194', '5636170')
# (('4834194', '5636170'), <MemoryView of 'NoneType' at 0x7fd5c7664e50>)
# 24986

# 5783668

# within function
# 16057417
# -2016694208
# in impute


# in dint to map
# within function
# 16057417
# -685801696
# in impute

# within function
# 16057417
# 1651884080
# in impute
# Bus error (core dumped)


# within function
# 16057417
# 2038592000
# in impute