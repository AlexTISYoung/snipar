import logging
# from impute_from_sibs import *
from cython_impute_from_sibs import *
import time
import argparse
import h5py
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
parser = argparse.ArgumentParser()
parser.add_argument('chr',type=int,help='Which chromosome (integer)')
parser.add_argument('ibd',type=str,help='IBD file')
parser.add_argument('genotypes',type=str,help='Genotypes in .bed format')
parser.add_argument('ped',type=str,help='Pedigree file with siblings sharing family ID')
parser.add_argument('out',type=str,help='Prefix of hdf5 output of imputed parental genotypes')
parser.add_argument('--bim',type=str,default = None, help='Bim file giving positions of SNPs in KING IBD file if different from Bim file of genotypes')
parser.add_argument('--start', type=int,
                    help='Start index of SNPs to perform imputation for in genotype file (starting at zero)',
                    default=0)
parser.add_argument('--end', type=int,
                    help='End index of SNPs to perform imputation for in genotype file (goes from 0 to (end-1)',
                    default=None)
args=parser.parse_args()
sibships, iid_to_bed_index, gts, ibd, pos = prepare_data(args.ped, args.genotypes, args.ibd, args.chr, args.start, args.end, args.bim)
start_time = time.time()
imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos)
end_time = time.time()
logging.info("imputation time: "+str(end_time-start_time))

with h5py.File('test_data/tmp.hdf5','w') as f:
    f.create_dataset('gts',imputed_par_gts.shape,dtype = 'f',chunks = True, compression = 'gzip', compression_opts=9, data = imputed_par_gts)
    f['fids'] = imputed_fids

# with h5py.File('test_data/tmp.hdf5','r') as f:
#     gts = np.array(f["gts"])
#     fids = np.array(f["fids"])