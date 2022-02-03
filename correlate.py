import snipar.preprocess as preprocess
import numpy as np
import argparse
from snipar.correlate import *

parser = argparse.ArgumentParser()
parser.add_argument('sumstats', type=str, help='Address of sumstats files in SNIPar sumstats.gz text format (without .sumstats.gz suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22')
parser.add_argument('outprefix',type=str,help='Prefix for output file(s)')
parser.add_argument('--ldscores',type=str,help='Address of ldscores as output by LDSC',default=None)
parser.add_argument('--bedfiles', type=str,
                    help='Address of observed genotype files in .bed format (without .bed suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22.',
                    default=None)
parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.05)', default=0.05)
parser.add_argument('--corr_filter',type=float,help='Filter out SNPs with outlying sampling correlations more than corr_filter SDs from mean (default 6)',default=6.0)
parser.add_argument('--n_blocks',type=int,help='Number of blocks to use for block-jacknife variance estimate (default 200)',default=200)
parser.add_argument('--save_delete',action='store_true',help='Save jacknife delete values',default=False)
parser.add_argument('--ld_wind',type=float,help='The window, in cM, within which LD scores are computed (default 1cM)',default=1.0)
parser.add_argument('--ld_out',type=str,help='Output LD scores in LDSC format to this address',default=None)
args=parser.parse_args()

if args.ldscores is None and args.bedfiles is None:
    raise(ValueError('Must provide either LD scores or genotypes to compute ld scores'))
if args.ldscores is not None and args.bedfiles is not None:
    raise(ValueError('Both LD scores and genotypes provided. Provided LD scores will be used'))

# Find sumstats files
sumstats_files, chroms = preprocess.parse_obsfiles(args.sumstats, obsformat='sumstats.gz')
# Read sumstats
s = read_sumstats_files(sumstats_files, chroms)

# Filter
print('Filtering on missing values')
s.filter_NAs()
print('Filtering out SNPs with MAF<'+str(args.min_maf))
s.filter_maf(args.min_maf)
print('Filtering out SNPs with sampling correlations more than '+str(args.corr_filter)+' SDs from mean')
s.filter_corrs(args.corr_filter)

# Get LD scores
if args.ldscores is not None:
    ld_files, chroms = preprocess.parse_obsfiles(args.ldscores, obsformat='l2.ldscore.gz')
    s.scores_from_ldsc(ld_files)
    s.filter_NAs()
elif args.bedfiles is not None:
    bedfiles, chroms = preprocess.parse_obsfiles(args.bedfiles, obsformat='bed')
    s.compute_ld_scores(bedfiles, chroms, args.ld_wind, args.ld_out)
    s.filter_NAs()

# Compute correlations 
print('Using '+str(s.sid.shape[0])+' SNPs to compute correlations')
r_dir_pop, r_dir_pop_SE, r_dir_pop_delete = s.cor_direct_pop(args.n_blocks)
print('Correlation between direct and population effects: '+str(round(r_dir_pop,4))+' (S.E. '+str(round(r_dir_pop_SE,4))+')')
r_dir_avg_par, r_dir_avg_par_SE, r_dir_avg_par_delete = s.cor_direct_avg_par(args.n_blocks)
print('Correlation between direct and average parental effects: '+str(round(r_dir_avg_par,4))+' (S.E. '+str(round(r_dir_avg_par_SE,4))+')')

outfile = str(args.outprefix)+'_corrs.txt'
print('Saving correlation estimates to '+outfile)
first_col = np.array(['r_direct_population','r_direct_avg_parental']).reshape((2,1))
header = np.array(['correlation','est','SE']).reshape((1,3))
cor_out = np.zeros((2,2))
cor_out[0,:] = np.array([r_dir_pop,r_dir_pop_SE])
cor_out[1,:] = np.array([r_dir_avg_par,r_dir_avg_par_SE])
outarray = np.vstack((header,np.hstack((first_col,cor_out))))
np.savetxt(outfile,outarray,fmt='%s')

if args.save_delete:
    delete_outfile = args.outprefix+'_delete.txt'
    print('Saving jacknife delete values to '+str(delete_outfile))
    delete_out = np.vstack((r_dir_pop_delete,r_dir_avg_par_delete)).T
    delete_out = np.vstack((np.array(['r_direct_population','r_direct_avg_parental']).reshape((1,2)),delete_out))
    np.savetxt(delete_outfile, delete_out,fmt='%s')