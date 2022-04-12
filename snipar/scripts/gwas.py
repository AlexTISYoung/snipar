#!/usr/bin/env python
import argparse, h5py
import snipar.read as read
import numpy as np
import snipar.lmm as lmm
from snipar.utilities import *
from snipar.gwas import *
from numba import set_num_threads
from numba import config as numba_config
from snipar.pedigree import get_sibpairs_from_ped

######### Command line arguments #########
parser=argparse.ArgumentParser()
parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
parser.add_argument('--bgen',
                    type=str,help='Address of the phased genotypes in .bgen format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--bed',
                    type=str,help='Address of the unphased genotypes in .bed format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--imp', type=str, help='Address of hdf5 files with imputed parental genotypes (without .hdf5 suffix). If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range (chr_range is an optional parameters for this script).', default = None)
parser.add_argument('--chr_range',
                    type=parseNumRange,
                    nargs='*',
                    action=NumRangeAction,
                    help='number of the chromosomes to be imputed. Should be a series of ranges with x-y format or integers.', default=None)
parser.add_argument('--out', type=str, help="The summary statistics will output to this path, one file for each chromosome. If the path contains '@', the '@' will be replaced with the chromosome number. Otherwise, the summary statistics will be output to the given path with file names chr_1.sumstats.gz, chr_2.sumstats.gz, etc. for the text sumstats, and chr_1.sumstats.hdf5, etc. for the HDF5 sumstats",default='./')
parser.add_argument('--pedigree',type=str,help='Address of pedigree file. Must be provided if not providing imputed parental genotypes.',default=None)
parser.add_argument('--parsum',action='store_true',help='Regress onto proband and sum of (imputed/observed) maternal and paternal genotypes. Default uses separate paternal and maternal genotypes when available.',default = False)
parser.add_argument('--fit_sib',action='store_true',help='Fit indirect effect from sibling ',default=False)
parser.add_argument('--covar',type=str,help='Path to file with covariates: plain text file with columns FID, IID, covar1, covar2, ..', default=None)
parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                    default=1)
parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)', default=0.01)
parser.add_argument('--threads',type=int,help='Number of threads to use for IBD inference. Uses all available by default.',default=None)
parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)', default=5)
parser.add_argument('--batch_size',type=int,help='Batch size of SNPs to load at a time (reduce to reduce memory requirements)',default=100000)
parser.add_argument('--no_hdf5_out',action='store_true',help='Suppress HDF5 output of summary statistics',default=False)
parser.add_argument('--no_txt_out',action='store_true',help='Suppress text output of summary statistics',default=False)
parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)', default='NA')
parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                    default=1)
args=parser.parse_args()

# Set number of threads
if args.threads is not None:
    if args.threads < numba_config.NUMBA_NUM_THREADS:
        set_num_threads(args.threads)
        print('Number of threads: '+str(args.threads))

# Check arguments
if args.bed is None and args.bgen is None:
    raise(ValueError('Must provide one of --bedfiles and --bgenfiles'))
if args.bed is not None and args.bgen is not None:
    raise(ValueError('Both bed files and bgen files provided. Please provide only one'))
if args.imp is None and args.pedigree is None:
    raise(ValueError('Must provide pedigree if not providing imputed parental genotypes file(s)'))

# Find observed and imputed files
if args.imp is None:
    print('Warning: no imputed parental genotypes provided. Will analyse only individuals with both parents genotyped.')
    if args.bed is not None:
        bedfiles, chroms = parse_obsfiles(args.bed, 'bed', chromosomes=args.chr_range)
        bgenfiles = [None for x in range(chroms.shape[0])]
    elif args.bgen is not None:
        bgenfiles, chroms = parse_obsfiles(args.bgen, 'bgen', chromosomes=args.chr_range)
        bedfiles = [None for x in range(chroms.shape[0])]
    pargts_list = [None for x in range(chroms.shape[0])]
else:
    if args.bed is not None:
        bedfiles, pargts_list, chroms = parse_filelist(args.bed, args.imp, 'bed', chromosomes=args.chr_range)
        bgenfiles = [None for x in range(chroms.shape[0])]
    elif args.bgen is not None:
        bgenfiles, pargts_list, chroms = parse_filelist(args.bgen, args.imp, 'bgen', chromosomes=args.chr_range)
        bedfiles = [None for x in range(chroms.shape[0])]
if chroms.shape[0]==0:
    raise(ValueError('No input genotype files found'))

# Read phenotype and covariates
######### Read Phenotype ########
y = read.phenotype.read_phenotype(args.phenofile, missing_char=args.missing_char, phen_index=args.phen_index)
######## Read covariates ########
if args.covar is not None:
    print('Reading covariates')
    covariates = read.phenotype.read_covariates(args.covar, pheno_ids=y.ids, missing_char=args.missing_char)
    # Match to pheno ids
    covariates.filter_ids(y.ids)
else:
    covariates = None

# Read pedigree
if args.imp is None:
    print('Reading pedigree from '+str(args.pedigree))
    ped = np.loadtxt(args.pedigree,dtype=str)
    if ped.shape[1] < 4:
        raise(ValueError('Not enough columns in pedigree file'))
    elif ped.shape[1] > 4:
        print('Warning: pedigree file has more than 4 columns. The first four columns only will be used')
    # Remove rows with missing parents
    sibpairs, ped = get_sibpairs_from_ped(ped)
    if sibpairs is not None:
        print('Found '+str(sibpairs.shape[0])+' sibling pairs in pedigree')
    else:
        print('Found 0 sibling pairs')
else:
    # Read pedigree
    par_gts_f = h5py.File(pargts_list[0],'r')
    ped = convert_str_array(par_gts_f['pedigree'])
    ped = ped[1:ped.shape[0]]
    # Remove control fams
    controls = np.array([x[0]=='_' for x in ped[:,0]])
    ped = ped[~controls,:]

####### Fit null model ######
# Match to pedigree
ped_dict = make_id_dict(ped,1)
y.filter_ids(ped[:,1])
print(str(y.shape[0])+' individuals with phenotype values found in pedigree')
ped_indices = np.array([ped_dict[x] for x in y.ids])
y.fams = ped[ped_indices,0]

# Fit variance components
print('Fitting variance components')
if args.covar is not None:
    # Match covariates
    covariates.filter_ids(y.ids)
    # Fit null model
    null_model, sigma2, tau, null_alpha, null_alpha_cov = lmm.fit_model(y.gts[:,0], covariates.gts, y.fams, add_intercept=True,
                                                                        tau_init=args.tau_init)
    # Adjust for covariates
    y.gts[:,0] = y.gts[:,0]-(null_alpha[0]+covariates.gts.dot(null_alpha[1:null_alpha.shape[0]]))
else:
    # Fit null model
    null_model, sigma2, tau = lmm.fit_model(y.gts[:,0], np.ones((y.shape[0], 1)), y.fams,
                                            tau_init = args.tau_init, return_fixed = False)
    y.gts[:,0] = y.gts[:,0]-np.mean(y.gts[:,0])
print('Family variance estimate: '+str(round(sigma2/tau,4)))
print('Residual variance estimate: ' + str(round(sigma2,4)))

# Diagonalize y
print('Transforming phenotype')
L = null_model.sigma_inv_root(tau, sigma2)
y.diagonalise(L)

for i in range(chroms.shape[0]):
    if args.bed is not None:
        print('Observed genotypes file: '+bedfiles[i])
    if args.bgen is not None:
        print('Observed genotypes file: '+bgenfiles[i])
    if args.imp is not None:
        print('Imputed genotypes file: '+pargts_list[i])
    if chroms.shape[0]>1:
        print('Estimating SNP effects for chromosome '+str(chroms[i]))
    else:
        print('Estimaing SNP effects')
    process_chromosome(chroms[i], y, ped, tau, sigma2, args.out, bedfile=bedfiles[i], bgenfile=bgenfiles[i], 
                        par_gts_f=pargts_list[i], fit_sib=args.fit_sib, parsum=args.parsum, 
                        max_missing=args.max_missing, min_maf=args.min_maf, batch_size=args.batch_size, 
                        no_hdf5_out=args.no_hdf5_out, no_txt_out=args.no_txt_out)