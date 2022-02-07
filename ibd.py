import argparse
from numba import set_num_threads
from numba import config as numba_config
import snipar.preprocess as preprocess
import snipar.ibd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('bedfiles', type=str,
                    help='Address of observed genotype files in .bed format (without .bed suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22.',
                    default=None)
parser.add_argument('--king',
                    type=str,
                    default = None,
                    help='Address of the king file')
parser.add_argument('--agesex',
                    type=str,
                    default=None,
                    help='Address of file with age and sex information')
parser.add_argument('--pedigree',
                    type=str,
                    default=None,
                    help='Address of pedigree file')
parser.add_argument('--map',type=str,default=None)
parser.add_argument('--outprefix',
                    type=str,
                    default = 'ibd',
                    help="Writes the result of IBD inference to outprefix.ibd.segments.gz")
parser.add_argument('--p_error',type=float,help='Probability of genotyping error. By default, this is estimated from genotyped parent-offspring pairs.',default=None)
parser.add_argument('--min_length',type=float,help='Smooth segments with length less than min_length (cM)',
                    default=0.01)
parser.add_argument('--threads',type=int,help='Number of threads to use for IBD inference. Uses all available by default.',default=None)
parser.add_argument('--min_maf',type=float,help='Minimum minor allele frequency',default=0.01)
parser.add_argument('--max_missing', type=float,
                    help='Ignore SNPs with greater percent missing calls than max_missing (default 5)', default=5)
parser.add_argument('--max_error', type=float, help='Maximum per-SNP genotyping error probability', default=0.01)
parser.add_argument('--ibdmatrix',action='store_true',default=False,help='Output a matrix of SNP IBD states (in addition to segments file)')
parser.add_argument('--ld_out',action='store_true',default=False,help='Output LD scores of SNPs (used internally for weighting).')
args = parser.parse_args()

# Find bed files
bedfiles, chroms = preprocess.parse_obsfiles(args.bedfiles, obsformat='bed')

# Set parameters
min_length = args.min_length
kinfile = args.king

if 0 <= args.min_maf <= 0.5:
    min_maf = args.min_maf
else:
    raise(ValueError('Min MAF must be between 0 and 0.5'))

mapfile = args.map
outprefix = args.outprefix

if 0 <= args.max_missing <= 100:
    max_missing = args.max_missing
else:
    raise(ValueError('Max missing % must be between 0 and 100'))

if 0 <= args.max_error <= 1:
    max_error = args.max_error
else:
    raise(ValueError('Max missing % must be between 0 and 100'))

#### Find sibling pairs ####
if args.pedigree is not None:
    print('Reading pedigree from '+str(args.pedigree))
    ped = np.loadtxt(args.pedigree,dtype=str)
    if ped.shape[1] < 4:
        raise(ValueError('Not enough columns in pedigree file'))
    elif ped.shape[1] > 4:
        print('Warning: pedigree file has more than 4 columns. The first four columns only will be used')
    # Remove rows with missing parents
    sibpairs, ped = preprocess.get_sibpairs_from_ped(ped)
    if sibpairs is None:
        raise(ValueError('No sibpairs found'))
elif kinfile is not None:
    print('Reading relationships from '+str(kinfile))
    sibpairs = preprocess.get_sibpairs_from_king(kinfile)
else:
    raise(ValueError('Must provide either KING kinship file or pedigree'))

if sibpairs.shape[0]==0:
    raise(ValueError('No sibling pairs found'))
print('Found '+str(sibpairs.shape[0])+' full sibling pairs')

# Set number of threads
if args.threads is not None:
    if args.threads < numba_config.NUMBA_NUM_THREADS:
        set_num_threads(args.threads)
        print('Number of threads: '+str(args.threads))

#### Get genotyping error probability ####
if args.p_error is None:
    print('No genotyping error probability provided. Will attempt to estimate from parent-offspring pairs.')
    if args.pedigree is None:
        if args.agesex is not None:
            print('Constructing pedigree from '+str(kinfile)+' and age and sex information from '+str(args.agesex))
            ped = np.array(preprocess.create_pedigree(kinfile, args.agesex), dtype=str)
        else:
            raise(ValueError('Must provide age and sex information (--agesex) in addition to KING kinship file, if estimating genotyping error probability'))
    else:
        ped = np.loadtxt(args.pedigree, dtype=str)
    error_prob, error_probs = preprocess.estimate_genotyping_error_rate(bedfiles, ped, min_maf)
    print('Estimated mean genotyping error probability: '+str(round(error_prob, 6)))
    if error_prob > 0.01:
        print('Warning: high genotyping error rate detected. Check pedigree and/or genotype data.')
else:
    error_prob = args.p_error
    error_probs = None

######### Infer IBD ###########
for i in range(bedfiles.shape[0]):
    if error_probs is None:
        error_probs_i = None
    else:
        error_probs_i = error_probs[i]
    snipar.ibd.infer_ibd_chr(bedfiles[i], sibpairs, error_prob, error_probs_i, outprefix,
                             min_length=min_length, mapfile=args.map,
                             ibdmatrix=args.ibdmatrix, ld_out=args.ld_out,
                             min_maf=min_maf, max_missing=max_missing, max_error=max_error)