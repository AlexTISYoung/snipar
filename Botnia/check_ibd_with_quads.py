#!/usr/bin/env python
import argparse, code, gzip
from numba import set_num_threads
from numba import config as numba_config
import snipar.ibd
import numpy as np
from snipar.errors import estimate_genotyping_error_rate
from snipar.utilities import *
from snipar.pedigree import *

def compute_ibd_matrix(true_ibdfile,ibdfile):
    true_ibdf = gzip.open(true_ibdfile, 'r')
    ibd_header = true_ibdf.readline()
    true_ibd_snps = np.array(str(ibd_header).split(' '))
    true_ibd_snps = true_ibd_snps[2:true_ibd_snps.shape[0]]
    true_ibd_snps[true_ibd_snps.shape[0] - 1] = true_ibd_snps[true_ibd_snps.shape[0] - 1].split('\\')[0]
    ibdcols = tuple(np.arange(2, 2 + true_ibd_snps.shape[0]))
    true_ibd = np.loadtxt(true_ibdfile, usecols=ibdcols, skiprows=1, dtype=float)
    true_ibd_sibs = np.loadtxt(true_ibdfile, usecols=(0, 1), skiprows=1, dtype=str)
    ## Read inferred IBD ##
    ibdf = gzip.open(ibdfile, 'r')
    ibd_header = ibdf.readline()
    ibd_snps = np.array(str(ibd_header).split(' '))
    ibd_snps = ibd_snps[2:ibd_snps.shape[0]]
    ibd_snps[ibd_snps.shape[0] - 1] = ibd_snps[ibd_snps.shape[0] - 1].split('\\')[0]
    ibdcols = tuple(np.arange(2, 2 + ibd_snps.shape[0]))
    ibd = np.loadtxt(ibdfile, usecols=ibdcols, skiprows=1, dtype=int)
    ibd_sibs = np.loadtxt(ibdfile, usecols=(0, 1), skiprows=1, dtype=str)
    ##### Compare inferred to true IBD #####
    # Match SNPs
    sid_dict = make_id_dict(true_ibd_snps)
    true_ibd = true_ibd[:, [sid_dict[x] for x in ibd_snps]]
    # Match sibpairs
    sibpair_indices = []
    missing_pairs = []
    for i in range(true_ibd_sibs.shape[0]):
        sibs = true_ibd_sibs[i, :]
        sibpair_index = np.where(np.logical_and(ibd_sibs[:, 0] == sibs[0], ibd_sibs[:, 1] == sibs[1]))[0]
        if len(sibpair_index) == 0:
            sibpair_index = np.where(np.logical_and(ibd_sibs[:, 0] == sibs[1], ibd_sibs[:, 1] == sibs[0]))[0]
        if len(sibpair_index) == 1:
            sibpair_indices.append(sibpair_index[0])
        else:
            missing_pairs.append(i)
    ibd = ibd[sibpair_indices,:]
    if len(missing_pairs) > 0:
        true_ibd = np.delete(true_ibd,missing_pairs,0)
    # Get correct %
    ibd_matrix = np.zeros((3, 3))
    for i in range(0, 3):
        for j in range(0, 3):
            ibd_matrix[i, j] = np.sum(ibd[true_ibd == i] == j)
    return ibd_matrix

parser = argparse.ArgumentParser()
parser.add_argument('ibd',type=str,help='Address of ibd.segments.gz files. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--bgen',
                    type=str,help='Address of the phased genotypes in .bgen format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--bed',
                    type=str,help='Address of the unphased genotypes in .bed format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--chr_range',
                    type=parseNumRange,
                    nargs='*',
                    action=NumRangeAction,
                    help='number of the chromosomes to be imputed. Should be a series of ranges with x-y format or integers.', default=None)
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
parser.add_argument('--out',
                    type=str,
                    default = 'ibd',
                    help="The IBD segments will output to this path, one file for each chromosome. If the path contains '#', the '#' will be replaced with the chromosome number. Otherwise, the segments will be output to the given path with file names chr_1.ibd.segments.gz, chr_2.segments.gz, etc.")
parser.add_argument('--threads',type=int,help='Number of threads to use for IBD inference. Uses all available by default.',default=None)
args = parser.parse_args()

if args.threads is not None:
    if args.threads < numba_config.NUMBA_NUM_THREADS:
        set_num_threads(args.threads)
        print('Number of threads: '+str(args.threads))

# Check arguments
if args.bed is None and args.bgen is None:
    raise(ValueError('Must provide one of --bedfiles and --bgenfiles'))
if args.bed is not None and args.bgen is not None:
    raise(ValueError('Both bedfiles and bgenfiles provided. Please provide only one'))

# Find bed files
if args.bed is not None:
    bedfiles, chroms = parse_obsfiles(args.bed, 'bed', chromosomes=args.chr_range)
    bgenfiles = [None for x in range(chroms.shape[0])]
elif args.bgen is not None:
    bgenfiles, chroms = parse_obsfiles(args.bgen, 'bgen', chromosomes=args.chr_range)
    bedfiles = [None for x in range(chroms.shape[0])]
if args.chrom is not None and chroms.shape[0]==1:
    chroms = np.array([args.chrom],dtype=int)
# Set parameters
min_length = args.min_length
kinfile = args.king

if 0 <= args.min_maf <= 0.5:
    min_maf = args.min_maf
else:
    raise(ValueError('Min MAF must be between 0 and 0.5'))

mapfile = args.map
outprefix = args.out

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
    sibpairs, ped = get_sibpairs_from_ped(ped)
    if sibpairs is None:
        raise(ValueError('No sibpairs found'))
elif kinfile is not None:
    print('Reading relationships from '+str(kinfile))
    sibpairs = get_sibpairs_from_king(kinfile)
else:
    raise(ValueError('Must provide either KING kinship file or pedigree'))

if sibpairs.shape[0]==0:
    raise(ValueError('No sibling pairs found'))
print('Found '+str(sibpairs.shape[0])+' full sibling pairs')

#### Get genotyping error probability ####
if args.p_error is None:
    print('No genotyping error probability provided. Will attempt to estimate from parent-offspring pairs.')
    if args.pedigree is None:
        if args.agesex is not None:
            print('Constructing pedigree from '+str(kinfile)+' and age and sex information from '+str(args.agesex))
            ped = np.array(create_pedigree(kinfile, args.agesex), dtype=str)
        else:
            raise(ValueError('Must provide age and sex information (--agesex) in addition to KING kinship file, if estimating genotyping error probability'))
    if args.bed is not None:
        error_prob, error_probs = estimate_genotyping_error_rate(ped, bedfiles=bedfiles, min_maf=min_maf)
    elif args.bgen:
        error_prob, error_probs = estimate_genotyping_error_rate(ped, bgenfiles=bgenfiles, min_maf=min_maf)
    print('Estimated mean genotyping error probability: '+str(round(error_prob, 6)))
    if error_prob > 0.01:
        print('Warning: high genotyping error rate detected. Check pedigree and/or genotype data.')
else:
    error_prob = args.p_error
    error_probs = None

#ibd_matrix = compute_ibd_matrix('quad_ibd/true/chr_22.ibd.gz','ibd/22_from_ped.ibd.gz')

ibd_matrix = np.zeros((3,3))
for chr in range(1,23):
    true_ibdfile = 'quad_ibd/true/chr_'+str(chr)+'.ibd.gz'
    ibdfile = 'quad_ibd/gloabl_p_estimated/chr_'+str(chr)+'.ibd.segments.gz'
    ibd_matrix += compute_ibd_matrix(true_ibdfile,ibdfile)
    print('Done chromosome '+str(chr))

pc_correct = np.zeros((3,3))
for i in range(0,3):
    pc_correct[i,:] = ibd_matrix[i,:]/np.sum(ibd_matrix[i,:])

prob_correct = np.sum(np.array([0.25,0.5,0.25]*np.diag(pc_correct)))

np.savetxt('quad_ibd/p_estimated/ibd_matrix.txt',pc_correct)