#!/usr/bin/env python
"""Infers identity-by-descent (IBD) segments shared between full-siblings.

Minimally: the script requires observed sibling genotypes in either .bed or .bgen format, along with information
on the relations present in the dataset, which can be provided using a pedigree file or the results
of KING kinship inference along with age and sex information (from which a pedigree can be constructed).

**Usage**

    ``ibd.py --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ld_out``

Args:
@parser@

Results:
    IBD segments
        For each chromosome, a gzipped text file containing the IBD segments for the siblings is output. 
        
"""
import argparse
from numba import set_num_threads
from numba import config as numba_config
import snipar.ibd
import numpy as np
from snipar.errors import estimate_genotyping_error_rate
from snipar.utilities import *
from snipar.pedigree import *
from snipar.numrange import parseNumRange, NumRangeAction
parser = argparse.ArgumentParser()
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
parser.add_argument('--map', type=str, help='Path to a SHAPEIT-format genetic map with pposition and gposition columns.')
parser.add_argument('--out',
                    type=str,
                    default = 'ibd',
                    help="Output prefix for per-chromosome IBD segments. An @ in the prefix is replaced with the chromosome number.")
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
parser.add_argument('--chrom',type=int,help='The chromosome of the input .bgen file. Helpful if inputting a single .bgen file without chromosome information.',default=None)
parser.add_argument('--batches',type=int,help='Number of batches to split the data (by sibpair) into for IBD inference. Useful for large datasets.',default=1)
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))
def main(args):
    """"Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    # Set number of threads
    if args.threads is not None:
        num_threads = min([args.threads, numba_config.NUMBA_NUM_THREADS])
    else:
        num_threads = numba_config.NUMBA_NUM_THREADS
    set_num_threads(num_threads)
    print('Number of threads: '+str(num_threads))

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
    if args.batches < 1:
        raise ValueError("--batches must be positive")
    min_length = args.min_length
    kinfile = args.king

    if 0 <= args.min_maf <= 0.5:
        min_maf = args.min_maf
    else:
        raise(ValueError('Min MAF must be between 0 and 0.5'))

    if 0 <= args.max_missing <= 100:
        max_missing = args.max_missing
    else:
        raise(ValueError('Max missing % must be between 0 and 100'))

    if 0 <= args.max_error <= 1:
        max_error = args.max_error
    else:
        raise ValueError('Max error probability must be between 0 and 1')

    #### Find sibling pairs ####
    if args.pedigree is not None:
        print('Reading pedigree from '+str(args.pedigree))
        ped = np.loadtxt(args.pedigree, dtype=str, ndmin=2)
        if ped.size and tuple(ped[0, :4]) == ("FID", "IID", "FATHER_ID", "MOTHER_ID"):
            ped = ped[1:]
        if not ped.size:
            raise ValueError("No sibling pairs found")
        if ped.shape[1] < 4:
            raise(ValueError('Not enough columns in pedigree file'))
        elif ped.shape[1] > 4:
            print('Warning: pedigree file has more than 4 columns. The first four columns only will be used')
        # Remove rows with missing parents
        sibpairs, _ = get_sibpairs_from_ped(ped)
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
        if not 0 <= error_prob <= 1:
            raise ValueError("--p_error must be between 0 and 1")
        error_probs = None

    batches = np.array_split(sibpairs, min(args.batches, len(sibpairs)))
    for i, chrom in enumerate(chroms):
        error_probs_i = None if error_probs is None else error_probs[i]
        for j, batch in enumerate(batches, 1):
            outprefix = args.out
            if args.batches > 1:
                print(f"Inferring IBD on batch {j} of {len(batches)}")
                outprefix += str(j)
            snipar.ibd.infer_ibd_chr(
                batch, error_prob, error_probs_i, outprefix, bedfile=bedfiles[i],
                bgenfile=bgenfiles[i], chrom=None if chrom == 0 else chrom,
                min_length=min_length, mapfile=args.map, ibdmatrix=args.ibdmatrix,
                ld_out=args.ld_out, min_maf=min_maf, max_missing=max_missing,
                max_error=max_error,
            )

if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
