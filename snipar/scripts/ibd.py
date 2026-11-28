#!/usr/bin/env python
"""Infers identity-by-descent (IBD) segments shared between full-siblings.

The script requires observed sibling genotypes in BED, BGEN or VCF format, along with information
on the relations present in the dataset, which can be provided using a pedigree file or the results
of KING kinship inference along with age and sex information (from which a pedigree can be constructed).

**Usage**

    ``ibd.py --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ld_out``

Args:
@parser@

Results:
    IBD segments
        Gzipped hard-call segments by chromosome, with optional per-marker state probabilities in HDF5.

    Genome-wide relatedness
        Posterior expected sharing for each sibling pair, weighted separately by genetic and physical distance.
"""
import argparse
from numba import set_num_threads
from snipar.resources import available_cpus
from numba import config as numba_config
import snipar.ibd
import numpy as np
from snipar.errors import estimate_genotyping_error_rate
from snipar.ibd_output import Relatedness
from snipar.utilities import get_parser_doc, parse_obsfiles
from snipar.pedigree import create_pedigree, get_sibpairs_from_king, get_sibpairs_from_ped
from snipar.numrange import parseNumRange, NumRangeAction
from snipar.genotype_cli import add_vcf_options, with_genotypes
parser = argparse.ArgumentParser()
parser.add_argument('--bgen',
                    type=str, help='Prefix of phased or unphased .bgen genotypes. Use @ to substitute chromosome numbers.')
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
parser.add_argument('--map', type=str, help='SHAPEIT-format genetic map with pposition (bp) and gposition (cM) '
                                           'columns. SNIPar linearly extends the map at both ends.')
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
parser.add_argument('--posterior', action='store_true',
                    help='Write marginal IBD0, IBD1, and IBD2 probabilities at retained SNPs to HDF5.')
parser.add_argument('--joint', action='store_true',
                    help='Jointly model complete sets of three or four genotyped siblings, using unphased dosages.')
parser.add_argument('--ld_out',action='store_true',default=False,help='Output LD scores of SNPs (used internally for weighting).')
parser.add_argument('--chrom',type=int,help='The chromosome of the input .bgen file. Helpful if inputting a single .bgen file without chromosome information.',default=None)
parser.add_argument('--batches',type=int,help='Number of batches to split the data (by sibpair) into for IBD inference. Useful for large datasets.',default=1)
add_vcf_options(parser, bgen_quality=True)
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))

def _load_relationships(args):
    ped = None
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
    elif args.king is not None:
        print('Reading relationships from '+str(args.king))
        sibpairs = get_sibpairs_from_king(args.king)
    else:
        raise(ValueError('Must provide either KING kinship file or pedigree'))
    if sibpairs is None or not len(sibpairs):
        raise ValueError("No sibling pairs found")
    print(f"Found {len(sibpairs)} full sibling pairs")
    if args.p_error is None and ped is None:
        if args.agesex is None:
            raise ValueError("Provide age and sex information (--agesex) with KING kinship "
                             "to estimate the genotyping error probability")
        print("Constructing pedigree from " + str(args.king) + " and " + str(args.agesex))
        ped = np.array(create_pedigree(args.king, args.agesex), dtype=str)
    return sibpairs, ped


def _vcf_samples(args, available_ids):
    sibpairs, ped = _load_relationships(args)
    available = set(sibpairs.flat)
    if args.p_error is None:
        available.update(ped[:, 1:4].flat)
    available.intersection_update(available_ids)
    selected = set()
    for first, second in sibpairs:
        if first in available and second in available:
            selected.update((first, second))
    if not selected:
        raise ValueError("No genotyped sibling pairs found")
    if args.p_error is None:
        for child, father, mother in ped[:, 1:4]:
            if child in available:
                if father in available:
                    selected.update((child, father))
                if mother in available:
                    selected.update((child, mother))
    # Preserve the original pair ordering and batch boundaries. Each batch
    # still tests sample availability as for a direct BED/BGEN invocation.
    args._vcf_relationships = sibpairs, ped
    return selected


def _pair_batches(sibpairs, count, *, joint_groups=()):
    batches = np.array_split(np.arange(len(sibpairs), dtype=np.intp), min(count, len(sibpairs)))
    if not joint_groups or len(batches) == 1:
        return batches

    assignments = np.repeat(np.arange(len(batches)), [len(batch) for batch in batches])
    for rows in joint_groups:
        assignments[rows] = assignments[rows[0]]
    boundaries = np.bincount(assignments, minlength=len(batches)).cumsum()[:-1]
    order = np.argsort(assignments, kind='stable')
    return [rows for rows in np.split(order, boundaries) if len(rows)]


@with_genotypes(vcf_samples=_vcf_samples)
def main(args):
    """"Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    # Set number of threads
    if args.threads is not None:
        num_threads = min(args.threads, numba_config.NUMBA_NUM_THREADS, available_cpus())
    else:
        num_threads = min(numba_config.NUMBA_NUM_THREADS, available_cpus())
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

    relationships = getattr(args, "_vcf_relationships", None)
    sibpairs, ped = _load_relationships(args) if relationships is None else relationships

    #### Get genotyping error probability ####
    if args.p_error is None:
        print('No genotyping error probability provided. Will attempt to estimate from parent-offspring pairs.')
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

    joint_groups = []
    if getattr(args, 'joint', False):
        from snipar.ibd_joint import joint_pair_groups

        joint_groups = joint_pair_groups(sibpairs)
    joint_mask = np.zeros(len(sibpairs), dtype=bool)
    for rows in joint_groups:
        joint_mask[rows] = True
    batches = _pair_batches(sibpairs, args.batches, joint_groups=joint_groups)
    relatedness = Relatedness(sibpairs)
    for i, chrom in enumerate(chroms):
        error_probs_i = None if error_probs is None else error_probs[i]
        for j, rows in enumerate(batches, 1):
            outprefix = args.out
            if args.batches > 1:
                print(f"Inferring IBD on batch {j} of {len(batches)}")
                outprefix += str(j)
            posterior = outprefix if getattr(args, 'posterior', False) else None
            with relatedness.chromosome_batch(rows, posterior) as store:
                joint_rows = joint_mask[rows]
                snipar.ibd.infer_ibd_chr(
                    sibpairs[rows], error_prob, error_probs_i, outprefix, bedfile=bedfiles[i],
                    bgenfile=bgenfiles[i], chrom=None if chrom == 0 else chrom,
                    min_length=min_length, mapfile=args.map, ibdmatrix=args.ibdmatrix,
                    ld_out=args.ld_out, min_maf=min_maf, max_missing=max_missing,
                    max_error=max_error, num_threads=num_threads, posterior_callback=store,
                    joint=joint_rows.any(), joint_pair_mask=joint_rows,
                )
    if relatedness.chromosomes.any():
        print(f'Wrote sibling relatedness to {relatedness.write(args.out)}')

if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
