#!/usr/bin/env python
"""Output multidimensional array with proband and (imputed/observed) parental genotypes

Args:
@parser@

Results:
    gts_matrix.hdf5
        For each chromosome, a HDF5 file containing proband and (imputed/observed) parental genotypes
        
"""
import argparse

import h5py
import snipar.read as read
import numpy as np
from snipar.utilities import encode_str_array, outfile_name, parse_filelist, parse_obsfiles
from snipar.docgen import get_parser_doc
from numba import set_num_threads
from numba import config as numba_config
from snipar.pedigree import get_sibpairs_from_ped
from snipar.numrange import parseNumRange, NumRangeAction

def write_genoarray(G, outfile, parsum, sib, ped=None):
    """Write family genotypes and matching sample and variant metadata."""
    print('Writing output to ' + outfile)
    outcols = ['proband']
    if sib:
        outcols.append('sib')
    if parsum:
        outcols.append('parent_sum')
    else:
        outcols.extend(['paternal', 'maternal'])
    if len(outcols) != G.shape[1]:
        raise ValueError('Genotype columns do not match the genotype matrix')
    with h5py.File(outfile, 'w') as output:
        output['bim'] = encode_str_array(np.column_stack((G.chrom, G.sid, G.pos, G.alleles)))
        genotypes = G.gts.astype('f', copy=False).filled(np.nan) if G.has_NAs else G.gts.data
        output.create_dataset('genotypes', data=genotypes, dtype='f',
                              compression='gzip', compression_opts=3, shuffle=True)
        output['ids'] = encode_str_array(G.ids)
        output['fams'] = encode_str_array(G.fams)
        output['genotype_cols'] = encode_str_array(outcols)
        if ped is not None:
            output['pedigree'] = encode_str_array(ped)


def gts_matrix_chromosome(chrom_out, pedigree, outprefix, bedfile=None, bgenfile=None, par_gts_f=None,
                        fit_sib=False, parsum=False, max_missing=5, min_maf=0.01, imp_fams=None):
    if imp_fams is None and par_gts_f is not None:
        imputed_pedigree, imp_fams = read.build_ped_from_par_gts(par_gts_f)
        if pedigree is None:
            pedigree = imputed_pedigree
    G = read.get_gts_matrix(ped=pedigree, imp_fams=imp_fams, bedfile=bedfile, bgenfile=bgenfile,
                            par_gts_f=par_gts_f, parsum=parsum, sib=fit_sib,
                            verbose=True, print_sample_info=True)
    if chrom_out:
        G.chrom = np.where(G.chrom == '', str(chrom_out), G.chrom)
    G.filter_maf(min_maf)
    G.filter_missingness(max_missing)
    hdf5_outfile = outfile_name(outprefix, '.genoarray.hdf5', chrom=chrom_out or None)
    # The genotype reader selects parental sums automatically when every parent is imputed.
    write_genoarray(G, hdf5_outfile, G.shape[1] == 2 + int(fit_sib), fit_sib, ped=pedigree)

######### Command line arguments #########
parser=argparse.ArgumentParser()
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
parser.add_argument('--out', type=str, help="Output prefix. Writes one .genoarray.hdf5 file per chromosome, replacing @ with the chromosome number if present.",default='./')
parser.add_argument('--pedigree',type=str,help='Address of pedigree file. Must be provided if not providing imputed parental genotypes.',default=None)
parser.add_argument('--parsum',action='store_true',help='Regress onto proband and sum of (imputed/observed) maternal and paternal genotypes. Default uses separate paternal and maternal genotypes when available.',default = False)
parser.add_argument('--fit_sib',action='store_true',help='Fit indirect effect from sibling ',default=False)
parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)', default=0.01)
parser.add_argument('--threads',type=int,help='Number of threads to use. Uses all available by default.',default=None)
parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)', default=5)
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))
# Set number of threads
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
    ped = None
    if args.imp is None:
        print('Reading pedigree from ' + args.pedigree)
        ped = np.loadtxt(args.pedigree, dtype=str, ndmin=2)
        if ped.shape[1] < 4:
            raise ValueError('Not enough columns in pedigree file')
        if ped.shape[1] > 4:
            print('Warning: only the first four pedigree columns will be used')
        if ped.size and tuple(ped[0, :2]) == ('FID', 'IID'):
            ped = ped[1:]
        _, ped = get_sibpairs_from_ped(ped[:, :4])

    for i in range(chroms.shape[0]):
        if args.bed is not None:
            print('Observed genotypes file: '+bedfiles[i])
        if args.bgen is not None:
            print('Observed genotypes file: '+bgenfiles[i])
        if args.imp is not None:
            print('Imputed genotypes file: '+pargts_list[i])
        if chroms.shape[0]>1:
            print('Obtaining genotype matrix for '+str(chroms[i]))
        else:
            print('Obtaining genotype matrix')
        gts_matrix_chromosome(chroms[i], ped, args.out, bedfile=bedfiles[i], bgenfile=bgenfiles[i], par_gts_f=pargts_list[i],
                        fit_sib=args.fit_sib, parsum=args.parsum, max_missing=args.max_missing, min_maf=args.min_maf)
        
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
