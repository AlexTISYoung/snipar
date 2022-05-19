#!/usr/bin/env python
"""Output multidimensional array with proband and (imputed/observed) parental genotypes

Args:
@parser@

Results:
    gts_matrix.hdf5
        For each chromosome, a HDF5 file containing proband and (imputed/observed) parental genotypes
        
"""
import argparse, h5py
import snipar.read as read
import numpy as np
import snipar.lmm as lmm
from snipar.utilities import *
from snipar.gwas import *
from numba import set_num_threads
from numba import config as numba_config
from snipar.pedigree import get_sibpairs_from_ped
from snipar.utilities import get_parser_doc

def write_genoarray(G, outfile, parsum, sib, ped=None):
    """
    Write fitted SNP effects and other parameters to output HDF5 file.
    """
    print('Writing output to ' + outfile)
    outfile = h5py.File(outfile, 'w')
    outbim = np.column_stack((G.chrom,G.sid,G.pos,G.alleles))
    outfile['bim'] = encode_str_array(outbim)
    X_length = 1
    outcols = ['proband']
    if sib:
        X_length += 1
        outcols.append('sib')
    if parsum:
        X_length += 1
        outcols.append('parent_sum')
    else:
        X_length += 2
        outcols = outcols + ['paternal','maternal']
    # Write genotypes
    outfile.create_dataset('genotypes', G.shape, dtype='f', chunks=True,
                           compression='gzip', compression_opts=9)
    outfile['genotypes'][:] = G.gts
    # Write ids
    outfile['ids'] = encode_str_array(G.ids)
    # Write families
    outfile['fams'] = encode_str_array(G.fams)
    # Write genotypes
    outfile['genotype_cols'] = encode_str_array(outcols)
    # Write pedigree
    if ped is not None:
        outfile['pedigree'] = encode_str_array(ped)
    outfile.close()


def gts_matrix_chromosome(chrom_out, pedigree, outprefix, bedfile=None, bgenfile=None, par_gts_f=None,
                        fit_sib=False, parsum=False, max_missing=5, min_maf=0.01):
    if bedfile is None and bgenfile is None:
        raise(ValueError('Must supply either bed or bgen file with observed genotypes'))
    if bedfile is not None and bgenfile is not None:
        raise(ValueError('Both --bed and --bgen specified. Please specify one only'))
    if bedfile is not None:
        bed = Bed(bedfile,count_A1 = True)
        snp_ids = bed.sid
        pos = np.array(bed.pos[:,2],dtype=int)
        alleles = np.loadtxt(bedfile.split('.bed')[0]+'.bim',dtype=str,usecols=(4,5))
        chrom = np.array(bed.pos[:,0],dtype=int)
    elif bgenfile is not None:
        bgen = open_bgen(bgenfile, verbose=False)
        snp_ids = bgen.ids
        # If SNP IDs are broken, try rsids
        if np.unique(snp_ids).shape[0] == 1:
            snp_ids = bgen.rsids
        pos = np.array(bgen.positions,dtype=int)
        alleles = np.array([x.split(',') for x in bgen.allele_ids])
        chrom = np.array(bgen.chromosomes,dtype='U2')
        # If chromosomse unknown, set to chromosome inferred from filename
        chrom[[len(x)==0 for x in chrom]] = chrom_out
    print('Found '+str(snp_ids.shape[0])+' SNPs')
    # Remove duplicates
    unique_snps, counts = np.unique(snp_ids, return_counts=True)
    non_duplicate = set(unique_snps[counts==1])
    if np.sum(counts>1)>0:
        print('Removing '+str(np.sum(counts>1))+' duplicate SNP ids')
        not_duplicated = np.array([x in non_duplicate for x in snp_ids])
        snp_ids = snp_ids[not_duplicated]
        pos = pos[not_duplicated]
        chrom = chrom[not_duplicated]
        alleles = alleles[not_duplicated,:]
    # Get genotype matrix
    G = read.get_gts_matrix(ped=pedigree, bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, snp_ids=snp_ids, 
                            parsum=parsum, sib=fit_sib, verbose=True, print_sample_info=True)
    # Filter
    G.filter_maf(min_maf)
    G.filter_missingness(max_missing)
    # Write to file 
    hdf5_outfile = outfile_name(outprefix, '.genoarray.hdf5')
    write_genoarray(G, hdf5_outfile, parsum, fit_sib, ped=pedigree)

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
parser.add_argument('--out', type=str, help="The summary statistics will output to this path, one file for each chromosome. If the path contains '@', the '@' will be replaced with the chromosome number. Otherwise, the summary statistics will be output to the given path with file names chr_1.sumstats.gz, chr_2.sumstats.gz, etc. for the text sumstats, and chr_1.sumstats.hdf5, etc. for the HDF5 sumstats",default='./')
parser.add_argument('--pedigree',type=str,help='Address of pedigree file. Must be provided if not providing imputed parental genotypes.',default=None)
parser.add_argument('--parsum',action='store_true',help='Regress onto proband and sum of (imputed/observed) maternal and paternal genotypes. Default uses separate paternal and maternal genotypes when available.',default = False)
parser.add_argument('--fit_sib',action='store_true',help='Fit indirect effect from sibling ',default=False)
parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)', default=0.01)
parser.add_argument('--threads',type=int,help='Number of threads to use for IBD inference. Uses all available by default.',default=None)
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
                        fit_sib=args.fit_sib, parsum=args.parsum, max_missing=5, min_maf=0.01)
        
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)