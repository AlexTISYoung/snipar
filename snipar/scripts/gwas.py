#!/usr/bin/env python
"""Infers direct effects, non-transmitted coefficients (NTCs), and population effects of genome-wide SNPs on a phenotype.

Minimally: the script requires observed genotypes on phenotyped individuals and their parents, and/or 
parental genotypes imputed by snipar's impute.py script, along with a phenotype file. 

Args:
@parser@

Results:
    sumstats.gz
        For each chromosome, a gzipped text file containing the SNP level summary statistics. 
        
"""
import argparse
import os
import time

# import logging

######### Command line arguments #########
parser=argparse.ArgumentParser()
parser.add_argument('--threads',type=int,help='Number of threads to use for IBD inference. Uses all available by default.',default=None)
args, extra_args = parser.parse_known_args()
num_threads = args.threads if args.threads is not None else 1
# export OMP_NUM_THREADS=...
os.environ['OMP_NUM_THREADS'] = str(num_threads)
# export OPENBLAS_NUM_THREADS=...
os.environ['OPENBLAS_NUM_THREADS'] = str(num_threads)
# export MKL_NUM_THREADS=...
os.environ['MKL_NUM_THREADS'] = str(num_threads)
# export VECLIB_MAXIMUM_THREADS=...
os.environ['VECLIB_MAXIMUM_THREADS'] = str(num_threads)
# export NUMEXPR_NUM_THREADS=...
os.environ['NUMEXPR_NUM_THREADS'] = str(num_threads)
print('Number of threads for numpy: '+str(num_threads))

import numpy as np
from numba import set_num_threads
from numba import config as numba_config
import snipar.read as read
import snipar.slmm as slmm
from snipar.gwas import process_chromosome
from snipar.pedigree import get_sibpairs_from_ped
from snipar.preprocess import remove_sibs
from snipar.utilities import get_parser_doc, parseNumRange, NumRangeAction, parse_obsfiles, parse_filelist, make_id_dict

parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
parser.add_argument('--bgen', type=str,help='Address of the phased genotypes in .bgen format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--bed', type=str,help='Address of the unphased genotypes in .bed format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--imp', type=str, help='Address of hdf5 files with imputed parental genotypes (without .hdf5 suffix). If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range (chr_range is an optional parameters for this script).', default = None)
parser.add_argument('--grm_path', type=str, help='Path to gcta grm gz file (without .gz prefix).', default=None)
parser.add_argument('--ibdrel_path', type=str, help='Path to KING IBD segment inference output (without .seg prefix).', default=None)
parser.add_argument('--sparse_thres', type=float, help='Threshold of GRM/IBD sparsity', default=0.05)
parser.add_argument('--chr_range', type=parseNumRange, nargs='*', action=NumRangeAction, help='number of the chromosomes to be imputed. Should be a series of ranges with x-y format or integers.', default=None)
parser.add_argument('--out', type=str, help="The summary statistics will output to this path, one file for each chromosome. If the path contains '@', the '@' will be replaced with the chromosome number. Otherwise, the summary statistics will be output to the given path with file names chr_1.sumstats.gz, chr_2.sumstats.gz, etc. for the text sumstats, and chr_1.sumstats.hdf5, etc. for the HDF5 sumstats",default='./')
parser.add_argument('--impute_unrel', action='store_true', default=False, help='Whether to include unrelated individuals and impute their parental genotypes or not.')
parser.add_argument('--robust', action='store_true', default=False, help='Whether to use the robust estimator')
parser.add_argument('--sib_diff', action='store_true', default=False, help='Use sibling difference method')
# parser.add_argument('--standard_gwas', action='store_true', default=False, help='Whether to fit standard gwas, i.e., without modelling parental genotypes.')
parser.add_argument('--pedigree',type=str,help='Address of pedigree file. Must be provided if not providing imputed parental genotypes.',default=None)
parser.add_argument('--parsum',action='store_true',help='Regress onto proband and sum of (imputed/observed) maternal and paternal genotypes. Default uses separate paternal and maternal genotypes when available.',default = False)
parser.add_argument('--fit_sib',action='store_true',help='Fit indirect effect from sibling ',default=False)
parser.add_argument('--covar',type=str,help='Path to file with covariates: plain text file with columns FID, IID, covar1, covar2, ..', default=None)
parser.add_argument('--fit_res', action='store_true', default=False,help='Use residualized phenotypes.')
parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)', default=1)
parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)', default=0.01)
parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)', default=5)
parser.add_argument('--batch_size',type=int,help='Batch size of SNPs to load at a time (reduce to reduce memory requirements)',default=100000)
parser.add_argument('--no_hdf5_out',action='store_true',help='Suppress HDF5 output of summary statistics',default=False)
parser.add_argument('--no_txt_out',action='store_true',help='Suppress text output of summary statistics',default=False)
parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)', default='NA')
parser.add_argument('--no_grm_var', action='store_true', default=False, help='whether to exclude grm variance component.')
parser.add_argument('--no_sib_var', action='store_true', default=False, help='whether to exclude sib variance component.')
parser.add_argument('--vc_out', type=str, help='Prefix of output filename for variance component array (without .npy).')
parser.add_argument('--vc_in', type=str, help='Prefix of input filename for variance component array (without .npy).')
parser.add_argument('--vc_list', type=float, nargs='+', default=None, help='Pass in variance components as a list of floats.')
parser.add_argument('--vc_only', action='store_true', help='Only perform variance component estimation.')
parser.add_argument('--keep', default=None, type=str, help='Filename of IDs to be kept for analysi (No header).')
parser.add_argument('--cpus', type=int, help='Number of cpus to distribute batches across', default=1)


parser.add_argument('--debug', action='store_true', default=False, help='Debug code in single process mode.')
# parser.add_argument('--loglevel', type=str, default='INFO', help='Case insensitive Logging level: INFO, DEBUG, ...')


__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))


def main(args):
    """"Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    # loglevel = args.loglevel
    # FORMAT = '%(process)d :: %(asctime)-10s :: %(levelname)s :: %(message)s'
    # numeric_level = getattr(logging, loglevel.upper(), None)
    # if not isinstance(numeric_level, int):
    #     raise ValueError('Invalid log level: %s' % loglevel)
    # logging.basicConfig(
    #     format=FORMAT, level=numeric_level)
    # logger = logging.getLogger(__name__)

    if int(args.robust) + (args.sib_diff) + (args.impute_unrel) > 1:
        raise argparse.ArgumentTypeError('Only one of --robust, --sib_diff and --impute_unrel.')
    if args.no_sib_var and not args.no_grm_var:
        raise argparse.ArgumentTypeError('Only one of --no_sib_var and --no_grm_var.')
    if not args.no_grm_var:
        if args.ibdrel_path is None and args.grm_path is None:
            raise argparse.ArgumentTypeError('Need to input GRM.')
    # if args.standard_gwas and (args.robust or args.sib_diff):
    #     raise argparse.ArgumentTypeError('Cannot set --robust or --sib_diff if --standard_gwas is set.')
    if (args.robust or args.sib_diff) and args.fit_sib:
        raise argparse.ArgumentTypeError('Cannot fit sib IGE while --robust or --sib_diff is supplied.')
    if args.sib_diff and args.parsum:
        raise argparse.ArgumentTypeError('--parsum is meaningless while --sib_diff is supplied.')
            
    # Set number of threads for numba
    if args.threads is not None:
        num_threads = min([args.threads, numba_config.NUMBA_NUM_THREADS])
    else:
        num_threads = numba_config.NUMBA_NUM_THREADS
    set_num_threads(num_threads)
    print('Number of threads for numba: '+str(num_threads))
    print('Number of processes: '+str(args.cpus))

    # logger.info(f'Output will be written to {args.out}.')
    if args.out == './':
        print(f'Output will be written to chr_@')
    elif '@' in args.out:
        print(f'Output will be written to {args.out}')
    else:
        print(f'Output will be written to {args.out}_chr_@')

    # Check arguments
    if args.bed is None and args.bgen is None:
        raise(ValueError('Must provide one of --bedfiles and --bgenfiles'))
    if args.bed is not None and args.bgen is not None:
        raise(ValueError('Both bed files and bgen files provided. Please provide only one'))
    if args.imp is None and args.pedigree is None:
        raise(ValueError('Must provide pedigree if not providing imputed parental genotypes file(s)'))

    # Find observed and imputed files
    if args.imp is None:
        print("No imputation provided.")
        if args.robust:
            raise ValueError("The robust estimator requires imputed parental genotypes. If these are not available, remove the --robust flag and a meta analysis of trios and siblings, using genetic differences between siblings will be performed.")
        trios_sibs = True if not args.sib_diff else False
        # if trios_sibs and not args.impute_unrel:
        #     print('Defaulting to meta-analyzing trios and siblings using genetic differences between siblings.')
        if args.bed is not None:
            bedfiles, chroms = parse_obsfiles(args.bed, 'bed', chromosomes=args.chr_range)
            bgenfiles = [None for x in range(chroms.shape[0])]
        elif args.bgen is not None:
            bgenfiles, chroms = parse_obsfiles(args.bgen, 'bgen', chromosomes=args.chr_range)
            bedfiles = [None for x in range(chroms.shape[0])]
        pargts_list = [None for x in range(chroms.shape[0])]
    else:
        trios_sibs = False
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
        # covariates.filter_ids(y.ids)
    else:
        covariates = None

    # if --robust is true, set impute_unrel to false
    if args.robust:
        args.impute_unrel = False

    # Read pedigree
    if args.imp is None:
        # logger.info('Reading pedigree from '+str(args.pedigree))
        print('Reading pedigree from '+str(args.pedigree))
        ped = np.loadtxt(args.pedigree,dtype=str)[1:,]
        if ped.shape[1] < 4:
            raise(ValueError('Not enough columns in pedigree file'))
        elif ped.shape[1] > 4:
            # logger.warning('Warning: pedigree file has more than 4 columns. The first four columns only will be used')
            print('WARNING: pedigree file has more than 4 columns. The first four columns only will be used')
        # Remove rows with missing parents
        sibpairs, ped = get_sibpairs_from_ped(ped)
        imp_fams = None
        if sibpairs is not None:
            # logger.info('Found '+str(sibpairs.shape[0])+' sibling pairs in pedigree')
            print('Found '+str(sibpairs.shape[0])+' sibling pairs in pedigree')
        else:
            # logger.info('Found 0 sibling pairs')
            print('Found 0 sibling pairs')
    else:
        # Read pedigree
        ped, imp_fams = read.build_ped_from_par_gts(pargts_list[0])
        if args.sib_diff:
            imp_fams = None # imputation not used

    if args.robust:
        print("== The robust estimator will be used.")
    elif args.sib_diff:
        print("== The sib-difference estimator will be used.")
    elif trios_sibs:
        if args.impute_unrel:
            print("== The unified estimator will be used.")
        else:
            print("== Defaulting to meta-analyzing trios and siblings using genetic differences between siblings.")
    elif args.impute_unrel:
        print("== The unified estimator will be used.")
    else:
        print("== The Young estimator will be used.")
    print("Check if the dataset is feasible for the selected estimator...")
    if args.sib_diff:
        # unrelated_inds = None
        ids = y.ids
        fam_labels = y.fams
        ped_dict = make_id_dict(ped,1)
        y.filter_ids(ped[:,1])
        print(str(y.shape[0])+' individuals with phenotype values found in pedigree')
        ped_indices = np.array([ped_dict[x] for x in y.ids])
        y.fams = ped[ped_indices,0]
        ids = y.ids
        fam_labels = y.fams
        ids, fam_labels = read.get_ids_with_sibs(bedfiles[0] if args.bed is not None else bgenfiles[0], ped, 
                                                 y.ids, return_info=False, ibdrel_path=args.ibdrel_path)
    elif trios_sibs:
        if args.impute_unrel:
            ids = remove_sibs(
                ped, bedfiles[0] if args.bed is not None else bgenfiles[0], y.ids
            )

            ids, fam_labels = read.get_ids_with_par(
                bedfiles[0] if args.bed is not None else bgenfiles[0], ped, imp_fams,
                ids, sib=args.fit_sib, include_unrel=args.impute_unrel, ibdrel_path=args.ibdrel_path,
                return_info=False
            )
            print(ids.shape)
            
            trios_sibs = False
        else:
            ids, fam_labels = read.get_ids_with_trios_sibs(
                bedfiles[0] if args.bed is not None else bgenfiles[0], ped, y.ids, ibdrel_path=args.ibdrel_path
            )
    else:
        # unrelated_inds = None
        ids, fam_labels = read.get_ids_with_par(
            bedfiles[0] if args.bed is not None else bgenfiles[0], ped, imp_fams,
            y.ids, sib=args.fit_sib, include_unrel=args.impute_unrel, ibdrel_path=args.ibdrel_path,
            return_info=False
        )
    y.filter_ids(ids)
    if args.sib_diff:
        ids = y.ids
        fam_labels = y.fams
    np.testing.assert_array_equal(ids, y.ids)
    if args.covar:
        if not args.sib_diff:
            covariates.filter_ids(ids)
            np.testing.assert_array_equal(ids, covariates.ids)
        else:
            covariates.filter_ids(y.ids)


    if args.ibdrel_path is not None:
        id_dict = make_id_dict(ids)
        grm_data, grm_row_ind, grm_col_ind = slmm.build_ibdrel_arr(
            args.ibdrel_path, id_dict=id_dict, keep=ids, thres=args.sparse_thres)
    elif args.grm_path is not None:
        ids, fam_labels = slmm.match_grm_ids(
            ids, fam_labels, grm_path=args.grm_path, grm_source='gcta')
        id_dict = make_id_dict(ids)
        grm_data, grm_row_ind, grm_col_ind = slmm.build_grm_arr(
            args.grm_path, id_dict=id_dict, thres=args.sparse_thres)

    if not args.no_sib_var:
        sib_data, sib_row_ind, sib_col_ind = slmm.build_sib_arr(fam_labels)
    
    if not args.no_grm_var and not args.no_sib_var:
        varcomp_lst = (
            (grm_data, grm_row_ind, grm_col_ind),
            (sib_data, sib_row_ind, sib_col_ind),
        )
    elif not args.no_grm_var:
        varcomp_lst = (
            (grm_data, grm_row_ind, grm_col_ind),
        )
    elif not args.no_sib_var:
        varcomp_lst = (
            (sib_data, sib_row_ind, sib_col_ind),
        )
    if args.vc_list and args.vc_in:
        raise ValueError('Only one of --vc_list and --vc_in is needed.')
    if args.vc_list:
        varcomps = tuple(args.vc_list)
        if len(varcomps) != len(varcomp_lst) + 1:
            raise ValueError(f'Supplied varcomps length {len(varcomps)} does not match the number of variance components.')
    elif args.vc_in:
        varcomps = tuple(np.loadtxt(args.vc_in))
        if len(varcomps) != len(varcomp_lst) + 1:
            raise ValueError(f'Supplied varcomps length {len(varcomps)} does not match the number of variance components.')
    else:
        varcomps = None
    if args.covar is None:
        y.gts -= y.gts.mean()
        model = slmm.LinearMixedModel(y.gts.reshape(-1).data, varcomps=varcomps, varcomp_arr_lst=varcomp_lst, covar_X=None, add_intercept=False, add_jitter=False)
    else:
        if args.fit_res:
            covar_1 = np.hstack((np.ones((y.shape[0], 1), dtype=y.dtype), covariates.gts.data))
            alpha_covar = np.linalg.solve(covar_1.T.dot(covar_1), covar_1.T.dot(y.gts))
            y.gts = y.gts -  alpha_covar[0] - covariates.gts.dot(alpha_covar[1:])
            # logger.info(f'--fit_res specified. Phenotypes residualized. Variance of y: {np.var(y.gts)}')
            print(f'--fit_res specified. Phenotypes residualized. Variance of y: {np.var(y.gts)}')
            covariates = None
            model = slmm.LinearMixedModel(y.gts.reshape(-1).data, varcomps=varcomps, varcomp_arr_lst=varcomp_lst, covar_X=None, add_intercept=False, add_jitter=False)
        else:
            model = slmm.LinearMixedModel(y.gts.reshape(-1).data, varcomps=varcomps, varcomp_arr_lst=varcomp_lst, covar_X=covariates.gts.data, add_intercept=True, add_jitter=False)
    if not varcomps:
        # logger.info(f'Optimizing variance components...')
        print(f'Optimizing variance components...')
        # print('Nonzero entries: ', model.V.nnz, 'Number of individuals: ', model.n, 'density: ', model.V.nnz / model.n / model.n)
        start = time.time()
        model.scipy_optimize()
        # logger.info(f'Time for variance component estimation: {time.time() - start}s.')
        print(f'Time for variance component estimation: {time.time() - start:.2f}s.')
    else:
        # logger.info('varcomps supplied.')
        print('varcomps supplied.')
    if args.vc_out:
        np.savetxt(f'{args.vc_out}', np.array(model.varcomps))
        # logger.info(f'varcomps saved to {args.vc_out}.')
        print(f'varcomps saved to {args.vc_out}.')
    # logger.info(f'Variance components: {list(i / y.gts.data.var() for i in model.varcomps)}')
    print(f'Variance components: {list(i / y.gts.data.var() for i in model.varcomps)}')
    if args.vc_only:
        exit('Variance component estimation finished.')
    sigmas = model.varcomps
    if not args.no_grm_var:
        if not args.no_sib_var:
            varcomp_lst = (
                (grm_data, grm_row_ind, grm_col_ind),
                (sib_data, sib_row_ind, sib_col_ind),
            )
        else:
            varcomp_lst = (
                (grm_data, grm_row_ind, grm_col_ind),
            )
    elif not args.no_sib_var:
        varcomp_lst = (
            (sib_data, sib_row_ind, sib_col_ind),
        )


    start = time.time()
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
        process_chromosome(chroms[i], y, varcomp_lst,
                           ped, imp_fams, sigmas, args.out, covariates, 
                           bedfile=bedfiles[i], bgenfile=bgenfiles[i],
                           par_gts_f=pargts_list[i], fit_sib=args.fit_sib, sib_diff=args.sib_diff, parsum=args.parsum, standard_gwas=False,
                           impute_unrel=args.impute_unrel, robust=args.robust, trios_sibs=trios_sibs,
                           max_missing=args.max_missing, min_maf=args.min_maf, batch_size=args.batch_size, 
                           no_hdf5_out=args.no_hdf5_out, no_txt_out=args.no_txt_out, cpus=args.cpus, add_jitter=False,
                           debug=args.debug)
    # logger.info(f'Time used: {time.time() - start}.')
    print(f'Time used: {time.time() - start:.2f}s.')

if __name__ == "__main__":
    args = parser.parse_args(extra_args)
    main(args)