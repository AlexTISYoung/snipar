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
import argparse, os, time
from snipar.numrange import parseNumRange, NumRangeAction
from snipar.docgen import get_parser_doc
######### Command line arguments #########
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
parser.add_argument('--bgen', type=str,help='Address of the phased genotypes in .bgen format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--bed', type=str,help='Address of the unphased genotypes in .bed format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--imp', type=str, help='Address of hdf5 files with imputed parental genotypes (without .hdf5 suffix). If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range (chr_range is an optional parameters for this script).', default = None)
parser.add_argument('--pedigree',type=str,help='Address of pedigree file. Must be provided if not providing imputed parental genotypes.',default=None)
parser.add_argument('--covar',type=str,help='Path to file with covariates: plain text file with columns FID, IID, covar1, covar2, ..', default=None)
parser.add_argument('--chr_range', type=parseNumRange, nargs='*', action=NumRangeAction, help='Chromosomes to analyse. Should be a series of ranges with x-y format (e.g. 1-22) or integers.', default=None)
parser.add_argument('--out', type=str, help="The summary statistics will output to this path, one file for each chromosome. If the path contains '@', the '@' will be replaced with the chromosome number. Otherwise, the summary statistics will be output to the given path with file names chr_1.sumstats.gz, chr_2.sumstats.gz, etc. for the text sumstats, and chr_1.sumstats.hdf5, etc. for the HDF5 sumstats",default='./')
parser.add_argument('--grm', type=str, help='Path to GRM file giving pairwise relatednsss information. Designed to work with KING IBD segment inference output (.seg file).', default=None)
parser.add_argument('--grmgz', type=str, help='Path to GRM in GCTA grm.gz format (without .grm.gz suffix). Assumes .grm.id file with same root path also available.', default=None)
parser.add_argument('--sparse_thresh', type=float, help='Threshold of GRM sparsity — elements below this value are set to zero', default=0.05)
parser.add_argument('--impute_unrel', action='store_true', default=False, help='Whether to include unrelated individuals and impute their parental genotypes lineary or not. See Unified estimator in Guan et al.')
parser.add_argument('--robust', action='store_true', default=False, help='Use the robust estimator')
parser.add_argument('--sib_diff', action='store_true', default=False, help='Use the sibling difference method')
parser.add_argument('--parsum',action='store_true',help='Regress onto proband and sum of (imputed/observed) maternal and paternal genotypes. Default uses separate paternal and maternal genotypes when available.',default = False)
parser.add_argument('--fit_sib',action='store_true',help='Fit indirect effect from sibling ',default=False)
parser.add_argument('--phen',type=str,help='Name of the phenotype to be analysed — case sensitive. Default uses first phenotype in file.', default=None)
parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)', default=1)
parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)', default='NA')
parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)', default=0.01)
parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)', default=5)
parser.add_argument('--vc_out', type=str, help='Prefix of output filename for variance component array (without .npy).')
parser.add_argument('--vc_list', type=float, nargs='+', default=None, help='Pass in variance components as a list of floats.')
parser.add_argument('--no_sib_var', action='store_true', default=False, help='Do not fit sibling variance component. Not recommended for family-GWAS.')
parser.add_argument('--keep', default=None, type=str, help='Filename of IDs to be kept for analysis (No header).')
parser.add_argument('--cpus', type=int, help='Number of cpus to distribute batches across', default=1)
parser.add_argument('--threads',type=int,help='Number of threads to use per CPU. Uses all available by default.',default=1)
parser.add_argument('--no_hdf5_out',action='store_true',help='Suppress HDF5 output of summary statistics',default=False)
parser.add_argument('--batch_size',type=int,help='Batch size of SNPs to load at a time (reduce to reduce memory requirements)',default=100000)
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))
args = parser.parse_args()
num_threads = args.threads
print('Number of threads for numpy: '+str(num_threads))
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
import numpy as np
from numba import set_num_threads
from numba import config as numba_config
import snipar.read as read
import snipar.slmm as slmm
from snipar.gwas import process_chromosome
from snipar.pedigree import get_sibpairs_from_ped
from snipar.preprocess import remove_sibs
from snipar.utilities import parse_obsfiles, parse_filelist, make_id_dict
def main(args):
    """"Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    # Check if mutually incompatible design arguments supplied
    if int(args.robust) + int(args.sib_diff) + int(args.impute_unrel) > 1:
        raise argparse.ArgumentTypeError('Only supply one of --robust, --sib_diff and --impute_unrel.')
    if (args.robust or args.sib_diff) and args.fit_sib:
        raise argparse.ArgumentTypeError('Cannot fit sib IGE while --robust or --sib_diff is supplied.')
    if (args.sib_diff or args.robust) and args.parsum:
        raise argparse.ArgumentTypeError('--parsum is ignored with --sib_diff and --robust estimator.')
    # Check if GRM provided
    if args.grm is None and args.grmgz is None:
        print('No GRM provided.')
        args.no_grm_var = True
    else:
        args.no_grm_var = False
        if args.grm:
            if os.path.exists(args.grm):
                print('Using pairwise relationships from: '+args.grm+' for GRM variance component')
            else:
                raise argparse.ArgumentTypeError('GRM not found.')
        elif args.grmgz:
            if os.path.exists(args.grmgz+'.grm.gz'):
                print('Using: '+args.grmgz+' for GRM variance component')
            else:
                raise argparse.ArgumentTypeError('GRM file not found.')
    # Check if fitting sibling variance component
    if args.no_sib_var:
        if args.no_grm_var:
            raise argparse.ArgumentTypeError('Need to fit at least one of sibling and GRM variance components.')
        else:
            print('Fitting GRM and residual variance components. Note this may result in statistically inefficient direct genetic effect estimates. It is advised to fit the sibling variance component.')
    elif args.no_grm_var:
        print('Fitting sibling and residual variance components')
    else:
        print('Fitting sibling, GRM, and residual variance components')
    # Set number of threads for numba
    if args.threads is not None:
        num_threads = min([args.threads, numba_config.NUMBA_NUM_THREADS])
    else:
        num_threads = numba_config.NUMBA_NUM_THREADS
    set_num_threads(num_threads)
    print('Number of threads for numba: '+str(num_threads))
    print('Number of processes: '+str(args.cpus))
    # Print where output written to
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
            raise ValueError("The robust estimator requires imputed parental genotypes. If these are not available, remove the --robust flag and a meta analysis of trios and siblings will be performed.")
        trios_sibs = True if not args.sib_diff else False
        trios_only = False
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
        trios_only = False
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
    y = read.phenotype.read_phenotype(args.phenofile, column=args.phen, column_index=args.phen_index, na_values=args.missing_char)
    ######## Read covariates ########
    if args.covar is not None:
        print('Reading covariates')
        covariates = read.phenotype.read_covariates(args.covar, pheno_ids=y.ids, missing_char=args.missing_char)
        # Match to pheno ids
        # covariates.filter_ids(y.ids)
    else:
        covariates = None
    ########## Read pedigree ########
    if args.imp is None:
        # Read pedigree from file
        print('Reading pedigree from '+str(args.pedigree))
        ped = np.loadtxt(args.pedigree,dtype=str)[1:,]
        if ped.shape[1] < 4:
            raise(ValueError('Not enough columns in pedigree file'))
        elif ped.shape[1] > 4:
            print('WARNING: pedigree file has more than 4 columns. The first four columns only will be used.')
        # Remove rows with missing parents
        sibpairs, ped = get_sibpairs_from_ped(ped)
        imp_fams = None
    else:
        ped, imp_fams = read.build_ped_from_par_gts(pargts_list[0])
        if args.sib_diff:
            imp_fams = None # imputation not used
    # # Check if sibling pairs are present
    # if sibpairs is not None:
    #     print('Found '+str(sibpairs.shape[0])+' sibling pairs in pedigree')
    # else:
    #     print('Found 0 sibling pairs')
    #     if args.sib_diff:
    #         raise argparse.ArgumentTypeError('No sibling pairs found in pedigree. Cannot use sibling difference estimator.')
    #     if not args.no_sib_var:
    #         print('No sibling pairs found in pedigree. Dropping sibling variance component')
    #         args.no_sib_var = True
    ####### Print which estimator will be used #########
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
    ## Find which individuals can be used for the selected estimator ##
    if args.sib_diff:
        # Find individuals with genotyped siblings
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
                                                 y.ids, return_info=False, ibdrel_path=args.grm)
    elif trios_sibs:
        if args.impute_unrel:
            # Remove siblings
            ids = remove_sibs(
                ped, bedfiles[0] if args.bed is not None else bgenfiles[0], y.ids
            )
            # Find trios
            ids, fam_labels = read.get_ids_with_par(
                bedfiles[0] if args.bed is not None else bgenfiles[0], ped, imp_fams,
                ids, sib=args.fit_sib, include_unrel=args.impute_unrel, ibdrel_path=args.grm,
                return_info=False
            )
            trios_sibs = False
        else:
            # Find individuals with genotyped siblings and/or both parents genotyped
            ids, fam_labels, trios_only = read.get_ids_with_trios_sibs(
                bedfiles[0] if args.bed is not None else bgenfiles[0], ped, y.ids, ibdrel_path=args.grm
            )
    else:
        # Find individuals with observed and/or imputed parental genotypes
        ids, fam_labels = read.get_ids_with_par(
            bedfiles[0] if args.bed is not None else bgenfiles[0], ped, imp_fams,
            y.ids, sib=args.fit_sib, include_unrel=args.impute_unrel, ibdrel_path=args.grm,
            return_info=False
        )
    #### Filter phenotype to valid IDs ###
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
    # Read GRM if provided 
    if args.grm is not None:
        id_dict = make_id_dict(ids)
        grm_data, grm_row_ind, grm_col_ind = slmm.build_ibdrel_arr(
            args.grm, id_dict=id_dict, keep=ids, thres=args.sparse_thresh)
    elif args.grmgz is not None:
        ids, fam_labels = slmm.match_grm_ids(
            ids, fam_labels, grm_path=args.grmgz, grm_source='gcta')
        id_dict = make_id_dict(ids)
        grm_data, grm_row_ind, grm_col_ind = slmm.build_grm_arr(
            args.grmgz, id_dict=id_dict, thres=args.sparse_thresh)
    # Build sparse silbship array if fitting sib variance component
    if not args.no_sib_var:
        sib_data, sib_row_ind, sib_col_ind = slmm.build_sib_arr(fam_labels)
    # Set variance components
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
    if args.vc_list:
        varcomps = tuple(args.vc_list)
        if len(varcomps) != len(varcomp_lst) + 1:
            raise ValueError(f'Supplied varcomps length {len(varcomps)} does not match the number of variance components.')
    else:
        varcomps = None
    # Define sparse linear mixed model for variance component estimation
    if args.covar is None:
        y.gts -= y.gts.mean()
        model = slmm.LinearMixedModel(y.gts.reshape(-1).data, varcomps=varcomps, varcomp_arr_lst=varcomp_lst, covar_X=None, add_intercept=False, add_jitter=False)
    else:
        model = slmm.LinearMixedModel(y.gts.reshape(-1).data, varcomps=varcomps, varcomp_arr_lst=varcomp_lst, covar_X=covariates.gts.data, add_intercept=True, add_jitter=False)
    if not varcomps:
        print(f'Optimizing variance components...')
        # print('Nonzero entries: ', model.V.nnz, 'Number of individuals: ', model.n, 'density: ', model.V.nnz / model.n / model.n)
        start = time.time()
        model.scipy_optimize()
        print(f'Time used for variance component estimation: {time.time() - start:.2f}s')
    else:
        print('Variance components supplied.')
    if args.vc_out:
        np.savetxt(f'{args.vc_out}', np.array(model.varcomps))
        print(f'Variance components saved to {args.vc_out}.')
    ## Print variance components
    varcomps = list(i for i in model.varcomps)
    vcomp_index=0
    if not args.no_grm_var:
        print(f'GRM variance component: {varcomps[vcomp_index]:.3f}')
        print(f'Variance explained by GRM: {varcomps[vcomp_index] / np.sum(varcomps) * 100:.1f}%')
        vcomp_index+=1
    if not args.no_sib_var:
        print(f'Sibling variance component: {varcomps[vcomp_index]:.3f}')
        print(f'Variance explained by sibling component: {varcomps[vcomp_index] / np.sum(varcomps) * 100:.1f}%')
    print(f'Total variance: {np.sum(varcomps):.3f}')
    # Print implied sibling correlation
    if args.no_grm_var:
        print(f'Implied sibling correlation: {varcomps[0] / np.sum(varcomps):.3f}')
    elif not args.no_sib_var:
        print(f'Implied sibling correlation: {(0.5*varcomps[0] + varcomps[1]) / np.sum(varcomps):.3f}')
    else:
        print(f'Implied sibling correlation: {0.5*varcomps[0] / np.sum(varcomps):.3f}')
    # Define variance components: does this need to be done twice?     
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
    ######### Process chromosomes ###########
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
                           impute_unrel=args.impute_unrel, robust=args.robust, trios_sibs=trios_sibs, trios_only=trios_only,
                           max_missing=args.max_missing, min_maf=args.min_maf, batch_size=args.batch_size, 
                           no_hdf5_out=args.no_hdf5_out, no_txt_out=False, cpus=args.cpus, add_jitter=False,
                           debug=False)
    print(f'Time used: {time.time() - start:.2f}s')

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)