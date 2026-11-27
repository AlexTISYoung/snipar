#!/usr/bin/env python
"""Infers direct effects, non-transmitted coefficients (NTCs), and population effects of genome-wide SNPs on a phenotype.

Minimally: the script requires observed genotypes on phenotyped individuals along with a phenotype file.
If no imputed parental genotypes are provided, a pedigree file is required, and the script will analyze
samples with siblings and/or both parents genotyped by default.  

**Usage**

    ``gwas.py phenotype.txt --bed chr_@ --imp chr_@ --chr_range 1 --cpus 1 --out chr_@_young``

Args:
@parser@

Results:
    ``OUT.sumstats.gz``, ``OUT.sumstats.hdf5``
        Per-chromosome SNP summary statistics in text and HDF5 formats. Use ``--no_hdf5_out`` to write only the text file.
        
"""
import argparse
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
parser.add_argument('--impute_unrel', action='store_true', default=False, help='Whether to include unrelated individuals and impute their parental genotypes linearly or not. See Unified estimator in Guan et al.')
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
parser.add_argument('--threads',type=int,help='Number of threads to use per CPU.',default=1)
parser.add_argument('--batch_size',type=int,help='Maximum SNPs per batch. For analyses using imputed parents, batches contain matched observed and imputed SNPs. Smaller batches use less memory.',default=None)
parser.add_argument('--no_hdf5_out',action='store_true',help='Suppress HDF5 output of summary statistics',default=False)
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))

def main(args):
    """Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    import os, time
    if args.cpus < 1 or args.threads < 1:
        raise ValueError('--cpus and --threads must be positive integers')
    if args.batch_size is not None and args.batch_size < 1:
        raise ValueError('--batch_size must be a positive integer')
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
    # Check if mutually incompatible design arguments supplied
    if int(args.robust) + int(args.sib_diff) + int(args.impute_unrel) > 1:
        raise argparse.ArgumentTypeError('Only supply one of --robust, --sib_diff and --impute_unrel.')
    if (args.robust or args.sib_diff or args.impute_unrel) and args.fit_sib:
        raise argparse.ArgumentTypeError('Cannot fit sib IGE with --robust, --sib_diff or --impute_unrel.')
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
        trios_sibs = not (args.sib_diff or args.fit_sib)
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
    y = read.phenotype.read_phenotype(args.phenofile, column=args.phen, column_index=args.phen_index, na_values=args.missing_char)
    if args.keep is not None:
        keep = np.loadtxt(args.keep, dtype=str, usecols=0, ndmin=1)
        y.filter_ids(y.ids[np.isin(y.ids, keep)])
    if args.imp is None:
        print('Reading pedigree from ' + args.pedigree)
        ped = np.loadtxt(args.pedigree, dtype=str, ndmin=2)
        if ped.shape[1] < 4:
            raise ValueError('Not enough columns in pedigree file')
        if ped.shape[1] > 4:
            print('WARNING: only the first four pedigree columns will be used.')
        if ped.size and tuple(ped[0, :2]) == ('FID', 'IID'):
            ped = ped[1:]
        _, ped = get_sibpairs_from_ped(ped[:, :4])
        pedigrees = [(ped, None)] * len(chroms)
    else:
        pedigrees = [read.build_ped_from_par_gts(path) for path in pargts_list]
    if args.robust:
        print('== The robust estimator will be used.')
    elif args.sib_diff:
        print('== The sib-difference estimator will be used.')
    elif args.impute_unrel:
        print('== The unified estimator will be used.')
    elif trios_sibs:
        print('== Defaulting to meta-analyzing trios and siblings using genetic differences between siblings.')
    else:
        print('== The Young estimator will be used.')

    observed_files = bedfiles if args.bed is not None else bgenfiles
    ids = y.ids
    trios_by_chromosome = []
    for observed, (ped, imp_fams) in zip(observed_files, pedigrees):
        trios = None
        if args.sib_diff:
            ids, fam_labels = read.get_ids_with_sibs(observed, ped, ids, ibdrel_path=args.grm)
        elif trios_sibs and not args.impute_unrel:
            ids, fam_labels, genotype_indices = read.get_ids_with_trios_sibs(
                observed, ped, ids, ibdrel_path=args.grm, return_info=True)
            trios = ids[np.all(genotype_indices >= 0, axis=1)]
        else:
            if trios_sibs:
                ids = remove_sibs(ped, observed, ids)
            ids, fam_labels = read.get_ids_with_par(observed, ped, imp_fams, ids,
                                                     sib=args.fit_sib, include_unrel=args.impute_unrel,
                                                     ibdrel_path=args.grm)
        trios_by_chromosome.append(trios)
    trios_sibs = trios_sibs and not args.impute_unrel
    print(f'Using {len(ids)} individuals eligible on all chromosomes.')
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
    y.filter_ids(ids)
    covariates = None
    if args.covar is not None:
        print('Reading covariates')
        covariates = read.phenotype.read_covariates(args.covar, pheno_ids=y.ids, missing_char=args.missing_char)
        covariates.filter_ids(y.ids)

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
    # Fit the null model on the same ordered cohort used for the SNP regressions.
    model = slmm.LinearMixedModel(y.gts.data.reshape(-1), varcomps=varcomps,
                                  varcomp_arr_lst=varcomp_lst,
                                  covar_X=None if covariates is None else covariates.gts.data,
                                  add_intercept=covariates is not None, drop_collinear_covariates=True)
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
    sigmas = model.varcomps
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
            print('Estimating SNP effects')
        ped, imp_fams = pedigrees[i]
        trios_only = trios_sibs and np.all(np.isin(y.ids, trios_by_chromosome[i]))
        process_chromosome(chroms[i], y, varcomp_lst,
                           ped, imp_fams, sigmas, args.out, covariates, 
                           bedfile=bedfiles[i], bgenfile=bgenfiles[i],
                           par_gts_f=pargts_list[i], fit_sib=args.fit_sib, sib_diff=args.sib_diff, parsum=args.parsum, standard_gwas=False,
                           impute_unrel=args.impute_unrel, robust=args.robust, trios_sibs=trios_sibs, trios_only=trios_only,
                           max_missing=args.max_missing, min_maf=args.min_maf,  
                           no_hdf5_out=args.no_hdf5_out, no_txt_out=False, cpus=args.cpus, add_jitter=False,
                           debug=False,batch_size=args.batch_size)
    print(f'Time used: {time.time() - start:.2f}s')

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
