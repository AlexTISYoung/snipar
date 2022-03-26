from sibreg.sibreg import *
from sibreg.model import *
import os
from os import path
import argparse
import h5py
from bgen_reader import open_bgen
from pysnptools.snpreader import Bed
from functools import partial
import ctypes

# TODO: set logging level for imported packages
# see https://stackoverflow.com/questions/47265898/how-to-set-the-log-level-for-an-imported-module
# FORMAT = '%(asctime)-15s :: %(levelname)s :: %(name)s :: %(funcName)s :: %(message)s'
# numeric_level = getattr(logging, loglevel.upper(), None)
# if not isinstance(numeric_level, int):
#     raise ValueError('Invalid log level: %s' % loglevel)
# logging.basicConfig(
#     format=FORMAT, level=logging.DEBUG if __debug__ else logging.INFO)
# logger = logging.getLogger(__name__)


def transform_phenotype(inv_root, y, fam_indices, null_mean=None):
    """
    Transform phenotype based on inverse square root of phenotypic covariance matrix.
    If the null model included covariates, the fitted mean is removed rather than the overall mean
    """
    # Mean normalise phenotype
    if null_mean is None:
        y = y - np.mean(y)
    else:
        y = y - null_mean
    # Transform by family
    for fam in fam_indices.keys():
        famsize = fam_indices[fam].shape[0]
        if famsize == 1:
            y[fam_indices[fam]] = inv_root[1] * y[fam_indices[fam]]
        else:
            y[fam_indices[fam]] = inv_root[famsize].dot(y[fam_indices[fam]])
    return y


def find_common_ind_ids(obsfiles, impfiles, pheno_ids, from_chr=None, covar=None, impute_unrel=True):
    if len(obsfiles) == 1:
        if from_chr is None:
            raise TypeError('from_chr should not be None.')
        ids, fam_labels = get_ids_with_par(impfiles[from_chr], obsfiles[from_chr], pheno_ids, impute_unrel=impute_unrel)
        return ids, fam_labels
    # for now impfiles and obsfiles contain all chromosome filenames
    # TODO: if from_chr is not 1, keyerror
    ids, fam_labels = get_ids_with_par(impfiles[1], obsfiles[1], pheno_ids, impute_unrel=impute_unrel)
    df = pd.DataFrame({'fam_labels_1': fam_labels}, index=ids)
    # for obs, imp in zip(obsfiles[1:], impfiles[1:]):
    for i in range(2, 23):
        if i not in impfiles or i not in obsfiles:
            continue
        imp = impfiles[i]
        obs = obsfiles[i]
        ids_, fam_labels_ = get_ids_with_par(imp, obs, pheno_ids, impute_unrel=impute_unrel)
        df_ = pd.DataFrame({f'fam_labels_{i}': fam_labels_}, index=ids_)
        # merge on index, i.e., ids
        df = df.merge(df_, how='inner', left_index=True, right_index=True)
    if len(df.index) == 0:
        raise ValueError('No commond ids among chromosome files.')
    logger.info(
        f'{len(df.index)} commond individuals among all chromosome files.')
    # for j in range(1, 22):
    #     if not df[f'fam_labels_{j}'].equals(df[f'fam_labels_{j + 1}']):
    #         raise ValueError(
    #             'fam_labels are not consistent across chromosome files.')
    if covar is not None:
        df_covar = pd.DataFrame(index=covar.ids)
        df = df.merge(df_covar, how='inner', left_index=True, right_index=True)
    df = df.reset_index()
    ids = df['index'].values
    fam_labels = df['fam_labels_1'].values
    logger.info(f'Found {len(ids)} common ids.')
    return ids, fam_labels


def parse_filelist(obsfile, impfile, obsformat, from_chr, to_chr, match_all=True):
    obs_files = {}  # []
    imp_files = {}  # []
    if ('~' in obsfile and '~' in impfile):
        bed_ixes = obsfile.split('~')
        imp_ixes = impfile.split('~')
        start = 1 if match_all else from_chr
        end = 23 if match_all else to_chr
        for i in range(start, end):  # from_chr, to_chr
            obsfile_ = bed_ixes[0]+str(i)+bed_ixes[1]+'.'+obsformat
            impfile_ = imp_ixes[0]+str(i)+imp_ixes[1]+'.hdf5'
            if path.exists(impfile_) and path.exists(obsfile_):
                obs_files[i] = obsfile_
                imp_files[i] = impfile_
                # obs_files.append(obsfile)
                # imp_files.append(impfile)
        logger.info(str(len(imp_files)) +
                    ' matched observed and imputed genotype files found')
    else:
        # [obsfiles+'.'+obsformat]
        obs_files = {from_chr: obsfile+'.'+obsformat}
        imp_files = {from_chr: impfile+'.hdf5'}  # [impfiles+'.hdf5']
    # return np.array(obs_files), np.array(imp_files)
    return obs_files, imp_files


def write_output(chrom, snp_ids, pos, alleles, outprefix, parsum, sib, alpha, alpha_ses, alpha_cov, sigmas, freqs):
    """
    Write fitted SNP effects and other parameters to output HDF5 file.
    """
    logger.info('Writing output to ' + outprefix + '.sumstats.hdf5')
    outfile = h5py.File(outprefix+'.sumstats.hdf5', 'w')
    outbim = np.column_stack((chrom, snp_ids, pos, alleles))
    outfile['bim'] = encode_str_array(outbim)
    X_length = 1
    outcols = ['direct']
    if sib:
        X_length += 1
        outcols.append('sib')
    if parsum:
        X_length += 1
        outcols.append('avg_parental')
    else:
        X_length += 2
        outcols = outcols + ['paternal', 'maternal']
    outfile.create_dataset('estimate_covariance', (snp_ids.shape[0], X_length, X_length), dtype='f', chunks=True,
                           compression='gzip', compression_opts=9)
    outfile.create_dataset('estimate', (snp_ids.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile.create_dataset('estimate_ses', (snp_ids.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile['estimate'][:] = alpha
    outfile['estimate_cols'] = encode_str_array(np.array(outcols))
    outfile['estimate_ses'][:] = alpha_ses
    outfile['estimate_covariance'][:] = alpha_cov
    for i in range(len(sigmas)):
        outfile[f'sigma_{i}'] = sigmas[i]
    outfile['freqs'] = freqs
    outfile.close()


def outarray_effect(est, ses, freqs, vy):
    N_effective = vy/(2*freqs*(1-freqs)*np.power(ses, 2))
    Z = est/ses
    P = -log10(np.exp(1))*chi2.logsf(np.power(Z, 2), 1)
    array_out = np.column_stack((N_effective, est, ses, Z, P))
    array_out = np.round(array_out, decimals=6)
    array_out[:, 0] = np.round(array_out[:, 0], 0)
    return array_out


def write_txt_output(chrom, snp_ids, pos, alleles, outprefix, parsum, sib, alpha, alpha_cov, sigmas, freqs):
    outbim = np.column_stack(
        (chrom, snp_ids, pos, alleles, np.round(freqs, 3)))
    header = ['chromosome', 'SNP', 'pos', 'A1', 'A2', 'freq']
    # Which effects to estimate
    effects = ['direct']
    if sib:
        effects.append('sib')
    if not parsum:
        effects += ['paternal', 'maternal']
    effects += ['avg_parental', 'population']
    effects = np.array(effects)
    if not parsum:
        paternal_index = np.where(effects == 'paternal')[0][0]
        maternal_index = np.where(effects == 'maternal')[0][0]
    avg_par_index = np.where(effects == 'avg_parental')[0][0]
    population_index = avg_par_index+1
    # Get transform matrix
    A = np.zeros((len(effects), alpha.shape[1]))
    A[0:alpha.shape[1], 0:alpha.shape[1]] = np.identity(alpha.shape[1])
    if not parsum:
        A[alpha.shape[1]:(alpha.shape[1]+2), :] = 0.5
        A[alpha.shape[1], 0] = 0
        A[alpha.shape[1]+1, 0] = 1
    else:
        A[alpha.shape[1], :] = 1
    # Transform effects
    alpha = alpha.dot(A.T)
    alpha_ses_out = np.zeros((alpha.shape[0], A.shape[0]))
    corrs = ['r_direct_avg_parental', 'r_direct_population']
    if sib:
        corrs.append('r_direct_sib')
    if not parsum:
        corrs.append('r_paternal_maternal')
    ncor = len(corrs)
    alpha_corr_out = np.zeros((alpha.shape[0], ncor))
    for i in range(alpha_cov.shape[0]):
        alpha_cov_i = A.dot(alpha_cov[i, :, :].dot(A.T))
        alpha_ses_out[i, :] = np.sqrt(np.diag(alpha_cov_i))
        # Direct to average parental
        alpha_corr_out[i, 0] = alpha_cov_i[0, avg_par_index] / \
            (alpha_ses_out[i, 0]*alpha_ses_out[i, avg_par_index])
        # Direct to population
        alpha_corr_out[i, 1] = alpha_cov_i[0, population_index] / \
            (alpha_ses_out[i, 0]*alpha_ses_out[i, population_index])
        # Direct to sib
        if sib:
            alpha_corr_out[i, 2] = alpha_cov_i[0, 1] / \
                (alpha_ses_out[i, 0]*alpha_ses_out[i, 1])
        # Paternal to maternal
        if not parsum:
            alpha_corr_out[i, ncor-1] = alpha_cov_i[paternal_index, maternal_index] / \
                (alpha_ses_out[i, maternal_index]
                 * alpha_ses_out[i, paternal_index])
    # Create output array
    vy = sum(sigmas)
    outstack = [outbim]
    for i in range(len(effects)):
        outstack.append(outarray_effect(
            alpha[:, i], alpha_ses_out[:, i], freqs, vy))
        header += [effects[i]+'_N', effects[i]+'_Beta',
                   effects[i]+'_SE', effects[i]+'_Z', effects[i]+'_log10_P']
    outstack.append(np.round(alpha_corr_out, 6))
    header += corrs
    # Output array
    outarray = np.row_stack((np.array(header), np.column_stack(outstack)))
    logger.info('Writing text output to '+outprefix+'.sumstats.gz')
    np.savetxt(outprefix+'.sumstats.gz', outarray, fmt='%s')


def split_batches(nsnp, niid, nbytes, num_cpus, parsum=False, fit_sib=False):
    max_nbytes = 2147483647
    if nbytes > max_nbytes:
        raise ValueError('Too many bytes to handle for multiprocessing.')
    n_tasks = num_cpus
    f = lambda x, y, z: x * y * 2 / z if parsum is False else x * y * 3 / z
    while int(nsnp / n_tasks) * nbytes > max_nbytes / 2 or f(nsnp, niid, n_tasks) > max_nbytes:
        n_tasks += num_cpus
    return np.array_split(np.arange(nsnp), n_tasks)


_var_dict = {}


def _init_worker(y_, varcomps_, covar_, covar_shape, **kwargs):
    _var_dict['y_'] = y_
    _var_dict['varcomps_'] = varcomps_
    for i in range(len(varcomps_) - 1):
        _var_dict[f'varcomp_data_{i}_'] = kwargs[f'varcomp_data_{i}_']
        _var_dict[f'varcomp_row_ind_{i}_'] = kwargs[f'varcomp_row_ind_{i}_']
        _var_dict[f'varcomp_col_ind_{i}_'] = kwargs[f'varcomp_col_ind_{i}_']
    if covar_ is not None and covar_shape is not None:
        _var_dict['covar_'] = covar_
        _var_dict['covar_shape'] = covar_shape


def process_batch(snp_ids, pheno_ids, pargts_f, gts_f, parsum=False,
                  fit_sib=False, max_missing=5, min_maf=0.01, min_info=0.9, verbose=False, print_sample_info=False,
                  impute_unrel=True, no_mean_normalise=False, common_freqs=False):
    y = np.frombuffer(_var_dict['y_'], dtype='float')
    varcomps = _var_dict['varcomps_']
    # no need to provide residual var arr
    varcomp_arr_lst = tuple(
        (
            np.frombuffer(_var_dict[f'varcomp_data_{i}_'], dtype='float'),
            np.frombuffer(_var_dict[f'varcomp_row_ind_{i}_'], dtype='uint32'),
            np.frombuffer(_var_dict[f'varcomp_col_ind_{i}_'], dtype='uint32'),
        ) for i in range(len(varcomps) - 1)
    )
    covar = np.frombuffer(_var_dict['covar_'], dtype='float').reshape(
        _var_dict['covar_shape']) \
            if 'covar_' in _var_dict else None
    lmm = LinearMixedModel(y, varcomp_arr_lst=varcomp_arr_lst,
                           varcomps=varcomps, covar_X=covar)
    ####### Construct family based genotype matrix #######
    G = get_gts_matrix(pargts_f, gts_f, snp_ids=snp_ids, ids=pheno_ids,
                       parsum=parsum, sib=fit_sib, print_sample_info=print_sample_info, impute_unrel=impute_unrel)
    # Check for empty fam labels
    no_fam = np.array([len(x) == 0 for x in G.fams])
    if np.sum(no_fam) > 0:
        ValueError('No family label from pedigree for some individuals')
    G.compute_freqs(ind=G.gts[:, 1, 0] > -1 if common_freqs else None)
    # impute parental genotypes of unrelated individuals
    if impute_unrel:
        logger.info('imputing unrelated ind parental geno.')
        G = impute_unrel_par_gts(G, sib=fit_sib, parsum=parsum)
    else:
        logger.info('Not imputing unrelated individuals parental geno.')
    #### Filter SNPs ####
    if verbose:
        logger.info('Filtering based on MAF')
    G.filter_maf(min_maf)
    gt_filetype = gts_f.split('.')[1]
    # if gt_filetype == 'bed':
    if verbose:
        logger.info('Filtering based on missingness')
    G.filter_missingness(max_missing)
    if gt_filetype == 'bgen':
        if verbose:
            logger.info('Filtering based on imputation INFO')
        G.filter_info(min_info)
    if verbose:
        logger.info(str(G.shape[2])+' SNPs that pass filters')
    #### Fill NAs ####
    if verbose:
        logger.info('Imputing missing values with population frequencies')
    if not no_mean_normalise:
        NAs = G.fill_NAs(ind=None)
    else:
        def _fill_NAs(G):
            # NAs = np.sum(G.gts.mask, axis=0)
            G.gts[G.gts.mask] = 0
        _fill_NAs(G)
    alpha, alpha_cov, alpha_ses = lmm.fit_snps_eff(G.gts)
    return G.chrom, G.pos, G.alleles, G.freqs, G.sid, alpha, alpha_cov, alpha_ses


######### Command line arguments #########
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'pargts', type=str, help='HDF5 file with imputed parental genotypes (without .hdf5 suffix)')
    parser.add_argument('phenofile', type=str,
                        help='Location of the phenotype file')
    parser.add_argument(
        '--bed', type=str, help='Bed file with observed genotypes (without .bed suffix).', default=None)
    parser.add_argument(
        '--bgen', type=str, help='bgen file with observed genotypes (without .bgen suffix).', default=None)
    parser.add_argument('--outprefix', type=str,
                        help='Location to output association statistic hdf5 file. Outputs text output to outprefix.sumstats.gz and HDF5 output to outprefix.sumstats.hdf5', default='')
    parser.add_argument('--parsum', action='store_true',
                        help='Regress onto proband and sum of parental genotypes (useful when parental genotypes imputed from sibs only)', default=False)
    parser.add_argument('--fit_sib', action='store_true',
                        help='Fit indirect effect from sibling ', default=False)
    parser.add_argument(
        '--covar', type=str, help='Path to file with covariates: plain text file with columns FID, IID, covar1, covar2, ..', default=None)
    parser.add_argument(
        '--fit_res', action='store_true', default=False,
        help='Use residualized phenotypes.'
    )
    parser.add_argument('--phen_index', type=int, help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                        default=1)
    parser.add_argument(
        '--start', type=int, help='Start index of the SNPs to use in the observed genotype file, counting from zero', default=0)
    parser.add_argument(
        '--end', type=int, help='End index of SNPs in the observed genotype file. The script will use SNPs with indices in the range [start,end-1], indexing from zero.', default=None)
    parser.add_argument('--min_maf', type=float,
                        help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)', default=0.01)
    parser.add_argument('--max_missing', type=float,
                        help='Ignore SNPs with greater percent missing calls than max_missing (default 5)', default=5)
    parser.add_argument('--min_info', type=float,
                        help='Ignore SNPs with imputation INFO score below this threshold (default 0.90)', default=0.90)
    parser.add_argument('--missing_char', type=str,
                        help='Missing value string in phenotype file (default NA)', default='NA')
    parser.add_argument('--tau_init', type=float, help='Initial value for ratio between shared family environmental variance and residual variance',
                        default=1)
    parser.add_argument('--output_covar_ests', action='store_true',
                        help='Output null model estimates of covariate fixed effects (default False)', default=False)
    parser.add_argument('--batch_size', type=int,
                        help='Batch size of SNPs to load at a time (reduce to reduce memory requirements)', default=100000)
    parser.add_argument('--no_hdf5_out', action='store_true',
                        help='Suppress HDF5 output of summary statistics', default=False)
    parser.add_argument('--no_txt_out', action='store_true',
                        help='Suppress text output of summary statistics', default=False)

    parser.add_argument('--hapmap_bed', type=str,
                        help='Bed file with observed hapmap3 snps (without suffix, chromosome number should be #).', default=None)
    parser.add_argument('--gcta_path', type=str,
                        help='Path to gcta executable.', default=None)
    parser.add_argument('--grm_path', type=str,
                        help='Path to gcta grm output (without prefix).', default=None)
    parser.add_argument('--plink_path', type=str,
                        help='Path to plink2 executable.', default=None)
    parser.add_argument('--grm_npz_path', type=str,
                        help='Path to subsetted grm npz file.', default=None)

    parser.add_argument('--ibdrel_path', type=str,
                        help='Path to KING IBD segment inference output (without .seg prefix).', default=None)
    
    parser.add_argument('--grm_only', action='store_true', default=False,
                        help='whether to only include grm variance component.')

    parser.add_argument('--sparse_thres', type=float,
                        help='Threshold of GRM/IBD sparsity', default=0.05)

    parser.add_argument('--num_threads', type=int,
                        help='Number of threads numpy uses.', default=None)

    parser.add_argument('--num_cpus', type=int,
                        help='Number of cpus to distribute batches across', default=1)

    parser.add_argument('--from_chr',
                        type=int,
                        help='Which chromosome (>=). Should be used with to_chr parameter.')
    parser.add_argument('--to_chr',
                        type=int,
                        help='Which chromosome (<). Should be used with from_chr parameter.')

    parser.add_argument('--match_all',
                        action='store_true', default=False,
                        help='Match individual IDs from all genotype files or not.')
    
    parser.add_argument('--vc_out',
                        type=str,
                        help='Prefix of output filename for variance component array (without .npy).')
    
    parser.add_argument('--vc_in',
                        type=str,
                        help='Prefix of input filename for variance component array (without .npy).')
    
    parser.add_argument('--impute_unrel',
                        action='store_true', default=False,
                        help='Whether to impute parental genotype of unrelated individuals or not.')

    parser.add_argument('--no_mean_normalise',
                        action='store_true', default=False,
                        help='Whether to mean-normalise genotype or not (for testin purposes).')
                    
    parser.add_argument('--common_freqs',
                        action='store_true', default=False,
                        help='Whether to use common allele frequencies.')

    parser.add_argument('--debug',
                     action='store_true', default=False,
                     help='Debug code in single process mode.')
    
    parser.add_argument('--loglevel',
                        type=str, default='INFO',
                        help='Case insensitive Logging level: INFO, DEBUG, ...')

    args = parser.parse_args()

    loglevel = args.loglevel
    FORMAT = '%(asctime)-15s :: %(levelname)s :: %(name)s :: %(funcName)s :: %(message)s'
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(
        format=FORMAT, level=numeric_level)
    logger = logging.getLogger(__name__)

    if (args.ibdrel_path is not None) + (args.grm_path is not None) + (args.gcta_path is not None) > 1:
        raise parser.error(
            'Only one of ibdrel_path and grm_path/gcta_path should be supplied.')
    if args.ibdrel_path is None and args.grm_path is None and args.gcta_path is None:
        # raise parser.error(
        #     'One of ibdrel_path and grm_path/gcta_path should be supplied.')
        if args.grm_only:
            logger.error(
                'If --grm_only is supplied, must also provide either --ibdrel_path, --grm_path or gcta_path.'
            )
        logger.info('Only sibship variance component is modelled.')
    if args.gcta_path is not None and (args.plink_path is None or args.hapmap_bed is None):
        raise parser.error(
            'Should provide plink_path and hapmap_bed is gcta_path is supplied.')

    if args.num_threads:
        # export OMP_NUM_THREADS=...
        os.environ["OMP_NUM_THREADS"] = str(args.num_threads)
        # export OPENBLAS_NUM_THREADS=...
        os.environ["OPENBLAS_NUM_THREADS"] = str(args.num_threads)
        # export MKL_NUM_THREADS=...
        os.environ["MKL_NUM_THREADS"] = str(args.num_threads)
        # export VECLIB_MAXIMUM_THREADS=...
        os.environ["VECLIB_MAXIMUM_THREADS"] = str(args.num_threads)
        # export NUMEXPR_NUM_THREADS=...
        os.environ["NUMEXPR_NUM_THREADS"] = str(args.num_threads)
    import numpy as np
    from scipy.stats import chi2
    from math import log10

    ######### Read Phenotype ########
    y, pheno_ids = read_phenotype(
        args.phenofile, missing_char=args.missing_char, phen_index=args.phen_index)
    ######## Read covariates ########
    if args.covar is not None:
        logger.info('Reading covariates')
        covariates = read_covariates(
            args.covar, missing_char=args.missing_char)
        # Match to pheno ids
        # covariates.filter_ids(pheno_ids)
    ######## Check for bed/bgen #######
    if args.bed is None and args.bgen is None:
        raise(ValueError('Must supply either bed or bgen file with observed genotypes'))
    if args.bed is not None and args.bgen is not None:
        raise(ValueError('Both --bed and --bgen specified. Please specify one only'))
    ######## Find common ind ids #####
    if args.bed is not None:
        gts_list, pargts_list = parse_filelist(
            args.bed, args.pargts, 'bed', args.from_chr, args.to_chr, args.match_all)
    elif args.bgen is not None:
        gts_list, pargts_list = parse_filelist(
            args.bgen, args.pargts, 'bgen', args.from_chr, args.to_chr, args.match_all)
    # if gts_list.shape[0] == 0:
    if len(gts_list) == 0:
        raise(ValueError('No input genotype files found'))
    # if not gts_list.shape[0] == pargts_list.shape[0]:
    if len(gts_list) != len(pargts_list):
        raise ValueError(
            f'Lists of imputed and observed genotype files not of same length {len(gts_list)} and {len(pargts_list)}')
    if args.covar:
        ids, fam_labels = find_common_ind_ids(gts_list, pargts_list, pheno_ids, from_chr=args.from_chr, covar=covariates, impute_unrel=args.impute_unrel)
    else:
        ids, fam_labels = find_common_ind_ids(gts_list, pargts_list, pheno_ids, from_chr=args.from_chr, covar=None, impute_unrel=args.impute_unrel)
    ########## Construct GRM/IBD and sibship matrix ##########
    # logger.info('Building GRM...')
    if args.ibdrel_path is not None:
        # ids, fam_labels = match_grm_ids(
        #     ids, fam_labels, grm_path=args.ibdrel_path, grm_source='ibdrel')
        id_dict = make_id_dict(ids)
        grm_data, grm_row_ind, grm_col_ind = build_ibdrel_arr(
            args.ibdrel_path, id_dict=id_dict, keep=ids, thres=args.sparse_thres)
    elif args.grm_path is not None or args.gcta_path is not None:
        if args.grm_path is None:
            run_gcta_grm(args.plink_path, args.gcta_path,
                         args.hapmap_bed, args.outprefix, ids)
            grm_path = args.outprefix
        else:
            grm_path = args.grm_path
        ids, fam_labels = match_grm_ids(
            ids, fam_labels, grm_path=grm_path, grm_source='gcta')
        id_dict = make_id_dict(ids)
        if args.grm_npz_path is not None:
            grm_data, grm_row_ind, grm_col_ind = build_grm_arr_from_npz(
                id_filepath=grm_path + '.grm.id', npz_path=args.grm_npz_path,
                ids=ids, id_dict=id_dict)
        else:
            grm_data, grm_row_ind, grm_col_ind = build_grm_arr(
                grm_path, id_dict=id_dict, thres=args.sparse_thres)
    if not args.grm_only:
        sib_data, sib_row_ind, sib_col_ind = build_sib_arr(fam_labels)
    ################## Optimize variance components ##################
    y = match_phenotype_(ids, y, pheno_ids)
    if 'grm_data' in locals() and not args.grm_only:
        varcomp_lst = (
            (grm_data, grm_row_ind, grm_col_ind),
            (sib_data, sib_row_ind, sib_col_ind),
        )
    elif args.grm_only:
        varcomp_lst = (
            (grm_data, grm_row_ind, grm_col_ind),
        )
    else:
        varcomp_lst = (
            (sib_data, sib_row_ind, sib_col_ind),
        )
    if args.vc_in:
        varcomps = tuple(np.load(f'{args.vc_in}.npz')['varcomps'])
        if len(varcomps) != len(varcomp_lst) + 1:
            raise ValueError('Supplied varcomps length does not match the length of ')
    else:
        varcomps = None
    if args.covar is None:
        lmm = LinearMixedModel(y, varcomps=varcomps, varcomp_arr_lst=varcomp_lst, covar_X=None)
    else:
        covariates.filter_ids(ids)
        if args.fit_res:
            alpha = np.linalg.solve(covariates.gts.T.dot(covariates.gts), covariates.gts.T.dot(y))
            aa = covariates.gts.dot(alpha)
            y -= aa
            logger.info(f'--fit_res specified. Phenotypes residualized. Variance of y: {np.var(y)}')
            lmm = LinearMixedModel(y, varcomps=varcomps, varcomp_arr_lst=varcomp_lst, covar_X=None)
        else:
            lmm = LinearMixedModel(y, varcomps=varcomps, varcomp_arr_lst=varcomp_lst, covar_X=covariates.gts)
    if not varcomps:
        logger.info('Optimizing variance components...')
        lmm.scipy_optimize()
    else:
        logger.info('varcomps supplied.')
    if args.vc_out:
        np.savez(f'{args.vc_out}', varcomps=np.array(lmm.varcomps))
        logger.info(f'varcomps saved to {args.vc_out}.npz.')
    logger.info(f'Variance components: {lmm.varcomps}')
    sigmas = lmm.varcomps
    ####### Iterate through all chromosomes #######
    chr_ = int(args.from_chr)
    while chr_ < args.to_chr:
        gts_f = gts_list[chr_]
        pargts_f = pargts_list[chr_]
        if args.bed is not None:
            bed = Bed(gts_f, count_A1=True)
            snp_ids = bed.sid
            pos = bed.pos[:, 2]
            chrom = bed.pos[:, 0]
            alleles = np.loadtxt(gts_f[:-4]+'.bim', dtype=str, usecols=(4, 5))
        elif args.bgen is not None:
            bgen = open_bgen(gts_f, verbose=False)
            snp_ids = bgen.rsids # ids
            pos = np.array(bgen.positions)
            chrom = np.array(bgen.chromosomes)
            # If chromosomse unknown, set to zero
            chrom[[len(x) == 0 for x in chrom]] = 0
            alleles = np.array([x.split(',') for x in bgen.allele_ids])
        ####### Compute batches #######
        logger.info('Found '+str(snp_ids.shape[0])+' variants in '+gts_f)
        if args.end is not None:
            snp_ids = snp_ids[args.start:args.end]
            pos = pos[args.start:args.end]
            chrom = chrom[args.start:args.end]
            alleles = alleles[args.start:args.end]
            logger.info('Using SNPs with indices from ' +
                        str(args.start)+' to '+str(args.end))
        # Remove duplicates
        unique_snps, counts = np.unique(snp_ids, return_counts=True)
        non_duplicate = set(unique_snps[counts == 1])
        if np.sum(counts > 1) > 0:
            logger.info('Removing '+str(np.sum(counts > 1)) +
                        ' duplicate SNP ids')
            not_duplicated = np.array([x in non_duplicate for x in snp_ids])
            snp_ids = snp_ids[not_duplicated]
            pos = pos[not_duplicated]
            chrom = chrom[not_duplicated]
            alleles = alleles[not_duplicated, :]
        snp_dict = make_id_dict(snp_ids)
        alpha_dim = 2
        if args.fit_sib:
            alpha_dim += 1
        if not args.parsum:
            alpha_dim += 1
        # Create output files
        alpha = np.zeros((snp_ids.shape[0], alpha_dim), dtype=np.float32)
        alpha[:] = np.nan
        alpha_cov = np.zeros(
            (snp_ids.shape[0], alpha_dim, alpha_dim), dtype=np.float32)
        alpha_cov[:] = np.nan
        alpha_ses = np.zeros((snp_ids.shape[0], alpha_dim), dtype=np.float32)
        alpha_ses[:] = np.nan
        freqs = np.zeros((snp_ids.shape[0]), dtype=np.float32)
        freqs[:] = np.nan
        nbytes = int((alpha.nbytes + alpha_cov.nbytes +
                     alpha_ses.nbytes + freqs.nbytes) / snp_ids.shape[0])
        # Compute batches
        batches = split_batches(
            snp_ids.shape[0], len(ids), nbytes, args.num_cpus, parsum=args.parsum, fit_sib=args.fit_sib)
        if len(batches) == 1:
            logger.info('Using 1 batch')
        else:
            logger.info('Using '+str(len(batches))+' batches')
        ##############  Process batches of SNPs ##############
        # make shared memory for multiprocessing
        y_ = RawArray('d', y.shape[0])
        # write to shared memory
        y_buffer = np.frombuffer(y_, dtype='float')
        np.copyto(y_buffer, y)
        # varcomp_mat_ = RawArray('f', int(2 * grm.shape[0]))
        # varcomp_mat_buffer = np.frombuffer(
        #     varcomp_mat_, dtype='float').reshape((2, grm.shape[0]))
        # np.copyto(varcomp_mat_buffer, np.array([grm, sib], dtype='float'))
        if args.covar and not args.fit_res:
            covar_ = RawArray(
                'd', int(covariates.gts.shape[0] * covariates.gts.shape[1]))
            covar_buffer = np.frombuffer(covar_, dtype='float').reshape(
                (covariates.gts.shape[0], covariates.gts.shape[1]))
            np.copyto(covar_buffer, covariates.gts)
            covar_shape = covariates.gts.shape
        else: covar_ = None; covar_shape = None
        varcomps_ = lmm.varcomps
        varcomp_dict = {}
        # put variance component data to shared memory
        for i in range(len(varcomps_) - 1):
            l = len(varcomp_lst[i][0])
            # data
            varcomp_dict[f'varcomp_data_{i}_'] = RawArray('d', l)
            varcomp_dict[f'varcomp_data_{i}_buffer'] = np.frombuffer(
                varcomp_dict[f'varcomp_data_{i}_'], dtype='float')
            np.copyto(
                varcomp_dict[f'varcomp_data_{i}_buffer'], varcomp_lst[i][0])
            # row indices
            varcomp_dict[f'varcomp_row_ind_{i}_'] = RawArray(ctypes.c_uint32, l)
            varcomp_dict[f'varcomp_row_ind_{i}_buffer'] = np.frombuffer(
                varcomp_dict[f'varcomp_row_ind_{i}_'], dtype='uint32')
            np.copyto(
                varcomp_dict[f'varcomp_row_ind_{i}_buffer'], varcomp_lst[i][1])
            # col indices
            varcomp_dict[f'varcomp_col_ind_{i}_'] = RawArray(ctypes.c_uint32, l)
            varcomp_dict[f'varcomp_col_ind_{i}_buffer'] = np.frombuffer(
                varcomp_dict[f'varcomp_col_ind_{i}_'], dtype='uint32')
            np.copyto(
                varcomp_dict[f'varcomp_col_ind_{i}_buffer'], varcomp_lst[i][2])
        process_batch_ = partial(process_batch, pheno_ids=ids, pargts_f=pargts_f, gts_f=gts_f,
                                 parsum=args.parsum, fit_sib=args.fit_sib,
                                 max_missing=args.max_missing, min_maf=args.min_maf, min_info=args.min_info,
                                 impute_unrel=args.impute_unrel,
                                 no_mean_normalise=args.no_mean_normalise,
                                 common_freqs=args.common_freqs)
        _init_worker_ = partial(_init_worker, **{k: v for k, v in varcomp_dict.items() if 'buffer' not in k})
        if args.debug:
            _init_worker_(y_, varcomps_, covar_, covar_shape)
            process_batch_(snp_ids[batches[0]])
            exit('Debug finishsed.')
        with Pool(
            processes=args.num_cpus,
            initializer=_init_worker_,
            initargs=(y_, varcomps_, covar_, covar_shape,)
        ) as pool:
            result = pool.map(
                process_batch_, [snp_ids[ind] for ind in batches])
        # varcomp_mat = np.array([grm, sib], dtype='float32').T
        # varcomps = varcomps_
        # covar = covariates.gts
        # result = []
        # for ind in batches:
        #     result.append(process_batch_(snp_ids[ind]))
        for i in range(0, len(result)):
            batch_chrom, batch_pos, batch_alleles, batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses = result[
                i]
            # Fill in fitted SNPs
            batch_indices = np.array([snp_dict[x] for x in batch_snps])
            alpha[batch_indices, :] = batch_alpha
            alpha_cov[batch_indices, :, :] = batch_alpha_cov
            alpha_ses[batch_indices, :] = batch_alpha_ses
            freqs[batch_indices] = batch_freqs
            logger.info('Done batch '+str(i+1)+' out of '+str(len(batches)))
        ######## Save output #########
        if not args.no_hdf5_out:
            write_output(chrom, snp_ids, pos, alleles, args.outprefix.replace('~', str(chr_)), args.parsum, args.fit_sib, alpha, alpha_ses, alpha_cov,
                         sigmas, freqs)
        if not args.no_txt_out:
            write_txt_output(chrom, snp_ids, pos, alleles, args.outprefix.replace('~', str(chr_)), args.parsum, args.fit_sib, alpha, alpha_cov,
                             sigmas, freqs)
        chr_ += 1
