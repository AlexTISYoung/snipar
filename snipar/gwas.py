import h5py
import numpy as np
import pandas as pd
from bgen_reader import open_bgen
from pysnptools.snpreader import Bed
from scipy.stats import chi2
from math import log10
import snipar.read as read
import snipar.lmm as lmm
from snipar.gtarray import impute_missing, impute_unrel_par_gts
from snipar.utilities import *
from numba import njit, prange
from snipar.preprocess import find_par_gts
from typing import List, Union, Dict, Hashable
from multiprocessing import Pool, RawArray
from functools import partial
import ctypes
import logging


logger = logging.getLogger(__name__)


def find_common_ind_ids(obsfiles, impfiles, pheno_ids, from_chr=None, covar=None, keep=None, include_unrel=False):
    if len(obsfiles) == 1:
        if from_chr is None:
            raise TypeError('from_chr should not be None.')
        ids, fam_labels = read.get_ids_with_par(impfiles[from_chr], obsfiles[from_chr], pheno_ids, include_unrel=include_unrel)
        if covar or keep:
            df = pd.DataFrame({f'fam_labels': fam_labels}, index=ids)
            if covar is not None:
                df_covar = pd.DataFrame(index=covar.ids)
                df = df.merge(df_covar, how='inner', left_index=True, right_index=True)
            if keep is not None:
                df_keep = pd.read_csv(keep, sep='\s+', header=None)
                df_keep = df_keep.astype(str)
                df_keep.index = df_keep[0]
                del df_keep[0]
                df = df.merge(df_keep, how='inner', left_index=True, right_index=True)
            df = df.reset_index()
            ids = df['index'].values
            fam_labels = df[f'fam_labels'].values
        # logger.info(f'Found {len(ids)} common ids.')
        return ids, fam_labels
    first = True
    j = from_chr
    for i in range(1, 23):
        if i not in impfiles or i not in obsfiles:
            continue
        imp = impfiles[i]
        obs = obsfiles[i]
        ids_, fam_labels_ = read.get_ids_with_par(imp, obs, pheno_ids, impute_unrel=include_unrel)
        if first:
            df = pd.DataFrame({f'fam_labels_{i}': fam_labels_}, index=ids_)
            first = False
            j = i
        else:
            df_ = pd.DataFrame({f'fam_labels_{i}': fam_labels_}, index=ids_)
            # merge on index, i.e., ids
            df = df.merge(df_, how='inner', left_index=True, right_index=True)
    if len(df.index) == 0:
        raise ValueError('No commond ids among chromosome files.')
    # for j in range(1, 22):
    #     if not df[f'fam_labels_{j}'].equals(df[f'fam_labels_{j + 1}']):
    #         raise ValueError(
    #             'fam_labels are not consistent across chromosome files.')
    if covar is not None:
        df_covar = pd.DataFrame(index=covar.ids)
        df = df.merge(df_covar, how='inner', left_index=True, right_index=True)
    if keep is not None:
        df_keep = pd.DataFrame(keep, header=None)
        df_keep.index = df_keep[0]
        del df_keep[0]
        df = df.merge(df_covar, how='inner', left_index=True, right_index=True)
    df = df.reset_index()
    ids = df['index'].values
    fam_labels = df[f'fam_labels_{j}'].values
    return ids, fam_labels


def transform_phenotype(inv_root, y, fam_indices, null_mean = None):
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

@njit(parallel=True)
def fit_models(y,G):
    alpha = np.zeros((G.shape[2],G.shape[1]),dtype=np.float_)
    alpha_cov = np.zeros((G.shape[2],G.shape[1],G.shape[1]),dtype=np.float_)
    for i in prange(G.shape[2]):
        not_na = np.sum(np.isnan(G[:,:,i]),axis=1)==0
        xtx = G[not_na,:,i].T @ (G[not_na,:,i])
        xty = G[not_na,:,i].T @ y[not_na]
        alpha[i,:] = np.linalg.solve(xtx,xty)
        alpha_cov[i,:,:] = np.linalg.inv(xtx)
    return alpha, alpha_cov

@njit(parallel=True)
def compute_ses(alpha_cov):
    alpha_ses = np.zeros((alpha_cov.shape[0],alpha_cov.shape[1]),dtype=np.float_)
    for i in prange(alpha_cov.shape[0]):
        alpha_ses[i,:] = np.sqrt(np.diag(alpha_cov[i,:,:]))
    return alpha_ses

def write_output(chrom, snp_ids, pos, alleles, outfile, parsum, sib, alpha, alpha_ses, alpha_cov, sigmas, freqs):
    """
    Write fitted SNP effects and other parameters to output HDF5 file.
    """
    print('Writing output to ' + outfile)
    outfile = h5py.File(outfile, 'w')
    outbim = np.column_stack((chrom,snp_ids,pos,alleles))
    outfile['bim'] = encode_str_array(outbim)
    X_length = 1
    outcols = ['direct']
    if sib:
        X_length += 1
        outcols.append('sib')
    if parsum:
        X_length += 1
        outcols.append('avg_NTC')
    else:
        X_length += 2
        outcols = outcols + ['paternal','maternal']
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
    N_effective = vy/(2*freqs*(1-freqs)*np.power(ses,2))
    Z = est/ses
    P = -log10(np.exp(1))*chi2.logsf(np.power(Z,2),1)
    array_out = np.column_stack((N_effective,est,ses,Z,P))
    array_out = np.round(array_out, decimals=6)
    array_out[:,0] = np.round(array_out[:,0], 0)
    return array_out

def write_txt_output(chrom, snp_ids, pos, alleles, outfile, parsum, sib, alpha, alpha_cov, sigmas, freqs):
    outbim = np.column_stack((chrom, snp_ids, pos, alleles,np.round(freqs,3)))
    header = ['chromosome','SNP','pos','A1','A2','freq']
    # Which effects to estimate
    effects = ['direct']
    if sib:
        effects.append('sib')
    if not parsum:
        effects += ['paternal','maternal']
    effects += ['avg_NTC','population']
    effects = np.array(effects)
    if not parsum:
        paternal_index = np.where(effects=='paternal')[0][0]
        maternal_index = np.where(effects=='maternal')[0][0]
    avg_NTC_index = np.where(effects=='avg_NTC')[0][0]
    population_index = avg_NTC_index+1
    # Get transform matrix
    A = np.zeros((len(effects),alpha.shape[1]))
    A[0:alpha.shape[1],0:alpha.shape[1]] = np.identity(alpha.shape[1])
    if not parsum:
        A[alpha.shape[1]:(alpha.shape[1]+2), :] = 0.5
        A[alpha.shape[1], 0] = 0
        A[alpha.shape[1]+1, 0] = 1
    else:
        A[alpha.shape[1], :] = 1
    # Transform effects
    alpha = alpha.dot(A.T)
    alpha_ses_out = np.zeros((alpha.shape[0],A.shape[0]))
    corrs = ['r_direct_avg_NTC','r_direct_population']
    if sib:
        corrs.append('r_direct_sib')
    if not parsum:
        corrs.append('r_paternal_maternal')
    ncor = len(corrs)
    alpha_corr_out = np.zeros((alpha.shape[0],ncor))
    for i in range(alpha_cov.shape[0]):
        alpha_cov_i = A.dot(alpha_cov[i,:,:].dot(A.T))
        alpha_ses_out[i,:] = np.sqrt(np.diag(alpha_cov_i))
        # Direct to average NTC
        alpha_corr_out[i,0] = alpha_cov_i[0,avg_NTC_index]/(alpha_ses_out[i,0]*alpha_ses_out[i,avg_NTC_index])
        # Direct to population
        alpha_corr_out[i,1] = alpha_cov_i[0,population_index]/(alpha_ses_out[i,0]*alpha_ses_out[i,population_index])
        # Direct to sib
        if sib:
            alpha_corr_out[i,2] = alpha_cov_i[0,1]/(alpha_ses_out[i,0]*alpha_ses_out[i,1])
        # Paternal to maternal
        if not parsum:
            alpha_corr_out[i,ncor-1] = alpha_cov_i[paternal_index,maternal_index]/(alpha_ses_out[i,maternal_index]*alpha_ses_out[i,paternal_index])
    # Create output array
    vy = sum(sigmas)
    outstack = [outbim]
    for i in range(len(effects)):
        outstack.append(outarray_effect(alpha[:,i],alpha_ses_out[:,i],freqs,vy))
        header += [effects[i]+'_N',effects[i]+'_Beta',effects[i]+'_SE',effects[i]+'_Z',effects[i]+'_log10_P']
    outstack.append(np.round(alpha_corr_out,6))
    header += corrs
    # Output array
    outarray = np.row_stack((np.array(header),np.column_stack(outstack)))
    print('Writing text output to '+outfile)
    np.savetxt(outfile, outarray, fmt='%s')

def compute_batch_boundaries(snp_ids,batch_size):
    nsnp = snp_ids.shape[0]
    n_blocks = int(np.ceil(float(nsnp)/float(batch_size)))
    block_bounds = np.zeros((n_blocks,2),dtype=int)
    start = 0
    for i in range(n_blocks-1):
        block_bounds[i,0] = start
        block_bounds[i,1] = start+batch_size
        start += batch_size
    block_bounds[n_blocks-1,:] = np.array([start,nsnp])
    return block_bounds


def split_batches(nsnp, niid, nbytes, num_cpus, parsum=False, fit_sib=False):
    max_nbytes = 2147483647
    if nbytes > max_nbytes:
        raise ValueError('Too many bytes to handle for multiprocessing.')
    n_tasks = num_cpus
    f = lambda x, y, z: x * y * 2 / z if parsum is False else x * y * 3 / z
    while int(nsnp / n_tasks) * nbytes > max_nbytes / 4 or f(nsnp, niid, n_tasks) > max_nbytes / 4:
        n_tasks += num_cpus
    return np.array_split(np.arange(nsnp), n_tasks)


_var_dict = {}


def _init_worker(y_, varcomps_, covar_, covar_shape, ped_, ped_shape, **kwargs):
    _var_dict['y_'] = y_
    _var_dict['varcomps_'] = varcomps_
    for i in range(len(varcomps_) - 1):
        _var_dict[f'varcomp_data_{i}_'] = kwargs[f'varcomp_data_{i}_']
        _var_dict[f'varcomp_row_ind_{i}_'] = kwargs[f'varcomp_row_ind_{i}_']
        _var_dict[f'varcomp_col_ind_{i}_'] = kwargs[f'varcomp_col_ind_{i}_']
    if covar_ is not None and covar_shape is not None:
        _var_dict['covar_'] = covar_
        _var_dict['covar_shape'] = covar_shape
    _var_dict['ped_'] = ped_
    _var_dict['ped_shape'] = ped_shape


# def process_batch(y, pedigree, tau, sigma2, snp_ids=None, bedfile=None, bgenfile=None, par_gts_f=None, parsum=False,
#                   fit_sib=False, max_missing=5, min_maf=0.01, verbose=False, print_sample_info=False):
#     ####### Construct family based genotype matrix #######
#     G = read.get_gts_matrix(ped=pedigree, bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, snp_ids=snp_ids, 
#                                 ids=y.ids, parsum=parsum, sib=fit_sib, verbose=verbose, print_sample_info=print_sample_info)
#     G.compute_freqs()
#     #### Filter SNPs ####
#     if verbose:
#         print('Filtering based on MAF')
#     G.filter_maf(min_maf)
#     if verbose:
#         print('Filtering based on missingness')
#     G.filter_missingness(max_missing)
#     if verbose:
#         print(str(G.shape[2])+' SNPs that pass filters')
#     #### Match phenotype ####
#     y.filter_ids(G.ids)
#     if G.ids.shape[0] > y.ids.shape[0]:
#         G.filter_ids(y.ids)
#     ##### Transform genotypes ######
#     if verbose:
#         print('Transforming genotypes')
#     null_model = lmm.model(y.gts[:,0], np.ones((y.shape[0], 1)), y.fams)
#     L = null_model.sigma_inv_root(tau, sigma2)
#     G.diagonalise(L)
#     ### Fit models for SNPs ###
#     if verbose:
#         print('Estimating SNP effects')
#     alpha, alpha_cov = fit_models(np.array(y.gts[:,0],dtype=np.float32),G.gts)
#     alpha_ses = compute_ses(alpha_cov)
#     return G.freqs, G.sid, alpha, alpha_cov, alpha_ses


def process_batch(snp_ids, pheno_ids=None, bedfile=None, bgenfile=None, par_gts_f=None, parsum=False,
                  fit_sib=False, max_missing=5, min_maf=0.01, verbose=False, print_sample_info=False,
                  impute_unrel=False, ignore_na_fams=False, ignore_na_rows=False):
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
    ped = np.frombuffer(_var_dict['ped_'], dtype='str').reshape(
        _var_dict['ped_shape'])
    ####### Construct family based genotype matrix #######
    G = read.get_gts_matrix(ped=ped, bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, snp_ids=snp_ids, 
                            ids=pheno_ids, parsum=parsum, sib=fit_sib,
                            include_unrel=impute_unrel,
                            verbose=verbose, print_sample_info=print_sample_info)
    # Check for empty fam labels
    no_fam = np.array([len(x) == 0 for x in G.fams])
    if np.sum(no_fam) > 0:
        ValueError('No family label from pedigree for some individuals')
    G.compute_freqs()
    # impute parental genotypes of unrelated individuals
    if impute_unrel:
        G = impute_unrel_par_gts(G, sib=fit_sib, parsum=parsum)
    else:
        logger.info('Not imputing unrelated individuals parental geno.')
    #### Filter SNPs ####
    if verbose:
        logger.info('Filtering based on MAF')
    G.filter_maf(min_maf)
    if verbose:
        logger.info()('Filtering based on missingness')
    G.filter_missingness(max_missing)
    if verbose:
        logger.info()(str(G.shape[2])+' SNPs that pass filters')
    #### Fill NAs ####
    if verbose:
        logger.info('Imputing missing values with population frequencies')
    if not ignore_na_fams and not ignore_na_rows:
        # NAs = G.fill_NAs()
        G = impute_missing(G)
    G.mean_normalise()
    ### Fit models for SNPs ###
    model = lmm.LinearMixedModel(y, varcomp_arr_lst=varcomp_arr_lst,
                           varcomps=varcomps, covar_X=covar, add_intercept=True)
    alpha, alpha_cov, alpha_ses = model.fit_snps_eff(G.gts, G.fams, ignore_na_fams=ignore_na_fams, ignore_na_rows=ignore_na_rows)
    return G.freqs, G.sid, alpha, alpha_cov, alpha_ses


def process_chromosome(chrom_out, y, varcomp_lst,
                       pedigree, sigmas, outprefix, covariates=None, bedfile=None, bgenfile=None, par_gts_f=None,
                       fit_sib=False, parsum=False, impute_unrel=False, max_missing=5, min_maf=0.01, batch_size=10000, 
                       no_hdf5_out=False, no_txt_out=False, num_cpus=1, debug=False, ignore_na_fams=False, ignore_na_rows=False):
    ######## Check for bed/bgen #######
    if bedfile is None and bgenfile is None:
        raise(ValueError('Must supply either bed or bgen file with observed genotypes'))
    if bedfile is not None and bgenfile is not None:
        raise(ValueError('Both --bed and --bgen specified. Please specify one only'))
    if bedfile is not None:
        bed = Bed(bedfile,count_A1 = True)
        gts_id_dict = make_id_dict(bed.iid,1)
        snp_ids = bed.sid
        pos = np.array(bed.pos[:,2],dtype=int)
        alleles = np.loadtxt(bedfile.split('.bed')[0]+'.bim',dtype=str,usecols=(4,5))
        chrom = np.array(bed.pos[:,0],dtype=int)
    elif bgenfile is not None:
        bgen = open_bgen(bgenfile, verbose=False)
        gts_id_dict = make_id_dict(bgen.samples)
        snp_ids = bgen.ids
        # If SNP IDs are broken, try rsids
        if np.unique(snp_ids).shape[0] == 1:
            snp_ids = bgen.rsids
        pos = np.array(bgen.positions,dtype=int)
        alleles = np.array([x.split(',') for x in bgen.allele_ids])
        chrom = np.array(bgen.chromosomes,dtype='U2')
        # If chromosomse unknown, set to chromosome inferred from filename
        chrom[[len(x)==0 for x in chrom]] = chrom_out
    # Check for observed parents if not using parsum
    if not parsum:
        par_status, gt_indices, fam_labels = find_par_gts(y.ids, pedigree, gts_id_dict)
        parcount = np.sum(par_status==0,axis=1)
        if np.sum(parcount>0)==0:
            logger.warning('No individuals with genotyped parents found. Using sum of imputed maternal and paternal genotypes to prevent collinearity.')
            parsum = True
        elif 100 > np.sum(parcount>0) > 0:
            logger.warning('Warning: low number of individuals with observed parental genotypes. Consider using the --parsum argument to prevent issues due to collinearity.')
    ####### Compute batches #######
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
    snp_dict = make_id_dict(snp_ids)
    # Compute batches
    batch_bounds = compute_batch_boundaries(snp_ids,batch_size)
    if batch_bounds.shape[0] == 1:
        print('Using 1 batch')
    else:
        print('Using '+str(batch_bounds.shape[0])+' batches')
    alpha_dim = 2
    if fit_sib:
        alpha_dim += 1
    if not parsum:
        alpha_dim += 1
    # Create output files
    alpha = np.zeros((snp_ids.shape[0],alpha_dim),dtype=np.float32)
    alpha[:] = np.nan
    alpha_cov = np.zeros((snp_ids.shape[0], alpha_dim, alpha_dim),dtype=np.float32)
    alpha_cov[:] = np.nan
    alpha_ses = np.zeros((snp_ids.shape[0],alpha_dim),dtype=np.float32)
    alpha_ses[:] = np.nan
    freqs = np.zeros((snp_ids.shape[0]),dtype=np.float32)
    freqs[:] = np.nan


    nbytes = int((alpha.nbytes + alpha_cov.nbytes +
                     alpha_ses.nbytes + freqs.nbytes) / snp_ids.shape[0])
    # Compute batches
    batches = split_batches(
        snp_ids.shape[0], len(y.ids), nbytes, num_cpus, parsum=parsum, fit_sib=fit_sib)
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
    if covariates is not None:
        covar_ = RawArray(
            'd', int(covariates.gts.shape[0] * covariates.gts.shape[1]))
        covar_buffer = np.frombuffer(covar_, dtype='float').reshape(
            (covariates.gts.shape[0], covariates.gts.shape[1]))
        np.copyto(covar_buffer, covariates.gts)
        covar_shape = covariates.gts.shape
    else: covar_ = None; covar_shape = None
    varcomps_ = sigmas
    ped_ = RawArray(
        'd', int(pedigree.shape[0] * pedigree.shape[1]))
    ped_buffer = np.frombuffer(ped_, dtype='str').reshape(
        (pedigree.shape[0], pedigree.shape[1]))
    np.copyto(ped_buffer, pedigree)
    ped_shape = pedigree.shape
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
    
    process_batch_ = partial(process_batch, pheno_ids=y.ids,
                             bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f,
                             parsum=parsum, fit_sib=fit_sib,
                             max_missing=max_missing, min_maf=min_maf,
                             impute_unrel=impute_unrel,
                             ignore_na_fams=ignore_na_fams,
                             ignore_na_rows=ignore_na_rows)
    _init_worker_ = partial(_init_worker, **{k: v for k, v in varcomp_dict.items() if 'buffer' not in k})
    if debug:
        _init_worker_(y_, varcomps_, covar_, covar_shape, ped_, ped_shape)
        batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses = process_batch_(snp_ids[batches[46]])
        exit('Debug finishsed.')
    with Pool(
        processes=num_cpus,
        initializer=_init_worker_,
        initargs=(y_, varcomps_, covar_, covar_shape, ped_, ped_shape)
    ) as pool:
        result = pool.map(
            process_batch_, [snp_ids[ind] for ind in batches], chunksize=1)
    
    for i in range(0, len(result)):
        batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses \
            = result[i]
        # Fill in fitted SNPs
        if len(batch_snps) == 0:
            logger.info('Done batch '+str(i+1)+' out of '+str(len(batches)) + f'#snps:{len(batch_snps)}')
            continue
        batch_indices = np.array([snp_dict[x] for x in batch_snps])
        alpha[batch_indices, :] = batch_alpha
        alpha_cov[batch_indices, :, :] = batch_alpha_cov
        alpha_ses[batch_indices, :] = batch_alpha_ses
        freqs[batch_indices] = batch_freqs
        logger.info('Done batch '+str(i+1)+' out of '+str(len(batches)) + f'#snps:{len(batch_snps)}')
    ##############  Process batches of SNPs ##############
    # for i in range(0,batch_bounds.shape[0]):
    #     if i==0:
    #         print_sample_info = True
    #         verbose = True
    #     else:
    #         print_sample_info = False
    #         verbose = False
    #     batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses = process_batch(y, pedigree, 
    #                 tau, sigma2, snp_ids=snp_ids[batch_bounds[i, 0]:batch_bounds[i, 1]], bedfile=bedfile, bgenfile=bgenfile, 
    #                 par_gts_f = par_gts_f, parsum=parsum, fit_sib=fit_sib, max_missing=max_missing, min_maf=min_maf,
    #                 print_sample_info=print_sample_info, verbose=verbose)
    #     # Fill in fitted SNPs
    #     batch_indices = np.array([snp_dict[x] for x in batch_snps])
    #     alpha[batch_indices, :] = batch_alpha
    #     alpha_cov[batch_indices, :, :] = batch_alpha_cov
    #     alpha_ses[batch_indices, :] = batch_alpha_ses
    #     freqs[batch_indices] = batch_freqs
    #     print('Done batch '+str(i+1)+' out of '+str(batch_bounds.shape[0]))
    ######## Save output #########
    if not no_hdf5_out:
        if chrom_out==0:
            hdf5_outfile = outfile_name(outprefix, '.sumstats.hdf5')
        else:
            hdf5_outfile = outfile_name(outprefix, '.sumstats.hdf5', chrom=chrom_out)
        write_output(chrom, snp_ids, pos, alleles, hdf5_outfile, parsum, fit_sib, alpha, alpha_ses, alpha_cov,
                     sigmas, freqs)
    if not no_txt_out:
        if chrom_out==0:
            txt_outfile = outfile_name(outprefix, '.sumstats.gz')
        else:
            txt_outfile = outfile_name(outprefix, '.sumstats.gz', chrom=chrom_out)
        write_txt_output(chrom, snp_ids, pos, alleles, txt_outfile, parsum, fit_sib, alpha, alpha_cov,
                     sigmas, freqs)
