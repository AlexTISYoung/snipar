import h5py
import numpy as np
from bgen_reader import open_bgen
from pysnptools.snpreader import Bed
from scipy.stats import chi2
from math import log10
import snipar.read as read
import snipar.lmm as lmm
from snipar.utilities import *
from numba import njit, prange
from snipar.preprocess import find_par_gts

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

def write_output(chrom, snp_ids, pos, alleles, outfile, parsum, sib, alpha, alpha_ses, alpha_cov, sigma2, tau, freqs):
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
    outfile['sigma2'] = sigma2
    outfile['tau'] = tau
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

def write_txt_output(chrom, snp_ids, pos, alleles, outfile, parsum, sib, alpha, alpha_cov, sigma2, tau, freqs):
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
    vy = (1+1/tau)*sigma2
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

def process_batch(y, pedigree, tau, sigma2, snp_ids=None, bedfile=None, bgenfile=None, par_gts_f=None, parsum=False,
                  fit_sib=False, max_missing=5, min_maf=0.01, verbose=False, print_sample_info=False):
    ####### Construct family based genotype matrix #######
    G = read.get_gts_matrix(ped=pedigree, bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, snp_ids=snp_ids, 
                                ids=y.ids, parsum=parsum, sib=fit_sib, verbose=verbose, print_sample_info=print_sample_info)
    G.compute_freqs()
    #### Filter SNPs ####
    if verbose:
        print('Filtering based on MAF')
    G.filter_maf(min_maf)
    if verbose:
        print('Filtering based on missingness')
    G.filter_missingness(max_missing)
    if verbose:
        print(str(G.shape[2])+' SNPs that pass filters')
    #### Match phenotype ####
    y.filter_ids(G.ids)
    if G.ids.shape[0] > y.ids.shape[0]:
        G.filter_ids(y.ids)
    ##### Transform genotypes ######
    if verbose:
        print('Transforming genotypes')
    null_model = lmm.model(y.gts[:,0], np.ones((y.shape[0], 1)), y.fams)
    L = null_model.sigma_inv_root(tau, sigma2)
    G.diagonalise(L)
    ### Fit models for SNPs ###
    if verbose:
        print('Estimating SNP effects')
    alpha, alpha_cov = fit_models(np.array(y.gts[:,0],dtype=np.float32),G.gts)
    alpha_ses = compute_ses(alpha_cov)
    return G.freqs, G.sid, alpha, alpha_cov, alpha_ses

def process_chromosome(chrom_out, y, pedigree, tau, sigma2, outprefix, bedfile=None, bgenfile=None, par_gts_f=None,
                        fit_sib=False, parsum=False, max_missing=5, min_maf=0.01, batch_size=10000, 
                        no_hdf5_out=False, no_txt_out=False):
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
            print('No individuals with genotyped parents found. Using sum of imputed maternal and paternal genotypes to prevent collinearity.')
            parsum = True
        elif 100 > np.sum(parcount>0) > 0:
            print('Warning: low number of individuals with observed parental genotypes. Consider using the --parsum argument to prevent issues due to collinearity.')
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
    ##############  Process batches of SNPs ##############
    for i in range(0,batch_bounds.shape[0]):
        if i==0:
            print_sample_info = True
            verbose = True
        else:
            print_sample_info = False
            verbose = False
        batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses = process_batch(y, pedigree, 
                    tau, sigma2, snp_ids=snp_ids[batch_bounds[i, 0]:batch_bounds[i, 1]], bedfile=bedfile, bgenfile=bgenfile, 
                    par_gts_f = par_gts_f, parsum=parsum, fit_sib=fit_sib, max_missing=max_missing, min_maf=min_maf,
                    print_sample_info=print_sample_info, verbose=verbose)
        # Fill in fitted SNPs
        batch_indices = np.array([snp_dict[x] for x in batch_snps])
        alpha[batch_indices, :] = batch_alpha
        alpha_cov[batch_indices, :, :] = batch_alpha_cov
        alpha_ses[batch_indices, :] = batch_alpha_ses
        freqs[batch_indices] = batch_freqs
        print('Done batch '+str(i+1)+' out of '+str(batch_bounds.shape[0]))
    ######## Save output #########
    if not no_hdf5_out:
        if chrom_out==0:
            hdf5_outfile = outfile_name(outprefix, '.sumstats.hdf5')
        else:
            hdf5_outfile = outfile_name(outprefix, '.sumstats.hdf5', chrom=chrom_out)
        write_output(chrom, snp_ids, pos, alleles, hdf5_outfile, parsum, fit_sib, alpha, alpha_ses, alpha_cov,
                     sigma2, tau, freqs)
    if not no_txt_out:
        if chrom_out==0:
            txt_outfile = outfile_name(outprefix, '.sumstats.gz')
        else:
            txt_outfile = outfile_name(outprefix, '.sumstats.gz', chrom=chrom_out)
        write_txt_output(chrom, snp_ids, pos, alleles, txt_outfile, parsum, fit_sib, alpha, alpha_cov,
                     sigma2, tau, freqs)
