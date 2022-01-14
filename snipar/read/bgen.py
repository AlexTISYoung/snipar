import snipar.preprocess as preprocess
import numpy as np
from snipar.gtarray import gtarray
from bgen_reader import open_bgen
from snipar.utilities import *

def match_observed_and_imputed_snps(gts_f, par_gts_f, snp_ids=None, start=0, end=None):
    """
    Used in get_gts_matrix_given_ped to match observed and imputed SNPs and return SNP information on shared SNPs.
    Removes SNPs that have duplicated SNP ids.
    in_obs_sid contains the SNPs in the imputed genotypes that are present in the observed SNPs
    obs_sid_index contains the index in the observed SNPs of the common SNPs
    """
    # Match SNPs from imputed and observed and restrict to those in list
    if snp_ids is None:
        snp_ids = gts_f.ids
        if end is None:
            end = snp_ids.shape[0]
        snp_ids = snp_ids[start:end]
    # Get bim info
    alleles = np.array([x.split(',') for x in gts_f.allele_ids])
    pos = np.array(gts_f.positions)
    chromosome = np.array(gts_f.chromosomes)
    # Remove duplicate ids
    unique_snps, snp_indices, snp_counts = np.unique(snp_ids, return_index=True, return_counts=True)
    snp_set = set(snp_ids[snp_indices[snp_counts == 1]])
    if len(snp_set) < snp_ids.shape[0]:
        print(str(snp_ids.shape[0]-len(snp_set))+' SNPs with duplicate IDs removed')
    ## Read and match SNP ids
    imp_bim = convert_str_array(np.array(par_gts_f['bim_values']))
    # Find relevant column for SNP ids in imputed data
    imp_bim_cols = convert_str_array(np.array(par_gts_f['bim_columns']))
    if 'rsid' in imp_bim_cols:
        id_col = np.where('rsid' == imp_bim_cols)[0][0]
    elif 'id' in imp_bim_cols:
        id_col = np.where('id' == imp_bim_cols)[0][0]
    else:
        raise(ValueError('Could not find SNP ids in imputed parental genotypes'))
    imp_sid = imp_bim[:, id_col]
    obs_sid = gts_f.ids
    if np.unique(obs_sid).shape[0] == 1:
        obs_sid = gts_f.rsids
    obs_sid_dict = make_id_dict(obs_sid)
    in_obs_sid = np.zeros((imp_sid.shape[0]), dtype=bool)
    obs_sid_index = np.zeros((imp_sid.shape[0]), dtype=int)
    for i in range(0, imp_sid.shape[0]):
        if imp_sid[i] in obs_sid_dict and imp_sid[i] in snp_set:
            in_obs_sid[i] = True
            obs_sid_index[i] = obs_sid_dict[imp_sid[i]]
    if np.sum(in_obs_sid) == 0:
        raise ValueError('No SNPs in common between imputed and observed data')
    obs_sid_index = obs_sid_index[in_obs_sid]
    sid = imp_sid[in_obs_sid]
    alleles = alleles[obs_sid_index, :]
    if 'Chr' in imp_bim_cols:
        chr_col = np.where('Chr' == imp_bim_cols)[0][0]
    else:
        chr_col = 0
    chromosome = imp_bim[in_obs_sid,chr_col]
    pos = pos[obs_sid_index]
    return chromosome, sid, pos, alleles, in_obs_sid, obs_sid_index

def get_gts_matrix_given_ped(ped, par_gts_f, gts_f, snp_ids=None, ids=None, sib=False, parsum=False, start=0, end=None, verbose=False, print_sample_info = False):
    """
    Used in get_gts_matrix: see get_gts_matrix for documentation
    """
    ### Genotype file ###
    gts_f = open_bgen(gts_f, verbose=verbose)
    # get ids of genotypes and make dict
    gts_ids = gts_f.samples
    # Get families with imputed parental genotypes
    fams = convert_str_array(np.array(par_gts_f['families']))
    ### Find ids with observed/imputed parents and indices of those in observed/imputed data
    ids, observed_indices, imp_indices = preprocess.get_indices_given_ped(ped, fams, gts_ids, ids=ids, sib=sib, verbose=print_sample_info)
    ### Match observed and imputed SNPs ###
    if verbose:
        print('Matching observed and imputed SNPs')
    chromosome, sid, pos, alleles, in_obs_sid, obs_sid_index = match_observed_and_imputed_snps(gts_f, par_gts_f, snp_ids=snp_ids, start=start, end=end)
    # Read imputed parental genotypes
    if verbose:
        print('Reading imputed parental genotypes')
    if (imp_indices.shape[0]*in_obs_sid.shape[0]) < (np.sum(in_obs_sid)*fams.shape[0]):
        imp_gts = np.array(par_gts_f['imputed_par_gts'][imp_indices, :])
        imp_gts = imp_gts[:,np.arange(in_obs_sid.shape[0])[in_obs_sid]]
    else:
        imp_gts = np.array(par_gts_f['imputed_par_gts'][:,np.arange(in_obs_sid.shape[0])[in_obs_sid]])
        imp_gts = imp_gts[imp_indices,:]
    fams = fams[imp_indices]
    # Read observed genotypes
    if verbose:
        print('Reading observed genotypes')
    gts = np.sum(gts_f.read((observed_indices,obs_sid_index), np.float32)[:,:,np.array([0,2])],axis=2)
    gts_ids = gts_ids[observed_indices]
    gts_id_dict = make_id_dict(gts_ids)
    # Find indices in reduced data
    par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, fams, gts_id_dict)
    if verbose:
        print('Constructing family based genotype matrix')
    ### Make genotype design matrix
    if sib:
        if parsum:
            G = np.zeros((ids.shape[0], 3, gts.shape[1]), dtype=np.float32)
            G[:, np.array([0, 2]), :] = preprocess.make_gts_matrix(gts, imp_gts, par_status, gt_indices, parsum=parsum)
        else:
            G = np.zeros((ids.shape[0],4,gts.shape[1]), dtype=np.float32)
            G[:,np.array([0,2,3]),:] = preprocess.make_gts_matrix(gts, imp_gts, par_status, gt_indices, parsum=parsum)
        G[:,1,:] = preprocess.get_fam_means(ids, ped, gts, gts_ids, remove_proband=True).gts
    else:
        G = preprocess.make_gts_matrix(gts, imp_gts, par_status, gt_indices, parsum=parsum)
    del gts
    del imp_gts
    return gtarray(G, ids, sid, alleles=alleles, pos=pos, chrom=chromosome, fams=fam_labels, par_status=par_status)


