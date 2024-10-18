import numpy as np
from pysnptools.snpreader import Bed
from snipar.utilities import make_id_dict, open_bgen
from snipar.pedigree import find_individuals_with_sibs
from snipar.gtarray import gtarray
# import logging

# logger = logging.getLogger(__name__)

def get_indices_given_ped(ped, gts_ids, imp_fams=None, ids=None, sib=False, include_unrel=False, verbose=True):
    """
    Used in get_gts_matrix_given_ped to get the ids of individuals with observed/imputed parental genotypes and, if sib=True, at least one genotyped sibling.
    It returns those ids along with the indices of the relevant individuals and their first degree relatives in the observed genotypes (observed indices),
    and the indices of the imputed parental genotypes for those individuals.
    """
    # Made dictionary for observed genotypes
    gts_id_dict = make_id_dict(gts_ids)
    # If IDs not provided, use all individuals with observed genotypes
    if ids is None:
        ids = gts_ids
    # Find individuals with genotyped siblings
    if sib:
        # Look in full genotype sample in case some genotyped sibs are not in ids
        # if sib_diff:
        #     all_sibs = True
        ids = gts_ids
        ids = find_individuals_with_sibs(ids, ped, gts_ids, return_ids_only=True)
        if verbose:
            print('Found ' + str(ids.shape[0]) + ' individuals with genotyped siblings')
    ### Find parental status
    if verbose:
        print('Checking for observed/imputed parental genotypes')
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, gts_id_dict, imp_fams=imp_fams)
    # Find which individuals can be used
    if not include_unrel:
        # logger.info('filtering out unrel inds.')
        # print('filtering out unrel inds.')
        none_missing = np.min(gt_indices, axis=1)
        none_missing = none_missing >= 0
        N = np.sum(none_missing)
        if N == 0:
            raise ValueError(
                'No individuals with phenotype observations and complete observed/imputed genotype observations')
        # print(str(N) + ' individuals with phenotype observations and complete observed/imputed genotype observations')
    else:
        none_missing = gt_indices[:, 0] >= 0
        if np.sum(np.min(gt_indices, axis=1) >= 0) == 0:
            raise ValueError(
                'No individuals with phenotype observations and complete observed/imputed genotype observations. The unified estimator cannot be performed due to collinearity.')
        N = np.sum(none_missing)
        if N == 0:
            raise ValueError(
                'No genotyped and phenotyped individual.')
        # print(str(N) + ' individuals with phenotype and genotype observations')
    # Take those that can be used
    gt_indices = gt_indices[none_missing, :]
    par_status = par_status[none_missing, :]
    ids = ids[none_missing]
    parcount = np.sum(par_status==0,axis=1)
    # if verbose:
        # print(str(N) + ' individuals with phenotype observations and complete observed/imputed genotype observations')
        # print(str(np.sum(parcount==0))+' individuals with imputed but no observed parental genotypes')
        # print(str(np.sum(parcount==1))+' individuals with one observed and one imputed parent')
        # print(str(np.sum(parcount==2))+' individuals with both parents observed')
    # Find indices of individuals and their parents in observed genotypes
    observed_indices = np.sort(np.unique(np.hstack((gt_indices[:, 0],
                                                    gt_indices[par_status[:, 0] == 0, 1],
                                                    gt_indices[par_status[:, 1] == 0, 2]))))
    # Get indices of imputed parents
    imp_indices = np.sort(np.unique(np.hstack((gt_indices[par_status[:, 0] == 1, 1],
                                               gt_indices[par_status[:, 1] == 1, 2]))))
    # Return ids with imputed/observed parents
    return ids, observed_indices, imp_indices, parcount

def get_indices_given_ped_sibs(ped, gts_ids, ids=None, verbose=True):
    """
    Used in get_gts_matrix_given_ped to get the ids of individuals with genotyped sibs (might not have phenotype observations).
    """
    # Made dictionary for observed genotypes
    gts_id_dict = make_id_dict(gts_ids)
    # If IDs not provided, use all individuals with observed genotypes
    if ids is None:
        ids = gts_ids
    # ids = gts_ids
    # Find individuals with genotyped siblings
    ids = find_individuals_with_sibs(ids, ped, gts_ids, return_ids_only=True)
    if verbose:
        print('Found ' + str(ids.shape[0]) + ' individuals with genotyped siblings.')
    gt_indices, fam_labels = find_gts(ids, ped, gts_id_dict)
    # Find which individuals can be used
    none_missing = gt_indices >= 0
    N = np.sum(none_missing)
    if N == 0:
        raise ValueError(
            'No individuals with phenotype observations and genotype observations')
    # Take those that can be used
    gt_indices = gt_indices[none_missing]
    ids = ids[none_missing]
    # Find indices of individuals
    observed_indices = np.sort(np.unique(gt_indices))
    # Return ids
    return ids, observed_indices

def get_indices_given_ped_trios(ped, gts_ids, ids=None, verbose=True):
    """
    Used in get_gts_matrix_given_ped to get the ids of individuals with both parents' genotypes.
    """
    # Made dictionary for observed genotypes
    gts_id_dict = make_id_dict(gts_ids)
    # If IDs not provided, use all individuals with observed genotypes
    if ids is None:
        ids = gts_ids
    ### Find parental status
    if verbose:
        print('Checking for observed parental genotypes.')
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, gts_id_dict)
    # Find which individuals can be used
    trios_indices = np.nonzero(np.min(gt_indices, axis=1) >= 0)[0]
    # possible sibs indices might include those with complete parental genotypes
    possible_sibs_indices = np.nonzero(np.logical_and(gt_indices[:, 0] >= 0, np.sum(par_status, axis=1) < 0))[0]
    if possible_sibs_indices.shape[0] > 0:
        ids_with_sibs = find_individuals_with_sibs(ids[possible_sibs_indices], ped, gts_ids, return_ids_only=True)
        sibs_indices = np.nonzero(np.isin(ids, ids_with_sibs))[0]
    else:
        sibs_indices = np.array([], dtype=trios_indices.dtype)
    assert not any(np.isin(sibs_indices, trios_indices))
    print(f'{trios_indices.shape[0]} individuals have complete parental genotypes.')
    if sibs_indices.shape[0] > 0:
        print(f"{sibs_indices.shape[0]} individuals have no complete parental genotypes, but have sib genotypes. The power of the unified estimator can be increased by adding imputed parental genotypes.")
    if trios_indices.shape[0] == 0:
        raise ValueError(
            'No individuals with phenotype observations, complete observed genotype observations, and complete parental genotypes.')
    # Take those that can be used
    gt_indices = gt_indices[trios_indices, :]
    par_status = par_status[trios_indices, :]
    ids = ids[trios_indices]
    parcount = np.sum(par_status==0,axis=1)
    if verbose:
        print(str(np.sum(parcount==0))+' individuals with imputed but no observed parental genotypes.')
        print(str(np.sum(parcount==1))+' individuals with one observed and one imputed parent.')
        print(str(np.sum(parcount==2))+' individuals with both parents observed.')
    # Find indices of individuals and their parents in observed genotypes
    observed_indices = gt_indices[:, 0]
    observed_indices = np.unique(np.hstack((observed_indices,
                                            gt_indices[par_status[:, 0] == 0, 1],
                                            gt_indices[par_status[:, 1] == 0, 2])))
    # Return ids with observed parents
    return ids, observed_indices

def get_indices_given_ped_trios_sibs(ped, gts_ids, ids=None, verbose=True):
    """
    Used in get_gts_matrix_given_ped to get the ids of individuals with complete parental genotypes or with at least one sibs.
    """
    # Made dictionary for observed genotypes
    gts_id_dict = make_id_dict(gts_ids)
    # If IDs not provided, use all individuals with observed genotypes
    if ids is None:
        ids = gts_ids
    ### Find parental status
    if verbose:
        print('Checking for observed/imputed parental genotypes')
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, gts_id_dict)
    # Find which individuals can be used
    trios_indices = np.nonzero(np.min(gt_indices, axis=1) >= 0)[0]
    # possible sibs indices might include those with complete parental genotypes
    possible_sibs_indices = np.nonzero(np.logical_and(gt_indices[:, 0] >= 0, np.sum(par_status, axis=1) < 0))[0]
    if possible_sibs_indices.shape[0] > 0:
        ids_with_sibs = find_individuals_with_sibs(ids[possible_sibs_indices], ped, gts_ids, return_ids_only=True)
        sibs_indices = np.nonzero(np.isin(ids, ids_with_sibs))[0]
    else:
        sibs_indices = np.array([], dtype=trios_indices.dtype)
    assert not any(np.isin(sibs_indices, trios_indices))
    print(f'{trios_indices.shape[0]} individuals have complete parental genotypes.')
    print(f'{sibs_indices.shape[0]} individuals have no complete parental genotypes, but have sib genotypes.')
    if trios_indices.shape[0] == 0 and sibs_indices.shape[0] == 0:
        raise ValueError(
            'No individuals with phenotype observations, complete observed genotype observations, and (complete parental genotypes/at least one sib genotyped).')
    elif trios_indices.shape[0] == 0:
        print('No individuals with phenotype observations, complete observed genotype observations, and complete parental genotype; performing the sib-difference estimator instead.')
    elif sibs_indices.shape[0] == 0:
        print('WARNING: no individuals with phenotype observations, complete observed genotype observations, and sib genotypes; performing complete trio analysis instead.')
    # Take those that can be used
    observed_indices = np.hstack((trios_indices, sibs_indices)) # avoid sorting to keep track of trios and sibs
    gt_indices = gt_indices[observed_indices, :]
    par_status = par_status[observed_indices, :]
    ids = ids[observed_indices]
    parcount = np.sum(par_status==0,axis=1)
    if verbose:
        print(str(len(trios_indices)) + ' individuals with phenotype observations and complete observed/imputed genotype observations')
        print(str(np.sum(parcount==0))+' individuals with imputed but no observed parental genotypes.')
        print(str(np.sum(parcount==1))+' individuals with one observed and one imputed parent.')
        print(str(np.sum(parcount==2))+' individuals with both parents observed.')
    # Find indices of individuals and their parents in observed genotypes
    observed_indices = gt_indices[:, 0]
    trios_indices = observed_indices[:trios_indices.shape[0]]
    sibs_indices = observed_indices[trios_indices.shape[0]:]
    # avoid sorting to allow easy identification of trios and sibs when constructing gts matrix
    observed_indices, idx = np.unique(np.hstack((observed_indices,
                                            gt_indices[par_status[:, 0] == 0, 1],
                                            gt_indices[par_status[:, 1] == 0, 2])), return_index=True)
    # Sort indices to get original order
    observed_indices = observed_indices[np.argsort(idx)]
    # Return ids with observed parents or sibs
    return ids, observed_indices, trios_indices, sibs_indices

def remove_sibs(ped, gts_f, ids, verbose=True):
    """
    Remove ids that have at least one sib but no observed parent.
    """
    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f_ = Bed(gts_f, count_A1=True)
        gts_ids = gts_f_.iid[:, 1]
    elif gts_f[(len(gts_f) - 5):len(gts_f)] == '.bgen':
        gts_f_ = open_bgen(gts_f)
        gts_ids = gts_f_.samples
    # Made dictionary for observed genotypes
    gts_id_dict = make_id_dict(gts_ids)
    # If IDs not provided, use all individuals with observed genotypes
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, gts_id_dict)
    # Find which individuals can be used
    trios_indices = np.nonzero(np.min(gt_indices, axis=1) >= 0)[0]
    # possible sibs indices might include those with complete parental genotypes
    possible_sibs_indices = np.nonzero(np.logical_and(gt_indices[:, 0] >= 0, np.sum(par_status, axis=1) < 0))[0]
    if possible_sibs_indices.shape[0] > 0:
        ids_with_sibs = find_individuals_with_sibs(ids[possible_sibs_indices], ped, gts_ids, return_ids_only=True)
        sibs_indices = np.nonzero(np.isin(ids, ids_with_sibs))[0]
    else:
        sibs_indices = np.array([], dtype=trios_indices.dtype)
    if sibs_indices.shape[0] > 0:
        print(f"{sibs_indices.shape[0]} individuals have no complete parental genotypes, but have sib genotypes. The power of the unified estimator can be increased by adding imputed parental genotypes.")
    # Take those that can be used
    mask = np.array([i not in sibs_indices for i in range(len(ids))])
    ids = ids[mask]
    return ids

def find_par_gts(pheno_ids, ped, gts_id_dict, imp_fams=None):
    """
    Used in get_gts_matrix to find whether individuals have imputed or observed parental genotypes, and to
    find the indices of the observed/imputed parents in the observed/imputed genotype arrays.
    'par_status' codes whether an individual has parents that are observed or imputed or neither.
    'gt_indices' records the relevant index of the parent in the observed/imputed genotype arrays
    'fam_labels' records the family of the individual based on the pedigree
    """
    # Whether mother and father have observed/imputed genotypes
    par_status = np.zeros((pheno_ids.shape[0],2),dtype=int)
    par_status[:] = -1
    # Indices of obsered/imputed genotypes in relevant arrays
    gt_indices = np.zeros((pheno_ids.shape[0],3),dtype=int)
    gt_indices[:] = -1
    ## Build dictionaries
    # Where each individual is in the pedigree
    ped_dict = make_id_dict(ped,1)
    # Where the imputed data is for each family
    if imp_fams is not None:
        fam_dict = make_id_dict(imp_fams)
    # Store family ID of each individual
    fam_labels = np.zeros((pheno_ids.shape[0]),dtype=ped.dtype)
    # Find status and find indices
    for i in range(0,pheno_ids.shape[0]):
        # Find index in genotypes
        if pheno_ids[i] in gts_id_dict:
            gt_indices[i,0] = gts_id_dict[pheno_ids[i]]
        # Find index in pedigree
        if pheno_ids[i] in ped_dict:
            ped_i = ped[ped_dict[pheno_ids[i]], :]
            fam_labels[i] = ped_i[0]
            # Check for observed father
            if ped_i[2] in gts_id_dict:
                gt_indices[i,1] = gts_id_dict[ped_i[2]]
                par_status[i,0] = 0
            # Check for observed mother
            if ped_i[3] in gts_id_dict:
                gt_indices[i, 2] = gts_id_dict[ped_i[3]]
                par_status[i,1] = 0
            # If parent not observed, look for imputation
            if imp_fams is not None:
                if ped_i[0] in fam_dict:
                    imp_index = fam_dict[ped_i[0]]
                    # Check if this is imputation of father, or mother, or both
                    if ped_i[4] == 'False' and not par_status[i,0] == 0:
                        gt_indices[i, 1] = imp_index
                        par_status[i, 0] = 1
                    if ped_i[5] == 'False' and not par_status[i,1] == 0:
                        gt_indices[i, 2] = imp_index
                        par_status[i, 1] = 1
    return par_status, gt_indices, fam_labels

def find_gts(pheno_ids, ped, gts_id_dict):
    """
    Used in get_gts_matrix to find whether individuals have observed genotypes, and to
    find the indices. (used for the sib-difference method)
    'gt_indices' records the relevant index of the proband in the genotype arrays
    'fam_labels' records the family of the individual based on the pedigree
    """
    # Indices of observed genotypes in relevant arrays
    gt_indices = np.zeros((pheno_ids.shape[0]),dtype=int)
    gt_indices[:] = -1
    ## Build dictionaries
    # Where each individual is in the pedigree
    ped_dict = make_id_dict(ped,1)
    # Where the imputed data is for each family
    # Store family ID of each individual
    fam_labels = np.zeros((pheno_ids.shape[0]),dtype=ped.dtype)
    # Find status and find indices
    for i in range(0,pheno_ids.shape[0]):
        # Find index in genotypes
        if pheno_ids[i] in gts_id_dict:
            gt_indices[i] = gts_id_dict[pheno_ids[i]]
        # Find index in pedigree
        if pheno_ids[i] in ped_dict:
            ped_i = ped[ped_dict[pheno_ids[i]], :]
            fam_labels[i] = ped_i[0]
    return gt_indices, fam_labels

def make_gts_matrix(gts, par_status, gt_indices, imp_gts=None, parsum = False):
    """
    Used in get_gts_matrix to construct the family based genotype matrix given
    observed/imputed genotypes. 'gt_indices' has the indices in the observed/imputed genotype arrays;
    and par_status codes whether the parents are observed (0) or imputed (1).
    """
    # if np.min(gt_indices)<0:
    #     raise(ValueError('Missing genotype index'))
    if imp_gts is None:
        if np.max(par_status)==1:
            raise(ValueError('No imputed parental genotypes provided'))
    N = gt_indices.shape[0]
    if parsum:
        gdim = 2
    else:
        gdim = 3
    G = np.full((N,gdim,gts.shape[1]), fill_value=-1, dtype=np.float32)
    # Proband genotypes
    G[:,0,:] = gts[gt_indices[:,0],:]
    # Paternal genotypes
    G[par_status[:, 0] == 0, 1 ,:] = gts[gt_indices[par_status[:, 0] == 0, 1], :]
    if imp_gts is not None:
        G[par_status[:, 0] == 1, 1, :] = imp_gts[gt_indices[par_status[:, 0] == 1, 1], :]
    # Maternal genotypes
    if parsum:
        G[par_status[:, 1] == 0, 1, :] += gts[gt_indices[par_status[:, 1] == 0, 2], :]
        if imp_gts is not None:
            G[par_status[:, 1] == 1, 1, :] += imp_gts[gt_indices[par_status[:, 1] == 1, 2], :]
    else:
        G[par_status[:, 1] == 0, 2, :] = gts[gt_indices[par_status[:, 1] == 0, 2], :]
        if imp_gts is not None:
            G[par_status[:, 1] == 1, 2, :] = imp_gts[gt_indices[par_status[:, 1] == 1, 2], :]
    return G

def make_num_obs_par_al_matrix(num_obs_par_al, par_status,gt_indices, N):
    """
    Used in get_gts_matrix to construct the family based genotype matrix given
    observed/imputed genotypes. 'gt_indices' has the indices in the observed/imputed genotype arrays;
    and par_status codes whether the parents are observed (0) or imputed (1).
    """
    G = np.full((N,num_obs_par_al.shape[1]), fill_value=2, dtype=np.float16)
    # Paternal genotypes
    G[(par_status[:, 0] == 0)&(par_status[:, 1] == 0), :] = 4
    G[(par_status[:, 0] == 1), :] = num_obs_par_al[gt_indices[par_status[:, 0] == 1, 1], :]
    G[(par_status[:, 1] == 1), :] = num_obs_par_al[gt_indices[par_status[:, 1] == 1, 2], :]
    return G

def get_fam_means(ids,ped,gts,gts_ids,remove_proband = True, return_famsizes = False):
    """
    Used in get_gts_matrix to find the mean genotype in each sibship (family) for each SNP or for a PGS.
    The gtarray that is returned is indexed based on the subset of ids provided from sibships of size 2 or greater.
    If remove_proband=True, then the genotype/PGS of the index individual is removed from the fam_mean given for that individual.
    """
    ids, ids_fams, gts_fams = find_individuals_with_sibs(ids, ped, gts_ids)
    fams = np.unique(ids_fams)
    fams_dict = make_id_dict(fams)
    # Compute sums of genotypes in each family
    fam_sums = np.zeros((fams.shape[0],gts.shape[1]),dtype=gts.dtype)
    fam_counts = np.zeros((fams.shape[0]),dtype=int)
    for i in range(0,fams.shape[0]):
        fam_indices = np.where(gts_fams==fams[i])[0]
        fam_sums[i,:] = np.sum(gts[fam_indices,:],axis=0)
        fam_counts[i] = fam_indices.shape[0]
    # Place in vector corresponding to IDs
    if remove_proband:
        gts_id_dict = make_id_dict(gts_ids)
    G_sib = np.zeros((ids.shape[0],gts.shape[1]),dtype = np.float32)
    for i in range(0,ids.shape[0]):
        fam_index = fams_dict[ids_fams[i]]
        G_sib[i,:] = fam_sums[fam_index,:]
        n_i = fam_counts[fam_index]
        if remove_proband:
            G_sib[i,:] = G_sib[i,:] - gts[gts_id_dict[ids[i]],:]
            n_i = n_i-1
        G_sib[i,:] = G_sib[i,:]/float(n_i)
    if return_famsizes:
        return [gtarray(G_sib, ids),fam_counts,fam_sums]
    else:
        return gtarray(G_sib,ids)