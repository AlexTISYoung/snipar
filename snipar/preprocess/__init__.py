import numpy as np
from snipar.gtarray import gtarray
from os import path

def make_id_dict(x,col=0):
    """
    Make a dictionary that maps from the values in the given column (col) to their row-index in the input array
    """
    if len(x.shape)>1:
        x = x[:,col]
    id_dict = {}
    for i in range(0,x.shape[0]):
        id_dict[x[i]] = i
    return id_dict

def convert_str_array(x):
    """
    Convert an ascii array to unicode array (UTF-8)
    """
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.decode('UTF-8') for y in x])
    return x_out.reshape(x_shape)

def encode_str_array(x):
    """
    Encode a unicode array as an ascii array
    """
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.encode('ascii') for y in x])
    return x_out.reshape(x_shape)

def parse_filelist(obsfiles, impfiles, obsformat):
    obs_files = []
    imp_files = []
    if '~' in obsfiles and impfiles:
        bed_ixes = obsfiles.split('~')
        imp_ixes = impfiles.split('~')
        for i in range(1,23):
            obsfile = bed_ixes[0]+str(i)+bed_ixes[1]+'.'+obsformat
            impfile = imp_ixes[0]+str(i)+imp_ixes[1]+'.hdf5'
            if path.exists(impfile) and path.exists(obsfile):
                obs_files.append(obsfile)
                imp_files.append(impfile)
        print(str(len(imp_files))+' matched observed and imputed genotype files found')
    else:
            obs_files = [obsfiles+'.'+obsformat]
            imp_files = [impfiles+'.hdf5']
    return np.array(obs_files), np.array(imp_files)

def find_individuals_with_sibs(ids,ped,gts_ids, return_ids_only = False):
    """
    Used in get_gts_matrix and get_fam_means to find the individuals in ids that have genotyped siblings.
    """
    # Find genotyped sibships of size > 1
    ped_dict = make_id_dict(ped, 1)
    ids_in_ped = np.array([x in ped_dict for x in gts_ids])
    gts_fams = np.zeros((gts_ids.shape[0]),dtype=gts_ids.dtype)
    gts_fams[ids_in_ped] = np.array([ped[ped_dict[x], 0] for x in gts_ids[ids_in_ped]])
    fams, counts = np.unique(gts_fams[ids_in_ped], return_counts=True)
    sibships = set(fams[counts > 1])
    # Find individuals with genotyped siblings
    ids_in_ped = np.array([x in ped_dict for x in ids])
    ids = ids[ids_in_ped]
    ids_fams = np.array([ped[ped_dict[x], 0] for x in ids])
    ids_with_sibs = np.array([x in sibships for x in ids_fams])
    ids = ids[ids_with_sibs]
    ids_fams = ids_fams[ids_with_sibs]
    if return_ids_only:
        return ids
    else:
        return ids, ids_fams, gts_fams

def get_fam_means(ids,ped,gts,gts_ids,remove_proband = True, return_famsizes = False):
    """
    Used in get_gts_matrix to find the mean genotype in each sibship (family) for each SNP or for a PGS.
    The gtarray that is returned is indexed based on the subset of ids provided from sibships of size 2 or greater.
    If remove_proband=True, then the genotype/PGS of the index individual is removed from the fam_mean given for that individual.
    """
    ids, ids_fams, gts_fams = find_individuals_with_sibs(ids,ped,gts_ids)
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

def find_par_gts(pheno_ids, ped, fams, gts_id_dict):
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
    fam_dict = make_id_dict(fams)
    # Store family ID of each individual
    fam_labels = np.zeros((pheno_ids.shape[0]),dtype=fams.dtype)
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

def make_gts_matrix(gts,imp_gts,par_status,gt_indices, parsum = False):
    """
    Used in get_gts_matrix to construct the family based genotype matrix given
    observed/imputed genotypes. 'gt_indices' has the indices in the observed/imputed genotype arrays;
    and par_status codes whether the parents are observed (0) or imputed (1).
    """
    if np.min(gt_indices)<0:
        raise(ValueError('Missing genotype index'))
    N = gt_indices.shape[0]
    if parsum:
        gdim = 2
    else:
        gdim = 3
    G = np.zeros((N,gdim,gts.shape[1]),np.float32)
    # Proband genotypes
    G[:,0,:] = gts[gt_indices[:,0],:]
    # Paternal genotypes
    G[par_status[:, 0] == 0, 1 ,:] = gts[gt_indices[par_status[:, 0] == 0, 1], :]
    G[par_status[:, 0] == 1, 1, :] = imp_gts[gt_indices[par_status[:, 0] == 1, 1], :]
    # Maternal genotypes
    if parsum:
        G[par_status[:, 1] == 0, 1, :] += gts[gt_indices[par_status[:, 1] == 0, 2], :]
        G[par_status[:, 1] == 1, 1, :] += imp_gts[gt_indices[par_status[:, 1] == 1, 2], :]
    else:
        G[par_status[:, 1] == 0, 2, :] = gts[gt_indices[par_status[:, 1] == 0, 2], :]
        G[par_status[:, 1] == 1, 2, :] = imp_gts[gt_indices[par_status[:, 1] == 1, 2], :]
    return G

def get_indices_given_ped(ped, fams, gts_ids, ids=None, sib=False, verbose = False):
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
    # Find individuals with siblings
    if sib:
        ids = find_individuals_with_sibs(ids, ped, gts_ids, return_ids_only=True)
        if verbose:
            print('Found ' + str(ids.shape[0]) + ' individuals with genotyped siblings')
    ### Find parental status
    if verbose:
        print('Checking for observed/imputed parental genotypes')
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, fams, gts_id_dict)
    # Find which individuals can be used
    none_missing = np.min(gt_indices, axis=1)
    none_missing = none_missing >= 0
    N = np.sum(none_missing)
    if N == 0:
        raise ValueError(
            'No individuals with phenotype observations and complete observed/imputed genotype observations')
    if verbose:
        print(str(N) + ' individuals with phenotype observations and complete observed/imputed genotypes observations')
    # Take those that can be used
    gt_indices = gt_indices[none_missing, :]
    par_status = par_status[none_missing, :]
    ids = ids[none_missing]
    # Find indices of individuals and their parents in observed genotypes
    observed_indices = np.sort(np.unique(np.hstack((gt_indices[:, 0],
                                                    gt_indices[par_status[:, 0] == 0, 1],
                                                    gt_indices[par_status[:, 1] == 0, 2]))))
    # Get indices of imputed parents
    imp_indices = np.sort(np.unique(np.hstack((gt_indices[par_status[:, 0] == 1, 1],
                                               gt_indices[par_status[:, 1] == 1, 2]))))
    # Return ids with imputed/observed parents
    return ids, observed_indices, imp_indices


