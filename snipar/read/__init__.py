import snipar.preprocess as preprocess
from pysnptools.snpreader import Bed
from bgen_reader import open_bgen
import snipar.read.bed as bed
import snipar.read.bgen as bgen
import snipar.read.phenotype as phenotype
from snipar.pedigree import get_sibpairs_from_ped
import h5py
import numpy as np
from snipar.utilities import convert_str_array, make_id_dict
from typing import Tuple
import logging


logger = logging.getLogger(__name__)


def get_gts_matrix(ped_f=None, bedfile=None, bgenfile=None, par_gts_f=None, snp_ids = None, ids = None, parsum=False, sib = False, compute_controls = False, include_unrel=False, verbose=False, print_sample_info=False):
    """Reads observed and imputed genotypes and constructs a family based genotype matrix for the individuals with
    observed/imputed parental genotypes, and if sib=True, at least one genotyped sibling.

    Args:
        par_gts_f : :class:`str`
            path to HDF5 file with imputed parental genotypes
        gts_f : :class:`str`
            path to bed file with observed genotypes
        snp_ids : :class:`numpy.ndarray`
            If provided, only obtains the subset of SNPs specificed that are present in both imputed and observed genotypes
        ids : :class:`numpy.ndarray`
            If provided, only obtains the ids with observed genotypes and imputed/observed parental genotypes (and observed sibling genotypes if sib=True)
        sib : :class:`bool`
            Retrieve genotypes for individuals with at least one genotyped sibling along with the average of their siblings' genotypes and observed/imputed parental genotypes. Default False.
        compute_controls : :class:`bool`
            Compute polygenic scores for control families (families with observed parental genotypes set to missing). Default False.
        parsum : :class:`bool`
            Return the sum of maternal and paternal observed/imputed genotypes rather than separate maternal/paternal genotypes. Default False.

    Returns:
        G : :class:`snipar.gtarray`
            Genotype array for the subset of genotyped individuals with complete imputed/obsereved parental genotypes. The array is [N x k x L], where
            N is the number of individuals; k depends on whether sib=True and whether parsum=True; and  L is the number of SNPs. If sib=False and parsum=False,
            then k=3 and this axis indexes individual's genotypes, individual's father's imputed/observed genotypes, individual's mother's imputed/observed genotypes.
            If sib=True and parsum=False, then k=4, and this axis indexes the individual, the sibling, the paternal, and maternal genotypes in that order. If parsum=True and sib=False,
            then k=2, and this axis indexes the individual and sum of paternal and maternal genotypes; etc.
            If compute_controls=True, then a list is returned, where the first element is as above, and the following elements give equivalent genotyping arrays for control families where the mother has been set
            to missing, the father has been set to missing, and both parents have been set to missing.

    """
    ####### Find parental status #######
    if ped_f is None and par_gts_f is None:
        raise(ValueError('Must provide one of pedigree and imputed parental genotypes file'))
    if bedfile is None and bgenfile is None:
        raise(ValueError('Must provide one bed file or one bgen file'))
    if bedfile is not None and bgenfile is not None:
        raise(ValueError('Must provide one bed file or one bgen file'))
    if par_gts_f is not None:
        ### Imputed parental file ###
        par_gts_f = h5py.File(par_gts_f,'r')
        # Get pedigree
        ped = convert_str_array(par_gts_f['pedigree'])
        ped = ped[1:ped.shape[0],:]
    else:
        ped = np.loadtxt(ped_f, dtype=str)
        # Remove rows with missing parents
        sibpairs, ped = get_sibpairs_from_ped(ped)
    # Control families
    controls = np.array([x[0]=='_' for x in ped[:,0]])
    # Compute genotype matrices
    if bedfile is not None:
        G = [bed.get_gts_matrix_given_ped(ped[np.logical_not(controls),:], bedfile, par_gts_f=par_gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, 
                                      parsum=parsum, include_unrel=include_unrel,
                                      verbose=verbose, print_sample_info = print_sample_info)]
        if compute_controls:
            G.append(bed.get_gts_matrix_given_ped(ped[np.array([x[0:3]=='_p_' for x in ped[:,0]]),], bedfile,
                                                par_gts_f=par_gts_f, snp_ids=snp_ids, ids=ids, sib=sib, 
                                                parsum=parsum, include_unrel=include_unrel,
                                                verbose=verbose, print_sample_info = print_sample_info))
            G.append(
                bed.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_m_' for x in ped[:, 0]]),], bedfile, 
                                            par_gts_f=par_gts_f, snp_ids=snp_ids, ids=ids, sib=sib, 
                                            parsum=parsum, include_unrel=include_unrel,
                                            verbose=verbose, print_sample_info = print_sample_info))
            G.append(
                bed.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_o_' for x in ped[:, 0]]),], bedfile, 
                                            par_gts_f=par_gts_f, snp_ids=snp_ids, ids=ids, sib=sib, 
                                            parsum=parsum, include_unrel=include_unrel,
                                            verbose=verbose, print_sample_info = print_sample_info))
            return G
        else:
            return G[0]
    elif bgenfile is not None:
        G = [bgen.get_gts_matrix_given_ped(ped[np.logical_not(controls),:], bgenfile,
                                                    par_gts_f=par_gts_f,snp_ids=snp_ids, ids=ids, sib=sib, 
                                                    parsum=parsum, include_unrel=include_unrel,
                                                    verbose=verbose, print_sample_info = print_sample_info)]
        if compute_controls:
            G.append(bgen.get_gts_matrix_given_ped(ped[np.array([x[0:3]=='_p_' for x in ped[:,0]]),],bgenfile,
                                                    par_gts_f=par_gts_f,snp_ids=snp_ids, ids=ids, sib=sib, 
                                                    parsum=parsum, include_unrel=include_unrel,
                                                    verbose=verbose, print_sample_info = print_sample_info))
            G.append(
                bgen.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_m_' for x in ped[:, 0]]),], bgenfile,
                                                    par_gts_f=par_gts_f,snp_ids=snp_ids, ids=ids, sib=sib, 
                                                    parsum=parsum, include_unrel=include_unrel,
                                                    verbose=verbose, print_sample_info = print_sample_info))
            G.append(
                bgen.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_o_' for x in ped[:, 0]]),], bgenfile,
                                                    par_gts_f=par_gts_f,snp_ids=snp_ids, ids=ids, sib=sib, 
                                                    parsum=parsum, include_unrel=include_unrel,
                                                    verbose=verbose, print_sample_info = print_sample_info))
            return G
        else:
            return G[0]


def get_ids_with_par(gts_f: str,
                     par_gts_f: str,
                     ids: np.ndarray = None,
                     sib: bool = False,
                     include_unrel: bool = False,
                     return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Find ids with observed/imputed parents and family labels.
    """
    # Imputed parental file
    par_gts_f_ = h5py.File(par_gts_f, 'r')
    # Genotype file
    ped = convert_str_array(np.array(par_gts_f_['pedigree']))
    ped = ped[1:ped.shape[0], :]
    # logger.debug(f'ped {len(ped)}')
    # Remove control families
    controls = np.array([x[0] == '_' for x in ped[:, 0]])
    ped = ped[np.logical_not(controls), :]
    # logger.debug(f'ped {len(ped)}')
    # Get families with imputed parental genotypes
    fams = convert_str_array(np.array(par_gts_f_['families']))
    # logger.debug(f'#ids present in both phen and ped: {len(np.intersect1d(ids, ped[:, 1]))}')

    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f_: Bed = Bed(gts_f, count_A1=True)
        gts_ids = gts_f_.iid[:, 1]
        # logger.debug(f'Length of geno: {len(gts_ids)}')
        # logger.debug('#iid' + len(gts_ids))
        # logger.debug(f'#ids present in both phen and geno: {len(np.intersect1d(ids, gts_ids))}')
        # logger.debug(f'#ids present in both ped and geno: {len(np.intersect1d(ped[:, 1], gts_ids))}')
        ids, observed_indices, imp_indices, parcount = preprocess.get_indices_given_ped(ped,                                
                                                                             gts_ids,
                                                                             imp_fams=fams,
                                                                             ids=ids,
                                                                             sib=sib,
                                                                             impute_unrel=include_unrel)
        gts_ids = gts_f_.iid[observed_indices, 1]
        # logger.debug(f'Length of observed geno: {len(gts_ids)}')
    elif gts_f[(len(gts_f) - 5):len(gts_f)] == '.bgen':
        gts_f_ = open_bgen(gts_f)
        gts_ids = gts_f_.samples
        # logger.debug(f'Length of geno: {len(gts_ids)}')
        # logger.debug(f'#ids present in both phen and geno: {len(np.intersect1d(ids, gts_ids))}')
        # logger.debug(f'#ids present in both ped and geno: {len(np.intersect1d(ped[:, 1], gts_ids))}')
        ids, observed_indices, imp_indices, parcount = preprocess.get_indices_given_ped(ped,
                                                                             gts_ids,
                                                                             imp_fams=fams,
                                                                             ids=ids,
                                                                             sib=sib,
                                                                             include_unrel=include_unrel)
        gts_ids = gts_ids[observed_indices]
        # logger.debug(f'Length of observed geno: {len(gts_ids)}')
    else:
        raise ValueError('Unknown filetype for observed genotypes file: ' +
                         str(gts_f))

    fams = fams[imp_indices]
    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict, fams)
    # print(len(par_status), len(observed_indices))
    # print(sum([i == j == -1 for i, j in par_status]))
    # print(sum([i == j == -1 for _,i, j in gt_indices]))
    # print(sum([1 for i in fam_labels if len(i) == 0]))
    # ind = np.array([i for i in range(6797) if par_status[i, 0] == par_status[i, 1] == -1 and len(fam_labels[i]) != 0])
    # print(par_status[ind], fam_labels[ind], gt_indices[ind])
    # f = fam_labels[ind]
    # print([ped[i, :] for i in range(len(ped)) if ped[i, 1] in f])

    # logger.debug(f'Length of par status: {len(par_status)}')
    # logger.debug(f'#missing [F, M]: {(par_status == -1).sum(axis=0)}')
    # logger.debug(f'#missing F and M (from par_status): {sum([i == j == -1 for i, j in par_status])}')
    # logger.debug(f'#missing F and M (fro gt_indices): {sum([i == j == -1 for _,i, j in gt_indices])}')
    # logger.debug(f'#both par geno: {sum([i == j == 0 for _,i, j in gt_indices])}')
    # logger.debug(f'#both par geno: {sum([i == j == 0 for i, j in par_status])}')
    # logger.debug(f'#one par geno: {sum([i + j == 1 for i, j in par_status])}')
    n_empty_fams = sum([1 for i in fam_labels if len(i) == 0])
    # logger.debug(f'#unique fam_labels: {np.unique(fam_labels).__len__()}')
    # logger.debug(f'#empty fam_labels: {n_empty_fams}')

    # logger.debug(fam_labels[fam_labels == ''])
    # ind = fam_labels == ''
    fam_labels = fam_labels.astype('<U20')
    fam_labels[fam_labels == ''] = np.array([f'_not_{i}_' for i in range(n_empty_fams)], dtype='<U20')
    # # ind = np.array([i[4] == 'False' and i[5] == 'False' for i in ped])
    # logger.debug(f'#empty fam_labels after assignment: {sum([1 for i in fam_labels if len(i) == 0])}')
    # logger.debug(f'#unique fam_labels: {np.unique(fam_labels).__len__()}')
    # logger.debug(f'#iids to use {len(ids)}')
    # xx = len([1 for i in fam_labels if '_not_' in i])
    # logger.debug(f'{xx}')

    imp_obs_ids = [i for i in range(len(par_status)) if par_status[i, 0] >= 0 and par_status[i, 1] >= 0]
    imp_only_ids = [i for i in range(len(par_status)) if par_status[i, 0] == par_status[i, 1] == 1]
    obs_only_ids = [i for i in range(len(par_status)) if par_status[i, 0] == par_status[i, 1] == 0]
    logger.info(f'{len(imp_obs_ids)} individuals with phenotype observations and complete observed/imputed genotype observations.')
    logger.info(f'{len(imp_only_ids)} individuals with imputed but no observed parental genotypes.')
    logger.info(f'{len(obs_only_ids)} individuals with both parents observed.')

    if include_unrel:
        # missing_ids = [i for i in range(len(par_status)) if par_status[i, 0] == -1 and par_status[i, 1] == -1]
        missing_ids = np.array([par_status[i, 0] == -1 and par_status[i, 1] == -1 for i in range(len(par_status))])
        # logger.info(f'{len(missing_ids)} individuals with both parents missing.')
        logger.info(f'{sum(missing_ids)} individuals with both parents missing.')
    one_missing_ids = [i for i in range(len(par_status)) if sum(par_status[i, :]) <= 0 and par_status[i, 0] * par_status[i, 1] < 0]
    if len(one_missing_ids) > 0:
        logger.error(f'{len(one_missing_ids)} individuals have one missing parent and one non-missing parent.')
        exit(-1)
    ##### testing code
    # print(np.where(ped[:, 1] == ids[0])[0])
    # ids_id = [np.where(ped[:, 1] == i)[0][0] for i in ids]
    # np.testing.assert_array_equal(ped[ids_id, 0], fam_labels); print('ok', len(ids))
    # print(len([1 for i,j in par_status if i >= 0 and j >= 0]), len(par_status))
    # print(sum([i == j for _,i,j in gt_indices[imp_only_ids]]), len(imp_only_ids))
    # print(sum([i != j for _,i,j in gt_indices[obs_only_ids]]), len(obs_only_ids))

    # uni, unique_fam_inds, cnt = np.unique(fam_labels, return_index=True, return_counts=True)
    # x = list(uni[:5])
    # fam_ids = [i for i in range(len(fam_labels)) if fam_labels[i] in x]
    # ids = ids[fam_ids]
    # fam_labels = fam_labels[fam_ids]
    nnz = 0
    _, unique_fam_inds, cnt = np.unique(fam_labels, return_index=True, return_counts=True)
    for i in cnt:
        nnz += i * (i + 1) /2
    logger.info(f'nnz should be {nnz}')

    # unique_fam_ids = unique_fam_ids[cnt == 1]
    # father_obs = [i for i in range(len(par_status)) if par_status[i, 0] == 0 and par_status[i, 1] == 1]
    # father_obs = [i for i in range(len(par_status)) if par_status[i, 0] == 1 and par_status[i, 1] == 0]
    # indices = np.intersect1d(father_obs, unique_fam_ids)
    # indices = imp_only_ids
    # ids = ids[indices]
    # fam_labels = fam_labels[indices]
    # ids_id = [np.where(ped[:, 1] == i)[0][0] for i in ids]
    # print(sum([i == 'False' and j == 'False' for _, _, _, _, i, j in ped[ids_id, :]]), 'test if the ids truly have imputed father and observed mother')
    
    if return_info:
        return ids, fam_labels, par_status, gt_indices
    return ids, fam_labels