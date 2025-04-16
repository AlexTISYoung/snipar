import snipar.preprocess as preprocess
from pysnptools.snpreader import Bed
from bgen_reader import open_bgen
import snipar.read.bed as bed
import snipar.read.bgen as bgen
import snipar.read.phenotype as phenotype
from snipar.pedigree import get_sibpairs_from_ped
import h5py
import numpy as np
import pandas as pd
import random
from snipar.utilities import convert_str_array, make_id_dict
from typing import Tuple
import logging


logger = logging.getLogger(__name__)

def build_ped_from_par_gts(par_gts_f):
    # Imputed parental file
    par_gts_f_ = h5py.File(par_gts_f, 'r')
    # Genotype file
    ped = convert_str_array(np.array(par_gts_f_['pedigree']))
    ped = ped[1:ped.shape[0], :]
    # Remove control families
    controls = np.array([x[0] == '_' for x in ped[:, 0]])
    return ped[np.logical_not(controls), :], convert_str_array(np.array(par_gts_f_['families']))

def get_gts_matrix(ped=None, imp_fams=None, bedfile=None, bgenfile=None, par_gts_f=None, snp_ids = None, ids = None, parsum=False, sib = False, sib_diff = False, compute_controls = False, include_unrel=False, robust=False, trios_sibs=False, match_snp_ids=False, verbose=False, print_sample_info=False):
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
    if ped is None and par_gts_f is None:
        raise(ValueError('Must provide one of pedigree and imputed parental genotypes file'))
    if bedfile is None and bgenfile is None:
        raise(ValueError('Must provide one bed file or one bgen file'))
    if bedfile is not None and bgenfile is not None:
        raise(ValueError('Must provide one bed file or one bgen file'))
    if ped is None:
        ped, imp_fams = build_ped_from_par_gts(par_gts_f)
    # Control families
    controls = np.array([x[0]=='_' for x in ped[:,0]])
    if par_gts_f is not None:
        par_gts_f = h5py.File(par_gts_f, 'r')
    # Compute genotype matrices
    if bedfile is not None:
        G = [bed.get_gts_matrix_given_ped(ped[np.logical_not(controls),:], imp_fams, bedfile, par_gts_f=par_gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, sib_diff=sib_diff,
                                      parsum=parsum, include_unrel=include_unrel, robust=robust, trios_sibs=trios_sibs,
                                      match_snp_ids=match_snp_ids,
                                      verbose=verbose, print_sample_info = print_sample_info)]
        if compute_controls:
            G.append(bed.get_gts_matrix_given_ped(ped[np.array([x[0:3]=='_p_' for x in ped[:,0]]),], imp_fams, bedfile,
                                                par_gts_f=par_gts_f, snp_ids=snp_ids, ids=ids, sib=sib, sib_diff=sib_diff,
                                                parsum=parsum, include_unrel=include_unrel, robust=robust, trios_sibs=trios_sibs,
                                                match_snp_ids=match_snp_ids,
                                                verbose=verbose, print_sample_info = print_sample_info))
            G.append(
                bed.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_m_' for x in ped[:, 0]]),], imp_fams, bedfile, 
                                            par_gts_f=par_gts_f, snp_ids=snp_ids, ids=ids, sib=sib, sib_diff=sib_diff,
                                            parsum=parsum, include_unrel=include_unrel, robust=robust, trios_sibs=trios_sibs,
                                            match_snp_ids=match_snp_ids,
                                            verbose=verbose, print_sample_info = print_sample_info))
            G.append(
                bed.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_o_' for x in ped[:, 0]]),], imp_fams, bedfile, 
                                            par_gts_f=par_gts_f, snp_ids=snp_ids, ids=ids, sib=sib, sib_diff=sib_diff,
                                            parsum=parsum, include_unrel=include_unrel, robust=robust, trios_sibs=trios_sibs,
                                            match_snp_ids=match_snp_ids,
                                            verbose=verbose, print_sample_info = print_sample_info))
            return G
        else:
            return G[0]
    elif bgenfile is not None:
        G = [bgen.get_gts_matrix_given_ped(ped[np.logical_not(controls),:], imp_fams, bgenfile,
                                                    par_gts_f=par_gts_f,snp_ids=snp_ids, ids=ids, sib=sib, sib_diff=sib_diff,
                                                    parsum=parsum, include_unrel=include_unrel, robust=robust, trios_sibs=trios_sibs,
                                                    match_snp_ids=match_snp_ids,
                                                    verbose=verbose, print_sample_info = print_sample_info)]
        if compute_controls:
            G.append(bgen.get_gts_matrix_given_ped(ped[np.array([x[0:3]=='_p_' for x in ped[:,0]]),],imp_fams, bgenfile,
                                                    par_gts_f=par_gts_f,snp_ids=snp_ids, ids=ids, sib=sib, sib_diff=sib_diff,
                                                    parsum=parsum, include_unrel=include_unrel, robust=robust, trios_sibs=trios_sibs,
                                                    match_snp_ids=match_snp_ids,
                                                    verbose=verbose, print_sample_info = print_sample_info))
            G.append(
                bgen.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_m_' for x in ped[:, 0]]),], imp_fams, bgenfile,
                                                    par_gts_f=par_gts_f,snp_ids=snp_ids, ids=ids, sib=sib, sib_diff=sib_diff,
                                                    parsum=parsum, include_unrel=include_unrel, robust=robust, trios_sibs=trios_sibs,
                                                    match_snp_ids=match_snp_ids,
                                                    verbose=verbose, print_sample_info = print_sample_info))
            G.append(
                bgen.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_o_' for x in ped[:, 0]]),], imp_fams, bgenfile,
                                                    par_gts_f=par_gts_f,snp_ids=snp_ids, ids=ids, sib=sib, sib_diff=sib_diff,
                                                    parsum=parsum, include_unrel=include_unrel, robust=robust, trios_sibs=trios_sibs,
                                                    match_snp_ids=match_snp_ids,
                                                    verbose=verbose, print_sample_info = print_sample_info))
            return G
        else:
            return G[0]


def get_ids_with_par(gts_f: str,
                     ped: np.ndarray,
                     imp_fams: np.ndarray,
                     ids: np.ndarray = None,
                     sib: bool = False,
                     include_unrel: bool = False,
                     ibdrel_path: str = None,
                     return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Find ids with observed/imputed parents and family labels.
    """
    if ibdrel_path is not None:
        king = pd.read_csv(ibdrel_path,
                        sep='\t')[['ID1', 'ID2', 'InfType']]
        king['ID1'] = king['ID1'].astype(str)
        king['ID2'] = king['ID2'].astype(str)

        # remove one individual from each pair of MZ twins
        rm_mz = king[king['InfType'] == 'Dup/MZ'].apply(lambda row: row.ID1, axis=1).to_numpy()
        l = len(ids)
        if len(rm_mz) > 0:
            ids = np.setdiff1d(ids, rm_mz)
            print(f'WARNING: {len(rm_mz)} pairs of Duplicates or MZ twins in pedigree. Removed {l - len(ids)} individuals.')
            
    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f_: Bed = Bed(gts_f, count_A1=True)
        gts_ids = gts_f_.iid[:, 1]
        ids, observed_indices, imp_indices, parcount = preprocess.get_indices_given_ped(ped,                                
                                                                             gts_ids,
                                                                             imp_fams=imp_fams,
                                                                             ids=ids,
                                                                             sib=sib,
                                                                             include_unrel=include_unrel)
        gts_ids = gts_f_.iid[observed_indices, 1]
    elif gts_f[(len(gts_f) - 5):len(gts_f)] == '.bgen':
        gts_f_ = open_bgen(gts_f, verbose=False)
        gts_ids = gts_f_.samples
        ids, observed_indices, imp_indices, parcount = preprocess.get_indices_given_ped(ped,
                                                                             gts_ids,
                                                                             imp_fams=imp_fams,
                                                                             ids=ids,
                                                                             sib=sib,
                                                                             include_unrel=include_unrel)
        gts_ids = gts_ids[observed_indices]
    else:
        raise ValueError('Unknown filetype for observed genotypes file: ' +
                         str(gts_f))


    if imp_fams is not None:
        imp_fams = imp_fams[imp_indices]
    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict, imp_fams)
    n_empty_fams = sum([1 for i in fam_labels if len(i) == 0])


    fam_labels = fam_labels.astype('<U20')
    fam_labels[fam_labels == ''] = np.array([f'_not_{i}_' for i in range(n_empty_fams)], dtype='<U20')

    imp_obs_ids = [i for i in range(len(par_status)) if par_status[i, 0] >= 0 and par_status[i, 1] >= 0]
    imp_only_ids = [i for i in range(len(par_status)) if par_status[i, 0] == par_status[i, 1] == 1]
    obs_only_ids = [i for i in range(len(par_status)) if par_status[i, 0] == par_status[i, 1] == 0]
    # logger.info(f'{len(imp_obs_ids)} individuals with phenotype observations and complete observed/imputed genotype observations.')
    # logger.info(f'{len(imp_only_ids)} individuals with imputed but no observed parental genotypes.')
    # logger.info(f'{sum(par_status.sum(axis=1) == 1)} individuals with one imputed and one observed parental genotypes.')
    # logger.info(f'{len(obs_only_ids)} individuals with both parents observed.')
    print(f'{len(imp_obs_ids)} individuals with phenotype observations and complete observed/imputed genotype observations.')
    print(f'{len(imp_only_ids)} individuals with imputed but no observed parental genotypes.')
    print(f'{sum(par_status.sum(axis=1) == 1)} individuals with one imputed and one observed parental genotypes.')
    print(f'{len(obs_only_ids)} individuals with both parents observed.')

    if include_unrel:
        unrelated_inds = np.where(np.array([par_status[i, 0] == -1 and par_status[i, 1] == -1 for i in range(len(par_status))]))[0]
        # logger.info(f'{len(unrelated_inds)} individuals with both parents missing.')
        print(f'{len(unrelated_inds)} samples without imputed or observed parental genotypes will be included through linear imputation.')
    one_missing_ids = [i for i in range(len(par_status)) if sum(par_status[i, :]) <= 0 and par_status[i, 0] * par_status[i, 1] < 0]
    if len(one_missing_ids) > 0:
        raise RuntimeError(f'{len(one_missing_ids)} individuals have one missing parent and one non-missing parent.')
    
    if return_info:
        if include_unrel:
            return ids, fam_labels, unrelated_inds
        return ids, fam_labels, par_status, gt_indices
    return ids, fam_labels


def get_ids_with_sibs(gts_f: str,
                      ped: np.ndarray, 
                      ids: np.ndarray = None,
                      ibdrel_path: str = None,
                      return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Find ids with sibs and family labels.
    """
    if ibdrel_path is not None:
        king = pd.read_csv(ibdrel_path,
                        sep='\t')[['ID1', 'ID2', 'InfType']]
        king['ID1'] = king['ID1'].astype(str)
        king['ID2'] = king['ID2'].astype(str)

        # remove one individual from each pair of MZ twins
        # rm_mz = king[king['InfType'] == 'Dup/MZ'].apply(lambda row: random.choice([row.ID1, row.ID2]), axis=1).to_numpy()
        rm_mz = king[king['InfType'] == 'Dup/MZ'].apply(lambda row: row.ID1, axis=1).to_numpy()
        l = len(ids)
        if len(rm_mz) > 0:
            ids = np.setdiff1d(ids, rm_mz)
            print(f'WARNING: {len(rm_mz)} pairs of Duplicates or MZ twins in pedigree. Removed {l - len(ids)} individuals.')

    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f_: Bed = Bed(gts_f, count_A1=True)
        gts_ids = gts_f_.iid[:, 1]
        ids, observed_indices = preprocess.get_indices_given_ped_sibs(ped, gts_ids, ids=ids, verbose=False)
        gts_ids = gts_f_.iid[observed_indices, 1]
    elif gts_f[(len(gts_f) - 5):len(gts_f)] == '.bgen':
        gts_f_ = open_bgen(gts_f)
        gts_ids = gts_f_.samples
        ids, observed_indices = preprocess.get_indices_given_ped_sibs(ped, gts_ids, ids=ids, verbose=False)
        gts_ids = gts_ids[observed_indices]
    else:
        raise ValueError('Unknown filetype for observed genotypes file: ' +
                         str(gts_f))

    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    gt_indices, fam_labels = preprocess.find_gts(ids, ped, gts_id_dict)
    # Find which individuals can be used
    none_missing = gt_indices >= 0
    N = np.sum(none_missing)
    if N == 0:
        raise ValueError(
            'No individuals with phenotype observations and genotypes.')
    # Take those that can be used
    # gt_indices = gt_indices[none_missing]
    ids = ids[none_missing]

    n_empty_fams = sum([1 for i in fam_labels if len(i) == 0])

    fam_labels = fam_labels.astype('<U20')
    fam_labels[fam_labels == ''] = np.array([f'_not_{i}_' for i in range(n_empty_fams)], dtype='<U20')

    print(f'Using {len(ids)} samples for the sib-difference estimator')
    if return_info:
        return ids, fam_labels, gt_indices
    return ids, fam_labels

def get_ids_with_trios_sibs(gts_f: str,
                            ped: np.ndarray, 
                            ids: np.ndarray = None,
                            ibdrel_path: str = None,
                            return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Find ids with sibs and both parents, and their family labels.
    """
    if ibdrel_path is not None:
        king = pd.read_csv(ibdrel_path,
                        sep='\t')[['ID1', 'ID2', 'InfType']]
        king['ID1'] = king['ID1'].astype(str)
        king['ID2'] = king['ID2'].astype(str)

        # remove one individual from each pair of MZ twins
        # rm_mz = king[king['InfType'] == 'Dup/MZ'].apply(lambda row: random.choice([row.ID1, row.ID2]), axis=1).to_numpy()
        rm_mz = king[king['InfType'] == 'Dup/MZ'].apply(lambda row: row.ID1, axis=1).to_numpy()
        l = len(ids)
        if len(rm_mz) > 0:
            ids = np.setdiff1d(ids, rm_mz)
            print(f'WARNING: {len(rm_mz)} pairs of Duplicates or MZ twins in pedigree. Removed {l - len(ids)} individuals.')

    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f_: Bed = Bed(gts_f, count_A1=True)
        gts_ids = gts_f_.iid[:, 1]
        ids, observed_indices, trios_indices, sibs_indices = preprocess.get_indices_given_ped_trios_sibs(ped, gts_ids, ids=ids, verbose=False)
        gts_ids = gts_f_.iid[observed_indices, 1]
    elif gts_f[(len(gts_f) - 5):len(gts_f)] == '.bgen':
        gts_f_ = open_bgen(gts_f)
        gts_ids = gts_f_.samples
        ids, observed_indices, trios_indices, sibs_indices = preprocess.get_indices_given_ped_trios_sibs(ped, gts_ids, ids=ids, verbose=False)
        gts_ids = gts_ids[observed_indices]
    else:
        raise ValueError('Unknown filetype for observed genotypes file: ' +
                         str(gts_f))

    if sibs_indices.shape[0] == 0:
        trios_only = True # used to decide dimension of output
    else:
        trios_only = False

    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict)
    # Find which individuals can be used
    none_missing = gt_indices[:, 0] >= 0
    N = np.sum(none_missing)
    if N == 0:
        raise ValueError(
            'No individuals with phenotype observations and genotypes.')
    # Take those that can be used
    gt_indices = gt_indices[none_missing]
    ids = ids[none_missing]

    n_empty_fams = sum([1 for i in fam_labels if len(i) == 0])

    fam_labels = fam_labels.astype('<U20')
    fam_labels[fam_labels == ''] = np.array([f'_not_{i}_' for i in range(n_empty_fams)], dtype='<U20')

    if return_info:
        return ids, fam_labels, gt_indices
    return ids, fam_labels, trios_only
    
def get_ids_with_trios(gts_f: str,
                       ped: np.ndarray, 
                       ids: np.ndarray = None,
                       ibdrel_path: str = None,
                       return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Remove ids with sibs only. Used for the unified estimator when only pedigree file is supplied.
    """
    if ibdrel_path is not None:
        king = pd.read_csv(ibdrel_path,
                        sep='\t')[['ID1', 'ID2', 'InfType']]
        king['ID1'] = king['ID1'].astype(str)
        king['ID2'] = king['ID2'].astype(str)

        # remove one individual from each pair of MZ twins
        # rm_mz = king[king['InfType'] == 'Dup/MZ'].apply(lambda row: random.choice([row.ID1, row.ID2]), axis=1).to_numpy()
        rm_mz = king[king['InfType'] == 'Dup/MZ'].apply(lambda row: row.ID1, axis=1).to_numpy()
        l = len(ids)
        if len(rm_mz) > 0:
            ids = np.setdiff1d(ids, rm_mz)
            print(f'WARNING: {len(rm_mz)} pairs of Duplicates or MZ twins in pedigree. Removed {l - len(ids)} individuals.')

    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f_: Bed = Bed(gts_f, count_A1=True)
        gts_ids = gts_f_.iid[:, 1]
        ids, observed_indices = preprocess.get_indices_given_ped_trios(ped, gts_ids, ids=ids, verbose=False)
        gts_ids = gts_f_.iid[observed_indices, 1]
    elif gts_f[(len(gts_f) - 5):len(gts_f)] == '.bgen':
        gts_f_ = open_bgen(gts_f)
        gts_ids = gts_f_.samples
        ids, observed_indices = preprocess.get_indices_given_ped_trios(ped, gts_ids, ids=ids, verbose=False)
        gts_ids = gts_ids[observed_indices]
    else:
        raise ValueError('Unknown filetype for observed genotypes file: ' +
                         str(gts_f))

    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict)
    # TODO: might suffice to only check if the first column == 0?
    # Find which individuals can be used
    none_missing = np.sum(gt_indices == -1, axis=1)
    none_missing = none_missing == 0
    N = sum(none_missing)
    if N == 0:
        raise ValueError(
            'No individuals with phenotype observations, genotypes and complete parental genotypes.')
    # Take those that can be used
    gt_indices = gt_indices[none_missing]
    ids = ids[none_missing]

    n_empty_fams = sum([1 for i in fam_labels if len(i) == 0])

    fam_labels = fam_labels.astype('<U20')
    fam_labels[fam_labels == ''] = np.array([f'_not_{i}_' for i in range(n_empty_fams)], dtype='<U20')

    if return_info:
        return ids, fam_labels, gt_indices
    return ids, fam_labels