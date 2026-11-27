import snipar.preprocess as preprocess
from pysnptools.snpreader import Bed
from snipar.utilities import open_bgen
import snipar.read.bed as bed
import snipar.read.bgen as bgen
import snipar.read.phenotype as phenotype
import h5py
import numpy as np
import pandas as pd
from snipar.utilities import make_id_dict
from typing import Tuple
from contextlib import nullcontext

def build_ped_from_par_gts(par_gts_f, *, include_controls=False):
    """Read the pedigree and family-row identifiers stored with imputed genotypes."""
    with h5py.File(par_gts_f, 'r') as source:
        ped = source['pedigree'].asstr(encoding='utf-8')[1:].astype(str)
        families = source['families'].asstr(encoding='utf-8')[:].astype(str)
    if not include_controls:
        ped = ped[~np.char.startswith(ped[:, 0], '_')]
    return ped, families


def get_gts_matrix(ped=None, imp_fams=None, bedfile=None, bgenfile=None, par_gts_f=None, snp_ids = None, ids = None, parsum=False, sib = False, sib_diff = False, compute_controls = False, include_unrel=False, robust=False, trios_sibs=False, match_snp_ids=False, verbose=False, print_sample_info=False, matched_variants=None):
    """Reads observed and imputed genotypes and constructs a family based genotype matrix for the individuals with
    observed/imputed parental genotypes, and if sib=True, at least one genotyped sibling.

    Args:
        par_gts_f : :class:`str`
            path to HDF5 file with imputed parental genotypes
        gts_f : :class:`str`
            path to bed file with observed genotypes
        snp_ids : :class:`numpy.ndarray`
            If provided, only obtains the subset of SNPs specified that are present in both imputed and observed genotypes
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
    if ped is None and par_gts_f is None:
        raise ValueError('Must provide a pedigree or an imputed parental genotypes file')
    if (bedfile is None) == (bgenfile is None):
        raise ValueError('Must provide exactly one BED or BGEN file')
    genotype_path = bedfile if bedfile is not None else bgenfile
    reader = bed if bedfile is not None else bgen

    with (nullcontext() if par_gts_f is None else h5py.File(par_gts_f, 'r')) as source:
        if ped is None or compute_controls:
            if source is None:
                raise ValueError('Control genotypes require an imputation file with control families')
            ped = source['pedigree'].asstr(encoding='utf-8')[1:].astype(str)
        if imp_fams is None and source is not None and not sib_diff:
            imp_fams = source['families'].asstr(encoding='utf-8')[:].astype(str)

        groups = [~np.char.startswith(ped[:, 0], '_')]
        if compute_controls:
            controls = [np.char.startswith(ped[:, 0], prefix) for prefix in ('_p_', '_m_', '_o_')]
            if not all(np.any(group) for group in controls):
                raise ValueError('Imputation file must contain paternal, maternal, and both-parent control families')
            groups.extend(controls)
        matrices = [
            reader.get_gts_matrix_given_ped(
                ped[group], imp_fams, genotype_path, par_gts_f=source, snp_ids=snp_ids,
                ids=ids, sib=sib, sib_diff=sib_diff, parsum=parsum, include_unrel=include_unrel,
                robust=robust, trios_sibs=trios_sibs, match_snp_ids=match_snp_ids,
                matched_variants=matched_variants, verbose=verbose, print_sample_info=print_sample_info,
            )
            for group in groups
        ]
    return matrices if compute_controls else matrices[0]


def _exclude_duplicate_twins(ids, ibdrel_path):
    if ibdrel_path is None:
        return ids
    king = pd.read_csv(ibdrel_path, sep='\t', usecols=['ID1', 'InfType'],
                       dtype=str, keep_default_na=False)
    remove = king.loc[king['InfType'] == 'Dup/MZ', 'ID1']
    keep = ~np.isin(ids, remove)
    if len(remove):
        print(f'WARNING: {len(remove)} pairs of Duplicates or MZ twins in pedigree. '
              f'Removed {np.count_nonzero(~keep)} individuals.')
    return ids[keep]


def _observed_ids(filename):
    filename = str(filename)
    if filename.endswith('.bed'):
        return Bed(filename, count_A1=True).iid[:, 1]
    if filename.endswith('.bgen'):
        with open_bgen(filename) as reader:
            return np.array(reader.samples, dtype=str)
    raise ValueError(f'Unknown filetype for observed genotypes file: {filename}')


def _fill_unrelated_labels(fam_labels):
    missing = fam_labels == ''
    if not np.any(missing):
        return fam_labels
    # Use labels that cannot collide with genuine pedigree families.
    used = set(fam_labels[~missing])
    replacements = []
    index = 0
    for _ in range(np.count_nonzero(missing)):
        label = f'_not_{index}_'
        while label in used:
            index += 1
            label = f'_not_{index}_'
        replacements.append(label)
        index += 1
    replacements = np.asarray(replacements)
    result = fam_labels.astype(np.promote_types(fam_labels.dtype, replacements.dtype))
    result[missing] = replacements
    return result


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
    gts_ids = _observed_ids(gts_f)
    ids = _exclude_duplicate_twins(gts_ids if ids is None else np.asarray(ids), ibdrel_path)
    ids, observed_indices, imp_indices, _ = preprocess.get_indices_given_ped(
        ped, gts_ids, imp_fams=imp_fams, ids=ids, sib=sib, include_unrel=include_unrel,
    )
    gts_ids = gts_ids[observed_indices]



    if imp_fams is not None:
        imp_fams = imp_fams[imp_indices]
    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict, imp_fams)
    fam_labels = _fill_unrelated_labels(fam_labels)

    missing = np.count_nonzero(par_status < 0, axis=1)
    with_both = np.count_nonzero(missing == 0)
    imputed = np.count_nonzero(par_status == 1, axis=1)
    print(f'{with_both} individuals with phenotype observations and complete observed/imputed genotype observations.')
    print(f'{np.count_nonzero(imputed == 2)} individuals with imputed but no observed parental genotypes.')
    print(f'{np.count_nonzero((missing == 0) & (imputed == 1))} individuals with one imputed and one observed parental genotypes.')
    print(f'{np.count_nonzero((missing == 0) & (imputed == 0))} individuals with both parents observed.')

    if include_unrel:
        unrelated_inds = np.flatnonzero(missing == 2)
        print(f'{len(unrelated_inds)} samples without imputed or observed parental genotypes will be included through linear imputation.')
    one_missing = np.count_nonzero(missing == 1)
    if one_missing:
        raise RuntimeError(f'{one_missing} individuals have one missing parent and one non-missing parent.')
    
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
    gts_ids = _observed_ids(gts_f)
    ids = _exclude_duplicate_twins(gts_ids if ids is None else np.asarray(ids), ibdrel_path)
    ids, observed_indices = preprocess.get_indices_given_ped_sibs(ped, gts_ids, ids=ids, verbose=False)
    gts_ids = gts_ids[observed_indices]


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
    fam_labels = _fill_unrelated_labels(fam_labels[none_missing])

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
    gts_ids = _observed_ids(gts_f)
    ids = _exclude_duplicate_twins(gts_ids if ids is None else np.asarray(ids), ibdrel_path)
    ids, observed_indices, _, sibs_indices = preprocess.get_indices_given_ped_trios_sibs(
        ped, gts_ids, ids=ids, verbose=False,
    )
    gts_ids = gts_ids[observed_indices]


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
    fam_labels = _fill_unrelated_labels(fam_labels[none_missing])

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
    gts_ids = _observed_ids(gts_f)
    ids = _exclude_duplicate_twins(gts_ids if ids is None else np.asarray(ids), ibdrel_path)
    ids, observed_indices = preprocess.get_indices_given_ped_trios(ped, gts_ids, ids=ids, verbose=False)
    gts_ids = gts_ids[observed_indices]


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
    fam_labels = _fill_unrelated_labels(fam_labels[none_missing])

    if return_info:
        return ids, fam_labels, gt_indices
    return ids, fam_labels
