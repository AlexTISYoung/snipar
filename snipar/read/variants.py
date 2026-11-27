"""Alignment of observed variants to the columns in an imputation file."""
from typing import NamedTuple

import numpy as np


class NoMatchingVariantsError(ValueError):
    """No compatible SNPs remain after aligning observed and imputed genotypes."""


class MatchedVariants(NamedTuple):
    chromosome: np.ndarray
    sid: np.ndarray
    pos: np.ndarray
    alleles: np.ndarray
    flipped: np.ndarray
    imputed: np.ndarray
    observed: np.ndarray

    def subset(self, indices):
        return MatchedVariants(*(column[indices] for column in self))

    def with_imputed_mask(self, size):
        # Preserve the return format of the BED and BGEN matching functions.
        mask = np.zeros(size, dtype=bool)
        mask[self.imputed] = True
        return (*self[:5], mask, self.observed)


def bgen_variant_ids(reader):
    """Choose the SNPIDs used in output, falling back to rsIDs for placeholder SNPIDs."""
    sid = np.asarray(reader.ids)
    rsid = np.asarray(reader.rsids)
    if sid.size and np.all(sid == sid[0]) and (sid.size > 1 or sid[0] in ('', '.')):
        return rsid, None
    return sid, rsid


def _unique_ids(ids):
    if ids.size < 2 or np.all(ids[1:] > ids[:-1]):
        keep = (ids != '') & (ids != '.')
        return ids[keep], np.flatnonzero(keep)
    values, indices, counts = np.unique(ids, return_index=True, return_counts=True)
    keep = (counts == 1) & (values != '') & (values != '.')
    return values[keep], indices[keep]


def _intersect_sorted(a, b):
    """Match sorted unique IDs without sorting a concatenated string array."""
    if len(a) > len(b):
        right, left = _intersect_sorted(b, a)
        return left, right
    pos = np.searchsorted(b, a)
    possible = np.flatnonzero(pos < len(b))
    left = possible[b[pos[possible]] == a[possible]]
    return left, pos[left]


def _match_ids(sid, imp_sid, aliases):
    observed_ids, observed_idx = _unique_ids(sid)
    imputed_ids, imputed_idx = _unique_ids(imp_sid)
    obs, imp = _intersect_sorted(observed_ids, imputed_ids)
    obs, imp = observed_idx[obs], imputed_idx[imp]
    if aliases is not None and obs.size < min(observed_idx.size, imputed_idx.size):
        # Phased BGEN imputation stores rsIDs. Only use an alias when neither
        # record has a competing primary ID match, including a duplicate ID.
        alias_ids, alias_idx = _unique_ids(aliases)
        ai, ii = _intersect_sorted(alias_ids, imputed_ids)
        shared = alias_ids[ai]
        rows = alias_idx[ai]
        unique_observed = np.zeros(len(sid), dtype=bool)
        unique_observed[observed_idx] = True
        keep = unique_observed[rows] & ~np.isin(shared, sid) & ~np.isin(sid[rows], imp_sid)
        obs = np.concatenate((obs, rows[keep]))
        imp = np.concatenate((imp, imputed_idx[ii[keep]]))
    return obs, imp


def match_variants(sid, chromosome, pos, alleles, imputation, snp_ids=None, start=0, end=None, aliases=None):
    """Match unique IDs and allele pairs, in increasing HDF5 column order."""
    columns = list(imputation['bim_columns'].asstr(encoding='utf-8')[:])
    bim = imputation['bim_values'].asstr(encoding='utf-8')
    if 'rsid' in columns:
        imp_sid = bim[:, columns.index('rsid')]
    if 'rsid' not in columns or not np.any((imp_sid != '') & (imp_sid != '.')):
        if 'id' not in columns:
            raise ValueError('Cannot find imputed SNP ids')
        imp_sid = bim[:, columns.index('id')]

    obs, imp = _match_ids(sid, imp_sid.astype(str), aliases)
    if snp_ids is None and (start != 0 or end is not None):
        snp_ids = sid[start:end]
    if snp_ids is not None:
        requested, counts = np.unique(snp_ids, return_counts=True)
        keep = np.isin(sid[obs], requested[counts == 1])
        obs, imp = obs[keep], imp[keep]
    if obs.size == 0:
        raise NoMatchingVariantsError('No SNPs with unique IDs in common between imputed and observed data')
    order = np.argsort(imp)
    obs, imp = obs[order], imp[order]

    if 'allele_ids' in columns:
        imputed_pairs = bim[:, columns.index('allele_ids')][imp]
        a, b = alleles[obs].T
        forward = imputed_pairs == np.strings.add(np.strings.add(a, ','), b)
        reverse = imputed_pairs == np.strings.add(np.strings.add(b, ','), a)
    elif 'allele1' in columns and 'allele2' in columns:
        a1 = bim[:, columns.index('allele1')][imp]
        a2 = bim[:, columns.index('allele2')][imp]
        forward = (a1 == alleles[obs, 0]) & (a2 == alleles[obs, 1])
        reverse = (a1 == alleles[obs, 1]) & (a2 == alleles[obs, 0])
    else:
        raise ValueError('Cannot find imputed allele identifiers')

    keep = forward | reverse
    if not np.any(keep):
        raise NoMatchingVariantsError('No alleles match between observed and imputed genotypes')
    if not np.all(keep):
        print(f'Removing {np.sum(~keep)} SNPs with mismatched alleles between imputed and observed')
    obs, imp = obs[keep], imp[keep]
    chrom = chromosome[obs]
    if 'Chr' in columns and np.any(chrom == ''):
        chrom = np.where(chrom == '', bim[:, columns.index('Chr')][imp], chrom)
    return MatchedVariants(chrom, sid[obs], pos[obs], alleles[obs], (reverse & ~forward)[keep], imp, obs)


def observed_variant_indices(observed_ids, requested, preserve_order):
    """Match only unambiguous IDs, optionally retaining the caller's requested order."""
    unique, indices = _unique_ids(observed_ids)
    if requested is None:
        return np.sort(indices) if preserve_order else indices
    requested = np.asarray(requested)
    if not preserve_order:
        requested, _ = _unique_ids(requested)
    locations = np.searchsorted(unique, requested)
    valid = locations < len(unique)
    locations = locations[valid]
    return indices[locations[unique[locations] == requested[valid]]]


def sibling_sample_indices(sibpairs, id_dict, format_name):
    """Use each observed sibling once even if the pedigree contains several pairs."""
    pairs = np.asarray(sibpairs, dtype=str).reshape(-1, 2) if sibpairs is not None else []
    present = [(first, second) for first, second in pairs if first in id_dict and second in id_dict]
    dropped = len(pairs) - len(present)
    if dropped:
        print(f'{dropped} sibling pairs do not both have genotypes in {format_name}')
    if not present:
        raise ValueError(f'No sibling pairs have both genotypes in {format_name}')
    return np.array(sorted({id_dict[iid] for pair in present for iid in pair}), dtype=int)
