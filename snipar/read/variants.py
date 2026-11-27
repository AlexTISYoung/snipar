"""Alignment of observed variants to the columns in an imputation file."""
from typing import NamedTuple

import numpy as np


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
    values, indices, counts = np.unique(ids, return_index=True, return_counts=True)
    keep = (counts == 1) & (values != '') & (values != '.')
    return values[keep], indices[keep]


def _match_ids(sid, imp_sid, aliases):
    observed_ids, observed_idx = _unique_ids(sid)
    imputed_ids, imputed_idx = _unique_ids(imp_sid)
    _, obs, imp = np.intersect1d(observed_ids, imputed_ids, assume_unique=True, return_indices=True)
    obs, imp = observed_idx[obs], imputed_idx[imp]
    if aliases is not None and obs.size < min(observed_idx.size, imputed_idx.size):
        # Phased BGEN imputation stores rsIDs. Only use an alias when neither
        # record has a competing primary ID match, including a duplicate ID.
        alias_ids, alias_idx = _unique_ids(aliases)
        shared, ai, ii = np.intersect1d(alias_ids, imputed_ids, assume_unique=True, return_indices=True)
        rows = alias_idx[ai]
        keep = np.isin(rows, observed_idx) & ~np.isin(shared, sid) & ~np.isin(sid[rows], imp_sid)
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
        raise ValueError('No SNPs with unique IDs in common between imputed and observed data')
    order = np.argsort(imp)
    obs, imp = obs[order], imp[order]

    if 'allele_ids' in columns:
        imputed_pairs = bim[:, columns.index('allele_ids')][imp]
        forward = np.array([pair == a + ',' + b for pair, (a, b) in zip(imputed_pairs, alleles[obs])])
        reverse = np.array([pair == b + ',' + a for pair, (a, b) in zip(imputed_pairs, alleles[obs])])
    elif 'allele1' in columns and 'allele2' in columns:
        a1 = bim[:, columns.index('allele1')][imp]
        a2 = bim[:, columns.index('allele2')][imp]
        forward = (a1 == alleles[obs, 0]) & (a2 == alleles[obs, 1])
        reverse = (a1 == alleles[obs, 1]) & (a2 == alleles[obs, 0])
    else:
        raise ValueError('Cannot find imputed allele identifiers')

    keep = forward | reverse
    if not np.any(keep):
        raise ValueError('No alleles match between observed and imputed genotypes')
    if not np.all(keep):
        print(f'Removing {np.sum(~keep)} SNPs with mismatched alleles between imputed and observed')
    obs, imp = obs[keep], imp[keep]
    chrom = chromosome[obs]
    if 'Chr' in columns and np.any(chrom == ''):
        chrom = np.where(chrom == '', bim[:, columns.index('Chr')][imp], chrom)
    return MatchedVariants(chrom, sid[obs], pos[obs], alleles[obs], (reverse & ~forward)[keep], imp, obs)
