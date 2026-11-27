"""Estimate genotyping error from Mendelian inconsistencies in observed pairs."""

import numpy as np
from numba import njit, prange

from snipar.read.bed import read_PO_pairs_from_bed
from snipar.read.bgen import read_PO_pairs_from_bgen
from snipar.utilities import make_id_dict


class g_error:
    def __init__(self, error_ests, ME, sum_het, sid):
        if not (len(error_ests) == len(ME) == len(sum_het) == len(sid)):
            raise ValueError('Genotyping error arrays must have the same number of SNPs')
        self.error_ests = error_ests
        self.ME = ME
        self.sum_het = sum_het
        self.sid = sid
        self.sid_dict = make_id_dict(sid)

    def bayes_shrink(self, alpha, beta):
        self.error_ests = (self.ME + alpha) / (self.sum_het + beta)


def estimate_genotyping_error_rate(ped, bedfiles=None, bgenfiles=None, min_maf=0.01):
    if (bedfiles is None) == (bgenfiles is None):
        raise ValueError('Must provide BED or BGEN files')
    if bedfiles is not None:
        genome_errors = [mendelian_errors(ped, bedfile=f, min_maf=min_maf) for f in bedfiles]
    else:
        genome_errors = [mendelian_errors(ped, bgenfile=f, min_maf=min_maf) for f in bgenfiles]
    if not genome_errors:
        raise ValueError('No genotype files were provided')
    rates = np.concatenate([chrom.error_ests for chrom in genome_errors])
    exposure = np.concatenate([chrom.sum_het for chrom in genome_errors])
    informative = exposure > 0
    if not np.any(informative):
        raise ValueError('No informative parent-offspring genotypes for genotyping error estimation')
    rates = rates[informative]
    mean = np.mean(rates)
    if mean == 0:
        return mean, None
    overdispersion = np.var(rates) / mean - np.mean(1 / exposure[informative])
    if overdispersion <= 0:
        return mean, None
    beta = 1 / overdispersion
    for chromosome in genome_errors:
        chromosome.bayes_shrink(mean * beta, beta)
    return mean, genome_errors


def mendelian_errors(ped, bedfile=None, bgenfile=None, min_maf=0.01):
    if (bedfile is None) == (bgenfile is None):
        raise ValueError('Must provide exactly one BED or BGEN file')
    if bedfile is not None:
        gts, pair_ped, _ = read_PO_pairs_from_bed(ped, bedfile)
    else:
        gts, pair_ped, _ = read_PO_pairs_from_bgen(ped, bgenfile)
    indices = gts.id_dict
    pairs = np.array([(indices[child], indices[parent])
                      for _, child, father, mother, *_ in pair_ped
                      for parent in (father, mother) if parent in indices], dtype=int)
    gts.filter_maf(min_maf)
    errors, called = _count_pairs(np.asfortranarray(gts.gts.filled(np.nan)), pairs)
    exposure = called * gts.freqs * (1 - gts.freqs)
    rate = np.full(exposure.shape, np.nan)
    np.divide(errors, exposure, out=rate, where=exposure > 0)
    return g_error(rate, errors, exposure, gts.sid)


@njit(parallel=True, cache=True)
def _count_pairs(gts, pair_indices):
    n = gts.shape[1]
    errors = np.zeros(n, dtype=np.int64)
    called = np.zeros(n, dtype=np.int64)
    for j in prange(n):
        for i in range(pair_indices.shape[0]):
            child, parent = pair_indices[i, 0], pair_indices[i, 1]
            a, b = gts[child, j], gts[parent, j]
            if np.isfinite(a) and np.isfinite(b):
                called[j] += 1
                errors[j] += abs(a - b) > 1
    return errors, called


def count_ME(gts, pair_indices):
    return _count_pairs(gts, pair_indices)[0]
