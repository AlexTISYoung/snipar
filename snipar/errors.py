"""Estimate genotyping error from Mendelian inconsistencies in observed pairs."""

import numpy as np
from numba import get_num_threads, njit, prange

from snipar.observed import allele_stats, open_observed
from snipar.read.family import IdIndex
from snipar.resources import plan_batches
from snipar.utilities import make_id_dict


class UninformativeParentOffspringError(ValueError):
    """There are no informative genotyped parent-offspring observations."""


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


def estimate_genotyping_error_rate(ped, bedfiles=None, bgenfiles=None, min_maf=0.01, *,
                                  allow_missing=False):
    if (bedfiles is None) == (bgenfiles is None):
        raise ValueError('Must provide BED or BGEN files')
    files, key = (bedfiles, 'bedfile') if bedfiles is not None else (bgenfiles, 'bgenfile')
    genome_errors = []
    for filename in files:
        try:
            chromosome = mendelian_errors(ped, min_maf=min_maf, **{key: filename})
        except UninformativeParentOffspringError:
            if not allow_missing:
                raise
            chromosome = None
        genome_errors.append(chromosome)
    if not genome_errors:
        raise ValueError('No genotype files were provided')
    informative_chromosomes = [chrom for chrom in genome_errors if chrom is not None]
    if not informative_chromosomes:
        raise UninformativeParentOffspringError('No genotyped parent-offspring pairs for genotyping error estimation')
    rates = np.concatenate([chrom.error_ests for chrom in informative_chromosomes])
    exposure = np.concatenate([chrom.sum_het for chrom in informative_chromosomes])
    informative = exposure > 0
    if not np.any(informative):
        raise UninformativeParentOffspringError('No informative parent-offspring genotypes for genotyping error estimation')
    rates = rates[informative]
    mean = np.mean(rates)
    if mean == 0:
        return mean, None
    overdispersion = np.var(rates) / mean - np.mean(1 / exposure[informative])
    if overdispersion <= 0:
        return mean, None
    beta = 1 / overdispersion
    for chromosome in informative_chromosomes:
        chromosome.bayes_shrink(mean * beta, beta)
    return mean, genome_errors


def mendelian_errors(ped, bedfile=None, bgenfile=None, min_maf=0.01):
    if (bedfile is None) == (bgenfile is None):
        raise ValueError('Must provide exactly one BED or BGEN file')
    with open_observed(bedfile, bgenfile, num_threads=get_num_threads()) as source:
        row = IdIndex(source.ids).lookup(ped[:, 1:4])
        found = (row[:, :1] >= 0) & (row[:, 1:] >= 0) & (row[:, :1] != row[:, 1:]) & (ped[:, 2:4] != '0')
        children, parents = np.where(found)
        if children.size == 0:
            raise UninformativeParentOffspringError('No genotyped parent-offspring pairs for genotyping error estimation')
        pairs = np.column_stack((row[children, 0], row[children, parents + 1]))
        samples = np.unique(pairs)
        pairs = np.searchsorted(samples, pairs)
        m = len(source.sid)
        if not m:
            return g_error(np.empty(0), np.empty(0, dtype=int), np.empty(0), source.sid.copy())
        print(f'{len(pairs)} parent-offspring pairs found in {bedfile if bedfile is not None else bgenfile}')
        if len(pairs) * m < 10 ** 5:
            print('Warning: limited information for estimation of genotyping error probability.')
        _, batch = plan_batches(m, len(samples) * (2 if source.bed else 24),
                                fixed_bytes=source.fixed_bytes)
        errors = np.empty(m, dtype=np.int64)
        exposure = np.empty(m)
        maf_keep = np.empty(m, dtype=bool)
        for start in range(0, m, batch):
            end = min(start + batch, m)
            genotypes = source.read(samples, np.arange(start, end), require_calls=True)
            total, observed = allele_stats(genotypes)
            freqs = total / np.maximum(2 * observed, 1)
            maf_keep[start:end] = (freqs > min_maf) & (freqs < 1 - min_maf)
            errors[start:end], called = _count_pairs(genotypes, pairs)
            exposure[start:end] = called * freqs * (1 - freqs)
            del genotypes
        exposure = exposure[maf_keep]
        errors = errors[maf_keep]
        rate = np.full(exposure.shape, np.nan)
        np.divide(errors, exposure, out=rate, where=exposure > 0)
        return g_error(rate, errors, exposure, source.sid[maf_keep])


@njit(parallel=True, cache=True)
def _count_pairs(gts, pair_indices):
    n = gts.shape[1]
    errors = np.zeros(n, dtype=np.int64)
    called = np.zeros(n, dtype=np.int64)
    for j in prange(n):
        for i in range(pair_indices.shape[0]):
            child, parent = pair_indices[i, 0], pair_indices[i, 1]
            a, b = gts[child, j], gts[parent, j]
            if a >= 0 and b >= 0 and np.isfinite(a) and np.isfinite(b):
                called[j] += 1
                errors[j] += abs(float(a) - float(b)) > 1
    return errors, called


def count_ME(gts, pair_indices):
    return _count_pairs(gts, pair_indices)[0]
