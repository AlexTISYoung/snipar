"""Bounded reads for commands that operate on raw observed genotypes."""

from contextlib import contextmanager

import numpy as np
from bed_reader import open_bed
from numba import njit, prange

from snipar.read.variants import bgen_variant_ids
from snipar.utilities import open_bgen


class _Observation:
    def __init__(self, reader, bed, threads):
        self.reader = reader
        self.bed = bed
        self.threads = threads
        self.ids = np.asarray(reader.iid if bed else reader.samples, dtype=str)
        self.sid = np.asarray(reader.sid if bed else bgen_variant_ids(reader)[0])
        self.fixed_bytes = 0 if bed else 32 * len(self.ids) * (threads or 1)

    def read(self, samples, variants):
        if self.bed:
            return self.reader.read(index=(samples, variants), dtype='int8', order='F',
                                    num_threads=self.threads)
        try:
            # bgen-reader's native library uses worker-local decode buffers.
            probabilities = self.reader.read((samples, variants), dtype=np.float32, num_threads=self.threads)
        except TypeError as error:
            if 'num_threads' not in str(error):
                raise
            probabilities = self.reader.read((samples, variants), dtype=np.float32)
        if probabilities.shape[2] == 4:
            return probabilities[..., 0] + probabilities[..., 2]
        if probabilities.shape[2] == 3:
            return 2 * probabilities[..., 0] + probabilities[..., 1]
        raise ValueError('Expected diploid, biallelic BGEN genotypes')


@contextmanager
def open_observed(bedfile=None, bgenfile=None, *, num_threads=1):
    if (bedfile is None) == (bgenfile is None):
        raise ValueError('Provide exactly one BED or BGEN file')
    if bedfile is not None:
        filename = str(bedfile)
        if not filename.endswith('.bed'):
            filename += '.bed'
        source = open_bed(filename, count_A1=True, num_threads=num_threads)
    else:
        source = open_bgen(bgenfile, verbose=False)
    with source as reader:
        yield _Observation(reader, bedfile is not None, num_threads)


@njit(cache=True, parallel=True)
def allele_stats(genotypes):
    counts = np.zeros(genotypes.shape[1], dtype=np.int64)
    total = np.zeros(genotypes.shape[1])
    for j in prange(genotypes.shape[1]):
        for row in range(genotypes.shape[0]):
            value = genotypes[row, j]
            if value >= 0 and np.isfinite(value):
                total[j] += value
                counts[j] += 1
    return total, counts
