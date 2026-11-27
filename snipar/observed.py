"""Bounded reads for commands that operate on raw observed genotypes."""

from contextlib import contextmanager
from functools import cached_property

import numpy as np
from bed_reader import open_bed
from numba import njit, prange

from snipar.config import nan_integer
from snipar.read.bgen import _first_allele_dosage, bgen_is_phased, read_hard_calls
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

    @cached_property
    def phased(self):
        return bgen_is_phased(self.reader, 0)

    def read(self, samples, variants, *, require_calls=False):
        if self.bed:
            return self.reader.read(index=(samples, variants), dtype='int8', order='F',
                                    num_threads=self.threads)
        if require_calls:
            return read_hard_calls(self.reader, samples, variants, self.phased, nan_integer,
                                   num_threads=self.threads)
        return _first_allele_dosage(self.reader, samples, variants, num_threads=self.threads)


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
