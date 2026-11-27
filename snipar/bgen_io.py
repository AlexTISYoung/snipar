"""Array access to BGEN files on interpreters without a binary cbgen wheel."""

from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd
from bgen import BgenReader

from snipar.utilities import read_oxford_sample


class BgenFile:
    """Expose the variant metadata and sample × variant × probability slices used by SNIPar."""

    def __init__(self, filename):
        self._reader = BgenReader(str(filename))
        try:
            if self._reader.header.has_sample_ids:
                self.samples = np.asarray(self._reader.samples, dtype=str)
                if not len(self.samples) or self.samples[0] != "sample_0":
                    return
            sample_file = Path(filename).with_suffix(".sample")
            self.samples = read_oxford_sample(sample_file)
            if len(self.samples) != self._reader.header.nsamples:
                raise ValueError("Sample file length differs from the BGEN sample count")
        except BaseException:
            self.close()
            raise

    @cached_property
    def ids(self):
        return np.asarray(self._reader.varids())

    @cached_property
    def rsids(self):
        return np.asarray(self._reader.rsids())

    @cached_property
    def chromosomes(self):
        return np.asarray(self._reader.chroms())

    @cached_property
    def positions(self):
        return np.asarray(self._reader.positions())

    @cached_property
    def allele_ids(self):
        return np.asarray([",".join(self._reader[i].alleles) for i in range(len(self._reader))])

    def variant_metadata(self):
        """Read the identifying fields without decompressing genotype probabilities."""
        records = []
        for i in range(len(self._reader)):
            variant = self._reader[i]
            alleles = variant.alleles
            # Offset of the genotype block after the variable-length BGEN identifier fields.
            strings = (variant.varid, variant.rsid, variant.chrom, *alleles)
            header_size = (12 if variant.layout == 2 else 14) + 4 * len(alleles)
            vaddr = variant.fileoffset + header_size + sum(len(s.encode()) for s in strings)
            records.append((variant.varid, variant.rsid, variant.chrom, variant.pos,
                            len(alleles), ",".join(alleles), vaddr))
        return pd.DataFrame.from_records(records, columns=["id", "rsid", "chrom", "pos",
                                                           "nalleles", "allele_ids", "vaddr"])

    def read(self, index=None, dtype=np.float64):
        """Decompress only the requested variants, retaining sample order and duplicates."""
        sample_index, variant_index = (None, None) if index is None else index
        if sample_index is None:
            sample_index = slice(None)
        if variant_index is None:
            variant_index = slice(None)
        if isinstance(sample_index, slice):
            sample_count = len(range(*sample_index.indices(len(self.samples))))
        else:
            sample_index = _array_index(sample_index, len(self.samples), 'sample')
            sample_count = len(sample_index)
        nvariants = len(self._reader)
        if isinstance(variant_index, slice):
            variants = range(*variant_index.indices(nvariants))
        else:
            variants = _array_index(variant_index, nvariants, 'variant')
        if not nvariants:
            raise ValueError("BGEN file has no variants")
        first = int(variants[0]) if len(variants) else 0
        probs = self._reader[first].probabilities
        result = np.empty((sample_count, len(variants), probs.shape[1]), dtype=dtype, order="F")
        previous = first
        for column, i in enumerate(variants):
            if i != previous:
                probs = self._reader[int(i)].probabilities
                previous = i
            if probs.shape[1] != result.shape[2]:
                raise ValueError("Selected BGEN variants have different probability dimensions")
            result[:, column, :] = probs[sample_index]
        return result

    def read_dosage(self, samples=None, variants=None, *, order='F'):
        """Return first-allele dosage without retaining a batch of probability planes."""
        if samples is None:
            samples = slice(None)
        if isinstance(samples, slice):
            count = len(range(*samples.indices(len(self.samples))))
        else:
            samples = _array_index(samples, len(self.samples), 'sample')
            count = len(samples)
        if variants is None:
            variants = slice(None)
        if isinstance(variants, slice):
            variants = range(*variants.indices(len(self._reader)))
        else:
            variants = _array_index(variants, len(self._reader), 'variant')
        output = np.empty((count, len(variants)), np.float32, order=order)
        previous = None
        for column, variant in enumerate(variants):
            result = output[:, column]
            if variant == previous:
                result[:] = output[:, column - 1]
                continue
            probabilities = self._reader[int(variant)].probabilities
            if probabilities.shape[1] == 4:
                np.add(probabilities[samples, 0], probabilities[samples, 2], out=result, dtype=np.float32)
            elif probabilities.shape[1] == 3:
                np.multiply(probabilities[samples, 0], 2, out=result, dtype=np.float32)
                np.add(result, probabilities[samples, 1], out=result, dtype=np.float32)
            else:
                raise ValueError('Expected diploid, biallelic BGEN genotypes')
            previous = variant
        return output

    def close(self):
        self._reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()


def _array_index(index, size, axis):
    index = np.atleast_1d(index)
    if index.ndim != 1:
        raise IndexError(f'BGEN {axis} index must be one dimensional')
    if index.dtype.kind == 'b':
        if len(index) != size:
            raise IndexError(f'Boolean {axis} mask length differs from BGEN {axis} count')
        return np.flatnonzero(index)
    if not index.size:
        return index.astype(np.intp)
    if index.dtype.kind not in 'iu':
        raise IndexError(f'BGEN {axis} indices must be integers or booleans')
    if np.any(index >= size) or np.any(index < -size):
        raise IndexError(f'BGEN {axis} index out of bounds')
    return np.where(index < 0, index + size, index).astype(np.intp)
