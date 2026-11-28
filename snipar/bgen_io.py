"""Array access to BGEN files on interpreters without a binary cbgen wheel."""

from contextlib import contextmanager
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd
from bgen import BgenReader, __version__ as _bgen_version

from snipar.bgen_packed import PackedBgenSubset, subset_samples
from snipar.utilities import read_oxford_sample


class BgenFile:
    """Expose the variant metadata and sample × variant × probability slices used by SNIPar."""

    def __init__(self, filename):
        self._filename = Path(filename).absolute()
        self._reader = BgenReader(str(self._filename))
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
        return self.read_allele_ids(range(len(self._reader)))

    def read_allele_ids(self, variants):
        return np.asarray([",".join(self._reader[int(i)].alleles) for i in variants])

    def variant_metadata(self, variants=None):
        """Read the identifying fields without decompressing genotype probabilities."""
        records = []
        if variants is None:
            variants = range(len(self._reader))
        for i in variants:
            variant = self._reader[int(i)]
            alleles = variant.alleles
            records.append((variant.varid, variant.rsid, variant.chrom, variant.pos,
                            len(alleles), ",".join(alleles), _probability_offset(variant, alleles)))
        return pd.DataFrame.from_records(records, columns=["id", "rsid", "chrom", "pos",
                                                           "nalleles", "allele_ids", "vaddr"])

    def is_phased(self, index):
        return self._reader[int(index)].is_phased

    @contextmanager
    def _packed_subset(self, samples, sample_count):
        count = len(self.samples)
        selected = subset_samples(samples, count, sample_count, minimum_count=32_768)
        # The native reader's zlib implementation outperforms gathering through
        # Python zlib. Native decoding also wins for smaller sample counts.
        if selected is None or self._reader.header.compression == 'zlib':
            yield None
            return
        # Require a known native implementation when reusing its rounded LUT.
        with self._filename.open('rb') as stream:
            yield PackedBgenSubset.from_stream(stream, count, selected, native=True,
                                              native_table=_bgen_version == '1.11.0')

    def read(self, index=None, dtype=np.float64, *, expected_phased=None):
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
        def probabilities(index, packed):
            variant = self._reader[int(index)]
            decoded = _read_packed(variant, packed)
            if decoded is not None:
                values, phased = decoded
                if expected_phased is not None and phased != expected_phased:
                    raise ValueError('BGEN hard calls require a consistently phased or unphased file')
                return values
            if expected_phased is not None:
                if variant.is_phased != expected_phased:
                    raise ValueError('BGEN hard calls require a consistently phased or unphased file')
                if len(variant.alleles) != 2 or np.any(variant.ploidy[sample_index] != 2):
                    raise ValueError('Expected diploid, biallelic BGEN genotypes')
                return variant.probabilities[sample_index, :4 if expected_phased else 3]
            return variant.probabilities[sample_index]

        with self._packed_subset(sample_index, sample_count) as packed:
            probs = probabilities(first, packed)
            result = np.empty((sample_count, len(variants), probs.shape[1]), dtype=dtype, order="F")
            previous = first
            for column, i in enumerate(variants):
                if i != previous:
                    probs = probabilities(i, packed)
                    previous = i
                if probs.shape[1] != result.shape[2]:
                    raise ValueError("Selected BGEN variants have different probability dimensions")
                result[:, column, :] = probs
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
        with self._packed_subset(samples, count) as packed:
            previous = None
            for column, variant in enumerate(variants):
                result = output[:, column]
                if variant == previous:
                    result[:] = output[:, column - 1]
                    continue
                record = self._reader[int(variant)]
                decoded = _read_packed(record, packed)
                selected = samples
                if decoded is None:
                    if len(record.alleles) != 2 or np.any(record.ploidy[samples] != 2):
                        raise ValueError('Expected diploid, biallelic BGEN genotypes')
                    probabilities, phased = record.probabilities, record.is_phased
                else:
                    probabilities, phased = decoded
                    selected = slice(None)
                if phased:
                    np.add(probabilities[selected, 0], probabilities[selected, 2], out=result, dtype=np.float32)
                else:
                    np.multiply(probabilities[selected, 0], 2, out=result, dtype=np.float32)
                    np.add(result, probabilities[selected, 1], out=result, dtype=np.float32)
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


def _probability_offset(variant, alleles):
    # The genotype block follows the variable-length BGEN identifier fields.
    fields = (variant.varid, variant.rsid, variant.chrom, *alleles)
    header_size = (12 if variant.layout == 2 else 14) + 4 * len(alleles)
    return variant.fileoffset + header_size + sum(len(field.encode()) for field in fields)


def _read_packed(record, packed):
    if packed is not None:
        alleles = record.alleles
        if len(alleles) == 2:
            return packed.read(_probability_offset(record, alleles))
    return None
