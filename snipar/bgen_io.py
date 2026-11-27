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
        if self._reader.header.has_sample_ids:
            self.samples = np.asarray(self._reader.samples, dtype=str)
            if self.samples.size and self.samples[0] != "sample_0":
                return
        sample_file = Path(filename).with_suffix(".sample")
        self.samples = read_oxford_sample(sample_file)
        if len(self.samples) != self._reader.header.nsamples:
            raise ValueError("Sample file length differs from the BGEN sample count")

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
            # Offset of the genotype block after the variable-length BGEN identifier fields.
            strings = [variant.varid, variant.rsid, variant.chrom, *variant.alleles]
            header_size = (12 if variant.layout == 2 else 14) + 4 * len(variant.alleles)
            vaddr = variant.fileoffset + header_size + sum(len(s.encode()) for s in strings)
            records.append((variant.varid, variant.rsid, variant.chrom, variant.pos,
                            len(variant.alleles), ",".join(variant.alleles), vaddr))
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
            sample_index = np.atleast_1d(sample_index)
            sample_count = (np.count_nonzero(sample_index) if sample_index.dtype == bool
                            else len(sample_index))
        if isinstance(variant_index, slice):
            variants = range(*variant_index.indices(len(self._reader)))
        else:
            variants = np.atleast_1d(variant_index)
            if variants.dtype == bool:
                if len(variants) != len(self._reader):
                    raise IndexError("Boolean variant mask length differs from BGEN variant count")
                variants = np.flatnonzero(variants)
        if not len(self._reader):
            raise ValueError("BGEN file has no variants")
        first = int(variants[0]) if len(variants) else 0
        probs = self._reader[first].probabilities
        result = np.empty((sample_count, len(variants), probs.shape[1]), dtype=dtype, order="F")
        for column, i in enumerate(variants):
            if column:
                probs = self._reader[int(i)].probabilities
            if probs.shape[1] != result.shape[2]:
                raise ValueError("Selected BGEN variants have different probability dimensions")
            result[:, column, :] = probs[sample_index]
        return result

    def close(self):
        self._reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()
