"""Keep a chromosome-level BGEN variant selection without copying genotypes."""

from contextlib import ExitStack
from dataclasses import dataclass
from functools import cached_property
from itertools import repeat
import json
from pathlib import Path

import numpy as np


BGEN_SUBSET_SUFFIX = '.snipar-subset.bgen'
_VARIANT_FIELDS = frozenset(('ids', 'rsids', 'chromosomes', 'positions', 'allele_ids',
                             'nalleles', 'phased', '_vaddr'))


@dataclass(frozen=True)
class BGENSubset:
    path: Path
    n_variants: int
    n_kept: int


def _selected_rows(ids, aliases, passing_ids, known_ids, variant_ids):
    known_ids = passing_ids if known_ids is None else known_ids
    # Hits start unclaimed, then record a primary match (1) or an RSID match (2).
    hits = dict.fromkeys(passing_ids, 0)

    def matches():
        for identifier, rsid in zip(ids, repeat(None) if aliases is None else aliases):
            primary = identifier in known_ids
            key = identifier if primary else rsid
            previous = hits.get(key)
            if previous is None:
                yield False
                continue
            if previous:
                raise ValueError(f'BGEN quality ID {key!r} matches multiple variants')
            hits[key] = 1 if primary else 2
            yield variant_ids is None or identifier in variant_ids or rsid in variant_ids

    rows = np.flatnonzero(np.fromiter(matches(), dtype=bool, count=len(ids)))
    if aliases is not None and 2 in hits.values():
        # An RSID used as a fallback must identify one record even if another
        # record with that RSID already has its own primary-ID score.
        for rsid in aliases:
            previous = hits.get(rsid)
            if previous == 2:
                hits[rsid] = 3
            elif previous == 3:
                raise ValueError(f'BGEN quality RSID {rsid!r} matches multiple variants')
    return rows


def filter_bgen(source, target, passing_ids, *, variant_ids=None, known_ids=None):
    """Write an internal BGEN selection that can be opened by analysis workers.

    ``passing_ids`` and ``known_ids`` are sets of identifiers from the quality
    table. Use the RSID when the primary BGEN ID is absent from that table.
    The descriptor references the original BGEN and its sample information.
    Selection uses variant identifiers from the source reader, and copies no
    genotype blocks.
    """
    from snipar.read.variants import bgen_variant_ids
    from snipar.utilities import open_bgen

    source = Path(source).absolute()
    target = Path(target)
    if not target.name.endswith(BGEN_SUBSET_SUFFIX):
        raise ValueError(f'BGEN selection filename must end with {BGEN_SUBSET_SUFFIX}')
    if source.name.endswith(BGEN_SUBSET_SUFFIX):
        raise ValueError('Supply the original BGEN when creating a quality selection')
    original = source.stat()
    with open_bgen(source) as reader:
        ids, aliases = bgen_variant_ids(reader)
        n_variants = len(ids)
        rows = _selected_rows(ids, aliases, passing_ids, known_ids, variant_ids)
    n_kept = len(rows)
    if not n_kept:
        raise ValueError(f'No variants passed the BGEN quality filters in {source}')

    metadata = {'source': str(source), 'size': original.st_size, 'mtime_ns': original.st_mtime_ns,
                'n_variants': n_variants, 'n_kept': n_kept}
    with ExitStack() as rollback:
        if n_kept != n_variants:
            indices = target.with_suffix('.npy')
            with indices.open('xb') as stream:
                rollback.callback(indices.unlink)
                np.save(stream, rows.astype(np.uint32, copy=False), allow_pickle=False)
        with target.open('x', encoding='utf-8') as descriptor:
            rollback.callback(target.unlink)
            json.dump(metadata, descriptor)
        rollback.pop_all()
    return BGENSubset(target, n_variants, n_kept)


def open_bgen_subset(filename, open_source, verbose=False):
    """Open a selection created alongside the current analysis's temporary files."""
    filename = Path(filename)
    with filename.open(encoding='utf-8') as stream:
        metadata = json.load(stream)
    source = Path(metadata['source'])
    current = source.stat()
    if (current.st_size, current.st_mtime_ns) != (metadata['size'], metadata['mtime_ns']):
        raise ValueError(f'BGEN source changed after quality filtering: {source}')
    reader = open_source(source, verbose=verbose)
    if metadata['n_kept'] == metadata['n_variants']:
        return reader
    try:
        indices = np.load(filename.with_suffix('.npy'), mmap_mode='r', allow_pickle=False)
        cls = _NativeBGENSubset if hasattr(reader, 'read_dosage') else _BGENSubsetReader
        return cls(reader, indices)
    except BaseException:
        reader.close()
        raise


class _BGENSubsetReader:
    def __init__(self, reader, indices):
        self._reader = reader
        self._indices = indices
        self.nvariants = len(indices)

    def __getattr__(self, name):
        value = getattr(self._reader, name)
        if name in _VARIANT_FIELDS:
            value = value[self._indices]
            setattr(self, name, value)
        return value

    @cached_property
    def analysis_variant_ids(self):
        from snipar.read.variants import bgen_variant_ids

        _, aliases = bgen_variant_ids(self._reader)
        return (self.rsids, None) if aliases is None else (self.ids, self.rsids)

    @property
    def shape(self):
        nsamples, _, *states = self._reader.shape
        return (nsamples, self.nvariants, *states)

    def _variants(self, variants):
        return self._indices if variants is None else self._indices[variants]

    def read(self, index=None, dtype=np.float64, **kwargs):
        samples, variants = (None, None) if index is None else index
        return self._reader.read((samples, self._variants(variants)), dtype=dtype, **kwargs)

    def close(self):
        try:
            self._reader.close()
        finally:
            self._indices._mmap.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()


class _NativeBGENSubset(_BGENSubsetReader):
    @cached_property
    def allele_ids(self):
        return self._reader.read_allele_ids(self._indices)

    def is_phased(self, variant):
        return self._reader.is_phased(self._indices[variant])

    def read_dosage(self, samples=None, variants=None, *, order='F'):
        return self._reader.read_dosage(samples, self._variants(variants), order=order)

    def variant_metadata(self):
        return self._reader.variant_metadata(self._indices)
