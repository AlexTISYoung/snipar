"""Export parent-of-origin haplotypes from chromosome-local duoHMM calls."""

from contextlib import ExitStack
import gzip
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import NamedTuple

import h5py
from numba import njit
import numpy as np

from snipar.imputation.origin_io import _chromosome_label
from snipar.read.family import IdIndex
from snipar.resources import plan_batches
from snipar.utilities import convert_str_array, encode_str_array
from snipar.vcf import _PhasedWriter


_DUO_TILE = 256
_HEADER_DATA = b'SNIPar haplotype order: paternal|maternal'


class _OriginBatch(NamedTuple):
    path: Path
    children: np.ndarray
    maternal: np.ndarray


def _read_batch(path, chromosome, pedigree, children, grid):
    with h5py.File(path, 'r') as file:
        keys = {'positions', 'variant_ids', 'alleles', 'pairs', 'origin', 'phase_evidence'}
        if file.attrs.get('format_version') != 1 or not keys.issubset(file):
            raise ValueError(f'{path} lacks parent-of-origin labels and marker allele metadata')
        if _chromosome_label(file.attrs.get('chromosome')) != chromosome:
            raise ValueError(f'Parent-of-origin sidecar belongs to a different chromosome: {path}')
        positions, ids, alleles = (file[key][:] for key in ('positions', 'variant_ids', 'alleles'))
        if (positions.ndim != 1 or not len(positions) or ids.shape != positions.shape
                or alleles.shape != (len(positions), 2)):
            raise ValueError(f'Invalid parent-of-origin marker metadata in {path}')
        if grid is None:
            grid = (positions, ids, alleles)
        elif not all(np.array_equal(left, right) for left, right in zip(grid, (positions, ids, alleles))):
            raise ValueError(f'Parent-of-origin sidecars have different marker grids: {path}')

        pairs = file['pairs']
        if pairs.ndim != 2 or pairs.shape[1] != 2 or not len(pairs):
            raise ValueError(f'Invalid parent-of-origin duo metadata in {path}')
        expected = (len(pairs), len(positions))
        if (file['origin'].shape != expected or file['origin'].dtype != np.int8
                or file['phase_evidence'].shape != expected or file['phase_evidence'].dtype != np.uint8):
            raise ValueError(f'Invalid parent-of-origin duo metadata in {path}')
        offspring = np.empty(len(pairs), dtype=np.intp)
        maternal = np.empty(len(pairs), dtype=bool)
        for first in range(0, len(pairs), 8192):
            rows = slice(first, min(first + 8192, len(pairs)))
            duo_ids = convert_str_array(pairs[rows])
            child_rows = children.lookup(duo_ids[:, 1])
            parents = duo_ids[:, :1] == pedigree[child_rows, 2:4]
            if (np.any(child_rows < 0) or np.any(duo_ids[:, 0] == '0')
                    or np.any(parents[:, 0] == parents[:, 1])):
                raise ValueError(f'Parent-of-origin sidecar has a duo without a matching parental role in {path}')
            offspring[rows], maternal[rows] = child_rows, parents[:, 1]
    return _OriginBatch(path, offspring, maternal), grid


@njit(cache=True)
def _merge_origin_tile(evidence, origin, children, maternal, calls, labels):
    for marker in range(origin.shape[1]):
        for row in range(origin.shape[0]):
            value = evidence[row, marker]
            if value == 255:
                continue
            child = children[row]
            alleles = value >> 4
            previous = calls[marker, child]
            if previous != 255 and previous != alleles:
                raise ValueError('Parent-of-origin sidecars disagree on original offspring haplotypes')
            paternal = origin[row, marker]
            if (alleles & 3) > 2 or (alleles >> 2) > 2 or not -1 <= paternal <= 1:
                raise ValueError('Invalid allele or orientation in parent-of-origin sidecar')
            calls[marker, child] = alleles
            if paternal >= 0:
                labels[marker, child] |= 1 << (paternal ^ maternal[row])


@njit(cache=True)
def _order_alleles(calls, labels, paternal, maternal):
    for i in range(len(calls)):
        first, second = calls[i] & 3, calls[i] >> 2
        paternal[i], maternal[i] = -1, -1
        if first > 1 or second > 1:
            continue
        # Sidecar hard calls count the first allele; BGEN and VCF writers
        # expect allele indices. Homozygotes have an unambiguous order.
        if first == second or labels[i] == 1:
            paternal[i], maternal[i] = 1 - first, 1 - second
        elif labels[i] == 2:
            paternal[i], maternal[i] = 1 - second, 1 - first


@njit(cache=True)
def _pack_vcf(paternal, maternal, buffer):
    for sample in range(len(paternal)):
        j = 4 * sample
        buffer[j + 1] = paternal[sample] + 48 if paternal[sample] >= 0 else 46
        buffer[j + 3] = maternal[sample] + 48 if maternal[sample] >= 0 else 46


class _VcfWriter:
    def __init__(self, stack, path, samples, chromosome):
        if path.suffix == '.gz':
            self.file = stack.enter_context(gzip.open(path, 'xb', compresslevel=3))
        else:
            self.file = stack.enter_context(path.open('xb'))
        self.file.write(b'##fileformat=VCFv4.2\n##source=SNIPar\n'
                        b'##SNIPar_haplotype_order=paternal|maternal\n'
                        b'##contig=<ID=' + chromosome + b'>\n'
                        b'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype ordered paternal|maternal">\n'
                        b'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
        self.file.write(b'\t'.join(samples) + b'\n')
        self.buffer = np.empty(4 * len(samples), dtype=np.uint8)
        self.buffer[::4], self.buffer[2::4] = 9, 124

    def write(self, chromosome, pos, marker, ref, alt, paternal, maternal):
        _pack_vcf(paternal, maternal, self.buffer)
        self.file.write(b'\t'.join((chromosome, str(pos).encode(), marker, ref, alt, b'.', b'.', b'.', b'GT')))
        self.file.write(self.buffer)
        self.file.write(b'\n')


def _write_haplotypes(writer, chromosome, samples, batches, positions, ids, alleles):
    n_samples, n_markers = len(samples), len(positions)
    bytes_per_marker = 2 * (n_samples + min(_DUO_TILE, max(len(batch.children) for batch in batches)))
    cap = max(1, min(2048, (512 * 1024 ** 2) // bytes_per_marker))
    _, width = plan_batches(n_markers, bytes_per_marker, batch_size=cap, fixed_bytes=12 * n_samples)
    if width >= 256 and width < n_markers:
        width -= width % 256
    paternal = np.empty(n_samples, dtype=np.int8)
    maternal = np.empty_like(paternal)
    for first in range(0, n_markers, width):
        stop = min(first + width, n_markers)
        columns = slice(first, stop)
        calls = np.full((stop - first, n_samples), 255, dtype=np.uint8)
        labels = np.zeros_like(calls)
        for batch in batches:
            with h5py.File(batch.path, 'r') as file:
                for row in range(0, len(batch.children), _DUO_TILE):
                    rows = slice(row, row + _DUO_TILE)
                    _merge_origin_tile(file['phase_evidence'][rows, columns], file['origin'][rows, columns],
                                       batch.children[rows], batch.maternal[rows], calls, labels)
        for j, marker in enumerate(range(first, stop)):
            _order_alleles(calls[j], labels[j], paternal, maternal)
            writer.write(chromosome, positions[marker], ids[marker], *alleles[marker], paternal, maternal)


def export_parental_haplotypes(path, chromosome, pedigree, origin_paths):
    """Write a chromosome of paternal-first, maternal-second hard calls.

    The files include offspring with at least one genotyped parent processed by
    duoHMM. Homozygotes retain both called alleles. Unresolved or conflicting
    heterozygotes and incomplete genotypes are written as fully missing.
    Sidecars from multiple sibling batches must use the same marker grid.
    """
    path = Path(path)
    bgen = path.suffix == '.bgen'
    if not bgen and not path.name.endswith(('.vcf', '.vcf.gz')):
        raise ValueError('Parent-of-origin haplotype output must end in .bgen, .vcf, or .vcf.gz')
    if isinstance(origin_paths, (str, Path)):
        origin_paths = (origin_paths,)
    chromosome = _chromosome_label(chromosome)
    pedigree = convert_str_array(pedigree)
    if pedigree.ndim != 2 or pedigree.shape[1] < 4 or not len(pedigree):
        raise ValueError('Parent-of-origin export requires a four-column pedigree')
    children = IdIndex(pedigree[:, 1])
    keep = np.zeros(len(pedigree), dtype=bool)
    batches, grid = [], None
    for origin_path in origin_paths:
        batch, grid = _read_batch(Path(origin_path), chromosome, pedigree, children, grid)
        keep[batch.children] = True
        batches.append(batch)
    if not batches:
        raise ValueError(f'No parent-of-origin calls are available to export for chromosome {chromosome}')
    samples = encode_str_array(pedigree[keep, 1])
    if any(not sample or any(char in sample for char in (b' ', b'\t', b'\n', b'\r')) for sample in samples):
        raise ValueError('Parent-of-origin export requires nonempty sample IDs without whitespace')
    destination = np.flatnonzero(keep)
    mapping = np.empty(len(keep), dtype=np.intp)
    mapping[destination] = np.arange(len(destination), dtype=np.intp)
    for batch in batches:
        batch.children[:] = mapping[batch.children]
    del destination, mapping, keep, children, pedigree

    chrom = chromosome.encode()
    with TemporaryDirectory(prefix='.snipar-haplotypes-', dir=path.parent) as directory:
        staged = Path(directory) / path.name
        with ExitStack() as stack:
            writer = (_PhasedWriter(stack, staged, samples, header_data=_HEADER_DATA) if bgen else
                      _VcfWriter(stack, staged, samples, chrom))
            _write_haplotypes(writer, chrom, samples, batches, *grid)
            if bgen:
                writer.finish()
        if bgen and staged.with_suffix('.sample').exists():
            staged.with_suffix('.sample').replace(path.with_suffix('.sample'))
        staged.replace(path)
    return path
