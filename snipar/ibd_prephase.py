"""Bounded phased input and duoHMM pre-phasing for sibling IBD inference."""

from concurrent.futures import ThreadPoolExecutor

import h5py
from numba import njit
import numpy as np

from snipar.config import nan_integer
from snipar.ibd_phased import infer_phased_pair_ibd
from snipar.imputation.duohmm import call_parental_origin, recombination_probabilities
from snipar.imputation.origin_io import OriginWriter
from snipar.read.bgen import read_diploid_probabilities
from snipar.read.family import IdIndex
from snipar.resources import plan_batches


def observed_parent_duos(pedigree, sample_ids, child_rows=None, *, sample_index=None):
    """Return directed sample rows and the parental role (father 0, mother 1)."""
    if pedigree is None:
        return np.empty((0, 2), dtype=np.intp), np.empty(0, dtype=np.intp)
    pedigree = np.asarray(pedigree, dtype=str)
    if child_rows is None:
        child_rows = np.arange(len(sample_ids), dtype=np.intp)
    children = IdIndex(pedigree[:, 1]).lookup(sample_ids[child_rows])
    found = children >= 0
    child_rows = child_rows[found]
    parent_ids = pedigree[children[found], 2:4]
    if sample_index is None:
        sample_index = IdIndex(sample_ids)
    parents = sample_index.lookup(parent_ids)
    parents[parent_ids == '0'] = -1
    valid = (parents >= 0) & (parents != child_rows[:, None])
    # A malformed pedigree naming the same parent twice cannot assign a role.
    valid[parents[:, 0] == parents[:, 1]] = False
    child, role = np.nonzero(valid)
    return np.column_stack((parents[child, role], child_rows[child])), role


def read_phased_calls(source, rows, variants):
    from snipar.imputation.impute_from_sibs import decode_phased_genotypes

    probabilities = read_diploid_probabilities(
        source.reader, rows, variants, True, num_threads=source.threads,
    )
    return decode_phased_genotypes(probabilities, len(rows), np.arange(len(rows), dtype=np.intp),
                                   np.zeros(len(variants), dtype=np.uint8))


@njit(cache=True, inline='always')
def _encode(allele):
    return allele if allele >= 0 else 2


@njit(cache=True)
def pack_phased(calls, packed):
    """Store two phased variants per byte, with independent missing alleles."""
    for row in range(calls.shape[0]):
        for j in range(0, calls.shape[1], 2):
            code = _encode(calls[row, j, 0]) | (_encode(calls[row, j, 1]) << 2)
            if j + 1 < calls.shape[1]:
                code |= (_encode(calls[row, j + 1, 0]) << 4) | (_encode(calls[row, j + 1, 1]) << 6)
            packed[row, j // 2] = code


@njit(cache=True, nogil=True)
def _unpack(packed, rows, markers, missing):
    calls = np.empty((len(rows), len(markers), 2), dtype=np.int8)
    for i in range(len(rows)):
        for j in range(len(markers)):
            marker = markers[j]
            code = packed[rows[i], marker // 2] >> (4 * (marker % 2))
            for allele in range(2):
                value = (code >> (2 * allele)) & 3
                calls[i, j, allele] = value if value < 2 else missing
    return calls


def origin_grid(reader, positions):
    """Select the same unique marker IDs used by BGEN imputation."""
    rsids = np.asarray(reader.rsids)
    ids = np.where(np.isin(rsids, ('', '.')), reader.ids, rsids)
    _, first = np.unique(ids, return_index=True)
    first.sort()
    markers = first[np.argsort(positions[first], kind='stable')]
    alleles = np.array([value.split(',') for value in np.asarray(reader.allele_ids)[markers]])
    return markers, ids[markers], alleles


def write_parent_origins(path, chrom, pairs, sample_ids, packed, markers, positions, ids, alleles,
                         genetic_map, threads):
    """Run forward/backward duoHMM on the full chromosome before IBD pruning."""
    recombination = recombination_probabilities(positions, genetic_map - min(0., genetic_map[0]))
    worker_count, _ = plan_batches(len(pairs), 16 * len(markers), workers=threads,
                                   fixed_bytes_per_worker=32 * len(markers))
    _, batch = plan_batches(len(pairs), 16 * len(markers), batch_size=256,
                            fixed_bytes=worker_count * 32 * len(markers))
    with OriginWriter(path, chrom, sample_ids[pairs], positions, ids, alleles) as writer:
        for start in range(0, len(pairs), batch):
            rows = slice(start, min(start + batch, len(pairs)))
            samples = np.unique(pairs[rows])
            local_pairs = np.searchsorted(samples, pairs[rows])
            phased = _unpack(packed, samples, markers, nan_integer)
            origin, copy, _ = call_parental_origin(phased, local_pairs, recombination, threads=worker_count)
            writer.write(rows, slice(None), origin, copy, phased, local_pairs)


class PhasedPairInference:
    """Reuse shared marker data and a bounded thread pool across IBD pair batches."""

    def __init__(self, stack, genotypes, markers, freqs, genetic_map, weights, errors, switch_rate,
                 n_pairs, threads, *, duos, duo_roles, origin_path, origin_markers, posterior):
        self.genotypes = genotypes
        self.markers = markers
        self.freqs = freqs
        self.genetic_map = genetic_map
        self.weights = weights
        self.errors = errors / 2
        self.switch_rate = switch_rate
        self.posterior = posterior
        self.origin = None
        per_pair = 24 * len(markers) + 1
        if origin_path is not None:
            self.origin = stack.enter_context(h5py.File(origin_path, 'r'))
            self.duos = duos
            self.duo_by_child = np.full((len(genotypes), 2), -1, dtype=np.intp)
            self.duo_by_child[duos[:, 1], duo_roles] = np.arange(len(duos), dtype=np.intp)
            source_to_origin = np.full(2 * genotypes.shape[1], -1, dtype=np.intp)
            source_to_origin[origin_markers] = np.arange(len(origin_markers), dtype=np.intp)
            origin_columns = source_to_origin[markers]
            self.has_origin_column = origin_columns >= 0
            self.origin_columns = origin_columns[self.has_origin_column]
            per_pair += 16 * (len(origin_markers) + len(markers))
        worker_bytes = 320 * len(markers) + 16 * 1024
        workers, _ = plan_batches(n_pairs, per_pair, workers=threads, fixed_bytes_per_worker=worker_bytes)
        cap = max(1, min(10000, (64 * 1024**2) // per_pair))
        _, self.batch_size = plan_batches(n_pairs, per_pair, batch_size=cap, fixed_bytes=workers * worker_bytes)
        self.executor = stack.enter_context(ThreadPoolExecutor(max_workers=workers)) if workers > 1 else None

    def _read_labels(self, pairs):
        if self.origin is None:
            return None, None, None
        duo_rows = self.duo_by_child[pairs]
        present = duo_rows >= 0
        if not present.any() or not self.has_origin_column.any():
            return None, None, None
        duos = np.unique(duo_rows[present])
        local = np.full_like(duo_rows, -1)
        local[present] = np.searchsorted(duos, duo_rows[present])
        labels = []
        for field in ('origin', 'parent_copy'):
            # The last row is an all-unknown sentinel addressable by local index -1.
            data = np.full((len(duos) + 1, len(self.markers)), -1, dtype=np.int8)
            data[:-1, self.has_origin_column] = self.origin[field][duos][:, self.origin_columns]
            labels.append(data)
        return local, *labels

    def _pair_labels(self, pair, local, origins, copies):
        father, mother = origins[local[:, 0]], origins[local[:, 1]]
        maternal_origin = np.where(mother >= 0, 1 - mother, -1)
        paternal_origin = np.where(father >= 0, father, maternal_origin)
        conflict = (father >= 0) & (maternal_origin >= 0) & (father != maternal_origin)
        paternal_origin[conflict] = -1
        labels = {'parent_origin': paternal_origin}
        pair_duos = self.duo_by_child[pair]
        for role, name in enumerate(('parent_haplotype', 'other_parent_haplotype')):
            left, right = pair_duos[:, role]
            if left >= 0 and right >= 0 and self.duos[left, 0] == self.duos[right, 0]:
                homologs = copies[local[:, role]]
                homologs[conflict] = -1
                labels[name] = homologs
        return labels

    def __call__(self, pairs):
        local, origins, copies = self._read_labels(pairs)
        nmarkers = len(self.markers)
        paths = np.empty((len(pairs), nmarkers), dtype=np.int8)
        posterior = np.empty((len(pairs), nmarkers, 3), dtype=np.float32) if self.posterior else None

        def infer(index):
            calls = _unpack(self.genotypes, pairs[index], self.markers, -1)
            following = calls[:, 1:]
            heterozygous = (following[:, :, 0] >= 0) & (following[:, :, 1] >= 0)
            heterozygous &= following[:, :, 0] != following[:, :, 1]
            labels = {}
            if local is not None and (local[index] >= 0).any():
                labels = self._pair_labels(pairs[index], local[index], origins, copies)
            paths[index], probabilities = infer_phased_pair_ibd(
                calls, self.freqs, self.genetic_map, self.switch_rate * heterozygous,
                self.errors, weights=self.weights, posterior=self.posterior, **labels,
            )
            if posterior is not None:
                posterior[index] = probabilities

        if self.executor is None:
            for index in range(len(pairs)):
                infer(index)
        else:
            for _ in self.executor.map(infer, range(len(pairs))):
                pass
        return paths, posterior
