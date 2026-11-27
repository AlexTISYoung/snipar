"""Reusable integer gathers for building a genotype design from variant blocks."""

from typing import NamedTuple

import numpy as np
from numba import njit, prange


class IdIndex:
    """Sorted ID lookup with compact row indices and last-occurrence semantics."""

    def __init__(self, ids):
        ids = np.asarray(ids, dtype=str)
        self.rows = None
        if ids.size > 1 and np.any(ids[:-1] > ids[1:]):
            rows = np.argsort(ids, kind='stable')
            self.keys = ids[rows]
            self.rows = rows.astype(_index_dtype(len(ids)), copy=False)
        else:
            self.keys = ids

    def lookup(self, ids):
        ids = np.asarray(ids, dtype=str)
        result = np.full(ids.shape, -1, dtype=_index_dtype(self.keys.size))
        if self.keys.size:
            pos = np.searchsorted(self.keys, ids, side='right') - 1
            np.maximum(pos, 0, out=pos)
            same = self.keys[pos] == ids
            result[same] = pos[same] if self.rows is None else self.rows[pos[same]]
        return result


def _index_dtype(size):
    return np.int32 if size <= np.iinfo(np.int32).max else np.int64


class SamplePlan(NamedTuple):
    ids: np.ndarray
    fams: np.ndarray
    observed: np.ndarray
    imputed: np.ndarray
    genotype_rows: np.ndarray
    parent_status: np.ndarray
    sibling_members: np.ndarray
    member_offsets: np.ndarray
    sibling_targets: np.ndarray
    target_offsets: np.ndarray
    trio_count: int
    parsum: bool
    sib: bool
    sib_diff: bool
    trios_sibs: bool
    file_n_observed: int

    @property
    def n_components(self):
        if self.sib_diff:
            return 2
        if self.trios_sibs:
            return 3
        return 3 + int(self.sib) - int(self.parsum)


def _siblings(pedigree, observed_ids, pedigree_index):
    """Assign eligible observed relatives to family groups, using the pedigree's sibling rules."""
    rows = pedigree_index.lookup(observed_ids)
    genotyped = np.flatnonzero(rows >= 0)
    rows = rows[genotyped]
    _, first, codes, sizes = np.unique(pedigree[rows, 0], return_index=True,
                                       return_inverse=True, return_counts=True)
    viable = sizes > 1
    anchor = rows[first[codes]]
    for parent in (2, 3):
        conflict = pedigree[rows, parent] != pedigree[anchor, parent]
        viable[codes[conflict]] = False
    groups = np.full(len(observed_ids), -1, dtype=_index_dtype(len(viable)))
    keep = viable[codes]
    groups[genotyped[keep]] = codes[keep]
    return groups, len(viable)


def _remap_rows(rows, sample_count):
    index = np.empty(sample_count, dtype=_index_dtype(sample_count))
    index[rows] = np.arange(len(rows), dtype=index.dtype)
    return index


def _sample_plan(ped, observed_ids, observed_index, imp_fams, ids, *, parsum=False, sib=False,
                 sib_diff=False, include_unrel=False, trios_sibs=False, verbose=False, select_only=False):
    if ped.ndim != 2 or ped.shape[1] < 4:
        raise ValueError('Pedigree must include family, individual, paternal and maternal IDs')
    ids = observed_ids if ids is None else np.asarray(ids)
    row_dtype = _index_dtype(max(len(observed_ids), 0 if imp_fams is None else len(imp_fams)))
    rows = np.full((len(ids), 3), -1, dtype=row_dtype)
    statuses = np.full((len(ids), 2), -1, dtype=np.int8)
    rows[:, 0] = observed_index.lookup(ids)
    pedigree_index = IdIndex(ped[:, 1])
    ped_rows = pedigree_index.lookup(ids)
    in_ped = np.flatnonzero(ped_rows >= 0)
    ped_rows = ped_rows[in_ped]
    families = np.full(len(ids), '', dtype=ped[:, 0].dtype)
    families[in_ped] = ped[ped_rows, 0]
    if not sib_diff:
        family_indices = None
        if imp_fams is not None:
            family_indices = IdIndex(imp_fams).lookup(families[in_ped])
        for parent in (0, 1):
            observed_parent = observed_index.lookup(ped[ped_rows, 2 + parent])
            present = observed_parent >= 0
            rows[in_ped[present], 1 + parent] = observed_parent[present]
            statuses[in_ped[present], parent] = 0
            if family_indices is not None:
                missing = (~present) & (family_indices >= 0)
                if ped.shape[1] >= 6:
                    missing &= ped[ped_rows, 4 + parent] == 'False'
                elif np.any(missing):
                    raise ValueError('Imputed genotypes require pedigree columns has_father and has_mother')
                rows[in_ped[missing], 1 + parent] = family_indices[missing]
                statuses[in_ped[missing], parent] = 1

    if sib or sib_diff or trios_sibs:
        family_groups, n_families = _siblings(ped, observed_ids, pedigree_index)
        with_sibling = np.zeros(len(ids), dtype=bool)
        called = rows[:, 0] >= 0
        with_sibling[called] = family_groups[rows[called, 0]] >= 0
    if sib_diff:
        cohort = np.flatnonzero(with_sibling)
        trio_count = 0
    elif trios_sibs:
        trios = np.flatnonzero((rows >= 0).all(axis=1))
        siblings = np.flatnonzero(with_sibling & (statuses < 0).any(axis=1))
        trio_count = len(trios)
        cohort = np.concatenate((trios, siblings))
        if verbose:
            print(f'{trio_count} individuals have complete parental genotypes.')
            print(f'{len(siblings)} individuals have no complete parental genotypes, but have sib genotypes.')
    else:
        complete = (rows >= 0).all(axis=1)
        if sib:
            complete &= with_sibling
        if not np.any(complete):
            detail = ' The unified estimator cannot be performed due to collinearity.' if include_unrel else ''
            raise ValueError('No individuals with phenotype observations and complete observed/imputed '
                             f'genotype observations.{detail}')
        eligible = complete
        if include_unrel:
            # The linear unrelated estimator predicts both parental dosages.
            eligible = eligible | ((rows[:, 0] >= 0) & (statuses < 0).all(axis=1))
        if sib:
            eligible &= with_sibling
        cohort = np.flatnonzero(eligible)
        trio_count = 0
        if not parsum:
            n_observed_parents = np.count_nonzero((statuses[cohort] == 0).any(axis=1))
            if not n_observed_parents:
                if verbose:
                    print('No individuals with observed parental genotypes. Using the sum of parental genotypes.')
                parsum = True
            elif n_observed_parents < 100 and verbose:
                print('Warning: few individuals have observed parental genotypes. Consider --parsum.')
    if not len(cohort):
        raise ValueError('No individuals with phenotype observations and genotype observations for this family design')

    if select_only:
        return ids[cohort], families[cohort], statuses[cohort], trio_count

    ids, families, rows, statuses = ids[cohort], families[cohort], rows[cohort], statuses[cohort]
    observed_needed = np.zeros(len(observed_ids), dtype=bool)
    observed_needed[rows[:, 0]] = True
    observed_needed[rows[:, 1:][statuses == 0]] = True
    empty = np.empty(0, dtype=row_dtype)
    sib_members, member_offsets, sib_targets, target_offsets = empty, empty, empty, empty
    if sib or sib_diff or (trios_sibs and trio_count < len(ids)):
        target_start = trio_count if trios_sibs else 0
        target_groups = family_groups[rows[target_start:, 0]]
        selected_groups = np.zeros(n_families + 1, dtype=bool)
        selected_groups[target_groups] = True
        members = np.flatnonzero(selected_groups[family_groups])
        observed_needed[members] = True
        order = np.argsort(family_groups[members], kind='stable')
        sib_members = members[order]
        member_sizes = np.bincount(family_groups[sib_members], minlength=n_families)[selected_groups[:-1]]
        member_offsets = np.insert(member_sizes.cumsum(), 0, 0).astype(row_dtype, copy=False)
        sib_targets = (target_start + np.argsort(target_groups, kind='stable')).astype(row_dtype, copy=False)
        target_sizes = np.bincount(target_groups, minlength=n_families)[selected_groups[:-1]]
        target_offsets = np.insert(target_sizes.cumsum(), 0, 0).astype(row_dtype, copy=False)

    observed = np.flatnonzero(observed_needed).astype(row_dtype, copy=False)
    obs_index = _remap_rows(observed, len(observed_needed))
    rows[:, 0] = obs_index[rows[:, 0]]
    sib_members = obs_index[sib_members]
    observed_parent = statuses == 0
    rows[:, 1:][observed_parent] = obs_index[rows[:, 1:][observed_parent]]
    imputed_parent = statuses == 1
    imputed = empty
    if np.any(imputed_parent):
        imputed_needed = np.zeros(len(imp_fams), dtype=bool)
        imputed_needed[rows[:, 1:][imputed_parent]] = True
        imputed = np.flatnonzero(imputed_needed).astype(row_dtype, copy=False)
        imp_index = _remap_rows(imputed, len(imp_fams))
        rows[:, 1:][imputed_parent] = imp_index[rows[:, 1:][imputed_parent]]
    plan = SamplePlan(ids, families, observed, imputed, rows, statuses, sib_members, member_offsets,
                      sib_targets, target_offsets, trio_count, parsum, sib, sib_diff, trios_sibs,
                      len(observed_ids))
    for value in plan:
        if isinstance(value, np.ndarray):
            value.setflags(write=False)
    return plan


@njit(inline='always')
def _observed_dosage(observed, row, variant):
    value = observed[row, variant]
    return np.float32(np.nan) if value == -127 else np.float32(value)


@njit(inline='always')
def _parent_dosage(observed, imputed, index, status, variant, flipped):
    if status == 0:
        return _observed_dosage(observed, index, variant)
    if status == 1:
        value = imputed[index, variant]
        return np.float32(2) - value if flipped else value
    return np.float32(-1)


@njit(inline='always')
def _family_row(out, row, observed, imputed, rows, status, flipped, parsum, sib, sib_diff, robust, counts, num_al):
    for j in range(out.shape[2]):
        out[row, 0, j] = _observed_dosage(observed, rows[row, 0], j)
    if sib_diff:
        return
    first_parent = 2 if sib else 1
    dad_status, mom_status = status[row, 0], status[row, 1]
    dad_index, mom_index = rows[row, 1], rows[row, 2]
    for j in range(out.shape[2]):
        dad = _parent_dosage(observed, imputed, dad_index, dad_status, j, flipped[j])
        mom = _parent_dosage(observed, imputed, mom_index, mom_status, j, flipped[j])
        if parsum:
            out[row, first_parent, j] = dad + mom if mom_status >= 0 else dad
        else:
            out[row, first_parent, j] = dad
            out[row, first_parent + 1, j] = mom
        if robust:
            count = np.float32(4) if dad_status == 0 and mom_status == 0 else np.float32(2)
            if dad_status == 1:
                count = num_al[dad_index, j]
            if mom_status == 1:
                count = num_al[mom_index, j]
            counts[row, j] = count


@njit(inline='always')
def _sibship(out, observed, members, member_offsets, targets, target_offsets, sib, trios_sibs, group):
    start, stop = member_offsets[group], member_offsets[group + 1]
    target_start, target_stop = target_offsets[group], target_offsets[group + 1]
    count = stop - start - int(sib)
    for j in range(out.shape[2]):
        total = np.float32(0)
        for i in range(start, stop):
            total += _observed_dosage(observed, members[i], j)
        for index in range(target_start, target_stop):
            row = targets[index]
            mean = (total - out[row, 0, j]) / count if sib else total / count
            out[row, 1, j] = mean
            if trios_sibs:
                out[row, 2, j] = mean


@njit(cache=True, nogil=True)
def _assemble(out, observed, imputed, rows, status, members, member_offsets, targets, target_offsets,
              parsum, sib, sib_diff, trios_sibs, flipped, robust, counts, num_al):
    for row in range(out.shape[0]):
        _family_row(out, row, observed, imputed, rows, status, flipped, parsum, sib, sib_diff,
                    robust, counts, num_al)
    for group in range(len(member_offsets) - 1):
        _sibship(out, observed, members, member_offsets, targets, target_offsets, sib, trios_sibs, group)


@njit(cache=True, nogil=True, parallel=True)
def _assemble_parallel(out, observed, imputed, rows, status, members, member_offsets, targets, target_offsets,
                       parsum, sib, sib_diff, trios_sibs, flipped, robust, counts, num_al):
    for row in prange(out.shape[0]):
        _family_row(out, row, observed, imputed, rows, status, flipped, parsum, sib, sib_diff,
                    robust, counts, num_al)
    for group in prange(max(0, len(member_offsets) - 1)):
        _sibship(out, observed, members, member_offsets, targets, target_offsets, sib, trios_sibs, group)


def _read_imputed(dataset, rows, variants):
    """Read a row/variant gather with bounded HDF5 transfer buffers."""
    result = np.empty((len(rows), len(variants)), np.float32)
    if not result.size:
        return result
    if not np.all(variants[1:] > variants[:-1]):
        unique, restore = np.unique(variants, return_inverse=True)
        return np.ascontiguousarray(_read_imputed(dataset, rows, unique)[:, restore])
    first, last = int(variants[0]), int(variants[-1]) + 1
    contiguous = last - first == len(variants)
    selection = slice(first, last) if contiguous else variants
    if rows[-1] - rows[0] + 1 == len(rows):
        dataset.read_direct(result, source_sel=np.s_[rows[0]:rows[-1] + 1, selection])
        return result
    if contiguous and len(rows) <= 64:
        dataset.read_direct(result, source_sel=np.s_[rows, selection])
        return result
    # HDF5 permits only one fancy index. Read dense row tiles and gather
    # within each tile when selected families are interleaved with others.
    tile_rows = max(1, 524_288 // len(variants))
    start = 0
    while start < len(rows):
        lo = int(rows[start])
        stop = np.searchsorted(rows, lo + tile_rows)
        hi = int(rows[stop - 1]) + 1
        source = np.s_[lo:hi, selection]
        if hi - lo == stop - start:
            dataset.read_direct(result, source_sel=source, dest_sel=np.s_[start:stop, :])
        else:
            block = np.empty((hi - lo, len(variants)), np.float32)
            dataset.read_direct(block, source_sel=source)
            result[start:stop] = block[rows[start:stop] - lo]
        start = stop
    return result
