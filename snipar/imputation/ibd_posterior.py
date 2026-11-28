"""Read sibling IBD posterior probabilities for the current genotype block."""

from itertools import combinations

import h5py
import numpy as np


def _pair_key(first, second):
    return (first, second) if first <= second else (second, first)


def _aligned_runs(indices):
    """Yield slices in a sorted subset and its source for contiguous runs."""
    if not len(indices):
        return []
    ends = np.r_[0, np.flatnonzero(np.diff(indices) != 1) + 1, len(indices)]
    return [(slice(first, stop), slice(indices[first], indices[stop - 1] + 1))
            for first, stop in zip(ends[:-1], ends[1:])]


class SiblingIBDPosteriors:
    """An indexed sidecar that only decodes selected pairs and marker columns."""

    def __init__(self, path, sibships, *, chromosome=None, include_larger_families=False):
        self.path = path
        self.chromosome = chromosome
        self.required = {
            _pair_key(first, second)
            for children in sibships["IID"]
            if len(children) == 2 or include_larger_families
            for first, second in combinations(children, 2)
        }

    def __enter__(self):
        self.file = h5py.File(self.path, "r")
        try:
            self._read_index()
        except Exception:
            self.file.close()
            raise
        return self

    def __exit__(self, *exc):
        self.file.close()

    def _read_index(self):
        if not {"pairs", "positions", "probabilities"}.issubset(self.file):
            raise ValueError("Sibling IBD posterior sidecar requires pairs, positions and probabilities")
        if self.file.attrs.get("format_version", 1) != 1:
            raise ValueError("Unsupported sibling IBD posterior format version")
        stored_chromosome = self.file.attrs.get("chromosome")
        if isinstance(stored_chromosome, bytes):
            stored_chromosome = stored_chromosome.decode()
        if (self.chromosome is not None and stored_chromosome is not None and
                str(stored_chromosome).removeprefix("chr").removesuffix(".0") != str(self.chromosome)):
            raise ValueError("Sibling IBD posterior chromosome differs from genotype data")
        pairs = self.file["pairs"]
        positions = self.file["positions"]
        self.probabilities = self.file["probabilities"]
        if (pairs.ndim != 2 or pairs.shape[1] != 2 or positions.ndim != 1 or positions.dtype.kind not in "iu"
                or self.probabilities.shape != (len(pairs), len(positions), 3)
                or self.probabilities.dtype.kind != "f"):
            raise ValueError("Invalid sibling IBD posterior dataset shapes or types")
        self.positions = positions[:].astype(np.int64, copy=False)
        if np.any(np.diff(self.positions) < 0):
            raise ValueError("Sibling IBD posterior positions must be sorted")
        self.variant_ids = None
        if "variant_ids" in self.file:
            identifiers = self.file["variant_ids"]
            if identifiers.shape != self.positions.shape:
                raise ValueError("Sibling IBD posterior variant IDs must match positions")
            self.variant_ids = identifiers[:]
        self.pair_index = {}
        rows = []
        if self.required:
            for first in range(0, len(pairs), 8192):
                for i, (left, right) in enumerate(pairs[first:first + 8192], first):
                    pair = _pair_key(left, right)
                    if pair in self.required:
                        if pair in self.pair_index:
                            raise ValueError("Duplicate sibling pair in IBD posterior sidecar")
                        self.pair_index[pair] = len(rows)
                        rows.append(i)
        self.row_runs = _aligned_runs(rows)
        self.required = None

    def align(self, positions, identifiers):
        """Match markers exactly where possible and interpolate between covered sites."""
        coords, first_at_coord = np.unique(self.positions, return_index=True)
        self.left = np.full(len(positions), -1, dtype=np.intp)
        self.right = self.left.copy()
        self.fraction = np.zeros(len(positions), dtype=np.float32)
        if not len(coords):
            return
        upper = np.searchsorted(coords, positions)
        bounded = upper < len(coords)
        exact = bounded & (positions == coords[np.minimum(upper, len(coords) - 1)])
        self.left[exact] = first_at_coord[upper[exact]]
        if self.variant_ids is not None:
            duplicate = np.flatnonzero(np.r_[False, np.diff(self.positions) == 0])
            if len(duplicate):
                same_position = np.unique(np.r_[duplicate, duplicate - 1])
                lookup = {(self.positions[row], self.variant_ids[row]): row for row in same_position}
                many = exact & np.isin(positions, self.positions[duplicate])
                for row in np.flatnonzero(many):
                    self.left[row] = lookup.get((positions[row], identifiers[row]), self.left[row])
        self.right[exact] = self.left[exact]
        interpolate = bounded & (~exact) & (upper > 0)
        high = upper[interpolate]
        low = high - 1
        self.left[interpolate] = first_at_coord[low]
        self.right[interpolate] = first_at_coord[high]
        self.fraction[interpolate] = (positions[interpolate] - coords[low]) / (coords[high] - coords[low])

    def read(self, first, stop):
        """Return local pair × marker × state weights, with zero mass outside coverage."""
        n_pairs = len(self.pair_index)
        if not n_pairs:
            return np.empty((0, stop - first, 3), dtype=np.float32)
        left, right = self.left[first:stop], self.right[first:stop]
        present = left >= 0
        source = np.union1d(left[present], right[present])
        values = np.empty((n_pairs, len(source), 3), dtype=np.float32)
        column_runs = _aligned_runs(source)
        for target_rows, source_rows in self.row_runs:
            for target_columns, source_columns in column_runs:
                self.probabilities.read_direct(values, np.s_[source_rows, source_columns, :],
                                               np.s_[target_rows, target_columns, :])
        if np.any(~np.isfinite(values)) or np.any(values < 0) or np.any(values > 1):
            raise ValueError("Sibling IBD posterior probabilities must be finite and between zero and one")
        if np.array_equal(source, left) and np.array_equal(left, right):
            return values
        result = np.zeros((n_pairs, stop - first, 3), dtype=np.float32)
        result[:, present] = values[:, np.searchsorted(source, left[present])]
        blend = present & (left != right)
        fraction = self.fraction[None, first:stop, None][:, blend]
        result[:, blend] *= 1 - fraction
        result[:, blend] += values[:, np.searchsorted(source, right[blend])] * fraction
        return result
