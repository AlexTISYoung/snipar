"""Store and reuse duoHMM labels in the original genotype haplotype order."""

from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
from numba import njit

from snipar.imputation.ibd_posterior import _aligned_runs
from snipar.utilities import encode_str_array


@njit(cache=True, inline="always")
def _phase_code(allele, flipped):
    return (1 - allele if flipped else allele) if allele == 0 or allele == 1 else 2


@njit(cache=True)
def _pack_duos(phased, pairs, markers, flipped):
    evidence = np.empty((len(pairs), len(markers)), dtype=np.uint8)
    for row in range(len(pairs)):
        parent, child = pairs[row]
        for column in range(len(markers)):
            marker = markers[column]
            flip = flipped[column]
            evidence[row, column] = (_phase_code(phased[parent, marker, 0], flip)
                                     | _phase_code(phased[parent, marker, 1], flip) << 2
                                     | _phase_code(phased[child, marker, 0], flip) << 4
                                     | _phase_code(phased[child, marker, 1], flip) << 6)
    return evidence


def _pack_phase_evidence(phased, pair_indices, markers, flipped=None):
    phased = np.asarray(phased)
    pairs = np.asarray(pair_indices, dtype=np.intp)
    if (phased.ndim != 3 or phased.shape[2] != 2 or pairs.ndim != 2 or pairs.shape[1] != 2
            or (pairs < 0).any() or (pairs >= len(phased)).any()):
        raise ValueError("Phase evidence requires original haplotypes and valid directed duo rows")
    if flipped is None:
        flipped = np.zeros(len(markers), dtype=bool)
    return _pack_duos(phased, pairs, markers, flipped)


def _check_labels(origin, copy):
    if ((origin < -1).any() or (origin > 1).any() or (copy < -1).any() or (copy > 1).any()
            or ((origin == -1) & (copy != -1)).any()):
        raise ValueError("Invalid parent-of-origin or parent-copy label in origin sidecar")


def _chromosome_label(chromosome):
    return str(chromosome).removeprefix("chr").removesuffix(".0")


class OriginWriter:
    """Write directed (parent, child) labels and phase evidence in bounded tiles."""

    def __init__(self, path, chromosome, pairs, positions, variant_ids, alleles=None):
        self.path = Path(path)
        self.chromosome = _chromosome_label(chromosome)
        self.pairs = encode_str_array(pairs)
        self.positions = np.asarray(positions, dtype=np.int64)
        self.variant_ids = encode_str_array(variant_ids)
        self.alleles = None if alleles is None else encode_str_array(alleles)
        if (self.pairs.ndim != 2 or self.pairs.shape[1] != 2 or not len(self.pairs)
                or self.positions.ndim != 1 or not len(self.positions)
                or self.variant_ids.shape != self.positions.shape
                or self.alleles is not None and self.alleles.shape != (len(self.positions), 2)):
            raise ValueError("Expected directed duos and aligned origin marker metadata")

    def __enter__(self):
        self.directory = TemporaryDirectory(prefix=".snipar-origin-", dir=self.path.parent)
        self.temporary = Path(self.directory.name) / self.path.name
        try:
            self.file = h5py.File(self.temporary, "w")
            self.file.attrs.update(format_version=1, chromosome=self.chromosome)
            self.file["pairs"] = self.pairs
            self.file["positions"] = self.positions
            self.file["variant_ids"] = self.variant_ids
            if self.alleles is not None:
                self.file["alleles"] = self.alleles
            shape = len(self.pairs), len(self.positions)
            chunks = min(8, shape[0]), min(256, shape[1])
            self.file.create_dataset("origin", shape, dtype=np.int8, chunks=chunks, compression="lzf", fillvalue=-1)
            self.file.create_dataset("parent_copy", shape, dtype=np.int8, chunks=chunks, compression="lzf",
                                     fillvalue=-1)
            self.file.create_dataset("phase_evidence", shape, dtype=np.uint8, chunks=chunks, compression="lzf",
                                     fillvalue=255)
        except BaseException:
            self.__exit__(True, None, None)
            raise
        return self

    def write(self, rows, markers, origin, parent_copy, phased, pair_indices):
        """Store a duo and marker slice; phased and pair indices use tile-local rows."""
        parent_copy = np.asarray(parent_copy)
        origin = np.asarray(origin)
        shape = (len(range(*rows.indices(len(self.pairs)))), len(range(*markers.indices(len(self.positions)))))
        if (origin.shape != shape or parent_copy.shape != shape or tuple(phased.shape[1:]) != (shape[1], 2)
                or len(pair_indices) != shape[0] or origin.dtype != np.int8 or parent_copy.dtype != np.int8):
            raise ValueError("Origin tile labels must be int8 and aligned with the selected duos and markers")
        _check_labels(origin, parent_copy)
        evidence = _pack_phase_evidence(phased, pair_indices, np.arange(shape[1], dtype=np.intp))
        self.file["origin"][rows, markers] = origin
        self.file["parent_copy"][rows, markers] = parent_copy
        self.file["phase_evidence"][rows, markers] = evidence

    def __exit__(self, exc_type, exc, traceback):
        try:
            if hasattr(self, "file"):
                self.file.close()
            if exc_type is None:
                self.temporary.replace(self.path)
        finally:
            self.directory.cleanup()


class OriginReader:
    """Match stored directed duos and markers without loading a chromosome of labels."""

    def __init__(self, path, chromosome, pairs):
        self.path = path
        self.chromosome = _chromosome_label(chromosome)
        self.pairs = encode_str_array(pairs)

    def __enter__(self):
        self.file = h5py.File(self.path, "r")
        try:
            self._load_duos()
        except BaseException:
            self.file.close()
            raise
        return self

    def __exit__(self, *exc):
        self.file.close()

    def _load_duos(self):
        required = {"pairs", "positions", "variant_ids", "origin", "parent_copy", "phase_evidence"}
        if self.file.attrs.get("format_version") != 1 or not required.issubset(self.file):
            raise ValueError(f"Unsupported origin sidecar format in {self.path}")
        if _chromosome_label(self.file.attrs.get("chromosome")) != self.chromosome:
            raise ValueError(f"Origin sidecar chromosome differs from genotype input in {self.path}")
        pairs, positions, variant_ids = (self.file[key] for key in ("pairs", "positions", "variant_ids"))
        if pairs.ndim != 2 or pairs.shape[1] != 2 or positions.ndim != 1 or variant_ids.shape != positions.shape:
            raise ValueError(f"Invalid origin sidecar duo or marker metadata in {self.path}")
        shape = len(pairs), len(positions)
        for key, dtype in (("origin", np.int8), ("parent_copy", np.int8), ("phase_evidence", np.uint8)):
            if self.file[key].shape != shape or self.file[key].dtype != dtype:
                raise ValueError(f"Invalid {key} dataset in origin sidecar {self.path}")
        requested = {tuple(pair): row for row, pair in enumerate(self.pairs)}
        selected = np.zeros(len(self.pairs), dtype=bool)
        source_rows, target_rows = [], []
        for first in range(0, len(pairs), 8192):
            for row, pair in enumerate(pairs[first:first + 8192], first):
                target = requested.get(tuple(pair))
                if target is None:
                    continue
                if selected[target]:
                    raise ValueError(f"Duplicate directed duo in origin sidecar {self.path}")
                source_rows.append(row)
                target_rows.append(target)
                selected[target] = True
        self.row_runs = _aligned_runs(source_rows)
        self.target_rows = np.array(target_rows, dtype=np.intp)
        self.missing_rows = np.flatnonzero(~selected)

    def align(self, positions, variant_ids, alleles=None):
        """Match each marker by both its physical position and variant ID."""
        source_pos, source_ids = self.file["positions"][:], self.file["variant_ids"][:]
        index = {}
        for row, key in enumerate(zip(source_pos, source_ids)):
            if key in index:
                raise ValueError(f"Duplicate marker position and ID in origin sidecar {self.path}")
            index[key] = row
        variant_ids = encode_str_array(variant_ids)
        self.columns = np.fromiter((index.get(key, -1) for key in zip(positions, variant_ids)),
                                   dtype=np.intp, count=len(positions))
        present = self.columns >= 0
        self.flipped = np.zeros(len(positions), dtype=bool)
        if alleles is not None and "alleles" in self.file:
            sidecar_alleles = self.file["alleles"][:]
            if sidecar_alleles.shape != (len(source_pos), 2):
                raise ValueError(f"Invalid marker allele metadata in origin sidecar {self.path}")
            expected = sidecar_alleles[self.columns[present]]
            observed = encode_str_array(alleles)[present]
            same = np.all(expected == observed, axis=1)
            reversed_order = np.all(expected == observed[:, ::-1], axis=1)
            if not np.all(same | reversed_order):
                raise ValueError(f"Origin sidecar alleles differ from phased genotype input in {self.path}")
            self.flipped[present] = ~same
        return int(present.sum())

    def _read_tile(self, key, source_columns, dtype):
        tile = np.empty((len(self.target_rows), len(source_columns)), dtype=dtype)
        columns = _aligned_runs(source_columns)
        dataset = self.file[key]
        for target_rows, source_rows in self.row_runs:
            for target_markers, source_markers in columns:
                dataset.read_direct(tile, np.s_[source_rows, source_markers], np.s_[target_rows, target_markers])
        return tile

    def read_into(self, first, stop, phased, pair_indices, origin, parent_copy):
        """Validate local phase before copying known labels into the imputation tile."""
        mapped = self.columns[first:stop]
        target_columns = np.flatnonzero(mapped >= 0)
        if not len(target_columns) or not len(self.target_rows):
            return
        target_columns = target_columns[np.argsort(mapped[target_columns])]
        source_columns = mapped[target_columns]
        saved = self._read_tile("phase_evidence", source_columns, np.uint8)
        available = saved != 255
        actual = _pack_phase_evidence(phased, pair_indices[self.target_rows], target_columns,
                                      self.flipped[first:stop][target_columns])
        differences = available & (saved != actual)
        if differences.any():
            row, column = np.unravel_index(np.argmax(differences), differences.shape)
            parent, child = self.pairs[self.target_rows[row]]
            marker = self.file["variant_ids"][source_columns[column]]
            raise ValueError(f"Origin sidecar phased genotypes differ for {parent!r}/{child!r} at {marker!r}. "
                             f"Use the original phased genotypes or regenerate {self.path}.")
        del saved, actual, differences
        origin_tile = self._read_tile("origin", source_columns, np.int8)
        copy_tile = self._read_tile("parent_copy", source_columns, np.int8)
        _check_labels(origin_tile, copy_tile)
        target = np.ix_(self.target_rows, target_columns)
        origin[target] = np.where(available, origin_tile, -1)
        parent_copy[target] = np.where(available, copy_tile, -1)
