from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd
import pytest
from pysnptools.snpreader import Bed, SnpData

from snipar.config import nan_integer
from snipar.imputation import preprocess_data
from snipar.imputation.impute_from_sibs import ImputationFamilies, impute
from snipar.resources import MIB
from snipar.scripts import impute as impute_script


def test_imputation_memory_plan_budgets_core_markers_and_halo():
    reader = SimpleNamespace(n_samples=1_000_000, bgen=SimpleNamespace(samples=range(1_000_000)),
                             matched_pcs=None)
    budget = 2_048 * MIB
    width = impute_script._plan_imputation(reader, 100_000, 500_000, 5, 8, 1, budget)
    assert 1 <= width < 175
    assert 1 <= reader.read_threads <= 8
    assert 16_000_000 <= reader.probability_bytes <= 128 * MIB
    with pytest.raises(MemoryError, match="minimum workspace"):
        impute_script._plan_imputation(reader, 100_000, 500_000, 5, 8, 1, 64 * MIB)


def test_small_budget_reuses_chromosome_plan_and_output_matches(tmp_path, monkeypatch):
    n_families, n_markers = 2_000, 40
    ids = np.array([f"{prefix}{i}" for i in range(n_families) for prefix in ("c", "f")])
    fids = np.repeat([f"fam{i}" for i in range(n_families)], 2)
    rng = np.random.default_rng(9124)
    a, b, c = rng.binomial(1, .3, (3, n_families, n_markers))
    calls = np.empty((len(ids), n_markers), dtype=float)
    calls[0::2], calls[1::2] = a + b, a + c
    calls[2:4, 10:12] = np.nan
    base = tmp_path / "calls"
    pos = np.column_stack((np.ones(n_markers), np.zeros(n_markers), np.arange(1, n_markers + 1)))
    Bed.write(str(base), SnpData(iid=np.column_stack((fids, ids)),
                                sid=[f"m{i}" for i in range(n_markers)], pos=pos, val=calls), count_A1=True)
    ped = tmp_path / "pedigree.txt"
    pd.DataFrame({"FID": fids[0::2], "IID": ids[0::2], "FATHER_ID": ids[1::2], "MOTHER_ID": "0"}).to_csv(
        ped, sep="\t", index=False,
    )
    variants = []
    plans = []
    original = impute_script.ImputationGenotypes
    read = original.read_calls

    def tracked_plan(*args, **kwargs):
        plan = original(*args, **kwargs)
        plans.append(plan)
        return plan

    def tracked_read(self, markers, **kwargs):
        variants.append(len(markers))
        return read(self, markers, **kwargs)

    monkeypatch.setattr(impute_script, "ImputationGenotypes", tracked_plan)
    monkeypatch.setattr(original, "read_calls", tracked_read)
    args = ["--bed", str(base), "--pedigree", str(ped), "--threads", "2",
            "--output_compression", "gzip", "-silent_progress"]
    for mode, budget in (("full", 256 * MIB), ("bounded", 17 * MIB)):
        monkeypatch.setattr(impute_script, "get_memory_budget", lambda: budget)
        impute_script.main(impute_script.parser.parse_args([*args, "--out", str(tmp_path / mode)]))
    assert len(plans) == 2
    assert variants[0] == n_markers and all(width < n_markers for width in variants[1:])
    assert sum(variants[1:]) == n_markers
    with h5py.File(tmp_path / "full.hdf5") as full, h5py.File(tmp_path / "bounded.hdf5") as bounded:
        assert set(full) == set(bounded)
        for key in full:
            np.testing.assert_array_equal(full[key][:], bounded[key][:], err_msg=key)
    assert not list(tmp_path.glob(".snipar-impute-*"))


def test_prepared_segments_match_unsorted_input_in_original_order():
    families = pd.DataFrame({"FID": [b"fam"], "FATHER_ID": [b"f"], "MOTHER_ID": [b"m"],
                             "IID": [[b"a", b"b"]], "has_father": [False], "has_mother": [False],
                             "single_parent": [False], "sib_count": [2]})
    ids = {b"a": 0, b"b": 1}
    segments = {("a", "b"): [0, 20, 0, 20, 25, 2, 30, 40, 1]}
    pos = np.array([10, 20, 22, 30, 50])
    genotypes = np.array([[1, 2, 1, 1, nan_integer], [0, 1, 1, 2, nan_integer]], dtype=np.int8)
    frequencies = np.broadcast_to(np.linspace(.2, .4, len(pos)), genotypes.shape)
    prepared = ImputationFamilies(families, ids, segments, len(ids))
    same, shuffled = {}, {}
    perm = np.array([2, 1, 0, 4, 3])
    impute(families, None, None, genotypes, None, pos, same, "1", frequencies,
           prepared=prepared, return_genotypes=False, silent_progress=True)
    impute(families, ids, None, genotypes[:, perm], segments, pos[perm], shuffled, "1", frequencies[:, perm],
           silent_progress=True)
    assert same["ratio_ibd0"][1] == 1
    for key in ("imputed_par_gts", "num_observed_parental_alleles", "ratio_ibd0"):
        np.testing.assert_array_equal(same[key], shuffled[key][..., np.argsort(perm)])


def test_bgen_reads_probability_tiles_and_preserves_missing_calls(monkeypatch):
    ids = np.array(["i0", "i1", "i2"])
    haplotypes = np.array([[[1, 1], [0, 0], [0, 1], [0, 0], [0, 1]],
                          [[0, 0], [0, 1], [1, 1], [0, 1], [1, 0]],
                          [[0, 1], [0, 1], [0, 0], [1, 1], [0, 0]]])
    probabilities = np.stack((haplotypes[..., 0], 1 - haplotypes[..., 0],
                              haplotypes[..., 1], 1 - haplotypes[..., 1]), axis=-1).astype(np.float32)
    probabilities[0, 3, :2] = (.99, .01)
    haplotypes[0, 3, 0] = nan_integer
    observed_tiles = []

    def read(indices, dtype):
        samples, variants = indices
        observed_tiles.append(len(variants))
        return np.asfortranarray(probabilities[samples][:, variants], dtype=dtype)

    monkeypatch.setattr(preprocess_data, "open_bgen",
                        lambda *_, **__: SimpleNamespace(samples=ids, read=read))
    bim = pd.DataFrame({"id": list("abcde"), "coordinate": range(1, 6)})
    with preprocess_data.ImputationGenotypes("sample", None, bim, set(ids.astype("S")),
                                             probability_bytes=16 * len(ids) * 2) as genotypes:
        calls, dosages = genotypes.read_calls(np.arange(5))
    assert observed_tiles == [2, 2, 1]
    np.testing.assert_array_equal(calls, haplotypes)
    expected = np.where((haplotypes < 0).any(axis=-1), nan_integer, haplotypes.sum(axis=-1))
    np.testing.assert_array_equal(dosages, expected)


@pytest.mark.parametrize("individual_frequencies", [False, True])
def test_parent_child_fast_path_matches_general_imputation(individual_frequencies):
    # Missing and inconsistent combinations compare diagnostics as well as dosages.
    children = [0, 1, 2, nan_integer]
    genotypes = np.array([(c, p) for c in children for p in children], dtype=np.int8).flatten()
    matrix = np.broadcast_to(genotypes[:, None], (len(genotypes), 7)).copy()
    ids = {f"s{i}".encode(): i for i in range(len(matrix))}
    families = pd.DataFrame({
        "FID": [f"f{i}".encode() for i in range(len(matrix) // 2)],
        "IID": [[iid] for iid in list(ids)[::2]], "FATHER_ID": list(ids)[1::2],
        "MOTHER_ID": b"missing", "has_father": True, "has_mother": False,
        "single_parent": True, "sib_count": 1,
    })
    frequencies = np.broadcast_to(np.linspace(0, 1, matrix.shape[1]), matrix.shape)
    if individual_frequencies:
        frequencies = np.clip(frequencies + np.linspace(-.35, .35, len(matrix))[:, None], 0, 1)
    fake_phases = np.full((*matrix.shape, 2), nan_integer, dtype=np.int8)
    outputs = [{}, {}]
    for phase, out in zip((None, fake_phases), outputs):
        impute(families, ids, phase, matrix, {}, np.arange(matrix.shape[1]), out, "1",
               frequencies, threads=2, silent_progress=True)
    for key in outputs[0]:
        np.testing.assert_array_equal(outputs[0][key], outputs[1][key], err_msg=key)
