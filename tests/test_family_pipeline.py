"""Regression tests for chromosome-wide family inference and simulation."""
import argparse
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from snipar.config import nan_integer
from snipar.ibd import find_segments, smooth_segments
from snipar.imputation import preprocess_data
from snipar.imputation.impute_from_sibs import impute
from snipar.map import get_map_positions, pos_to_cM
from snipar.numrange import parseNumRange
from snipar.pedigree import create_pedigree, find_individuals_with_sibs, get_sibpairs_from_ped
from snipar.preprocess import get_indices_given_ped
from snipar.scripts import ibd, simulate
from snipar.simulate import am_indices, simulate_freqs


def test_map_lookup_supports_unsorted_positions_and_boundaries(tmp_path):
    positions = np.array([39, 10, 20, 40, 9])
    expected = np.array([2.1, 0.7, 2.1, np.nan, np.nan])
    np.testing.assert_array_equal(pos_to_cM(positions, np.array([10, 20, 40]), np.array([.7, 2.1])),
                                  expected)
    source = tmp_path / "genetic_map.txt"
    source.write_text("pposition\tgposition\n10\t0.7\n20\t2.1\n")
    numeric = SimpleNamespace(pos=np.array([20, 15, 10]))
    np.testing.assert_allclose(get_map_positions(source, numeric), [2.1, 1.4, .7])


def test_segments_cover_every_marker_including_terminal_and_single_snp():
    snps = np.array(["a", "b", "c", "d"])
    points = np.arange(4)
    single = find_segments(np.array([2]), points[:1], snps[:1], points[:1])
    assert [(s.start, s.end, s.start_snp, s.end_snp) for s in single] == [(0, 1, "a", "a")]
    path, segments = smooth_segments(np.array([1, 1, 2, 2]), points, snps, points, 2)
    assert list(path) == [2, 2, 1, 1]
    assert [(s.start, s.end, s.end_bp) for s in segments] == [(0, 2, 1), (2, 4, 3)]


@pytest.mark.parametrize("expression", ["1-2junk", "3-1", "2.3", "2 4"])
def test_ranges_reject_partial_or_descending_values(expression):
    with pytest.raises(argparse.ArgumentTypeError):
        parseNumRange(expression)


def test_pedigree_uses_parent_pairs_without_ambiguous_concatenation():
    ped = np.array([["a", "a1", "1", "23"], ["a", "a2", "1", "23"],
                    ["b", "b1", "12", "3"], ["b", "b2", "12", "3"]])
    pairs, _ = get_sibpairs_from_ped(ped)
    np.testing.assert_array_equal(pairs, [["a1", "a2"], ["b1", "b2"]])
    empty = find_individuals_with_sibs(np.array(["z"]), ped, np.array(["z"]), True)
    assert empty.shape == (0,)


def test_king_pedigree_unifies_sibling_groups_and_preserves_numeric_ids(tmp_path):
    king = tmp_path / "in.kin"
    king.write_text("FID1\tID1\tFID2\tID2\tInfType\n"
                    "f\t01\t02\t02\tFS\n"
                    "f\t03\t04\t04\tFS\n"
                    "f\t02\t03\t03\tFS\n"
                    "p\t05\tf\t01\tPO\n")
    age = tmp_path / "agesex"
    age.write_text("FID IID age sex\nf 01 20 M\nf 02 20 F\nf 03 20 M\nf 04 20 F\np 05 32 M\n")
    ped = create_pedigree(king, age, same_parents_in_ped=False).set_index("IID")
    assert ped.loc[["01", "02", "03", "04"], "FID"].nunique() == 1
    assert ped.loc["01", "FATHER_ID"] == "05"


def test_sibling_genotype_selection_keeps_requested_cohort():
    ped = np.array([["f", "1", "p", "m"], ["f", "2", "p", "m"]])
    gts_ids = np.array(["1", "2", "p", "m"])
    selected, rows, _, _ = get_indices_given_ped(ped, gts_ids, ids=np.array(["2"]), sib=True)
    np.testing.assert_array_equal(selected, ["2"])
    np.testing.assert_array_equal(rows, [1, 2, 3])


def test_ibd_batches_reused_for_each_chromosome(tmp_path, monkeypatch):
    for chrom in (1, 2):
        (tmp_path / f"g{chrom}.bed").touch()
    king = tmp_path / "in.kin"
    king.write_text("ID1\tID2\tInfType\n01\t02\tFS\n03\t04\tFS\n")
    calls = []
    monkeypatch.setattr(ibd.snipar.ibd, "infer_ibd_chr",
                        lambda pairs, *args, **kwargs: calls.append((kwargs["chrom"], tuple(pairs[0]))))
    args = ibd.parser.parse_args(["--bed", str(tmp_path / "g@"), "--king", str(king),
                                  "--p_error", "0.001", "--threads", "1", "--batches", "2"])
    ibd.main(args)
    assert calls == [(1, ("01", "02")), (1, ("03", "04")),
                     (2, ("01", "02")), (2, ("03", "04"))]


def test_native_imputation_diagnostics_match_with_multiple_threads():
    n = 250
    groups = [[f"sib{i}a".encode(), f"sib{i}b".encode()] for i in range(n)]
    sibships = pd.DataFrame({"FID": [str(i).encode() for i in range(n)], "FATHER_ID": b"father",
                            "MOTHER_ID": b"mother", "IID": groups, "has_father": False,
                            "has_mother": False, "single_parent": False, "sib_count": 2})
    ids = {iid: i for i, iid in enumerate(iid for pair in groups for iid in pair)}
    ibd = {(a.decode(), b.decode()): [10, 40, 0] for a, b in groups}
    dosages = np.ones((2 * n, 2), dtype=np.int8)
    dosages[2:, 1] = nan_integer
    freqs = np.full(dosages.shape, .3)
    outputs = [{}, {}]
    for threads, output in enumerate(outputs, 1):
        impute(sibships, ids, None, dosages, ibd, np.array([10, 20]), output, "1",
               freqs, threads=threads, silent_progress=True)
    for output in outputs:
        np.testing.assert_array_equal(output["ratio_ibd0"], [1, 1])
        np.testing.assert_array_equal(np.asarray(output["imputed_par_gts"])[0], [1, 1])
        assert np.isnan(np.asarray(output["imputed_par_gts"])[1:, 1]).all()
    with pytest.raises(ValueError, match="Missing genotype for sibling"):
        impute(sibships, {}, None, dosages, ibd, np.array([10, 20]), {}, "1", freqs)


def test_mixed_imputation_aligns_bgen_rows_markers_alleles_and_pcs(tmp_path, monkeypatch):
    from pysnptools.snpreader import Bed, SnpData
    base = str(tmp_path / "geno")
    known = np.array([[0, 1], [1, 2], [2, 0]])
    ids = np.array(["001", "002", "003"])
    Bed.write(base, SnpData(iid=np.column_stack((ids, ids)), sid=["a", "b"],
                            pos=[[1, 0, 10], [1, 0, 20]], val=known.astype(float)), count_A1=True)
    desired = np.array([[[0, 0], [1, 0]], [[0, 1], [1, 1]], [[1, 1], [0, 0]]])
    calls = desired[[2, 0, 1]][:, [1, 0]].copy()
    calls[:, 0] = 1 - calls[:, 0]
    probs = np.stack([calls[..., 0], 1 - calls[..., 0], calls[..., 1], 1 - calls[..., 1]], axis=2)
    bgen = SimpleNamespace(samples=ids[[2, 0, 1]], ids=np.array(["b", "a"]), rsids=np.array(["b", "a"]),
                           phased=[True, True],
                           positions=np.array([20, 10]), allele_ids=np.array(["G,A", "A,G"]),
                           read=lambda indexes, dtype: probs[indexes[0]][:, indexes[1]].astype(dtype))
    monkeypatch.setattr(preprocess_data, "open_bgen", lambda *args, **kwargs: bgen)
    bim = pd.DataFrame({"id": ["a", "b"], "coordinate": [10, 20],
                        "allele1": ["A", "A"], "allele2": ["G", "G"]})
    phased, unphased, rows, _, freqs, output = preprocess_data.prepare_gts(
        base, base, bim, np.empty((0, 0)), set(ids.astype("S")), ["1"],
        pcs=np.array([[.8], [.2], [.4]]), pc_ids=[b"003", b"001", b"002"])
    assert list(rows) == list(ids.astype("S"))
    np.testing.assert_array_equal(unphased, known)
    np.testing.assert_array_equal(phased, desired)
    np.testing.assert_array_equal(output["maf_x"][:, 0], [.2, .4, .8])
    expected, _ = preprocess_data.estimate_f(known, np.array([[.2], [.4], [.8]]))
    np.testing.assert_allclose(freqs, expected)


@pytest.mark.parametrize("r", [-.8, 0, .8])
def test_mating_handles_zero_negative_and_positive_correlation(r):
    np.random.seed(515)
    rng = np.random.default_rng(20)
    a, b = rng.normal(size=(2, 4096))
    male, female = am_indices(12 * a, b, r)
    assert abs(np.corrcoef(a[male], b[female])[0, 1] - r) < .08


def test_vertical_transmission_without_assortative_mating(tmp_path):
    simulate.main(simulate.parser.parse_args([
        "12", ".5", str(tmp_path) + "/", "--nfam", "20", "--beta_vert", ".1", "--min_maf", ".5",
    ]))
    assert (tmp_path / "phenotype.txt").is_file()
    np.testing.assert_array_equal(simulate_freqs(4, min_maf=.5), .5)


def test_hmm_matches_markerwise_reference_with_missing_calls():
    from snipar.ibd import infer_ibd, make_dynamic, viterbi

    rng = np.random.default_rng(282)
    for n in (1, 31, 251):
        gts = rng.integers(0, 3, (20, n)).astype(float)
        gts[rng.random(gts.shape) < .1] = np.nan
        pairs = np.vstack((np.arange(0, 20, 2), np.arange(1, 20, 2))).T.copy()
        freqs = rng.uniform(.01, .99, n)
        genetic_map = np.cumsum(rng.uniform(0, 2, n))
        weights = rng.uniform(.2, 1, n)
        errors = rng.uniform(0, .05, n)
        expected = np.vstack([viterbi(*make_dynamic(gts[a], gts[b], freqs, genetic_map, weights, errors))
                              for a, b in pairs])
        np.testing.assert_array_equal(infer_ibd(pairs, gts, freqs, genetic_map, weights, errors), expected)


def test_ld_moments_match_pairwise_complete_reference():
    from snipar.ld import compute_ld_scores, r2_est

    rng = np.random.default_rng(5)
    gts = rng.integers(0, 3, (500, 30)).astype(float)
    gts[rng.random(gts.shape) < .05] = np.nan
    genetic_map = np.arange(30) / 10
    expected = np.ones(30)
    for i in range(30):
        for j in range(30):
            if i == j or abs(genetic_map[i] - genetic_map[j]) >= 1:
                continue
            valid = np.isfinite(gts[:, i]) & np.isfinite(gts[:, j])
            rho2 = np.corrcoef(gts[valid, i], gts[valid, j])[0, 1] ** 2
            expected[i] += rho2 - (1 - rho2) / (valid.sum() - 2)
    np.testing.assert_allclose(compute_ld_scores(gts, genetic_map), expected, atol=1e-12)
    assert r2_est(gts[:, 0], np.repeat(.25, 500)) == 0
    assert r2_est(np.array([0, 1.]), np.array([1, 2.])) == 0


def test_pc_regression_excludes_missing_calls_instead_of_treating_as_zero():
    pc_scores = np.array([[-1.], [0.], [1.], [2.]])
    calls = np.array([[0], [1], [2], [nan_integer]])
    freqs, metrics = preprocess_data.estimate_f(calls, pc_scores)
    np.testing.assert_allclose(freqs[:, 0], [0, .5, 1, 1], atol=1e-14)
    np.testing.assert_allclose(metrics["coefs"][:, 0], [.5, .5])
    assert metrics["R2_1"][0] == 1


def test_phased_imputation_chunks_match_single_slice(tmp_path):
    from pathlib import Path

    import h5py
    import snipar
    from snipar.scripts.impute import main, parser

    data = Path(snipar.__file__).parent / "tests/test_data"
    common = ["--bgen", str(data / "sample_reduced@"), "--chr_range", "1",
              "--ibd", str(data / "sample@.our"), "--pedigree", str(data / "sample.ped"),
              "--start", "7", "--end", "151", "--threads", "2", "-silent_progress"]
    whole = tmp_path / "whole"
    chunked = tmp_path / "chunked"
    main(parser.parse_args(common + ["--out", str(whole)]))
    main(parser.parse_args(common + ["--out", str(chunked), "--chunks", "4",
                                    "--output_compression", "gzip"]))
    with h5py.File(f"{whole}.hdf5") as a, h5py.File(f"{chunked}.hdf5") as b:
        assert a.keys() == b.keys()
        for key in a:
            if a[key].dtype.kind in "fc":
                np.testing.assert_allclose(a[key][:], b[key][:], atol=1e-14, equal_nan=True, err_msg=key)
            else:
                np.testing.assert_array_equal(a[key][:], b[key][:], err_msg=key)
        assert b["imputed_par_gts"].compression == "gzip"
        assert b["num_observed_parental_alleles"].compression == "gzip"


def test_simulated_transmissions_have_expected_ibd_and_alleles():
    from snipar.simulate import produce_next_gen_unlinked

    parents = np.broadcast_to([False, True], (100, 33, 2)).copy()
    rows = np.arange(len(parents))
    offspring, ibd = produce_next_gen_unlinked(rows, rows, parents, parents)
    assert abs(offspring.mean() - .5) < .025
    np.testing.assert_array_equal(ibd, offspring[:, 0] == offspring[:, 1])
    proportions = np.bincount(ibd.sum(axis=2).ravel().astype(int)) / ibd[:, :, 0].size
    np.testing.assert_allclose(proportions, [.25, .5, .25], atol=.04)


def test_streamed_segments_preserve_text_output(tmp_path):
    import gzip

    from snipar.ibd import write_segs_from_matrix

    pairs = np.array([["a", "b"], ["c", "d"]])
    snps = np.array(["123", "rs298", "a.d", "tt"])
    gmap = np.array([.25, 1.0, 1.125, 1.75])
    pos = np.array([1, 2, 20, 300])
    ibd = np.array([[1, 1, 0, 2], [0, 1, 1, 1]], dtype=np.int8)
    in_memory, streamed = tmp_path / "array.gz", tmp_path / "stream.gz"
    write_segs_from_matrix(ibd, pairs, snps, pos, gmap, 1, in_memory)
    write_segs_from_matrix(ibd, pairs, snps, pos, gmap, 1, streamed, return_segments=False)
    assert gzip.decompress(in_memory.read_bytes()) == gzip.decompress(streamed.read_bytes())


@pytest.mark.parametrize("indirect", [False, True])
def test_simulator_exports_founder_genotypes_without_additional_generations(tmp_path, indirect):
    from pysnptools.snpreader import Bed

    command = ["20", ".5", str(tmp_path) + "/", "--nfam", "25", "--save_par_gts", "--maf", ".3"]
    if indirect:
        command.extend(["--v_indir", ".2", "--r_dir_indir", "0"])
    simulate.main(simulate.parser.parse_args(command))
    parents = Bed(str(tmp_path / "chr_1_par"), count_A1=True).read()
    pedigree = pd.read_csv(tmp_path / "pedigree.txt", sep=r"\s+", dtype=str)
    children = pedigree.loc[pedigree.FID.str.startswith("1_")]
    assert len(parents.iid) == 50
    assert set(parents.iid[:, 1]) == set(children.FATHER_ID) | set(children.MOTHER_ID)
    founders = pedigree.loc[pedigree.FID.str.startswith("0_")].set_index("IID")
    weights = pd.read_csv(tmp_path / "causal_effects.txt", sep=r"\s+")
    np.testing.assert_allclose(parents.val @ weights.direct, founders.loc[parents.iid[:, 1], "DIRECT"].astype(float),
                               atol=1e-10)
    relationship = pd.read_csv(tmp_path / "last_gen_rel.seg", sep=r"\s+")
    assert list(relationship.InfType.unique()) == ["FS"]
