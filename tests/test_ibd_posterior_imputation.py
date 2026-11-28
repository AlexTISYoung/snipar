import gzip

import h5py
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed, SnpData

from snipar.config import nan_integer
from snipar.imputation.ibd_posterior import SiblingIBDPosteriors
from snipar.imputation.impute_from_sibs import ImputationFamilies, impute
from snipar.scripts import ibd as ibd_script
from snipar.scripts import impute as impute_script


def _sidecar(path, pairs, positions, variant_ids, probabilities):
    with h5py.File(path, "w") as file:
        file["pairs"] = np.array(pairs, dtype="S")
        file["positions"] = positions
        file["variant_ids"] = np.array(variant_ids, dtype="S")
        file["probabilities"] = np.array(probabilities, dtype=np.float32)


def test_sidecar_reorders_pair_ids_and_resolves_duplicate_positions_by_marker(tmp_path):
    path = tmp_path / "ibd.hdf5"
    probs = np.eye(3, dtype=np.float32)
    _sidecar(path, [["unneeded", "outside"], ["b", "a"]], [100, 100, 200], ["x", "y", "z"], [probs, probs])
    sibships = pd.DataFrame({"IID": [[b"a", b"b"], [b"c", b"d"]]})
    with SiblingIBDPosteriors(path, sibships) as source:
        assert source.pair_index == {(b"a", b"b"): 0}
        source.align(np.array([90, 100, 100, 150, 200, 300]), [b"u", b"x", b"y", b"m", b"z", b"t"])
        np.testing.assert_array_equal(source.read(0, 6)[0],
                                      [[0, 0, 0], [1, 0, 0], [0, 1, 0], [.5, 0, .5], [0, 0, 1], [0, 0, 0]])


def _impute_test_family(genotypes, probabilities, *, parent=False, segments=None, phased=None, **kwargs):
    table = pd.DataFrame({"FID": [b"fam"], "FATHER_ID": [b"dad"], "MOTHER_ID": [b"mom"],
                          "IID": [[b"a", b"b"]], "has_father": [parent], "has_mother": [False],
                          "single_parent": [parent], "sib_count": [2]})
    sample_index = {b"a": 0, b"b": 1, b"dad": 2} if parent else {b"a": 0, b"b": 1}
    rows = {(b"a", b"b"): 0}
    genotypes = np.array(genotypes, dtype=np.int8)
    positions = np.arange(1, genotypes.shape[1] + 1)
    freqs = np.broadcast_to(.25, genotypes.shape)
    prepared = ImputationFamilies(table, sample_index, segments or {}, len(sample_index), posterior_index=rows)
    result = {}
    impute(table, None, phased, genotypes, None, positions, result, "1", freqs, prepared=prepared,
           return_genotypes=False, silent_progress=True, ibd_posteriors=np.array([probabilities], dtype=np.float32),
           posterior_hard_fallback=segments is None, **kwargs)
    return result


def test_two_sibling_posteriors_weight_only_compatible_states_and_keep_legacy_counts():
    result = _impute_test_family([[0, 0, 1, 0], [0, 2, nan_integer, 0]], [[.2, .5, .3]] * 3 + [[0, 0, 1]],
                                 segments={("a", "b"): [1, 4, 0]})
    np.testing.assert_allclose(result["imputed_par_gts"][0], [.1375, 1, .75, .25])
    np.testing.assert_array_equal(result["robust_imputed_par_gts"][0], [0, 1, .75, 0])
    np.testing.assert_array_equal(result["num_observed_parental_alleles"][0], [4, 4, 2, 4])
    np.testing.assert_allclose(result["num_expected_observed_parental_alleles"][0], [2.9, 4, 2, 2])
    np.testing.assert_array_equal(result["posterior_num_sibships"], [1, 1, 0, 1])


def test_known_parent_posterior_counts_average_uncertain_shared_parent():
    result = _impute_test_family([[0], [0], [0]], [[.2, .5, .3]], parent=True,
                                 segments={("a", "b"): [1, 1, 1]})
    np.testing.assert_allclose(result["imputed_par_gts"][0], [.5 / 7 + .3 / 4])
    np.testing.assert_allclose(result["robust_imputed_par_gts"][0], [1 / 7])
    assert np.isnan(result["num_observed_parental_alleles"][0, 0])
    np.testing.assert_allclose(result["num_expected_observed_parental_alleles"][0], [.2 * 4 + .5 * (3 + 3 / 7) + .3 * 3])


def test_posterior_uses_phased_duo_transmissions_within_each_state():
    phased = np.array([[[0, 0]], [[0, 1]], [[0, 1]]], dtype=np.int8)
    origins = np.zeros((2, 1), dtype=np.int8)
    result = _impute_test_family(phased.sum(axis=-1), [[.5, .5, 0]], phased=phased, parent=True,
                                 segments={("a", "b"): [1, 1, 0]}, parent_origin=origins,
                                 parent_haplotype=origins, parent_origin_index={(2, 0): 0, (2, 1): 1})
    np.testing.assert_array_equal(result["imputed_par_gts"][0], [.5])
    np.testing.assert_array_equal(result["robust_imputed_par_gts"][0], [0])
    np.testing.assert_array_equal(result["num_expected_observed_parental_alleles"][0], [4])


def test_impute_auto_discovers_sidecar_and_writes_chunked_outputs(tmp_path):
    positions = [100, 150, 200, 300]
    ids = ["a", "b", "frequency_only"]
    names = ["m1", "m2", "m3", "m4"]
    genotypes = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [2., 2., 2., 2.]])
    bed = tmp_path / "test"
    Bed.write(str(bed), SnpData(iid=np.column_stack((ids, ids)), sid=names,
                                pos=np.column_stack((np.ones(4), np.zeros(4), positions)), val=genotypes),
              count_A1=True)
    pedigree = tmp_path / "pedigree.txt"
    pedigree.write_text("FID IID FATHER_ID MOTHER_ID\nf a dad mom\nf b dad mom\nu frequency_only absent1 absent2\n")
    ibd = tmp_path / "test_ibd"
    with gzip.open(f"{ibd}.segments.gz", "wt") as file:
        file.write("ID1 ID2 Chr start_coordinate stop_coordinate IBDType\na b 1 100 300 0\n")
    args = ["--bed", str(bed), "--pedigree", str(pedigree), "--threads", "1", "-silent_progress"]
    hard = tmp_path / "hard"
    impute_script.main(impute_script.parser.parse_args([*args, "--ibd", str(ibd), "--out", str(hard)]))
    _sidecar(f"{ibd}.posterior.hdf5", [["a", "b"]], [100, 200], ["m1", "m3"], [[[1, 0, 0], [0, 0, 1]]])
    for chunks in ("1", "3"):
        output = tmp_path / f"soft_{chunks}"
        impute_script.main(impute_script.parser.parse_args([*args, "--ibd", str(ibd), "--out", str(output),
                                                          "--chunks", chunks]))
        with h5py.File(f"{output}.hdf5", "r") as source, h5py.File(f"{hard}.hdf5", "r") as legacy:
            np.testing.assert_array_equal(source["imputed_par_gts"][0], np.array([0, 1 / 6, 1 / 3, 0], dtype=np.float16))
            np.testing.assert_array_equal(source["robust_imputed_par_gts"][:], legacy["imputed_par_gts"][:])
            np.testing.assert_array_equal(source["num_observed_parental_alleles"][:], legacy["num_observed_parental_alleles"][:])
            np.testing.assert_array_equal(source["posterior_num_sibships"][:], [1, 1, 1, 0])
    posterior_only = tmp_path / "posterior_only"
    impute_script.main(impute_script.parser.parse_args([*args, "--ibd-posterior", f"{ibd}.posterior.hdf5",
                                                      "--out", str(posterior_only)]))
    with h5py.File(f"{posterior_only}.hdf5", "r") as source:
        np.testing.assert_array_equal(source["num_observed_parental_alleles"][0], [4, 4, 2, np.nan])


def test_ibd_cli_sidecar_is_consumed_and_larger_families_keep_hard_estimates(tmp_path):
    rng = np.random.default_rng(221)
    positions = np.arange(101, 141)
    ids = ["a", "b", "c", "d", "e"]
    genotypes = rng.binomial(2, .5, size=(len(ids), len(positions)))
    bed = tmp_path / "test"
    Bed.write(str(bed), SnpData(iid=np.column_stack((ids, ids)), sid=[f"m{i}" for i in positions],
                                pos=np.column_stack((np.ones(len(positions)), np.zeros(len(positions)), positions)),
                                val=genotypes.astype(float)), count_A1=True)
    pedigree = tmp_path / "pedigree.txt"
    pedigree.write_text("FID IID FATHER_ID MOTHER_ID\nsmall a f1 m1\nsmall b f1 m1\n"
                        "larger c f2 m2\nlarger d f2 m2\nlarger e f2 m2\n")
    genetic_map = tmp_path / "map.txt"
    genetic_map.write_text("pposition gposition\n" + "".join(f"{p} {i * .03}\n" for i, p in enumerate(positions)))
    base = tmp_path / "related"
    common = ["--bed", str(bed), "--pedigree", str(pedigree), "--threads", "1"]
    ibd_script.main(ibd_script.parser.parse_args([*common, "--p_error", ".001", "--posterior", "--map",
                                                 str(genetic_map), "--out", str(base)]))
    segments = str(next(tmp_path.glob("related_chr_*.ibd.segments.gz"))).removesuffix(".segments.gz")
    output = tmp_path / "result"
    impute_script.main(impute_script.parser.parse_args([*common, "--ibd", str(segments), "--out", str(output),
                                                       "--chunks", "3", "-silent_progress"]))
    with h5py.File(f"{output}.hdf5", "r") as source:
        row = np.flatnonzero(source["families"][:] == b"larger")[0]
        np.testing.assert_array_equal(source["imputed_par_gts"][row], source["robust_imputed_par_gts"][row])
        assert 0 < source["posterior_num_sibships"][:].sum()
        assert source["posterior_num_sibships"][:].max() == 1
