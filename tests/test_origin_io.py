"""Check cached duoHMM labels against real marker and haplotype order."""

import gzip

import h5py
import numpy as np
import pytest

from snipar.imputation import duohmm
from snipar.imputation.origin_io import OriginReader, OriginWriter
from snipar.scripts import impute


def test_origin_reader_matches_directed_pairs_and_retains_haplotype_slots(tmp_path):
    path = tmp_path / "parental.hdf5"
    source = np.array([[[0, 1], [1, 0], [1, -11]], [[0, 1], [0, 1], [1, 0]],
                       [[1, 0], [1, 0], [0, 1]]], dtype=np.int8)
    positions = np.array([100, 100, 103])
    ids = ["m0", "m1", "m2"]
    alleles = np.array([["A", "G"], ["A", "T"], ["G", "C"]])
    labels = np.array([[0, 1, 1], [1, 0, -1]], dtype=np.int8)
    pairs = np.array([[0, 1], [0, 2]])
    with OriginWriter(path, 1, [["p", "a"], ["p", "b"]], positions, ids, alleles=alleles) as writer:
        for row in range(len(pairs)):
            block = slice(row, row + 1)
            writer.write(block, slice(None), labels[block], labels[block], source, pairs[block])

    marker_order = np.array([2, 0, 1])
    sample_order = np.array([2, 0, 1])
    phased = source[sample_order][:, marker_order]
    valid = phased[:, 1] >= 0
    phased[:, 1][valid] = 1 - phased[:, 1][valid]
    current_alleles = alleles[marker_order]
    current_alleles[1] = current_alleles[1, ::-1]
    names = [["p", "b"], ["p", "a"], ["a", "p"]]
    new_pairs = np.array([[1, 0], [1, 2], [2, 1]])
    with OriginReader(path, "chr1", names) as reader:
        assert reader.align(positions[marker_order], ["m2", "m0", "absent"], current_alleles) == 2
        np.testing.assert_array_equal(reader.missing_rows, [2])
        origin = np.full((3, 3), -1, dtype=np.int8)
        parent_copy = origin.copy()
        reader.read_into(0, 3, phased, new_pairs, origin, parent_copy)
        np.testing.assert_array_equal(origin, [[-1, 1, -1], [1, 0, -1], [-1, -1, -1]])
        np.testing.assert_array_equal(origin, parent_copy)

        incorrect = phased.copy()
        incorrect[2, 0] = incorrect[2, 0, ::-1]
        with pytest.raises(ValueError, match="phased genotypes differ"):
            reader.read_into(0, 3, incorrect, new_pairs, origin, parent_copy)

    with OriginReader(path, 1, names) as reader:
        incorrect_alleles = alleles.copy()
        incorrect_alleles[0] = ["C", "G"]
        with pytest.raises(ValueError, match="alleles differ"):
            reader.align(positions, ids, incorrect_alleles)


@pytest.mark.parametrize("uncached_duo", [False, True])
def test_impute_reuses_origin_without_refitting_and_rejects_rephased_genotypes(tmp_path, monkeypatch, uncached_duo):
    positions = np.arange(101, 293)
    midpoint = len(positions) // 2
    phased = np.empty((2, len(positions), 2), dtype=np.int8)
    phased[0] = [1, 0]
    phased[0, ::4] = [1, 1]
    phased[1, :midpoint] = [1, 0]
    phased[1, midpoint:] = [0, 1]
    names = ["p", "c"]
    if uncached_duo:
        phased = np.concatenate((phased, phased))
        names.extend(("q", "d"))
    sample_order = [1, 0, 3, 2] if uncached_duo else [1, 0]
    vcf = tmp_path / "offspring.vcf"

    def write_vcf(genotypes, order):
        header = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        with vcf.open("w") as output:
            output.write(header + "\t".join(names[i] for i in sample_order) + "\n")
            for i in order:
                calls = ["|".join(str(1 - call) for call in genotypes[j, i]) for j in sample_order]
                output.write(f"1\t{positions[i]}\ts{i}\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")

    pedigree = tmp_path / "pedigree.txt"
    pedigree.write_text("FID IID FATHER_ID MOTHER_ID\nf c p 0\n" + ("g d q 0\n" if uncached_duo else ""))
    prefix = tmp_path / "segments"
    with gzip.open(f"{prefix}.segments.gz", "wt") as output:
        output.write("ID1 ID2 Chr start_coordinate stop_coordinate IBDType\n")
    common = ["--vcf", str(vcf), "--pedigree", str(pedigree), "--threads", "1", "--chunks", "5",
              "-silent_progress", "--ibd", str(prefix)]
    write_vcf(phased, np.arange(len(positions)))
    impute.main(impute.parser.parse_args([*common, "--out", str(tmp_path / "local")]))

    pair_indices = np.array([[0, 1]])
    recombination = duohmm.recombination_probabilities(positions)
    origins, copies, _ = duohmm.call_parental_origin(phased, pair_indices, recombination)
    assert set(np.unique(origins)) == {0, 1}
    sidecar = tmp_path / "segments.origin.hdf5"
    with OriginWriter(sidecar, 1, [["p", "c"]], positions, [f"s{i}" for i in range(len(positions))],
                      alleles=np.tile(["A", "G"], (len(positions), 1))) as writer:
        writer.write(slice(None), slice(None), origins, copies, phased, pair_indices)

    def check_refit(function):
        def checked(haplotypes, fitted_pairs, *args, **kwargs):
            assert uncached_duo, "An origin sidecar covering the cohort should bypass local duoHMM"
            np.testing.assert_array_equal(fitted_pairs, [[3, 2]])
            return function(haplotypes, fitted_pairs, *args, **kwargs)
        return checked

    for method in ("call_parental_origin", "advance_parental_origin"):
        monkeypatch.setattr(duohmm, method, check_refit(getattr(duohmm, method)))
    order = np.random.default_rng(2).permutation(len(positions))
    write_vcf(phased, order)
    for name, option in (("discovered", []), ("explicit", ["--parent-origin", str(sidecar)])):
        impute.main(impute.parser.parse_args([*common, *option, "--out", str(tmp_path / name)]))
        with h5py.File(tmp_path / "local.hdf5") as local, h5py.File(tmp_path / f"{name}.hdf5") as reused:
            for key in ("imputed_par_gts", "num_observed_parental_alleles", "standard_f", "pos"):
                np.testing.assert_array_equal(reused[key][:], local[key][:][..., order], err_msg=key)
            np.testing.assert_array_equal(reused["imputed_par_gts"][0, np.isin(order, [midpoint - 1, midpoint + 1])],
                                          [.5, .5])

    changed = phased.copy()
    changed[1, midpoint] = changed[1, midpoint, ::-1]
    write_vcf(changed, order)
    args = impute.parser.parse_args([*common, "--out", str(tmp_path / "rephased")])
    with pytest.raises(ValueError, match="phased genotypes differ"):
        impute.main(args)
    assert not (tmp_path / "rephased.hdf5").exists()
    # The opt-out disables both cached and local origin calls.
    impute.main(impute.parser.parse_args([*common, "--no_parent_origin", "--out", str(tmp_path / "legacy")]))
    assert (tmp_path / "legacy.hdf5").is_file()
