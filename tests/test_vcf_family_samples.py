"""Check VCF pruning against the same legacy conversion with every sample."""

import gzip

import h5py
import numpy as np
import pytest

from snipar.scripts import ibd, impute
from snipar.vcf import convert_vcf


def _write_vcf(path, positions, haplotypes, excluded_call=None):
    with path.open("w") as output:
        output.write("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        output.write("\t".join(haplotypes) + "\n")
        for index, position in enumerate(positions):
            calls = (excluded_call if name == "extra" and excluded_call is not None
                     else f"{sample[0, index]}|{sample[1, index]}" for name, sample in haplotypes.items())
            output.write(f"1\t{position}\ts{index}\tA\tG\t.\tPASS\t.\tGT\t")
            output.write("\t".join(calls) + "\n")


def _read_gzip(path):
    with gzip.open(path, "rt") as source:
        return source.read()


@pytest.mark.parametrize("source", ["pedigree", "king"])
def test_vcf_ibd_preserves_sibling_and_error_calibration_pairs(tmp_path, capsys, source):
    rng = np.random.default_rng(116)
    positions = np.arange(101, 141)
    samples = rng.integers(0, 2, size=(8, 2, len(positions)))
    # A separate parent/offspring pair calibrates genotyping errors.
    samples[5, 0] = samples[4, 0]
    names = ("a", "b", "c", "d", "err_parent", "err_child", "unmatched", "extra")
    haplotypes = dict(zip(names, samples))
    vcf = tmp_path / "input.vcf"
    _write_vcf(vcf, positions, haplotypes)
    all_samples = convert_vcf(vcf, tmp_path / "full")
    _write_vcf(vcf, positions, haplotypes, excluded_call="invalid")
    genetic_map = tmp_path / "map.txt"
    genetic_map.write_text("pposition gposition\n" +
                           "".join(f"{p} {i * .03}\n" for i, p in enumerate(positions)))
    if source == "pedigree":
        path = tmp_path / "pedigree.txt"
        path.write_text("FID IID FATHER_ID MOTHER_ID\n"
                        "a a f1 m1\na b f1 m1\nb c f2 m2\nb d f2 m2\n"
                        "err err_child err_parent 0\nmissing unmatched f3 m3\nmissing absent f3 m3\n")
        relatedness = ["--pedigree", str(path)]
    else:
        king, ages = tmp_path / "king.kin0", tmp_path / "ages.txt"
        king.write_text("ID1\tID2\tInfType\na\tb\tFS\nc\td\tFS\n"
                        "err_parent\terr_child\tPO\nunmatched\tabsent\tFS\n")
        ages.write_text("FID IID sex age\n" + "".join(
            f"f {iid} M {55 if iid == 'err_parent' else 25}\n" for iid in names))
        relatedness = ["--king", str(king), "--agesex", str(ages)]
    common = [*relatedness, "--threads", "1", "--map", str(genetic_map), "--ibdmatrix", "--ld_out", "--unphased-ibd"]
    for flag, path, prefix in [("--bed", all_samples.bed, "full"), ("--vcf", vcf, "selected")]:
        # LD estimation draws from the selected sibling rows.
        np.random.seed(116)
        ibd.main(ibd.parser.parse_args([*common, flag, str(path), "--out", str(tmp_path / prefix)]))
    output = capsys.readouterr().out
    assert "40 of 40 VCF records and 6 samples as BED" in output
    assert output.count("1 parent-offspring pairs found") == 2
    for suffix in (".ibd.segments.gz", ".ibdmatrix.gz", ".l2.ldscore.gz"):
        assert _read_gzip(tmp_path / f"full_chr_1{suffix}") == _read_gzip(tmp_path / f"selected_chr_1{suffix}")


def test_vcf_ibd_with_fixed_error_does_not_read_calibration_samples(tmp_path):
    pedigree = tmp_path / "pedigree.txt"
    pedigree.write_text("FID IID FATHER_ID MOTHER_ID\nf a father mother\nf b father mother\n")
    args = ibd.parser.parse_args(["--pedigree", str(pedigree), "--p_error", ".001", "--unphased-ibd"])
    assert ibd._vcf_samples(args, ("a", "b", "father", "mother", "extra")) == {"a", "b"}


@pytest.mark.parametrize("controls", [False, True])
def test_vcf_imputation_preserves_frequency_support_and_full_pedigree(tmp_path, capsys, controls):
    rng = np.random.default_rng(994)
    positions = np.arange(101, 165)
    sample_names = ("c", "p", "f_only", "s1", "s2", "trio", "trio_p", "trio_m", "unrelated", "extra")
    samples = dict(zip(sample_names, rng.integers(0, 2, size=(len(sample_names), 2, len(positions)))))
    samples["c"][0] = samples["p"][0]
    samples["trio"][0] = samples["trio_p"][0]
    samples["trio"][1] = samples["trio_m"][0]
    # Both siblings inherited the same first haplotype, consistent with the segment below.
    samples["s2"][0] = samples["s1"][0]
    vcf = tmp_path / "input.vcf"
    _write_vcf(vcf, positions, samples)
    all_samples = convert_vcf(vcf, tmp_path / "full", output="phased")
    # A call outside the needed set would invalidate all-sample phased conversion.
    _write_vcf(vcf, positions, samples, excluded_call="0/1")
    pedigree = tmp_path / "pedigree.txt"
    pedigree.write_text("FID IID FATHER_ID MOTHER_ID\n"
                        "duo c p absent_m\nfreq f_only absent_p absent_f\n"
                        "siblings s1 absent_sire absent_dam\nsiblings s2 absent_sire absent_dam\n"
                        "complete trio trio_p trio_m\n"
                        "parents trio_p 0 0\nparents trio_m 0 0\n"
                        "u unrelated 0 0\n")
    with gzip.open(tmp_path / "truth.segments.gz", "wt") as output:
        output.write("ID1 ID2 Chr start_coordinate stop_coordinate IBDType\n"
                     "s1 s2 1 101 164 1\n")
    options = ["--pedigree", str(pedigree), "--threads", "1", "--ibd", str(tmp_path / "truth"),
               "--chunks", "3", "-silent_progress"]
    if controls:
        options.append("-c")
    for flag, path, prefix in [("--bgen", all_samples.bgen, "full"), ("--vcf", vcf, "selected")]:
        impute.main(impute.parser.parse_args([*options, flag, str(path), "--out", str(tmp_path / prefix)]))
    assert f"64 of 64 VCF records and {8 if controls else 5} samples as phased BGEN" in capsys.readouterr().out
    with h5py.File(tmp_path / "full.hdf5") as reference, h5py.File(tmp_path / "selected.hdf5") as selected:
        assert set(reference) == set(selected)
        for name in reference:
            if name == "bim_values":
                # BGEN offsets correspond to the temporary files, whose
                # chromosome variant payloads have different byte lengths.
                columns = reference["bim_columns"].asstr()[:]
                keep = columns != "vaddr"
                np.testing.assert_array_equal(selected[name][:, keep], reference[name][:, keep])
                continue
            np.testing.assert_array_equal(selected[name][:], reference[name][:], err_msg=name)
        pedigree_ids = selected["pedigree"].asstr()[1:, 1]
        assert {"trio_p", "trio_m", "unrelated"} <= set(pedigree_ids)
        used = ["c", "p", "f_only", "s1", "s2"]
        if controls:
            used.extend(("trio", "trio_p", "trio_m"))
        # VCF-derived genotypes count the reference allele.
        np.testing.assert_allclose(selected["standard_f"][:], 1 - np.array([samples[i] for i in used]).mean(axis=(0, 1)))


def test_vcf_imputation_panel_maf_uses_all_samples(tmp_path, capsys):
    # The child's genotype and observed parent's genotype are identical.
    # A full-panel allele count retains the variant; the two imputation
    # samples alone would be monomorphic.
    vcf = tmp_path / "panel.vcf"
    _write_vcf(vcf, [101, 102], {"c": np.zeros((2, 2), dtype=int), "p": np.zeros((2, 2), dtype=int),
                                "extra": np.ones((2, 2), dtype=int)})
    pedigree = tmp_path / "pedigree.txt"
    pedigree.write_text("FID IID FATHER_ID MOTHER_ID\nf c p 0\n")
    out = tmp_path / "parents"
    impute.main(impute.parser.parse_args(["--vcf", str(vcf), "--pedigree", str(pedigree),
                                        "--vcf-min-maf", "0.2", "--threads", "1", "--out", str(out)]))
    assert "2 of 2 VCF records and 3 samples" in capsys.readouterr().out
    with h5py.File(out.with_suffix(".hdf5")) as output:
        np.testing.assert_array_equal(output["standard_f"][:], [1., 1.])
