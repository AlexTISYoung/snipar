"""End-to-end checks for the combined family genotype command."""

import gzip
from pathlib import Path

from bed_reader import to_bed
import h5py
import numpy as np
import pytest

from snipar import vcf
from snipar.read.bgen import bgen_is_phased, read_diploid_probabilities
from snipar.scripts import family
from snipar.utilities import open_bgen


def _pedigree(path, rows):
    path.write_text("FID IID FATHER_ID MOTHER_ID\n" + "\n".join(rows) + "\n")
    return path


def _genetic_map(path):
    path.write_text("pposition gposition\n100 0\n200 0.001\n")
    return path


def _write_vcf(path, samples):
    with path.open("w") as stream:
        stream.write("##fileformat=VCFv4.3\n")
        stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        stream.write("\t".join(samples) + "\n")
        for pos, calls in enumerate(zip(*samples.values()), 101):
            stream.write(f"1\t{pos}\ts{pos}\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")
    return path


def _write_bed(path, chromosome, samples):
    calls = np.array(list(samples.values()), dtype=np.float32)
    count = calls.shape[1]
    positions = np.arange(101, count + 101)
    to_bed(path, calls, properties={
        "iid": list(samples), "sid": [f"s{pos}" for pos in positions],
        "bp_position": positions, "chromosome": [str(chromosome)] * count,
    })


def _run(*args):
    family.main(family.parser.parse_args([str(arg) for arg in args]))


def test_combined_phased_vcf_reads_once_and_retains_complete_trio(tmp_path, monkeypatch):
    count = 64
    pedigree = _pedigree(tmp_path / "pedigree.txt", [
        "sibship a father absent", "sibship b father absent", "trio child sire dam",
    ])
    source = _write_vcf(tmp_path / "input.vcf", {
        "a": ["0|1"] * count,
        "b": ["0|1"] * count,
        "father": ["0|0" if i % 4 == 0 else "0|1" for i in range(count)],
        "child": ["1|0"] * count,
        "sire": ["0|0"] * count,
        "dam": ["1|1"] * count,
        "outside_pedigree": ["0/1"] * count,
    })
    genetic_map = _genetic_map(tmp_path / "map.txt")
    convert_vcf = vcf.convert_vcf
    prepared_samples = []

    def record_conversion(*args, **kwargs):
        result = convert_vcf(*args, **kwargs)
        prefix = result.bgen.replace("@", "1")
        with open_bgen(f"{prefix}.bgen", verbose=False) as reader:
            prepared_samples.append(set(reader.samples))
        return result

    monkeypatch.setattr(vcf, "convert_vcf", record_conversion)
    phased = tmp_path / "origin.vcf.gz"
    output = tmp_path / "result"
    _run("--vcf", source, "--pedigree", pedigree, "--map", genetic_map, "--threads", 1,
         "--p_error", ".001", "--out", output, "--ibd", "--impute", "--parent-origin-out", phased)

    assert prepared_samples == [{"a", "b", "father", "child", "sire", "dam"}]
    assert (tmp_path / "result_chr_1.ibd.segments.gz").is_file()
    with h5py.File(tmp_path / "result_chr_1.ibd.posterior.hdf5") as posterior:
        np.testing.assert_array_equal(posterior["pairs"].asstr()[:], [["a", "b"]])
    with h5py.File(output.with_suffix(".hdf5")) as imputed:
        assert np.isfinite(imputed["imputed_par_gts"][:]).any()
    with gzip.open(phased, "rt") as stream:
        for line in stream:
            if line.startswith("#CHROM"):
                export_ids = line.rstrip().split("\t")[9:]
                break
        for line in stream:
            if not line.startswith("#"):
                genotypes = line.rstrip().split("\t")[9:]
                break
    assert set(export_ids) == {"a", "b", "child"}
    assert genotypes[export_ids.index("child")].split(":")[0] == "0|1"


def test_singleton_duo_imputes_and_exports_phased_bgen_without_ibd(tmp_path):
    count = 32
    source = _write_vcf(tmp_path / "input.vcf", {
        "child": ["0|1"] * count, "father": ["1|1"] * count,
    })
    prepared = vcf.convert_vcf(source, tmp_path / "prepared", output="phased")
    pedigree = _pedigree(tmp_path / "pedigree.txt", ["duo child father 0"])
    output = tmp_path / "result"
    phased = tmp_path / "origin.bgen"
    _run("--bgen", f"{prepared.bgen}.bgen", "--pedigree", pedigree, "--map", _genetic_map(tmp_path / "map.txt"),
         "--threads", 1, "--out", output, "--impute", "--parent-origin-out", phased)

    assert not (tmp_path / "result_chr_1.ibd.segments.gz").exists()
    with h5py.File(output.with_suffix(".hdf5")) as imputed:
        assert np.isfinite(imputed["imputed_par_gts"][:]).all()
    with open_bgen(phased) as reader:
        np.testing.assert_array_equal(reader.samples, ["child"])
        assert bgen_is_phased(reader, 0)
        probabilities = read_diploid_probabilities(reader, np.array([0]), np.arange(count), True)
        np.testing.assert_array_equal(probabilities[0], np.tile([0, 1, 1, 0], (count, 1)))


def test_origin_only_export_protects_vcf_input_resolved_from_prefix(tmp_path):
    source = _write_vcf(tmp_path / "input.vcf", {"child": ["1|0"] * 32, "father": ["0|0"] * 32})
    pedigree = _pedigree(tmp_path / "pedigree.txt", ["duo child father 0"])
    original = source.read_bytes()
    common = ("--vcf", source.with_suffix(""), "--pedigree", pedigree, "--threads", 1,
              "--map", _genetic_map(tmp_path / "map.txt"), "--out", tmp_path / "origin")
    with pytest.raises(ValueError, match="different path"):
        _run(*common, "--parent-origin-out", source)
    assert source.read_bytes() == original

    exported = tmp_path / "offspring.vcf"
    _run(*common, "--parent-origin-out", exported)
    assert "\tGT\t0|1\n" in exported.read_text()
    assert not (tmp_path / "origin.hdf5").exists()
    assert not (tmp_path / "origin_chr_1.ibd.segments.gz").exists()


def test_impute_allows_sibships_without_genotyped_parents_or_explicit_error(tmp_path):
    source = _write_vcf(tmp_path / "input.vcf", {
        "a": ["0/1"] * 32, "b": ["0/1"] * 32,
    })
    pedigree = _pedigree(tmp_path / "pedigree.txt", ["f a absent_f absent_m", "f b absent_f absent_m"])
    output = tmp_path / "result"
    _run("--vcf", source, "--pedigree", pedigree, "--map", _genetic_map(tmp_path / "map.txt"),
         "--threads", 1, "--out", output, "--impute")

    with h5py.File(output.with_suffix(".hdf5")) as imputed:
        assert np.isfinite(imputed["imputed_par_gts"][:]).all()
    assert (tmp_path / "result_chr_1.ibd.posterior.hdf5").is_file()


def test_imputation_uses_chromosome_specific_genotyped_pairs(tmp_path):
    prefix = tmp_path / "input_@"
    columns = 32
    _write_bed(tmp_path / "input_1.bed", 1, {
        "a": [1] * columns, "b": [1] * columns, "singleton": [1] * columns, "father": [0] * columns,
    })
    _write_bed(tmp_path / "input_2.bed", 2, {
        "b": [1] * columns, "singleton": [1] * columns, "father": [0] * columns,
    })
    pedigree = _pedigree(tmp_path / "pedigree.txt", [
        "siblings a missing_f missing_m", "siblings b missing_f missing_m", "duo singleton father 0",
    ])
    output = tmp_path / "result_@"
    _run("--bed", prefix, "--chr_range", 1, 2, "--pedigree", pedigree,
         "--map", _genetic_map(tmp_path / "map.txt"), "--threads", 1, "--out", output, "--impute")

    assert (tmp_path / "result_1.ibd.segments.gz").is_file()
    assert not (tmp_path / "result_2.ibd.segments.gz").exists()
    for chromosome in [1, 2]:
        with h5py.File(Path(str(output).replace("@", str(chromosome)) + ".hdf5")) as imputed:
            assert "duo" in imputed["families"].asstr()[:]
            assert np.isfinite(imputed["imputed_par_gts"][:]).any()


def test_imputation_accepts_existing_ibd_input(tmp_path):
    source = _write_vcf(tmp_path / "input.vcf", {"a": ["0/1"] * 32, "b": ["0/1"] * 32})
    pedigree = _pedigree(tmp_path / "pedigree.txt", ["f a absent_f absent_m", "f b absent_f absent_m"])
    supplied = tmp_path / "supplied"
    with gzip.open(str(supplied) + ".segments.gz", "wt") as output:
        output.write("ID1 ID2 Chr start_coordinate stop_coordinate IBDType\na b 1 101 132 0\n")
    output = tmp_path / "result"
    _run("--vcf", source, "--pedigree", pedigree, "--threads", 1, "--out", output,
         "--impute", "--ibd-input", supplied)
    assert not (tmp_path / "result_chr_1.ibd.segments.gz").exists()
    with h5py.File(output.with_suffix(".hdf5")) as imputed:
        np.testing.assert_array_equal(imputed["imputed_par_gts"][0], np.ones(32))
