"""Integration checks for exact blockwise duo inference and legacy phased mode."""

import h5py
import numpy as np

from snipar.scripts import impute


def test_phased_imputation_across_a_phase_switch(tmp_path):
    pedigree = tmp_path / "pedigree.txt"
    pedigree.write_text("FID IID FATHER_ID MOTHER_ID\nf c p 0\n")
    vcf = tmp_path / "offspring.vcf"
    positions = np.arange(101, 293)
    switch = len(positions) // 2
    header = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tp\tc\n"
    records = []
    for marker, position in enumerate(positions):
        parent = "0|0" if marker % 4 == 0 else "0|1"
        child = "0|1" if marker < switch else "1|0"
        records.append(f"1\t{position}\ts{marker}\tA\tG\t.\tPASS\t.\tGT\t{parent}\t{child}\n")
    vcf.write_text(header + "".join(records))
    genetic_map = tmp_path / "map.txt"
    genetic_map.write_text("pposition gposition\n100 0\n300 0.001\n")
    common = ["--vcf", str(vcf), "--pedigree", str(pedigree), "--threads", "1",
              "--parent_origin_map", str(genetic_map), "-silent_progress"]
    for name, options in (("whole", []), ("blocks", ["--chunks", "4"]),
                          ("legacy", ["--no_parent_origin", "--chunks", "4"])):
        impute.main(impute.parser.parse_args([*common, *options, "--out", str(tmp_path / name)]))

    perm = np.random.default_rng(42).permutation(len(positions))
    vcf.write_text(header + "".join(records[index] for index in perm))
    impute.main(impute.parser.parse_args([*common, "--chunks", "4", "--out", str(tmp_path / "shuffled")]))

    with (h5py.File(tmp_path / "whole.hdf5") as whole, h5py.File(tmp_path / "blocks.hdf5") as blocks,
          h5py.File(tmp_path / "legacy.hdf5") as legacy, h5py.File(tmp_path / "shuffled.hdf5") as shuffled):
        assert set(whole) == set(blocks) == set(legacy)
        for key in whole:
            np.testing.assert_array_equal(whole[key][:], blocks[key][:], err_msg=key)
        # The VCF converter counts REF alleles. The missing parent transmitted ALT,
        # and the other allele's expectation is the sample frequency of REF (0.5).
        near_switch = [switch - 1, switch + 1]
        np.testing.assert_array_equal(whole["imputed_par_gts"][0, near_switch], [.5, .5])
        np.testing.assert_array_equal(legacy["imputed_par_gts"][0, near_switch], [1., 1.])
        np.testing.assert_array_equal(shuffled["pos"][:], positions[perm])
        np.testing.assert_array_equal(shuffled["imputed_par_gts"][:], whole["imputed_par_gts"][:][:, perm])
    assert not list(tmp_path.glob(".snipar-impute-*"))
