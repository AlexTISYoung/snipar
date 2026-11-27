"""Check the installed BGEN backend against the packaged PLINK genotypes."""

from pathlib import Path
import shutil

import numpy as np
import pytest
from pysnptools.snpreader import Bed

from snipar import example
from snipar.utilities import open_bgen, read_oxford_sample


DATA = Path(example.__file__).parent / "example_data"


def test_phased_bgen_agrees_with_bed(tmp_path):
    # Give backends that build a metadata index a writable directory.
    filename = tmp_path / "chr_1.bgen"
    shutil.copyfile(DATA / filename.name, filename)
    reader = open_bgen(str(filename))
    bed = Bed(str(DATA / "chr_1"), count_A1=True)
    try:
        np.testing.assert_array_equal(reader.samples, bed.iid[:, 1])
        np.testing.assert_array_equal(reader.rsids, bed.sid)
        np.testing.assert_array_equal(reader.positions, bed.pos[:, 2])
        samples = [7, 29, 19, 7, 2, 1, 0]
        variants = [8, 5, 5, 2, 0]
        probs = reader.read((samples, variants), dtype=np.float32)
        assert probs.dtype == np.float32
        assert probs.shape == (len(samples), len(variants), 4)
        expected = bed[samples, variants].read().val
        dosage = probs[:, :, [0, 2]].sum(axis=2)
        # The shipped BED deliberately masks some known phased genotypes.
        observed = np.isfinite(expected)
        np.testing.assert_array_equal(dosage[observed], expected[observed])
    finally:
        reader.close()


def test_native_bgen_metadata():
    pytest.importorskip("bgen")
    from snipar.bgen_io import BgenFile

    with BgenFile(DATA / "chr_1.bgen") as reader:
        metadata = reader.variant_metadata()
        assert list(metadata.columns) == ["id", "rsid", "chrom", "pos", "nalleles", "allele_ids", "vaddr"]
        assert metadata.iloc[0].tolist() == ["1", "1_rs0", "", 1000, 2, "A,G", 61952]
        np.testing.assert_array_equal(metadata.rsid, reader.rsids)
        np.testing.assert_array_equal(metadata.pos, reader.positions)
        np.testing.assert_array_equal(metadata.allele_ids, reader.allele_ids)
        np.testing.assert_array_equal(reader.read((slice(0, 3), [3, 0])),
                                      reader.read((np.array([0, 1, 2]), [0, 3]))[:, ::-1])


def test_oxford_sample_preserves_individual_ids(tmp_path):
    filename = tmp_path / "genotypes.sample"
    filename.write_text("ID_1 ID_2 missing\n0 0 0\nfamily 0001 0\nfamily NA 0\n")
    np.testing.assert_array_equal(read_oxford_sample(filename), ["0001", "NA"])
