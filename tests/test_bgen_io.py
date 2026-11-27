"""Check the installed BGEN backend against the packaged PLINK genotypes."""

import os
from pathlib import Path
import shutil
import subprocess
import sys
from unittest.mock import Mock

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


@pytest.mark.skipif(os.name != "posix", reason="test requires POSIX directory permissions")
@pytest.mark.parametrize("stale", [False, True], ids=["missing", "stale"])
def test_legacy_indexes_in_readonly_directory(tmp_path, stale):
    pytest.importorskip("bgen_reader")
    source = (DATA / "chr_1.bgen").read_bytes()
    cache = tmp_path / "cache"
    paths = []
    for label in ["a", "b"]:
        directory = tmp_path / label
        directory.mkdir()
        filename = directory / "chr.bgen"
        filename.write_bytes(source.replace(b"1_rs0", f"{label}_rs0".encode(), 1))
        if stale:
            for suffix in [".metafile", ".metadata2.mmm"]:
                index = directory / (filename.name + suffix)
                index.write_bytes(b"stale")
                os.utime(index, (1, 1))
        directory.chmod(0o555)
        paths.append(filename)

    try:
        if os.access(paths[0].parent, os.W_OK):
            pytest.skip("this account can write read-only directories")
        environment = dict(os.environ, BGEN_READER_CACHE_HOME=str(cache))
        code = """
import sys
from pathlib import Path
from snipar.utilities import open_bgen, read_bgen

source = Path(sys.argv[1])
with open_bgen(str(source)) as reader:
    assert reader.rsids[0] == source.parent.name + '_rs0'
assert read_bgen(str(source))['variants'].compute().rsid.iloc[0] == source.parent.name + '_rs0'
"""
        # Exercise two workers on one input as well as an identically named
        # chromosome from a second dataset, using the same initially empty cache.
        workers = [subprocess.Popen([sys.executable, "-c", code, str(source)], env=environment)
                   for source in [paths[0], paths[0], paths[1]]]
        assert [worker.wait(timeout=120) for worker in workers] == [0, 0, 0]
        assert len(list((cache / "metafile").glob("*.metadata2.mmm"))) == 2
        assert len(list((cache / "metafile").glob("*.metafile"))) == 2
        for filename in paths:
            for suffix in [".metafile", ".metadata2.mmm"]:
                sidecar = filename.with_name(filename.name + suffix)
                if stale:
                    assert sidecar.read_bytes() == b"stale"
                else:
                    assert not sidecar.exists()
    finally:
        for filename in paths:
            filename.parent.chmod(0o755)


def test_legacy_writable_indexes_refresh_and_reuse(tmp_path):
    pytest.importorskip("bgen_reader")
    from snipar.utilities import read_bgen

    filename = tmp_path / "chr.bgen"
    shutil.copyfile(DATA / "chr_1.bgen", filename)
    indexes = [tmp_path / "chr.bgen.metadata2.mmm", tmp_path / "chr.bgen.metafile"]
    def read_and_check(expected):
        with open_bgen(str(filename)) as reader:
            assert reader.rsids[0] == expected
        assert read_bgen(str(filename))["variants"].compute().rsid.iloc[0] == expected
        return [index.stat().st_mtime_ns for index in indexes]

    initial = read_and_check("1_rs0")
    assert read_and_check("1_rs0") == initial
    filename.write_bytes(filename.read_bytes().replace(b"1_rs0", b"2_rs0", 1))
    new_time = max(initial) + 2_000_000_000
    os.utime(filename, ns=(new_time, new_time))
    assert min(read_and_check("2_rs0")) >= new_time


def test_legacy_external_sample_validated_by_both_readers(tmp_path):
    pytest.importorskip("bgen_reader")
    from snipar.utilities import _bgen, read_bgen

    # Omit the BGEN sample-ID block to exercise the Oxford .sample fallback.
    data = bytearray((DATA / "chr_1.bgen").read_bytes())
    offset = int.from_bytes(data[:4], "little")
    header_size = int.from_bytes(data[4:8], "little")
    data[header_size + 3] &= 0x7F
    data[:4] = header_size.to_bytes(4, "little")
    filename = tmp_path / "chr.bgen"
    filename.write_bytes(data[:header_size + 4] + data[offset + 4:])
    sample = filename.with_suffix(".sample")
    sample.write_text((DATA / "chr_1.sample").read_text().replace("0_0 0_0", "family 00001", 1))
    with open_bgen(str(filename)) as reader:
        assert reader.samples[0] == "00001"
    assert read_bgen(str(filename))["samples"][0] == "00001"
    sample.write_text("ID_1 ID_2 missing\n0 0 0\nfamily 00001 0\n")
    for read in [open_bgen, read_bgen]:
        with pytest.raises(ValueError, match="sample file length"):
            read(str(filename))
    opened = Mock(nsamples=5)
    with pytest.raises(FileNotFoundError):
        _bgen(filename, tmp_path / "missing.sample", reader=opened)
    opened.close.assert_called_once_with()


def test_legacy_rebuilds_unreadable_metadata_with_current_mtime(tmp_path):
    pytest.importorskip('bgen_reader')
    filename = tmp_path / 'chr.bgen'
    shutil.copyfile(DATA / 'chr_1.bgen', filename)
    index = filename.with_name(filename.name + '.metadata2.mmm')
    index.write_bytes(b'partial metadata mmap')
    newer = filename.stat().st_mtime_ns + 2_000_000_000
    os.utime(index, ns=(newer, newer))
    with open_bgen(filename) as reader:
        assert reader.rsids[0] == '1_rs0'
        assert reader.read(([0], [0])).shape == (1, 1, 4)
    assert index.stat().st_size > 1000
