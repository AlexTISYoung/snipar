"""Check the installed BGEN backend against the packaged PLINK genotypes."""

import gzip
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
from unittest.mock import Mock
import zlib

import h5py
import numpy as np
import pandas as pd
import pytest
from pysnptools.snpreader import Bed, SnpData

from snipar import example
from snipar.config import nan_integer
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


def _write_bgen(filename, samples, variants, *, compression=1, bits=8):
    """Write Layout 2 rows with integer probabilities."""
    def string(value, length="H"):
        field = value.encode("ascii")
        return struct.pack("<" + length, len(field)) + field

    ids = b"".join(string(sample) for sample in samples)
    sample_block = struct.pack("<II", len(ids) + 8, len(samples)) + ids
    if compression == 2:
        if sys.version_info >= (3, 14):
            from compression.zstd import compress
        else:
            compress = pytest.importorskip('zstandard').compress
    else:
        compress = zlib.compress

    flags = (1 << 31) | (2 << 2) | compression
    with filename.open("wb") as bgen:
        bgen.write(struct.pack("<IIII4sI", 20 + len(sample_block), 20, len(variants), len(samples), b"bgen", flags))
        bgen.write(sample_block)
        for index, (phased, ploidy, probabilities) in enumerate(variants, 1):
            bgen.write(string(f"id{index}") + string(f"rs{index}") + string("1"))
            bgen.write(struct.pack("<IH", index * 100, 2) + string("A", "I") + string("G", "I"))
            copy_numbers = [value & 63 for value in ploidy]
            payload = struct.pack("<IHBB", len(samples), 2, min(copy_numbers), max(copy_numbers))
            payload += np.asarray(ploidy, dtype=np.uint8).tobytes() + bytes((phased, bits))
            if bits < 8:
                words = np.asarray(probabilities, dtype=np.uint8)
                padding = -len(words) % (8 // bits)
                values = np.pad(words, (0, padding)).reshape(-1, 8 // bits)
                payload += np.bitwise_or.reduce(values << np.arange(0, 8, bits, dtype=np.uint8), axis=1).tobytes()
            else:
                payload += np.asarray(probabilities, dtype=f'<u{bits // 8}').tobytes()
            if compression:
                block = compress(payload)
                bgen.write(struct.pack("<II", len(block) + 4, len(payload)) + block)
            else:
                bgen.write(struct.pack("<I", len(payload)) + payload)


@pytest.mark.parametrize('compression,bits,phased', [
    (0, 1, False), (1, 1, True), (1, 2, False), (1, 4, True),
    (1, 8, False), (1, 8, True), (1, 16, False), (2, 8, False),
    (2, 1, True), (2, 8, True), (2, 16, False),
])
@pytest.mark.parametrize('backend', ['legacy', 'native'])
def test_packed_bgen_subset_matches_full_decoder(tmp_path, monkeypatch, compression, bits, phased, backend):
    from snipar import bgen_packed
    from snipar.read.bgen import _first_allele_dosage, read_diploid_probabilities, read_hard_calls

    if backend == 'native':
        library = pytest.importorskip('bgen')
        from snipar.bgen_io import BgenFile
        opener = BgenFile
        packed_expected = compression != 1 and (bits == 1 or (
            bits == 8 and library.__version__ == '1.11.0' and (not phased or bgen_packed._NATIVE_PHASED_8)
        ))
    else:
        if compression == 0:
            pytest.skip('cbgen cannot index the uncompressed Layout 2 fixture')
        library = pytest.importorskip('bgen_reader')

        def opener(path):
            return library.open_bgen(path, verbose=False)

        packed_expected = True
    filename = tmp_path / 'packed.bgen'
    # Native decoding is faster on small files. Use a partial final byte for both readers.
    count = 32_769 if backend == 'native' else 2049
    ids = [f'sample{i}' for i in range(count)]
    rng = np.random.default_rng(2026)
    maximum = (1 << bits) - 1
    records = []
    for row in range(4):
        first = rng.integers(maximum + 1, size=count)
        second = rng.integers(maximum + 1 if phased else maximum + 1 - first, size=count)
        missing = (np.arange(count) + row) % 7 == 0
        states = np.column_stack((first, second))
        states[missing] = 0
        records.append((phased, np.where(missing, 130, 2), states.ravel()))
    _write_bgen(filename, ids, records, compression=compression, bits=bits)
    original_read = bgen_packed.PackedBgenSubset.read
    used = []

    def track_read(self, offset):
        result = original_read(self, offset)
        used.append(result is not None)
        return result

    monkeypatch.setattr(bgen_packed.PackedBgenSubset, 'read', track_read)
    variants = [3, 0, 2, 2, -3]
    samples = np.array([-1, 19, 0, 1, 6, 7, 2, 9, 8, 0, -7])
    boolean_mask = np.isin(np.arange(count), [0, 1, 6, 7, count - 1])
    with opener(str(filename)) as reader:
        probabilities = read_diploid_probabilities(reader, None, variants, phased)
        dosage = _first_allele_dosage(reader, None, variants)
        calls = read_hard_calls(reader, None, variants, phased, nan_integer)
        assert not used
        for selected in (samples, boolean_mask, slice(50, 3, -2)):
            np.testing.assert_array_equal(
                read_diploid_probabilities(reader, selected, variants, phased, num_threads=2),
                probabilities[selected],
            )
            np.testing.assert_array_equal(_first_allele_dosage(reader, selected, variants, num_threads=2),
                                          dosage[selected])
            np.testing.assert_array_equal(read_hard_calls(reader, selected, variants, phased, nan_integer,
                                                          num_threads=2), calls[selected])
    if compression == 2 and bits == 8 and not phased:
        from snipar.bgen_filter import BGEN_SUBSET_SUFFIX, filter_bgen, open_bgen_subset

        target = tmp_path / ('view' + BGEN_SUBSET_SUFFIX)
        filter_bgen(filename, target, {'id1', 'id4'})
        with open_bgen_subset(target, lambda path, verbose=False: opener(path)) as reader:
            np.testing.assert_array_equal(_first_allele_dosage(reader, samples, [1, 0, 1]),
                                          dosage[samples][:, [0, 1, 0]])
    assert any(used) == packed_expected
    assert all(value == packed_expected for value in used)


@pytest.fixture
def unphased_bgen(tmp_path):
    filename = tmp_path / "calls.bgen"
    ids = np.array(["child0", "father0", "child1", "father1"])
    # Each sample stores P(AA), P(AG), and implies P(GG) from one minus their sum.
    _write_bgen(filename, ids, [
        (0, [2, 2, 2, 2], [0, 254, 255, 0, 0, 0, 0, 255]),
        (0, [2, 2, 130, 2], [50, 153, 0, 255, 0, 0, 0, 0]),
        (0, [2, 2, 2, 2], [255, 0, 0, 255, 0, 255, 255, 0]),
    ])
    expected_calls = np.array([[1, nan_integer, 2], [2, 1, 1],
                               [0, nan_integer, 1], [1, 0, 2]], dtype=np.int8)
    expected_doses = expected_calls.astype(np.float32)
    expected_doses[0, 0] = 254 / 255
    expected_doses[0, 1] = 253 / 255
    expected_doses[2, 1] = np.nan
    return filename, ids, expected_calls, expected_doses


@pytest.fixture(params=["default", "native"])
def bgen_backend(request, monkeypatch):
    if request.param == "native":
        pytest.importorskip("bgen")
        from snipar import observed
        from snipar.bgen_io import BgenFile
        from snipar.imputation import preprocess_data

        def open_native(filename, verbose=False):
            return BgenFile(filename)

        monkeypatch.setattr(observed, "open_bgen", open_native)
        monkeypatch.setattr(preprocess_data, "open_bgen", open_native)


def test_unphased_bgen_expected_dosage_and_confident_calls(unphased_bgen, bgen_backend):
    from snipar.observed import open_observed

    filename, _, calls, doses = unphased_bgen
    samples, variants = [3, 0, 2, 2], [2, 0, 1, 0]
    with open_observed(bgenfile=filename, num_threads=2) as reader:
        np.testing.assert_allclose(reader.read(samples, variants), doses[np.ix_(samples, variants)], rtol=1e-6)
        np.testing.assert_array_equal(reader.read(samples, variants, require_calls=True), calls[np.ix_(samples, variants)])


@pytest.mark.parametrize("layout", ["mixed phase", "triploid"])
def test_hard_calls_validate_each_bgen_variant(tmp_path, bgen_backend, layout):
    from snipar.observed import open_observed

    filename = tmp_path / "layout.bgen"
    plain = (0, [2, 2], [0, 255, 0, 255])
    if layout == "mixed phase":
        # The legacy metadata probes rows 0, 1, 3, and misses the phase of row 2.
        variants = [plain, plain, (1, [2, 2], [255, 0, 0, 255]), plain]
        error = "consistently phased or unphased"
    else:
        variants = [(0, [3, 2], [0, 0, 255, 0, 255])]
        error = "diploid, biallelic"
    _write_bgen(filename, ["a", "b"], variants)
    with open_observed(bgenfile=filename, num_threads=2) as reader:
        with pytest.raises(ValueError, match=error):
            reader.read(np.arange(2), np.arange(len(variants)), require_calls=True)
        if layout == "triploid":
            with pytest.raises(ValueError, match=error):
                reader.read(np.arange(2), np.arange(len(variants)))
        else:
            np.testing.assert_array_equal(reader.read(np.arange(2), np.arange(len(variants))), 1)


def test_unphased_bgen_imputation_matches_bed(unphased_bgen, bgen_backend, tmp_path):
    from snipar.scripts import impute

    filename, ids, calls, _ = unphased_bgen
    fids = np.repeat(["family0", "family1"], 2)
    pedigree = tmp_path / "pedigree.txt"
    pd.DataFrame({"FID": fids[::2], "IID": ids[::2], "FATHER_ID": ids[1::2], "MOTHER_ID": "0"}).to_csv(
        pedigree, sep="\t", index=False)
    bed = tmp_path / "plink"
    genotypes = np.where(calls == nan_integer, np.nan, calls)
    Bed.write(str(bed), SnpData(iid=np.column_stack((fids, ids)), sid=["rs1", "rs2", "rs3"],
                               pos=[[1, 0, 100], [1, 0, 200], [1, 0, 300]], val=genotypes), count_A1=True)
    args = ["--pedigree", str(pedigree), "--threads", "2", "--chunks", "2", "-silent_progress"]
    for mode, source in (("bgen", filename.with_suffix("")), ("bed", bed)):
        impute.main(impute.parser.parse_args([*args, "--" + mode, str(source), "--out", str(tmp_path / mode)]))
    with h5py.File(tmp_path / "bgen.hdf5") as bgen, h5py.File(tmp_path / "bed.hdf5") as bed:
        for key in ("families", "standard_f", "imputed_par_gts", "num_observed_parental_alleles"):
            np.testing.assert_array_equal(bgen[key][:], bed[key][:], err_msg=key)


def test_unphased_bgen_ibd_matches_bed(unphased_bgen, bgen_backend, tmp_path):
    from snipar.ibd import infer_ibd_chr

    filename, ids, calls, _ = unphased_bgen
    bed = tmp_path / "plink"
    genotypes = np.where(calls == nan_integer, np.nan, calls)
    Bed.write(str(bed), SnpData(iid=np.column_stack((ids, ids)), sid=["id1", "id2", "id3"],
                               pos=[[1, 0, 100], [1, 0, 200], [1, 0, 300]], val=genotypes), count_A1=True)
    genetic_map = tmp_path / "map.txt"
    genetic_map.write_text("pposition gposition\n100 .1\n200 .2\n300 .3\n")
    pairs = ids.reshape(-1, 2)
    state = np.random.get_state()
    try:
        for label, source in (("bgen", {"bgenfile": filename}), ("bed", {"bedfile": bed.with_suffix(".bed")})):
            np.random.seed(17)
            infer_ibd_chr(pairs, .001, None, str(tmp_path / label), chrom=1, mapfile=genetic_map,
                          max_missing=50, ibdmatrix=True, num_threads=2, **source)
    finally:
        np.random.set_state(state)
    for suffix in (".ibd.segments.gz", ".ibdmatrix.gz"):
        with gzip.open(tmp_path / ("bgen_chr_1" + suffix), "rt") as bgen, gzip.open(tmp_path / ("bed_chr_1" + suffix), "rt") as bed:
            assert bgen.read() == bed.read()


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
