"""Round-trip paternal-first, maternal-second haplotype exports."""

import gzip
import struct

import numpy as np
import pytest

from snipar.config import nan_integer
from snipar.imputation.origin_export import export_parental_haplotypes
from snipar.imputation.origin_io import OriginWriter
from snipar.utilities import open_bgen


_PEDIGREE = [['fam', 'son', 'dad', 'mom'], ['fam2', 'daughter', '0', 'other-mom'],
             ['fam3', 'unrelated', '0', '0']]
_MARKERS = ['m1', 'm2', 'm3', 'm4', 'm5']
_POSITIONS = [100, 120, 140, 160, 180]
_ALLELES = [['A', 'G'], ['T', 'C'], ['A', 'G'], ['A', 'G'], ['A', 'G']]


def _write_origin(path, pairs, calls, origin, chromosome=1):
    with OriginWriter(path, chromosome, pairs, _POSITIONS, _MARKERS, _ALLELES) as writer:
        writer.write(slice(None), slice(None), np.asarray(origin, dtype=np.int8),
                     np.full((len(pairs), len(_MARKERS)), -1, dtype=np.int8),
                     calls, np.arange(2 * len(pairs)).reshape(-1, 2))


@pytest.fixture
def sidecars(tmp_path):
    son = np.array([[1, 0], [1, 0], [0, 0], [1, nan_integer], [0, 1]], dtype=np.int8)
    daughter = np.array([[1, 0], [0, 1], [1, 1], [1, 0], [0, 1]], dtype=np.int8)
    path1 = tmp_path / 'batch1.hdf5'
    path2 = tmp_path / 'batch2.hdf5'
    _write_origin(path1, [('other-mom', 'daughter'), ('dad', 'son')],
                  np.stack([np.ones_like(daughter), daughter, np.ones_like(son), son]),
                  [[0, -1, -1, 1, 0], [0, 0, -1, 0, -1]])
    _write_origin(path2, [('mom', 'son')], np.stack([np.ones_like(son), son]), [[1, 0, -1, -1, 0]])
    return path1, path2


_EXPECTED = [['0|1', '1|0'], ['.|.', '.|.'], ['1|1', '0|0'], ['.|.', '0|1'], ['0|1', '0|1']]


@pytest.mark.parametrize('suffix', ['.vcf', '.vcf.gz', '.bgen'])
def test_origin_export_formats(tmp_path, sidecars, suffix):
    path = tmp_path / ('ordered.haps' + suffix)
    result = export_parental_haplotypes(path, 'chr1', _PEDIGREE, sidecars)
    assert result == path
    if suffix == '.bgen':
        with path.open('rb') as f:
            _, size, _, _, magic = struct.unpack('<4I4s', f.read(20))
            assert magic == b'bgen'
            assert f.read(size - 20).endswith(b'paternal|maternal')
        with open_bgen(path) as reader:
            np.testing.assert_array_equal(reader.samples, ['son', 'daughter'])
            np.testing.assert_array_equal(reader.positions, _POSITIONS)
            np.testing.assert_array_equal(reader.rsids, _MARKERS)
            np.testing.assert_array_equal(reader.allele_ids, ['A,G', 'T,C', 'A,G', 'A,G', 'A,G'])
            probabilities = reader.read().transpose(1, 0, 2)
            for marker, row in enumerate(_EXPECTED):
                for sample, gt in enumerate(row):
                    probs = probabilities[marker, sample].reshape(2, 2)
                    if gt == '.|.':
                        assert np.isnan(probs).all()
                    else:
                        np.testing.assert_array_equal(probs.argmax(axis=1), [int(gt[0]), int(gt[2])])
        return
    opener = gzip.open if suffix == '.vcf.gz' else open
    with opener(path, 'rt') as f:
        lines = f.read().splitlines()
    assert '##SNIPar_haplotype_order=paternal|maternal' in lines
    headers = [line.split('\t') for line in lines if not line.startswith('##')]
    assert headers[0][-2:] == ['son', 'daughter']
    assert [line[3:5] for line in headers[1:]] == _ALLELES
    assert [line[-2:] for line in headers[1:]] == _EXPECTED


def test_export_rejects_mismatched_haplotypes_without_replacing_output(tmp_path, sidecars):
    path = tmp_path / 'existing.vcf'
    path.write_text('previous output\n')
    conflict = tmp_path / 'rephased.hdf5'
    swapped = np.ones((2, len(_MARKERS), 2), dtype=np.int8)
    swapped[1, 0] = [0, 1]
    _write_origin(conflict, [('mom', 'son')], swapped, [[1, -1, -1, -1, -1]])
    with pytest.raises(ValueError, match='disagree on original offspring haplotypes'):
        export_parental_haplotypes(path, 1, _PEDIGREE, [sidecars[0], conflict])
    assert path.read_text() == 'previous output\n'
