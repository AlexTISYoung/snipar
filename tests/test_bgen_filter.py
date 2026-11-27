from pathlib import Path
import shutil

import numpy as np
import pytest

from snipar import example, utilities
from snipar.bgen_filter import BGEN_SUBSET_SUFFIX, _selected_rows, filter_bgen
from snipar.config import nan_integer
from snipar.read.bgen import _first_allele_dosage, bgen_is_phased, read_hard_calls
from snipar.read.variants import bgen_variant_ids


DATA = Path(example.__file__).parent / 'example_data'


@pytest.fixture
def source(tmp_path):
    for filename in ('chr_1.bgen', 'chr_1.sample'):
        shutil.copyfile(DATA / filename, tmp_path / filename)
    return tmp_path / 'chr_1.bgen'


@pytest.mark.parametrize('backend', ['legacy', 'native'])
def test_bgen_selection_reuses_source_probabilities_and_sample_order(source, tmp_path, monkeypatch, backend):
    if backend == 'legacy':
        pytest.importorskip('bgen_reader')
    else:
        pytest.importorskip('bgen')
        from snipar.bgen_io import BgenFile

        monkeypatch.setattr(utilities, '_open_bgen', lambda filename, verbose: BgenFile(filename))

    target = tmp_path / ('filtered' + BGEN_SUBSET_SUFFIX)
    result = filter_bgen(source, target, {'1_rs90', '1_rs3', '1_rs5'}, variant_ids={'1_rs90', '1_rs3'})
    assert result.path == target
    assert (result.n_variants, result.n_kept) == (1000, 2)
    assert target.stat().st_size < 1024
    assert target.with_suffix('.npy').stat().st_size < 1024
    samples = np.array([29, 0, 4, 0])
    variants = np.array([1, 0, 1])

    with utilities.open_bgen(source) as original, utilities.open_bgen(target) as selected:
        np.testing.assert_array_equal(selected.samples, original.samples)
        np.testing.assert_array_equal(selected.ids, original.ids[[3, 90]])
        np.testing.assert_array_equal(selected.positions, original.positions[[3, 90]])
        np.testing.assert_array_equal(selected.allele_ids, original.allele_ids[[3, 90]])
        np.testing.assert_array_equal(bgen_variant_ids(selected)[0], ['1_rs3', '1_rs90'])
        assert bgen_is_phased(selected, 0)
        if backend == 'legacy':
            assert selected.shape == (len(original.samples), 2, 4)
        else:
            np.testing.assert_array_equal(selected.variant_metadata()['rsid'], selected.rsids)
        actual = selected.read((samples, variants), np.float32)
        expected = original.read((samples, [90, 3, 90]), np.float32)
        np.testing.assert_array_equal(actual, expected)
        dose = _first_allele_dosage(selected, samples, variants)
        calls = read_hard_calls(selected, samples, variants, True, nan_integer)
        np.testing.assert_array_equal(dose, expected[..., 0] + expected[..., 2])
        np.testing.assert_array_equal(calls, dose)


def test_bgen_single_variant_keeps_source_analysis_id(source, tmp_path):
    target = tmp_path / ('one' + BGEN_SUBSET_SUFFIX)
    filter_bgen(source, target, {'1_rs10'})
    with utilities.open_bgen(target) as selected:
        ids, aliases = bgen_variant_ids(selected)
        np.testing.assert_array_equal(ids, ['1_rs10'])
        assert aliases is None


def test_bgen_quality_primary_id_precedes_rsid():
    ids = ['s1', 's2', 's3']
    aliases = ['r1', 'r2', 'r3']
    selected = _selected_rows(ids, aliases, {'r1', 'r2', 's3'}, {'s1', 'r1', 'r2', 's3'}, None)
    np.testing.assert_array_equal(selected, [1, 2])


def test_bgen_quality_rejects_an_ambiguous_rsid_match():
    with pytest.raises(ValueError, match="RSID 'r1' matches multiple variants"):
        _selected_rows(['s1', 's2'], ['r1', 'r1'], {'s1', 'r1'}, None, None)


def test_bgen_quality_empty_and_all_kept_needs_no_mapping(source, tmp_path):
    none = tmp_path / ('none' + BGEN_SUBSET_SUFFIX)
    with pytest.raises(ValueError, match='No variants passed'):
        filter_bgen(source, none, {'1'})
    assert not none.exists()
    assert not none.with_suffix('.npy').exists()

    target = tmp_path / ('all' + BGEN_SUBSET_SUFFIX)
    result = filter_bgen(source, target, {f'1_rs{i}' for i in range(1000)})
    assert result.n_kept == result.n_variants == 1000
    assert not target.with_suffix('.npy').exists()
    with utilities.open_bgen(target) as selected:
        assert len(selected.ids) == 1000


def test_bgen_selection_preserves_an_existing_destination(source, tmp_path):
    target = tmp_path / ('existing' + BGEN_SUBSET_SUFFIX)
    target.write_text('original')
    with pytest.raises(FileExistsError):
        filter_bgen(source, target, {'1_rs10'})
    assert target.read_text() == 'original'
    assert not target.with_suffix('.npy').exists()
