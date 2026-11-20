"""Genotype LD input preserves SNP alignment and chromosome output conventions."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from snipar import coherent_ld


def _bed(tmp_path, name, records):
    prefix = tmp_path / name
    prefix.with_suffix('.bed').write_bytes(b'')
    prefix.with_suffix('.fam').write_text('family person 0 0 0 -9\n')
    prefix.with_suffix('.bim').write_text(''.join(
        '{} {} {} {} A G\n'.format(chrom, snp, i / 10, bp)
        for i, (chrom, snp, bp) in enumerate(records)))
    return prefix


@pytest.mark.parametrize('output_pattern', ['ld_', 'ld_chr@', 'ld_chr#'])
@pytest.mark.parametrize('bed_pattern', ['chr@', 'chr#'])
def test_chromosome_alignment_and_ldsc_exports(tmp_path, monkeypatch, output_pattern, bed_pattern):
    _bed(tmp_path, 'chr1', [(1, 'b', 20), (1, 'a', 10)])
    _bed(tmp_path, 'chr2', [(2, 'a', 30), (2, 'c', 40)])
    calls = []

    def compute(filename, chromosome, window, output):
        calls.append((Path(filename).name, chromosome, window, output))
        return ((np.array([2., 3.]), np.array(['a', 'b'])) if chromosome == 1 else
                (np.array([7., 11.]), np.array(['c', 'a'])))

    monkeypatch.setattr(coherent_ld, 'ldscores_from_bed', compute)
    values, metadata = coherent_ld.compute_ldscores(
        str(tmp_path / bed_pattern), ['c', 'a', 'a', 'b', 'absent', 'missing_chr'],
        [2, 1, 2, 1, 1, 3], ld_wind=.75, ld_out=tmp_path / output_pattern)
    np.testing.assert_allclose(values, [7., 2., 11., 3., np.nan, np.nan], equal_nan=True)
    assert calls == [('chr1.bed', 1, .75, None), ('chr2.bed', 2, .75, None)]
    assert len(metadata['bed_files']) == 2
    assert len(metadata['ldscore_files']) == 2
    output = pd.read_csv(metadata['ldscore_files'][0], sep='\t')
    assert list(output.columns) == ['CHR', 'SNP', 'BP', 'L2']
    assert output.to_dict('list') == {'CHR': [1, 1], 'SNP': ['a', 'b'], 'BP': [10, 20], 'L2': [2., 3.]}
    output_pattern = output_pattern.replace('#', '@')
    expected = output_pattern.replace('@', '1') if '@' in output_pattern else output_pattern + '1'
    assert Path(metadata['ldscore_files'][0]).name == expected + '.l2.ldscore.gz'


def test_single_bed_infers_chromosome_from_bim(tmp_path, monkeypatch):
    prefix = _bed(tmp_path, 'single', [(5, 'a', 10), (5, 'b', 20)])
    calls = []

    def compute(filename, chromosome, window, output):
        calls.append(chromosome)
        return np.array([4., 6.]), np.array(['b', 'a'])

    monkeypatch.setattr(coherent_ld, 'ldscores_from_bed', compute)
    values, _ = coherent_ld.compute_ldscores(prefix.with_suffix('.bed'), ['a', 'b'], [5, 5])
    np.testing.assert_allclose(values, [6., 4.])
    assert calls == [5]


def test_unknown_summary_chromosomes_align_by_unique_snp(tmp_path, monkeypatch):
    _bed(tmp_path, 'chr1', [(1, 'a', 10), (1, 'b', 20)])
    _bed(tmp_path, 'chr2', [(2, 'c', 30), (2, 'd', 40)])

    def compute(filename, chromosome, window, output):
        return ((np.array([2., 3.]), np.array(['b', 'a'])) if chromosome == 1 else
                (np.array([7., 11.]), np.array(['d', 'c'])))

    monkeypatch.setattr(coherent_ld, 'ldscores_from_bed', compute)
    values, _ = coherent_ld.compute_ldscores(tmp_path / 'chr@', ['c', 'b', 'absent', 'a', 'd'], None)
    np.testing.assert_allclose(values, [11., 2., np.nan, 3., 7.], equal_nan=True)
    _bed(tmp_path, 'chr2', [(2, 'a', 30), (2, 'b', 40)])
    monkeypatch.setattr(coherent_ld, 'ldscores_from_bed',
                        lambda *args: (np.array([2., 3.]), np.array(['a', 'b'])))
    with pytest.raises(ValueError, match='unique across BED files'):
        coherent_ld.compute_ldscores(tmp_path / 'chr@', ['a', 'b'], None)


def test_chromosome_mismatch_and_mixed_bim_are_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(coherent_ld, 'ldscores_from_bed', lambda *args: pytest.fail('backend must not run'))
    _bed(tmp_path, 'chr1', [(2, 'a', 10), (2, 'b', 20)])
    with pytest.raises(ValueError, match='does not match'):
        coherent_ld.compute_ldscores(tmp_path / 'chr@', ['a'], [1])
    prefix = _bed(tmp_path, 'mixed', [(1, 'a', 10), (2, 'b', 20)])
    with pytest.raises(ValueError, match='one chromosome per file'):
        coherent_ld.compute_ldscores(prefix, ['a', 'b'], [1, 2])


def test_duplicate_backend_ids_and_no_overlap_are_rejected(tmp_path, monkeypatch):
    prefix = _bed(tmp_path, 'single', [(1, 'a', 10), (1, 'b', 20)])
    monkeypatch.setattr(coherent_ld, 'ldscores_from_bed',
                        lambda *args: (np.array([1., 2.]), np.array(['a', 'a'])))
    with pytest.raises(ValueError, match='duplicate SNP IDs'):
        coherent_ld.compute_ldscores(prefix, ['a'], [1])
    monkeypatch.setattr(coherent_ld, 'ldscores_from_bed',
                        lambda *args: (np.array([1.]), np.array(['b'])))
    with pytest.raises(ValueError, match='No SNPs in common'):
        coherent_ld.compute_ldscores(prefix, ['a'], [1])
