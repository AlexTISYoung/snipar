import subprocess
import sys

from bed_reader import to_bed
import h5py
import numpy as np
import pandas as pd
import pytest


def _run(module, *args):
    result = subprocess.run([sys.executable, '-m', f'snipar.scripts.{module}', *map(str, args)],
                            text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def _write_observed(path, kind, ids, variants, values, invalid_id=None):
    if kind == 'bed':
        to_bed(path.with_suffix('.bed'), values.astype(np.float32), count_A1=True,
               properties=dict(iid=ids, sid=variants, chromosome=['1'] * len(variants),
                               bp_position=np.arange(1, len(variants) + 1),
                               allele_1=['A'] * len(variants), allele_2=['G'] * len(variants)))
        return '--bed', path
    calls = np.array(['1/1', '0/1', '0/0'])[values].astype('<U3')
    if kind == 'bgen':
        calls = np.char.replace(calls, '/', '|')
    elif invalid_id is not None:
        calls[ids.index(invalid_id)] = 'bad'
    lines = ['##fileformat=VCFv4.2',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *ids])]
    for position, (variant, genotypes) in enumerate(zip(variants, calls.T), 1):
        lines.append('\t'.join(['1', str(position), variant, 'A', 'G', '.', 'PASS', '.', 'GT', *genotypes]))
    vcf = path.with_suffix('.vcf')
    vcf.write_text('\n'.join(lines) + '\n')
    if kind == 'bgen':
        from snipar.vcf import convert_vcf

        return '--bgen', convert_vcf(vcf, path.parent / 'prepared', output='phased').bgen
    return '--vcf', vcf


@pytest.fixture
def sibling_cohort(tmp_path):
    ids = ['001', 'sib1', 'dad1', 'mom1', '002', 'sib2', 'dad2', 'mom2',
           '003', 'sib3', 'dad3', 'mom3', '004', 'sib4', 'dad4', 'mom4', 'unrelated']
    calls = np.array([
        [0, 2, 1], [1, 1, 1], [0, 1, 2], [2, 2, 0],
        [1, 0, 0], [2, 1, 0], [1, 0, 1], [2, 1, 0],
        [2, 0, 0], [1, 1, 0], [1, 1, 0], [2, 0, 1],
        [0, 1, 2], [1, 2, 1], [2, 2, 1], [0, 1, 1],
        [2, 0, 0],
    ])
    ped = tmp_path / 'pedigree.txt'
    ped.write_text('FID IID FATHER_ID MOTHER_ID\n' + ''.join(
        f'f{i} 00{i} dad{i} mom{i}\nf{i} sib{i} dad{i} mom{i}\n' for i in range(1, 5)))
    keep = tmp_path / 'keep.txt'
    keep.write_text('003\n001\n003\nmissing\n')
    weights = tmp_path / 'weights.txt'
    weights.write_text('SNP A1 A2 b\na A G 0.2\nb A G -0.6\nc A G 0.3\n')
    return ids, calls, ped, keep, weights


@pytest.mark.parametrize('kind', ['bed', 'bgen', 'vcf'])
def test_keep_scores_and_matrix_read_relatives_outside_selected_ids(sibling_cohort, tmp_path, kind):
    ids, calls, ped, keep, weights = sibling_cohort
    observed = _write_observed(tmp_path / 'observed', kind, ids, ['a', 'b', 'c'], calls,
                               invalid_id='unrelated')
    shared = [*observed, '--pedigree', ped, '--keep', keep, '--fit_sib', '--threads', '1']
    score_out = tmp_path / 'scores'
    score_stdout = _run('pgs', score_out, *shared, '--weights', weights, '--no_am_adj', '--batch_size', '1')
    matrix_out = tmp_path / 'matrix'
    matrix_stdout = _run('gts_matrix', '--out', matrix_out, *shared)
    if kind == 'vcf':
        assert 'and 8 samples as BED' in score_stdout
        assert 'and 8 samples as BED' in matrix_stdout
    expected = calls[np.array([[0, 1, 2, 3], [8, 9, 10, 11]])]
    with h5py.File(matrix_out.with_suffix('.genoarray.hdf5')) as matrix:
        np.testing.assert_array_equal(matrix['ids'].asstr()[:], ['001', '003'])
        np.testing.assert_array_equal(matrix['genotype_cols'].asstr()[:], ['proband', 'sib', 'paternal', 'maternal'])
        np.testing.assert_array_equal(matrix['genotypes'][:], expected)
    scores = pd.read_csv(score_out.with_suffix('.pgs.txt'), sep=r'\s+', dtype={'IID': str})
    assert scores.IID.tolist() == ['001', '003']
    allelic_weights = np.array([.2, -.6, .3])
    centered = expected - expected[:, :1].mean(axis=0)
    np.testing.assert_allclose(scores[['proband', 'sibling', 'paternal', 'maternal']],
                               centered @ allelic_weights, atol=2e-7)


def test_phenotype_alone_retains_full_score_file(sibling_cohort, tmp_path):
    ids, calls, ped, _, weights = sibling_cohort
    observed = _write_observed(tmp_path / 'observed', 'bed', ids, ['a', 'b', 'c'], calls)
    phenotypes = tmp_path / 'phenotype.txt'
    phenotypes.write_text('IID outcome\n001 1.0\n002 -1.0\n003 0.2\n004 0.7\n')
    output = tmp_path / 'uncropped'
    _run('pgs', output, *observed, '--weights', weights, '--pedigree', ped,
         '--phenofile', phenotypes, '--gen_models', '1', '--no_am_adj', '--threads', '1')
    scores = pd.read_csv(output.with_suffix('.pgs.txt'), sep=r'\s+', dtype={'IID': str})
    assert scores.IID.tolist() == [iid for iid in ids if iid.startswith(('00', 'sib'))]
    assert output.with_suffix('.effects.txt').is_file()


def test_keep_subset_uses_precomputed_scores_for_regression(tmp_path):
    scores = tmp_path / 'scores.txt'
    phenotypes = tmp_path / 'phenotype.txt'
    keep = tmp_path / 'keep.txt'
    dose = np.arange(10)
    y = np.r_[3 * dose[:6] + [0.2, -0.2, 0.3, 0.1, -0.25, 0.2], -3 * dose[6:]]
    ids = [f'{i:03}' for i in range(10)]
    scores.write_text('FID IID proband\n' + ''.join(
        f'f{i} {iid} {dose[i]}\n' for i, iid in enumerate(ids)))
    phenotypes.write_text('IID outcome\n' + ''.join(f'{iid} {value}\n' for iid, value in zip(ids, y)))
    keep.write_text('\n'.join(ids[:6]) + '\n')
    output = tmp_path / 'selected'
    _run('pgs', output, '--pgs', scores, '--phenofile', phenotypes, '--keep', keep,
         '--gen_models', '1', '--threads', '1')
    effects = pd.read_csv(output.with_suffix('.effects.txt'), sep=r'\s+').set_index('statistic')
    slope = np.linalg.lstsq(np.column_stack((np.ones(6), dose[:6])), y[:6], rcond=None)[0][1]
    np.testing.assert_allclose(effects.loc['1.proband', 'est'], slope, atol=1e-7)
    assert not output.with_suffix('.pgs.txt').exists()


@pytest.mark.parametrize('kind', ['bed', 'vcf'])
def test_grandparent_score_uses_auxiliary_parents_without_emitting_them(tmp_path, kind):
    ids = ['001', 'dad', 'mom', 'paternal_gp', 'paternal_gm', 'maternal_gp', 'maternal_gm', 'unrelated']
    values = np.array([[1], [1], [2], [2], [0], [2], [1], [0]])
    observed = _write_observed(tmp_path / 'observed', kind, ids, ['a'], values,
                               invalid_id='unrelated')
    ped = tmp_path / 'pedigree.txt'
    ped.write_text('FID IID FATHER_ID MOTHER_ID\nf0 001 dad mom\n'
                   'f1 dad paternal_gp paternal_gm\nf2 mom maternal_gp maternal_gm\n')
    weights = tmp_path / 'weights.txt'
    weights.write_text('SNP A1 A2 b\na A G 1\n')
    keep = tmp_path / 'keep.txt'
    keep.write_text('001\n')
    output = tmp_path / 'grandparents'
    stdout = _run('pgs', output, *observed, '--weights', weights, '--pedigree', ped, '--keep', keep,
                  '--grandpar', '--no_am_adj', '--threads', '1')
    if kind == 'vcf':
        assert 'and 7 samples as BED' in stdout
    scores = pd.read_csv(output.with_suffix('.pgs.txt'), sep=r'\s+', dtype={'IID': str})
    assert scores.IID.tolist() == ['001']
    np.testing.assert_allclose(scores[['gpp', 'gmp', 'gpm', 'gmm']],
                               np.array([[2., 0., 2., 1.]]) - 4 / 3, atol=2e-7)


def test_keep_reads_imputed_families_and_keeps_observed_relatives(sibling_cohort, tmp_path):
    ids, calls, _, keep, weights = sibling_cohort
    observed_rows = ~np.isin(ids, ['dad1', 'dad3'])
    observed = _write_observed(tmp_path / 'observed', 'vcf', list(np.array(ids)[observed_rows]),
                               ['a', 'b', 'c'], calls[observed_rows], invalid_id='unrelated')
    imp = tmp_path / 'parents.hdf5'
    with h5py.File(imp, 'w') as hdf:
        hdf['pedigree'] = np.array([
            ['FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'has_father', 'has_mother'],
            ['f1', '001', 'dad1', 'mom1', 'False', 'True'],
            ['f1', 'sib1', 'dad1', 'mom1', 'False', 'True'],
            ['f3', '003', 'dad3', 'mom3', 'False', 'True'],
            ['f3', 'sib3', 'dad3', 'mom3', 'False', 'True'],
        ], dtype='S')
        hdf['families'] = np.array(['f1', 'f3'], dtype='S')
        hdf['imputed_par_gts'] = np.array([[.5, 1.5, .75], [2., .25, .5]])
        hdf['bim_columns'] = np.array(['id', 'allele1', 'allele2'], dtype='S')
        hdf['bim_values'] = np.array([['a', 'A', 'G'], ['b', 'A', 'G'], ['c', 'A', 'G']], dtype='S')
    args = [*observed, '--imp', imp.with_suffix(''), '--keep', keep, '--fit_sib', '--threads', '1']
    scores = tmp_path / 'scores'
    stdout = _run('pgs', scores, *args, '--weights', weights, '--no_am_adj')
    assert 'and 6 samples as BED' in stdout
    matrix_out = tmp_path / 'matrix'
    _run('gts_matrix', '--out', matrix_out, *args)
    expected = calls[np.array([[0, 1, 2, 3], [8, 9, 10, 11]])].astype(float)
    expected[:, 2] = [[.5, 1.5, .75], [2., .25, .5]]
    with h5py.File(matrix_out.with_suffix('.genoarray.hdf5')) as matrix:
        np.testing.assert_array_equal(matrix['ids'].asstr()[:], ['001', '003'])
        np.testing.assert_allclose(matrix['genotypes'][:], expected)
    score_table = pd.read_csv(scores.with_suffix('.pgs.txt'), sep=r'\s+', dtype={'IID': str})
    assert score_table.IID.tolist() == ['001', '003']
    centered = expected - expected[:, :1].mean(axis=0)
    np.testing.assert_allclose(score_table[['proband', 'sibling', 'paternal', 'maternal']],
                               centered @ np.array([.2, -.6, .3]), atol=2e-7)
