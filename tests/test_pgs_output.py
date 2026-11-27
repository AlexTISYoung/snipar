"""The Python PGS command writes one complete estimate table per run."""

import csv
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from snipar import pgs


ROOT = Path(__file__).resolve().parents[1]


def _table(prefix):
    with Path(str(prefix) + '.effects.txt').open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        assert reader.fieldnames == ['statistic', 'est', 'SE']
        rows = list(reader)
    assert len({row['statistic'] for row in rows}) == len(rows)
    return {row['statistic']: row for row in rows}


def _fits():
    rng = np.random.default_rng(71)
    x = rng.normal(size=80)
    parent = .4 * x + rng.normal(size=80)
    response = .7*x + .3*parent + rng.normal(size=80)
    designs = [np.column_stack([np.ones(80), x]),
               np.column_stack([np.ones(80), x, parent])]
    fits = []
    for design, variance, cols in zip(designs, [1.5, 1.],
                                     [['intercept', 'proband'],
                                      ['intercept', 'proband', 'parental']]):
        inverse_design = design / variance
        covariance = np.linalg.inv(design.T @ inverse_design)
        fits.append(dict(alpha=[covariance @ inverse_design.T @ response, covariance],
                         cols=np.array(cols), ids=np.arange(80).astype(str),
                         model=SimpleNamespace(Z=design, Vinv_Z=inverse_design, y=response)))
    return fits


def test_ratio_api_writes_both_models_once_and_preserves_covariances(tmp_path):
    fit1, fit2 = _fits()
    prefix = tmp_path / 'nested' / 'analysis'
    result = pgs.direct_population_ratio(fit1, fit2, outprefix=prefix)
    rows = _table(prefix)
    assert list(rows) == ['1.intercept', '1.proband', '2.intercept',
                          '2.proband', '2.parental', 'direct_population_ratio']
    for generation, fit in [(1, fit1), (2, fit2)]:
        for index, name in enumerate(fit['cols']):
            row = rows['{}.{}'.format(generation, name)]
            np.testing.assert_equal(float(row['est']), fit['alpha'][0][index])
            np.testing.assert_equal(float(row['SE']), np.sqrt(fit['alpha'][1][index,index]))
    np.testing.assert_equal(float(rows['direct_population_ratio']['est']), result['ratio'])
    np.testing.assert_equal(float(rows['direct_population_ratio']['SE']), result['se'])
    joint = np.loadtxt(str(prefix)+'.direct_population.vcov.txt', dtype=str, delimiter='\t')
    assert list(joint[0]) == ['term', 'population', 'direct']
    np.testing.assert_allclose(joint[1:, 1:].astype(float),
                               [[result['population_se']**2, result['covariance']],
                                [result['covariance'], result['direct_se']**2]])
    assert Path(str(prefix)+'.1_2.cross_vcov.txt').is_file()
    for suffix in ('.1.effects.txt', '.2.effects.txt', '.direct_population_ratio.txt'):
        assert not Path(str(prefix)+suffix).exists()


def test_undefined_ratio_is_na_in_unified_table(tmp_path):
    fit1, fit2 = _fits()
    fit1['alpha'][0][1] = 0.
    prefix = tmp_path / 'undefined'
    with pytest.warns(RuntimeWarning, match='undefined_zero_denominator'):
        result = pgs.direct_population_ratio(fit1, fit2, outprefix=prefix)
    assert result['status'] == 'undefined_zero_denominator'
    rows = _table(prefix)
    assert rows['direct_population_ratio']['est'] == 'NA'
    assert rows['direct_population_ratio']['SE'] == 'NA'
    assert float(rows['1.proband']['est']) == 0


def test_adjustment_rows_are_combined_and_missing_se_is_explicit(tmp_path):
    fit = _fits()[1]
    parameters = ['delta', 'alpha', 'rk', 'k', 'r', 'h2_eq', 'rho',
                  'alpha_delta', 'v_eta_delta']
    estimates = dict(zip(parameters, np.arange(9) / 10))
    errors = dict.fromkeys(parameters, .05)
    errors['rho'] = np.nan
    prefix = tmp_path / 'adjusted'
    pgs.write_pgs_effects(prefix, {2: (fit['alpha'], fit['cols'])},
                          adjusted=(estimates, errors))
    rows = _table(prefix)
    assert list(rows)[3:] == ['am_adj.'+name for name in parameters]
    assert rows['am_adj.rho']['SE'] == 'NA'
    assert float(rows['am_adj.rho']['est']) == estimates['rho']
    assert not Path(str(prefix)+'.am_adj_pars.txt').exists()


@pytest.fixture
def inputs(tmp_path):
    rng = np.random.default_rng(681)
    families = np.repeat(np.arange(80), 2)
    n = len(families)
    parental = rng.normal(size=n)
    grandparent1 = rng.normal(size=n)
    grandparent2 = rng.normal(size=n)
    proband = .5*parental + rng.normal(size=n)
    age = rng.normal(size=n)
    y = (.8*proband + .3*parental + .1*grandparent1 + .1*grandparent2
         + .2*age + rng.normal(size=80)[families] + rng.normal(size=n))
    ids = ['id'+str(i) for i in range(n)]
    score_file = tmp_path / 'scores.txt'
    phenotype_file = tmp_path / 'phenotype.txt'
    covar_file = tmp_path / 'covariates.txt'
    with score_file.open('w') as handle:
        writer = csv.writer(handle, delimiter=' ', lineterminator='\n')
        writer.writerow(['FID','IID','proband','parental','gp','gm'])
        writer.writerows(zip(families,ids,proband,parental,grandparent1,grandparent2))
    with phenotype_file.open('w') as handle:
        writer = csv.writer(handle, delimiter=' ', lineterminator='\n')
        writer.writerow(['FID','IID','phenotype'])
        writer.writerows(zip(families,ids,y))
    with covar_file.open('w') as handle:
        writer = csv.writer(handle, delimiter=' ', lineterminator='\n')
        writer.writerow(['FID','IID','age'])
        writer.writerows(zip(families,ids,age))
    return score_file, phenotype_file, covar_file


@pytest.mark.parametrize('models', [None, '1', '2', '3', '1-3', '2-adjusted'])
def test_cli_default_and_selected_generations_share_one_table(tmp_path, inputs, models):
    scores, phenotype, covariates = inputs
    prefix = tmp_path / 'fit'
    command = [sys.executable, str(ROOT/'snipar/scripts/pgs.py'), str(prefix),
               '--pgs',str(scores),'--phenofile',str(phenotype),'--covar',str(covariates),
               '--gparsum']
    adjusted = models == '2-adjusted'
    if adjusted:
        command += ['--gen_models', '2', '--h2f', '.42,.02', '--rk', '.1,.01']
    elif models is not None:
        command += ['--gen_models', models]
    environment = os.environ.copy()
    environment.update(PYTHONPATH=str(ROOT), NUMBA_THREADING_LAYER='workqueue',
                       NUMBA_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    completed = subprocess.run(command, env=environment, cwd=tmp_path,
                               capture_output=True, text=True, timeout=120)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    rows = _table(prefix)
    generations = ([2] if adjusted else [1,2] if models is None else
                   [1,2,3] if models=='1-3' else [int(models)])
    expected = []
    for generation in generations:
        terms = ['intercept','age','proband']
        if generation >= 2:
            terms += ['parental']
        if generation == 3:
            terms += ['gp','gm']
        expected += ['{}.{}'.format(generation, term) for term in terms]
        covariance_path = Path(str(prefix)+'.'+str(generation)+'.vcov.txt')
        covariance = np.loadtxt(covariance_path,dtype=str,skiprows=1)
        assert covariance_path.read_text().splitlines()[0].split() == terms
        for index, term in enumerate(terms):
            actual = rows['{}.{}'.format(generation,term)]
            np.testing.assert_allclose(float(actual['SE'])**2,
                                       float(covariance[index,index+1]),rtol=1e-13)
    if {1,2}.issubset(generations):
        expected += ['direct_population_ratio']
        np.testing.assert_allclose(float(rows['direct_population_ratio']['est']),
                                   float(rows['2.proband']['est'])/float(rows['1.proband']['est']))
        assert Path(str(prefix)+'.1_2.cross_vcov.txt').is_file()
        assert Path(str(prefix)+'.direct_population.vcov.txt').is_file()
    if adjusted:
        expected += ['am_adj.'+name for name in
                     ['delta','alpha','rk','k','r','h2_eq','rho','alpha_delta','v_eta_delta']]
        proband_scores = np.loadtxt(scores, skiprows=1, usecols=2)
        outcomes = np.loadtxt(phenotype, skiprows=1, usecols=2)
        scale = np.std(proband_scores) / np.std(outcomes)
        np.testing.assert_allclose(float(rows['am_adj.alpha']['SE']),
                                   scale * float(rows['2.parental']['SE']))
        assert not Path(str(prefix)+'.am_adj_pars.txt').exists()
    assert list(rows) == expected
    assert len(list(tmp_path.glob('*.effects.txt'))) == 1
    assert not Path(str(prefix)+'.direct_population_ratio.txt').exists()
