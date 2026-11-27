import gc
import gzip
from pathlib import Path
import shutil
import subprocess
import sys
import weakref

import h5py
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pandas as pd
import pytest
from scipy.sparse import csc_matrix

from snipar import example, gwas
from snipar.gtarray import gtarray, impute_missing, impute_unrel_par_gts
from snipar.slmm import LinearMixedModel, build_grm_arr, build_ibdrel_arr, match_grm_ids


DATA = Path(example.__file__).parent / 'example_data'


def _dense_gls(y, covariates, covariance, genotypes):
    design = np.column_stack((genotypes, covariates))
    precision_design = np.linalg.solve(covariance, design)
    coefficient_covariance = np.linalg.inv(design.T @ precision_design)
    weights = coefficient_covariance @ precision_design.T
    count = genotypes.shape[1]
    return (weights[:count] @ y, coefficient_covariance[:count, :count], weights[0])


def _sparse_component(matrix):
    rows, cols = np.tril_indices(len(matrix))
    values = matrix[rows, cols]
    keep = values != 0
    return [(values[keep], rows[keep], cols[keep])]


def test_gls_covariates_and_intercept_match_joint_dense_regression():
    rng = np.random.default_rng(192)
    n = 34
    a = rng.normal(size=(n, n))
    relationship = a @ a.T / n
    residual = 0.7
    covariance = relationship + residual * np.eye(n)
    gts = rng.normal(size=(n, 2, 6)).astype(np.float32)
    covariate = gts[:, 0, 1] + rng.normal(size=n)
    y = 2 + covariate + rng.normal(size=n)
    model = LinearMixedModel(y, _sparse_component(relationship), (1, residual),
                             covar_X=np.column_stack((covariate, covariate)), add_intercept=True)
    gts[:, 1, 2] = gts[:, 0, 2]  # Unidentifiable parental and direct coefficients.
    gts[:, 1, 4] = 1  # The constant predictor is itself collinear with the intercept.
    gts[1, 1, 5] = np.nan
    alpha, alpha_cov, ses = model.fit_snps_eff(gts)
    for i in (0, 1, 3):
        beta, beta_cov, _ = _dense_gls(y, np.column_stack((np.ones(n), covariate)), covariance,
                                      gts[:, :, i])
        assert_allclose(alpha[i], beta, rtol=2e-12)
        assert_allclose(alpha_cov[i], beta_cov, rtol=2e-12)
        assert_allclose(ses[i], np.sqrt(np.diag(beta_cov)), rtol=2e-12)
    assert np.isnan(alpha[[2, 4, 5]]).all()
    assert np.isnan(alpha_cov[[2, 4, 5]]).all()


def test_trios_sibs_cross_covariance_uses_full_genotype_weights():
    rng = np.random.default_rng(212)
    n = 28
    a = rng.normal(size=(n, n))
    relationship = a @ a.T / n
    cov = relationship + np.eye(n)
    y = rng.normal(size=n)
    Z = rng.normal(size=(n, 1))
    gts = rng.normal(size=(n, 3, 4))
    trios = np.arange(0, n, 2)
    sibs = np.arange(1, n, 2)
    gts[sibs, 1, 2] = gts[sibs, 0, 2]
    model = LinearMixedModel(y, _sparse_component(relationship), (1, 1),
                             covar_X=Z, add_intercept=True)
    beta, variance, ses = model.trios_sibs_est(gts, trios, sibs)
    for i in range(gts.shape[2]):
        groups = [(trios, 3)] if i == 2 else [(trios, 3), (sibs, 2)]
        estimates = []
        weights = np.zeros((n, len(groups)))
        for j, (rows, predictors) in enumerate(groups):
            coefficients, _, first_weights = _dense_gls(
                y[rows], np.column_stack((np.ones(len(rows)), Z[rows])),
                cov[np.ix_(rows, rows)], gts[rows, :predictors, i])
            estimates.append(coefficients[0])
            weights[rows, j] = first_weights
        combination = np.linalg.solve(weights.T @ cov @ weights, np.ones(len(groups)))
        assert_allclose(beta[i, 0], combination @ estimates / combination.sum(), rtol=1e-12)
        assert_allclose(variance[i, 0, 0], 1 / combination.sum(), rtol=1e-12)
        assert_allclose(ses[i, 0] ** 2, variance[i, 0, 0], rtol=1e-12)


@pytest.mark.parametrize('covariates', [False, True])
def test_reml_methods_gradient_and_cache_lifetime_agree(covariates):
    rng = np.random.default_rng(210)
    n = 12
    a = rng.normal(size=(n, n))
    relationship = a @ a.T / n
    z = rng.normal(size=n) if covariates else None
    model = LinearMixedModel(rng.normal(size=n), _sparse_component(relationship), (1, 1),
                             covar_X=z)
    likelihood = model.reml_loglik
    factor = model.V_lu
    assert factor is model.V_lu
    for method in ('chol', 'cg'):
        assert_allclose(model.dense_reml_loglik(method), likelihood, atol=1e-11)
    gradient = model.grad.copy()
    step = 1e-5
    for i in range(2):
        perturbed = np.ones(2)
        perturbed[i] += step
        model._varcomps = tuple(perturbed)
        upper = model.reml_loglik
        perturbed[i] -= 2 * step
        model._varcomps = tuple(perturbed)
        lower = model.reml_loglik
        assert_allclose((upper - lower) / (2 * step), gradient[i], atol=1e-8)
    assert factor is not model.V_lu
    reference = weakref.ref(model)
    del model
    gc.collect()
    assert reference() is None


def test_sparse_dense_solver_handles_single_predictor_and_float32():
    v = csc_matrix([[2., 0.25], [0.25, 3.]])
    x = np.arange(6, dtype=np.float32).reshape(3, 2, 1)
    expected = [np.linalg.solve(v.toarray(), snp) for snp in x]
    assert_allclose(LinearMixedModel.sp_solve_dense3d(v, x), expected)
    assert_allclose(LinearMixedModel.sp_solve_dense3d_lu(v, x), expected)


def test_grid_reml_skips_singular_covariance_boundaries():
    from snipar.slmm import build_sib_arr

    sibling = build_sib_arr(np.repeat(['f1', 'f2', 'f3'], 2))
    model = LinearMixedModel(np.array([0., 0.2, 2., 2.1, 4.2, 4.]), [sibling])
    model.grid_search_reml()
    assert model.varcomps[1] > 0
    assert np.isfinite(model.reml_loglik)


def test_king_and_gcta_grm_keep_leading_zeros_and_sample_order(tmp_path):
    ids = np.array(['010', '007', '003'])
    pairs = tmp_path / 'pairs.seg'
    pairs.write_text('ID1  ID2  relatedness\n010  003  0.2\n007  003  0.01\n')
    values, row, col = build_ibdrel_arr(pairs, dict(zip(ids, range(3))), ids, thres=0.05)
    assert_allclose(values, [0.2, 1, 1, 1])
    assert_array_equal(row, [2, 0, 1, 2])
    assert_array_equal(col, [0, 0, 1, 2])
    prefix = tmp_path / 'sample'
    Path(f'{prefix}.grm.id').write_text('f1 003\nf2 010\n')
    with gzip.open(f'{prefix}.grm.gz', 'wt') as f:
        f.write('1 1 100 1.02\n2 1 100 0.25\n2 2 100 0.98\n')
    matched_ids, families = match_grm_ids(ids, np.array(['a', 'b', 'c']), prefix, 'gcta')
    assert_array_equal(matched_ids, ['010', '003'])
    assert_array_equal(families, ['a', 'c'])
    values, row, col = build_grm_arr(prefix, {'010': 0, '003': 1}, 0.05)
    assert_allclose(values, [1.02, 0.25, 0.98])
    assert_array_equal(row, [1, 1, 0])
    assert_array_equal(col, [1, 0, 0])


def test_genotype_filters_preserve_metadata_and_allow_zero_missingness(capsys):
    data = gtarray(np.ma.array([[0., 1.], [2., 1.]], mask=[[False, False], [False, True]]),
                   ids=np.array(['a', 'b']), sid=np.array(['ok', 'missing']))
    data.compute_info()
    assert_array_equal(data.filter_missingness(0, verbose=True), [True, False])
    assert '1 SNPs' in capsys.readouterr().out
    assert data.info.shape == data.variances.shape == (1,)
    assert not data.has_NAs
    data.filter_ids(['b'])
    data.compute_freqs()
    assert_allclose(data.freqs, [1.0])
    no_metadata = gtarray(np.ma.array([[1.], [3.]]), np.array(['a', 'b']))
    assert_array_equal(no_metadata.filter_missingness(0), [True])
    assert_array_equal(no_metadata.fill_NAs(), [0])
    assert_allclose(no_metadata.gts.data.ravel(), [-1, 1])
    with pytest.raises(ValueError, match='Shape of SNP ids'):
        gtarray(np.zeros((2, 2)), np.array(['a', 'b']), sid=np.array(['a']))


def test_missing_parent_imputation_uses_valid_allele_dosages():
    separate = gtarray(np.array([[[1.], [1.], [np.nan]], [[0.], [1.], [1.]]]),
                       np.array(['a', 'b']))
    impute_missing(separate)
    assert_allclose(separate.gts.data[0, 2], [0.75])
    summed = gtarray(np.array([[[1.], [2.]], [[np.nan], [np.nan]]]), np.array(['a', 'b']))
    impute_missing(summed)
    assert_allclose(summed.gts.data[1, :, 0], [1, 2])
    assert not summed.has_NAs


def test_conditional_gaussian_independent_proband_matches_linear_imputation():
    ids = np.array(['a', 'b'])
    genotypes = np.array([[[0., 2.], [-1., -1.], [-1., -1.]],
                          [[2., 0.], [1., 0.], [1., 2.]]])
    expected = impute_unrel_par_gts(gtarray(genotypes.copy(), ids)).gts.data
    unrelated_indices = np.array([0])
    result = impute_unrel_par_gts(gtarray(genotypes.copy(), ids),
                                 ped=np.array([['f', 'a', '0', '0'], ['g', 'b', '0', '0']]),
                                 grm=(np.ones(2), np.arange(2), np.arange(2)),
                                 unrelated_inds=unrelated_indices)
    assert_allclose(result.gts.data, expected)


def test_text_sumstats_uses_parental_average_and_stable_normal_tail(tmp_path):
    assert gwas.outarray_effect(np.array([50.]), np.array([1.]), np.array([0.5]), 1)[0, 4] > 500
    output = tmp_path / 'summary.gz'
    gwas.write_txt_output(np.array([1]), np.array(['s']), np.array([100]), np.array([['A', 'C']]),
                          str(output), False, True, False, np.array([[1., 2., 3., 5.]]),
                          np.eye(4)[None, :, :], [1.], np.array([0.5]))
    summary = pd.read_csv(output, sep=r'\s+').iloc[0]
    assert summary['avg_NTC_Beta'] == 4.
    assert summary['population_Beta'] == 6.
    assert_allclose(summary['avg_NTC_SE'], np.sqrt(0.5), atol=1e-6)
    assert_allclose(summary['population_SE'], np.sqrt(1.75), atol=1e-6)


def test_genotype_matrix_cli_writes_each_chromosome_and_applies_maf(tmp_path):
    for suffix in (1, 2):
        shutil.copyfile(DATA / 'chr_1_trios.bgen', tmp_path / f'input_{suffix}.bgen')
    command = [sys.executable, '-m', 'snipar.scripts.gts_matrix', '--bgen', str(tmp_path / 'input_@'),
               '--chr_range', '1-2', '--pedigree', str(DATA / 'pedigree.txt'), '--threads', '1',
               '--min_maf', '0.4', '--max_missing', '0', '--out', str(tmp_path / 'output_@')]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    for chrom in (1, 2):
        with h5py.File(tmp_path / f'output_{chrom}.genoarray.hdf5') as f:
            assert_array_equal(f['genotype_cols'].asstr()[:], ['proband', 'paternal', 'maternal'])
            calls = f['genotypes'][:]
            assert 0 < calls.shape[2] < 1000
            frequency = calls[:, 0].mean(axis=0, dtype=float) / 2
            assert np.all((frequency > 0.4) & (frequency < 0.6))


def test_gwas_keep_integer_phenotype_and_gcta_id_subset(tmp_path):
    original = pd.read_csv(DATA / 'phenotype.txt', sep=r'\s+')
    eligible = [f'{family}_{child}' for family in range(10) for child in (0, 1)]
    phenotype = original[original['IID'].isin(eligible)].copy()
    phenotype['phenotype'] = np.arange(len(phenotype))
    phenofile = tmp_path / 'integer.tsv'
    phenotype.to_csv(phenofile, sep='\t', index=False)
    keep = tmp_path / 'keep'
    keep.write_text('\n'.join(eligible[:-2]) + '\n')
    grm_prefix = tmp_path / 'input'
    Path(f'{grm_prefix}.grm.id').write_text(''.join(f'f {iid}\n' for iid in eligible[6:]))
    with gzip.open(f'{grm_prefix}.grm.gz', 'wt') as f:
        for i in range(1, len(eligible[6:]) + 1):
            f.write(f'{i} {i} 100 1\n')
    out = tmp_path / 'gwas'
    command = [sys.executable, '-m', 'snipar.scripts.gwas', str(phenofile),
               '--bed', str(DATA / 'chr_1'), '--pedigree', str(DATA / 'pedigree.txt'),
               '--sib_diff', '--vc_list', '0.25', '0.25', '0.5', '--keep', str(keep),
               '--grmgz', str(grm_prefix), '--out', str(out)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'Using 18 individuals eligible on all chromosomes.' in result.stdout
    with h5py.File(out.with_suffix('.sumstats.hdf5')) as f:
        assert np.isfinite(f['estimate'][:]).any()


def test_robust_shared_covariance_matches_dense_group_regressions():
    rng = np.random.default_rng(193)
    n = 40
    ids = np.array([f'id{i}' for i in range(n)])
    snps = np.array(['a', 'b', 'c'])
    a = rng.normal(size=(n, n))
    kinship = a @ a.T / n
    covariance = kinship + np.eye(n)
    phenotype = rng.normal(size=n)
    covariate = rng.normal(size=n)
    gts = rng.uniform(0, 2, (n, 3, 3))
    sib_gts = rng.uniform(0, 2, (10, 2, 3))
    status = np.repeat([[0, 0], [1, 1], [0, 1], [1, 0]], 10, axis=0)
    alleles = np.full((n, 3), 3)
    alleles[:10, :] = 4
    alleles[20:30, 2] = 2  # Missing group makes this design unavailable.
    model = LinearMixedModel(phenotype, _sparse_component(kinship), (1, 1), covar_X=covariate,
                             add_intercept=True)
    g = gtarray(gts, ids, sid=snps, par_status=status, num_obs_par_al=alleles)
    sib = gtarray(sib_gts[::-1], ids[10:20][::-1], sid=snps)
    estimates, variances, _ = model.robust_est(g, sib)
    for locus in (0, 1):
        paternal = gts[20:30, :, locus].copy()
        maternal = gts[30:40, :, locus].copy()
        paternal[:, 2] = paternal[:, 0] - np.floor(paternal[:, 2])
        maternal[:, 1] = maternal[:, 0] - np.floor(maternal[:, 1])
        designs = (gts[:10, :, locus], sib_gts[:, :, locus], paternal, maternal)
        weights = np.zeros((n, 4))
        effects = []
        for i, design in enumerate(designs):
            rows = slice(10 * i, 10 * (i + 1))
            coef, _, w = _dense_gls(phenotype[rows],
                                     np.column_stack((np.ones(10), covariate[rows])),
                                     covariance[rows, rows], design)
            weights[rows, i] = w
            effects.append(coef[0])
        combination = np.linalg.solve(weights.T @ covariance @ weights, np.ones(4))
        assert_allclose(estimates[locus, 0], combination @ effects / combination.sum(), rtol=1e-12)
        assert_allclose(variances[locus, 0, 0], 1 / combination.sum(), rtol=1e-12)
    assert np.isnan(estimates[2, 0])


def test_gwas_aligns_shared_cohort_across_chromosomes(tmp_path):
    from pysnptools.snpreader import Bed, SnpData

    rng = np.random.default_rng(194)
    ids = np.array([f'child{i:02}' for i in range(20)])
    families = np.repeat([f'family{i:02}' for i in range(10)], 2)
    variants = np.array(['s1', 's2', 's3'])
    matrix = rng.integers(0, 3, size=(len(ids), len(variants))).astype(float)
    for chromosome, indices in ((1, np.arange(20)), (2, np.arange(17, -1, -1))):
        pos = np.column_stack((np.full(3, chromosome), np.zeros(3), np.arange(1, 4)))
        bed = SnpData(iid=np.column_stack((families, ids))[indices], sid=variants,
                      val=matrix[indices], pos=pos)
        Bed.write(str(tmp_path / f'input_{chromosome}.bed'), bed, count_A1=True)
    pedigree = tmp_path / 'pedigree.txt'
    rows = [[family, iid, f'{family}_father', f'{family}_mother'] for family, iid in zip(families, ids)]
    np.savetxt(pedigree, rows, fmt='%s', header='FID IID FATHER_ID MOTHER_ID', comments='')
    phenotype = tmp_path / 'phenotype.txt'
    phenotype.write_text('IID phenotype\n' + ''.join(f'{iid} {value}\n' for iid, value in
                                                  zip(ids, rng.normal(size=len(ids)))))
    command = [sys.executable, '-m', 'snipar.scripts.gwas', str(phenotype),
               '--bed', str(tmp_path / 'input_@'), '--chr_range', '1-2',
               '--pedigree', str(pedigree), '--sib_diff',
               '--vc_list', '0.5', '0.5', '--out', str(tmp_path / 'output_@')]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'Using 18 individuals eligible on all chromosomes.' in result.stdout
    with h5py.File(tmp_path / 'output_1.sumstats.hdf5') as left, h5py.File(tmp_path / 'output_2.sumstats.hdf5') as right:
        assert np.isfinite(left['estimate'][:]).all()
        for field in ('estimate', 'estimate_ses', 'estimate_covariance', 'freqs'):
            assert_allclose(left[field][:], right[field][:], rtol=1e-7)


def test_genotype_matrix_uses_imputed_pedigree_and_reports_actual_parental_sum(tmp_path):
    variants = np.arange(10)
    bim = np.loadtxt(DATA / 'chr_1.bim', dtype=str)[variants]
    pedigree = np.loadtxt(DATA / 'pedigree.txt', dtype=str)
    families = np.arange(16).astype(str)
    children = np.char.endswith(pedigree[:, 1], '_0') | np.char.endswith(pedigree[:, 1], '_1')
    rows = pedigree[np.isin(pedigree[:, 0], families) & children]
    imputation = tmp_path / 'parents.hdf5'
    with h5py.File(imputation, 'w') as f:
        f['families'] = families.astype('S')
        f['pedigree'] = np.vstack(([*pedigree[0], 'has_father', 'has_mother'],
                                   np.column_stack((rows, np.full((len(rows), 2), 'False'))))).astype('S')
        f['bim_columns'] = np.array(['Chr', 'id', 'morgans', 'coordinate', 'allele1', 'allele2'], dtype='S')
        f['bim_values'] = bim.astype('S')
        f['imputed_par_gts'] = np.full((len(families), len(variants)), 0.75)
    from snipar.scripts.gts_matrix import gts_matrix_chromosome

    gts_matrix_chromosome(0, None, str(tmp_path / 'export'),
                          bedfile=str(DATA / 'chr_1.bed'), par_gts_f=imputation)
    with h5py.File(tmp_path / 'export.genoarray.hdf5') as f:
        assert_array_equal(f['genotype_cols'].asstr()[:], ['proband', 'parent_sum'])
        assert 0 < f['genotypes'].shape[2] <= len(variants)
        assert_array_equal(f['genotypes'][:, 1, :], 1.5)
