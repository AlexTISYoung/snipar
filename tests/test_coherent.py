"""Independent numerical checks of the joint, noise-corrected moment estimator."""

import numpy as np
import pytest


from snipar import coherent as estimator


def _exact_effects(second_moment, repetitions=20):
    """A balanced finite sample with the prescribed uncentered second moment."""
    p = len(second_moment)
    rows = np.sqrt(p) * np.linalg.cholesky(second_moment).T
    return np.tile(np.concatenate([rows, -rows]), (repetitions, 1))


def _heteroskedastic_sample(n=600, p=2, seed=891):
    rng = np.random.default_rng(seed)
    signal = np.array([[1.4, 0.23, 0.12],
                       [0.23, 0.8, -0.08],
                       [0.12, -0.08, 0.65]])[:p, :p]
    sampling_correlation = np.array([[1.0, -0.55, 0.2],
                                     [-0.55, 1.0, -0.1],
                                     [0.2, -0.1, 1.0]])[:p, :p]
    sampling_sd = rng.uniform(0.12, 0.6, size=(n, p))
    sampling = (sampling_sd[:, :, None] * sampling_correlation
                * sampling_sd[:, None, :])
    effects = np.einsum("nij,nj->ni", np.linalg.cholesky(sampling + signal),
                        rng.normal(size=(n, p)))
    ld = rng.uniform(1, 4, size=n)
    return effects, sampling, ld


def _metric_map(result, field="estimates"):
    return dict(zip(result["metric_names"], result[field]))


def test_gaussian_products_include_all_cross_moment_covariances():
    # For the order XX, XY, YY, Isserlis' identity gives this explicit matrix.
    total = np.array([[[2.0, -0.6], [-0.6, 3.0]],
                      [[1.1, 0.2], [0.2, 0.8]]])
    expected = []
    for a, b, c in [(2.0, -0.6, 3.0), (1.1, 0.2, 0.8)]:
        expected.append([[2*a*a, 2*a*b, 2*b*b],
                         [2*a*b, a*c+b*b, 2*b*c],
                         [2*b*b, 2*b*c, 2*c*c]])
    np.testing.assert_allclose(
        estimator.gaussian_product_covariance(total), expected)


def test_three_component_precision_reduction_matches_explicit_product_inverse():
    rng = np.random.default_rng(606)
    roots = rng.normal(size=(8, 3, 3))
    total = roots @ roots.transpose(0, 2, 1) + 0.3*np.eye(3)
    total *= np.geomspace(1e-6, 1e3, len(total))[:, None, None]
    expected = np.linalg.inv(estimator.gaussian_product_covariance(total))
    np.testing.assert_allclose(estimator._product_precision(total), expected,
                               rtol=2e-12, atol=1e-12)


def test_homoskedastic_solution_matches_exact_corrected_second_moment():
    signal = np.array([[0.7, -0.11], [-0.11, 0.4]])
    sampling = np.array([[0.15, -0.075], [-0.075, 0.2]])
    effects = _exact_effects(signal + sampling)
    fitted, _ = estimator.fit_covariance(
        effects, np.broadcast_to(sampling, (len(effects), 2, 2)).copy())
    np.testing.assert_allclose(fitted, signal, rtol=2e-6, atol=1e-9)


def test_joint_fit_satisfies_independent_heteroskedastic_gls_equations():
    effects, sampling, ld = _heteroskedastic_sample()
    fitted, _ = estimator.fit_covariance(
        effects, sampling, ldscores=ld, tol=1e-9, max_iter=200)
    assert np.linalg.eigvalsh(fitted).min() > 0.1
    pair = list(zip(*np.triu_indices(2)))
    information = np.zeros((3, 3))
    rhs = np.zeros(3)
    for y, s, l in zip(effects, sampling, ld):
        omega = fitted + s
        product_covariance = np.array([
            [omega[a, c]*omega[b, d] + omega[a, d]*omega[b, c]
             for c, d in pair] for a, b in pair])
        precision = np.linalg.inv(product_covariance) / l
        observed = np.array([y[a]*y[b] - s[a, b] for a, b in pair])
        information += precision
        rhs += precision @ observed
    expected = np.linalg.solve(information, rhs)
    np.testing.assert_allclose(fitted[np.triu_indices(2)], expected,
                               rtol=2e-6, atol=1e-8)


def test_psd_projection_uses_information_metric():
    # Symmetry reduces the constrained optimum to [[t,t],[t,t]], t=11/6.
    # Simply clipping the unconstrained matrix's eigenvalues gives t=1.5.
    information = np.diag([1.0, 10.0, 1.0])
    unconstrained = np.array([1.0, 2.0, 1.0])
    fitted = estimator.psd_gls(information, information @ unconstrained, 2)
    np.testing.assert_allclose(fitted, np.full((2, 2), 11.0/6),
                               rtol=3e-6, atol=1e-8)
    assert np.linalg.eigvalsh(fitted).min() >= -1e-10


@pytest.mark.parametrize("rank", [0, 1, 2])
def test_three_component_psd_solver_recovers_known_boundary_optimum(rank):
    # Construct primal and dual PSD matrices with complementary ranges.
    # The resulting KKT conditions certify the answer independently of a solver.
    rng = np.random.default_rng(764 + rank)
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    expected = rotation @ np.diag([0.4]*rank + [0]*(3-rank)) @ rotation.T
    dual = rotation @ np.diag([0]*rank + [0.7]*(3-rank)) @ rotation.T
    root = rng.normal(size=(6, 6))
    information = root.T @ root + 0.4*np.eye(6)
    i, j = np.triu_indices(3)
    rhs = (information @ expected[i, j]
           - dual[i, j]*np.where(i == j, 1, 2))
    fitted = estimator.psd_gls(information, rhs, 3)
    np.testing.assert_allclose(fitted, expected, rtol=1e-5, atol=1e-6)
    assert np.linalg.eigvalsh(fitted).min() >= -1e-10


def test_parental_matrix_and_population_outputs_are_one_linear_transform():
    signal = np.array([[0.7, 0.1, -0.03],
                       [0.1, 0.4, 0.12],
                       [-0.03, 0.12, 0.3]])
    sampling = np.array([[0.1, -0.02, -0.01],
                         [-0.02, 0.08, 0.01],
                         [-0.01, 0.01, 0.12]])
    effects = _exact_effects(signal + sampling)
    result = estimator.fit_coherent(
        effects, np.broadcast_to(sampling, (len(effects), 3, 3)).copy(),
        n_blocks=10, effect_names=["direct", "paternal", "maternal"])
    transform = np.array([[1, 0, 0], [0, 0.5, 0.5], [1, 0.5, 0.5]])
    expected = transform @ signal @ transform.T
    np.testing.assert_allclose(result["covariance"], signal,
                               rtol=3e-6, atol=1e-8)
    np.testing.assert_allclose(result["expanded_covariance"], expected,
                               rtol=3e-6, atol=1e-8)
    metric = _metric_map(result)
    dd, nn, pp = np.diag(expected)
    np.testing.assert_allclose(metric["r_direct_population"],
                               expected[0, 2]/np.sqrt(dd*pp), rtol=3e-6)
    np.testing.assert_allclose(metric["r_direct_avg_NTC"],
                               expected[0, 1]/np.sqrt(dd*nn), rtol=3e-6)
    for numerator, denominator, expected_ratio in [
        ("direct", "population", dd/pp),
        ("avg_NTC", "population", nn/pp),
        ("direct", "avg_NTC", dd/nn),
    ]:
        forward = metric[f"sqrt_var_ratio_{numerator}_{denominator}"]
        reverse = metric[f"sqrt_var_ratio_{denominator}_{numerator}"]
        np.testing.assert_allclose(forward, np.sqrt(expected_ratio), rtol=3e-6)
        np.testing.assert_allclose(forward*reverse, 1.0, rtol=1e-10)
    # P-D-N has identically zero variance, not a third fitted effect.
    contrast = np.array([-1.0, -1.0, 1.0])
    np.testing.assert_allclose(result["expanded_covariance"] @ contrast,
                               0, atol=1e-12)


def test_jackknife_reestimates_precision_and_propagates_joint_uncertainty():
    effects, sampling, ld = _heteroskedastic_sample(n=240)
    blocks = np.repeat(np.arange(6), 40)
    result = estimator.fit_coherent(
        effects, sampling, ldscores=ld, blocks=blocks, n_blocks=6,
        effect_names=["direct", "avg_NTC"], jackknife="refit", tol=1e-9,
        max_iter=200)
    expected_delete = []
    expected_correlation = []
    expected_ratio = []
    for block in range(6):
        keep = blocks != block
        fitted, _ = estimator.fit_covariance(
            effects[keep], sampling[keep], ldscores=ld[keep], tol=1e-9,
            max_iter=200)
        expected_delete.append(fitted[np.triu_indices(2)])
        dd, dn, nn = fitted[0, 0], fitted[0, 1], fitted[1, 1]
        pp = dd + 2*dn + nn
        expected_correlation.append((dd+dn)/np.sqrt(dd*pp))
        expected_ratio.append(np.sqrt(dd/pp))
    expected_delete = np.array(expected_delete)
    np.testing.assert_allclose(result["moment_delete"], expected_delete,
                               rtol=1e-5, atol=1e-8)
    centered = expected_delete - expected_delete.mean(axis=0)
    expected_covariance = 5/6 * centered.T @ centered
    np.testing.assert_allclose(result["moment_covariance"], expected_covariance,
                               rtol=2e-5, atol=1e-9)
    np.testing.assert_allclose(result["moment_se"]**2,
                               np.diag(expected_covariance), rtol=2e-5)
    for name, expected in [
        ("r_direct_population", expected_correlation),
        ("sqrt_var_ratio_direct_population", expected_ratio),
    ]:
        column = list(result["metric_names"]).index(name)
        np.testing.assert_allclose(result["delete_values"][:, column], expected,
                                   rtol=1e-5, atol=1e-8)
        expected = np.asarray(expected)
        variance = 5/6 * np.sum((expected - expected.mean())**2)
        np.testing.assert_allclose(result["se"][column]**2, variance,
                                   rtol=3e-5, atol=1e-10)


def test_units_and_snp_specific_allele_coding_do_not_change_ratios():
    effects, sampling, ld = _heteroskedastic_sample(n=160)
    kwargs = dict(ldscores=ld, n_blocks=4,
                  effect_names=["direct", "avg_NTC"])
    result = estimator.fit_coherent(effects, sampling, **kwargs)
    signs = np.where(np.arange(len(effects)) % 3 == 0, -1, 1)
    altered = estimator.fit_coherent(
        1e-4 * effects * signs[:, None], 1e-8 * sampling, **kwargs)
    np.testing.assert_allclose(altered["covariance"] / 1e-8,
                               result["covariance"], rtol=2e-5, atol=1e-8)
    moment_columns = np.array([name.startswith(("var_", "cov_"))
                               for name in result["metric_names"]])
    for field, tolerance in [("estimates", 2e-5), ("se", 2e-4)]:
        converted = np.asarray(altered[field]).copy()
        converted[moment_columns] /= 1e-8
        np.testing.assert_allclose(converted, result[field], rtol=tolerance,
                                   atol=1e-8, equal_nan=True)


def test_zero_signal_has_undefined_correlations_and_ratios():
    result = estimator.fit_coherent(
        np.zeros((40, 2)), np.tile(np.eye(2), (40, 1, 1)), n_blocks=4,
        effect_names=["direct", "avg_NTC"])
    np.testing.assert_allclose(result["covariance"], 0, atol=1e-10)
    metrics = _metric_map(result)
    errors = _metric_map(result, "se")
    for name in ["r_direct_population", "r_direct_avg_NTC",
                 "sqrt_var_ratio_direct_population",
                 "sqrt_var_ratio_avg_NTC_population",
                 "sqrt_var_ratio_direct_avg_NTC"]:
        assert np.isnan(metrics[name])
        assert np.isnan(errors[name])


def test_common_weight_diagnostic_preserves_the_fixed_finite_panel_target():
    effects, sampling, ld = _heteroskedastic_sample(n=120)
    result = estimator.fit_coherent(effects, sampling, ldscores=ld, n_blocks=4)
    corrected = effects[:, :, None] * effects[:, None, :] - sampling
    expected = np.einsum('n,nij->ij', 1 / ld, corrected) / np.sum(1 / ld)
    np.testing.assert_allclose(result['common_weight_covariance'], expected,
                               rtol=1e-12, atol=1e-14)
    expected_delete = []
    for excluded in np.array_split(np.arange(len(effects)), 4):
        keep = np.ones(len(effects), dtype=bool)
        keep[excluded] = False
        matrix = (np.einsum('n,nij->ij', 1 / ld[keep], corrected[keep])
                  / np.sum(1 / ld[keep]))
        expected_delete.append(matrix[np.triu_indices(2)])
    np.testing.assert_allclose(result['common_weight_moment_delete'], expected_delete)
    centered = expected_delete - np.mean(expected_delete, axis=0)
    np.testing.assert_allclose(result['common_weight_moment_covariance'],
                               .75 * centered.T @ centered)


def test_frozen_jackknife_uses_explicit_full_fit_product_weights():
    effects, sampling, ld = _heteroskedastic_sample(n=120)
    result = estimator.fit_coherent(effects, sampling, ldscores=ld, n_blocks=4,
                                    jackknife='frozen')
    precision = np.linalg.inv(estimator.gaussian_product_covariance(
        sampling + result['covariance'])) / ld[:, None, None]
    i, j = np.triu_indices(2)
    products = effects[:, i] * effects[:, j] - sampling[:, i, j]
    expected = []
    for excluded in np.array_split(np.arange(len(effects)), 4):
        keep = np.ones(len(effects), dtype=bool)
        keep[excluded] = False
        info = precision[keep].sum(axis=0)
        rhs = np.einsum('nij,nj->i', precision[keep], products[keep])
        expected.append(estimator.psd_gls(info, rhs, 2)[i, j])
    np.testing.assert_allclose(result['moment_delete'], expected, rtol=1e-7, atol=1e-9)
    assert 'conditional' in ' '.join(result['diagnostics']['warnings'])


@pytest.mark.parametrize("fault", ["nan_effect", "asymmetric", "indefinite",
                                    "singular", "zero_ld", "nan_ld"])
def test_invalid_inputs_are_rejected(fault):
    effects, sampling, ld = _heteroskedastic_sample(n=40)
    if fault == "nan_effect":
        effects[0, 0] = np.nan
    elif fault == "asymmetric":
        sampling[0, 0, 1] += 0.1
    elif fault == "indefinite":
        sampling[0] = [[1, 1.2], [1.2, 1]]
    elif fault == "singular":
        sampling[0] = [[1, 1], [1, 1]]
    elif fault == "zero_ld":
        ld[0] = 0
    elif fault == "nan_ld":
        ld[0] = np.nan
    with pytest.raises(ValueError):
        estimator.fit_covariance(effects, sampling, ldscores=ld)


def test_homogeneous_covariance_refits_weighted_psd_solution():
    # A correlated sampling covariance makes clipping eigenvalues in the
    # original units differ from clipping after whitening.
    cov = np.array([[1.2, .35], [.35, .4]])
    rng = np.random.default_rng(31)
    effects = rng.normal(size=(120, 2)) @ np.array([[1.8, .2], [.1, .4]])
    ld = rng.uniform(1, 3, len(effects))
    blocks = np.repeat(np.arange(6), 20)
    sampling = np.broadcast_to(cov, (len(effects), 2, 2))
    fit = estimator.fit_coherent(effects, sampling, ld, blocks=blocks)
    l = np.linalg.cholesky(cov)

    def weighted_oracle(keep):
        whitened = np.linalg.solve(l, effects[keep].T).T
        moment = whitened.T @ (whitened / ld[keep, None]) / np.sum(1 / ld[keep])
        spectrum, eigenvectors = np.linalg.eigh(moment)
        return l @ ((eigenvectors * np.maximum(spectrum - 1, 0)) @ eigenvectors.T) @ l.T

    np.testing.assert_allclose(fit['covariance'], weighted_oracle(np.ones(len(effects), dtype=bool)),
                               rtol=2e-7, atol=2e-9)
    i, j = np.triu_indices(2)
    np.testing.assert_allclose(fit['moment_delete'],
                               [weighted_oracle(blocks != block)[i, j] for block in np.unique(blocks)],
                               rtol=2e-7, atol=2e-9)
