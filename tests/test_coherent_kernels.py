"""Fused kernels must reproduce the retained NumPy reference equations."""

import importlib.util

import numpy as np
import pytest

from snipar import coherent, coherent_kernels as fast


# Load a separate estimator module so disabling acceleration cannot affect
# production imports or the other tests. Relative imports retain package scope.
_spec = importlib.util.spec_from_file_location(
    'snipar._coherent_numpy_reference', coherent.__file__)
original = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(original)
original.fast_kernels = None


def _sample(n, p, seed=951):
    rng = np.random.default_rng(seed)
    roots = rng.normal(size=(n, p, p))
    s = roots @ roots.transpose(0, 2, 1) + np.eye(p)*.5
    effect_roots = rng.normal(size=(p, p))
    sigma = effect_roots @ effect_roots.T * .2
    y = np.einsum('nij,nj->ni', np.linalg.cholesky(s+sigma), rng.normal(size=(n, p)))
    weight = rng.uniform(.1, 3., n)
    return y, s, weight, sigma


@pytest.mark.parametrize('p', [2, 3])
@pytest.mark.parametrize('chunk_size', [13, 20000])
def test_fused_equations_match_original_for_heteroskedastic_data(p, chunk_size):
    y, s, w, sigma = _sample(777, p)
    (information, rhs), info = fast.normal_equations(y, s, w, sigma, chunk_size, return_info=True)
    expected_information, expected_rhs = original._normal_equations(y, s, w, sigma, chunk_size)
    np.testing.assert_allclose(information, expected_information, rtol=2e-12, atol=1e-12)
    np.testing.assert_allclose(rhs, expected_rhs, rtol=2e-12, atol=1e-12)
    np.testing.assert_allclose(fast.criterion(y, s, w, sigma, chunk_size),
                               original._criterion(y, s, w, sigma, chunk_size), rtol=2e-13, atol=1e-13)
    assert info['n_fallback'] == 0 if fast.AVAILABLE else info['n_fallback'] == len(y)


@pytest.mark.parametrize('p', [2, 3])
def test_ill_conditioned_rows_use_original_numpy_fallback(p):
    y, s, w, sigma = _sample(100, p)
    sigma[:] = 0
    correlations = np.r_[np.linspace(0, .999, 90), np.repeat(1-1e-10, 10)]
    for k, r in enumerate(correlations):
        s[k] = np.eye(p)*(1-r) + np.ones((p, p))*r
    (information, rhs), info = fast.normal_equations(y, s, w, sigma, 17, return_info=True)
    expected_information, expected_rhs = original._normal_equations(y, s, w, sigma, 17)
    np.testing.assert_allclose(information, expected_information, rtol=2e-12, atol=1e-12)
    np.testing.assert_allclose(rhs, expected_rhs, rtol=2e-12, atol=1e-12)
    objective, criterion_info = fast.criterion(y, s, w, sigma, 17, return_info=True)
    np.testing.assert_allclose(objective, original._criterion(y, s, w, sigma, 17), rtol=2e-12)
    assert info['n_fallback'] >= 10
    assert criterion_info['n_fallback'] >= 10


@pytest.mark.parametrize('p', [2, 3])
def test_component_and_global_scales_and_condition_range(p):
    rng = np.random.default_rng(348)
    condition = np.geomspace(1, 1e12, 101)
    s = []
    for value in condition:
        rotation = np.linalg.qr(rng.normal(size=(p, p)))[0]
        s.append((rotation*np.geomspace(1, 1/value, p)) @ rotation.T)
    s = np.array(s)*np.geomspace(1e-5, 1e5, len(condition))[:, None, None]
    y = np.einsum('nij,nj->ni', np.linalg.cholesky(s), rng.normal(size=(len(s), p)))
    w, sigma = np.ones(len(s)), np.zeros((p, p))
    (actual_i, actual_h), info = fast.normal_equations(y, s, w, sigma, return_info=True)
    expected_i, expected_h = original._normal_equations(y, s, w, sigma, 20000)
    assert np.max(abs(actual_i-expected_i))/np.max(abs(expected_i)) < 2e-12
    assert np.max(abs(actual_h-expected_h))/np.max(abs(expected_h)) < 2e-12
    np.testing.assert_allclose(fast.criterion(y, s, w, sigma),
                               original._criterion(y, s, w, sigma, 20000), rtol=2e-12)
    assert info['n_fallback'] > 50


@pytest.mark.parametrize('p', [2, 3])
def test_without_numba_uses_original_numpy(monkeypatch, p):
    y, s, w, sigma = _sample(127, p)
    monkeypatch.setattr(fast, 'AVAILABLE', False)
    actual, info = fast.normal_equations(y, s, w, sigma, 21, return_info=True)
    expected = original._normal_equations(y, s, w, sigma, 21)
    np.testing.assert_array_equal(actual[0], expected[0])
    np.testing.assert_array_equal(actual[1], expected[1])
    np.testing.assert_array_equal(fast.criterion(y, s, w, sigma, 21),
                                  original._criterion(y, s, w, sigma, 21))
    assert info['n_fallback'] == len(y)


def test_invalid_shapes_raise():
    y, s, w, sigma = _sample(8, 3)
    with pytest.raises(ValueError, match='match'):
        fast.normal_equations(y, s, w[:-1], sigma)
    with pytest.raises(ValueError, match='chunk_size'):
        fast.criterion(y, s, w, sigma, chunk_size=0)


@pytest.mark.skipif(not fast.AVAILABLE, reason='Regression guards the fused accumulation')
def test_million_row_objective_keeps_backtracking_precision():
    # Naively summing this constant contribution over 1.5M rows produced a
    # 2.9e-11 objective error, larger than the estimator's 1e-13 relative slack.
    n = 1_500_000
    y = np.ones((n, 3))
    s = np.broadcast_to(np.eye(3), (n, 3, 3)).copy()
    w, sigma = np.ones(n), np.eye(3)
    expected = .5*(3*np.log(2)+1.5)
    actual = fast.criterion(y, s, w, sigma)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=3e-15)
