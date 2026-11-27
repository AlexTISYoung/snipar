"""Numerical and scheduling checks for memory-bounded family GWAS."""

from concurrent.futures import ThreadPoolExecutor
import threading

import numpy as np
from numpy.testing import assert_allclose
import pytest

from snipar.gwas import _batch_slices, _bounded_results, split_batches
from snipar.slmm import LinearMixedModel, build_sib_arr


@pytest.mark.parametrize('fixed_effects', [False, True])
def test_structured_covariance_matches_sparse_reml_and_snp_gls(fixed_effects):
    rng = np.random.default_rng(239)
    labels = rng.permutation(np.repeat(['a', 'b', 'c', 'd'], [4, 1, 8, 17]))
    y = rng.normal(size=len(labels))
    covariates = rng.normal(size=(len(labels), 2)) if fixed_effects else None
    sigma = (.48, .27)
    sparse = LinearMixedModel(y, [build_sib_arr(labels)], sigma, covariates, fixed_effects)
    blocked = LinearMixedModel(y, (), sigma, covariates, fixed_effects, family_labels=labels)
    assert_allclose(blocked.V_logdet, sparse.V_logdet)
    assert_allclose(blocked.reml_loglik, sparse.reml_loglik)
    assert_allclose(blocked.grad, sparse.grad, atol=1e-12)
    assert_allclose(blocked.hessian, sparse.hessian, atol=1e-12)

    gts = rng.integers(0, 3, (len(y), 3, 5)).astype(np.float32)
    gts[:, 1, 0] = gts[:, 0, 0]  # The rank check must reject a single singular SNP.
    gts[2, 0, 2] = np.nan
    assert blocked._component_mats is None
    for actual, expected in zip(blocked.fit_snps_eff(gts), sparse.fit_snps_eff(gts)):
        assert_allclose(actual, expected, atol=1e-12, equal_nan=True)
    assert blocked._component_mats is None


def test_preconditioned_sparse_gls_matches_lu_on_common_outputs():
    rng = np.random.default_rng(933)
    rows = np.arange(20, dtype='uint32')
    labels = np.repeat(np.arange(10), 2)
    diagonal = (np.ones(20), rows, rows)
    adjacency = (np.full(18, .1), rows[2:], rows[:-2])
    component = tuple(np.concatenate(parts) for parts in zip(diagonal, adjacency))
    components = [component, build_sib_arr(labels)]
    y = rng.normal(size=20)
    covariates = rng.normal(size=(20, 2))
    effects = rng.normal(size=(20, 2, 8)).astype(np.float32)
    dense = LinearMixedModel(y, components, (.41, .31, .5), covariates, True)
    iterative = LinearMixedModel(y, components, (.41, .31, .5), covariates, True,
                                  family_labels=labels, snp_solver='cg')
    for actual, expected in zip(iterative.fit_snps_eff(effects), dense.fit_snps_eff(effects)):
        assert_allclose(actual, expected, rtol=2e-8, atol=1e-9)
    assert '_cached_V_lu' not in iterative.__dict__
    assert_allclose(iterative.reml_loglik, dense.reml_loglik, atol=1e-12)


@pytest.mark.parametrize('loci,workers,maximum', [(1, 4, None), (37, 3, 5), (29, 4, 1), (0, 1, None)])
def test_batch_slices_cover_each_variant_once(loci, workers, maximum):
    blocks = list(_batch_slices(loci, workers, maximum))
    known = split_batches(loci, workers, maximum)
    for block, indices in zip(blocks, known, strict=True):
        np.testing.assert_array_equal(np.arange(loci)[block], indices)
    where = [i for block in blocks for i in range(block.start, block.stop)]
    assert where == list(range(loci))
    if maximum is not None:
        assert max(block.stop - block.start for block in blocks) <= maximum


def test_bounded_submission_does_not_consume_future_batches():
    released = threading.Event()
    taken = []

    def data():
        for i in range(20):
            taken.append(i)
            yield i

    def run(value):
        if value == 1:
            assert released.wait(timeout=5)
        return value

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = _bounded_results(executor, run, data(), workers=2)
        try:
            assert next(results) == 0
            assert taken == [0, 1]
        finally:
            released.set()
        assert sorted(results) == list(range(1, 20))


def test_structured_trio_sibling_covariance_matches_sparse_without_materializing_v():
    rng = np.random.default_rng(947)
    labels = np.repeat(np.arange(17), 2)
    y = rng.normal(size=len(labels))
    gts = rng.integers(0, 3, size=(len(y), 3, 6)).astype(np.float32)
    gts[1::2, 2] = gts[1::2, 1]
    first, second = np.arange(0, len(y), 2), np.arange(1, len(y), 2)
    covariance = (.4, .8)
    sparse = LinearMixedModel(y, (build_sib_arr(labels),), covariance, add_intercept=True)
    blocked = LinearMixedModel(y, (), covariance, add_intercept=True, family_labels=labels)
    reference = sparse.trios_sibs_est(gts, first, second)
    for result in (blocked.trios_sibs_est(gts, first, second), blocked.trios_sibs_est(gts, first, second)):
        for actual, expected in zip(result, reference, strict=True):
            assert_allclose(actual, expected, atol=1e-12, rtol=1e-12)
    assert blocked._component_mats is None
