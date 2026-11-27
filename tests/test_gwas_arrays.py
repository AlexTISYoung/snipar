import numpy as np
import pytest

from snipar.gtarray import gtarray
from snipar.slmm import _solve_snp_effects


@pytest.mark.parametrize("family_data", [False, True])
def test_frequencies_ignore_missing_calls_and_filter_unobserved_loci(family_data):
    genotypes = np.array([[0, np.nan, np.nan], [2, 1, np.nan]])
    if family_data:
        genotypes = np.stack([genotypes, np.zeros_like(genotypes)], axis=1)
    data = gtarray(genotypes, ids=np.array(["a", "b"]), sid=np.array(["x", "y", "z"]))
    data.compute_freqs()
    assert type(data.freqs) is np.ndarray
    np.testing.assert_allclose(data.freqs, [0.5, 0.5, np.nan])
    data.filter_maf()
    np.testing.assert_array_equal(data.sid, ["x", "y"])


@pytest.mark.parametrize("n_snps", [1, 2, 5])
def test_batched_snp_coefficients_solve_each_locus_independently(n_snps):
    rng = np.random.default_rng(132)
    predictors = rng.normal(size=(n_snps, 10, 2))
    information = predictors.transpose(0, 2, 1) @ predictors
    scores = rng.normal(size=(n_snps, 2))
    expected = [np.linalg.solve(matrix, score) for matrix, score in zip(information, scores)]
    np.testing.assert_allclose(_solve_snp_effects(information, scores), expected)
