import numpy as np
import pytest
from scipy.sparse import csc_matrix, csr_matrix

from snipar.gwas_cg import SibshipCG


@pytest.mark.parametrize('sparse', [csc_matrix, csr_matrix])
def test_cg_uses_exact_principal_family_blocks_and_matches_dense(sparse):
    labels = np.array(['a', 'b', 'c', 'a', 'b', 'b'])
    rng = np.random.default_rng(624)
    relatedness = rng.normal(size=(6, 6))
    covariance = relatedness @ relatedness.T + (labels[:, None] == labels[None, :]) + np.eye(6)
    partitioned = covariance * (labels[:, None] == labels[None, :])
    rhs = rng.normal(size=(6, 2))
    solver = SibshipCG(sparse(covariance), labels)
    np.testing.assert_allclose(solver.preconditioner @ rhs, np.linalg.solve(partitioned, rhs), atol=1e-12)
    np.testing.assert_allclose(solver.solve(rhs), np.linalg.solve(covariance, rhs), atol=2e-10, rtol=1e-9)
    np.testing.assert_array_equal(solver.solve(np.zeros(6)), 0)


def test_cg_rejects_unconverged_solves_unless_lu_fallback_requested():
    matrix = np.array([[4., 1., 0.], [1., 3., 1.], [0., 1., 4.]])
    rhs = np.array([2., 3., 5.])
    labels = np.arange(3)
    covariance = csc_matrix(matrix)
    with pytest.raises(RuntimeError, match='relative residual'):
        SibshipCG(covariance, labels, maxiter=1).solve(rhs)
    solver = SibshipCG(covariance, labels, maxiter=1, fallback=True)
    np.testing.assert_allclose(solver.solve(rhs), np.linalg.solve(matrix, rhs))
    np.testing.assert_allclose(solver.solve(rhs[:, None]), np.linalg.solve(matrix, rhs[:, None]))
