"""Preconditioned generalized least squares for large sparse covariances."""

import numpy as np
import pandas as pd
from numba import njit
from scipy.sparse import issparse
from scipy.sparse.linalg import LinearOperator, cg, splu


@njit(cache=True, nogil=True)
def _principal_blocks(indptr, indices, values, codes, positions, rows):
    count, size = rows.shape
    blocks = np.zeros((count, size, size))
    for i in range(count):
        for j in range(size):
            row = rows[i, j]
            for index in range(indptr[row], indptr[row + 1]):
                column = indices[index]
                if codes[column] == codes[row]:
                    blocks[i, j, positions[column]] += values[index]
    return blocks


@njit(cache=True, nogil=True)
def _apply_blocks(vector, inverse, rows, out):
    count, size = rows.shape
    for i in range(count):
        for j in range(size):
            total = 0.
            for k in range(size):
                total += inverse[i, j, k] * vector[rows[i, k]]
            out[rows[i, j]] = total


class SibshipCG:
    """Solve an SPD sparse system with a sibship-block Jacobi preconditioner.

    Principal blocks are extracted from the full covariance, including any
    within-family entries from the relationship matrix. Larger families use
    diagonal preconditioning to bound the setup workspace. This class supplies
    linear solves only. It does not approximate a REML log determinant.
    """

    def __init__(self, covariance, family_labels, *, rtol=1e-9, maxiter=None, fallback=False):
        if not issparse(covariance) or covariance.shape[0] != covariance.shape[1]:
            raise ValueError('CG requires a square sparse covariance matrix')
        n = covariance.shape[0]
        if not n or not np.isfinite(rtol) or rtol <= 0 or (maxiter is not None and maxiter < 1):
            raise ValueError('CG requires a nonempty covariance and positive tolerance and iteration limit')
        if covariance.format not in ('csc', 'csr'):
            covariance = covariance.tocsr()
        self.covariance = covariance.astype(np.float64, copy=False)
        self.rtol = rtol
        self.maxiter = min(n, 1000) if maxiter is None else maxiter
        self.fallback = fallback
        self._lu = None
        diag = self.covariance.diagonal()
        if np.any(~np.isfinite(diag)) or np.any(diag <= 0):
            raise ValueError('Covariance diagonal must be finite and positive')
        self._inverse_diagonal = 1 / diag
        self._blocks = self._factor_blocks(family_labels)
        self.preconditioner = LinearOperator((n, n), matvec=self._precondition, dtype=np.float64)

    def _factor_blocks(self, labels):
        labels = np.asarray(labels)
        if labels.shape != self._inverse_diagonal.shape:
            raise ValueError('Provide one family label per covariance row')
        codes, families = pd.factorize(labels, sort=False)
        if np.any(codes < 0):
            raise ValueError('Family labels must be nonmissing')
        sizes = np.bincount(codes, minlength=len(families))
        offsets = np.cumsum(sizes) - sizes
        order = np.argsort(codes, kind='stable')
        position = np.empty(len(labels), dtype=np.int32)
        position[order] = np.arange(len(labels)) - np.repeat(offsets, sizes)
        covariance = self.covariance
        blocks = []
        # Using the CSC columns as rows is valid for symmetric covariance.
        for size in np.unique(sizes[(sizes > 1) & (sizes <= 16)]):
            starts = offsets[sizes == size]
            rows = order[starts[:, None] + np.arange(size)]
            matrices = _principal_blocks(covariance.indptr, covariance.indices, covariance.data,
                                         codes, position, rows)
            # Confirm that the preconditioner is positive definite. A caller
            # still needs to supply a positive-definite full covariance.
            chol = np.linalg.cholesky(matrices)
            inverse_chol = np.linalg.inv(chol)
            inverse = inverse_chol.transpose(0, 2, 1) @ inverse_chol
            blocks.append((rows, inverse))
        return blocks

    def _precondition(self, rhs):
        vector = rhs.reshape(-1)
        result = vector * self._inverse_diagonal
        for rows, inverse in self._blocks:
            _apply_blocks(vector, inverse, rows, result)
        return result

    def _solve_one(self, rhs):
        norm = np.linalg.norm(rhs)
        if norm == 0:
            return np.zeros_like(rhs)
        if not np.isfinite(norm):
            raise ValueError('CG right-hand side must have finite norm')
        result, status = cg(self.covariance, rhs, M=self.preconditioner, rtol=self.rtol * 0.5,
                            atol=0, maxiter=self.maxiter)
        # Check the true residual because the recursive CG residual can drift.
        relative_residual = np.linalg.norm(self.covariance @ result - rhs) / norm
        if status >= 0 and relative_residual <= self.rtol:
            return result
        if self.fallback:
            if self._lu is None:
                self._lu = splu(self.covariance.tocsc())
            return self._lu.solve(rhs)
        raise RuntimeError(f'Preconditioned CG did not reach rtol={self.rtol:g} '
                           f'(relative residual {relative_residual:.3g}, status {status})')

    def solve(self, rhs):
        """Apply the inverse to a vector or dense columns using bounded scratch space."""
        rhs = np.asarray(rhs, dtype=np.float64)
        if rhs.ndim not in (1, 2) or rhs.shape[0] != self.covariance.shape[0]:
            raise ValueError('Right-hand side rows must match the covariance')
        if self._lu is not None:
            return self._lu.solve(rhs)
        if rhs.ndim == 1:
            return self._solve_one(rhs)
        result = np.empty(rhs.shape, dtype=np.float64, order='F')
        for i in range(rhs.shape[1]):
            result[:, i] = self._solve_one(rhs[:, i]) if self._lu is None else self._lu.solve(rhs[:, i])
        return result
