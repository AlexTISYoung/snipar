"""Optional fused kernels for the unchanged Gaussian covariance equations.

Only row arithmetic and summation are accelerated. PSD fitting, weights,
convergence tolerances and jackknife inference remain the caller's responsibility.
Inputs must already have passed the estimator's validation. Difficult total
covariances use the original NumPy expressions, with no regularization.
"""

from __future__ import annotations

import numpy as np

try:
    from numba import njit as _numba_njit
    AVAILABLE = True
except ImportError:
    _numba_njit = None
    AVAILABLE = False


def njit(*args, **kwargs):
    """Use the package-qualified module name for CLI and API Numba caches."""
    def decorate(function):
        if _numba_njit is None:
            return function
        try:
            return _numba_njit(*args, **kwargs)(function)
        except RuntimeError as error:
            if not kwargs.get('cache') or 'no locator available for file' not in str(error):
                raise
            # Shared installations may have no writable Numba cache.
            return _numba_njit(*args, **dict(kwargs, cache=False))(function)
    return decorate


@njit(inline='always', cache=True, nogil=True)
def _row_algebra(y, s, sigma, row, p, need_criterion):
    a = s[row, 0, 0] + sigma[0, 0]
    b = s[row, 1, 0] + sigma[1, 0]
    d = s[row, 1, 1] + sigma[1, 1]
    c, e, f = 0.0, 0.0, 0.0
    if p == 3:
        c = s[row, 2, 0] + sigma[2, 0]
        e = s[row, 2, 1] + sigma[2, 1]
        f = s[row, 2, 2] + sigma[2, 2]
    scale = max(abs(a), abs(d), abs(f))
    invalid = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False)
    if not np.isfinite(scale) or scale <= 0:
        return invalid
    a, b, c, d, e, f = a/scale, b/scale, c/scale, d/scale, e/scale, f/scale
    if a <= 1.4210854715202004e-14:
        return invalid
    l00 = np.sqrt(a)
    l10 = b/l00
    pivot11 = d-l10*l10
    if pivot11 <= 1.4210854715202004e-14:
        return invalid
    l11 = np.sqrt(pivot11)
    b00, b11 = 1/l00, 1/l11
    b10 = -l10*b00*b11
    ic, ie, iff = 0.0, 0.0, 0.0
    l20, l21, l22 = 0.0, 0.0, 1.0
    if p == 2:
        ia, ib, idd = b00*b00+b10*b10, b10*b11, b11*b11
        matrix_norm = max(a+abs(b), d+abs(b))
        inverse_norm = max(ia+abs(ib), idd+abs(ib))
        residual = max(abs(a*ia+b*ib-1), abs(b*ib+d*idd-1),
                       abs(a*ib+b*idd), abs(b*ia+d*ib))
    else:
        l20 = c/l00
        l21 = (e-l20*l10)/l11
        pivot22 = f-l20*l20-l21*l21
        if pivot22 <= 1.4210854715202004e-14:
            return invalid
        l22 = np.sqrt(pivot22)
        b22 = 1/l22
        b21 = -l21*b11*b22
        b20 = -(l20*b00+l21*b10)*b22
        ia = b00*b00+b10*b10+b20*b20
        ib, ic = b10*b11+b20*b21, b20*b22
        idd, ie, iff = b11*b11+b21*b21, b21*b22, b22*b22
        matrix_norm = max(a+abs(b)+abs(c), d+abs(b)+abs(e), f+abs(c)+abs(e))
        inverse_norm = max(ia+abs(ib)+abs(ic), idd+abs(ib)+abs(ie), iff+abs(ic)+abs(ie))
        residual = max(abs(a*ia+b*ib+c*ic-1), abs(b*ib+d*idd+e*ie-1),
                       abs(c*ic+e*ie+f*iff-1), abs(a*ib+b*idd+c*ie),
                       abs(a*ic+b*ie+c*iff), abs(b*ia+d*ib+e*ic),
                       abs(b*ic+d*ie+e*iff), abs(c*ia+e*ib+f*ic),
                       abs(c*ib+e*idd+f*ie))
    condition = matrix_norm*inverse_norm
    if not np.isfinite(condition) or condition > 1e4 or not np.isfinite(residual) or residual > 5e-13:
        return invalid
    logdet, quadratic = 0.0, 0.0
    if need_criterion:
        logdet = p*np.log(scale)+2*(np.log(l00)+np.log(l11))
        root_scale = np.sqrt(scale)
        z0 = (y[row, 0]/root_scale)/l00
        z1 = (y[row, 1]/root_scale-l10*z0)/l11
        quadratic = z0*z0+z1*z1
        if p == 3:
            logdet += 2*np.log(l22)
            z2 = (y[row, 2]/root_scale-l20*z0-l21*z1)/l22
            quadratic += z2*z2
        if not np.isfinite(logdet) or not np.isfinite(quadratic):
            return invalid
    ia, ib, ic, idd, ie, iff = ia/scale, ib/scale, ic/scale, idd/scale, ie/scale, iff/scale
    if not np.isfinite(ia) or not np.isfinite(idd) or not np.isfinite(iff):
        return invalid
    return ia, ib, ic, idd, ie, iff, logdet, quadratic, True


@njit(cache=True, nogil=True)
def _normal_equations_jit(y, s, weight, sigma):
    n, p = y.shape
    q = p*(p+1)//2
    left, right = np.empty(q, np.int64), np.empty(q, np.int64)
    mult = np.empty(q)
    k = 0
    for a in range(p):
        for b in range(a, p):
            left[k], right[k], mult[k] = a, b, 1.0 if a == b else 2.0
            k += 1
    information, rhs = np.zeros((q, q)), np.zeros(q)
    inverse, products, contributions = np.empty((3, 3)), np.empty(q), np.empty(q)
    bad = np.zeros(n, np.bool_)
    for row in range(n):
        ia, ib, ic, idd, ie, iff, _, _, good = _row_algebra(y, s, sigma, row, p, False)
        if not good:
            bad[row] = True
            continue
        inverse[0, 0], inverse[0, 1], inverse[0, 2] = ia, ib, ic
        inverse[1, 0], inverse[1, 1], inverse[1, 2] = ib, idd, ie
        inverse[2, 0], inverse[2, 1], inverse[2, 2] = ic, ie, iff
        for k in range(q):
            products[k] = y[row, left[k]]*y[row, right[k]]-s[row, left[k], right[k]]
            contributions[k] = 0.
        # Product precision is symmetric. Accumulate each off-diagonal
        # once and use it in both coordinates of the right-hand side.
        for k in range(q):
            a, b = left[k], right[k]
            for ell in range(k, q):
                c, d = left[ell], right[ell]
                precision = ((inverse[a, c]*inverse[b, d]+inverse[a, d]*inverse[b, c])
                             * mult[k]*mult[ell]/4.0)
                precision *= weight[row]
                information[k, ell] += precision
                contributions[k] += precision*products[ell]
                if ell != k:
                    contributions[ell] += precision*products[k]
            rhs[k] += contributions[k]
    for k in range(q):
        for ell in range(k):
            information[k, ell] = information[ell, k]
    return information, rhs, bad


@njit(cache=True, nogil=True)
def _criterion_jit(y, s, weight, sigma):
    n, p = y.shape
    bad = np.zeros(n, np.bool_)
    value = 0.0
    compensation = 0.0
    for row in range(n):
        _, _, _, _, _, _, logdet, quadratic, good = _row_algebra(y, s, sigma, row, p, True)
        if good:
            # Backtracking can compare objectives within 1e-13. Compensated
            # summation prevents SNP-count-dependent drift at that tolerance.
            adjusted = .5*weight[row]*(logdet+quadratic)-compensation
            updated = value+adjusted
            compensation = (updated-value)-adjusted
            value = updated
        else:
            bad[row] = True
    return value, bad


def _check(y, s, weight, sigma, chunk_size):
    y, s, weight, sigma = (np.asarray(value, dtype=float) for value in (y, s, weight, sigma))
    if y.ndim != 2 or y.shape[1] not in (2, 3):
        raise ValueError('y must have shape (n,2) or (n,3)')
    n, p = y.shape
    if s.shape != (n, p, p) or sigma.shape != (p, p) or weight.shape != (n,):
        raise ValueError('s, weight and sigma must match y')
    if int(chunk_size) != chunk_size or chunk_size < 1:
        raise ValueError('chunk_size must be a positive integer')
    return y, s, weight, sigma, int(chunk_size)


def normal_equations(y, s, weight, sigma, chunk_size=20000, return_info=False):
    """Return the original Gaussian GLS information and RHS on validated data.

    A fused serial kernel replaces per-SNP 6-by-6 allocations. No fastmath,
    parallel reduction, reweighting or statistical approximation is enabled.
    Rejected rows, or every row without Numba, use the original NumPy algebra.
    """
    y, s, weight, sigma, chunk_size = _check(y, s, weight, sigma, chunk_size)
    p = y.shape[1]
    i, j = np.triu_indices(p)
    information, rhs = np.zeros((len(i), len(i))), np.zeros(len(i))
    n_fallback = 0
    mult = np.where(i == j, 1.0, 2.0)
    a, b, c, d = i[:, None], j[:, None], i[None, :], j[None, :]
    for start in range(0, len(y), chunk_size):
        stop = min(start+chunk_size, len(y))
        if AVAILABLE:
            chunk_information, chunk_rhs, bad = _normal_equations_jit(
                y[start:stop], s[start:stop], weight[start:stop], sigma)
        else:
            chunk_information, chunk_rhs = np.zeros_like(information), np.zeros_like(rhs)
            bad = np.ones(stop-start, bool)
        rows = start+np.flatnonzero(bad)
        n_fallback += len(rows)
        if len(rows):
            y_rows, s_rows = y[rows], s[rows]
            inverse = np.linalg.inv(s_rows+sigma)
            precision = ((inverse[:, a, c]*inverse[:, b, d]+inverse[:, a, d]*inverse[:, b, c])
                         * mult[:, None]*mult[None, :]/4.0)
            precision *= weight[rows, None, None]
            # Index rows before selecting traits, as in the NumPy estimator.
            # Joint advanced indexing makes the trait axis contiguous instead,
            # which changes einsum's reduction for near-singular covariances.
            products = y_rows[:, i]*y_rows[:, j]-s_rows[:, i, j]
            chunk_information += precision.sum(axis=0)
            chunk_rhs += np.einsum('nij,nj->ni', precision, products).sum(axis=0)
        information += chunk_information
        rhs += chunk_rhs
    result = information, rhs
    if return_info:
        return result, dict(n_matrices=len(y), n_fallback=n_fallback, compiled=AVAILABLE)
    return result


def criterion(y, s, weight, sigma, chunk_size=20000, return_info=False):
    """Return the same weighted mean Gaussian objective on validated data."""
    y, s, weight, sigma, chunk_size = _check(y, s, weight, sigma, chunk_size)
    value, compensation, n_fallback = 0.0, 0.0, 0
    for start in range(0, len(y), chunk_size):
        stop = min(start+chunk_size, len(y))
        if AVAILABLE:
            chunk_value, bad = _criterion_jit(y[start:stop], s[start:stop], weight[start:stop], sigma)
        else:
            chunk_value, bad = 0.0, np.ones(stop-start, bool)
        rows = start+np.flatnonzero(bad)
        n_fallback += len(rows)
        if len(rows):
            total = s[rows]+sigma
            sign, logdet = np.linalg.slogdet(total)
            if np.any(sign <= 0):
                value = np.inf
                break
            solved = np.linalg.solve(total, y[rows, :, None])[..., 0]
            chunk_value += .5*np.sum(weight[rows]*(logdet+np.einsum('ni,ni->n', y[rows], solved)))
        if AVAILABLE:
            adjusted = chunk_value-compensation
            updated = value+adjusted
            compensation = (updated-value)-adjusted
            value = updated
        else:
            # Preserve the original NumPy operation sequence without Numba.
            value += chunk_value
    result = value/weight.sum()
    if return_info:
        return result, dict(n_matrices=len(y), n_fallback=n_fallback, compiled=AVAILABLE)
    return result
