"""Joint, sampling-error-corrected covariance estimation for snipar effects.

The zero-mean working model is E[y_j y_j.T] = Sigma + S_j. Gaussian
product covariances provide feasible GLS weights; inverse LD scores multiply
each SNP's estimating equation. A PSD-constrained Fisher-scoring fit solves
the weighted Gaussian covariance score, with likelihood backtracking.
No empirical centering, LD deconvolution, or SE calibration is performed.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from . import coherent_kernels as fast_kernels


def _pairs(p):
    return np.triu_indices(p)


def _matrix(theta, p):
    out = np.zeros((p, p), dtype=float)
    i, j = _pairs(p)
    out[i, j] = out[j, i] = theta
    return out


def _project_psd(matrix):
    values, vectors = np.linalg.eigh((matrix + matrix.T) / 2)
    return (vectors * np.maximum(values, 0)) @ vectors.T


def gaussian_product_covariance(total_covariance):
    """Cov(vech(yy')) by Isserlis, in numpy.triu_indices order.

    Accepts (p,p) or (...,p,p). Subtracting known sampling covariance from
    the products does not change this working covariance.
    """
    total = np.asarray(total_covariance, dtype=float)
    i, j = _pairs(total.shape[-1])
    a, b, c, d = i[:, None], j[:, None], i[None, :], j[None, :]
    return total[..., a, c] * total[..., b, d] + total[..., a, d] * total[..., b, c]


def _product_precision(total):
    """Inverse product covariance without an inversion in p(p+1)/2 space."""
    precision = np.linalg.inv(total)
    i, j = _pairs(total.shape[-1])
    multiplicity = np.where(i == j, 1.0, 2.0)
    return (gaussian_product_covariance(precision)
            * multiplicity[:, None] * multiplicity[None, :] / 4.0)


def psd_gls(information, rhs, p):
    """Minimize theta'I theta/2 - theta'rhs subject to Sigma(theta) >= 0.

    This is a projection in the GLS information metric. Eigenvalue clipping
    alone would optimize the wrong metric for heterogeneous sampling errors.
    """
    information = np.asarray(information, dtype=float)
    rhs = np.asarray(rhs, dtype=float)
    q = p * (p + 1) // 2
    if information.shape != (q, q) or rhs.shape != (q,):
        raise ValueError('GLS information/rhs have incompatible dimensions')
    information = (information + information.T) / 2
    np.linalg.cholesky(information)
    unconstrained = np.linalg.solve(information, rhs)
    raw = _matrix(unconstrained, p)
    if np.linalg.eigvalsh(raw)[0] >= 0:
        return raw
    # Scale both parameter and criterion: GWAS moments may be much smaller
    # than scipy's absolute optimizer tolerances.
    size = max(float(np.max(np.abs(unconstrained))), np.finfo(float).tiny)
    target = unconstrained / size
    info = information / np.linalg.eigvalsh(information)[-1]
    i, j = _pairs(p)
    start = _project_psd(raw / size)[i, j]

    def objective(x):
        residual = x - target
        return .5 * residual @ info @ residual, info @ residual

    def cone_jacobian(x):
        _, vectors = np.linalg.eigh(_matrix(x, p))
        return (vectors[i, :] * vectors[j, :] * np.where(i == j, 1, 2)[:, None]).T

    result = minimize(
        objective, start, jac=True, method='SLSQP',
        constraints={'type': 'ineq', 'fun': lambda x: np.linalg.eigvalsh(_matrix(x, p)),
                     'jac': cone_jacobian},
        options={'ftol': 1e-12, 'maxiter': 500},
    )
    # Eigenvalue constraints are nonsmooth at repeated zero eigenvalues.
    # SLSQP can report failure there, or stop on an objective tolerance before
    # the constrained score is accurate enough for an outer Fisher step.
    # Certify and, when necessary, refine its solution by projected accelerated
    # gradient in orthonormal svec coordinates. In these coordinates Euclidean
    # projection is exactly the Frobenius PSD projection, while the objective
    # retains the full GLS information metric.
    multiplicity = np.where(i == j, 1.0, np.sqrt(2.0))
    hessian = info / np.outer(multiplicity, multiplicity)
    spectrum = np.linalg.eigvalsh(hessian)
    lipschitz, curvature = spectrum[-1], spectrum[0]
    if curvature <= 0:
        raise RuntimeError('PSD GLS information is numerically singular')
    target_svec = multiplicity * target

    def project_svec(x):
        return _project_psd(_matrix(x / multiplicity, p))[i, j] * multiplicity

    def gradient_step(x):
        return project_svec(x - hessian @ (x - target_svec) / lipschitz)

    starting_point = result.x if np.all(np.isfinite(result.x)) else start
    current = project_svec(multiplicity * starting_point)
    # Strong convexity converts a gradient-mapping residual to an error bound.
    # The floor only accommodates floating-point projection roundoff.
    tolerance = (max(5e-14, 2e-10 * curvature / lipschitz)
                 * max(1.0, np.linalg.norm(target_svec)))
    acceleration = ((np.sqrt(lipschitz) - np.sqrt(curvature))
                    / (np.sqrt(lipschitz) + np.sqrt(curvature)))
    extrapolated = current.copy()
    for _ in range(20000):
        residual = gradient_step(current) - current
        if np.linalg.norm(residual) <= tolerance:
            return _project_psd(_matrix(current / multiplicity, p)) * size
        updated = gradient_step(extrapolated)
        extrapolated = updated + acceleration * (updated - current)
        current = updated
    raise RuntimeError('PSD GLS optimization failed to satisfy the constrained '
                       'score (projected step %.3g)' % np.linalg.norm(residual))


def _validate(effects, sampling_covariance, ldscores):
    y = np.asarray(np.ma.filled(effects, np.nan), dtype=float)
    s = np.asarray(np.ma.filled(sampling_covariance, np.nan), dtype=float)
    if y.ndim != 2 or y.shape[1] not in (2, 3) or y.shape[0] < 2:
        raise ValueError('effects must have shape (n,2) or (n,3), with at least two SNPs')
    n, p = y.shape
    if s.shape != (n, p, p):
        raise ValueError('sampling_covariance must have shape (n,p,p) matching effects')
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(s)):
        raise ValueError('effects and sampling covariances must be finite; filter complete SNP vectors first')
    magnitudes = np.max(np.abs(s), axis=(1, 2))
    if np.any(np.max(np.abs(s - s.transpose(0, 2, 1)), axis=(1, 2)) > 1e-10 * magnitudes):
        raise ValueError('sampling covariance matrices must be symmetric')
    s = (s + s.transpose(0, 2, 1)) / 2
    eigenvalues = np.linalg.eigvalsh(s)
    if np.any(eigenvalues[:, 0] <= 1e-12 * eigenvalues[:, -1]) or np.any(eigenvalues[:, -1] <= 0):
        raise ValueError('sampling covariance matrices must be positive definite and numerically nonsingular; use a nonredundant effect basis')
    ld = np.ones(n) if ldscores is None else np.asarray(np.ma.filled(ldscores, np.nan), dtype=float)
    if ld.shape != (n,) or not np.all(np.isfinite(ld)) or np.any(ld <= 0):
        raise ValueError('LD scores must be a finite, positive vector of length n')
    return y, s, ld


def _normal_equations(y, s, weight, sigma, chunk_size, block_ids=None):
    if block_ids is None and fast_kernels is not None and fast_kernels.AVAILABLE:
        return fast_kernels.normal_equations(y, s, weight, sigma, chunk_size)
    n, p = y.shape
    i, j = _pairs(p)
    q = len(i)
    groups = 1 if block_ids is None else int(block_ids.max()) + 1
    info = np.zeros((groups, q, q))
    rhs = np.zeros((groups, q))
    for start in range(0, n, chunk_size):
        sl = slice(start, min(start + chunk_size, n))
        products = y[sl, i] * y[sl, j] - s[sl, i, j]
        precision = _product_precision(s[sl] + sigma)
        precision *= weight[sl, None, None]
        contributions = np.einsum('nij,nj->ni', precision, products)
        if block_ids is None:
            info[0] += precision.sum(axis=0)
            rhs[0] += contributions.sum(axis=0)
        else:
            # Only block totals are retained, not n x q x q matrices.
            for block in np.unique(block_ids[sl]):
                mask = block_ids[sl] == block
                info[block] += precision[mask].sum(axis=0)
                rhs[block] += contributions[mask].sum(axis=0)
    return (info[0], rhs[0]) if block_ids is None else (info, rhs)


def _criterion(y, s, weight, sigma, chunk_size):
    if fast_kernels is not None and fast_kernels.AVAILABLE:
        return fast_kernels.criterion(y, s, weight, sigma, chunk_size)
    value = 0.0
    for start in range(0, len(y), chunk_size):
        sl = slice(start, min(start + chunk_size, len(y)))
        total = s[sl] + sigma
        sign, logdet = np.linalg.slogdet(total)
        if np.any(sign <= 0):
            return np.inf
        solved = np.linalg.solve(total, y[sl, :, None])[..., 0]
        value += .5 * np.sum(weight[sl] * (logdet + np.einsum('ni,ni->n', y[sl], solved)))
    return value / weight.sum()


def fit_covariance(effects, sampling_covariance, ldscores=None, max_iter=100,
                   tol=1e-7, chunk_size=20000, initial=None):
    """Fit a PSD signal covariance using iterated Gaussian product GLS.

    PSD Fisher steps with objective backtracking solve the constrained
    covariance score. Optimality refers to the common-covariance Gaussian
    working model and not to arbitrary LD or heterogeneous true covariances.
    Failure to converge raises instead of returning a partial fit.
    """
    y, s, ld = _validate(effects, sampling_covariance, ldscores)
    return _fit_covariance_validated(y, s, ld, max_iter=max_iter, tol=tol,
                                     chunk_size=chunk_size, initial=initial)


def _fit_covariance_validated(y, s, ld, max_iter=100, tol=1e-7,
                              chunk_size=20000, initial=None):
    """Private fit for arrays already checked by _validate; subsets stay valid."""
    if int(max_iter) != max_iter or max_iter < 1 or not np.isfinite(tol) or tol <= 0:
        raise ValueError('max_iter must be positive integer and tol must be finite and positive')
    if int(chunk_size) != chunk_size or chunk_size < 1:
        raise ValueError('chunk_size must be a positive integer')
    chunk_size = int(chunk_size)
    scale = np.sqrt(np.mean(y * y + np.diagonal(s, axis1=1, axis2=2), axis=0))
    outer_scale = np.outer(scale, scale)
    y = y / scale
    s = s / outer_scale
    weight = 1.0 / ld
    weight /= weight.mean()
    p = y.shape[1]
    if initial is None:
        raw = (np.einsum('n,ni,nj->ij', weight, y, y)
               - np.einsum('n,nij->ij', weight, s)) / weight.sum()
        sigma = _project_psd(raw)
    else:
        initial = np.asarray(initial, dtype=float)
        if initial.shape != (p, p) or not np.all(np.isfinite(initial)):
            raise ValueError('initial must be a finite (p,p) covariance matrix')
        sigma = _project_psd(initial / outer_scale)
    criterion = _criterion(y, s, weight, sigma, chunk_size)
    for iteration in range(1, int(max_iter) + 1):
        info, rhs = _normal_equations(y, s, weight, sigma, chunk_size)
        candidate = psd_gls(info, rhs, p)
        change = np.linalg.norm(candidate - sigma) / max(np.linalg.norm(sigma), 1e-5)
        if change <= tol:
            # Current iterate already satisfies the constrained score.
            break
        step = 1.0
        while step >= 2.0 ** -24:
            trial = sigma + step * (candidate - sigma)
            trial_criterion = _criterion(y, s, weight, trial, chunk_size)
            if trial_criterion <= criterion + 1e-13 * max(1.0, abs(criterion)):
                sigma, criterion = trial, trial_criterion
                break
            step *= .5
        else:
            raise RuntimeError('Covariance fit stalled before satisfying the constrained score')
    else:
        raise RuntimeError('Covariance fit did not converge within %d iterations (relative score step %.3g)' % (max_iter, change))
    covariance = _project_psd(sigma) * outer_scale
    eigenvalues = np.linalg.eigvalsh(sigma)
    return covariance, {
        'iterations': iteration, 'converged': True, 'relative_score_step': float(change),
        'linear_algebra_backend': 'numba_fused' if fast_kernels is not None and fast_kernels.AVAILABLE else 'numpy',
        'scaled_objective': float(criterion),
        'boundary': bool(eigenvalues[0] <= 1e-8 * max(eigenvalues[-1], 1e-12)),
        'min_scaled_eigenvalue': float(eigenvalues[0]),
    }


def _moment_names(effect_names):
    return [('var_' + effect_names[i] if i == j else
             'cov_' + effect_names[i] + '_' + effect_names[j])
            for i, j in zip(*_pairs(len(effect_names)))]


def _derived(covariance):
    p = len(covariance)
    transform = (np.array([[1, 0], [0, 1], [1, 1]]) if p == 2 else
                 np.array([[1, 0, 0], [0, .5, .5], [1, .5, .5]]))
    expanded = transform @ covariance @ transform.T
    names = ['direct', 'avg_NTC', 'population']
    keys = _moment_names(names)
    values = list(expanded[_pairs(3)])
    variance = np.diag(expanded)
    # Treat numerical cancellation at a zero-variance boundary as zero.
    threshold = 1e-10 * max(float(np.max(np.abs(expanded))), np.finfo(float).tiny)
    positive = variance > threshold
    for a, b in [(0, 2), (0, 1)]:
        keys.append('r_' + names[a] + '_' + names[b])
        value = expanded[a, b] / np.sqrt(variance[a] * variance[b]) if positive[a] and positive[b] else np.nan
        values.append(float(np.clip(value, -1, 1)))
    # Both directions make the numerator/denominator convention explicit.
    for a, b in [(0, 2), (1, 2), (0, 1), (2, 0), (2, 1), (1, 0)]:
        keys.append('sqrt_var_ratio_' + names[a] + '_' + names[b])
        values.append(np.sqrt(max(variance[a], 0) / variance[b]) if positive[b] else np.nan)
    return keys, np.asarray(values), expanded


def _jackknife_covariance(delete):
    """No nan-dropping: an undefined deletion makes the affected row/col NA."""
    b, k = delete.shape
    covariance = np.full((k, k), np.nan)
    valid = np.all(np.isfinite(delete), axis=0)
    delta = delete[:, valid] - delete[:, valid].mean(axis=0)
    covariance[np.ix_(valid, valid)] = (b - 1) / b * delta.T @ delta
    return covariance


def fit_coherent(effects, sampling_covariance, ldscores=None, n_blocks=200,
                 blocks=None, effect_names=None, max_iter=100, tol=1e-7,
                 chunk_size=20000, jackknife='refit'):
    """Estimate moments, correlations, SD ratios and their joint jackknife SEs.

    Columns must be [D, average NTC] or [D, paternal NTC, maternal NTC].
    Rows must be genomically ordered if explicit block labels are not supplied.
    ``refit`` refits all weights and constraints per deletion. ``frozen`` is a
    faster conditional approximation retaining full-fit product weights; it
    is not the full refitted jackknife. Boundary SEs are descriptive and do
    not justify ordinary Wald intervals for zero variances or |r|=1.
    """
    y, s, ld = _validate(effects, sampling_covariance, ldscores)
    n, p = y.shape
    if jackknife not in ('refit', 'frozen'):
        raise ValueError('jackknife must be refit or frozen')
    if blocks is None:
        if int(n_blocks) != n_blocks or not 2 <= n_blocks <= n:
            raise ValueError('n_blocks must be an integer between 2 and the number of retained SNPs')
        n_blocks = int(n_blocks)
        # Array-split sized contiguous blocks, with at most one SNP difference.
        block_ids = np.repeat(np.arange(n_blocks), [len(x) for x in np.array_split(np.arange(n), n_blocks)])
        block_labels = np.arange(n_blocks)
    else:
        labels = np.asarray(blocks)
        if labels.shape != (n,) or any(x is None or str(x).lower() in ('nan', 'nat', '') for x in labels):
            raise ValueError('blocks must supply one nonmissing label per retained SNP')
        block_labels, block_ids = np.unique(labels, return_inverse=True)
    b = len(block_labels)
    if b < 2 or n - np.bincount(block_ids).max() < 2:
        raise ValueError('At least two blocks are needed, each leaving at least two SNPs after deletion')
    if effect_names is None:
        effect_names = ['direct', 'avg_NTC'] if p == 2 else ['direct', 'paternal', 'maternal']
    if len(effect_names) != p or len(set(effect_names)) != p:
        raise ValueError('effect_names must contain p distinct names in direct/NTC order')
    options = dict(max_iter=max_iter, tol=tol, chunk_size=chunk_size)
    covariance, diagnostics = _fit_covariance_validated(y, s, ld, **options)
    chunk_size = int(chunk_size)
    moment_names = _moment_names(effect_names)
    metric_names, estimates, expanded = _derived(covariance)
    moment_delete = np.empty((b, len(moment_names)))
    delete_values = np.empty((b, len(metric_names)))
    boundary_deletes = 0
    delete_iterations = []
    if jackknife == 'frozen':
        # Use a fixed diagonal change of units for well-conditioned cone solves.
        scale = np.sqrt(np.mean(y * y + np.diagonal(s, axis1=1, axis2=2), axis=0))
        outer_scale = np.outer(scale, scale)
        block_info, block_rhs = _normal_equations(y / scale, s / outer_scale,
            1 / ld, covariance / outer_scale, chunk_size, block_ids)
        total_info, total_rhs = block_info.sum(axis=0), block_rhs.sum(axis=0)
    for block in range(b):
        if jackknife == 'refit':
            keep = block_ids != block
            try:
                deleted, detail = _fit_covariance_validated(y[keep], s[keep], ld[keep], initial=covariance, **options)
            except RuntimeError as exc:
                raise RuntimeError('Jackknife block %s failed: %s' % (block_labels[block], exc)) from exc
            delete_iterations.append(detail['iterations'])
            boundary_deletes += int(detail['boundary'])
        else:
            deleted = psd_gls(total_info - block_info[block], total_rhs - block_rhs[block], p) * outer_scale
            eig = np.linalg.eigvalsh(deleted / outer_scale)
            boundary_deletes += int(eig[0] <= 1e-8 * max(eig[-1], 1e-12))
        moment_delete[block] = deleted[_pairs(p)]
        delete_values[block] = _derived(deleted)[1]
    moment_covariance = _jackknife_covariance(moment_delete)
    metric_covariance = _jackknife_covariance(delete_values)
    undefined = ~np.isfinite(estimates)
    metric_covariance[undefined, :] = metric_covariance[:, undefined] = np.nan
    # A separate finite-panel diagnostic uses ONE scalar weight per SNP for
    # every matrix entry. Unlike the model-based GLS fit, it is not constrained
    # or reweighted using the estimated covariance, and can be indefinite.
    scalar = 1 / ld
    scalar /= scalar.sum()
    i, j = _pairs(p)
    common_totals = np.zeros((b, len(i)))
    common_mass = np.bincount(block_ids, weights=scalar, minlength=b)
    for start in range(0, n, chunk_size):
        sl = slice(start, min(start + chunk_size, n))
        products = (y[sl, i] * y[sl, j] - s[sl, i, j]) * scalar[sl, None]
        for block in np.unique(block_ids[sl]):
            common_totals[block] += products[block_ids[sl] == block].sum(axis=0)
    common_moments = common_totals.sum(axis=0)
    common_delete = (common_moments - common_totals) / (1 - common_mass[:, None])
    common_sampling_covariance = _jackknife_covariance(common_delete)
    common_covariance = _matrix(common_moments, p)
    common_eigenvalues = np.linalg.eigvalsh(common_covariance)
    block_sizes = np.bincount(block_ids)
    diagnostics.update({
        'n_snps': n, 'n_blocks': b, 'block_labels': [str(x) for x in block_labels],
        'block_sizes': block_sizes.tolist(), 'jackknife': jackknife,
        'jackknife_boundary_count': boundary_deletes,
        'jackknife_max_iterations': max(delete_iterations) if delete_iterations else None,
        'common_weight_min_eigenvalue': float(common_eigenvalues[0]),
        'common_weight_is_psd': bool(common_eigenvalues[0] >= 0),
        'undefined_se_metrics': [metric_names[k] for k in range(len(metric_names))
                                 if not np.isfinite(metric_covariance[k, k])],
        'warnings': [],
    })
    warnings = diagnostics['warnings']
    if diagnostics['boundary'] or boundary_deletes:
        warnings.append('PSD boundary encountered; jackknife SEs are descriptive, ordinary Wald intervals can be unreliable.')
    if jackknife == 'frozen':
        warnings.append('Frozen full-fit weights: SEs are conditional approximations, not a full refitted jackknife.')
    if block_sizes.max() > 2 * block_sizes.min():
        warnings.append('Highly unequal block sizes: the ordinary equal-block jackknife may be poorly calibrated.')
    if diagnostics['undefined_se_metrics']:
        warnings.append('An undefined full-sample or deletion statistic has SE=NA; invalid deletions were not dropped.')
    return dict(
        effect_names=list(effect_names), covariance=covariance,
        moment_names=moment_names, moments=covariance[_pairs(p)],
        moment_covariance=moment_covariance, moment_se=np.sqrt(np.maximum(np.diag(moment_covariance), 0)),
        moment_delete=moment_delete, metric_names=metric_names, estimates=estimates,
        se=np.sqrt(np.maximum(np.diag(metric_covariance), 0)), delete_values=delete_values,
        metric_covariance=metric_covariance, expanded_covariance=expanded,
        common_weight_covariance=common_covariance, common_weight_moments=common_moments,
        common_weight_moment_covariance=common_sampling_covariance,
        common_weight_moment_se=np.sqrt(np.maximum(np.diag(common_sampling_covariance), 0)),
        common_weight_moment_delete=common_delete,
        expanded_names=['direct', 'avg_NTC', 'population'], diagnostics=diagnostics,
    )
