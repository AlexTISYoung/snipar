from operator import attrgetter
from typing import Tuple, NamedTuple, Sequence
from typing import Literal
from itertools import product

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from numpy.linalg import slogdet, solve
from scipy.linalg import cho_factor, cho_solve, qr
from scipy.sparse import csc_matrix, tril
from scipy.sparse.linalg import splu, SuperLU, cg
from numba import njit
from snipar.types import FamLabels, SparseGRMRepr, Ids, IdDict
from snipar.gtarray import gtarray




def _solve_snp_effects(information, scores):
    """Solve for one coefficient vector per SNP using an explicit RHS column."""
    return solve(information, scores[..., None])[..., 0]


def build_sib_arr(fam_labels: FamLabels) -> SparseGRMRepr:    
    """Build lower-triangular nonzero entries of sibship matrix.

    Args:
        fam_labels (FamLabels): ndarray of family id strings corresponding to each individual.

    Returns:
        SparseGRMRepr: sparse GRM representation.
    """
    partition = _FamilyPartition(fam_labels)
    entries = (partition.sizes * (partition.sizes + 1) // 2).sum()
    rows, columns = _family_coordinates(partition.rows, partition.offsets, entries)
    return np.ones(entries, dtype=np.uint8), rows, columns


def build_ibdrel_arr(ibdrel_path: str, id_dict: IdDict,
                     keep: Ids,
                     thres: float = 0.05) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build ibd relatedness array (lower triangular entries) and corresponding indices from KING ibdseg output.

    Args:
        ibdrel_path (str): Path to ibdseg output.
        id_dict (IdDict): dictionary of id-index pairs.
        keep (Ids): list of ids to keep.
        ignore_sib (bool): whether to set sibling entries to 0.
        thres (float, optional): sparsity threshold. Defaults to 0.0205.

    Returns:
        SparseGRMRepr: sparse GRM representation.
    """
    with open(ibdrel_path) as grm:
        header = grm.readline().split()
    if not {'ID1', 'ID2'}.issubset(header):
        raise ValueError('GRM should have ID1 and ID2 columns.')
    if 'PropIBD' in header:
        relcol = 'PropIBD'
    elif 'relatedness' in header:
        relcol = 'relatedness'
    else:
        raise ValueError('GRM should have PropIBD or relatedness column.')
    cohort = pd.Index(keep)
    row_for_sample = np.fromiter((id_dict[iid] for iid in keep), dtype=np.uint32, count=len(keep))
    values, rows, cols = [], [], []
    with pd.read_csv(ibdrel_path, sep=r'\s+', usecols=['ID1', 'ID2', relcol],
                     dtype={'ID1': str, 'ID2': str, relcol: float},
                     keep_default_na=False, chunksize=250_000) as chunks:
        for table in chunks:
            r = table[relcol].to_numpy()
            relevant = r > thres
            r = r[relevant]
            left = cohort.get_indexer(table.loc[relevant, 'ID1'])
            right = cohort.get_indexer(table.loc[relevant, 'ID2'])
            include = (left >= 0) & (right >= 0)
            invalid = np.flatnonzero(include & (r > 0.8))
            if invalid.size:
                result = table.loc[relevant].iloc[invalid[0]]
                raise ValueError(f"Impossible ibd relatedness {result[relcol]} for pair "
                                 f"({result['ID1']}, {result['ID2']}). MZ pairs should be removed.")
            left, right = row_for_sample[left[include]], row_for_sample[right[include]]
            if np.any(left == right):
                raise ValueError('GRM contains a self-pair; diagonal entries are added automatically.')
            values.append(r[include])
            rows.append(np.maximum(left, right))
            cols.append(np.minimum(left, right))
    return (np.concatenate((*values, np.ones(len(keep)))),
            np.concatenate((*rows, row_for_sample)), np.concatenate((*cols, row_for_sample)))


def build_grm_arr(grm_path: str, id_dict: IdDict,
                  thres: float) -> SparseGRMRepr:
    """Build GRM data array and corresponding indices.
    """
    grm_ids = np.loadtxt(f'{grm_path}.grm.id', dtype=str, usecols=1, ndmin=1)
    # GRM indices are one-based. Preserve -1 as the marker for off-cohort IDs.
    index = np.full(grm_ids.size + 1, -1, dtype=np.int64)
    index[1:] = [id_dict.get(iid, -1) for iid in grm_ids]
    arrays = ([], [], [])
    with pd.read_csv(f'{grm_path}.grm.gz', sep=r'\s+', header=None, usecols=[0, 1, 3],
                     dtype={0: np.int64, 1: np.int64, 3: np.float64}, chunksize=250_000) as chunks:
        for chunk in chunks:
            relevant = chunk[3].to_numpy() >= thres
            ids = chunk.loc[relevant, [0, 1]].to_numpy()
            if np.any(ids < 1) or np.any(ids >= len(index)):
                raise ValueError('GCTA GRM index lies outside the .grm.id file')
            row, col = index[ids].T
            valid = (row >= 0) & (col >= 0)
            arrays[0].append(chunk[3].to_numpy()[relevant][valid])
            arrays[1].append(np.maximum(row[valid], col[valid]).astype('uint32'))
            arrays[2].append(np.minimum(row[valid], col[valid]).astype('uint32'))
    return tuple(np.concatenate(chunks) if chunks else np.empty(0, dtype)
                 for chunks, dtype in zip(arrays, ('float64', 'uint32', 'uint32')))


def match_grm_ids(ids: Ids,
                  fam_labels: FamLabels,
                  grm_path: str,
                  grm_source: Literal['gcta', 'ibdrel'],
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """Match ids with GRM individual ids.
    """
    if grm_source == 'gcta':
        grm_ids = np.loadtxt(f'{grm_path}.grm.id', dtype=str, usecols=1, ndmin=1)
    elif grm_source == 'ibdrel':
        pairs = pd.read_csv(grm_path, sep=r'\s+', usecols=['ID1', 'ID2'],
                            dtype=str, keep_default_na=False)
        grm_ids = np.unique(pairs.to_numpy().ravel())
    else:
        raise ValueError('Incorrect source.')
    common = np.isin(ids, grm_ids)
    if not np.any(common):
        raise ValueError('No analysis individuals are present in the GRM.')
    return ids[common], fam_labels[common]


class GradHessComponents(NamedTuple):
    P_y: np.ndarray
    P_varcomp_mats: Tuple[np.ndarray, ...]


def cached_property_depends_on(*args):
    attrs = attrgetter(*args)

    def decorator(func):
        name = f'_cached_{func.__name__}'

        def get(self):
            key = attrs(self)
            cached = self.__dict__.get(name)
            if cached is None or cached[0] != key:
                cached = self.__dict__[name] = (key, func(self))
            return cached[1]

        return property(get, doc=func.__doc__)

    return decorator


def _precision_terms(lu, y, covariates):
    """Precompute the covariate projection in inverse covariance geometry."""
    precision_y = lu.solve(np.asarray(y, dtype=float))
    if covariates is None or not covariates.shape[1]:
        return precision_y, None
    covariates = np.asarray(covariates, dtype=float)
    precision_z = lu.solve(covariates)
    # A covariate may be constant or absent within a particular sample group.
    projected_z = precision_z @ np.linalg.pinv(covariates.T @ precision_z, hermitian=True)
    return precision_y - projected_z @ (covariates.T @ precision_y), (covariates, projected_z)


class _FamilyPartition:
    def __init__(self, labels):
        labels = np.asarray(labels)
        if labels.ndim != 1:
            raise ValueError('Family labels must be one-dimensional')
        self.codes, families = pd.factorize(labels, sort=False)
        if np.any(self.codes < 0):
            raise ValueError('Family labels must be nonmissing')
        self.sizes = np.bincount(self.codes, minlength=len(families))
        self.rows = np.argsort(self.codes, kind='stable')
        self.offsets = np.concatenate(([0], np.cumsum(self.sizes)))


@njit(cache=True, nogil=True)
def _family_coordinates(rows, offsets, entries):
    row_ind, col_ind = np.empty(entries, np.uint32), np.empty(entries, np.uint32)
    entry = 0
    for family in range(len(offsets) - 1):
        start, end = offsets[family], offsets[family + 1]
        for i in range(start, end):
            for j in range(i, end):
                row_ind[entry], col_ind[entry] = rows[j], rows[i]
                entry += 1
    return row_ind, col_ind


@njit(cache=True, nogil=True)
def _group_columns(values, codes, groups):
    sums = np.zeros((groups, values.shape[1]))
    for row in range(len(codes)):
        code = codes[row]
        if code >= 0:
            for j in range(values.shape[1]):
                sums[code, j] += values[row, j]
    return sums


@njit(cache=True, nogil=True)
def _solve_family_covariance(rhs, rows, offsets, within_precision, mean_precision):
    result = np.empty(rhs.shape)
    for f in range(len(mean_precision)):
        start, end = offsets[f], offsets[f + 1]
        if end - start == 1:
            row = rows[start]
            for column in range(rhs.shape[1]):
                result[row, column] = rhs[row, column] * mean_precision[f]
            continue
        for column in range(rhs.shape[1]):
            total = 0.
            for j in range(start, end):
                total += rhs[rows[j], column]
            mean = total / (end - start)
            for j in range(start, end):
                row = rows[j]
                result[row, column] = (rhs[row, column] - mean) * within_precision + mean * mean_precision[f]
    return result


class _FamilyFactor:
    """Exact inverse for a covariance with one shared effect per family."""

    def __init__(self, partition, sibling, residual):
        self.partition = partition
        self.mean_variance = residual + sibling * partition.sizes
        self.within_count = len(partition.rows) - len(partition.sizes)
        self.residual = residual
        if (self.within_count and residual <= 0) or np.any(self.mean_variance <= 0):
            raise RuntimeError('Family covariance is singular')
        self.within_precision = 1 / residual if self.within_count else 0.
        self.mean_precision = 1 / self.mean_variance

    def solve(self, rhs):
        rhs = np.asarray(rhs, dtype=float)
        if rhs.ndim == 1:
            return self.solve(rhs[:, None])[:, 0]
        return _solve_family_covariance(rhs, self.partition.rows, self.partition.offsets,
                                        self.within_precision, self.mean_precision)

    @property
    def logdet(self):
        within = self.within_count * np.log(self.residual) if self.within_count else 0.
        return np.log(self.mean_variance).sum() + within


@njit(cache=True, nogil=True)
def _family_information(gts, rows, offsets, within_precision, mean_precision,
                        precision_y, precision_z):
    """Accumulate small SNP information matrices without solving for every genotype."""
    _, p, n_loci = gts.shape
    n_covariates = precision_z.shape[1]
    information = np.zeros((n_loci, p, p))
    scores = np.zeros((n_loci, p))
    cross = np.zeros((n_loci, p, n_covariates))
    family_sum = np.empty((p, n_loci))
    for family in range(len(mean_precision)):
        start, end = offsets[family], offsets[family + 1]
        size = end - start
        family_sum[:] = 0.
        for j in range(start, end):
            row = rows[j]
            for a in range(p):
                for locus in range(n_loci):
                    value = np.float64(gts[row, a, locus])
                    family_sum[a, locus] += value
                    scores[locus, a] += value * precision_y[row]
                    for c in range(n_covariates):
                        cross[locus, a, c] += value * precision_z[row, c]
        if size > 1:
            family_sum /= size
        for a in range(p):
            for b in range(a + 1):
                for locus in range(n_loci):
                    mean_a, mean_b = family_sum[a, locus], family_sum[b, locus]
                    value = size * mean_precision[family] * mean_a * mean_b
                    if size > 1:
                        within = 0.
                        for j in range(start, end):
                            row = rows[j]
                            within += ((gts[row, a, locus] - mean_a)
                                       * (gts[row, b, locus] - mean_b))
                        value += within_precision * within
                    information[locus, a, b] += value
    for a in range(p):
        for b in range(a):
            information[:, b, a] = information[:, a, b]
    return information, scores, cross


def _identified_information(information, unadjusted_variance, n):
    # Rank in correlation units, independent of predictor scale.
    k = information.shape[1]
    diagonal = np.diagonal(information, axis1=1, axis2=2)
    tolerance = k * np.finfo(float).eps
    projection_tolerance = max(n, k, 16) * np.finfo(float).eps
    valid = np.flatnonzero(np.all(np.isfinite(information), axis=(1, 2)) &
                           np.all(diagonal > projection_tolerance * unadjusted_variance, axis=1))
    scales = np.sqrt(diagonal[valid])
    normalized = information[valid] / scales[:, :, None] / scales[:, None, :]
    identified = np.linalg.eigvalsh(normalized)[:, 0] > tolerance
    scales = scales[identified]
    covariance = np.linalg.inv(normalized[identified]) / scales[:, :, None] / scales[:, None, :]
    return valid[identified], covariance


def _fit_gls_batch(gts, lu, precision_y, projection, return_weights=False):
    """Fit individual SNPs using shared covariance and isolate unidentifiable loci."""
    n, k, loci = gts.shape
    alpha = np.full((loci, k), np.nan)
    covariance = np.full((loci, k, k), np.nan)
    weights = np.full((n, loci), np.nan) if return_weights else None
    if n < k or not loci:
        return alpha, covariance, weights
    # Limit the work arrays for the multi-right-hand-side SuperLU solve.
    block_size = max(1, min(256, (32 * 1024 ** 2) // (8 * n * k)))
    for start in range(0, loci, block_size):
        block = gts[..., start:start + block_size]
        # Fill the solver's column-major RHS directly, including strided SNP slices.
        x = np.empty((n, k * block.shape[2]), dtype=float, order='F')
        np.copyto(x.reshape(block.shape), block)
        precision_x = lu.solve(x)
        unadjusted_variance = np.einsum('ij,ij->j', x, precision_x).reshape(k, -1).T
        if projection is not None:
            z, projected_z = projection
            precision_x -= projected_z @ (z.T @ precision_x)
        x = x.reshape(n, k, -1).transpose(2, 1, 0)
        precision_x = precision_x.reshape(n, k, -1).transpose(2, 0, 1)
        information = x @ precision_x
        information = (information + information.transpose(0, 2, 1)) * 0.5
        scores = x @ precision_y
        valid, vc = _identified_information(information, unadjusted_variance, n)
        covariance[start + valid] = vc
        alpha[start + valid] = (vc @ scores[valid, :, None])[..., 0]
        if return_weights:
            weights[:, start + valid] = np.einsum('lnk,lk->nl', precision_x[valid], vc[:, :, 0])
    return alpha, covariance, weights


def _fit_family_gls(gts, factor, precision_y, precision_z, covariate_inverse):
    n, p, loci = gts.shape
    alpha = np.full((loci, p), np.nan)
    covariance = np.full((loci, p, p), np.nan)
    if n < p or not loci:
        return alpha, covariance
    partition = factor.partition
    information, scores, cross = _family_information(
        gts, partition.rows, partition.offsets, factor.within_precision,
        factor.mean_precision, precision_y, precision_z,
    )
    unadjusted = np.diagonal(information, axis1=1, axis2=2).copy()
    if precision_z.shape[1]:
        information -= cross @ covariate_inverse @ cross.transpose(0, 2, 1)
    valid, fitted = _identified_information(information, unadjusted, n)
    covariance[valid] = fitted
    alpha[valid] = (fitted @ scores[valid, :, None])[:, :, 0]
    return alpha, covariance


def _combine_direct_effects(effects, covariance):
    """Minimum-variance combination of direct-effect estimates with a common target."""
    try:
        weights = solve(covariance, np.ones(len(effects)))
    except np.linalg.LinAlgError:
        return np.nan, np.nan
    total = weights.sum()
    if total <= 0 or not np.isfinite(total):
        return np.nan, np.nan
    return weights @ effects / total, 1 / total


class LinearMixedModel:
    """
    Wrapper of data and functions that compute estimates of variance components or SNP effects.
    Core functions: fit_snp_eff, robust_est, sib_diff_est
    """
    # logger = logger.getChild(__qualname__)


    def __init__(self, y: np.ndarray, 
                 varcomp_arr_lst: Sequence[Tuple[np.ndarray, np.ndarray, np.ndarray]],
                 varcomps: Sequence[float] = None,
                 covar_X: np.ndarray = None,
                 add_intercept: bool = False,
                 add_jitter: bool = False,
                 drop_collinear_covariates: bool = False, *,
                 family_labels: np.ndarray | None = None) -> None:
        """Initialize a LinearMixedModel instance.

        Args:
            y (np.ndarray): 1-d array phenotpye
            varcomp_arr_lst (Sequence[np.ndarray, ...]): a Sequence of arbitrary length, holding variance component matrices, excluding the residual variance matrix.
            varcomps (Sequence[float, ...]): a sequence of variance component parameters.
            covar_X (np.ndarray, optional): 2-d array holding fixed covariates. Defaults to None.
            add_intercept: (bool): a boolean indicating whether to add intercept or not.
            add_jitter: (bool): a boolean indicating whether to add jitter to the diagonal of the covariance matrix or not.
            drop_collinear_covariates: (bool): reduce nuisance covariates to independent columns before fitting.
            family_labels: represent the family covariance implicitly when no
                sparse components are provided.

        Raises:
            ValueError: y should be a 1-d array.
            ValueError: covar_X should be a 2-d array.
            ValueError: Dimensions of y and covar_X do not match.
        """
        if y.ndim != 1:
            raise ValueError('y should be a 1-d array.')
        self.y = np.asarray(y, dtype=float)
        self.n = len(y)
        self.has_covar = add_intercept or covar_X is not None
        if covar_X is not None:
            covar_X = np.asarray(covar_X, dtype=float)
            # logger.info('has covar')
            if covar_X.ndim == 1:
                covar_X = covar_X[:, None]
            if covar_X.ndim != 2:
                raise ValueError('covar_X should be a 1-d or 2-d array.')
            if covar_X.shape[0] != y.shape[0]:
                raise ValueError(
                    f'Dimensions of y and covar_X do not match ({y.shape} and {covar_X.shape}).')
            self.Z  = np.hstack((np.ones((self.n, 1),dtype=covar_X.dtype), covar_X)) if add_intercept else covar_X
            if drop_collinear_covariates and self.Z.shape[1]:
                _, triangular, pivot = qr(self.Z, mode='economic', pivoting=True)
                diagonal = np.abs(np.diag(triangular))
                tolerance = diagonal.max(initial=0) * max(self.Z.shape) * np.finfo(float).eps
                independent = np.sort(pivot[:np.count_nonzero(diagonal > tolerance)])
                self.Z = self.Z[:, independent]
            self.has_covar = self.Z.shape[1] > 0
        elif add_intercept:
            self.Z = np.ones((self.n, 1))
        # self.logger.info(f'#individuals: {self.n}.')
            
        if family_labels is not None and len(family_labels) != self.n:
            raise ValueError('Provide one family label per sample')
        self._family = None
        if family_labels is not None and not len(varcomp_arr_lst):
            self._family = _FamilyPartition(family_labels)
            self._component_mats = None
            self.n_varcomps = 2
            if not self.has_covar:
                self._family_y_sum = np.bincount(self._family.codes, weights=self.y)
                group_means = self._family_y_sum / self._family.sizes
                centered = self.y - group_means[self._family.codes]
                self._family_y_within = centered @ centered
        else:
            self._component_mats = self._build_components(varcomp_arr_lst)
            self.n_varcomps = len(varcomp_arr_lst) + 1
        self.y_var  = self.y.var()
        # self.logger.info(f'Phenotypic variance: {self.y_var}.')
        self.jitter = self.y_var / 1000.0 if add_jitter else 0.0
        if varcomps is not None:
            if len(varcomps) != self.n_varcomps:
                raise ValueError('varcomps and varcomps_arr_lst must have the same length.')
            # _varcomps should be tuple for cached_dependent_property to work
            self._varcomps = tuple(float(value) for value in varcomps)
            if not all(np.isfinite(value) and value >= 0 for value in self._varcomps):
                raise ValueError('Variance components must be finite and nonnegative.')
            self.optimized = True
        else:
            self._varcomps = tuple(
                # self.y.var() if i == self.n_varcomps - 1 else 0. for i in range(self.n_varcomps))
                self.y_var / self.n_varcomps for i in range(self.n_varcomps))
            # self.logger.info(f'init varcomps: {self._varcomps}.')
            self.optimized = False

    @staticmethod
    def _empty_solver_result(sps_mat: csc_matrix, dense_mat: np.ndarray) -> np.ndarray:
        if dense_mat.ndim != 3:
            raise ValueError('dense_mat must be a 3-d array.')
        n1, n2 = sps_mat.shape
        m, n, c = dense_mat.shape
        if n != n2:
            raise ValueError(f'Input dims do not match: {n2} and {n}.')
        out = np.empty((m, n1, c))
        return out
    
    @staticmethod
    def sp_solve_dense3d(sps_mat: csc_matrix,
                         dense_mat: np.ndarray) -> np.ndarray:
        """Compute product of the inverse of given sparse matrix and given 3-d dense array; used for repeated OLS.

        Args:
            sps_mat (csc_matrix): 2-d sparse matrix.
            dense_mat (np.ndarray): 3-d dense array.

        Raises:
            ValueError: dense_mat should be 3-d.
            ValueError: 2nd dimensions of both inputs should match.

        Returns:
            np.ndarray: 3-d dense array.
        """
        return __class__.sp_solve_dense3d_lu(sps_mat, dense_mat)

    @staticmethod
    def sp_solve_dense3d_lu(sps_mat: csc_matrix,
                            dense_mat: np.ndarray) -> np.ndarray:
        """Compute product of the inverse of given sparse matrix and given 3-d dense array using LU; used for repeated OLS.

        Args:
            sps_mat (csc_matrix): 2-d sparse matrix.
            dense_mat (np.ndarray): 3-d dense array.

        Raises:
            ValueError: dense_mat should be 3-d.
            ValueError: 2nd dimensions of both inputs should match.

        Returns:
            np.ndarray: 3-d dense array.
        """
        out = __class__._empty_solver_result(sps_mat, dense_mat)
        lu = splu(sps_mat)
        for i in range(dense_mat.shape[0]):
            out[i, :, :] = lu.solve(np.asarray(dense_mat[i], dtype=float))
        return out
    
    @staticmethod
    def sp_mul_dense3d(sps_mat: csc_matrix,
                       dense_mat: np.ndarray) -> np.ndarray:
        """Compute product of given sparse matrix and given 3-d dense array; used for repeated OLS.

        Args:
            sps_mat (csc_matrix): 2-d sparse matrix.
            dense_mat (np.ndarray): 3-d dense array.

        Raises:
            ValueError: dense_mat should be 3-d.
            ValueError: 2nd dimensions of both inputs should match.

        Returns:
            np.ndarray: 3-d dense array.
        """
        out = __class__._empty_solver_result(sps_mat, dense_mat)
        for i in range(dense_mat.shape[0]):
            out[i, :, :] = sps_mat.dot(dense_mat[i, :, :])
        return out
    
    def _build_sym_mat(self, data: np.ndarray, row_ind: np.ndarray, col_ind: np.ndarray) -> csc_matrix:
        """Build sparse matrix using given lower triangular entries.

        Args:
            data (np.ndarray): data array.
            row_ind (np.ndarray): row indices.
            col_ind (np.ndarray): column indices.

        Returns:
            csc_matrix: symmetric sparse matrix in csc format.
        """
        tril_mat = csc_matrix((data, (row_ind, col_ind)), shape=(self.n, self.n))
        triu_mat = tril(tril_mat, k=-1, format='csc')
        return tril_mat + triu_mat.T

    def _build_components(self, components):
        return tuple(self._build_sym_mat(data, row_ind, col_ind)
                     for data, row_ind, col_ind in components) + (
            csc_matrix((np.ones(self.n), (np.arange(self.n), np.arange(self.n))),
                       shape=(self.n, self.n)),
        )

    @property
    def _varcomp_mats(self):
        if self._component_mats is None:
            partition = self._family
            entries = (partition.sizes * (partition.sizes + 1) // 2).sum()
            rows, cols = _family_coordinates(partition.rows, partition.offsets, entries)
            self._component_mats = self._build_components(((np.ones(entries), rows, cols),))
        return self._component_mats

    @property
    def varcomps(self) -> Tuple[float, ...]:
        """Return optimized variance components.

        Raises:
            ValueError: self.optimized should be True

        Returns:
            Tuple[float, ...]: a tuple of variance components
        """
        if not self.optimized:
            raise ValueError('Variance components are not optimized.')
        return self._varcomps

    @cached_property_depends_on('_varcomps')
    def V(self) -> csc_matrix:
        """Compute V.

        Returns:
            csc_matrix: V in sparse csc format
        """
        # self.logger.debug('Calculating V...')
        varcomps = np.array(self._varcomps, dtype=float)
        varcomps[-1] += self.jitter
        V: csc_matrix = sum(
            varcomps[i] * self._varcomp_mats[i] for i in range(self.n_varcomps)
        )
        return V

    @cached_property_depends_on('_varcomps')
    def V_lu(self) -> SuperLU | _FamilyFactor:
        """Compute sparse LU factorization of V.

        Returns:
            SuperLU: wrapper object holding LU factorization of V
        """
        if self._family is not None:
            sibling, residual = self._varcomps
            return _FamilyFactor(self._family, sibling, residual + self.jitter)
        return splu(self.V)

    @cached_property_depends_on('_varcomps')
    def V_logdet(self) -> float:
        """Compute log determinant of V using LU.

        Returns:
            float: log determinant of V
        """
        factor = self.V_lu
        if isinstance(factor, _FamilyFactor):
            return factor.logdet
        return np.log(np.abs(factor.U.diagonal())).sum()

    @cached_property_depends_on('_varcomps')
    def Vinv_y(self) -> np.ndarray:
        """Compute matrix-vector product of inverse of V and y

        Returns:
            np.ndarray: Vinv_e
        """
        # self.logger.debug('Calculating Vinv_y')
        return self.V_lu.solve(self.y)

    def Vinv_mat(self, dense_mat: np.ndarray) -> np.ndarray:
        """Calculate matrix-matrix product of inverse of V and a dense matrix

        Args:
            dense_mat (np.ndarray): 2-d array in the dense format

        Raises:
            ValueError: dense_mat must be a 2-d array
            ValueError: dimensions of V and dense_mat must match

        Returns:
            np.ndarray: matrix product in the dense format
        """        
        if dense_mat.ndim != 2:
            raise ValueError('dense_mat must be a 2-d array.')
        n, c = dense_mat.shape
        if n != self.n:
            raise ValueError(f'Input dims do not match: {self.n} and {n}.')
        return self.V_lu.solve(np.asarray(dense_mat, dtype=float))

    @cached_property_depends_on('_varcomps')
    def Vinv_Z(self):
        return self.Vinv_mat(self.Z)

    @cached_property_depends_on('_varcomps')
    def Vinv_e(self) -> np.ndarray:
        """Compute matrix-vector product of inverse of V and one-vector

        Returns:
            np.ndarray: Vinv_e
        """
        # self.logger.debug('Calculating V_inv_e')
        e = np.ones(self.n)
        return self.V_lu.solve(e)

    @cached_property_depends_on('_varcomps')
    def Vinv_varcomp_mats(self) -> Tuple[csc_matrix, ...]:
        """Compute matrix multiplications of inverse of V and all variance component matrices.

        Returns:
            Tuple[csc_matrix, ...]: a tuple holding all Vinv_varcomp_mat
        """
        # self.logger.debug('Calculating V_inv_varcomp_mats...')
        return tuple(self.Vinv_mat(matrix.toarray()) for matrix in self._varcomp_mats)

    @cached_property_depends_on('_varcomps')
    def P_attrs(self) -> GradHessComponents:
        """Compute ingredients for gradient and hessian (Vinv_y, Vinv_varcomp_mats).

        Returns:
            GradHessComponents: a NameTuple holding the computation results
        """
        # self.logger.debug('Calculating P_mats...')
        if self.has_covar:
            P_comp = self.Vinv_Z @ solve(self.Z.T @ self.Vinv_Z, self.Vinv_Z.T)
        else:
            P_comp = np.outer(
                self.Vinv_e, self.Vinv_e) / (np.ones(self.n) @ self.Vinv_e)
        P_y = self.Vinv_y - P_comp @ self.y
        P_varcomp_mats = tuple(
            self.Vinv_varcomp_mats[i] - P_comp @ self._varcomp_mats[i] for i in range(self.n_varcomps))
        return GradHessComponents(P_y=P_y, P_varcomp_mats=P_varcomp_mats)

    @cached_property_depends_on('_varcomps')
    def grad(self) -> np.ndarray:
        """Compute gradient.

        Returns:
            np.ndarray: 1-d array of gradient
        """
        if self._family is not None:
            return self._family_reml[1]
        P_y = self.P_attrs.P_y
        grad = -0.5 * np.array(
            [
                P_mat.diagonal().sum() - self.y @ P_mat @ P_y for P_mat in self.P_attrs.P_varcomp_mats
            ]
        )
        return grad if self.has_covar else grad / self.n

    @cached_property_depends_on('_varcomps')
    def hessian(self) -> np.ndarray:
        """Compute hessian.

        Returns:
            np.ndarray: 2-d n_varcomps-by-n_varcomps array
        """
        if self._family is not None:
            design = self.Z if self.has_covar else np.ones((self.n, 1))
            py, projection = _precision_terms(self.V_lu, self.y, design)
            sums = np.bincount(self._family.codes, weights=py)
            effects = np.column_stack((sums[self._family.codes], py))
            precision_effects = self.V_lu.solve(effects)
            z, projected_z = projection
            precision_effects -= projected_z @ (z.T @ precision_effects)
            information = .5 * (effects.T @ precision_effects)
            return information if self.has_covar else information / self.n
        hessian = np.empty((self.n_varcomps, self.n_varcomps))
        P_y = self.P_attrs.P_y
        for i in range(self.n_varcomps):
            for j in range(self.n_varcomps):
                hessian[i, j] = self.y @ self.P_attrs.P_varcomp_mats[i] @ self.P_attrs.P_varcomp_mats[j] @ P_y
        return 0.5 * hessian if self.has_covar else 0.5 * hessian / self.n

    @cached_property_depends_on('_varcomps')
    def reml_loglik(self) -> float:
        """Compute REML log likelihood

        Returns:
            float: REML log likelihood
        """
        if self._family is not None:
            return self._family_reml[0]
        if self.has_covar:
            ZT_Vinv_Z = self.Z.T @ self.Vinv_Z
            xx_ = solve(ZT_Vinv_Z, self.Vinv_Z.T.dot(self.y)) # @ self.y
            # P_comp: np.ndarray = self.Vinv_Z @ xx_
            # P_y: np.ndarray = self.Vinv_y - P_comp @ self.y
            P_y = self.Vinv_y - self.Vinv_Z @ xx_
            _, logdet_ZT_Vinv_Z = slogdet(ZT_Vinv_Z)
            return -0.5 * (self.V_logdet + logdet_ZT_Vinv_Z + self.y @ P_y)
        else:
            e = np.ones(self.n)
            yT_Vinv_y = self.y @ self.Vinv_y
            eT_Vinv_e = e @ self.Vinv_e
            eT_Vinv_y = e @ self.Vinv_y
            logdet_eT_Vinv_e = np.log(np.ones(self.n) @ self.Vinv_e)
            return -0.5 * (self.V_logdet + logdet_eT_Vinv_e + yT_Vinv_y - eT_Vinv_y ** 2 / eT_Vinv_e) / self.n

    @cached_property_depends_on('_varcomps')
    def _family_reml(self):
        factor = self.V_lu
        partition, sizes, weight = self._family, self._family.sizes, factor.mean_precision
        mean_trace = weight.sum() + factor.within_count * factor.within_precision
        family_trace = sizes @ weight
        if not self.has_covar:
            totals, within = self._family_y_sum, self._family_y_within
            mean = (weight @ totals) / family_trace
            centered_sum = totals - sizes * mean
            between = (centered_sum ** 2 / sizes) @ weight
            likelihood = -.5 * (factor.logdet + np.log(family_trace)
                               + between + factor.within_precision * within) / self.n
            square_weight = np.square(weight)
            quadratic = centered_sum ** 2 * square_weight
            trace = np.array([family_trace - (sizes ** 2) @ square_weight / family_trace,
                              mean_trace - sizes @ square_weight / family_trace])
            quad = np.array([quadratic.sum(),
                             np.sum(quadratic / sizes) + factor.within_precision ** 2 * within])
            return likelihood, .5 * (quad - trace) / self.n

        precision_z = factor.solve(self.Z)
        information = self.Z.T @ precision_z
        logdet_information = slogdet(information)[1]
        inverse = np.linalg.inv(information)
        precision_y = factor.solve(self.y)
        py = precision_y - precision_z @ (inverse @ (self.Z.T @ precision_y))
        grouped_z = np.column_stack([
            np.bincount(partition.codes, weights=z) for z in precision_z.T
        ])
        grouped_y = np.bincount(partition.codes, weights=py)
        trace = np.array([family_trace - np.trace(inverse @ (grouped_z.T @ grouped_z)),
                          mean_trace - np.trace(inverse @ (precision_z.T @ precision_z))])
        quad = np.array([grouped_y @ grouped_y, py @ py])
        likelihood = -.5 * (factor.logdet + logdet_information + self.y @ py)
        return likelihood, .5 * (quad - trace)

    def dense_reml_loglik(self, method: Literal['chol', 'cg'] = 'cg') -> float:
        """Evaluate the REML log likelihood using dense Cholesky or CG solves."""
        design = self.Z if self.has_covar else np.ones((self.n, 1))
        V = self.V.toarray()
        rhs = np.column_stack((self.y, design))
        if method == 'chol':
            factor = cho_factor(V)
            logdet_V = 2 * np.log(np.diag(factor[0])).sum()
            precision = cho_solve(factor, rhs)
        elif method == 'cg':
            _, logdet_V = slogdet(V)
            precision = np.empty_like(rhs)
            for i, column in enumerate(rhs.T):
                precision[:, i], status = cg(self.V, column, atol=0, rtol=1e-10)
                if status != 0:
                    raise RuntimeError(f'Conjugate gradient failed with status {status}.')
        else:
            raise ValueError(f'Unknown REML solver {method!r}')
        precision_y, precision_design = precision[:, 0], precision[:, 1:]
        information = design.T @ precision_design
        _, logdet_information = slogdet(information)
        projection = precision_y - precision_design @ solve(information, design.T @ precision_y)
        likelihood = -0.5 * (logdet_V + logdet_information + self.y @ projection)
        return likelihood if self.has_covar else likelihood / self.n

    def grid_search_reml(self) -> None:
        """Perform grid search on variance component parameters.
        """
        if self.optimized:
            print('WARNING: variance components are already optimized.')
        upper = self.y_var - self.jitter
        if not np.isfinite(upper) or upper <= 0:
            raise ValueError('Variance component estimation requires a nonconstant, finite phenotype.')
        max_loglik = float('-inf')
        max_sigma = None
        possible_values = np.linspace(0, upper, 101)
        for sigma in product(possible_values, repeat=self.n_varcomps - 1):
            if sum(sigma) > upper:
                continue
            self._varcomps = sigma + (upper - sum(sigma),)
            try:
                likelihood = self.reml_loglik
            except RuntimeError:
                # A zero residual component may make the covariance singular.
                continue
            if likelihood > max_loglik:
                max_loglik = likelihood
                max_sigma = self._varcomps
        if max_sigma is None:
            raise ValueError('No estimable variance components on the REML grid.')
        self._varcomps = max_sigma
        self.optimized = True

    def ai_reml(self) -> None:
        """Perform AI-REML algorithm to obtain maximum likelihood estimates of variance components.
        """
        if self.optimized:
            # self.logger.warning('Variance components are already optimized.')
            print('WARNING: variance components are already optimized.')
        previous = -np.inf
        best_likelihood = -np.inf
        best_components = self._varcomps
        for _ in range(50):
            likelihood = self.reml_loglik
            if likelihood > best_likelihood:
                best_likelihood, best_components = likelihood, self._varcomps
            if abs(likelihood - previous) < 1e-4:
                break
            previous = likelihood
            step = solve(self.hessian, self.grad)
            self._varcomps = tuple(np.maximum(np.asarray(self._varcomps) + step, self.y_var * 1e-5))
        self._varcomps = best_components
        self.optimized = True

    def scipy_optimize(self) -> None:
        """Perform LBFGS-B to optimize variance components.
        """
        if self.optimized:
            print('WARNING: variance components are already optimized.')
        yvar = self.y_var
        if not np.isfinite(yvar) or yvar <= 0:
            raise ValueError('Variance component estimation requires a nonconstant, finite phenotype.')

        def nll(components):
            self._varcomps = tuple(components)
            if self._family is not None:
                value, gradient = self._family_reml
                return -value, -gradient
            return -self.reml_loglik

        minimum = yvar * 1e-5
        bounds = [(minimum, yvar)] * self.n_varcomps
        initial = self._varcomps
        result = minimize(nll, x0=initial, method='L-BFGS-B', bounds=bounds,
                          jac=self._family is not None)
        if not result.success and self.n_varcomps > 2:
            # Keep the component array aligned with the inputs when a boundary term is dropped.
            eligible = np.flatnonzero(result.x[:-1] <= minimum * 1.01)
            for index in eligible:
                retry_bounds = bounds.copy()
                retry_bounds[index] = (0., 0.)
                candidate = minimize(nll, x0=result.x, method='L-BFGS-B', bounds=retry_bounds)
                if candidate.success:
                    print(f'WARNING: variance component {index} was fixed at zero.')
                    result = candidate
                    break
        if not result.success:
            self._varcomps = initial
            raise ValueError(f'Scipy minimize failed: {result.message}')
        self._varcomps = tuple(result.x)
        self.optimized = True

    def fit_snps_eff(self, gts: np.ndarray, standard_gwas: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Perform repeated OLS to estimate SNP effects and sampling variance-covariance.

        Args:
            gts (np.ndarray): 3-d array of genetic data.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: 3 arrays of SNP effects, covarinaces and standard errors.
        """
        if gts.ndim != 3 or gts.shape[0] != self.n:
            raise ValueError('Size of genotype matrix does not match phenotype.')
        if standard_gwas:
            gts = gts[:, :1, :]
        if self._family is None:
            precision_y, projection = self._snp_precision_terms
            alpha, covariance, _ = _fit_gls_batch(gts, self.V_lu, precision_y, projection)
        else:
            alpha, covariance = _fit_family_gls(gts, self.V_lu, *self._family_snp_terms)
        return alpha, covariance, np.sqrt(np.diagonal(covariance, axis1=1, axis2=2))

    @cached_property_depends_on('_varcomps')
    def _family_snp_terms(self):
        factor = self.V_lu
        precision_y = factor.solve(self.y)
        if not self.has_covar:
            return precision_y, np.empty((self.n, 0)), None
        precision_z = factor.solve(self.Z)
        inverse = np.linalg.pinv(self.Z.T @ precision_z, hermitian=True)
        precision_y -= precision_z @ (inverse @ (self.Z.T @ precision_y))
        return precision_y, precision_z, inverse

    @cached_property_depends_on('_varcomps')
    def _snp_precision_terms(self):
        return _precision_terms(self.V_lu, self.y, self.Z if self.has_covar else None)

    def _factor_for_rows(self, rows):
        if self._family is not None:
            sibling, residual = self._varcomps
            partition = _FamilyPartition(self._family.codes[rows])
            return _FamilyFactor(partition, sibling, residual + self.jitter)
        return splu(self.V[rows, :][:, rows])

    def robust_est(self, G: gtarray, G_sib: gtarray):
        """Combine design-specific direct effects using their sampling covariance."""
        if G.shape[0] != self.n or not np.array_equal(G.sid, G_sib.sid):
            raise ValueError('Genotype matrices must match the phenotype and SNPs.')
        alpha = np.full(G.shape[-1], np.nan)
        variance = np.full_like(alpha, np.nan)
        num_observed = G.num_obs_par_al
        par_status = G.par_status
        sibling_rows = np.array([G_sib.id_dict.get(iid, -1) for iid in G.ids])
        both_imputed = np.all(par_status == 1, axis=1)
        observed_father = par_status[:, 0] == 0
        observed_mother = par_status[:, 1] == 0
        group_models = [None] * 4
        for s in range(G.shape[-1]):
            observed = num_observed[:, s]
            three = observed == 3
            groups = [np.flatnonzero(mask) for mask in (
                observed == 4, three & both_imputed,
                three & observed_father, three & observed_mother)]
            if any(not len(group) for group in groups):
                continue
            both, one, paternal, maternal = groups
            if np.any(sibling_rows[one] < 0):
                continue
            x_pat = G.gts.data[paternal, :, s].copy()
            x_mat = G.gts.data[maternal, :, s].copy()
            x_pat[:, 2] = x_pat[:, 0] - np.floor(x_pat[:, 2])
            x_mat[:, 1] = x_mat[:, 0] - np.floor(x_mat[:, 1])
            matrices = (G.gts.data[both, :, s], G_sib.gts.data[sibling_rows[one], :, s], x_pat, x_mat)
            estimates = np.empty(4)
            weights = np.zeros((self.n, 4))
            for i, (rows, design) in enumerate(zip(groups, matrices)):
                cached = group_models[i]
                if cached is None or not np.array_equal(rows, cached[0]):
                    covariates = self.Z[rows] if self.has_covar else np.ones((len(rows), 1))
                    lu = self._factor_for_rows(rows)
                    precision_y, projection = _precision_terms(lu, self.y[rows], covariates)
                    cached = group_models[i] = (rows, lu, precision_y, projection)
                _, lu, precision_y, projection = cached
                coef, _, first_weights = _fit_gls_batch(
                    design[..., None], lu, precision_y, projection, return_weights=True)
                estimates[i] = coef[0, 0]
                weights[rows, i] = first_weights[:, 0]
            if np.all(np.isfinite(estimates)):
                if self._family is None:
                    covariance = weights.T @ (self.V @ weights)
                else:
                    sibling, residual = self._varcomps
                    totals = _group_columns(weights, self._family.codes, len(self._family.sizes))
                    covariance = (residual + self.jitter) * (weights.T @ weights) + sibling * (totals.T @ totals)
                alpha[s], variance[s] = _combine_direct_effects(estimates, covariance)
        return alpha[:, None], variance[:, None, None], np.sqrt(variance[:, None])

    def sib_diff_est(self, gts: np.ndarray):
        """
        Perform the sib-difference estimator
        """
        if gts.ndim != 3 or gts.shape[1] != 2:
            raise ValueError('Second dimension of the genotype matrix for the sib-difference estimator should be of size 2.')
        alpha, covariance, ses = self.fit_snps_eff(gts)
        return alpha[:, :1], covariance[:, :1, :1], ses[:, :1]

    def trios_sibs_est(self, gts: np.ndarray,
                       complete_trios_inds: np.ndarray,
                       sibs_inds: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Perform default regression on trios and sib-pairs.

        Args:
            gts (np.ndarray): 3-d array of genetic data.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: 3 arrays of SNP effects, covarinaces and standard errors.
        """
        if gts.ndim != 3 or gts.shape[0] != self.n:
            raise ValueError('Size of genotype matrix does not match phenotype.')
        loci = gts.shape[2]
        trios = np.asarray(complete_trios_inds)
        siblings = np.asarray(sibs_inds)
        groups, cross_covariance = self._trio_group_factors(trios, siblings, gts.shape[1])
        alpha = np.full(loci, np.nan)
        variance = np.full(loci, np.nan)
        block_size = max(1, min(256, (32 * 1024 ** 2) // (8 * self.n * gts.shape[1])))
        for start in range(0, loci, block_size):
            block = slice(start, start + block_size)
            estimates, variances, weights = [], [], []
            for indices, predictors, lu, precision_y, projection in groups:
                coef, covariance, first_weights = _fit_gls_batch(
                    gts[indices, :predictors, block], lu, precision_y, projection, return_weights=True)
                estimates.append(coef[:, 0])
                variances.append(covariance[:, 0, 0])
                weights.append(first_weights)
            a, b = variances
            only_trios = np.isfinite(a) & ~np.isfinite(b)
            only_sibs = ~np.isfinite(a) & np.isfinite(b)
            alpha[block][only_trios], variance[block][only_trios] = estimates[0][only_trios], a[only_trios]
            alpha[block][only_sibs], variance[block][only_sibs] = estimates[1][only_sibs], b[only_sibs]
            both = np.flatnonzero(np.isfinite(a) & np.isfinite(b))
            a, b = a[both], b[both]
            if self._family is None:
                cross = np.einsum('ij,ij->j', weights[0][:, both], cross_covariance @ weights[1][:, both])
            else:
                first, second, n_shared = cross_covariance
                sib_component = self._varcomps[0]
                cross = np.zeros(len(both))
                if n_shared and sib_component > 0:
                    t_sum = _group_columns(weights[0], first, n_shared)
                    s_sum = _group_columns(weights[1], second, n_shared)
                    cross = sib_component * np.einsum('ij,ij->j', t_sum[:, both], s_sum[:, both])
            # Closed-form minimum-variance weights for a two-estimate meta-analysis.
            determinant = a * b - cross * cross
            denominator = a + b - 2 * cross
            valid = (determinant > 0) & (denominator > 0)
            combined = both[valid]
            w = (b[valid] - cross[valid]) / denominator[valid]
            alpha[start + combined] = w * estimates[0][combined] + (1 - w) * estimates[1][combined]
            variance[start + combined] = determinant[valid] / denominator[valid]
        return alpha[:, None], variance[:, None, None], np.sqrt(variance[:, None])

    def _trio_group_factors(self, trios, siblings, predictors):
        cached = self.__dict__.get('_last_trio_groups')
        key = self._varcomps, predictors
        if cached is not None:
            previous, previous_trios, previous_sibs, result = cached
            if previous == key and np.array_equal(previous_trios, trios) and np.array_equal(previous_sibs, siblings):
                return result
        groups = []
        for indices, columns in ((trios.copy(), predictors), (siblings.copy(), 2)):
            z = self.Z[indices] if self.has_covar else np.ones((len(indices), 1))
            factor = self._factor_for_rows(indices)
            py, projection = _precision_terms(factor, self.y[indices], z)
            groups.append((indices, columns, factor, py, projection))
        if self._family is None:
            cross = self.V[trios, :][:, siblings]
        else:
            codes = self._family.codes
            n_groups = len(self._family.sizes)
            first, second = codes[trios], codes[siblings]
            shared = np.flatnonzero((np.bincount(first, minlength=n_groups) > 0) &
                                   (np.bincount(second, minlength=n_groups) > 0))
            mapping = np.full(n_groups, -1, dtype=np.intp)
            mapping[shared] = np.arange(len(shared))
            cross = mapping[first], mapping[second], len(shared)
        result = groups, cross
        self._last_trio_groups = key, groups[0][0], groups[1][0], result
        return result
