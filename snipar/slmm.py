import gzip
# import logging
from functools import lru_cache
from operator import attrgetter
from typing import Tuple, NamedTuple, Sequence
from typing_extensions import Literal
from itertools import product, combinations_with_replacement
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.optimize import minimize, OptimizeResult
from numpy.linalg import slogdet, solve, inv
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse import csc_matrix, tril
from scipy.sparse.linalg import splu, SuperLU, spsolve, cg
from snipar.types import FamLabels, SparseGRMRepr, Ids, IdDict
from snipar.gtarray import gtarray

# logger = logging.getLogger(__name__)



def build_sib_arr(fam_labels: FamLabels) -> SparseGRMRepr:    
    """Build lower-triangular nonzero entries of sibship matrix.

    Args:
        fam_labels (FamLabels): ndarray of family id strings corresponding to each individual.

    Returns:
        SparseGRMRepr: sparse GRM representation.
    """
    # logger.info('Building sibship arr...')
    data = []
    row_ind = []
    col_ind = []

    label_indices = defaultdict(list)
    for l, f in enumerate(fam_labels):
        label_indices[f].append(l)
    f = lambda lst: len(lst) * (len(lst) + 1) / 2
    n = sum(map(f, label_indices.values()))
    # logger.info('Done creating label_indices. ' + str(len(label_indices)) + ' ' + str(int(n)))

    for f, indices_lst in label_indices.items():
        for pair in combinations_with_replacement(indices_lst, 2):
            ind1, ind2 = max(pair[0], pair[1]), min(pair[0], pair[1])
            data.append(1)
            row_ind.append(ind1)
            col_ind.append(ind2)
    # logger.info(f'Done building sibship arr. nnz={len(data)}')
    return np.array(data), np.array(row_ind, dtype='uint32'), np.array(col_ind, dtype='uint32')


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
    # Determine delimiter and read header
    with open(ibdrel_path, 'r') as file:
        first_line = file.readline()
        delimiter = '\t' if '\t' in first_line else ' '  
        header = first_line.split(delimiter)
        header[-1] = header[-1].strip()  # Remove newline character
    if not 'ID1' in header and 'ID2' in header:
        raise ValueError('GRM should have ID1 and ID2 columns.')
    if 'PropIBD' in header:
        relcol = 'PropIBD'
    elif 'relatedness' in header:
        relcol = 'relatedness'
    else:
        raise ValueError('GRM should have PropIBD or relatedness column.')
    king = pd.read_csv(ibdrel_path,
                       sep=delimiter)[['ID1', 'ID2',relcol]]
    king['ID1'] = king['ID1'].astype(str)
    king['ID2'] = king['ID2'].astype(str)

    # filter out IDs that are not in keep
    # logger.info('Filtering...')
    king = king.query('ID1 in @keep & ID2 in @keep & PropIBD > @thres')
    king = king.reset_index(drop=True)
    n = king.shape[0]
    # the additional len(keep) entries are for diagonal terms
    data = np.zeros(n + len(keep))
    row_ind = np.zeros(n + len(keep), dtype=int)
    col_ind = np.zeros(n + len(keep), dtype=int)
    for row in king.itertuples():
        id1 = row.ID1
        id2 = row.ID2
        if row.PropIBD > 0.8:
            raise ValueError(f'Impossible ibd relatedness {row.PropIBD} for pair ({id1}, {id2}). MZ pairs should be removed.')
        ind1, ind2 = id_dict[id1], id_dict[id2]
        ind1, ind2 = max(ind1, ind2), min(ind1, ind2)
        # data[row.Index] = 0. if ignore_sib and row.InfType == 'FS' else row.PropIBD # if row.PropIBD > thres else 0.
        data[row.Index] = row.PropIBD # if row.PropIBD > thres else 0.
        if ind1 == ind2: raise RuntimeError('ind1 and ind2 cannot equal')
        # data[row.Index] = min(0.5, row.PropIBD) # if row.PropIBD > thres else 0.
        row_ind[row.Index] = ind1
        col_ind[row.Index] = ind2
    for i in range(len(keep)):
        data[n + i] = 1.
        row_ind[n + i] = i
        col_ind[n + i] = i
    # logger.info('Done building ibd arr.')
    return np.array(data), np.array(row_ind, dtype='uint32'), np.array(col_ind, dtype='uint32')


def build_grm_arr(grm_path: str, id_dict: IdDict,
                  thres: float) -> SparseGRMRepr:
    """Build GRM data array and corresponding indices.
    """
    gcta_id = dict()

    row_n = 1  # index starts from 1 in gcta grm output
    with open(f'{grm_path}.grm.id', 'r') as f:
        for line in f:
            id = line.strip().split('\t')[1]
            gcta_id[row_n] = id
            row_n += 1
    n = len(id_dict)
    # grm_arr = np.full(int(n * (n + 1) / 2), np.nan)
    data = []
    row_ind = []
    col_ind = []

    # logger.info('Building grm...')
    with gzip.open(f'{grm_path}.grm.gz', 'rt') as f:
        for line in f:
            id1, id2, _, gr = line.strip().split('\t')
            gr_: float = float(gr)
            if gr_ < thres:
                # ignore relatedness coeffs less than thres
                continue
            id1, id2 = gcta_id[int(id1)], gcta_id[int(id2)]
            if id1 not in id_dict or id2 not in id_dict:
                # ignore individuals that are not in id_dict
                continue
            ind1, ind2 = id_dict[id1], id_dict[id2]
            ind1, ind2 = max(ind1, ind2), min(ind1, ind2)
            data.append(gr_)
            row_ind.append(ind1)
            col_ind.append(ind2)
    # logger.info(f'Done building grm. nnz={len(data)}')
    return np.array(data), np.array(row_ind, dtype='uint32'), np.array(col_ind, dtype='uint32')


def match_grm_ids(ids: Ids,
                  fam_labels: FamLabels,
                  grm_path: str,
                  grm_source: Literal['gcta', 'ibdrel'],
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """Match ids with GRM individual ids.
    """
    orig_ids = pd.DataFrame({'ids': ids, 'fam_labels': fam_labels})
    if grm_source == 'gcta':
        grm_ids = pd.read_csv(f'{grm_path}.grm.id', sep='\s+', names=[
                              '_', 'ids_'], dtype={'_': str, 'ids_': str})[['ids_']]
    elif grm_source == 'ibdrel':
        grm_ids = pd.read_csv(grm_path,
                              sep='\t')[['ID1', 'ID2']]
        id1 = grm_ids['ID1'].to_numpy(dtype=str)
        id2 = grm_ids['ID2'].to_numpy(dtype=str)
        grm_ids = pd.DataFrame({'ids_': np.union1d(id1, id2)})
    else:
        raise ValueError('Incorrect source.')
    orig_ids = orig_ids.merge(grm_ids, how='inner',
                              left_on='ids', right_on='ids_')
    return orig_ids['ids'].to_numpy(), orig_ids['fam_labels'].to_numpy()


class GradHessComponents(NamedTuple):
    P_y: np.ndarray
    P_varcomp_mats: Tuple[np.ndarray, ...]


# cache class attribute calculation dependent on another attribute
# https://stackoverflow.com/questions/48262273/python-bookkeeping-dependencies-in-cached-attributes-that-might-change#answer-48264126
def cached_property_depends_on(*args):
    attrs = attrgetter(*args)
    def decorator(func):
        _cache = lru_cache(maxsize=None)(lambda self, _: func(self))
        def _with_tracked(self):
            return _cache(self, attrs(self))
        return property(_with_tracked, doc=func.__doc__)
    return decorator


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
                 add_jitter: bool = False) -> None:
        """Initialize a LinearMixedModel instance.

        Args:
            y (np.ndarray): 1-d array phenotpye
            varcomp_arr_lst (Sequence[np.ndarray, ...]): a Sequence of arbitrary length, holding variance component matrices, excluding the residual variance matrix.
            varcomps (Sequence[float, ...]): a sequence of variance component parameters.
            covar_X (np.ndarray, optional): 2-d array holding fixed covariates. Defaults to None.
            add_intercept: (bool): a boolean indicating whether to add intercept or not.
            add_jitter: (bool): a boolean indicating whether to add jitter to the diagonal of the covariance matrix or not.

        Raises:
            ValueError: y should be a 1-d array.
            ValueError: covar_X should be a 2-d array.
            ValueError: Dimensions of y and covar_X do not match.
        """
        if y.ndim != 1:
            raise ValueError('y should be a 1-d array.')
        self.y = y
        self.n = len(y)
        self.has_covar  = False
        if covar_X is not None:
            # logger.info('has covar')
            if covar_X.ndim == 1:
                covar_X = covar_X.reshape((self.n,1))
            if covar_X.ndim > 2:
                raise ValueError('covar_X should be a 1-d or 2-d array.')
            if covar_X.shape[0] != y.shape[0]:
                raise ValueError(
                    f'Dimensions of y and covar_X do not match ({y.shape} and {covar_X.shape}).')
            self.Z  = np.hstack((np.ones((self.n, 1),dtype=covar_X.dtype), covar_X)) if add_intercept else covar_X
            self.has_covar = True
        else:
            pass
            # self.logger.info('no covar')
        # self.logger.info(f'#individuals: {self.n}.')
            
        self._varcomp_mats: Tuple[csc_matrix, ...] = \
            tuple(self._build_sym_mat(data, row_ind, col_ind) for data, row_ind, col_ind in varcomp_arr_lst) \
            + (
                csc_matrix((np.ones(self.n), (np.arange(0, self.n), np.arange(0, self.n))), shape=(self.n, self.n)),
            )
        self.n_varcomps  = len(varcomp_arr_lst) + 1
        self.y_var  = self.y.var()
        # self.logger.info(f'Phenotypic variance: {self.y_var}.')
        self.jitter = self.y_var / 1000.0 if add_jitter else 0.0
        if varcomps is not None:
            if len(varcomps) != self.n_varcomps:
                raise ValueError('varcomps and varcomps_arr_lst must have the same length.')
            # _varcomps should be tuple for cached_dependent_property to work
            self._varcomps = tuple(varcomps)
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
        out  = __class__._empty_solver_result(sps_mat, dense_mat)
        for i in range(dense_mat.shape[0]):
            out[i, :, :] = spsolve(sps_mat, dense_mat[i, :, :])
        return out

    @staticmethod
    def sp_solve_dense3d_lu(sps_mat: csc_matrix,
                            dense_mat: np.ndarray) -> np.ndarray:
        """Compute product of the inverse of given sparse matrix and given 3-d dense array uisng LU; used for repeated OLS.

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
            for j in range(dense_mat.shape[-1]):
                out[i, :, j] = lu.solve(dense_mat[i, :, j])
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
        varcomps = np.array(self._varcomps)
        varcomps[-1] += self.jitter
        V: csc_matrix = sum(
            varcomps[i] * self._varcomp_mats[i] for i in range(self.n_varcomps)
        )
        return V

    @cached_property_depends_on('_varcomps')
    def V_lu(self) -> SuperLU:
        """Compute sparse LU factorization of V.

        Returns:
            SuperLU: wrapper object holding LU factorization of V
        """
        # self.logger.debug('Calculating V_lu')
        lu = splu(self.V)
        return lu

    @cached_property_depends_on('_varcomps')
    def V_logdet(self) -> float:
        """Compute log determinant of V using LU.

        Returns:
            float: log determinant of V
        """
        diag_l = self.V_lu.L.diagonal()
        diag_u = self.V_lu.U.diagonal()
        V_logdet = np.log(diag_l).sum() + np.log(diag_u).sum()
        return V_logdet

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
        out = np.empty((self.n, c))
        for i in range(c):
            out[:, i] = self.V_lu.solve(dense_mat[:, i])
        return out

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
        return tuple(self.Vinv_mat(self._varcomp_mats[i]) for i in range(self.n_varcomps))

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
            (self.Vinv_varcomp_mats[i] - P_comp @ self._varcomp_mats[i]).A for i in range(self.n_varcomps))
        return GradHessComponents(P_y=P_y, P_varcomp_mats=P_varcomp_mats)

    @cached_property_depends_on('_varcomps')
    def grad(self) -> np.ndarray:
        """Compute gradient.

        Returns:
            np.ndarray: 1-d array of gradient
        """
        # self.logger.debug('Calculating grad...')
        P_y = self.P_attrs.P_y
        grad = -0.5 * np.array(
            [
                P_mat.diagonal().sum() - self.y @ P_mat @ P_y for P_mat in self.P_attrs.P_varcomp_mats
            ]
        )
        return grad

    @cached_property_depends_on('_varcomps')
    def hessian(self) -> np.ndarray:
        """Compute hessian.

        Returns:
            np.ndarray: 2-d n_varcomps-by-n_varcomps array
        """
        # self.logger.debug('Calculating hessian...')
        hessian = np.empty((self.n_varcomps, self.n_varcomps))
        P_y = self.P_attrs.P_y
        for i in range(self.n_varcomps):
            for j in range(self.n_varcomps):
                hessian[i, j] = self.y @ self.P_attrs.P_varcomp_mats[i] @ self.P_attrs.P_varcomp_mats[j] @ P_y
        return 0.5 * hessian

    @cached_property_depends_on('_varcomps')
    def reml_loglik(self) -> float:
        """Compute REML log likelihood

        Returns:
            float: REML log likelihood
        """
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

    # TODO: add support for covariates; current version doesn't allow controlling for covariates
    # TODO: check sparsity of V and apply dense version if not sparse enough
    def dense_reml_loglik(self, method: Literal['chol', 'cg'] = 'cg') -> float:
        """Dense version of reml_loglik

        Args:
            method (Literal['chol', 'cg'], optional): a string specifying method for computing V inverse. Defaults to 'cg'.

        Raises:
            RuntimeError: Conjugate gradient did not converge
            RuntimeError: Illigal input for conjugate gradient

        Returns:
            float: REML log likelihood
        """
        varcomps = np.array(self._varcomps)
        varcomps[-1] += self.jitter
        V: csc_matrix = sum(
            varcomps[i] * self._varcomp_mats[i] for i in range(self.n_varcomps)
        )
        V_ = V.toarray()
        if method == 'chol':
            del V
        _, logdet_V = slogdet(V_)
        e = np.ones(self.n)
        if method == 'chol':
            c, low = cho_factor(V_)
            Vinv_y = cho_solve((c, low), self.y)
            Vinv_e = cho_solve((c, low), e)
        elif method == 'cg':
            Vinv_y, info1 = cg(V, self.y, atol='legacy')
            if info1 == 0:
                pass
            elif info1 > 0:
                self.logger.warning('Conjugate gradient did not converge.')
            elif info1 < 1:
                raise RuntimeError(
                    'Illegal input or breakdown for conjugate gradient.')
            Vinv_e, info2 = cg(V, e, atol='legacy')
            if info2 == 0:
                pass
            elif info2 > 0:
                self.logger.warning('Conjugate gradient did not converge.')
            elif info2 < 1:
                raise RuntimeError(
                    'Illegal input or breakdown for conjugate gradient.')
        y_T_V_inv_y = self.y @ Vinv_y
        e_T_V_inv_e = e @ Vinv_e
        e_T_V_inv_y = e @ Vinv_y
        logdet_e_T_V_inv_e = np.log(e_T_V_inv_e)
        return -0.5 * (logdet_V + logdet_e_T_V_inv_e + y_T_V_inv_y - e_T_V_inv_y ** 2 / e_T_V_inv_e)

    def grid_search_reml(self) -> None:
        """Perform grid search on variance component parameters.
        """
        if self.optimized:
            self.logger.warning('Variance components are already optimized.')
            print('WARNING: variance components are already optimized.')
        # reserve jitter for sibma_sib
        upper = np.var(self.y) - self.jitter
        step = upper / 100.
        max_loglik = float('-inf')
        max_sigma = (0.,)
        # logging.info('Starting grid search...')
        i = 1
        possible_values = np.arange(0.0, upper + step, step)
        for sigma in product(possible_values, repeat=self.n_varcomps - 1):
            if sum(sigma) > upper:
                continue
            self._varcomps = sigma + (upper - sum(sigma),)
            # logging.info(f'{i}: {self.reml_loglik}')
            i += 1
            if self.reml_loglik > max_loglik:
                max_loglik = self.reml_loglik
                max_sigma = self._varcomps
        # logging.info('Finished grid search.')
        self._varcomps = max_sigma
        self.optimized = True

    def ai_reml(self) -> None:
        """Perform AI-REML algorithm to obtain maximum likelihood estimates of variance components.
        """
        if self.optimized:
            # self.logger.warning('Variance components are already optimized.')
            print('WARNING: variance components are already optimized.')
        ll_ = float('inf')
        max_ll = float('-inf')
        iter = 1
        # self.logger.info('Starting AI-REML...')
        while True:
            # self.logger.info(f'Iter: {iter}\tloglik: {self.reml_loglik}')
            if abs(self.reml_loglik - ll_) < 1e-4:
                break
            if iter > 50:
                self._varcomps = max_sigma
            ll_ = self.reml_loglik
            sigma = np.array(
                self._varcomps) + solve(self.hessian, self.grad)
            sigma[sigma <= 0.] = self.y.var() * 1e-5
            self._varcomps = tuple(sigma)
            iter += 1
            if ll_ > max_ll:
                max_ll = ll_
                max_sigma = self._varcomps
        self.optimized = True
        # logger.info('Finished AI-REML.')

    def scipy_optimize(self) -> None:
        """Perform LBFGS-B to optimize variance components.
        """
        if self.optimized:
            # self.logger.warning('Variance components are already optimized.')
            print('WARNING: variance components are already optimized.')
        # self.logger.info('Starting LBFGS-B...')
        yvar = self.y.var()
        def nll(x):
            self._varcomps = tuple(x)
            return -1 * self.reml_loglik
        # def c(x):
        #     self.logger.info(inspect.currentframe().f_back.f_locals)  
        def callbackF(Xi):
            print(Xi, nll(Xi))
        res: OptimizeResult = minimize(nll, x0=np.array(self._varcomps), # options={'gtol': 1e-5, 'eps': 1e-5},
                                    method='L-BFGS-B', bounds=[(0.00001 * yvar, yvar) for i in range(self.n_varcomps)])
                                    # callback=callbackF)
        __varcomp_mats = self._varcomp_mats
        if not res.success:
            if self.n_varcomps == 2:
                raise ValueError(f'Scipy minimize failed: {res.message}')
            else:
                # self.logger.warning(f'Scipy minimize failed: {res.message}. Trying a reduced model...')
                print(f'WARNING: Scipy minimize failed: {res.message}. Trying a reduced model...')
                self.n_varcomps -= 1
                if (0.00001 * self.y_var >= res.x[0] or res.x[0] >= self.y_var) and \
                    (0.00001 * self.y_var >= res.x[1] or res.x[1] >= self.y_var):
                    raise ValueError(f'Scipy minimize failed: {res.message}).')
                elif 0.00001 * self.y_var >= res.x[0] or res.x[0] >= self.y_var:
                    # self.logger.warning(f'Dropping variance component 0: {res.x[0]}, and re-estimate...')
                    print(f'WARNING: dropping variance component 0: {res.x[0]}, and re-estimate...')
                    self._varcomp_mats = __varcomp_mats[1:]
                    self._varcomps = tuple(
                        self.y_var / self.n_varcomps for _ in range(self.n_varcomps))
                    res: OptimizeResult = minimize(nll, x0=np.array(self._varcomps), # options={'gtol': 1e-5, 'eps': 1e-5},
                                                method='L-BFGS-B', bounds=[(0.00001 * yvar, yvar) for i in range(self.n_varcomps)],
                                                callback=callbackF)
                    if not res.success:
                        raise ValueError(f'Scipy minimize failed: {res.message}') 
                elif 0.00001 * self.y_var >= res.x[1] or res.x[1] >= self.y_var:
                    # self.logger.warning(f'Dropping variance component 1: {res.x[1]}, and re-estimate...')
                    print(f'WARNING: dropping variance component 1: {res.x[1]}, and re-estimate...')
                    self._varcomp_mats = __varcomp_mats[::2]
                    self._varcomps = tuple(
                        self.y_var / self.n_varcomps for _ in range(self.n_varcomps))
                    res: OptimizeResult = minimize(nll, x0=np.array(self._varcomps), # options={'gtol': 1e-5, 'eps': 1e-5},
                                                method='L-BFGS-B', bounds=[(0.00001 * yvar, yvar) for i in range(self.n_varcomps)],
                                                callback=callbackF)
                    if not res.success:
                        raise ValueError(f'Scipy minimize failed: {res.message}')
                else: raise RuntimeError()
        self.optimized = True
        self._varcomps = tuple(res.x)
        # self.logger.info(f'Finished LBFGS-B. # of Likelihood evaluations: {res.nfev}')
    
    def fit_snps_eff(self, gts: np.ndarray, standard_gwas: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Perform repeated OLS to estimate SNP effects and sampling variance-covariance.

        Args:
            gts (np.ndarray): 3-d array of genetic data.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: 3 arrays of SNP effects, covarinaces and standard errors.
        """
        if standard_gwas:
            gts = gts[:, 0:1, :]
        # np.testing.assert_array_almost_equal(gts[:, 1, :], gts[:, 2, :])
        n, k, l = gts.shape
        if n != self.n:
            raise ValueError(f'Size of genotype matrix does not match pheno: {n},{self.n}.')
        if self.has_covar:
            gts_ = gts.reshape((gts.shape[0], int(k * l)))
            M_X = gts_ - self.Z.dot(solve(self.Z.T @ self.Z, self.Z.T.dot(gts_)))
            X_ = M_X.reshape((gts_.shape[0], k, l)).transpose(2, 0, 1)
            y = self.y - self.Z @ solve(self.Z.T @ self.Z, self.Z.T.dot(self.y))
            Vinv_X = self.sp_solve_dense3d_lu(self.V, X_)
            XT_Vinv_X = np.einsum('...ij,...ik', X_, Vinv_X)
            XT_Vinv_y = np.einsum('...ij,i', Vinv_X, y)
            alpha = solve(XT_Vinv_X, XT_Vinv_y)
            alpha_cov = np.linalg.inv(XT_Vinv_X)
        else:
            gts = gts.transpose(2, 0, 1)
            Vinv_X = self.sp_solve_dense3d_lu(self.V, gts)
            XT_Vinv_X = np.einsum('...ij,...ik', gts, Vinv_X)
            XT_Vinv_y = np.einsum('...ij,i', gts, self.Vinv_y)
            alpha = solve(XT_Vinv_X, XT_Vinv_y)
            alpha_cov = np.linalg.inv(XT_Vinv_X)

        alpha_ses: np.ndarray = np.sqrt(
            np.diagonal(alpha_cov, axis1=1, axis2=2))
        return alpha, alpha_cov, alpha_ses
    
    def robust_est(self, G: gtarray, G_sib: gtarray):
        num_obs_par_al = G.num_obs_par_al
        par_status = G.par_status
        ped = G.ped
        X = G.gts.data
        X_sib = G_sib.gts.data
        n, k, l = G.shape
        assert n == self.n
        assert l == G_sib.shape[2]

        alpha = np.full((l,1), fill_value=np.nan)
        alpha_ses = np.full((l,1), fill_value=np.nan)
        alpha_cov = np.full((l,1,1), fill_value=np.nan)
        ct = 0
        for s in range(l):
            if np.isnan(num_obs_par_al[:, s]).all():
                ct += 1
                continue
            notnan = np.isfinite(num_obs_par_al[:, s])
            both = (num_obs_par_al[:, s] == 4) * notnan
            pat = (num_obs_par_al[:, s] == 3) * (par_status[:, 0] == 0) * notnan
            mat = (num_obs_par_al[:, s] == 3) * (par_status[:, 1] == 0)  * notnan
            one = (num_obs_par_al[:, s] == 3) * ((par_status[:, 0] == 1) & (par_status[:, 1] == 1)) * notnan
            if both.sum() == 0 or one.sum() == 0 or pat.sum() == 0 or mat.sum() == 0:
                ct += 1
                # TODO: handle the case where one or more groups are missing
                continue
            X_both = X[both, :, s]
            X_pat = X[pat, :, s]
            X_pat[:, 2] = X_pat[:, 0] - X_pat[:, 2] // 1 # see np.divmod or np.modf
            X_mat = X[mat, :, s]
            X_mat[:, 1] = X_mat[:, 0] - X_mat[:, 1] // 1
            # TODO: check G_sib.ids in G.ids
            idx_one = [G_sib.id_dict[i] for i in G.ids[one]]
            X_one = X_sib[idx_one, :, s]

            y_both = self.y[both]
            y_pat = self.y[pat]
            y_mat = self.y[mat]
            y_one = self.y[one]
            V_both = self.V[both, :][:, both]
            V_pat = self.V[pat, :][:, pat]
            V_mat = self.V[mat, :][:, mat]
            V_one = self.V[one, :][:, one]

            V_both_one = self.V[both, :][:, one]
            V_both_pat = self.V[both, :][:, pat]
            V_both_mat = self.V[both, :][:, mat]
            V_one_pat = self.V[one, :][:, pat]
            V_one_mat = self.V[one, :][:, mat]
            V_pat_mat = self.V[pat, :][:, mat]

            
            if self.has_covar:
                X_both -= self.Z[both] @ solve(self.Z[both].T @ self.Z[both], self.Z[both].T.dot(X_both))
                y_both -= self.Z[both] @ solve(self.Z[both].T @ self.Z[both], self.Z[both].T.dot(y_both))

                X_pat -= self.Z[pat] @ solve(self.Z[pat].T @ self.Z[pat], self.Z[pat].T.dot(X_pat))
                y_pat -= self.Z[pat] @ solve(self.Z[pat].T @ self.Z[pat], self.Z[pat].T.dot(y_pat))

                X_mat -= self.Z[mat] @ solve(self.Z[mat].T @ self.Z[mat], self.Z[mat].T.dot(X_mat))
                y_mat -= self.Z[mat] @ solve(self.Z[mat].T @ self.Z[mat], self.Z[mat].T.dot(y_mat))

                X_one -= self.Z[one] @ solve(self.Z[one].T @ self.Z[one], self.Z[one].T.dot(X_one))
                y_one -= self.Z[one] @ solve(self.Z[one].T @ self.Z[one], self.Z[one].T.dot(y_one))


            X_both -= np.mean(X_both, axis=0)
            X_pat -= np.mean(X_pat, axis=0)
            X_mat -= np.mean(X_mat, axis=0)
            X_one -= np.mean(X_one, axis=0)
            y_both -= y_both.mean()
            y_pat -= y_pat.mean()
            y_mat -= y_mat.mean()
            y_one -= y_one.mean()
            
            try:
                Vinv_X_both = spsolve(V_both, X_both)
                XT_Vinv_X_both = X_both.T @ Vinv_X_both
                XT_Vinv_y_both = Vinv_X_both.T @ y_both
                alpha_both = solve(XT_Vinv_X_both, XT_Vinv_y_both)[0]
                alpha_cov_both = np.linalg.inv(XT_Vinv_X_both)[0,0]

                Vinv_X_one = spsolve(V_one, X_one)
                XT_Vinv_X_one = X_one.T @ Vinv_X_one
                XT_Vinv_y_one = Vinv_X_one.T @ y_one
                alpha_one = solve(XT_Vinv_X_one, XT_Vinv_y_one)[0]
                alpha_cov_one = np.linalg.inv(XT_Vinv_X_one)[0,0]

                Vinv_X_mat = spsolve(V_mat, X_mat)
                XT_Vinv_X_mat = X_mat.T @ Vinv_X_mat
                XT_Vinv_y_mat = Vinv_X_mat.T @ y_mat
                alpha_mat = solve(XT_Vinv_X_mat, XT_Vinv_y_mat)[0]
                alpha_cov_mat = np.linalg.inv(XT_Vinv_X_mat)[0,0]

                Vinv_X_pat = spsolve(V_pat, X_pat)
                XT_Vinv_X_pat = X_pat.T @ Vinv_X_pat
                XT_Vinv_y_pat = Vinv_X_pat.T @ y_pat
                alpha_pat = solve(XT_Vinv_X_pat, XT_Vinv_y_pat)[0]
                alpha_cov_pat = np.linalg.inv(XT_Vinv_X_pat)[0,0]
            except np.linalg.LinAlgError:
                continue
            

            cov_both_one = alpha_cov_both * Vinv_X_both[:, 0].T @ V_both_one @ Vinv_X_one[:, 0] * alpha_cov_one
            cov_both_pat = alpha_cov_both * Vinv_X_both[:, 0].T @ V_both_pat @ Vinv_X_pat[:, 0] * alpha_cov_pat
            cov_both_mat = alpha_cov_both * Vinv_X_both[:, 0].T @ V_both_mat @ Vinv_X_mat[:, 0] * alpha_cov_mat

            cov_one_mat = alpha_cov_one * Vinv_X_one[:, 0].T @ V_one_mat @ Vinv_X_mat[:, 0] * alpha_cov_mat
            cov_one_pat = alpha_cov_one * Vinv_X_one[:, 0].T @ V_one_pat @ Vinv_X_pat[:, 0] * alpha_cov_pat
            cov_pat_mat = alpha_cov_pat * Vinv_X_pat[:, 0].T @ V_pat_mat @ Vinv_X_mat[:, 0] * alpha_cov_mat

            S = np.block(
                [[alpha_cov_both, cov_both_one, cov_both_pat, cov_both_mat],
                [cov_both_one.T, alpha_cov_one, cov_one_pat, cov_one_mat],
                [cov_both_pat.T, cov_one_pat.T, alpha_cov_pat, cov_pat_mat],
                [cov_both_mat.T, cov_one_mat.T, cov_pat_mat.T, alpha_cov_mat]]
            )

            A = np.ones(4)
            robust_var = np.power(A @ solve(S, A), -1)
            alpha[s,0] = robust_var * (A @ solve(S, np.array([alpha_both, alpha_one, alpha_pat, alpha_mat])))
            alpha_cov[s,0,0] = robust_var
            alpha_ses[s,0] = robust_var ** 0.5
        return alpha, alpha_cov, alpha_ses

    def sib_diff_est(self, gts: np.ndarray):
        """
        Perform the sib-difference estimator
        """
        n, k, l = gts.shape
        if k != 2:
            raise ValueError('Second dimension of the genotype matrix for the sib-difference estimator should be of size 2.')
        if n != self.n:
            raise ValueError(f'Size of genotype matrix does not match pheno: {n},{self.n}.')
        alpha = np.full((l,1), fill_value=np.nan)
        alpha_ses = np.full((l,1), fill_value=np.nan)
        alpha_cov = np.full((l,1,1), fill_value=np.nan)
        if self.has_covar:
            gts_ = gts.reshape((gts.shape[0], int(k * l)))
            M_X = gts_ - self.Z.dot(solve(self.Z.T @ self.Z, self.Z.T.dot(gts_)))
            X_ = M_X.reshape((gts_.shape[0], k, l)).transpose(2, 0, 1)
            y = self.y - self.Z @ solve(self.Z.T @ self.Z, self.Z.T.dot(self.y))
            Vinv_X = self.sp_solve_dense3d_lu(self.V, X_)
            XT_Vinv_X = np.einsum('...ij,...ik', X_, Vinv_X)
            XT_Vinv_y = np.einsum('...ij,i', Vinv_X, y)
            alpha[:, 0] = solve(XT_Vinv_X, XT_Vinv_y)[:, 0]
            alpha_cov[:, 0, 0] = np.linalg.inv(XT_Vinv_X)[:, 0, 0]

        else:
            gts_ = gts.transpose(2, 0, 1)
            Vinv_X = self.sp_solve_dense3d_lu(self.V, gts_)
            XT_Vinv_X = np.einsum('...ij,...ik', gts_, Vinv_X)
            XT_Vinv_y = np.einsum('...ij,i', gts_, self.Vinv_y)
            alpha[:, 0] = solve(XT_Vinv_X, XT_Vinv_y)[:,0]
            alpha_cov[:, 0, 0] = np.linalg.inv(XT_Vinv_X)[:,0,0]
        alpha_ses[:, 0] = np.sqrt(
            alpha_cov)[:, 0, 0]
        return alpha, alpha_cov, alpha_ses
    
    def trios_sibs_est(self, gts: np.ndarray,
                       complete_trios_inds: np.ndarray,
                       sibs_inds: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Perform default regression on trios and sib-pairs.

        Args:
            gts (np.ndarray): 3-d array of genetic data.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: 3 arrays of SNP effects, covarinaces and standard errors.
        """
        n, k, l = gts.shape
        assert n == self.n
        if self.has_covar:
            gts_ = gts.reshape((gts.shape[0], int(k * l)))
            M_X = gts_ - self.Z.dot(solve(self.Z.T @ self.Z, self.Z.T.dot(gts_)))
            X = M_X.reshape((gts_.shape[0], k, l)).transpose(2, 0, 1)
            y = self.y - self.Z @ solve(self.Z.T @ self.Z, self.Z.T.dot(self.y))
        else:
            X = gts.transpose(2, 0, 1)
            # y = self.y - self.y.mean()
            y = self.y

        alpha = np.full((l,1), fill_value=np.nan)
        alpha_ses = np.full((l,1), fill_value=np.nan)
        alpha_cov = np.full((l,1,1), fill_value=np.nan)
        ct = 0

        X_trios = X[:, complete_trios_inds, :]
        X_sibs = X[:, sibs_inds, :2]
        y_trios = y[complete_trios_inds]
        y_sibs = y[sibs_inds]
        V_trios = self.V[complete_trios_inds, :][:, complete_trios_inds]
        V_sibs = self.V[sibs_inds, :][:, sibs_inds]
        V_trios_sibs = self.V[complete_trios_inds, :][:, sibs_inds]

        X_trios -= X_trios.mean(axis=1, keepdims=True)
        X_sibs -= X_sibs.mean(axis=1, keepdims=True)

        y_trios -= y_trios.mean()
        y_sibs -= y_sibs.mean()
        
        Vinv_X_trios = self.sp_solve_dense3d_lu(V_trios, X_trios)
        XT_Vinv_X_trios = np.einsum('...ij,...ik', X_trios, Vinv_X_trios)
        XT_Vinv_y_trios = np.einsum('...ij,i', Vinv_X_trios, y_trios)
        alpha_trios = solve(XT_Vinv_X_trios, XT_Vinv_y_trios)[:, 0]
        alpha_cov_trios = np.linalg.inv(XT_Vinv_X_trios)[:, 0, 0]

        Vinv_X_sibs = self.sp_solve_dense3d_lu(V_sibs, X_sibs)
        XT_Vinv_X_sibs = np.einsum('...ij,...ik', X_sibs, Vinv_X_sibs)
        XT_Vinv_y_sibs = np.einsum('...ij,i', Vinv_X_sibs, y_sibs)
        alpha_sibs = solve(XT_Vinv_X_sibs, XT_Vinv_y_sibs)[:, 0]
        alpha_cov_sibs = np.linalg.inv(XT_Vinv_X_sibs)[:, 0, 0]
        # numpy cannot do the following einsum
        # print(np.einsum('ji,ik->jk', Vinv_X_trios[:, :, 0], V_trios_sibs))
        for i in range(l):
            cov_trios_sibs = alpha_cov_trios[i] * Vinv_X_trios[i, :, 0].T @ V_trios_sibs @ Vinv_X_sibs[i, :, 0] * alpha_cov_sibs[i]

            S = np.block(
                [[alpha_cov_trios[i], cov_trios_sibs],
                [cov_trios_sibs.T, alpha_cov_sibs[i]]]
            )

            A = np.ones(2)
            robust_var = np.power(A @ solve(S, A), -1)
            alpha[i, 0] = robust_var * (A @ solve(S, np.array([alpha_trios[i], alpha_sibs[i]])))
            alpha_cov[i, 0, 0] = robust_var
            alpha_ses[i, 0] = robust_var ** 0.5
        return alpha, alpha_cov, alpha_ses