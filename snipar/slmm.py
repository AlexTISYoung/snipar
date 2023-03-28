import subprocess
import tempfile
import gzip
import numpy as np
import pandas as pd
from scipy.optimize import fmin_l_bfgs_b, minimize, OptimizeResult
from numpy.linalg import slogdet, solve, inv
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse import csc_matrix, tril
from scipy.sparse.linalg import splu, SuperLU, spsolve, cg
from typing import List, Dict, Tuple, Callable, NamedTuple, Optional
from typing_extensions import Literal
from collections import defaultdict
from itertools import combinations_with_replacement, product
import inspect
import logging
from functools import lru_cache
from operator import attrgetter
import snipar.utilities as utilities


logger = logging.getLogger(__name__)

FamLabels = np.ndarray
Ids = np.ndarray
IdDict = Dict[str, int]
SparseGRMRepr = Tuple[np.ndarray, np.ndarray, np.ndarray]



def build_sib_arr(fam_labels: FamLabels) -> SparseGRMRepr:    
    """Build lower-triangular nonzero entries of sibship matrix.

    Args:
        fam_labels (FamLabels): ndarray of family id strings corresponding to each individual.

    Returns:
        SparseGRMRepr: sparse GRM representation.
    """
    logger.info('Building sibship arr...')
    data = []
    row_ind = []
    col_ind = []

    label_indices = defaultdict(list)
    for l, f in enumerate(fam_labels):
        label_indices[f].append(l)
    f = lambda lst: len(lst) * (len(lst) + 1) / 2
    n = sum(map(f, label_indices.values()))
    logger.info('Done creating label_indices. ' + str(len(label_indices)) + ' ' + str(int(n)))

    for f, indices_lst in label_indices.items():
        for pair in combinations_with_replacement(indices_lst, 2):
            ind1, ind2 = max(pair[0], pair[1]), min(pair[0], pair[1])
            data.append(1)
            row_ind.append(ind1)
            col_ind.append(ind2)
    logger.info(f'Done building sibship arr. nnz={len(data)}')
    return np.array(data), np.array(row_ind, dtype='uint32'), np.array(col_ind, dtype='uint32')


def build_ibdrel_arr(ibdrel_path: str, id_dict: IdDict,
                     keep: Ids,
                     ignore_sib: bool = False,
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
    logger.info(f'Reading {ibdrel_path}.seg...')
    king = pd.read_csv(ibdrel_path + '.seg',
                       sep='\t')[['ID1', 'ID2', 'PropIBD', 'InfType']]
    king['ID1'] = king['ID1'].astype(str)
    king['ID2'] = king['ID2'].astype(str)

    # filter out IDs that are not in keep
    logger.info('Filtering...')
    king = king.query('ID1 in @keep & ID2 in @keep & PropIBD > @thres')
    king = king.reset_index(drop=True)
    n = king.shape[0]
    # the additional len(keep) entries are for diagonal terms
    data = np.zeros(n + len(keep))
    row_ind = np.zeros(n + len(keep), dtype=int)
    col_ind = np.zeros(n + len(keep), dtype=int)
    logger.info(f'Building ibd arr... zeroing out sib entries: {ignore_sib}.')
    for row in king.itertuples():
        id1 = row.ID1
        id2 = row.ID2
        ind1, ind2 = id_dict[id1], id_dict[id2]
        ind1, ind2 = max(ind1, ind2), min(ind1, ind2)
        data[row.Index] = 0. if ignore_sib and row.InfType == 'FS' else row.PropIBD # if row.PropIBD > thres else 0.
        row_ind[row.Index] = ind1
        col_ind[row.Index] = ind2
    for i in range(len(keep)):
        data[n + i] = 1.
        row_ind[n + i] = i
        col_ind[n + i] = i
    logger.info('Done building ibd arr.')
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

    logger.info('Building grm...')
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
    logger.info(f'Done building grm. nnz={len(data)}')
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
        grm_ids = pd.read_csv(grm_path + '.seg',
                              sep='\t')[['ID1', 'ID2']]
        id1 = grm_ids['ID1'].to_numpy(dtype=str)
        id2 = grm_ids['ID2'].to_numpy(dtype=str)
        grm_ids = pd.DataFrame({'ids_': np.union1d(id1, id2)})
    else:
        raise ValueError('Incorrect source.')
    orig_ids = orig_ids.merge(grm_ids, how='inner',
                              left_on='ids', right_on='ids_')
    return orig_ids['ids'].to_numpy(), orig_ids['fam_labels'].to_numpy()


def run_gcta_grm(plink_path: str,
                 gcta_path: str,
                 filename: str,
                 output_path: str,
                 keep: Optional[List[str]] = None) -> None:
    """
    Build GRM using GCTM.

    Args:
        gcta_path : str
            path of gcta64 executable.
        filename : str
            prefix of bed files; if '#' is in it, create a file containing file names of 22 chromosomes.
        output_path : str
            prefix of output path.
        keep : Optional[List]
            if List, create a txt file with each row being containing one IID to keep.
    """
    logger.info('Start running gcta grm...')
    args = [
        gcta_path,
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        if keep is not None:
            keep_file = f'{tmpdir}/keep.txt'
            with open(keep_file, 'w') as f:
                # f.write(f'#IID\n')
                for k in keep:
                    f.write(f'{str(k)}\t{str(k)}\n')

        if '#' in filename:
            args.append('--mbfile')
            chr_file = f'{tmpdir}/chr.txt'
            with open(chr_file, 'w') as f:
                for i in range(1, 23):
                    c = filename.replace('#', str(i))
                    plink_args = [
                        plink_path, '--bfile', c, '--rm-dup',
                        'exclude-mismatch', '--make-bed', '--out',
                        f'{tmpdir}/chr{str(i)}'
                    ]
                    if keep is not None:
                        plink_args += ['--keep', keep_file]
                    subprocess.run(plink_args)
                    f.write(f'{tmpdir}/chr{str(i)}' + '\n')
            args.append(chr_file)
        else:
            plink_args = [
                plink_path, '--bfile', filename, '--rm-dup',
                'exclude-mismatch', '--make-bed', '--out', f'{tmpdir}/chr'
            ]
            if keep is not None:
                plink_args += ['--keep', keep_file]
            subprocess.run(plink_args)
            args += ['--bfile', f'{tmpdir}/chr']

        args += ['--make-grm-gz', '--out', output_path]
        subprocess.run(args)
    logger.info('Finished running gcta grm.')


'''
class model(object):
    """Define a linear model with within-class correlations.

    Args:
        y : :class:`~numpy:numpy.array`
            1D array of phenotype observations
        X : :class:`~numpy:numpy.array`
            Design matrix for the fixed mean effects.
        labels : :class:`~numpy:numpy.array`
            1D array of sample labels

    Returns:
        model : :class:`snipar.model`

    """

    def __init__(self, y, X, labels, add_intercept=False):
        if y.shape[0] == X.shape[0] and X.shape[0] == labels.shape[0]:
            pass
        else:
            raise (ValueError('inconsistent sample sizes of response, covariates, and labels'))
        # Get sample size
        self.n = X.shape[0]
        if X.ndim == 1:
            X = X.reshape((self.n, 1))
        if add_intercept:
            X = np.hstack((np.ones((self.n, 1), dtype=X.dtype), X))
        self.X = X
        # Label mapping
        self.label_counts = dict()
        self.label_indices = dict()
        for l in range(0, labels.shape[0]):
            if labels[l] not in self.label_counts:
                self.label_counts[labels[l]] = 1
                self.label_indices[labels[l]] = [l]
            else:
                self.label_counts[labels[l]] += 1
                self.label_indices[labels[l]].append(l)
        self.y_lab = dict()
        self.X_lab = dict()
        for label in self.label_indices.keys():
            self.y_lab[label] = y[self.label_indices[label]]
            self.X_lab[label] = X[self.label_indices[label], :]
        self.n_labels = len(self.y_lab.keys())
        # response
        self.y = y
        self.labels = labels

    def alpha_mle(self, tau, sigma2, compute_cov=False, xtx_out=False):
        """
        Compute the MLE of alpha given variance parameters

        Args:
            sigma2 : :class:`float`
                variance of model residuals
            tau : :class:`float`
                ratio of variance of model residuals to variance explained by mean differences between classes

        Returns:
            alpha : :class:`~numpy:numpy.array`
                MLE of alpha

        """
        X_T_X = np.zeros((self.X.shape[1], self.X.shape[1]), dtype=np.float64)
        X_T_y = np.zeros((self.X.shape[1]), dtype=np.float64)

        for label in self.y_lab.keys():
            sigma_u = sigma2 / tau
            Sigma_lab = sigma_u * np.ones((self.label_counts[label], self.label_counts[label]))
            np.fill_diagonal(Sigma_lab, sigma_u + sigma2)
            Sigma_lab_inv = np.linalg.inv(Sigma_lab)
            X_T_X = X_T_X + np.dot(self.X_lab[label].T, Sigma_lab_inv.dot(self.X_lab[label]))
            X_T_y = X_T_y + np.dot(self.X_lab[label].T, Sigma_lab_inv.dot(self.y_lab[label]))

        if xtx_out:
            return [X_T_X, X_T_y.reshape((self.X.shape[1]))]
        else:
            alpha = np.linalg.solve(X_T_X, X_T_y)
            alpha = alpha.reshape((alpha.shape[0],))

            if compute_cov:
                alpha_cov = np.linalg.inv(X_T_X)
                return [alpha, alpha_cov]
            else:
                return alpha

    # Compute likelihood of data given beta, alpha
    def likelihood_and_gradient(self, sigma2, tau):
        """
        Compute the loss function, which is -2 times the likelihood along with its gradient

        Args:
            sigma2 : :class:`float`
                variance of model residuals
            tau : :class:`float`
                ratio of variance of model residuals to variance explained by mean differences between classes

        Returns:
            L, grad : :class:`float`
                loss function and gradient, divided by sample size

        """
        ## Likelihood
        alpha = self.alpha_mle(tau, sigma2)
        resid = self.y - self.X.dot(alpha)
        RSS = np.sum(np.square(resid))

        L = self.n * np.log(sigma2) + RSS / sigma2

        ## Gradient with respect to sigma2
        grad_sigma2 = self.n / sigma2 - RSS / np.square(sigma2)

        ## Gradient with respect to tau
        grad_tau = 0

        for label in self.y_lab.keys():
            resid_label = resid[self.label_indices[label]]
            resid_sum = np.sum(resid_label)
            resid_square_sum = np.square(resid_sum)
            # Add to likelihood
            L = L - resid_square_sum / (sigma2 * (tau + self.label_counts[label])) + np.log(
                1 + self.label_counts[label] / tau)
            # Add to grad sigma2
            grad_sigma2 += resid_square_sum / (np.square(sigma2) * (tau + self.label_counts[label]))
            # Add to grad tau
            grad_tau += (resid_square_sum / sigma2 - self.label_counts[label] * (
                        1 + self.label_counts[label] / tau)) / np.square(tau + self.label_counts[label])

        # Overall gradient vector
        grad = np.hstack((grad_sigma2, grad_tau))

        return L / self.n, grad / self.n

    def optimize_model(self, init_params):
        """
        Find the parameters that minimise the loss function for a given regularisation parameter

        Args:
            init_param : :class:`array`
                initial values for residual variance (sigma^2_epsilon) followed by ratio
                of residual variance to within-class variance (tau)

        Returns:
            optim : :class:`dict`
                dictionary with keys: 'success', whether optimisation was successful (bool);
                'warnflag', output of L-BFGS-B algorithm giving warnings; 'sigma2', MLE of
                residual variance; 'tau', MLE of ratio of residual variance to within-class variance;
                'likelihood', maximum of likelihood.

        """
        # Paramtere boundaries
        parbounds = [(0.00001, None), (0.00001, None)]
        # Optimize
        optimized = fmin_l_bfgs_b(func=lik_and_grad, x0=init_params,
                                  args=(self.y, self.X, self.labels),
                                  bounds=parbounds)

        # Get MLE
        optim = {}
        optim['success'] = True
        optim['warnflag'] = optimized[2]['warnflag']
        if optim['warnflag'] != 0:
            print('Optimization unsuccessful.')
            optim['success'] = False
        optim['sigma2'] = optimized[0][0]
        optim['tau'] = optimized[0][1]
        # Get parameter covariance
        optim['likelihood'] = -0.5 * np.float64(self.n) * (optimized[1] + np.log(2 * np.pi))

        return optim

    def sigma_inv_root(self, tau, sigma2):
        sigma_u = sigma2 / tau
        sigma2_nsqrt = dict()
        famsizes = np.unique(list(self.label_counts.values()))
        sigma2_nsqrt[1] = np.power(sigma_u + sigma2, -0.5)
        famsizes = famsizes[famsizes > 1]
        for famsize in famsizes:
            Sigma_lab = sigma_u * np.ones((famsize, famsize))
            np.fill_diagonal(Sigma_lab, sigma_u + sigma2)
            vals, vectors = np.linalg.eigh(Sigma_lab)
            vals = np.power(vals, 0.25)
            vectors = vectors / vals
            sigma2_nsqrt[famsize] = vectors.dot(vectors.T)
        return sigma2_nsqrt

    def predict(self, X):
        """
        Predict new observations based on model regression coefficients

        Args:
            X : :class:`array`
                matrix of covariates to predict from

        Returns:
            y : :class:`array`
                predicted values

        """
        if hasattr(self, 'alpha'):
            return X.dot(self.alpha)
        else:
            raise (AttributeError('Model does not have known regression coefficients. Try optimizing model first'))

    def set_alpha(self, alpha):
        self.alpha = alpha


def lik_and_grad(pars, *args):
    # Wrapper for function to pass to L-BFGS-B
    y, X, labels = args
    mod = model(y, X, labels)
    return mod.likelihood_and_gradient(pars[0], pars[1])

def simulate(n, alpha, sigma2, tau):
    """Simulate from a linear model with correlated observations within-class. The mean for each class
     is drawn from a normal distribution.

    Args:
        n : :class:`int`
            sample size
        alpha : :class:`~numpy:numpy.array`
            value of regression coefficeints
        sigma2 : :class:`float`
            variance of residuals
        tau : :class:`float`
            ratio of variance of residuals to variance of distribution of between individual means

    Returns:
        model : :class:`regrnd.model`
            linear model with repeated observations

    """
    c = alpha.shape[0]
    # X = np.random.randn((n * c)).reshape((n, c))
    X_cov = np.ones((c, c))
    np.fill_diagonal(X_cov, 1.2)
    X = np.random.multivariate_normal(np.zeros((c)), X_cov, n).reshape((n, c))
    labels = np.random.choice(n // 10, n)
    random_effects = np.sqrt(sigma2 // tau) * np.random.randn(n)
    y = X.dot(alpha) + random_effects[labels - 1] + np.random.randn(n) * np.sqrt(sigma2)
    return model(y, X, labels)

def fit_model(y, X, fam_labels, add_intercept=False, tau_init=1, return_model=True, return_vcomps=True, return_fixed=True):
    """Compute the MLE for the fixed effects in a family-based linear mixed model.

    Args:
        y : :class:`~numpy:numpy.array`
            vector of phenotype values
        X: :class:`~numpy:numpy.array`
            regression design matrix for fixed effects
        fam_labels : :class:`~numpy:numpy.array`
            vector of family labels: residual correlations in y are modelled between family members (that share a fam_label)
        add_intercept : :class:`bool`
            whether to add an intercept to the fixed effect design matrix

    Returns:
       model : :class:`snipar.model`
            the snipar model object, if return_model=True
       vcomps: :class:`float`
            the MLEs for the variance parameters: sigma2 (residual variance) and tau (ratio between sigma2 and family variance), if return_vcomps=True
       alpha : :class:`~numpy:numpy.array`
            MLE of fixed effects, if return_fixed=True
       alpha_cov : :class:`~numpy:numpy.array`
            sampling variance-covariance matrix for MLE of fixed effects, if return_fixed=True

    """
    # Optimize  model
    sigma_2_init = np.var(y)*tau_init/(1+tau_init)
    null_model = model(y, X, fam_labels, add_intercept=add_intercept)
    null_optim = null_model.optimize_model(np.array([sigma_2_init,tau_init]))
    # Create return list
    return_list = []
    if return_model:
        return_list.append(null_model)
    if return_vcomps:
        return_list += [null_optim['sigma2'], null_optim['tau']]
    if return_fixed:
        return_list += null_model.alpha_mle(null_optim['tau'], null_optim['sigma2'], compute_cov=True)
    return return_list
'''


class GradHessComponents(NamedTuple):
    P_y: np.ndarray
    P_varcomp_mats: Tuple[np.ndarray, ...]


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
    """Wrapper of data and functions that compute estimates of variance components or SNP effects.
    """
    logger = logger.getChild(__qualname__)

    _vec_linear2coord: Callable[..., Tuple[np.ndarray, np.ndarray]] = np.vectorize(
        utilities.linear2coord)
    _vec_coord2linear: Callable[..., np.ndarray] = np.vectorize(utilities.coord2linear)

    def __init__(self, y: np.ndarray, 
                 varcomp_arr_lst: Tuple[Tuple[np.ndarray, np.ndarray, np.ndarray], ...],
                 varcomps: Tuple[float, ...] = None,
                 covar_X: np.ndarray = None,
                 add_intercept: bool = False,
                 add_jitter: bool = False) -> None:
        """Initilize a LinearMixedModel instance.

        Args:
            y (np.ndarray): 1-d array phenotpye
            varcomp_arr_lst (Tuple[np.ndarray, ...]): a tuple of arbitrary length, holding variance component matrices, excluding the residual variance matrix
            covar_X (np.ndarray, optional): 2-d array holding fixed covariates. Defaults to None.

        Raises:
            ValueError: y should be a 1-d array.
            ValueError: covar_X should be a 2-d array.
            ValueError: Dimensions of y and covar_X do not match.
        """
        if y.ndim != 1:
            raise ValueError('y should be a 1-d array.')
        self.y: np.ndarray = y
        self.n: int = len(y)
        self.has_covar: bool = False
        if covar_X is not None:
            logger.info('has covar')
            if covar_X.ndim == 1:
                covar_X = covar_X.reshape((self.n,1))
            if covar_X.ndim > 2:
                raise ValueError('covar_X should be a 1-d or 2-d array.')
            if covar_X.shape[0] != y.shape[0]:
                raise ValueError(
                    f'Dimensions of y and covar_X do not match ({y.shape} and {covar_X.shape}).')
            self.Z: np.ndarray = covar_X
            self.has_covar = True
        else:
            logger.info('no covar')
        self.logger.info(f'#individuals: {self.n}.')

        if add_intercept and self.has_covar:
             self.Z = np.hstack((np.ones((self.n, 1),dtype=self.Z.dtype), self.Z))
            
        # def to_nz(arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        #     nz = arr.nonzero()[0]
        #     return arr[nz], nz
        # self._varcomp_mats: Tuple[csc_matrix, ...] = \
            # tuple(self._build_sp_mat_(*to_nz(arr), self.n) for arr in varcomp_arr_lst) \
            # + (
                # self._build_sp_mat_(np.ones(self.n), self._vec_coord2linear(
                    # np.arange(self.n), np.arange(self.n)), self.n),
        # )
        self._varcomp_mats: Tuple[csc_matrix, ...] = \
            tuple(self._build_sym_mat(data, row_ind, col_ind) for data, row_ind, col_ind in varcomp_arr_lst) \
            + (
                csc_matrix((np.ones(self.n), (np.arange(0, self.n), np.arange(0, self.n))), shape=(self.n, self.n)),
            )
        self.n_varcomps: int = len(varcomp_arr_lst) + 1
        self._varcomps: Tuple[float, ...]
        self.optimized: bool
        self.jitter: float
        self.y_var: float = self.y.var()
        self.logger.info(f'Phenotypic variance: {self.y_var}.')
        if add_jitter:
            self.jitter = self.y_var / 1000.
        else:
            self.jitter = 0.
        if varcomps is not None:
            if len(varcomps) != self.n_varcomps:
                raise ValueError('varcomps and varcomps_arr_lst must have the same length.')
            self._varcomps = varcomps 
            self.optimized = True
        else:
            self._varcomps = tuple(
                # self.y.var() if i == self.n_varcomps - 1 else 0. for i in range(self.n_varcomps))
                self.y_var / self.n_varcomps for i in range(self.n_varcomps))
            self.logger.info(f'init varcomps: {self._varcomps}.')
            # self._varcomps = tuple(
            #     y_var / self.n_varcomps for i in range(self.n_varcomps))
            self.optimized = False

    # deprecated
    @staticmethod
    def fit_covar(y: np.ndarray, covar_X: np.ndarray) -> np.ndarray:
        if covar_X.ndim != 2:
            raise TypeError('Wrong input dimension.')
        if y.shape[0] != covar_X.shape[0]:
            raise TypeError(
                f'Input dimensions do not match. y: {y.shape}. covar_X: {covar_X.shape}.'
            )
        y_adjusted: np.ndarray = y - covar_X.dot(solve(
            covar_X.T.dot(covar_X), covar_X.T.dot(y)))
        return y_adjusted

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
        if dense_mat.ndim != 3:
            raise ValueError('dense_mat must be a 3-d array.')
        n1: int
        n2: int
        m: int
        n: int
        c: int
        n1, n2 = sps_mat.shape
        m, n, c = dense_mat.shape
        if n != n2:
            raise ValueError(f'Input dims do not match: {n2} and {n}.')
        out: np.ndarray = np.empty((m, n1, c))
        for i in range(m):
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
        if dense_mat.ndim != 3:
            raise ValueError('dense_mat must be a 3-d array.')
        n1: int
        n2: int
        m: int
        n: int
        c: int
        n1, n2 = sps_mat.shape
        m, n, c = dense_mat.shape
        if n != n2:
            raise ValueError(f'Input dims do not match: {n2} and {n}.')
        lu = splu(sps_mat)
        out: np.ndarray = np.empty((m, n1, c))
        for i in range(m):
            for j in range(c):
                out[i, :, j] = lu.solve(dense_mat[i, :, j])
        return out

    

    # testing
    def fit_snps_eff(self, gts: np.ndarray, fam_labels: np.ndarray, ignore_na_fams: bool = True, ignore_na_rows: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Perform repeated OLS to estimate SNP effects and sampling variance-covariance.

        Args:
            gts (np.ndarray): 3-d array of genetic data.

        Raises:
            RuntimeError: should adjust for covariates if not yet.
            ValueError: gts should be 3-d.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: 3 arrays of SNP effects, covarinaces and standard errors.
        """
        n, k, l = gts.shape
        assert n == self.n
        if self.has_covar:
            logger.critical(f'starting... ignore_na_fams: {ignore_na_fams}. ignore_na_rows: {ignore_na_rows}.')
            if not ignore_na_fams and not ignore_na_rows:
                gts_ = gts.reshape((gts.shape[0], int(k * l)))
                M_X: np.ndarray = gts_ - self.Z.dot(solve(self.Z.T @ self.Z, self.Z.T.dot(gts_)))
                logger.info('Projecting genotype...')
                X_: np.ndarray = M_X.reshape((gts_.shape[0], k, l)).transpose(2, 0, 1)
                # if __debug__:
                #     M: np.ndarray = - self.Z @ solve(self.Z.T @ self.Z, self.Z.T)
                #     M.flat[::self.n + 1] += 1
                #     X__ = M.dot(gts_).reshape((self.n, k, l)).transpose(2, 0, 1)
                #     np.testing.assert_array_almost_equal(X__, X_)
                # self.logger.info('Done projecting genotype.')
                # self.logger.info('Projecting phenotype...')
                # y: np.ndarray = M @ self.y
                y: np.ndarray = self.y - self.Z @ solve(self.Z.T @ self.Z, self.Z.T.dot(self.y))
                self.logger.info('Start estimating snp effects...')
                Vinv_X: np.ndarray = self.sp_solve_dense3d_lu(self.V, X_)
                Vinv_y: np.ndarray = self.V_lu.solve(y)
                XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', X_, Vinv_X)
                XT_Vinv_y: np.ndarray = np.einsum('...ij,i', X_, Vinv_y)
                alpha: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
                alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
            else:
                alpha = np.zeros((gts.shape[2],gts.shape[1]),dtype=np.float_)
                alpha_cov = np.zeros((gts.shape[2],gts.shape[1],gts.shape[1]),dtype=np.float_) 
                for i in range(gts.shape[2]):
                    if ignore_na_fams:
                        is_na = np.sum(gts[:,:,i].mask,axis=1) > 0
                        is_na_fams = np.unique(fam_labels[is_na])
                        is_na = np.in1d(fam_labels, is_na_fams)
                        not_na = ~is_na
                    else: not_na = np.sum(gts[:,:,i].mask,axis=1) == 0
                    gts_ = gts[not_na, :, i]
                    Z = self.Z[not_na, :]
                    # gts_ = gts_.reshape((gts_.shape[0], int(k * l)))
                    # self.logger.info('Calculating projecting matrix...')
                    X_: np.ndarray = gts_ - Z.dot(solve(Z.T @ Z, Z.T.dot(gts_)))
                    # logger.info('Projecting genotype...')
                    # if __debug__:
                    #     M: np.ndarray = - self.Z @ solve(self.Z.T @ self.Z, self.Z.T)
                    #     M.flat[::self.n + 1] += 1
                    #     X__ = M.dot(gts_).reshape((self.n, k, l)).transpose(2, 0, 1)
                    #     np.testing.assert_array_almost_equal(X__, X_)
                    # self.logger.info('Done projecting genotype.')
                    # self.logger.info('Projecting phenotype...')
                    # y: np.ndarray = M @ self.y
                    y: np.ndarray = self.y[not_na] - Z @ solve(Z.T @ Z, Z.T.dot(self.y[not_na]))
                    # self.logger.info('Start estimating snp effects...')
                    v = self.V[not_na, :][:, not_na]
                    vinvx = spsolve(v, X_)
                    xtx = X_.T.dot(vinvx)
                    xty = vinvx.T.dot(y)
                    alpha[i,:] = np.linalg.solve(xtx,xty)
                    alpha_cov[i,:,:] = np.linalg.inv(xtx)
        else:
            gts = gts.transpose(2, 0, 1)
            # Vinv_X: np.ndarray = self.sp_solve_dense3d_lu(self.V, gts)
            # XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', gts, Vinv_X)
            # XT_Vinv_y: np.ndarray = np.einsum('...ij,i', gts, self.Vinv_y)

            logger.critical(f'starting... ignore_na_fams: {ignore_na_fams}. ignore_na_rows: {ignore_na_rows}.')
            if not ignore_na_fams and not ignore_na_rows:
                Vinv_X: np.ndarray = self.sp_solve_dense3d_lu(self.V, gts)
                XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', gts, Vinv_X)
                XT_Vinv_y: np.ndarray = np.einsum('...ij,i', gts, self.Vinv_y)
                alpha: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
                alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
            else:
                alpha = np.zeros((gts.shape[0],gts.shape[2]),dtype=np.float_)
                alpha_cov = np.zeros((gts.shape[0],gts.shape[2],gts.shape[2]),dtype=np.float_) 
                for i in range(gts.shape[0]):
        
                    ### debug ###
                    # is_na = np.where(np.sum(gts[i,:,:].mask,axis=1)>0)[0]
                    # print(is_na)
                    # print(len(is_na)) 
                    # for a in is_na:
                    #     f = np.array([fam_labels[a], ])
                    #     particular = np.in1d(fam_labels, f)
                    #     print(freqs[i])
                    #     print(fam_labels[a])
                    #     print(par_status[a, :])
                    #     print(gts[i, particular, :])
                    #     print()
                    # exit()
                    ############

                    if ignore_na_fams:
                        ## ignore whole fam
                        is_na = np.sum(gts[i,:,:].mask,axis=1) > 0
                        is_na_fams = np.unique(fam_labels[is_na])
                        is_na = np.in1d(fam_labels, is_na_fams)
                        not_na = ~is_na
                    else: not_na = np.sum(gts[:,:,i].mask,axis=1) == 0
                    
                    # not_na = np.sum(np.isnan(gts[:,:,i]),axis=1)==0
                    # unq, count = np.unique(gts[not_na,:,i], axis=0, return_counts=True)
            
                    v = self.V[not_na, :][:, not_na]
                    vinvx = spsolve(v, gts[i, not_na, :])
                    xtx = gts[i,not_na,:].T.dot(vinvx)
                    xty = vinvx.T.dot(self.y[not_na])
                    alpha[i,:] = np.linalg.solve(xtx,xty)
                    alpha_cov[i,:,:] = np.linalg.inv(xtx)

        # alpha: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
        # alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
        alpha_ses: np.ndarray = np.sqrt(
            np.diagonal(alpha_cov, axis1=1, axis2=2))
        return alpha, alpha_cov, alpha_ses

    # deprecated
    def _build_sp_mat_(self, arr: np.ndarray, nonzero_ind: np.ndarray, n: int) -> csc_matrix:
        """Build sparse matrix using given lower triangular entries.

        Args:
            arr (np.ndarray): nonzero lower triangular entries.
            nonzero_ind (np.ndarray): linear indices of lower triangular entries.
            n (int): length of matrix dimension.

        Returns:
            csc_matrix: symmetric sparse matrix.
        """
        rows: np.ndarray
        cols: np.ndarray
        rows, cols = self._vec_linear2coord(nonzero_ind)
        tril_mat: csc_matrix = csc_matrix((arr, (rows, cols)), shape=(n, n))
        triu_mat: csc_matrix = tril(tril_mat, k=-1, format='csc')
        return tril_mat + triu_mat.T
    
    def _build_sym_mat(self, data: np.ndarray, row_ind: np.ndarray, col_ind: np.ndarray) -> csc_matrix:
        """Build sparse matrix using given lower triangular entries.

        Args:
            data (np.ndarray): data array.
            row_ind (np.ndarray): row indices.
            col_ind (np.ndarray): column indices.

        Returns:
            csc_matrix: symmetric sparse matrix in csc format.
        """
        tril_mat: csc_matrix = csc_matrix((data, (row_ind, col_ind)), shape=(self.n, self.n))
        triu_mat: csc_matrix = tril(tril_mat, k=-1, format='csc')
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
        # return tuple(
        #     self._varcomps[i] if i < self.n_varcomps - 1 else self._varcomps[i] + self.jitter 
        #     for i in range(self.n_varcomps)
        # )
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
        lu: SuperLU = splu(self.V)
        return lu

    @cached_property_depends_on('_varcomps')
    def V_logdet(self) -> float:
        """Compute log determinant of V using LU.

        Returns:
            float: log determinant of V
        """
        diag_l: np.ndarray = self.V_lu.L.diagonal()
        diag_u: np.ndarray = self.V_lu.U.diagonal()
        V_logdet: float = np.log(diag_l).sum() + np.log(diag_u).sum()
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
        out: np.ndarray = np.empty((self.n, c))
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
        logger.debug('Calculating V_inv_e')
        e = np.ones(self.n)
        return self.V_lu.solve(e)

    @cached_property_depends_on('_varcomps')
    def Vinv_varcomp_mats(self) -> Tuple[csc_matrix, ...]:
        """Compute matrix multiplications of inverse of V and all variance component matrices.

        Returns:
            Tuple[csc_matrix, ...]: a tuple holding all Vinv_varcomp_mat
        """
        self.logger.debug('Calculating V_inv_varcomp_mats...')
        return tuple(self.Vinv_mat(self._varcomp_mats[i]) for i in range(self.n_varcomps))

    # TODO: test
    @cached_property_depends_on('_varcomps')
    def P_attrs(self) -> GradHessComponents:
        """Compute ingredients for gradient and hessian (Vinv_y, Vinv_varcomp_mats).

        Returns:
            GradHessComponents: a NameTuple holding the computation results
        """
        self.logger.debug('Calculating P_mats...')
        P_comp: np.ndarray
        if self.has_covar:
            P_comp = self.Vinv_Z @ solve(self.Z.T @ self.Vinv_Z, self.Vinv_Z.T)
        else:
            P_comp = np.outer(
                self.Vinv_e, self.Vinv_e) / (np.ones(self.n) @ self.Vinv_e)
        P_y: np.ndarray = self.Vinv_y - P_comp @ self.y
        P_varcomp_mats: Tuple[np.ndarray, ...] = tuple(
            (self.Vinv_varcomp_mats[i] - P_comp @ self._varcomp_mats[i]).A for i in range(self.n_varcomps))
        return GradHessComponents(P_y=P_y, P_varcomp_mats=P_varcomp_mats)

    @cached_property_depends_on('_varcomps')
    def grad(self) -> np.ndarray:
        """Compute gradient.

        Returns:
            np.ndarray: 1-d array of gradient
        """
        self.logger.debug('Calculating grad...')
        P_y: np.ndarray = self.P_attrs.P_y
        grad: np.ndarray = -0.5 * np.array(
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
        self.logger.debug('Calculating hessian...')
        hessian: np.ndarray = np.empty((self.n_varcomps, self.n_varcomps))
        P_y: np.ndarray = self.P_attrs.P_y
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
            ZT_Vinv_Z: np.ndarray = self.Z.T @ self.Vinv_Z
            xx_ = solve(ZT_Vinv_Z, self.Vinv_Z.T.dot(self.y)) # @ self.y
            # P_comp: np.ndarray = self.Vinv_Z @ xx_
            # P_y: np.ndarray = self.Vinv_y - P_comp @ self.y
            P_y: np.ndarray = self.Vinv_y - self.Vinv_Z @ xx_
            logdet_ZT_Vinv_Z: float
            _, logdet_ZT_Vinv_Z = slogdet(ZT_Vinv_Z)
            return -0.5 * (self.V_logdet + logdet_ZT_Vinv_Z + self.y @ P_y)
        else:
            e: np.ndarray = np.ones(self.n)
            yT_Vinv_y: float = self.y @ self.Vinv_y
            eT_Vinv_e: float = e @ self.Vinv_e
            eT_Vinv_y: float = e @ self.Vinv_y
            logdet_eT_Vinv_e: float = np.log(np.ones(self.n) @ self.Vinv_e)
            return -0.5 * (self.V_logdet + logdet_eT_Vinv_e + yT_Vinv_y - eT_Vinv_y ** 2 / eT_Vinv_e) / self.n

    # TODO: add Vinv_Z
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
        varcomps: np.ndarray = np.array(self._varcomps)
        varcomps[-1] += self.jitter
        V: csc_matrix = sum(
            varcomps[i] * self._varcomp_mats[i] for i in range(self.n_varcomps)
        )
        V_: np.ndarray = V.toarray()
        if method == 'chol':
            del V
        _: int
        logdet_V: float
        _, logdet_V = slogdet(V_)
        e: np.ndarray = np.ones(self.n)
        Vinv_y: np.ndarray
        Vinv_e: np.ndarray
        if method == 'chol':
            c, low = cho_factor(V_)
            Vinv_y = cho_solve((c, low), self.y)
            Vinv_e = cho_solve((c, low), e)
        elif method == 'cg':
            info1: int
            info2: int
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
        y_T_V_inv_y: float = self.y @ Vinv_y
        e_T_V_inv_e: float = e @ Vinv_e
        e_T_V_inv_y: float = e @ Vinv_y
        logdet_e_T_V_inv_e: float = np.log(e_T_V_inv_e)
        return -0.5 * (logdet_V + logdet_e_T_V_inv_e + y_T_V_inv_y - e_T_V_inv_y ** 2 / e_T_V_inv_e)

    def grid_search_reml(self) -> None:
        """Perform grid search on variance component parameters.
        """
        if self.optimized:
            self.logger.warning('Variance components are already optimized.')
        # reserve jitter for sibma_sib
        upper: float = np.var(self.y) - self.jitter
        step: float = upper / 100.
        max_loglik: float = float('-inf')
        max_sigma: Tuple[float, ...] = (0.,)
        logging.info('Starting grid search...')
        i: int = 1
        possible_values: np.ndarray = np.arange(0.0, upper + step, step)
        for sigma in product(possible_values, repeat=self.n_varcomps - 1):
            if sum(sigma) > upper:
                continue
            self._varcomps = sigma + (upper - sum(sigma),)
            # logging.info(f'{i}: {self.reml_loglik}')
            i += 1
            if self.reml_loglik > max_loglik:
                max_loglik = self.reml_loglik
                max_sigma = self._varcomps
        logging.info('Finished grid search.')
        self._varcomps = max_sigma
        self.optimized = True

    def ai_reml(self) -> None:
        """Perform AI-REML algorithm to obtain maximum likelihood estimates of variance components.
        """
        if self.optimized:
            logger.warning('Variance components are already optimized.')
        ll_: float = float('inf')
        max_ll: float = float('-inf')
        iter: int = 1
        max_sigma: Tuple[float, ...]
        self.logger.info('Starting AI-REML...')
        while True:
            self.logger.info(f'Iter: {iter}\tloglik: {self.reml_loglik}')
            if abs(self.reml_loglik - ll_) < 1e-4:
                break
            if iter > 50:
                self._varcomps = max_sigma
            ll_ = self.reml_loglik
            sigma: np.ndarray = np.array(
                self._varcomps) + solve(self.hessian, self.grad)
            sigma[sigma <= 0.] = self.y.var() * 1e-5
            self._varcomps = tuple(sigma)
            iter += 1
            if ll_ > max_ll:
                max_ll = ll_
                max_sigma = self._varcomps
        self.optimized = True
        logger.info('Finished AI-REML.')

    def scipy_optimize(self) -> None:
        """Perform LBFGS-B to optimize variance components.
        """
        if self.optimized:
            self.logger.warning('Variance components are already optimized.')
        self.logger.info('Starting LBFGS-B...')
        yvar = self.y.var()
        def nll(x):
            self._varcomps = tuple(x)
            return -1 * self.reml_loglik
        def c(x):
            self.logger.info(inspect.currentframe().f_back.f_locals)
        
        def callbackF(Xi):
            print(Xi, nll(Xi))
        res: OptimizeResult = minimize(nll, x0=np.array(self._varcomps), # options={'gtol': 1e-5, 'eps': 1e-5},
                                    method='L-BFGS-B', bounds=[(0.00001 * yvar, yvar) for i in range(self.n_varcomps)],)
                                    # callback=callbackF)
        __varcomp_mats = self._varcomp_mats
        if not res.success:
            if self.n_varcomps == 2:
                raise ValueError(f'Scipy minimize failed: {res.message}')
            else:
                self.logger.warning(f'Scipy minimize failed: {res.message}. Trying a reduced model...')
                self.n_varcomps -= 1
                if (0.00001 * self.y_var >= res.x[0] or res.x[0] >= self.y_var) and \
                    (0.00001 * self.y_var >= res.x[1] or res.x[1] >= self.y_var):
                    raise ValueError(f'Scipy minimize failed: {res.message}).')
                elif 0.00001 * self.y_var >= res.x[0] or res.x[0] >= self.y_var:
                    self.logger.warning(f'Dropping variance component 0: {res.x[0]}, and re-estimate...')
                    self._varcomp_mats = __varcomp_mats[1:]
                    self._varcomps = tuple(
                        self.y_var / self.n_varcomps for _ in range(self.n_varcomps))
                    res: OptimizeResult = minimize(nll, x0=np.array(self._varcomps), # options={'gtol': 1e-5, 'eps': 1e-5},
                                                method='L-BFGS-B', bounds=[(0.00001 * yvar, yvar) for i in range(self.n_varcomps)],
                                                callback=callbackF)
                    if not res.success:
                        raise ValueError(f'Scipy minimize failed: {res.message}') 
                elif 0.00001 * self.y_var >= res.x[1] or res.x[1] >= self.y_var:
                    self.logger.warning(f'Dropping variance component 1: {res.x[1]}, and re-estimate...')
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
        self.logger.info(f'Finished LBFGS-B. # of Likelihood evaluations: {res.nfev}')

    # TODO: add Z
    def test_grad_hessian(self) -> Tuple[float, np.ndarray, np.ndarray]:
        """Calculate REML loglikelihood, grad and hessian in dense format for testing purposes.

        Returns:
            Tuple[float, np.ndarray, np.ndarray]: a tuple containing the three results
        """
        varcomps: np.ndarray = np.array(self._varcomps)
        varcomps[-1] += self.jitter
        V_: np.ndarray = sum(
            varcomps[i] * self._varcomp_mats[i] for i in range(self.n_varcomps)
        )
        V: np.ndarray = V_.toarray()
        _: int
        logdet_V: float
        _, logdet_V = slogdet(V)
        e: np.ndarray = np.ones(self.n)
        Vinv_y: np.ndarray = solve(V, self.y)
        Vinv_e: np.ndarray = solve(V, e)
        y_T_V_inv_y: float = self.y @ Vinv_y
        e_T_V_inv_e: float = e @ Vinv_e
        e_T_V_inv_y: float = e @ Vinv_y
        logdet_e_T_V_inv_e: float = np.log(e_T_V_inv_e)
        ll: float = -0.5 * (logdet_V + logdet_e_T_V_inv_e +
                            y_T_V_inv_y - e_T_V_inv_y ** 2 / e_T_V_inv_e)
        P_comp: np.ndarray = np.outer(
            Vinv_e, Vinv_e) / (np.ones(self.n) @ Vinv_e)
        P_y: np.ndarray = Vinv_y - P_comp @ self.y
        Vinv_varcomp_mats = tuple(
            solve(V, self._varcomp_mats[i].toarray()) for i in range(self.n_varcomps))
        P_varcomp_mats: Tuple[np.ndarray, ...] = tuple(
            Vinv_varcomp_mats[i] - P_comp @ self._varcomp_mats[i].toarray() for i in range(self.n_varcomps))
        grad: np.ndarray = -0.5 * np.array(
            [
                P_mat.diagonal().sum() - self.y @ P_mat @ P_y for P_mat in P_varcomp_mats
            ]
        )
        hessian: np.ndarray = np.empty((self.n_varcomps, self.n_varcomps))
        for i in range(self.n_varcomps):
            for j in range(self.n_varcomps):
                hessian[i, j] = 0.5 * \
                    self.y @ P_varcomp_mats[i] @ P_varcomp_mats[j] @ P_y
        return ll, grad, hessian