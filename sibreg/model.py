from numpy.linalg.linalg import solve
from .sibreg import convert_str_array, get_indices_given_ped, make_id_dict, find_par_gts
import numpy as np
import numpy.ma as ma
from numpy.linalg import slogdet
from scipy.sparse import csc_matrix, csr_matrix, tril
from scipy.sparse.linalg import splu, SuperLU, spsolve, cg
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize, OptimizeResult
from bgen_reader import open_bgen
from pysnptools.snpreader import Bed
import h5py
import subprocess
import logging
import tempfile
import gzip
from functools import wraps
from time import time
from typing import List, NamedTuple, Optional, Dict, Tuple, Callable, Hashable, Union
from typing_extensions import Literal
from collections import defaultdict
from itertools import product, combinations_with_replacement
from functools import lru_cache
import pandas as pd
from operator import attrgetter


FORMAT = '%(asctime)-15s :: %(levelname)s :: %(filename)s :: %(funcName)s :: %(message)s'
# numeric_level = getattr(logging, loglevel.upper(), None)
# if not isinstance(numeric_level, int):
#     raise ValueError('Invalid log level: %s' % loglevel)
logging.basicConfig(format=FORMAT, level=logging.DEBUG if __debug__ else logging.INFO)
logger = logging.getLogger(__name__)


def cached_property_depends_on(*args):
    attrs = attrgetter(*args)
    def decorator(func):
        _cache = lru_cache(maxsize=None)(lambda self, _: func(self))
        def _with_tracked(self):
            return _cache(self, attrs(self))
        return property(_with_tracked, doc=func.__doc__)
    return decorator


def timethis(f: Callable):
    @wraps(f)
    def wrap(*args, **kwargs):
        ts = time()
        result = f(*args, **kwargs)
        te = time()
        logger.info(f'{f.__name__} took: {te - ts:2.4f} sec.')
        return result
    return wrap


def match_phenotype_(ids: List[str], y: np.ndarray, pheno_ids: List[str]):
    """Match a phenotype to a genotype array by individual IDs.
    """
    id_dict = make_id_dict(ids)
    in_G_dict = np.array([x in id_dict for x in pheno_ids])
    y = y[in_G_dict]
    pheno_ids = pheno_ids[in_G_dict]
    pheno_id_dict = make_id_dict(pheno_ids)
    y = y[[pheno_id_dict[x] for x in ids]]
    return y


def coord2linear(ind1: int, ind2: int) -> int:
    row_ind, col_ind = max(ind1, ind2), min(ind1, ind2)
    return int(row_ind * (row_ind + 1) / 2 + col_ind)


def linear2coord(ind: int) -> Tuple[int, int]:
    ind += 1
    ind1 = np.ceil((2 * ind + 0.25) ** 0.5 - 0.5)
    ind2 = ind - (ind1 - 1) * ind1 / 2
    return int(ind1 - 1), int(ind2 - 1)


n_tril: Callable[[int], int] = lambda n: int(n * (n + 1) / 2)


def get_ids_with_par(par_gts_f: str,
                     gts_f: str,
                     ids: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray]:
    """Find ids with observed/imputed parents and family labels
    """
    # Imputed parental file
    par_gts_f_ = h5py.File(par_gts_f + '.hdf5', 'r')
    # Genotype file
    ped = convert_str_array(np.array(par_gts_f_['pedigree']))
    ped = ped[1:ped.shape[0], :]
    # Remove control families
    controls = np.array([x[0] == '_' for x in ped[:, 0]])
    ped = ped[np.logical_not(controls), :]
    # Get families with imputed parental genotypes
    fams = convert_str_array(np.array(par_gts_f_['families']))

    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f_: Bed = Bed(gts_f, count_A1=True)
        gts_ids = gts_f_.iid[:, 1]
        ids, observed_indices, imp_indices = get_indices_given_ped(ped,
                                                                   fams,
                                                                   gts_ids,
                                                                   ids=ids)
        gts_ids = gts_f_.iid[observed_indices, 1]
    elif gts_f[(len(gts_f) - 5):len(gts_f)] == '.bgen':
        gts_f_ = open_bgen(gts_f)
        gts_ids = gts_f_.samples
        ids, observed_indices, imp_indices = get_indices_given_ped(ped,
                                                                   fams,
                                                                   gts_ids,
                                                                   ids=ids)
        gts_ids = gts_ids[observed_indices]
    else:
        raise ValueError('Unknown filetype for observed genotypes file: ' +
                         str(gts_f))

    fams = fams[imp_indices]
    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, fams,
                                                      gts_id_dict)

    return ids, fam_labels


def run_gcta_grm(plink_path: str,
                 gcta_path: str,
                 filename: str,
                 output_path: str,
                 keep: Optional[List[str]] = None) -> None:
    """
    Build GRM using GCTM.

    Args:
        gcta_path : str
            path of gcta64 executable
        filename : str
            prefix of bed files; if '#' is in it, create a file containing file names of 22 chromosomes
        output_path : str
            prefix of output path
        keep : Optional[List]
            if List, create a txt file with each row being containing one IID to keep
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

def build_grm_arr(grm_path: str, id_dict: Dict[Hashable, int],
                  thres: float) -> np.ndarray:
    """Build GRM data array for HE regression.
    """
    gcta_id = dict()

    row_n = 1  # index starts from 1 in gcta grm output
    with open(f'{grm_path}.grm.id', 'r') as f:
        for line in f:
            id = line.strip().split('\t')[1]
            gcta_id[row_n] = id
            row_n += 1
    n = len(id_dict)
    grm_arr = np.full(int(n * (n + 1) / 2), np.nan)

    with gzip.open(f'{grm_path}.grm.gz', 'rt') as f:
        for line in f:
            id1, id2, _, gr = line.strip().split('\t')
            id1, id2 = gcta_id[int(id1)], gcta_id[int(id2)]
            ind1, ind2 = id_dict[id1], id_dict[id2]
            arr_ind = coord2linear(ind1, ind2)
            gr_: float = float(gr)
            grm_arr[arr_ind] = gr_ if gr_ >= thres else 0.

    if np.isnan(grm_arr).any():
        raise ValueError('Fewer pairs in GCTA GRM outout than expected.')
    return grm_arr


def build_ibdseg_arr(ibdseg_path: str, id_dict: Dict[Hashable, int],
                     keep: List[Union[str, int]], thres: float = 0.05) -> np.ndarray:
    """Build ibdseg array (lower triangular entries) from KING ibdseg output.

    Args:
        ibdseg_path (str): Path to ibdseg output.
        id_dict (Dict[Hashable, int]): dictionary of id-index pairs.
        keep (List[Union[str, int]]): list of ids to keep.
        thres (float, optional): sparsity threshold. Defaults to 0.0205.

    Returns:
        np.ndarray: 1-d ibdseg array.
    """
    logger.info(f'Reading {ibdseg_path}.seg...')
    king = pd.read_csv(ibdseg_path + '.seg', sep='\t')[['ID1', 'ID2', 'PropIBD']]
    king['ID1'] = king['ID1'].astype(str)
    king['ID2'] = king['ID2'].astype(str)

    # filter out IDs that are not in keep
    logger.info('Filtering...')
    king = king.query('ID1 in @keep & ID2 in @keep')

    n = len(id_dict)
    ibd_arr = np.zeros(n_tril(n))
    logger.info('Building ibd arr...')
    for row in king.itertuples():
        id1 = row.ID1
        id2 = row.ID2
        ind1, ind2 = id_dict[id1], id_dict[id2]
        arr_ind = coord2linear(ind1, ind2)
        ibd_arr[arr_ind] = row.PropIBD if row.PropIBD > thres else 0.
    for i in range(n):
        diag_ind = coord2linear(i, i)
        ibd_arr[diag_ind] = 1.
    logger.info('Done building ibd arr.')
    return ibd_arr


def write_ibdseg_to_grm(ibdseg_path: str, output_path: str,
                        ids: List[str], thres: float = 0.0205) -> None:
    """Write ibdseg from KING to GCTA GRM format.
    """
    logger.info(f'Reading {ibdseg_path}...')
    king = pd.read_csv(ibdseg_path, sep='\t')[['ID1', 'ID2', 'PropIBD']]
    king['ID1'] = king['ID1'].astype(str)
    king['ID2'] = king['ID2'].astype(str)
    logger.info('Filtering...')
    king = king.query('ID1 in @ids & ID2 in @ids')
    id_dict = make_id_dict(ids)
    logger.info('making full indices...')
    x = np.tril_indices(len(ids))
    full_king = pd.DataFrame({
        'row_ind': x[0] + 1,
        'col_ind': x[1] + 1,
        'non_missing': np.ones(len(x[0]), dtype=np.int)
    })

    def make_ind(df):
        rows, cols = [], []
        for row in df.itertuples():
            id1 = row.ID1
            id2 = row.ID2
            ind1, ind2 = id_dict[id1], id_dict[id2]
            row_ind, col_ind = max(ind1, ind2), min(ind1, ind2)
            rows.append(row_ind + 1)  # GCTA indices start from 1
            cols.append(col_ind + 1)
        return pd.Series({'row_ind': rows, 'col_ind': cols})

    logger.info('Making ibd row and col indices...')
    king[['row_ind', 'col_ind']] = make_ind(king)
    king_ = king.merge(full_king, on=['row_ind', 'col_ind'], how='right')[[
        'row_ind', 'col_ind', 'non_missing', 'PropIBD'
    ]]
    del king, full_king
    king_ = king_.fillna(0.)
    king_.loc[king_.row_ind == king_.col_ind, 'PropIBD'] = 1.
    # sparsification
    king_.loc[king_['PropIBD'] < thres, 'PropIBD'] = 0.
    logger.info(f'Writting ibd matrix to {output_path}.grm.gz...')
    king_.to_csv(output_path + '.grm.gz',
                 index=False,
                 header=None,
                 sep='\t',
                 compression='gzip')
    logger.info(f'Writing id info to {output_path}.grm.id...')
    pd.DataFrame({
        'FID': np.array(ids), 'ID': np.array(ids)
    }).to_csv(output_path + '.grm.id', index=False, header=None, sep='\t')
    logger.info('Finished writing.')


def build_sib_arr(fam_labels: List[str]) -> np.ndarray:
    """Build sibship array for HE regression.

    Args:
        fam_labels (List[str]): List of family id strings corresponding to each individual.

    Returns:
        np.ndarray: lower triangular entries of sibship matrix.
    """
    n = len(fam_labels)
    sib_arr = np.zeros(n_tril(n))

    label_indices = defaultdict(list)
    for l, f in enumerate(fam_labels):
        label_indices[f].append(l)

    for f, indices_lst in label_indices.items():
        for pair in combinations_with_replacement(indices_lst, 2):
            arr_ind = coord2linear(*pair)
            sib_arr[arr_ind] = 1.

    return sib_arr


def build_res_arr(nonzero_ind: np.ndarray) -> np.ndarray:
    """Build residual array for HE regression.
    """
    def is_diag(x, y): return int(x == y)
    x = [is_diag(*linear2coord(i)) for i in nonzero_ind]
    return np.array(x, dtype=np.float)


def build_pheno_prod_arr(y: np.ndarray, nonzero_ind: np.ndarray) -> np.ndarray:
    """Build phenotype product array (lower triangular entries).

    Args:
        y (np.ndarray): 1-d array of phenotype.
        nonzero_ind (np.ndarray): 1-d array of nonzeo indices.

    Returns:
        np.ndarray: pheno product array.
    """
    y = y - y.mean()
    vec_linear2coord = np.vectorize(linear2coord)
    logger.info('Building index vector...')
    ind1_vec, ind2_vec = vec_linear2coord(nonzero_ind)
    logger.info('Done.')
    return y[ind1_vec] * y[ind2_vec]


def he_reg(grm_arr: np.ndarray, sib_arr: np.ndarray, pheno_prod_arr: np.ndarray) -> Tuple[float, float, float]:
    """Perform HE regression.
    """
    X = np.vstack([np.ones(len(grm_arr)), grm_arr, sib_arr]).T
    results = np.linalg.lstsq(X, pheno_prod_arr, rcond=None)
    beta0, beta1, beta2 = results[0]
    err = results[1][0] / (pheno_prod_arr.shape[0] - 2)
    logger.info(beta0, beta1, beta2, err)
    # beta = np.linalg.solve(X.T.dot(X), X.T.dot(pheno_prod_arr))
    # res = pheno_prod_arr - X @ beta
    # err = res.T @ res / (pheno_prod_arr.shape[0] - 2)
    beta1 = 0 if beta1 < 0 else beta1
    beta2 = 0 if beta2 < 0 else beta2
    return beta1, beta2, err


class GradHessComponents(NamedTuple):
    P_y: np.ndarray
    P_varcomp_mats: Tuple[np.ndarray, ...]


class LinearMixedModel:
    """Wrapper of data and functions that compute estimates of variance components or SNP effects.
    """

    jitter = 1e-3
    _vec_linear2coord: Callable[..., Tuple[np.ndarray, np.ndarray]] = np.vectorize(
        linear2coord)
    _vec_coord2linear: Callable[..., np.ndarray] = np.vectorize(coord2linear)

    def __init__(self, y: np.ndarray, varcomp_arr_lst: Tuple[np.ndarray, ...],
                 covar_X: np.ndarray = None) -> None:
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
        self.y: np.ndarray
        if covar_X is not None:
            if covar_X.ndim != 2:
                raise ValueError('covar_X should be a 2-d array.')
            if covar_X.shape[0] != y.shape[0]:
                raise ValueError('Dimensions of y and covar_X do not match.')
            logger.info('Adjusting phenotype for fixed covariates...')
            self.y = self.fit_covar(y, covar_X)
            logger.info('Finished adjusting.')
        else:
            self.y = y
        self.n: int = len(y)
        def to_nz(arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]: 
            nz = arr.nonzero()[0]; return arr[nz], nz
        self._varcomp_mats: Tuple[csc_matrix, ...] = \
            tuple(self._build_sp_mat(*to_nz(arr), self.n) for arr in varcomp_arr_lst) \
            + (
                self._build_sp_mat(np.ones(self.n), self._vec_coord2linear(
                    np.arange(self.n), np.arange(self.n)), self.n),
        )
        self.n_varcomps: int = len(varcomp_arr_lst) + 1
        logger.debug(f'Number of variance components: {self.n_varcomps}.')
        self._varcomps: Tuple[float, ...] = tuple(
            self.y.var() if i == self.n_varcomps - 1 else 0. for i in range(self.n_varcomps))
        self.optimized: bool = False

    @staticmethod
    def fit_covar(y: np.ndarray, covar_X: np.ndarray) -> np.ndarray:
        if covar_X.ndim != 2:
            raise TypeError('Wrong input dimension.')
        if y.shape[0] != covar_X.shape[0]:
            raise TypeError(
                f'Input dimensions do not match. y: {y.shape}. covar_X: {covar_X.shape}.'
            )
        y_adjusted: np.ndarray = y - covar_X @ np.linalg.inv(
            covar_X.T @ covar_X) @ covar_X.T @ y
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
            raise ValueError('Input dims do not match.')
        out: np.ndarray = np.empty((m, n1, c))
        for i in range(m):
            out[i, :, :] = spsolve(sps_mat, dense_mat[i, :, :])
        return out

    def fit_snps_eff(
            self,
            gts: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Perform repeated OLS to estimate SNP effects and sampling variance-covariance in transformed model.

        Args:
            gts (np.ndarray): 3-d array of genetic data.

        Raises:
            RuntimeError: should adjust for covariates if not yet.
            ValueError: gts should be 3-d.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: 3 arrays of SNP effects, covarinaces and standard errors.
        """
        if not self.optimized:
            raise RuntimeError('Variance components are not optimized yer.')
        if gts.ndim != 3:
            raise ValueError('gts is expected to be of dimension 3.')
        gts = gts.transpose(2, 0, 1)
        # TODO: computing self._Vinv is inefficient; use LU decomp on self._grm_mat, self._sib_mat... instead
        Vinv_X: np.ndarray = self.sp_solve_dense3d(self.V, gts)
        XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', gts, Vinv_X)
        XT_Vinv_y: np.ndarray = np.einsum('...ij,i', gts, self.Vinv_y)
        alpha: np.ndarray = np.linalg.solve(XT_Vinv_X, XT_Vinv_y)
        alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
        alpha_ses: np.ndarray = np.sqrt(
            np.diagonal(alpha_cov, axis1=1, axis2=2))
        return alpha, alpha_cov, alpha_ses

    def _build_sp_mat(self, arr: np.ndarray, nonzero_ind: np.ndarray, n: int) -> csr_matrix:
        """Build sparse matrix using given lower triangular entries.

        Args:
            arr (np.ndarray): nonzero lower triangular entries.
            nonzero_ind (np.ndarray): linear indices of lower triangular entries.
            n (int): length of matrix dimension.

        Returns:
            csr_matrix: symmetric sparse matrix.
        """
        rows: np.ndarray
        cols: np.ndarray
        rows, cols = self._vec_linear2coord(nonzero_ind)
        tril_mat: csc_matrix = csc_matrix((arr, (rows, cols)), shape=(n, n))
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
        return self._varcomps
    
    @cached_property_depends_on('_varcomps')
    def V(self) -> csc_matrix:
        """Compute V.

        Returns:
            csc_matrix: V in sparse csc format
        """
        logger.debug('Calculating V...')
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
        logger.debug('Calculating V_lu')
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
        logger.debug('Calculating Vinv_y')
        return self.V_lu.solve(self.y)

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
        logger.debug('Calculating V_inv_varcomp_mats...')
        return tuple(spsolve(self.V, self._varcomp_mats[i]) for i in range(self.n_varcomps))

    @cached_property_depends_on('_varcomps')
    def P_attrs(self) -> GradHessComponents:
        """Compute ingredients for gradient and hessian (Vinv_y, Vinv_varcomp_mats).

        Returns:
            GradHessComponents: a NameTuple holding the computation results
        """
        logger.debug('Calculating P_mats...')
        P_comp: np.ndarray = np.outer(
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
        logger.debug('Calculating grad...')
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
        logger.debug('Calculating hessian...')
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
        e: np.ndarray = np.ones(self.n)
        yT_Vinv_y: float = self.y @ self.Vinv_y
        eT_Vinv_e: float = e @ self.Vinv_e
        eT_Vinv_y: float = e @ self.Vinv_y
        logdet_eT_Vinv_e: float = np.log(np.ones(self.n) @ self.Vinv_e)
        return -0.5 * (self.V_logdet + logdet_eT_Vinv_e + yT_Vinv_y - eT_Vinv_y ** 2 / eT_Vinv_e)

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
        Vinv_y: np.ndarray; Vinv_e: np.ndarray
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
                logger.warning('Conjugate gradient did not converge.')
            elif info1 < 1:
                raise RuntimeError(
                    'Illegal input or breakdown for conjugate gradient.')
            Vinv_e, info2 = cg(V, e, atol='legacy')
            if info2 == 0:
                pass
            elif info2 > 0:
                logger.warning('Conjugate gradient did not converge.')
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
            logger.warning('Variance components are already optimized.')
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
            logging.info(f'{i}: {self.reml_loglik}')
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
        logger.info('Starting AI-REML...')
        while True:
            logger.info(f'Iter: {iter}\tloglik: {self.reml_loglik}')
            if abs(self.reml_loglik - ll_) < 1e-4:
                break
            if iter > 50:
                self._varcomps = max_sigma
            ll_ = self.reml_loglik
            sigma: np.ndarray = np.array(
                self._varcomps) + np.linalg.solve(self.hessian, self.grad)
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
            logger.warning('Variance components are already optimized.')
        logger.info('Starting LBFGS-B...')

        def nll(x):
            self._varcomps = tuple(x)
            return -1 * self.reml_loglik
        res: OptimizeResult = minimize(nll, x0=np.array(self._varcomps),
                                       method='L-BFGS-B', bounds=[(0, self.y.var()) for i in range(self.n_varcomps)])
        if res.success:
            self.optimized = True
            self._varcomps = tuple(res.x)
            logger.info('Finished LBFGS-B.')
        else:
            raise ValueError('Scipy minimize failed.')
    
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
        ll: float = -0.5 * (logdet_V + logdet_e_T_V_inv_e + y_T_V_inv_y - e_T_V_inv_y ** 2 / e_T_V_inv_e)
        P_comp: np.ndarray = np.outer(
            Vinv_e, Vinv_e) / (np.ones(self.n) @ Vinv_e)
        P_y: np.ndarray = Vinv_y - P_comp @ self.y
        Vinv_varcomp_mats = tuple(solve(V, self._varcomp_mats[i].toarray()) for i in range(self.n_varcomps))
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
                hessian[i, j] = 0.5 * self.y @ P_varcomp_mats[i] @ P_varcomp_mats[j] @ P_y
        return ll, grad, hessian
