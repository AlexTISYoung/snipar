from .sibreg import convert_str_array, get_indices_given_ped, make_id_dict, find_par_gts, model
import numpy as np
import numpy.ma as ma
from numpy.linalg import slogdet
from scipy.sparse import csc, csc_matrix, csr_matrix, tril
from scipy.sparse.linalg import inv, splu, SuperLU, spsolve, eigsh, cg
from scipy.linalg import cho_factor, cho_solve
from bgen_reader import open_bgen
from pysnptools.snpreader import Bed
import h5py
import subprocess
import logging
import tempfile
import gzip
from functools import wraps
from time import time
from typing import List, Optional, Dict, Tuple, Callable, Any, Iterable, Hashable, Union, Concatenate
from typing_extensions import Literal
from collections import defaultdict
from itertools import combinations
import pandas as pd


FORMAT = '%(asctime)-15s :: %(levelname)-8s :: %(filename)-8s :: %(funcName)-25s :: %(message)s'
logging.basicConfig(format=FORMAT, level=logging.INFO)

logger = logging.getLogger(__name__)


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


def combo_with_rep(it: Iterable[Any]) -> List[Tuple[Any, Any]]:
    """Return all combinations of length 2 plus duplicates (i.e., [a, a])
    """
    return list(combinations(it, 2)) + [(l, l) for l in it]


def coord2linear(ind1: int, ind2: int) -> int:
    row_ind, col_ind = max(ind1, ind2), min(ind1, ind2)
    return int(row_ind * (row_ind + 1) / 2 + col_ind)


def linear2coord(ind: int) -> Tuple[int, int]:
    ind += 1
    ind1 = np.ceil((2 * ind + 0.25) ** 0.5 - 0.5)
    ind2 = ind - (ind1 - 1) * ind1 / 2
    return int(ind1 - 1), int(ind2 - 1)


n_tril: Callable[[int], int] = lambda n: int(n * (n + 1) / 2)


class OrderedIdPair(tuple):
    """Container to store ordered id pair.
    NOTE: not used now
    """
    def __new__(cls, lst):
        if len(lst) != 2:
            raise ValueError('Input must of length 2.')
        if any(not isinstance(i, str) for i in lst):
            raise TypeError('Ids must be string.')
        id1, id2 = sorted(lst)
        return tuple.__new__(cls, [id1, id2])


def get_ids_with_par(par_gts_f: str,
                     gts_f: str,
                     ids: List[str] = None) -> Tuple[List[str], List[str]]:
    """Find ids with observed/imputed parents and family labels

    Args:
        ids: (List[str]) list of samples ids that we want to analyse; possibly from pheno_ids
    """
    # Imputed parental file
    par_gts_f = h5py.File(par_gts_f + '.hdf5', 'r')
    # Genotype file
    ped = convert_str_array(np.array(par_gts_f['pedigree']))
    ped = ped[1:ped.shape[0], :]
    # Remove control families
    controls = np.array([x[0] == '_' for x in ped[:, 0]])
    ped = ped[np.logical_not(controls), :]
    # Get families with imputed parental genotypes
    fams = convert_str_array(np.array(par_gts_f['families']))

    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f = Bed(gts_f, count_A1=True)
        gts_ids = gts_f.iid[:, 1]
        ids, observed_indices, imp_indices = get_indices_given_ped(ped,
                                                                   fams,
                                                                   gts_ids,
                                                                   ids=ids)
        gts_ids = gts_f.iid[observed_indices, 1]
    elif gts_f[(len(gts_f) - 5):len(gts_f)] == '.bgen':
        gts_f = open_bgen(gts_f)
        gts_ids = gts_f.samples
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


def make_id_pair_dict(ids: List[str]) -> Dict[OrderedIdPair, int]:
    """Make ordered ID pair dict from list of IDs.
    NOTE: not used now
    """
    logger.info('Making ID pairs...')
    # convert ids to str
    ids = map(str, ids)
    pairs = map(OrderedIdPair, combo_with_rep(ids))

    logger.info('Building ibd arr...')
    id_pair_dict = {p: ind for ind, p in enumerate(pairs)}
    # for p in map(OrderedIdPair, pairs):
    #     id_pair_dict[p] = i
    #     i += 1
    return id_pair_dict


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
            gr = float(gr)
            grm_arr[arr_ind] = gr if gr >= thres else 0.

    if np.isnan(grm_arr).any():
        raise ValueError('Fewer pairs in GCTA GRM outout than expected.')
    return grm_arr


def build_ibdseg_arr(ibdseg_path: str, id_dict: Dict[Hashable, int],
                     keep: List[Union[str, int]], thres: float = 0.0205) -> np.ndarray:
    """Build ibdseg array (lower triangular entries) from KING ibdseg output.

    Args:
        ibdseg_path (str): Path to ibdseg output.
        id_dict (Dict[Hashable, int]): dictionary of id-index pairs.
        keep (List[Union[str, int]]): list of ids to keep.
        thres (float, optional): sparsity threshold. Defaults to 0.0205.

    Returns:
        np.ndarray: 1-d ibdseg array.
    """                     
    """Build ibdseg coeff array from king ibdseg output.
    """
    logger.info(f'Reading {ibdseg_path}...')
    king = pd.read_csv(ibdseg_path, sep='\t')[['ID1', 'ID2', 'PropIBD']]
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
        ibd_arr[arr_ind] = row.PropIBD if row.PropIBD > thres else 0.   # 0.0205
    for i in range(n):
        diag_ind = coord2linear(i, i)
        ibd_arr[diag_ind] = 1.
    logger.info('Done building ibd arr')
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
        for pair in combo_with_rep(indices_lst):
            arr_ind = coord2linear(*pair)
            sib_arr[arr_ind] = 1.

    return sib_arr


def build_res_arr(nonzero_ind: np.ndarray) -> np.ndarray:
    """Build residual array for HE regression.
    """
    def is_diag(x, y): return int(x == y)
    n_ = len(nonzero_ind)
    x = [is_diag(*linear2coord(i)) for i in range(n_)]
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


class LinearMixedModel:
    """Wrapper of data and functions that compute estimates of variance components or SNP effects.
    """

    jitter = 1e-5

    def __init__(self, y: np.ndarray,
                 grm_arr: np.ndarray, sib_arr: np.ndarray,
                 include_covar: bool = True) -> None:
        """Initialize a linear model with within-class correlations.

        Args:
            y (np.ndarray): 1-d array of phenotpye observations.
            grm_arr (np.ndarray): 1-d array of GRM coefficients.
            sib_arr (np.ndarray): 1-d array of sibship coefficients.
            include_covar (bool, optional): whether to include fixed covariates. Defaults to True.

        Raises:
            ValueError: y should be 1-d.
        """
        if y.ndim != 1:
            raise ValueError('y should be a 1-d array.')

        self.y: np.ndarray = y
        self.include_covar: bool = include_covar
        if include_covar:
            self._covar_adjusted: bool = False
        self.n: bool = len(y)
        grm_nz: np.ndarray = grm_arr.nonzero()[0]
        grm_arr: np.ndarray = grm_arr[grm_nz]
        sib_nz: np.ndarray = sib_arr.nonzero()[0]
        sib_arr: np.ndarray = sib_arr[sib_nz]
        self._vec_linear2coord: Callable[..., Tuple[np.ndarray, np.ndarray]] = np.vectorize(
            linear2coord)
        self._vec_coord2linear: Callable[...,
                                         np.ndarray] = np.vectorize(coord2linear)
        self._grm_mat: csc_matrix = self._build_sp_mat(grm_arr, grm_nz, self.n)
        self._sib_mat: csc_matrix = self._build_sp_mat(sib_arr, sib_nz, self.n)
        self._res_mat: csc_matrix = self._build_sp_mat(np.ones(self.n), self._vec_coord2linear(
            np.arange(self.n), np.arange(self.n)), self.n)
        self.grm_nnz: int = grm_arr.shape[0]
        assert self._grm_mat.diagonal().sum() == self.n
        assert self._sib_mat.diagonal().sum() == self.n
        assert self._res_mat.diagonal().sum() == self.n

    def fit_covar(self, covar_X: np.ndarray) -> None:
        """Fit and adjust for covariates.

        Args:
            covar_X (np.ndarray): array of fixed covariates.

        Raises:
            TypeError: covar_X should be of dimension 3.
            TypeError: dimensioin 1 of covar_X and length of self.y should match.
        """
        if not self.include_covar:
            logger.warning('This model does not include covariates.')
            return
        if getattr(self, '_covar_adjusted'):
            logger.warning('Covariates have been adjusted for.')
            return

        if covar_X.ndim != 2:
            raise TypeError('Wrong input dimension.')
        if self.y.shape[0] != covar_X.shape[0]:
            raise TypeError(
                f'Input dimensions do not match. y: {self.y.shape}. covar_X: {covar_X.shape}.'
            )
        self.y: np.ndarray = self.y - covar_X @ np.linalg.inv(
            covar_X.T @ covar_X) @ covar_X.T @ self.y
        self.covar_adjusted: bool = True

    @staticmethod
    def sp_mul_dense3d(sps_mat: csc_matrix,
                       dense_mat: np.ndarray) -> np.ndarray:
        """Compute product of given sparse matrix and 3-d dense array; used for repeated OLS.

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
            out[i, :, :] = sps_mat.dot(dense_mat[i, :, :])
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
        if self.include_covar and not getattr(self, '_covar_adjusted'):
            raise RuntimeError('Please adjust phenotype for covariates first.')
        if gts.ndim != 3:
            raise ValueError('gts is expected to be of dimension 3.')
        if not hasattr(self, '_Vinv'):
            raise RuntimeError('V inverse is not computed.')
        gts = gts.transpose(2, 0, 1)
        # TODO: computing self._Vinv is inefficient; use LU decomp on self._grm_mat, self._sib_mat... instead
        Vinv_X: np.ndarray = self.sp_mul_dense3d(self._Vinv, gts)
        X_T_Vinv_X: np.ndarray = np.einsum('...ij,...ik', gts, Vinv_X)
        X_T_Vinv_y: np.ndarray = np.einsum('...ij,i', gts, self._Vinv_y)
        alpha: np.ndarray = np.linalg.solve(X_T_Vinv_X, X_T_Vinv_y)
        alpha_cov: np.ndarray = np.linalg.inv(X_T_Vinv_X)
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

    def _get_V(self, sigma_grm: float, sigma_sib: float, sigma_res: float) -> csc_matrix:
        """Get V matrix using given variance component parameters.

        Args:
            sigma_grm (float): variance component contributed by GRM.
            sigma_sib (float): variance component contributed by sibship.
            sigma_res (float): residual variance component (need to be >= self.jitter).

        Returns:
            csc_matrix: sparse V matrix.
        """
        if sigma_res < 1e-5:
            logging.warning('sigma_res should at >= self.jitter')
        return sigma_grm * self._grm_mat + sigma_sib * self._sib_mat + sigma_res * self._res_mat

    def reml_loglik(self, sigma_grm: float, sigma_sib: float, sigma_res: float) -> float:
        """Compute REML log likelihood

        Args:
            sigma_grm (float): variance component contributed by GRM
            sigma_sib (float): variance component contributed by sibship
            sigma_res (float): residual variance component

        Returns:
            float: REML log likelihood
        """
        V: csc_matrix = self._get_V(sigma_grm, sigma_sib, sigma_res)
        lu: SuperLU = splu(V)
        diag_l: np.ndarray = lu.L.diagonal()
        diag_u: np.ndarray = lu.U.diagonal()
        logdet_V: float = np.log(diag_l).sum() + np.log(diag_u).sum()
        e: np.ndarray = np.ones(self.n)
        V_inv_y: np.ndarray = lu.solve(self.y)
        V_inv_e: np.ndarray = lu.solve(e)
        y_T_V_inv_y: float = self.y @ V_inv_y
        e_T_V_inv_e: float = e @ V_inv_e
        e_T_V_inv_y: float = e @ V_inv_y
        logdet_e_T_V_inv_e: float = np.log(e_T_V_inv_e)
        return -0.5 * (logdet_V + logdet_e_T_V_inv_e + y_T_V_inv_y - e_T_V_inv_y ** 2 / e_T_V_inv_e)

    def dense_reml_loglik(self, sigma_grm: float, sigma_sib: float, sigma_res: float,
                          method: Literal['chol', 'cg']) -> float:
        """Dense version of self.reml_loglik.
        
        Support cholesky decomposition and conjugate gradient for matrix inversion.
        """
        logger.info(f'Using {method}.')
        V: csc_matrix = self._get_V(sigma_grm, sigma_sib, sigma_res)
        V_: np.ndarray = V.toarray()
        if method == 'chol':
            del V
        _: int
        logdet_V: float
        V_: np.ndarray
        _, logdet_V = np.linalg.slogdet(V_)
        e: np.ndarray = np.ones(self.n)
        if method == 'chol':
            c, low = cho_factor(V_)
            V_inv_y: np.ndarray = cho_solve((c, low), self.y)
            V_inv_e: np.ndarray = cho_solve((c, low), e)
        elif method == 'cg':
            info1: int
            info2: int
            V_inv_y: float
            V_inv_e: float
            V_inv_y, info1 = cg(V, self.y, atol='legacy')
            if info1 == 0:
                pass
            elif info1 > 0:
                logger.warning('Conjugate gradient did not converge.')
            elif info1 < 1:
                raise RuntimeError(
                    'Illegal input or breakdown for conjugate gradient.')
            V_inv_e, info2 = cg(V, e, atol='legacy')
            if info2 == 0:
                pass
            elif info2 > 0:
                logger.warning('Conjugate gradient did not converge.')
            elif info2 < 1:
                raise RuntimeError(
                    'Illegal input or breakdown for conjugate gradient.')
        y_T_V_inv_y: float = self.y @ V_inv_y
        e_T_V_inv_e: float = e @ V_inv_e
        e_T_V_inv_y: float = e @ V_inv_y
        logdet_e_T_V_inv_e: float = np.log(e_T_V_inv_e)
        return -0.5 * (logdet_V + logdet_e_T_V_inv_e + y_T_V_inv_y - e_T_V_inv_y ** 2 / e_T_V_inv_e)

    def grid_search_reml(self) -> Tuple[float, float, float]:
        """Perform grid search on variance component parameters.

        Returns:
            Tuple[float, float, float]: Tuple of the 3 optimal variance components.
        """
        # reserve jitter for sibma_sib
        upper: float = np.var(self.y) - self.jitter
        step: float = upper / 100.
        max_loglik: float = float('-inf')
        max_sigma: List[float] = [0., 0., 0.]
        logging.info('Starting grid search...')
        i = 1
        for sigma_grm in np.arange(0, upper, step):
            for sigma_sib in np.arange(0, upper, step):
                if sigma_grm + sigma_sib > upper:
                    continue
                sigma_res: float = upper - sigma_grm - sigma_sib
                loglik: float = self.reml_loglik(sigma_grm=sigma_grm, sigma_sib=sigma_sib,
                                                 sigma_res=sigma_sib)
                logging.info(f'{i}: {loglik}')
                i += 1
                if loglik > max_loglik:
                    max_loglik = loglik
                    max_sigma[0] = sigma_grm
                    max_sigma[1] = sigma_sib
                    max_sigma[2] = sigma_res + self.jitter
        logging.info('Finished grid search.')
        return tuple(max_sigma)

