from .sibreg import convert_str_array, get_indices_given_ped, make_id_dict, find_par_gts
import numpy as np
import numpy.ma as ma
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import inv
from bgen_reader import open_bgen
from pysnptools.snpreader import Bed
import h5py
import subprocess
import sys
import tempfile
import gzip
from functools import wraps
from time import time
from typing import List, Optional, Dict, Tuple, Callable, Any, Iterable
from collections import defaultdict
from itertools import combinations


def timethis(f: Callable):
    @wraps(f)
    def wrap(*args, **kwargs):
        ts = time()
        result = f(*args, **kwargs)
        te = time()
        print(
            f'{f.__name__} took: {te - ts:2.4f} sec.'
        )
        return result
    return wrap


def match_phenotype_(
    ids: List[str], y: np.ndarray, pheno_ids: List[str]
):
    """
    Match a phenotype to a genotype array by individual IDs.
    """
    id_dict = make_id_dict(ids)
    in_G_dict = np.array([x in id_dict for x in pheno_ids])
    y = y[in_G_dict]
    pheno_ids = pheno_ids[in_G_dict]
    pheno_id_dict = make_id_dict(pheno_ids)
    y = y[[pheno_id_dict[x] for x in ids]]
    return y


def combo_with_rep(it: Iterable[Any]) -> List[Tuple[Any, Any]]: 
    """
    Return all combinations of length 2 plus duplicates (i.e., [a, a])
    """
    return list(combinations(it, 2)) + [(l, l) for l in it]


class OrderedIdPair(tuple):
    """
    Container to store ordered id pair.
    """
    def __new__(cls, lst):
        if len(lst) != 2:
            raise ValueError('Input must of length 2.')
        if any(not isinstance(i, str) for i in lst):
            raise TypeError('Ids must be string.')
        id1, id2 = sorted(lst)
        return tuple.__new__(cls, [id1, id2])


def get_ids_with_par(
    par_gts_f: str, gts_f: str, ids: List[str] = None
) -> Tuple[List[str],Dict]:
    """
    Find ids with observed/imputed parents and family labels

    Args:
        ids: List[str]
            list of samples ids that we want to analyse; possibly from pheno_ids
    """
    ### Imputed parental file ### 
    par_gts_f = h5py.File(par_gts_f + '.hdf5', 'r')
    ### Genotype file ###
    ped = convert_str_array(np.array(par_gts_f['pedigree']))
    ped = ped[1:ped.shape[0], :]
    # Remove control families
    controls = np.array([x[0] == '_' for x in ped[:, 0]])
    ped = ped[np.logical_not(controls),:]
    # Get families with imputed parental genotypes
    fams = convert_str_array(np.array(par_gts_f['families']))

    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f = Bed(gts_f, count_A1=True)
        gts_ids = gts_f.iid[:, 1]
        ids, observed_indices, imp_indices = get_indices_given_ped(ped, fams, gts_ids, ids=ids)
        gts_ids = gts_f.iid[observed_indices, 1]
    elif gts_f[(len(gts_f) - 5):len(gts_f)]  == '.bgen':
        gts_f = open_bgen(gts_f)
        gts_ids = gts_f.samples
        ids, observed_indices, imp_indices = get_indices_given_ped(ped, fams, gts_ids, ids=ids)
        gts_ids = gts_ids[observed_indices]
    else:
        raise ValueError('Unknown filetype for observed genotypes file: ' + str(gts_f))
    
    fams = fams[imp_indices]
    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, fams, gts_id_dict)

    return ids, fam_labels


def run_gcta_grm(
    plink_path: str, gcta_path: str, filename: str, output_path: str, keep: Optional[List] = None
) -> None:
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
    print('------------- start running gcta grm ----------------')
    args = [gcta_path, ]

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
                        plink_path, '--bfile', c, '--rm-dup', 'exclude-mismatch', 
                        '--make-bed', '--out', f'{tmpdir}/chr{str(i)}'
                    ]
                    if keep is not None:
                        plink_args += ['--keep', keep_file]
                    subprocess.run(plink_args)
                    f.write(f'{tmpdir}/chr{str(i)}' + '\n')
            args.append(chr_file)
        else:
            plink_args = [
                plink_path, '--bfile', filename, '--rm-dup', 'exclude-mismatch', 
                '--make-bed', '--out', f'{tmpdir}/chr'
            ]
            if keep is not None:
                plink_args += ['--keep', keep_file]
            subprocess.run(plink_args)
            args += ['--bfile', f'{tmpdir}/chr']
        
        args += ['--make-grm-gz', '--out', output_path]
        subprocess.run(args)
    print('------------- finished running gcta grm --------------')


def make_id_pair_dict(ids: List[str]) -> Dict[OrderedIdPair, int]:
    """
    Make ordered ID pair dict from list of IDs.
    """
    pairs = combo_with_rep(ids)
    id_pair_dict = dict()
    i = 0
    for p in map(OrderedIdPair, pairs):
       id_pair_dict[p] = i
       i += 1
    return id_pair_dict


def build_grm_arr(
    grm_path: str, he_id_dict: Dict[OrderedIdPair, int], thres: float = 0.05
) -> np.ndarray:
    """
    Build GRM data array for HE regression.
    """
    gcta_id = dict()

    row_n = 1  # index starts from 1 in gcta grm output
    with open(f'{grm_path}.id', 'r') as f:
        for line in f:
            id = line.strip().split('\t')[1]
            gcta_id[row_n] = id
            row_n += 1
    grm_arr = np.full(len(he_id_dict), np.nan)

    with gzip.open(f'{grm_path}.gz', 'rt') as f:
        for line in f:
            id1, id2, _, gr = line.strip().split('\t')
            id1, id2 = gcta_id[int(id1)], gcta_id[int(id2)]
            gr = float(gr)
            id_pair = OrderedIdPair([id1, id2])
            try:
                grm_arr[he_id_dict[id_pair]] = gr if gr >= thres else 0.
            except IndexError:
                sys.exit('Indices in he_id_dict should probably start from 0.')
            except KeyError:
                print(f'ID pair {id_pair} in GCTA GRM output is not requested in he_id_dict.')
                # sys.exit(f'ID pair {id_pair} in GCTA GRM output is not requested in he_id_dict.')

    if np.isnan(grm_arr).any():
        raise ValueError('Fewer pairs in GCTA GRM outout than expected.')
    return grm_arr


def build_sibship_arr(
    fam_labels: List[str], ids: List[str], he_id_dict: Dict[OrderedIdPair, int], 
) -> np.ndarray:
    """
    Build sibship array for HE regression.

    Args:
        fam_labels: List[str]
            List of family id strings corresponding to each individual
    """
    sibship_arr = np.full(len(he_id_dict), 0.)

    label_indices = defaultdict(list)
    for l, f in enumerate(fam_labels):
        label_indices[f].append(l)
    
    for f, indices_lst in label_indices.items():
        if len(indices_lst) == 1:
            continue
        for pair in combo_with_rep(indices_lst):
            id_pair = list(map(ids.__getitem__, pair))
            try:
               sibship_arr[he_id_dict[OrderedIdPair(id_pair)]] = 1.
            except KeyError:
               sys.exit(f'ID pair {id_pair} from fam_indices is not repuested in he_id_dict.')
    
    return sibship_arr


def build_res_arr(
    he_id_dict: Dict[OrderedIdPair, int], 
) -> np.ndarray:
    """
    Build residual array for HE regression.
    """
    return np.array([int(i == j) for i, j in he_id_dict], dtype=np.float)


def build_pheno_prod_arr(
    y: np.ndarray, ids: List[str], he_id_dict: Dict[OrderedIdPair, int]
) -> np.ndarray:
    """
    Build phenotype product array.
    """
    y = y - y.mean()
    pheno_prod_arr = np.full(len(he_id_dict), np.nan)
    id_dict = make_id_dict(ids)

    for id_pair, ind in he_id_dict.items():
        id1, id2 = id_pair
        prod = y[id_dict[id1]] * y[id_dict[id2]]
        pheno_prod_arr[ind] = prod
    
    if np.isnan(pheno_prod_arr).any():
        raise ValueError('Fewer pairs than expected.')
    return pheno_prod_arr


def he_reg(
    grm_arr: np.ndarray, sibship_arr: np.ndarray, res_arr: np.ndarray, pheno_prod_arr: np.ndarray
) -> Tuple[float, float, float]:
    """
    Perform HE regression.
    """
    if grm_arr.ndim > 1 or sibship_arr.ndim > 1 or pheno_prod_arr.ndim > 1:
        raise TypeError('Wrong input dimension. Expecting ndarrays of dim 1.')
    if len(grm_arr) != len(sibship_arr) or len(sibship_arr) != len(pheno_prod_arr):
        raise ValueError('Lengths of inputs do not match.')

    X = np.vstack([np.ones(len(grm_arr)), grm_arr, sibship_arr]).T
    beta = np.linalg.solve(X.T.dot(X), X.T.dot(pheno_prod_arr))
    res = pheno_prod_arr - X @ beta
    err = res.T @ res / (pheno_prod_arr.shape[0] - 2)
    beta1 = 0 if beta[1] < 0 else beta[1]
    beta2 = 0 if beta[2] < 0 else beta[2]
    return beta1, beta2, err


class LinearMixedModel(object):
    """Define a linear model with within-class correlations.

    Args:
        y: np.ndarray
            1D array of phenotype observations
        sigma_grm: float
            variance component contributed by GRM
        sigma_sib: float
            variance component contributed by sibship
        sigma_res: float
            redidual variance component

    Returns:
        model : :class:`sibreg.model`

    """

    def __init__(
        self,
        y: np.ndarray,
        # ids: List[str], fam_labels: List[str], 
        sigma_grm: float, sigma_sib: float, sigma_res: float,
        include_covar: bool = True
    ) -> None:
        # if len(ids) != len(fam_labels):
        #     raise ValueError('ids and fam_labels should have equal lengths.')
        if y.ndim != 1:
            raise ValueError('y should be a 1-d array.')

        self.y = y
        self.include_covar = include_covar
        if include_covar:
           self._covar_adjusted = False
        # self.ids = ids
        # self.fam_labels = fam_labels
        self.n = len(self.y)
        self._sigma_grm = sigma_grm
        self._sigma_sib = sigma_sib
        self._sigma_res = sigma_res

     
    def fit_covar(self, covar_X: np.ndarray) -> np.ndarray:
        if not self.include_covar:
            print('This model does not include covariates.')
            return
        if getattr(self, '_covar_adjusted'):
            print('Covariates have been adjusted for.')
            return

        if covar_X.ndim != 2:
            raise TypeError(
                'Wrong input dimension.'
            )
        if self.y.shape[0] != covar_X.shape[0]:
            raise TypeError(
                f'Input dimensions do not match. y: {self.y.shape}. covar_X: {covar_X.shape}.'
            )
        self.y = self.y - covar_X @ np.linalg.inv(covar_X.T @ covar_X) @ covar_X.T @ self.y
        self.covar_adjusted = True
    

    def compute_V_inv(
        self,
        ids: Dict[str, int],
        id_pair_dict: Dict[OrderedIdPair, int],
        grm_arr: np.ndarray, sibship_arr: np.ndarray
    ) -> None:
        """
        Construct covariance matrix V, compute its inverse and its product with self.y
        """
        if self.include_covar and not getattr(self, '_covar_adjusted'):
            raise RuntimeError(
                'Please adjust phenotype for covariates first.'
            )
        if grm_arr.ndim > 1 or sibship_arr.ndim > 1:
            raise TypeError(
                'Wrong input dimension. grm_arr and sibship_arr should have dimension (1, ).'
            )
        if len(id_pair_dict) != len(grm_arr) or len(grm_arr) != len(sibship_arr):
            raise ValueError(
                'Lengths of inputs should match.'
            )
        
        id_dict = make_id_dict(ids)

        # split into diag and off-diag
        data = self._sigma_grm * grm_arr + self._sigma_sib * sibship_arr
        off_data = data[np.array([ind for (i, j), ind in id_pair_dict.items() if i != j])]
        off_row_ind = np.array([id_dict[i] for i, _ in id_pair_dict if i != _])
        off_col_ind = np.array([id_dict[j] for _, j in id_pair_dict if _ != j])
        diag_data = data[np.array([ind for (i, j), ind in id_pair_dict.items() if i == j])]
        diag_row_ind = np.array([id_dict[i] for i, _ in id_pair_dict if i == _])
        diag_col_ind = np.array([id_dict[j] for _, j in id_pair_dict if _ == j])

        # keep only nonzero entries
        off_nonzero = off_data.nonzero()
        off_data, off_row_ind, off_col_ind = \
            off_data[off_nonzero], off_row_ind[off_nonzero], off_col_ind[off_nonzero]
        
        # combine diag and off-diag
        data = np.concatenate([off_data, off_data, diag_data + self._sigma_res])
        row_ind, col_ind = \
            np.concatenate([off_row_ind, off_col_ind, diag_row_ind]), np.concatenate([off_col_ind, off_row_ind, diag_col_ind])
        
        # construct sparse V
        V = csc_matrix((data, (row_ind, col_ind)), shape=(self.n, self.n))
        self._Vinv = inv(V)
        self._Vinv_y = self._Vinv.dot(self.y)


    @staticmethod
    def sp_mul_dense3d(
        sps_mat: csc_matrix, dense_mat: np.ndarray
    ) -> np.ndarray:
        """
        Compute product of given sparse matrix and 3d dense matrix
        """
        if dense_mat.ndim != 3:
            raise ValueError(
                'dense_mat must be a 3-d array.'
            )
        n1, n2 = sps_mat.shape
        m, n, c = dense_mat.shape
        if n != n2:
            raise ValueError(
                'Input dims do not match.'
            )

        out = np.empty((m, n1, c))
        for i in range(m):
            out[i, :, :] = sps_mat.dot(dense_mat[i, :, :])
        return out


    def fit_snps_eff(
        self,
        gts: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Perform repeated OLS to estimate SNP effects and sampling variance-covariance in transformed model
        """
        if self.include_covar and not getattr(self, '_covar_adjusted'):
            raise RuntimeError(
                'Please adjust phenotype for covariates first.'
            )
        if gts.ndim != 3:
            raise ValueError(
                'gts is expected to be of dimension 3.'
            )
        if not hasattr(self, '_Vinv'):
            raise RuntimeError(
                'V inverse is not computed.'
            )

        gts = gts.transpose(2, 0, 1)
        Vinv_X = self.sp_mul_dense3d(self._Vinv, gts)
        X_T_Vinv_X = np.einsum('...ij,...ik', gts, Vinv_X)
        X_T_Vinv_y = np.einsum('...ij,i', gts, self._Vinv_y)
        alpha = np.linalg.solve(X_T_Vinv_X, X_T_Vinv_y)
        alpha_cov = np.linalg.inv(X_T_Vinv_X)
        alpha_ses = np.sqrt(np.diagonal(alpha_cov, axis1=1, axis2=2))
        return alpha, alpha_cov, alpha_ses