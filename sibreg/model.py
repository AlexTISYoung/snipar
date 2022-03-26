import code
from numpy.linalg import solve
from .sibreg import convert_str_array, get_indices_given_ped, gtarray, make_id_dict, find_par_gts
import numpy as np
import numpy.ma as ma
from numpy.linalg import slogdet
from scipy.sparse import csc_matrix, tril
from scipy.sparse.linalg import splu, SuperLU, spsolve, cg, spsolve_triangular
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
from itertools import product, combinations_with_replacement, combinations
from functools import lru_cache
import pandas as pd
from operator import attrgetter
from multiprocessing import Pool, RawArray
import ctypes


# FORMAT = '%(asctime)-15s :: %(levelname)s :: %(filename)s :: %(funcName)s :: %(message)s'
# numeric_level = getattr(logging, loglevel.upper(), None)
# if not isinstance(numeric_level, int):
#     raise ValueError('Invalid log level: %s' % loglevel)
# logging.basicConfig(
#     format=FORMAT, level=logging.DEBUG if __debug__ else logging.INFO)
logger = logging.getLogger(__name__)


# https://stackoverflow.com/questions/48262273/python-bookkeeping-dependencies-in-cached-attributes-that-might-change#answer-48264126
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


def match_phenotype_(ids: np.ndarray, y: np.ndarray, pheno_ids: List[str]):
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
                     ids: np.ndarray = None,
                     sib: bool = False,
                     impute_unrel: bool = True,
                     return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Find ids with observed/imputed parents and family labels
    """
    # Imputed parental file
    par_gts_f_ = h5py.File(par_gts_f, 'r')
    # Genotype file
    ped = convert_str_array(np.array(par_gts_f_['pedigree']))
    ped = ped[1:ped.shape[0], :]
    logger.debug(f'ped {len(ped)}')
    # Remove control families
    controls = np.array([x[0] == '_' for x in ped[:, 0]])
    ped = ped[np.logical_not(controls), :]
    logger.debug(f'ped {len(ped)}')
    # Get families with imputed parental genotypes
    fams = convert_str_array(np.array(par_gts_f_['families']))

    # logger.debug(f'#ids present in both phen and ped: {len(np.intersect1d(ids, ped[:, 1]))}')

    if gts_f[(len(gts_f) - 4):len(gts_f)] == '.bed':
        gts_f_: Bed = Bed(gts_f, count_A1=True)
        gts_ids = gts_f_.iid[:, 1]
        logger.debug(f'Length of geno: {len(gts_ids)}')
        # logger.debug('#iid' + len(gts_ids))
        # logger.debug(f'#ids present in both phen and geno: {len(np.intersect1d(ids, gts_ids))}')
        logger.debug(f'#ids present in both ped and geno: {len(np.intersect1d(ped[:, 1], gts_ids))}')
        ids, observed_indices, imp_indices = get_indices_given_ped(ped,
                                                                   fams,
                                                                   gts_ids,
                                                                   ids=ids,
                                                                   sib=sib,
                                                                   impute_unrel=impute_unrel)
        gts_ids = gts_f_.iid[observed_indices, 1]
        logger.debug(f'Length of observed geno: {len(gts_ids)}')
    elif gts_f[(len(gts_f) - 5):len(gts_f)] == '.bgen':
        gts_f_ = open_bgen(gts_f)
        gts_ids = gts_f_.samples
        logger.debug(f'Length of geno: {len(gts_ids)}')
        # logger.debug(f'#ids present in both phen and geno: {len(np.intersect1d(ids, gts_ids))}')
        logger.debug(f'#ids present in both ped and geno: {len(np.intersect1d(ped[:, 1], gts_ids))}')
        ids, observed_indices, imp_indices = get_indices_given_ped(ped,
                                                                   fams,
                                                                   gts_ids,
                                                                   ids=ids,
                                                                   sib=sib,
                                                                   impute_unrel=impute_unrel)
        gts_ids = gts_ids[observed_indices]
        logger.debug(f'Length of observed geno: {len(gts_ids)}')
    else:
        raise ValueError('Unknown filetype for observed genotypes file: ' +
                         str(gts_f))

    fams = fams[imp_indices]
    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, fams,
                                                      gts_id_dict)
    # print(len(par_status), len(observed_indices))
    # print(sum([i == j == -1 for i, j in par_status]))
    # print(sum([i == j == -1 for _,i, j in gt_indices]))
    # print(sum([1 for i in fam_labels if len(i) == 0]))
    # ind = np.array([i for i in range(6797) if par_status[i, 0] == par_status[i, 1] == -1 and len(fam_labels[i]) != 0])
    # print(par_status[ind], fam_labels[ind], gt_indices[ind])
    # f = fam_labels[ind]
    # print([ped[i, :] for i in range(len(ped)) if ped[i, 1] in f])
    logger.debug(f'Length of par status: {len(par_status)}')
    logger.debug(f'#missing [F, M]: {(par_status == -1).sum(axis=0)}')
    logger.debug(f'#missing F and M (from par_status): {sum([i == j == -1 for i, j in par_status])}')
    logger.debug(f'#missing F and M (fro gt_indices): {sum([i == j == -1 for _,i, j in gt_indices])}')
    logger.debug(f'#both par geno: {sum([i == j == 0 for _,i, j in gt_indices])}')
    logger.debug(f'#both par geno: {sum([i == j == 0 for i, j in par_status])}')
    logger.debug(f'#one par geno: {sum([i + j == 1 for i, j in par_status])}')
    n_empty_fams = sum([1 for i in fam_labels if len(i) == 0])
    logger.debug(f'#unique fam_labels: {np.unique(fam_labels).__len__()}')
    logger.debug(f'#empty fam_labels: {n_empty_fams}')

    logger.debug(fam_labels[fam_labels == ''])
    # ind = fam_labels == ''
    fam_labels = fam_labels.astype('<U20')
    fam_labels[fam_labels == ''] = np.array([f'_not_{i}_' for i in range(n_empty_fams)], dtype='<U20')
    # ind = np.array([i[4] == 'False' and i[5] == 'False' for i in ped])
    logger.debug(f'#empty fam_labels after assignment: {sum([1 for i in fam_labels if len(i) == 0])}')
    logger.debug(f'#unique fam_labels: {np.unique(fam_labels).__len__()}')
    logger.debug(f'#iids to use {len(ids)}')
    xx = len([1 for i in fam_labels if '_not_' in i])
    logger.debug(f'{xx}')
    if return_info:
        return ids, fam_labels, par_status, gt_indices
    return ids, fam_labels


def match_grm_ids(ids: np.ndarray,
                  fam_labels: np.ndarray,
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


def build_grm_arr_(grm_path: str, id_dict: Dict[Hashable, int],
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
            if id1 not in id_dict or id2 not in id_dict:
                # ignore individuals that are not in id_dict
                continue
            ind1, ind2 = id_dict[id1], id_dict[id2]
            arr_ind = coord2linear(ind1, ind2)
            gr_: float = float(gr)
            grm_arr[arr_ind] = gr_ if gr_ >= thres else 0.

    if np.isnan(grm_arr).any():
        raise ValueError('Fewer pairs in GCTA GRM outout than expected.')
    return grm_arr


def build_ibdrel_arr_(ibdrel_path: str, id_dict: Dict[Hashable, int],
                     keep: List[Union[str, int]], thres: float = 0.05) -> np.ndarray:
    """Build ibd relatedness array (lower triangular entries) from KING ibdseg output.

    Args:
        ibdrel_path (str): Path to ibdseg output.
        id_dict (Dict[Hashable, int]): dictionary of id-index pairs.
        keep (List[Union[str, int]]): list of ids to keep.
        thres (float, optional): sparsity threshold. Defaults to 0.0205.

    Returns:
        np.ndarray: 1-d ibdrel array.
    """
    logger.info(f'Reading {ibdrel_path}.seg...')
    king = pd.read_csv(ibdrel_path + '.seg',
                       sep='\t')[['ID1', 'ID2', 'PropIBD']]
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


def build_grm_arr(grm_path: str, id_dict: Dict[Hashable, int],
                  thres: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    grm_arr = np.full(int(n * (n + 1) / 2), np.nan)
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


def build_grm_arr_from_npz(id_filepath: str, npz_path: str, ids: np.ndarray, id_dict: Dict[Hashable, int]):
    logger.info('Building grm from npz...')
    _ids = pd.read_csv(id_filepath, sep="\s+", header=None)[1].astype(str).to_numpy()
    f = np.load(npz_path)
    orig = pd.DataFrame({'data': np.float64(f['data']), 'row_ids': _ids[f['row_ind']], 'col_ids': _ids[f['col_ind']]})
    _ids = pd.read_csv(id_filepath, sep="\s+", header=None)[1].astype(str)
    _ids = _ids[~_ids.isin(ids)]
    orig = orig[orig['row_ids'].isin(ids)]
    orig = orig[orig['col_ids'].isin(ids)]
    orig['row_ind'] = orig['row_ids'].apply(lambda x: id_dict[x])
    orig['col_ind'] = orig['col_ids'].apply(lambda x: id_dict[x])
    # https://stackoverflow.com/questions/45504391/swapping-column-values-based-on-column-conditions-pandas-dataframe
    cond = orig['row_ind'] < orig['col_ind']
    orig.loc[cond, ['row_ind', 'col_ind']] = orig.loc[cond, ['col_ind', 'row_ind']].values
    data = orig['data'].to_numpy()
    row_ind = orig['row_ind'].to_numpy(dtype='uint32')
    col_ind = orig['col_ind'].to_numpy(dtype='uint32')
    logger.info('Done building grm from npz.')
    return data, row_ind, col_ind


def build_ibdrel_arr(ibdrel_path: str, id_dict: Dict[Hashable, int],
                     keep: List[Union[str, int]], 
                     thres: float = 0.05) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build ibd relatedness array (lower triangular entries) and corresponding indices from KING ibdseg output.

    Args:
        ibdrel_path (str): Path to ibdseg output.
        id_dict (Dict[Hashable, int]): dictionary of id-index pairs.
        keep (List[Union[str, int]]): list of ids to keep.
        thres (float, optional): sparsity threshold. Defaults to 0.0205.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: a tuple of ndarrays containing data and indices.
    """
    logger.info(f'Reading {ibdrel_path}.seg...')
    king = pd.read_csv(ibdrel_path + '.seg',
                       sep='\t')[['ID1', 'ID2', 'PropIBD']]
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
    logger.info('Building ibd arr...')
    for row in king.itertuples():
        id1 = row.ID1
        id2 = row.ID2
        ind1, ind2 = id_dict[id1], id_dict[id2]
        ind1, ind2 = max(ind1, ind2), min(ind1, ind2)
        data[row.Index] = row.PropIBD  # if row.PropIBD > thres else 0.
        row_ind[row.Index] = ind1
        col_ind[row.Index] = ind2
    for i in range(len(keep)):
        data[n + i] = 1.
        row_ind[n + i] = i
        col_ind[n + i] = i
    logger.info('Done building ibd arr.')
    return np.array(data), np.array(row_ind, dtype='uint32'), np.array(col_ind, dtype='uint32')


def write_ibdrel_to_grm(ibdrel_path: str, output_path: str,
                        ids: np.ndarray, thres: float = 0.0205) -> None:
    """Write ibd relatedness from KING to GCTA GRM format.
    """
    logger.info(f'Reading {ibdrel_path}...')
    king = pd.read_csv(ibdrel_path, sep='\t')[['ID1', 'ID2', 'PropIBD']]
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
    # print(full_king)
    def make_ind(df):
        rows, cols = [], []
        for row in df.itertuples():
            id1 = row.ID1
            id2 = row.ID2
            ind1, ind2 = id_dict[id1], id_dict[id2]
            row_ind, col_ind = max(ind1, ind2), min(ind1, ind2)
            rows.append(row_ind + 1)  # GCTA indices start from 1
            cols.append(col_ind + 1)
        return rows, cols  # pd.Series({'row_ind': rows, 'col_ind': cols})

    logger.info('Making ibd row and col indices...')
    # king[['row_ind', 'col_ind']] = make_ind(king)
    rows, cols = make_ind(king)
    logger.info('Finished making ibd row and col indices...')
    # king['row_ind'] = np.array(rows, dtype=np.int)
    # king['col_ind'] = np.array(cols, dtype=np.int)
    king = king.assign(row_ind=np.array(rows, dtype=np.int), col_ind=np.array(cols, dtype=np.int))
    king_ = king.merge(full_king, on=['row_ind', 'col_ind'], how='right')[[
        'row_ind', 'col_ind', 'non_missing', 'PropIBD'
    ]]
    del king, full_king
    king_ = king_.fillna(0.)
    king_.loc[king_.row_ind == king_.col_ind, 'PropIBD'] = 1.
    # sparsification
    # king_.loc[king_['PropIBD'] < thres, 'PropIBD'] = 0.
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


def build_sib_arr_(fam_labels: List[str], y) -> np.ndarray:
    """Build sibship array for HE regression.

    Args:
        fam_labels (List[str]): List of family id strings corresponding to each individual.

    Returns:
        np.ndarray: lower triangular entries of sibship matrix.
    """
    logger.info('Building sibship arr...')
    n = len(fam_labels)
    sib_arr = np.zeros(n_tril(n))
    # sib_arr = np.zeros(n_tril(n) - n)
    # # y_arr_ = np.zeros(n_tril(n))
    # y = y - y.mean()
    # y_i = np.tril_indices(n, k=-1)
    # y_arr = y[y_i[0]] * y[y_i[1]]
    label_indices = defaultdict(list)
    for l, f in enumerate(fam_labels):
        label_indices[f].append(l)

    for f, indices_lst in label_indices.items():
        # for pair in combinations(indices_lst, 2):
        for pair in combinations_with_replacement(indices_lst, 2):
            arr_ind = coord2linear(*pair)
            sib_arr[arr_ind] = 1.
            # y_arr_[arr_ind] = y_arr[arr_ind]
    logger.info('Done building sibship arr.')
    # v = y_arr.var()
    # print(v, y_arr.dot(sib_arr).sum() / y_arr.std() / sib_arr.std() / len(y_arr))
    # return sib_arr, y_arr
    return sib_arr


def build_sib_arr(fam_labels: List[str]) -> np.ndarray:
    """Build lower-triangular nonzero entries of sibship matrix.

    Args:
        fam_labels (List[str]): List of family id strings corresponding to each individual.

    Returns:
        np.ndarray: lower triangular entries of sibship matrix.
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


def impute_unrel_par_gts(G: gtarray, sib: bool = False, parsum: bool = False) -> gtarray:
    """Impute parental genotypes of unrelated individuals.

    Args:
        G (gtarray): gtarray object holding genetic data.
        sib (bool, optional): Whether sib effect is modelled. Defaults to False.
        parsum (bool, optional): Whether sum of parental genotypes is modelled. Defaults to False.

    Returns:
        gtarray: gtarray object with imputed parental genotypes.
    """    
    if sib:
        logger.warning(
            'No need to impute parental genotypes of unrelated individuals if sib effect is modelled.'
        )
        return G
    # if G.freqs is None:
    #     G.compute_freqs()
    # missing_ind = np.min(gt_indices, axis=1)
    # missing_ind = none_missing_ind == -1
    missing_ind = G.gts[:, 1, 0] == -1
    # G.compute_freqs(G.gts[:, 1, 0] > -1)
    np.testing.assert_array_equal(missing_ind, G.gts[:, 1, 1] == -1, err_msg='Missing indices don\'t match across snps.')
    np.testing.assert_array_equal(missing_ind, G.gts[:, 1, 2] == -1, err_msg='Missing indices don\'t match across snps.')
    logger.debug(f'Number of unrelated individuals that need linear imputations: {missing_ind.sum()}.')
    imp = 0.5 * G.gts[missing_ind, 0, :] + G.freqs
    G.gts[missing_ind, 1, :] = imp
    if parsum:
        G.gts[missing_ind, 1, :] += imp
    else:
        if __debug__:
            np.testing.assert_array_equal(missing_ind, G.gts[:, 2, 0] == -1)
        G.gts[missing_ind, 2, :] = imp
    logger.debug(f'Number of missingness left after linear imputation: {(G.gts[:, 1, 0] == -1).sum()}.')
    return G


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


def compute_neff_ratio(sib_corr: float, f: float, n0: int, n1: int, phased: bool = False) -> float:
    """Compute theoretical effective sample size ratio of combined sample analysis versus related sample analysis.

    Args:
        sib_corr (float): sibling correlation.
        f (float): allele frequency.
        n0 (int): number of sib pairs in related sample.
        n1 (int): unrelated sample size.

    Returns:
        float: effective sample size ratio.
    """
    if phased:
        nu = 3 / 4
    else:
        nu = 3 / 4 - f * (1 - f) / 8.
    var_direct = np.linalg.inv(
        n0 / (1 - sib_corr ** 2) * np.array(
            [[2 - sib_corr, 2 * (1 - sib_corr)],
            [2 * (1 - sib_corr), 4 * nu * (1 - sib_corr)]]
        ) + n1
    )[0, 0]
    var_direct0 = 4 * nu * (1 - sib_corr ** 2) / (4 * n0 * (nu * (2 - sib_corr) + sib_corr - 1))
    return var_direct0 / var_direct


class GradHessComponents(NamedTuple):
    P_y: np.ndarray
    P_varcomp_mats: Tuple[np.ndarray, ...]


_var_dict = {}


NUM_CPU = 8


def init_worker(L_data_, L_indices_, L_indptr_,
                U_data_, U_indices_, U_indptr_,
                dense_mat_, dense_mat_shape,
                perm_c_, perm_r_):
    _var_dict['L_data_'] = L_data_
    _var_dict['L_indices_'] = L_indices_
    _var_dict['L_indptr_'] = L_indptr_
    _var_dict['U_data_'] = U_data_
    _var_dict['U_indices_'] = U_indices_
    _var_dict['U_indptr_'] = U_indptr_
    _var_dict['dense_mat_'] = dense_mat_
    _var_dict['dense_mat_shape'] = dense_mat_shape
    _var_dict['perm_c_'] = perm_c_
    _var_dict['perm_r_'] = perm_r_


def worker_func(ind_lst):
    def spsolve_lu(L, U, b, perm_c=None, perm_r=None):
        if perm_r is not None:
            b_old = b.copy()
            for old_ndx, new_ndx in enumerate(perm_r):
                b[new_ndx] = b_old[old_ndx]
        try:
            c = spsolve_triangular(L, b, lower=True, unit_diagonal=True)
        except TypeError:
            c = spsolve_triangular(L, b, lower=True)
        px = spsolve_triangular(U, c, lower=False)
        if perm_c is None:
            return px
        return px[perm_c]
    L_data = np.frombuffer(_var_dict['L_data_'])
    L_indices = np.frombuffer(_var_dict['L_indices_'], dtype='int32')
    L_indptr = np.frombuffer(_var_dict['L_indptr_'], dtype='int32')
    U_data = np.frombuffer(_var_dict['U_data_'])
    U_indices = np.frombuffer(_var_dict['U_indices_'], dtype='int32')
    U_indptr = np.frombuffer(_var_dict['U_indptr_'], dtype='int32')
    dense_mat = np.frombuffer(_var_dict['dense_mat_']).reshape(
        _var_dict['dense_mat_shape'])
    perm_c = np.frombuffer(_var_dict['perm_c_'], dtype='int32')
    perm_r = np.frombuffer(_var_dict['perm_r_'], dtype='int32')
    L = csc_matrix((L_data, L_indices, L_indptr), shape=(
        dense_mat.shape[1], dense_mat.shape[1])).tocsr()
    U = csc_matrix((U_data, U_indices, U_indptr), shape=(
        dense_mat.shape[1], dense_mat.shape[1])).tocsr()
    out = np.empty((len(ind_lst), L.shape[0], dense_mat.shape[2]))
    for i, ind in enumerate(ind_lst):
        start = time()
        out[i, :, :] = spsolve_lu(
            L, U, dense_mat[ind, :, :], perm_c=perm_c, perm_r=perm_r)
        print(time() - start)
    return out


def init_worker_(data_, indices_, indptr_, shape,
                 dense_mat_, dense_mat_shape):
    _var_dict['data_'] = data_
    _var_dict['indices_'] = indices_
    _var_dict['indptr_'] = indptr_
    _var_dict['dense_mat_'] = dense_mat_
    _var_dict['dense_mat_shape'] = dense_mat_shape
    _var_dict['shape'] = shape


def worker_func_(ind_lst):
    data = np.frombuffer(_var_dict['data_'])
    indices = np.frombuffer(_var_dict['indices_'], dtype='int32')
    indptr = np.frombuffer(_var_dict['indptr_'], dtype='int32')
    dense_mat = np.frombuffer(_var_dict['dense_mat_']).reshape(
        _var_dict['dense_mat_shape'])
    V = csc_matrix((data, indices, indptr), shape=_var_dict['shape'])
    lu = splu(V)
    out = np.empty((len(ind_lst), V.shape[0], dense_mat.shape[2]))
    for i, ind in enumerate(ind_lst):
        for j in range(dense_mat.shape[2]):
            out[i, :, j] = lu.solve(dense_mat[ind, :, j])
    logger.debug(f'Start sending result: {out.nbytes} bytes...')
    return out


class LinearMixedModel:
    """Wrapper of data and functions that compute estimates of variance components or SNP effects.
    """
    logger = logger.getChild(__qualname__)

    # jitter = 1e-6  # 1e-3
    _vec_linear2coord: Callable[..., Tuple[np.ndarray, np.ndarray]] = np.vectorize(
        linear2coord)
    _vec_coord2linear: Callable[..., np.ndarray] = np.vectorize(coord2linear)

    def __init__(self, y: np.ndarray, 
                 varcomp_arr_lst: Tuple[Tuple[np.ndarray, np.ndarray, np.ndarray], ...],
                 varcomps: Tuple[float, ...] = None,
                 covar_X: np.ndarray = None,
                 add_intercept: bool = False) -> None:
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
            if covar_X.ndim == 1:
                covar_X = covar_X.reshape((self.n,1))
            if covar_X.ndim > 2:
                raise ValueError('covar_X should be a 1-d or 2-d array.')
            if covar_X.shape[0] != y.shape[0]:
                raise ValueError(
                    f'Dimensions of y and covar_X do not match ({y.shape} and {covar_X.shape}).')
            self.Z: np.ndarray = ma.getdata(covar_X)
            self.has_covar = True
            # self.logger.info('Adjusting phenotype for fixed covariates...')
            # self.y = self.fit_covar(y, covar_X)
            # self.logger.info('Finished adjusting.')
        # self.n: int = len(y)
        self.logger.info(f'#individuals: {self.n}.')

        if add_intercept:
            self.Z = np.hstack((np.ones((self.n, 1),dtype=self.Z.dtype), self.Z))
            
        def to_nz(arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            nz = arr.nonzero()[0]
            return arr[nz], nz
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
        y_var: float = self.y.var()
        self.logger.info(f'Phenotypic variance: {y_var}.')
        self.jitter = y_var / 1000.
        if varcomps is not None:
            if len(varcomps) != self.n_varcomps:
                raise ValueError('varcomps and varcomps_arr_lst must have the same length.')
            self._varcomps = varcomps 
            self.optimized = True
        else:
            self._varcomps = tuple(
                # self.y.var() if i == self.n_varcomps - 1 else 0. for i in range(self.n_varcomps))
                self.y.var() / self.n_varcomps for i in range(self.n_varcomps))
            self.logger.info(f'init varcomps: {self._varcomps}.')
            # self._varcomps = tuple(
            #     y_var / self.n_varcomps for i in range(self.n_varcomps))
            self.optimized = False
        self.logger.info(f'Number of variance components: {self.n_varcomps}.')

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
    @staticmethod
    def sp_solve_dense3d_lu_mp(sps_mat: csc_matrix,
                               dense_mat: np.ndarray) -> np.ndarray:
        """Compute product of the inverse of given sparse matrix and given 3-d dense array using LU in parallel; used for repeated OLS.
        Each SNP requires one LU factorization.

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
        lu = timethis(splu)(sps_mat)
        L, U = lu.L, lu.U

        L_data, L_indices, L_indptr = L.data, L.indices, L.indptr
        U_data, U_indices, U_indptr = U.data, U.indices, U.indptr
        L_data_shape = L_data.shape
        L_indices_shape = L_indices.shape
        L_indptr_shape = L_indptr.shape
        U_data_shape = U_data.shape
        U_indices_shape = U_indices.shape
        U_indptr_shape = U_indptr.shape
        L_data_ = RawArray('d', L_data_shape[0])
        L_data_buffer = np.frombuffer(L_data_)
        np.copyto(L_data_buffer, L_data)
        L_indices_ = RawArray(ctypes.c_int32, L_indices_shape[0])
        L_indices_buffer = np.frombuffer(L_indices_, dtype='int32')
        np.copyto(L_indices_buffer, L_indices)
        L_indptr_ = RawArray(ctypes.c_int32, L_indptr_shape[0])
        L_indptr_buffer = np.frombuffer(L_indptr_, dtype='int32')
        np.copyto(L_indptr_buffer, L_indptr)

        U_data_ = RawArray('d', U_data_shape[0])
        U_data_buffer = np.frombuffer(U_data_)
        np.copyto(U_data_buffer, U_data)
        U_indices_ = RawArray(ctypes.c_int32, U_indices_shape[0])
        U_indices_buffer = np.frombuffer(U_indices_, dtype='int32')
        np.copyto(U_indices_buffer, U_indices)
        U_indptr_ = RawArray(ctypes.c_int32, U_indptr_shape[0])
        U_indptr_buffer = np.frombuffer(U_indptr_, dtype='int32')
        np.copyto(U_indptr_buffer, U_indptr)

        perm_c_ = RawArray(ctypes.c_int32, lu.perm_c.shape[0])
        perm_c_buffer = np.frombuffer(perm_c_, dtype='int32')
        np.copyto(perm_c_buffer, lu.perm_c)
        perm_r_ = RawArray(ctypes.c_int32, lu.perm_r.shape[0])
        perm_r_buffer = np.frombuffer(perm_r_, dtype='int32')
        np.copyto(perm_r_buffer, lu.perm_r)

        dense_mat_shape = dense_mat.shape
        dense_mat_ = RawArray(
            'd', dense_mat_shape[0] * dense_mat_shape[1] * dense_mat_shape[2])
        dense_mat_buffer = np.frombuffer(dense_mat_).reshape(*dense_mat_shape)
        np.copyto(dense_mat_buffer, dense_mat)

        ind_lst = np.array_split(np.arange(dense_mat.shape[0]), NUM_CPU)
        with Pool(processes=NUM_CPU, initializer=init_worker, initargs=(L_data_, L_indices_, L_indptr_,
                                                                        U_data_, U_indices_, U_indptr_,
                                                                        dense_mat_, dense_mat_shape,
                                                                        perm_c_, perm_r_)) as pool:
            result = pool.map(worker_func, ind_lst)
        return np.vstack(result)

    # testing
    @staticmethod
    def sp_solve_dense3d_mp(sps_mat: csc_matrix,
                            dense_mat: np.ndarray, tasks: int = 8) -> np.ndarray:
        """Compute product of the inverse of given sparse matrix and given 3-d dense array in parallel; used for repeated OLS.
        Each batch of SNPs require one LU factorization.

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
        data, indices, indptr = sps_mat.data, sps_mat.indices, sps_mat.indptr
        data_shape = data.shape
        indices_shape = indices.shape
        indptr_shape = indptr.shape

        data_ = RawArray('d', data_shape[0])
        data_buffer = np.frombuffer(data_)
        np.copyto(data_buffer, data)
        indices_ = RawArray(ctypes.c_int32, indices_shape[0])
        indices_buffer = np.frombuffer(indices_, dtype='int32')
        np.copyto(indices_buffer, indices)
        indptr_ = RawArray(ctypes.c_int32, indptr_shape[0])
        indptr_buffer = np.frombuffer(indptr_, dtype='int32')
        np.copyto(indptr_buffer, indptr)

        dense_mat_shape = dense_mat.shape
        dense_mat_ = RawArray(
            'd', dense_mat_shape[0] * dense_mat_shape[1] * dense_mat_shape[2])
        dense_mat_buffer = np.frombuffer(dense_mat_).reshape(*dense_mat_shape)
        np.copyto(dense_mat_buffer, dense_mat)
        shape = sps_mat.shape

        nbytes = int(dense_mat.shape[1] * dense_mat.shape[2] * 8)
        max_nbytes = 2147483647
        if nbytes > max_nbytes:
            raise ValueError('Too many bytes to handle for multiprocessing.')
        oldtasks = tasks
        while int(m / tasks) * nbytes > max_nbytes:
            tasks += NUM_CPU
        logger.info(f'Original number of tasks: {oldtasks}; Final number of tasks: {tasks}')
        ind_lst = np.array_split(np.arange(dense_mat.shape[0]), tasks)
        with Pool(processes=tasks if NUM_CPU > tasks else NUM_CPU, initializer=init_worker_, initargs=(data_, indices_, indptr_, shape,
                                                                         dense_mat_, dense_mat_shape)) as pool:
            result = pool.map(worker_func_, ind_lst)
            logger.info('Got result back.')
        return np.vstack(result)

    # testing
    def fit_snps_eff(self, gts: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
            gts_ = gts.reshape((gts.shape[0], int(k * l)))
            self.logger.info('Calculating projecting matrix...')
            M_X: np.ndarray = gts_ - self.Z.dot(solve(self.Z.T @ self.Z, self.Z.T.dot(gts_)))
            logger.info('Projecting genotype...')
            X_: np.ndarray = M_X.reshape((self.n, k, l)).transpose(2, 0, 1)
            if __debug__:
                M: np.ndarray = - self.Z @ solve(self.Z.T @ self.Z, self.Z.T)
                M.flat[::self.n + 1] += 1
                X__ = M.dot(gts_).reshape((self.n, k, l)).transpose(2, 0, 1)
                np.testing.assert_array_almost_equal(X__, X_)
            self.logger.info('Done projecting genotype.')
            self.logger.info('Projecting phenotype...')
            # y: np.ndarray = M @ self.y
            y: np.ndarray = self.y - self.Z @ solve(self.Z.T @ self.Z, self.Z.T.dot(self.y))
            self.logger.info('Start estimating snp effects...')
            Vinv_X: np.ndarray = self.sp_solve_dense3d_lu(self.V, X_)
            Vinv_y: np.ndarray = self.V_lu.solve(y)
            XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', X_, Vinv_X)
            XT_Vinv_y: np.ndarray = np.einsum('...ij,i', X_, Vinv_y)
        else:
            gts = gts.transpose(2, 0, 1)
            Vinv_X: np.ndarray = self.sp_solve_dense3d_lu(self.V, gts)
            XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', gts, Vinv_X)
            XT_Vinv_y: np.ndarray = np.einsum('...ij,i', gts, self.Vinv_y)
        alpha: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
        alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
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

    # testing
    def Vinv_X(self, gts) -> np.ndarray:
        return self.sp_solve_dense3d(self.V, gts.transpose(2, 0, 1))

    # testing
    def Vinv_X_lu(self, gts) -> np.ndarray:
        return self.sp_solve_dense3d_lu(self.V, gts.transpose(2, 0, 1))

    # testing
    def Vinv_X_lu_mp(self, gts) -> np.ndarray:
        return self.sp_solve_dense3d_lu_mp(self.V, gts.transpose(2, 0, 1))

    # testing
    def Vinv_X_mp(self, gts) -> np.ndarray:
        return self.sp_solve_dense3d_mp(self.V, gts.transpose(2, 0, 1), gts.shape[2])

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
            xx_ = solve(ZT_Vinv_Z, self.Vinv_Z.T) @ self.y
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
            return -0.5 * (self.V_logdet + logdet_eT_Vinv_e + yT_Vinv_y - eT_Vinv_y ** 2 / eT_Vinv_e)

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

        def nll(x):
            self._varcomps = tuple(x)
            return -1 * self.reml_loglik
        res: OptimizeResult = minimize(nll, x0=np.array(self._varcomps), options={'gtol': 1e-8, 'eps': 1e-10},
                                       method='L-BFGS-B', bounds=[(0.0001 * self.y.var(), self.y.var()) for i in range(self.n_varcomps)])
        if res.success:
            self.optimized = True
            self._varcomps = tuple(res.x)
            self.logger.info(f'Finished LBFGS-B. # of Likelihood evaluations: {res.nfev}')
        else:
            raise ValueError(f'Scipy minimize failed: {res.message}')

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
