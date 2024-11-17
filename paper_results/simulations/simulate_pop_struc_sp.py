import numpy as np
from scipy.stats import multivariate_normal
from itertools import combinations_with_replacement, product
from scipy.optimize import fmin_l_bfgs_b, minimize, OptimizeResult
from numpy.linalg import slogdet, solve, inv
from scipy.sparse import csc_matrix, tril
from scipy.sparse.linalg import splu
from collections import defaultdict
import os
import code

num_threads = 1
os.environ["OMP_NUM_THREADS"] = str(num_threads)
# export OPENBLAS_NUM_THREADS=...
os.environ["OPENBLAS_NUM_THREADS"] = str(num_threads)
# export MKL_NUM_THREADS=...
os.environ["MKL_NUM_THREADS"] = str(num_threads)
# export VECLIB_MAXIMUM_THREADS=...
os.environ["VECLIB_MAXIMUM_THREADS"] = str(num_threads)
# export NUMEXPR_NUM_THREADS=...
os.environ["NUMEXPR_NUM_THREADS"] = str(num_threads)

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

def build_sib_arr(fam_labels):   
    """Build lower-triangular nonzero entries of sibship matrix.

    Args:
        fam_labels (FamLabels): ndarray of family id strings corresponding to each individual.

    Returns:
        SparseGRMRepr: sparse GRM representation.
    """
    data = []
    row_ind = []
    col_ind = []

    label_indices = defaultdict(list)
    for l, f in enumerate(fam_labels):
        label_indices[f].append(l)
    f = lambda lst: len(lst) * (len(lst) + 1) / 2
    n = sum(map(f, label_indices.values()))

    for f, indices_lst in label_indices.items():
        for pair in combinations_with_replacement(indices_lst, 2):
            ind1, ind2 = max(pair[0], pair[1]), min(pair[0], pair[1])
            data.append(1)
            row_ind.append(ind1)
            col_ind.append(ind2)
    return np.array(data), np.array(row_ind, dtype='uint32'), np.array(col_ind, dtype='uint32')

def simulate_ind(nsnp,f):
    return np.random.binomial(1,f,nsnp*2).reshape((nsnp,2))

def simulate_sibs(father,mother, blocksize=200):
    # Compute blocks without recombination
    n_blocks = np.int(np.ceil(father.shape[0]/float(blocksize)))
    blocksizes = np.zeros((n_blocks),dtype=int)
    for i in range(n_blocks-1):
        blocksizes[i] = blocksize
    blocksizes[n_blocks-1] = father.shape[0]-np.sum(blocksizes)
    meioses = [[np.zeros((father.shape[0]),dtype=np.int8),np.zeros((father.shape[0]),dtype=np.int8)],  # sib1_f, sib1_m
               [np.zeros((father.shape[0]),dtype=np.int8),np.zeros((father.shape[0]),dtype=np.int8)]]  # sib2_f, sib2_m
    for i in range(2):
        for j in range(2):
            block_start=0
            for k in range(n_blocks):
                block_end = block_start+blocksizes[k]
                meioses[i][j][block_start:block_end] = np.random.binomial(1,0.5)
                block_start = block_end
    ibd = np.array(meioses[0][0]==meioses[1][0],dtype=np.int8)+np.array(meioses[0][1]==meioses[1][1],dtype=np.int8)
    gts = np.zeros((2,father.shape[0],2),dtype=np.int8)
    snp_index = np.array([x for x in range(0, father.shape[0])])
    for i in range(0,2):
        gts[i, :, 0] = father[snp_index,meioses[i][0]]
        gts[i, :, 1] = mother[snp_index, meioses[i][1]]
    return [gts,ibd]

def simulate_sibs_phased(father,mother, blocksize=200):
    # Compute blocks without recombination
    n_blocks = np.int(np.ceil(father.shape[0]/float(blocksize)))
    blocksizes = np.zeros((n_blocks),dtype=int)
    for i in range(n_blocks-1):
        blocksizes[i] = blocksize
    blocksizes[n_blocks-1] = father.shape[0]-np.sum(blocksizes)
    meioses = [[np.zeros((father.shape[0]),dtype=np.int8),np.zeros((father.shape[0]),dtype=np.int8)],  # sib1_f, sib1_m
               [np.zeros((father.shape[0]),dtype=np.int8),np.zeros((father.shape[0]),dtype=np.int8)]]  # sib2_f, sib2_m
    for i in range(2):
        for j in range(2):
            block_start=0
            for k in range(n_blocks):
                block_end = block_start+blocksizes[k]
                meioses[i][j][block_start:block_end] = np.random.binomial(1,0.5) # the whole block passes to sibling j
                block_start = block_end
    # (nfams, 2) ibd: entry in first column is 1 (0) if (not) sharing allele from father, entry in second column is 1 if sharing allele from mother
    # sum of two column give IBD state
    ibd = np.vstack([meioses[0][0]==meioses[1][0], meioses[0][1]==meioses[1][1]],)
    gts = np.zeros((2,father.shape[0],2),dtype=np.int8)
    snp_index = np.array([x for x in range(0, father.shape[0])])
    for i in range(0,2):
        gts[i, :, 0] = father[snp_index,meioses[i][0]] # get father allele from meiosis
        gts[i, :, 1] = mother[snp_index, meioses[i][1]] # get mother allele from meiosis
    return [gts,ibd]

def gt_convert(g):
    if g==0:
        return 'AA'
    elif g==1:
        return 'AG'
    elif g==2:
        return 'GG'

def pedline(ped,gts):
    pline = ped[0]+' '
    for i in range(1,4):
        pline += ped[i]+' '
    pline += '1 -9 '
    for i in range(0,gts.shape[0]-1):
        pline += gt_convert(gts[i])+' '
    pline += gt_convert(gts[i+1])+'\n'
    return pline

def impute_from_sibs(g1, g2, ibd, f):
    if ibd==2:
        return g1+2*f
    elif ibd==0:
        return g1+g2
    elif ibd==1:
        gsum = g1+g2
        if gsum==0:
            return f
        elif gsum==1:
            return 1+f
        elif gsum==2:
            return 1+2*f
        elif gsum==3:
            return 2+f
        elif gsum==4:
            return 3+f

def impute_from_sibs_phased(g1,g2,ibd,f):
    imp = 0.0
    if ibd[0]:
        imp += int(g1[0])+f
    else:
        imp += int(g1[0])+int(g2[0])
    if ibd[1]:
        imp += int(g1[1]) + f
    else:
        imp += int(g1[1]) + int(g2[1])
    return imp

def impute_from_sibs_po_phased(g1,g2,gm, ibd,f):
    # imp = 0.0
    # if ibd[0]:
    #     imp += int(g1[0])+f
    # else:
    #     imp += int(g1[0])+int(g2[0])
    # if ibd[1]:
    #     imp += int(g1[1]) + f
    # else:
    #     imp += int(g1[1]) + int(g2[1])
    # return imp
    if sum(ibd) == 0:
        return sum(g1) + sum(g2) - sum(gm)
    if sum(ibd) == 1:
        shared_sibs = 0 if ibd[0] == 1 else 1 # allele shared between sib pair
        assert g1[shared_sibs] == g2[shared_sibs]
        if shared_sibs == 0:  # shared between sib pair and unobserved parent (father)
            return g1[shared_sibs] + f
        else: # shared between sib pair and observed parent (mother)
            unshared_sibs = 0 if shared_sibs == 1 else 1
            return g1[unshared_sibs] + g2[unshared_sibs]
    if sum(ibd) == 2:
        # take the allele unshared with observed parent (mother), i.e., shared with father
        assert g1[0] == g2[0]
        assert g1[1] == g2[1]
        return g1[0] + f


def loglik(sf, se, Vf, y):
    n = Vf.shape[0]
    V = sf * Vf + se * csc_matrix((np.ones(n), (np.arange(0, n), np.arange(0, n))), shape=(n, n))
    e: np.ndarray = np.ones(n)
    lu = splu(V)
    Vinv_y = lu.solve(y)
    Vinv_e = lu.solve(e)
    yT_Vinv_y: float = y @ Vinv_y
    eT_Vinv_e: float = e @ Vinv_e
    eT_Vinv_y: float = e @ Vinv_y
    logdet_eT_Vinv_e: float = np.log(e @ Vinv_e)
    diag_l: np.ndarray = lu.L.diagonal()
    diag_u: np.ndarray = lu.U.diagonal()
    V_logdet: float = np.log(diag_l).sum() + np.log(diag_u).sum()
    # _, V_logdet = np.linalg.slogdet(V)
    return -0.5 * (V_logdet + logdet_eT_Vinv_e + yT_Vinv_y - eT_Vinv_y ** 2 / eT_Vinv_e)




def optimize(Vf, y):
    yvar = np.var(y)
    sf = se = yvar / 2
    def nll(x):
        sf, se = x
        return -1 * loglik(sf, se, Vf, y)
    res: OptimizeResult = minimize(nll, x0=np.array([sf, se]),
                                       method='L-BFGS-B', bounds=[(0.0001 * yvar, yvar) for i in range(2)])
    if res.success:
        optimized = True
        sf, se = tuple(res.x)
    else:
        raise ValueError(f'Scipy minimize failed: {res.message}')
    return sf, se


def _simple_fit(V, gts, inds, y):
    gts = gts.transpose(2, 0, 1)
    Vinv_X: np.ndarray = sp_solve_dense3d_lu(V[inds, :][:, inds], gts)
    XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', gts, Vinv_X)
    XT_Vinv_y: np.ndarray = Vinv_X.transpose(0, 2, 1).dot(y[inds])
    alpha: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
    alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
    return alpha, alpha_cov


def gwas(V_, X_, X_comple_, y_, n_rel=None, r=0, Fst=0, compare_with_theo=False):
    if n_rel is not None:
        V = V_[:(2*n_rel), :(2*n_rel)].copy()
        X = X_[:(2*n_rel)].copy()
        X_comple = X_comple_[:(2*n_rel)].copy()
        y = y_[:(2*n_rel)].copy()
    else:
        V = V_.copy()
        X = X_.copy()
        X_comple = X_comple_.copy()
        y = y_.copy()
    lu = splu(V)
    y -= y.mean()
    for i in range(0, X.shape[1]):
        X[:, i, :] = X[:, i, :] - np.mean(X[:, i, :], axis=0)
        X_comple[:, i, :] = X_comple[:, i, :] - np.mean(X_comple[:, i, :], axis=0)
    Vinv_X: np.ndarray = sp_solve_dense3d_lu(V, X.transpose(2,0,1))
    XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', X.transpose(2,0,1), Vinv_X)
    XT_Vinv_y: np.ndarray = np.einsum('...ij,i', X.transpose(2,0,1), lu.solve(y))
    alpha: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
    alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
    alpha_ses =  np.sqrt(
        np.diagonal(alpha_cov, axis1=1, axis2=2))

    Vinv_X: np.ndarray = sp_solve_dense3d_lu(V, X_comple.transpose(2,0,1))
    XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', X_comple.transpose(2,0,1), Vinv_X)
    XT_Vinv_y: np.ndarray = np.einsum('...ij,i', X_comple.transpose(2,0,1), lu.solve(y))
    alpha_comple: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
    alpha_cov_comple: np.ndarray = np.linalg.inv(XT_Vinv_X)
    alpha_ses_comple =  np.sqrt(
        np.diagonal(alpha_cov_comple, axis1=1, axis2=2))
    if n_rel is not None or compare_with_theo:
        compute_theoretical = lambda f,n: 6*(1-r**2)/(7+5*r)/(1-Fst)/n/f/(1-f)
        print(Fst, 'young et al. empirical', alpha_cov[:, 0, 0].mean(), compute_theoretical(0.5, 0.5*X_.shape[0]))
    return alpha, alpha_ses, alpha_cov, alpha_comple, alpha_ses_comple, alpha_cov_comple

def sib_diff(X, y, sf, se, Fst=0):
    compute_theoretical = lambda f,n: (1-sf)/(1-Fst)/n/f/(1-f)
    X_diff = X[:, 0,:] - X[:,1,:]
    y_diff = y[:, 0] - y[:, 1]
    X_diff -= np.mean(X_diff, axis=0)
    y_diff -= y_diff.mean()
    XX = np.zeros((X.shape[2],1))
    Xy = np.zeros((X.shape[2],1))
    for i in range(XX.shape[0]):
        XX[i,0] = (X_diff[:,i]**2).sum()
        Xy[i,0] = (X_diff[:,i]*y_diff[:]).sum()
    print(XX.shape, Xy.shape)
    alpha: np.ndarray = Xy / XX
    alpha_cov: np.ndarray = (2-2*sf)/XX
    alpha_ses =  np.sqrt(alpha_cov)
    print(Fst, 'sib difference empirical var', alpha_cov.mean(), compute_theoretical(0.5, X.shape[0]))
    return alpha, alpha_ses,alpha_cov,

def get_VC(X_, y_):
    X = X_[:,0].copy()
    y = y_.copy()
    y -= y.mean()
    X -= np.mean(X, axis=0)
    XT_Vinv_X: np.ndarray = np.einsum('...i,...i', X.T, X.T)
    XT_Vinv_y: np.ndarray = np.einsum('...i,i', X.T, y)
    print(XT_Vinv_X.shape, XT_Vinv_y.shape)
    alpha: np.ndarray = XT_Vinv_y / XT_Vinv_X
    print(alpha.shape, X.shape)
    # res = np.einsum('...ij', X.T, alpha) - y
    res = alpha[0] * X[:, 0] - y
    print(1 - np.sum(res ** 2) / np.sum(y ** 2))
    print(alpha.shape, res.shape)
    exit()
    alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
    alpha_ses =  np.sqrt(alpha_cov)

    return alpha, alpha_ses, alpha_cov,


def standard_gwas(V_, X_, y_, n_rel=None):
    if n_rel is not None:
        V = V_[:(2*n_rel), :(2*n_rel)].copy()
        X = X_[:(2*n_rel),0:1].copy()
        y = y_[:(2*n_rel)].copy()
    else:
        V = V_.copy()
        X = X_[:,0:1].copy()
        y = y_.copy()
    lu = splu(V)
    y -= y.mean()
    for i in range(0, X.shape[1]):
        X[:, i, :] = X[:, i, :] - np.mean(X[:, i, :], axis=0)
    Vinv_X: np.ndarray = sp_solve_dense3d_lu(V, X.transpose(2,0,1))
    XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', X.transpose(2,0,1), Vinv_X)
    XT_Vinv_y: np.ndarray = np.einsum('...ij,i', X.transpose(2,0,1), lu.solve(y))
    alpha: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
    alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
    alpha_ses =  np.sqrt(alpha_cov)

    return alpha, alpha_ses, alpha_cov,


def standard_gwas_with_no_V(X_, y_, n_rel=None):
    if n_rel is not None:
        X = X_[:(2*n_rel),0:1].copy()
        y = y_[:(2*n_rel)].copy()
    else:
        X = X_[:,0:1].copy()
        y = y_.copy()
    y -= y.mean()
    for i in range(0, X.shape[1]):
        X[:, i, :] = X[:, i, :] - np.mean(X[:, i, :], axis=0)
    # Vinv_X: np.ndarray = sp_solve_dense3d_lu(V, X.transpose(2,0,1))
    XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', X.transpose(2,0,1), X.transpose(2,0,1))
    XT_Vinv_y: np.ndarray = np.einsum('...ij,i', X.transpose(2,0,1), y)
    alpha: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
    alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
    alpha_ses =  np.sqrt(alpha_cov)

    return alpha, alpha_ses, alpha_cov,

def robust_est(ibd_, V_, X_, X_comple_, y_, n_rel=None, separate=False,sf=0.5, se=0.5, npop=2, nfam=20000, Fst=0.1):
    if n_rel is not None:
        ibd = ibd_[:(2*n_rel),:].copy()
        V = V_[:(2*n_rel), :(2*n_rel)].copy()
        X = X_[:(2*n_rel)].copy()
        X_comple = X_comple_[:(2*n_rel)].copy()
        y = y_[:(2*n_rel)].copy()
    else:
        ibd = ibd_.copy()
        V = V_.copy()
        X = X_.copy()
        X_comple = X_comple_.copy()
        y = y_.copy()
    y -= y.mean()
    alpha = []
    alpha_complete = []
    alpha_ses = []
    alpha_ses_complete = []
    alpha_cov = []
    alpha_cov_complete = []
    if separate:
        alpha0 = []
        alpha0_complete = []
        alpha0_ses = []
        alpha0_ses_complete = []
        alpha0_cov = []
        alpha0_cov_complete = []
        alpha1 = []
        alpha1_complete = []
        alpha1_ses = []
        alpha1_ses_complete = []
        alpha1_cov = []
        alpha1_cov_complete = []
        bias0 = []
        bias1 = []
        bias1_2ndrow = []
        bias = []
        theoretical_direct_var0 = []
        theoretical_direct_var1 = []
        compute_theoretical0 = lambda f,n: (1-sf)/2/(1-Fst)/n/f/(1-f)
        compute_theoretical1 = lambda f,n: 3*(1-sf**2)/2/(2+sf)/(1-Fst)/n/f/(1-f)
    theoretical_direct_var = []
    compute_theoretical = lambda f,n: 6*(1-sf**2)/(7+5*sf)/(1-Fst)/n/f/(1-f)
    for s in range(X.shape[2]):
        if np.isnan(X[:, :, s]).sum() > 0:
            raise ValueError('Missing entry.')
        ibd0_inds = ibd[:, s] == 0
        ibd1_inds = ibd[:, s] == 1
        gts0 = X[ibd0_inds, :, s:s+1]
        gts0_complete = X_comple[ibd0_inds, :, s:s+1]
        gts1 = X[ibd1_inds, :, s:s+1]
        gts1_complete = X_comple[ibd1_inds, :, s:s+1]
        if np.mean(gts0[:, 0, :], axis=0)/2 > 0.99 or np.mean(gts0[:, 0, :], axis=0)/2 < 0.01 or \
            np.mean(gts1[:, 0, :], axis=0)/2 < 0.01 or np.mean(gts1[:, 0, :], axis=0)/2 > 0.99:
            print(s, 'low maf', np.mean(gts0[:, 0, :], axis=0)/2, np.mean(gts1[:, 0, :], axis=0)/2)
            continue
        for d in range(0, gts0.shape[1]):
            gts0[:, d, :] = gts0[:, d, :] - np.mean(gts0[:, d, :], axis=0)
            gts0_complete[:, d, :] = gts0_complete[:, d, :] - np.mean(gts0_complete[:, d, :], axis=0)
        try:
            alpha0_, alpha0_cov_, = _simple_fit(V, gts0, ibd0_inds, y)
            alpha0_complete_, alpha0_cov_complete_, = _simple_fit(V, gts0_complete, ibd0_inds, y)
        except np.linalg.LinAlgError:
            code.interact(local=locals())
            raise ValueError(str(s) + ' linalg')
            continue
        for d in range(0, gts1.shape[1]):
            gts1[:, d, :] = gts1[:, d, :] - np.mean(gts1[:, d, :], axis=0)
            gts1_complete[:, d, :] = gts1_complete[:, d, :] - np.mean(gts1_complete[:, d, :], axis=0)
        
        try:
            alpha1_, alpha1_cov_, = _simple_fit(V, gts1, ibd1_inds, y)
            alpha1_complete_, alpha1_cov_complete_, = _simple_fit(V, gts1_complete, ibd1_inds, y)
        except np.linalg.LinAlgError:
            code.interact(local=locals())
            raise ValueError(str(s) + ' linalg')
            continue
        if separate:
            lu = splu(V[ibd0_inds, :][:, ibd0_inds])
            bias0.append(np.linalg.inv(gts0.squeeze().T.dot(lu.solve(gts0.squeeze()))).dot(gts0.squeeze().T.dot(lu.solve(gts0_complete.squeeze())))[0,:])
            lu = splu(V[ibd1_inds, :][:, ibd1_inds])
            bias1.append(np.linalg.inv(gts1.squeeze().T.dot(lu.solve(gts1.squeeze()))).dot(gts1.squeeze().T.dot(lu.solve(gts1_complete.squeeze())))[0,:])
            bias1_2ndrow.append(np.linalg.inv(gts1.squeeze().T.dot(gts1.squeeze())).dot(gts1.squeeze().T.dot(gts1_complete.squeeze()))[1,:])

        robust_var = (alpha0_cov_[0, 0, 0]**(-1) + alpha1_cov_[0, 0, 0]**(-1))**(-1)
        robust_var_complete = (alpha0_cov_complete_[0,0,0]**(-1) + alpha1_cov_complete_[0,0,0]**(-1))**(-1)
        robust_alpha = robust_var * (alpha0_cov_[0,0,0]**(-1) * alpha0_[0,0] + alpha1_cov_[0,0,0]**(-1) * alpha1_[0,0])
        robust_alpha_complete = robust_var_complete * (alpha0_cov_complete_[0,0,0]**(-1) * alpha0_complete_[0,0] + alpha1_cov_complete_[0,0,0]**(-1) * alpha1_complete_[0,0])
        alpha.append(robust_alpha)
        alpha_complete.append(robust_alpha_complete)
        alpha_cov.append(robust_var)
        alpha_cov_complete.append(robust_var_complete)
        alpha_ses.append(np.sqrt(robust_var))
        alpha_ses_complete.append(np.sqrt(robust_var_complete))
        if separate:
            alpha0.append(alpha0_)
            alpha0_complete.append(alpha0_complete_)
            alpha0_cov.append(alpha0_cov_.squeeze())
            alpha0_cov_complete.append(alpha0_cov_complete_.squeeze())
            alpha0_ses.append(np.sqrt(np.diag(alpha0_cov_.squeeze())))
            alpha0_ses_complete.append(np.sqrt(np.diag(alpha0_cov_complete_.squeeze())))

            alpha1.append(alpha1_)
            alpha1_complete.append(alpha1_complete_)
            alpha1_cov.append(alpha1_cov_.squeeze())
            alpha1_cov_complete.append(alpha1_cov_complete_.squeeze())
            alpha1_ses.append(np.sqrt(np.diag(alpha1_cov_.squeeze())))
            alpha1_ses_complete.append(np.sqrt(np.diag(alpha1_cov_complete_.squeeze())))
            theoretical_direct_var0.append(compute_theoretical0(X[:, 0, s].mean()/2, gts0.shape[0]/2))
            theoretical_direct_var1.append(compute_theoretical1(X[:, 0, s].mean()/2, gts1.shape[0]/2))
        theoretical_direct_var.append(compute_theoretical(X[:, 0, s].mean()/2, X.shape[0]/2))

    alpha = np.vstack(alpha)
    alpha_complete = np.vstack(alpha_complete)
    alpha_ses = np.vstack(alpha_ses)
    alpha_ses_complete = np.vstack(alpha_ses_complete)
    alpha_cov = np.stack(alpha_cov, axis=0)
    alpha_cov_complete = np.stack(alpha_cov_complete, axis=0)
    if separate:
        alpha0 = np.vstack(alpha0)
        alpha0_complete = np.vstack(alpha0_complete)
        alpha0_ses = np.vstack(alpha0_ses)
        alpha0_ses_complete = np.vstack(alpha0_ses_complete)
        alpha0_cov = np.stack(alpha0_cov, axis=0)
        alpha0_cov_complete = np.stack(alpha0_cov_complete, axis=0)

        alpha1 = np.vstack(alpha1)
        alpha1_complete = np.vstack(alpha1_complete)
        alpha1_ses = np.vstack(alpha1_ses)
        alpha1_ses_complete = np.vstack(alpha1_ses_complete)
        alpha1_cov = np.stack(alpha1_cov, axis=0)
        alpha1_cov_complete = np.stack(alpha1_cov_complete, axis=0)
        bias0 = np.vstack(bias0)
        bias1 = np.vstack(bias1)
        bias1_2ndrow = np.vstack(bias1_2ndrow)
        # bias = np.vstack(bias)
        def true_bias1(Fst, r):
            return 3*(Fst*(2-3*r)-r+2)/(1-Fst)/(2+r) - 4*(2*Fst+1)*(1-r)/(1-Fst)/(r+2)
        def true_bias1_2ndrow(Fst, r):
            return 3*(3*Fst+1)*(Fst*(2-3*r)-r+2)/(1-Fst)/(2+r)/(2*Fst+1) - 4*(3*Fst+1)*(1-r)/(1-Fst)/(r+2)
        print(np.mean(bias0, axis=0), np.var(bias0, axis=0), )
        print(np.mean(bias1, axis=0), np.var(bias1, axis=0), true_bias1(Fst, 0.5) )
        print(np.mean(bias1_2ndrow, axis=0), np.var(bias1_2ndrow, axis=0), true_bias1_2ndrow(Fst, 0.5))
        # print(np.mean(bias, axis=0), np.var(bias, axis=0), )

        print(alpha0_cov[:,0,0].mean(), alpha0_cov[:,0,0].var(), np.mean(theoretical_direct_var0))
        print(alpha1_cov[:,0,0].mean(), alpha1_cov[:,0,0].var(), np.mean(theoretical_direct_var1))
    print(alpha_cov.mean(), alpha_cov.var(), np.mean(theoretical_direct_var))
    
    
    
    if separate:
        return alpha, alpha_ses, alpha_cov, alpha_complete, alpha_ses_complete, alpha_cov_complete, \
            alpha0, alpha0_ses, alpha0_cov, alpha0_complete, alpha0_ses_complete, alpha0_cov_complete, \
            alpha1, alpha1_ses, alpha1_cov, alpha1_complete, alpha1_ses_complete, alpha1_cov_complete
    return alpha, alpha_ses, alpha_cov, alpha_complete, alpha_ses_complete, alpha_cov_complete

def savez(*args, **kwds):
    if save:
        print('SAVING')
        np.savez(*args, **kwds)
    else:
        print('NOT SAVING')


def main(nsnp=1500, nfam=1500, v=1, phased=True, blocksize=1, f=0.5, h2quad=0.5, h2_pop=0.1, unrel=False, npop=4, Fst=0.01, balding_nicols=True, save=True, diff=0.05, robust=False):
    gens = 1 
    fsize = 2



    direct_var = 1/100 #1/nsnp
    indirect_var = 0#0.25/nsnp
    direct_indirect_corr = 0.
    direct_indirect_cov = direct_indirect_corr*np.sqrt(indirect_var*direct_var)
    direct_indirect_cov_matrix = [[direct_var, direct_indirect_cov],[direct_indirect_cov, indirect_var]]
    direct_indirect = multivariate_normal.rvs(cov = direct_indirect_cov_matrix, size=nsnp)
    if nsnp == 1:
        direct_indirect = direct_indirect.reshape((-1, 2))
    direct_effects = direct_indirect[:,0].reshape((-1, 1))
    direct_effects[np.random.choice([i for i in range(nsnp)], size=int(nsnp*0.9), replace=True), 0] = 0.
    indirect_effects = direct_indirect[:,1].reshape((-1, 1))

    sib_gts_lst = []
    sib_phen_lst = []
    ibd_lst = []
    father_gts_lst = []
    mother_gts_lst = []

    if Fst == 0 and balding_nicols:
        pop_fs = np.ones((npop, nsnp), dtype=float) * 0.5
    else:
        if balding_nicols:
            a = (1-Fst)/Fst*0.5
            b = (1-Fst)/Fst*0.5
        pop_fs = np.zeros((npop, nsnp), dtype=float)
        for snp in range(nsnp):
            if nsnp == 1 and npop == 2:
                pop_fs[0, snp] = 0.5 + np.sqrt(Fst) / 2
                pop_fs[1, snp] = 0.5 - np.sqrt(Fst) / 2
            else:
                if balding_nicols:
                    pop_fs[:, snp] = np.random.beta(a, b, size=npop)
                else:
                    pop_fs[0, snp] = 0.5 + np.sqrt(Fst) / 2
                    pop_fs[1, snp] = 0.5 - np.sqrt(Fst) / 2
              

    pop_effs = np.random.randn(npop)

    pop_indiv = np.hstack([[e for i in range(nfam * 2)] for e in pop_effs])
    scale_fctr = np.sqrt(h2_pop/np.var(pop_indiv))
    pop_effs = pop_effs * scale_fctr
    for s, pf in enumerate(pop_fs):
        if phased: ibd = np.zeros((nfam, 2, nsnp),dtype=np.int8)
        else: ibd = np.zeros((nfam,fsize*(fsize-1)//2,nsnp),dtype=np.int8)
        father_gts = np.full((nfam,nsnp,2), fill_value=np.nan, dtype=int)
        mother_gts = np.full((nfam,nsnp,2), fill_value=np.nan, dtype=int)
        for i,f in enumerate(pf):
            father_gts[:,i,:] = np.random.binomial(1, f, size=(nfam,2))
            mother_gts[:,i,:] = np.random.binomial(1, f, size=(nfam,2))
        sib_gts = np.zeros((nfam,fsize,nsnp,2))
        # second fpgs can be replaced with normal joint regression 
        father_sum = np.sum(father_gts, axis=2)
        father_indirect = np.zeros(nfam)
        father_direct = (father_sum@direct_effects).flatten()
        father_noise = np.random.normal(size=father_direct.shape)
        father_phen = father_direct+father_indirect+father_noise

        mother_sum = np.sum(mother_gts, axis=2)
        mother_indirect = np.zeros(nfam)
        mother_direct = (mother_sum@direct_effects).flatten()
        mother_noise = np.random.normal(size=mother_direct.shape)
        mother_phen = mother_direct+mother_indirect


        father_indexes = np.argsort(father_phen) # + np.random.normal(scale=am_noise*np.std(father_phen), size=father_phen.shape))
        father_direct = father_direct[father_indexes]
        father_indirect = father_indirect[father_indexes]
        father_phen = father_phen[father_indexes]
        father_gts = father_gts[father_indexes,:,:]
        father_sum = father_sum[father_indexes,:]

        mother_indexes = np.argsort(mother_phen) # + np.random.normal(scale=am_noise*np.std(mother_phen), size=mother_phen.shape))
        mother_direct = mother_direct[mother_indexes]
        mother_indirect = mother_indirect[mother_indexes]
        mother_phen = mother_phen[mother_indexes]
        mother_gts = mother_gts[mother_indexes,:,:]
        mother_sum = mother_sum[mother_indexes,:]
        for i in range(0,nfam):
            if i % 2000 == 0:
                print(i)
            father = father_gts[i, :, :]
            mother = mother_gts[i, :, :]
            if phased:
                sibs = simulate_sibs_phased(father,mother, blocksize = blocksize)
            else:
                sibs = simulate_sibs(father,mother, blocksize = blocksize)
            sib_gts[i, :, :, :] = sibs[0]
            ibd[i,:,:] = sibs[1]
        # tmp=((father_sum+mother_sum)@indirect_effects).flatten()

        # sib_indirect = np.tile(tmp, (2,1)).T
        sib_direct = np.zeros((nfam, fsize))
        sib_phen = np.zeros((nfam, fsize))
        for i in range(fsize):
            sib_direct[:, i] = (np.sum(sib_gts[:, i, :, :], axis=2)@direct_effects).flatten()
            sib_noise = np.random.normal(0, v, size=sib_direct[:, i].shape)
            sib_phen[:, i] = sib_direct[:, i] # +sib_indirect[:, i]
            a_factor = np.sqrt(h2quad)*np.power(np.std(sib_phen[:, i],axis=0),-1)
  
            sib_phen[:, i] = sib_phen[:, i]*a_factor+pop_effs[s]+np.sqrt(1-h2quad-h2_pop)*sib_noise
        sib_gts_lst.append(sib_gts)
        sib_phen_lst.append(sib_phen)
        ibd_lst.append(ibd)
        father_gts_lst.append(father_gts)
        mother_gts_lst.append(mother_gts)
        

    
    sib_gts = np.vstack(sib_gts_lst)
    sib_phen = np.vstack(sib_phen_lst)
    ibd = np.vstack(ibd_lst)
    father_gts = np.vstack(father_gts_lst)
    mother_gts = np.vstack(mother_gts_lst)
    pop_ids = np.hstack([np.ones(nfam, dtype=int) * i for i in range(npop)])
    
    nfam = nfam * npop
    del father_gts_lst, mother_gts_lst, sib_gts_lst, ibd_lst, sib_phen_lst
    import gc
    gc.collect()


    ### imputation
    ### snipar imputation on first half of fams, linear imputation on second half of fams

    # if not phased, get sum of alleles
    if not phased:
        sib_gts = sib_gts.sum(axis=3)
        father_gts = father_gts.sum(axis=2)
        mother_gts = mother_gts.sum(axis=2)
        
    # vec_impute_from_sibs_phased = np.vectorize(impute_from_sibs_phased)
    imp = np.zeros((sib_gts.shape[0],sib_gts.shape[2]),dtype=np.float_)
    fam_ids = np.array([str(i) for i in range(nfam)])
    unrel_inds = np.sort(np.random.choice([i for i in range(nfam)], int(nfam/10*9), replace=False))
    rel_inds = np.setdiff1d([i for i in range(nfam)], unrel_inds)

    freqs = np.mean(sib_gts[:, 0, :, 0], axis=0)
    if not unrel:
        for i in range(nfam): #rel_inds:
            for j in range(sib_gts.shape[2]):
                if phased:
                    imp[i, j] = impute_from_sibs_phased(sib_gts[i, 0, j, :], sib_gts[i, 1, j, :], ibd[i, :, j], freqs[j])
                else:
                    imp[i, j] = impute_from_sibs(sib_gts[i, 0, j], sib_gts[i, 1, j], ibd[i, 0, j], freqs[j])
        # for i in unrel_inds:
        #     for j in range(sib_gts.shape[2]):
        #         if phased:
        #             imp[i, j] = sib_gts[i, 0, j].sum() + 2 * freqs[j]
        #         else:
        #             imp[i, j] = sib_gts[i, 0, j] + 2 * freqs[j]
    else:
        for i in rel_inds:
            for j in range(sib_gts.shape[2]):
                if phased:
                    imp[i, j] = impute_from_sibs_phased(sib_gts[i, 0, j, :], sib_gts[i, 1, j, :], ibd[i, :, j], freqs[j])
                else:
                    imp[i, j] = impute_from_sibs(sib_gts[i, 0, j], sib_gts[i, 1, j], ibd[i, 0, j], freqs[j])
        for i in unrel_inds:
            for j in range(sib_gts.shape[2]):
                if phased:
                    imp[i, j] = sib_gts[i, 0, j].sum() + 2 * freqs[j]
                else:
                    imp[i, j] = sib_gts[i, 0, j] + 2 * freqs[j]

    if phased:
        sib_gts = sib_gts.sum(axis=3)
        father_gts = father_gts.sum(axis=2)
        mother_gts = mother_gts.sum(axis=2)
    ibd = ibd.sum(axis=1)

    if not unrel:
        fam_labels = np.hstack([fam_ids, fam_ids])
        ibd = np.vstack([ibd, ibd])
        X = np.full((sib_gts.shape[0]*2, 2, sib_gts.shape[2]), fill_value=np.nan)
        X_comple = np.full((sib_gts.shape[0]*2, 2, sib_gts.shape[2]), fill_value=np.nan)
        y = np.full(sib_gts.shape[0]*2, np.nan)
        X[:nfam, 0, :] = sib_gts[:, 0, :]
        X[nfam:, 0, :] = sib_gts[:, 1, :]
        X_comple[:nfam, 0, :] = sib_gts[:, 0, :]
        X_comple[nfam:, 0, :] = sib_gts[:, 1, :]
        X[:nfam, 1, :] = imp
        X[nfam:, 1, :] = imp
        X_comple[:nfam, 1, :] = father_gts + mother_gts
        X_comple[nfam:, 1, :] = father_gts + mother_gts
        y[:nfam] = sib_phen[:, 0]
        y[nfam:] = sib_phen[:, 1]

        pop_effs = np.hstack([[e for i in range(int(nfam / len(pop_effs)))] for e in pop_effs])
        pop_effs = np.hstack([pop_effs, pop_effs])

        
        np.testing.assert_equal(0, np.isnan(X).sum())
        np.testing.assert_equal(0, np.isnan(X_comple).sum())
        np.testing.assert_equal(0, np.isnan(y).sum())
    else:
        fam_labels = np.empty(2 * len(rel_inds) + len(unrel_inds), dtype=fam_ids.dtype)
        fam_labels[:(2*len(rel_inds)):2] = fam_ids[rel_inds]
        fam_labels[1:(2*len(rel_inds)):2] = fam_ids[rel_inds]
        fam_labels[(2*len(rel_inds)):] = fam_ids[unrel_inds]
        ibd_ = np.empty((2 * len(rel_inds) + len(unrel_inds), ibd.shape[1]), dtype=ibd.dtype)
        ibd_[:(2*len(rel_inds)):2,:] = ibd[rel_inds,:]
        ibd_[1:(2*len(rel_inds)):2,:] = ibd[rel_inds,:]
        ibd_[(2*len(rel_inds)):,:] = ibd[unrel_inds,:]
        ibd = ibd_
        X = np.full((2 * len(rel_inds) + len(unrel_inds), 2, sib_gts.shape[2]), fill_value=np.nan)
        X_comple = np.full((2 * len(rel_inds) + len(unrel_inds), 2, sib_gts.shape[2]), fill_value=np.nan)
        y = np.full(2 * len(rel_inds) + len(unrel_inds), np.nan)
        pop_effs_ = np.hstack([[e for i in range(int(nfam / len(pop_effs)))] for e in pop_effs])
        pop_effs = np.full(2 * len(rel_inds) + len(unrel_inds), np.nan)

        
        X[:(2*len(rel_inds)):2, 0, :] = sib_gts[rel_inds, 0, :]
        X[1:(2*len(rel_inds)):2, 0, :] = sib_gts[rel_inds, 1, :]
        X[(2*len(rel_inds)):, 0, :] = sib_gts[unrel_inds, 0, :]
        X_comple[:(2*len(rel_inds)):2, 0, :] = sib_gts[rel_inds, 0, :]
        X_comple[1:(2*len(rel_inds)):2, 0, :] = sib_gts[rel_inds, 1, :]
        X_comple[(2*len(rel_inds)):, 0, :] = sib_gts[unrel_inds, 0, :]

        X[:(2*len(rel_inds)):2, 1, :] = imp[rel_inds, :]
        X[1:(2*len(rel_inds)):2, 1, :] = imp[rel_inds, :]
        X[(2*len(rel_inds)):, 1, :] = imp[unrel_inds, :]
        X_comple[:(2*len(rel_inds)):2, 1, :] = father_gts[rel_inds, :] + mother_gts[rel_inds, :]
        X_comple[1:(2*len(rel_inds)):2, 1, :] = father_gts[rel_inds, :] + mother_gts[rel_inds, :]
        X_comple[(2*len(rel_inds)):, 1, :] = father_gts[unrel_inds, :] + mother_gts[unrel_inds, :]
        y[:(2*len(rel_inds)):2] = sib_phen[rel_inds, 0]
        y[1:(2*len(rel_inds)):2] = sib_phen[rel_inds, 1]
        y[(2*len(rel_inds)):] = sib_phen[unrel_inds, 0]
        pop_effs[:(2*len(rel_inds)):2] = pop_effs_[rel_inds]
        pop_effs[1:(2*len(rel_inds)):2] = pop_effs_[rel_inds]
        pop_effs[(2*len(rel_inds)):] = pop_effs_[unrel_inds]

        np.testing.assert_equal(0, np.isnan(X).sum())
        np.testing.assert_equal(0, np.isnan(X_comple).sum())
        np.testing.assert_equal(0, np.isnan(y).sum())

 
    sib_data, sib_row_ind, sib_col_ind = build_sib_arr(fam_labels)
    tril_mat: csc_matrix = csc_matrix((sib_data, (sib_row_ind, sib_col_ind)), shape=(y.shape[0], y.shape[0]))
    triu_mat: csc_matrix = tril(tril_mat, k=-1, format='csc')
    V=tril_mat + triu_mat.T
    sf,se = optimize(V, y - y.mean())
    sf = h2_pop
    se = 1 - sf
    V = sf * V + se * csc_matrix((np.ones(y.shape[0]), (np.arange(0, y.shape[0]), np.arange(0, y.shape[0]))), shape=(y.shape[0], y.shape[0]))




    

    


    if unrel:
        alpha, alpha_ses, alpha_cov, = standard_gwas(V, X, y, len(rel_inds))
        np.savez(f'sp_pairs/Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}_standard_no_unrel', alpha=alpha, alpha_cov=alpha_cov)
        alpha, alpha_ses, alpha_cov, = standard_gwas(V, X, y)
        np.savez(f'sp_pairs/Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}_standard_with_unrel', alpha=alpha, alpha_cov=alpha_cov)
        alpha, alpha_ses, alpha_cov, alpha_comple, alpha_ses_comple, alpha_cov_comple = gwas(V, X, X_comple, y, len(rel_inds), r=sf)
        np.savez(f'sp_pairs/Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}_no_unrel', alpha=alpha, alpha_cov=alpha_cov, alpha_comple=alpha_comple, alpha_cov_comple=alpha_cov_comple)
        alpha1, alpha_ses1, alpha_cov1 = sib_diff(sib_gts[rel_inds], sib_phen[rel_inds], sf, se, Fst=Fst)
        np.savez(f'sp_pairs/Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}_sibdiff', alpha=alpha1, alpha_cov=alpha_cov1)
        alpha, alpha_ses, alpha_cov, alpha_comple, alpha_ses_comple, alpha_cov_comple = gwas(V, X, X_comple, y)
        np.savez(f'sp_pairs/Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}_with_unrel', alpha=alpha, alpha_cov=alpha_cov, alpha_comple=alpha_comple, alpha_cov_comple=alpha_cov_comple)
        alpha, alpha_ses, alpha_cov, alpha_comple, alpha_ses_comple, alpha_cov_comple = robust_est(ibd, V, X, X_comple, y, len(rel_inds))
        np.savez(f'sp_pairs/Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}_robust', alpha=alpha, alpha_cov=alpha_cov, alpha_comple=alpha_comple, alpha_cov_comple=alpha_cov_comple)
    else:
        alpha, alpha_ses, alpha_cov, = standard_gwas(V, X, y)
        if Fst == 0:
            tmp = 'tmp/'
        else:
            tmp = ''
        savez(f'sp_pairs/allsibs/{tmp}Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}_standard', alpha=alpha, alpha_cov=alpha_cov)
        alpha, alpha_ses, alpha_cov, alpha_comple, alpha_ses_comple, alpha_cov_comple = gwas(V, X, X_comple, y, r=sf, Fst=Fst)
        savez(f'sp_pairs/allsibs/{tmp}Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}', alpha=alpha, alpha_cov=alpha_cov, alpha_comple=alpha_comple, alpha_cov_comple=alpha_cov_comple)
        alpha, alpha_ses, alpha_cov, alpha_comple, alpha_ses_comple, alpha_cov_comple = robust_est(ibd, V, X, X_comple, y)
        savez(f'sp_pairs/allsibs/{tmp}Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}_robust', alpha=alpha, alpha_cov=alpha_cov, alpha_comple=alpha_comple, alpha_cov_comple=alpha_cov_comple)
        alpha1, alpha_ses1, alpha_cov1 = sib_diff(sib_gts, sib_phen, sf, se, Fst=Fst)
        savez(f'sp_pairs/allsibs/{tmp}Fst_{Fst}_h2pop_{h2_pop}_h2quad_{h2quad}_sibdiff', alpha=alpha1, alpha_cov=alpha_cov1)







save = True

main(nsnp=20_000, nfam=20_000, h2quad=0.0, h2_pop=0.5, unrel=False, Fst=0, npop=2, save=save, balding_nicols=True)
main(nsnp=20_000, nfam=20_000, h2quad=0.0, h2_pop=0.5, unrel=False, Fst=0.001, npop=2, save=save, balding_nicols=True)
main(nsnp=20_000, nfam=20_000, h2quad=0.0, h2_pop=0.5, unrel=False, Fst=0.01, npop=2, save=save, balding_nicols=True)
main(nsnp=20_000, nfam=20_000, h2quad=0.0, h2_pop=0.5, unrel=False, Fst=0.1, npop=2, save=save, balding_nicols=True)










