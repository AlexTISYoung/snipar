
import numpy as np
from scipy.stats import multivariate_normal
from scipy.linalg import block_diag
from itertools import combinations_with_replacement, product
from scipy.optimize import fmin_l_bfgs_b, minimize, OptimizeResult
from numpy.linalg import slogdet, solve, inv, norm
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse import csc_matrix, tril
from scipy.sparse.linalg import splu, SuperLU, spsolve, cg
from collections import defaultdict
import scipy
import os


save = False
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

def build_sib_arr(fam_labels, x=1, y=1):   
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
            if ind1 == ind2:
                data.append(1)
            elif '_' in f:
                data.append(y)
            else:
                data.append(x)
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





def impute_cond_gau_gwas(sib_gts, f, true_par=None): # true_father, true_mother, true_grandpar_f=None, true_grandpar_m=None):
    R12 = scipy.linalg.block_diag(*[np.array([[1/2,1/8],[1/8,1/2]]) for i in range(int(sib_gts.shape[0]/2))])
    R22 = scipy.linalg.block_diag(*[np.array([[1,1/8],[1/8,1]]) for i in range(int(sib_gts.shape[0]/2))])
    # R12 = 2*f*(1-f)*R12
    # R22 = 2*f*(1-f)*R22
    factor = R12.T @ inv(R22)
    a_minus_mu = sib_gts.T - 2*np.expand_dims(f, axis=-1)
    imp = np.einsum('jk,...k', factor, a_minus_mu) + 2*np.expand_dims(f, axis=-1)
    imp[imp > 2] = 2
    imp[imp < 0] = 0
    linimp = sib_gts + 2*f
    father_imp = imp[:, ::2]
    mother_imp = imp[:, 1::2]

    return 2*father_imp, 2*mother_imp
    

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


from scipy.optimize import minimize, OptimizeResult
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




def savez(*args, **kwds):
    if save:
        print('SAVING')
        np.savez(*args, **kwds)
    else:
        print('NOT SAVING')



def main(nsnp=1500, nfam=1500, v=1, phased=True, blocksize=1, f=0.5, h2quad=0.5, h2_pop=0.1, unrel=False, npop=4, Fst=0.01, balding_nicols=True, save=True, diff=0.05, robust=False):

    fsize = 2

    direct_var = 1/100 #1/nsnp
    indirect_var = 1/100#0.25/nsnp
    direct_indirect_corr = 0.
    direct_indirect_cov = direct_indirect_corr*np.sqrt(indirect_var*direct_var)
    direct_indirect_cov_matrix = [[direct_var, direct_indirect_cov],[direct_indirect_cov, indirect_var]]
    direct_indirect = multivariate_normal.rvs(cov = direct_indirect_cov_matrix, size=nsnp)
    if nsnp == 1:
        direct_indirect = direct_indirect.reshape((-1, 2))
    direct_effects = direct_indirect[:,0].reshape((-1, 1))
    direct_effects[np.random.choice([i for i in range(nsnp)], size=int(nsnp*0.9), replace=True), 0] = 0.
    direct_effects = [1]


    

    cousin_gts_lst = []
    cousin_phen_lst = []
    father_gts_lst = []


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
        for i in range(0,nfam):

            father = father_gts[i, :, :]
            mother = mother_gts[i, :, :]
            if phased:
                sibs = simulate_sibs_phased(father,mother, blocksize = blocksize)
            else:
                sibs = simulate_sibs(father,mother, blocksize = blocksize)
            sib_gts[i, :, :, :] = sibs[0]
            # ibd[i,:,:] = sibs[1]
        


        father_i = sum([[0+3*i, 2+3*i] for i in range(int(sib_gts.shape[0] / 3))],[])
        mother_i = sum([[1+3*i, 0+3*i] for i in range(int(sib_gts.shape[0] / 3))],[])
        

        father_gts = sib_gts[father_i, 0, :, :].copy()
        mother_gts = sib_gts[mother_i, 1, :, :].copy()
        sib_gts = np.zeros((father_gts.shape[0],fsize,nsnp,2))
        for i in range(0,father_gts.shape[0]):
            father = father_gts[i, :, :]
            mother = mother_gts[i, :, :]
            if phased:
                sibs = simulate_sibs_phased(father,mother, blocksize = blocksize)
            else:
                sibs = simulate_sibs(father,mother, blocksize = blocksize)
            sib_gts[i, :, :, :] = sibs[0]
            # ibd[i,:,:] = sibs[1]
        # tmp=((father_sum+mother_sum)@indirect_effects).flatten()

        # n = sib_gts.shape[0]
        cousin_gts = sib_gts[:, 0, :, :].sum(axis=2)
       
        cousin_noise = np.random.normal(0, v, size=(cousin_gts.shape[0],))
        cousin_phen = pop_effs[s] + np.sqrt(1-h2_pop)* cousin_noise
        cousin_gts_lst.append(cousin_gts)
        cousin_phen_lst.append(cousin_phen)
        father_gts_lst.append(father_gts.sum(axis=2))


    
   
    cousin_gts = np.vstack(cousin_gts_lst)
    cousin_phen = np.hstack(cousin_phen_lst)
    father_gts = np.vstack(father_gts_lst)
    pop_ids = np.hstack([np.ones(cousin_gts.shape[0], dtype=int) * i for i in range(npop)])


    


    

    X = np.full((cousin_gts.shape[0], 2, cousin_gts.shape[1]), fill_value=np.nan)
    X[:, 0, :] = cousin_gts
    par1, par2 = impute_cond_gau_gwas(cousin_gts, cousin_gts.mean(axis=0))
    X[::2,1,:] = par1.T
    X[1::2,1,:] = par2.T


    fam_labels = sum([[str(i), str(i)] for i in range(int(cousin_gts.shape[0])//2)],[])
    sib_data, sib_row_ind, sib_col_ind = build_sib_arr(fam_labels, x=1, y=1)
    n = cousin_gts.shape[0]
    tril_mat: csc_matrix = csc_matrix((sib_data, (sib_row_ind, sib_col_ind)), shape=(n,n))
    triu_mat: csc_matrix = tril(tril_mat, k=-1, format='csc')
    V=tril_mat + triu_mat.T
    sf_cousin,se = optimize(V, cousin_phen-cousin_phen.mean())
    V_cousin = sf_cousin * V + se * csc_matrix((np.ones(n), (np.arange(0, n), np.arange(0, n))), shape=(n,n))
    alpha_cousin,alpha_ses_cousin, alpha_cov_cousin = standard_gwas(V_cousin, X, cousin_phen)
    Z = alpha_cousin/alpha_ses_cousin
    return alpha_cousin,alpha_cov_cousin


def standard_gwas(V_, X_, y_, n_rel=None):
    X = X_.copy()
    y = y_.copy()
    lu = splu(V_)
    y -= y.mean()
    for i in range(0, X.shape[1]):
        X[:, i, :] = X[:, i, :] - np.mean(X[:, i, :], axis=0)
    Vinv_X: np.ndarray = sp_solve_dense3d_lu(V_, X.transpose(2,0,1))
    XT_Vinv_X: np.ndarray = np.einsum('...ij,...ik', X.transpose(2,0,1), Vinv_X)
    XT_Vinv_y: np.ndarray = np.einsum('...ij,i', X.transpose(2,0,1), lu.solve(y))
    alpha: np.ndarray = solve(XT_Vinv_X, XT_Vinv_y)
    alpha_cov: np.ndarray = np.linalg.inv(XT_Vinv_X)
    alpha_ses =  np.zeros((alpha_cov.shape[0], 2))
    for i in range(alpha_cov.shape[0]):
        alpha_ses[i, :] = np.sqrt(np.diag(alpha_cov[i,:,:]))
    return alpha, alpha_ses, alpha_cov

save = True

import sys
alpha, cov = main(nsnp=3000, nfam=5000, h2quad=0.9, h2_pop=0.5, unrel=False, Fst=0, npop=2, save=save, balding_nicols=True)
np.savez(f'impute_cond_gau_results/3000SNPs_50000cousinpairs_0.5h2pop_0Fst.npz', alpha=alpha, cov=cov)