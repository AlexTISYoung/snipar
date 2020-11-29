import numpy as np
import glob
import pandas as pd
from scipy.optimize import minimize
import scipy.optimize
from scipy.special import comb
from scipy.misc import derivative
import scipy.stats
from numba import jit, njit, prange, int64, float64
import multiprocessing as mp
from functools import partial

def delete_obs_jk(var, start_idx, end_idx, end_cond):

    '''
    var: numpy array
    end_cond : boolean
    Function helps take out observations
    for a jackknife routine
    '''
    
    if var is not None: # in case we pass a none given f value

        if end_cond:

            var_jk = np.delete(var, range(start_idx, end_idx), 
                                axis = 0)

        else:
            
            var_jk = np.delete(var, range(start_idx, var.shape[0]), 
                                axis = 0)
            # var_jk = np.delete(var_jk, range(end_idx - var.shape[0]))
    else:
        
        var_jk = None
        
    return var_jk


@njit
def normalize_S(S, norm):
    '''
    A function which normalizes a vector of S matrices
    '''
    
    N = S.shape[0]
    S_norm = np.zeros_like(S)
    for idx in range(N):
        
        Si = S[idx]
        normi = norm[idx]
        S_norm[idx] = Si * normi
    
    return S_norm

@njit
def calc_inv_root(S):
    '''
    A stable solver for S^{-1/2}
    '''
    
    if ~np.any(np.isnan(S)):
        S_eig = np.linalg.eig(S)
        l = np.zeros(S.shape)
        np.fill_diagonal(l,np.power(S_eig[0],-0.5))
        S_inv_root = S_eig[1].dot(np.dot(l,S_eig[1].T))
    else:
        S_inv_root =  np.empty_like(S)
        S_inv_root[:] = np.nan
    return S_inv_root


@njit
def makeDmat(S, M):
    '''
    Makes D matrix = M [[1/sigma_1, 0], [0, 1/sigma_2]]
    '''

    sigma1 = np.sqrt(M * S[0, 0])
    sigma2 = np.sqrt(M * S[1, 1])
    
    Dmat = np.sqrt(M) * np.array([[1/sigma1, 0], 
                                [0, 1/sigma2]])
    
    return Dmat

@njit(parallel=True)
def makeSnew_vec(S, M):

    '''
    Makes a vector of S and a scalar value M
    into a vector of Snew = Dmat @ Si @ Dmat
    for each Si in S
    '''
    
    Snew_mat = np.zeros_like(S)

    for idx in prange(len(S)):
        Si = S[idx]
        Dmat = makeDmat(Si, M)
        Snew = Dmat @ Si @ Dmat
        Snew_mat[idx] = Snew
        
    return Snew_mat

@njit
def standardize_mat(V, S, M):
    '''
    Standardizes V and S matrices by constructing
    D = M [[1/sigma_1, 0], [0, 1/sigma_2]]. sigma_1 and sigma_2
    come from the S matrix provided.
    
    Then Vnew = Dmat @ V @ Dmat
    and Snew = Dmat @ S @ Dmat
    '''
    
    Dmat = makeDmat(S, M)
    Snew = Dmat @ S @ Dmat
    Vnew = Dmat @ V @ Dmat
    
    return Vnew, Snew

@njit
def V2Vmat(V, M):
    '''
    Transforms a 1 dimensional V array
    with 3 elements v1, v2 and r
    into a V matrix
    '''
    # getting important scalars
    v1 = V[0]
    v2 = V[1]
    r = V[2]

    Vmat = (1/M) * np.array([[v1, r * np.sqrt(v1 * v2)], [r * np.sqrt(v1 * v2), v2]])

    return Vmat

@njit 
def Vmat2V(Vmat, M):
    '''
    Makes a 2x2 V matrix
    into a 1 dimensional V
    array containing v1, v2 and
    r
    '''

    v1 = M * Vmat[0, 0]
    v2 = M * Vmat[1, 1]
    r = M * Vmat[0, 1]/np.sqrt(v1 * v2)
    r_check = M * Vmat[1, 0]/np.sqrt(v1 * v2)

    assert np.abs(r - r_check) < 10 ** -6

    return np.array([v1, v2, r])


@njit
def _log_ll(V, z, S, l, N):

    """
    Returns the log likelihood matrix for a given SNP i as formulated by:

    .. math::
        l_i = -\frac{d}{2} log (2 \pi) - \frac{1}{2} log ( |\Sigma| ) -
                \frac{1}{2} z_i^T (\Sigma) ^{-1} z_i

    Inputs:
    V = dxd numpy matrix
    z = dx1 numpy matrix
    S = dxd numpy matrix
    l = scalar
    f = scalar

    Outputs:
    logll = scalar
    """
    Vmat = V2Vmat(V, N)

    Vnew, Snew = standardize_mat(Vmat, S, N)
    Sigma = Snew + l * Vnew
    logdet = np.linalg.slogdet(Sigma)

    det = np.linalg.det(Sigma)
    if det > 1e-6 or det < -1e-6:
        Sigma_inv = np.linalg.inv(Sigma)
    else:
        Sigma_inv = np.linalg.pinv(Sigma)
    
    Sigma_inv = np.ascontiguousarray(Sigma_inv)

    d = Vmat.shape[0]
    z = z.reshape(d,1)

    L = - (d/2) * np.log(2 * np.pi) \
        - (1/2) * logdet[0]*logdet[1] \
        - (1/2) * z.T @ Sigma_inv @ z

    return L[0, 0]


@njit
def _grad_ll_v(V, z, S, l, N):

    """
    """

    Vmat = V2Vmat(V, N)
    d = S.shape[0]

    Vnew, Snew = standardize_mat(Vmat, S, N)
    Sigma = Snew + l * Vnew
    
    # getting important scalars
    v1 = V[0]
    v2 = V[1]
    r = V[2]

    rs = Snew[0, 1]

    sigma1 = np.sqrt(N * S[0, 0])
    sigma2 = np.sqrt(N * S[1, 1])

    sigma1sq = sigma1 ** 2
    sigma2sq = sigma2 ** 2
    
    assert len(z.shape) == 1
    
    z1 = z[0]
    z2 = z[1]
    
    z1sq = z1 ** 2
    z2sq = z2 ** 2

    det = np.linalg.det(Sigma)

    T = z1sq * (1 + (l * v2/sigma2sq)) - \
            2 * z1 * z2 * (rs + r * l * np.sqrt(v1 * v2)/(sigma1 * sigma2)) + \
            z2sq * (1 + l * v1/sigma1sq)

    # gradient wrt v1
    ddet_dv1 = l/sigma1sq * (1 + l * v2/sigma2sq) - (l * r/(sigma1 * sigma2)) * \
            (rs * np.sqrt(v2/v1) + r * l * v2/(sigma1 * sigma2)) 

    dT_dv1 = (z2sq * l/sigma1sq) - (r * l * z1 * z2)/(sigma1 * sigma2) * np.sqrt(v2/v1)

    dl_dv1 = -(1/2) * 1/det * (ddet_dv1 * (1 - T/det) + dT_dv1)

    # gradient wrt v2
    ddet_dv2 = l/sigma2sq * (1 + l * v1/sigma1sq) - (l * r/(sigma1 * sigma2)) * \
            (rs * np.sqrt(v1/v2) + r * l * v1/(sigma1 * sigma2)) 

    dT_dv2 = (z1sq * l/sigma2sq) - (r * l * z1 * z2)/(sigma1 * sigma2) * np.sqrt(v1/v2)

    dl_dv2 = -(1/2) * 1/det * (ddet_dv2 * (1 - T/det) + dT_dv2)

    # gradient wrt r
    ddet_dr = -2 * l * np.sqrt(v1 * v2)/(sigma1 * sigma2) * \
         (rs + (r * l * np.sqrt(v1 * v2))/(sigma1 * sigma2))

    dT_dr = -2 * l * np.sqrt(v1 * v2)/(sigma1 * sigma2) * z1 * z2

    dl_dr = -(1/2) * 1/det * (ddet_dr * (1 - T/det) + dT_dr)

    return np.array([dl_dv1, dl_dv2, dl_dr])


@njit(parallel = True, fastmath = True)
def neg_logll_grad(V, z, S, l, u, M):

    """
    Returns the loglikelihood and its gradient wrt V for a given SNP i as formulated by:

    .. math::
        l_i = -\frac{d}{2} log (2 \pi) - \frac{1}{2} log ( |I + r_i S_i^{-1/2} V S_i^{-1/2}| ) -
                \frac{1}{2} z_i^T (I + r_i S_i^{-1/2} V S_i^{-1/2}) ^{-1} z_i

    and

    .. math::
        \frac{dl}{dV} = S^{-1/2} \Sigma_i^{-1} (\Sigma - z_i z_i^T) \Sigma_i^{-1} S^{-1/2}

    Inputs:
    V = dxd numpy matrix
    z = dxN numpy matrix
    S = dxd numpy matrix
    u = 1 numpy matrix
    r = 1 numpy matrix
    f = 1 numpy matrix
    M = Scalar

    Outputs:
    -log_ll = 1x1 scalar
    -Gvec = dxd numpy matrix
    """

    # Unflatten V into a matrix
    N = len(S)
    
    G = np.zeros((3))
    log_ll = 0.0

    for i in prange(N):

        Si = S[i]
        zi = z[i, :]
        ui = u[i]
        li = l[i]

        log_ll += (1/ui) * _log_ll(V, zi, Si, li, M)
        G += (1/ui) * _grad_ll_v(V, zi, Si, li, M) 

    return -log_ll , -G

def neglike_wrapper(V, z, S, l, u, f, M):
    
    '''
    Wrapper for neg_logll_grad to convert V from an
    array of individual parameters to a symmetric
    matrix and solve for the negative log likelihood
    '''

    N = S.shape[0]
    
    normalizer = 2 * f * (1 - f) if f is not None else np.ones(N)
    S = normalize_S(S, normalizer)
    
    logll, Gvec = neg_logll_grad(V, 
                               z, S, 
                               l, u, M)
    
    # print(f"Logll : {logll}")
    # print(f"V : {V}")
    
    return logll, Gvec


@njit
def _num_grad_V(V, z, S, l, M):
    """
    Returns numerical gradient vector of _log_ll
    Mostly meant to check if _grad_ll_v is working
    properly
        
    Inputs:
    V = dx1 numpy matrix
    z = dx1 numpy matrix
    S = dxd numpy matrix
    u = 1 numpy matrix
    r = 1 numpy matrix
    f = 1 numpy matrix
        
    Outputs:
    g = dxd matrix 
    """
    
    g = np.zeros(V.shape)

    for i in range(0,V.shape[0]):
        dV = np.zeros(V.shape)
        dV[i] = 10 ** (-6)
        V_upper = V+dV
        V_lower = V-dV
        g[i] = (_log_ll(V_upper, z, S, l, M) - \
                    _log_ll(V_lower, z, S, l, M)) / (2 * 10 ** (-6))
    return g

@njit
def _num_grad2_V(V, z, S, l, M):
    """
    Calculates second derivative matrix (the Hessian) of 
    the log likelihood at a particular observation

    Used to calculate the standard errors of the estimates.
    """

    h = np.zeros((V.shape[0], V.shape[0]))

    for i in range(V.shape[0]):
        dV = np.zeros_like(V)
        dV[i] = 10 ** (-6)
        V_upper = V+dV
        V_lower = V-dV
        h[i, :] = (_grad_ll_v(V_upper, z, S, l, M) - \
                    _grad_ll_v(V_lower, z, S, l, M)) / (2 * 10 ** (-6))
    
    return h

@njit(parallel = True, fastmath = True)
def _data_hessian(V, z, S, l, u, M):
    """
    Get hessian matrix at a particular value
    of V across all data points
    """

    # Unflatten V into a matrix
    N = len(S)
    H = np.zeros((3, 3))

    for i in prange(N):

        Si = S[i]
        zi = z[i, :]
        ui = u[i]
        li = l[i]

        H += (1/ui) * _num_grad2_V(V, zi, Si, li, M) 

    return -H

def get_hessian(V, z, S, l, u, f, M):
    """
    Get Hessian Matrix for dataset
    """

    N = S.shape[0]
    
    normalizer = 2 * f * (1 - f) if f is not None else np.ones(N)
    S = normalize_S(S, normalizer)
    
    H = _data_hessian(V, z, S, l, u, M)

    return H


    
def Vinit(z, S, l, M):
    '''
    Get initial estimate to start
    solver from
    '''
    
    # compile function
    _ = makeSnew_vec(S[0:2], M)
    
    # actual run should be much faster
    Snew_mat = makeSnew_vec(S, M)
        
    Snew = np.average(Snew_mat, axis=0, weights = 1/l)
    z_var = np.cov(z.T, aweights = 1/l)
    l_bar = np.mean(l)
    v1 = (z_var[0, 0] - 1)/l_bar
    v2 = (z_var[1, 1] - 1)/l_bar
    r = (z_var[0, 1] - Snew[0, 1])/(l_bar * np.sqrt(v1 * v2))

    est_init = np.array([v1, v2, r])
    return est_init

class sibreg():
    
    def __init__(self, S, z = None, l = None, u = None,  f = None, M = None):
        
        if S.ndim > 1:
            for s in S:
                n, m = s.shape
                assert n == m

        if z is None:
            print("Warning there is no value for z. Maybe consider simulating it")
        if u is None:
            print("No value for U given. Generating a vector of ones (all SNPs weighted equally)")
            u = np.ones(S.shape[0])
        if l is None:
            print("No value for LD Scores given. Generating a vector of ones for l")
            l = np.ones(S.shape[0])
        if f is None:
            print("Warning: No value given for allele frequencies. Some parameters won't be normalized.")
        if M is None:
            print("No value for effective number of loci is given. Using total number of loci instead")
            M = len(S[~np.any(np.isnan(S), axis = (1, 2))])
        
        self.z = None if z is None else z[~np.any(np.isnan(z), axis = 1)]
        self.S = S[~np.any(np.isnan(S), axis = (1, 2))]
        self.u = u[~np.isnan(u)]
        self.l = l[~np.isnan(l)]
        self.f = None if f is None else f[~np.isnan(f)]
        self.M = M
    

    def simdata(self, V,  N, simld = False):
        """
        Simulates data for z scores.
        
        Inputs:
        V = varcov matrix of true effects
        N = Number of obs/SNPs to generate
        simr = boolean indicating if we want
                to simulate ldscores
                
        
        Outputs:
        None
        
        - It creates an object within the class
        called z
        """
        
        
        S = self.S
        
        if simld:
            self.l = np.random.uniform(low=1, high=5, size=N)
            print("Simulated LD scores!")
        
        l = self.l

        zhat_vec = np.empty((N, V.shape[1]))
        for i in range(N):
            
            Si = S[i]
            li = l[i]
            
            V = np.array(V)
            Si = np.array(Si)

            # get shape of V
            d = V.shape[0]
            zeromat = np.zeros(d)
            Vnew, Snew = standardize_mat(V, Si, N)

            # generate true effect vector
            sim = np.random.multivariate_normal(zeromat, Snew + li * Vnew)
            
            # Append to vector of effects
            zhat_vec[i, :] = sim
        
        self.snp = np.arange(1, N+1, 1)
        self.pos = np.arange(1, N+1, 1)
        self.z = zhat_vec



    def solve(self,
              z = None, 
              S = None,
              l = None,
              u = None,
              f = None,
              M = None,
              neg_logll_grad_func = neglike_wrapper,
              est_init = None,
              printout = True,
              rbounds = True):
        
        """
        Solves the ldsc problem of infering the V matrix
                
        Inputs:
        z = Nx1 numpy matrix
        S = dxd numpy matrix
        u = 1 numpy matrix
        r = 1 numpy matrix
        f = 1 numpy matrix
        
        Outputs:
        output_matrix = dxd numpy matrix
        result = result of scipy solver 
        """
        
        # inherit parameters from the class if they aren't defined
        z = self.z if z is None else z
        S = self.S if S is None else S
        l = self.l if l is None else l
        u = self.u if u is None else u
        f = self.f if f is None else f
        M = self.M if M is None else M

        # == Solves our MLE problem == #
        m = 3
        
        if est_init is None:
            if printout == True:
                print("No initial guess provided.")
                print("Making Method of Moments Guess")
                
            est_init = Vinit(z, S, u, M)
            print(f"Initial estimate: {est_init}")
        
        # exporting for potential later reference
        self.est_init = est_init
        
        rlimit = (-1, 1) if rbounds else (None, None)

        result = minimize(
            neg_logll_grad_func, 
            est_init,
            jac = True,
            args = (z, S, l, u, f, M),
            bounds = [(1e-6, None), (1e-6, None), rlimit],
            method = 'L-BFGS-B'
            # options = {'ftol' : 1e-20}
            
        )

        output = dict(
            v1 = result.x[0],
            v2 = result.x[1],
            r = result.x[2]
        )

        
        # re-normnalizing output matrix 
        self.output = output

        # Getting Inverse Hessian
        H = get_hessian(result.x, z, S, l, u, f, M)
        invH = np.linalg.inv(H)
        std_err_mat = np.sqrt(invH)
        
        output["std_err_mat"] = std_err_mat
        
        return output, result 

# == Parallelized Jackknife == #

def jkse_core(indices,
             model,
             full_est):
    
    '''
    This runs the core estimation
    for each jackknife iteration
    '''
    
    mask = np.ones(model.z.shape[0], dtype=bool)
    mask[indices] = False
    
    z = model.z[mask]
    S = model.S[mask]
    l = model.l[mask]
    u = model.u[mask]
    f = model.f[mask] if model.f is not None else model.f
    M = model.M
    
    output, _ = model.solve(z = z,
                    S = S,
                    l = l,
                    u = u,
                    f = f,
                    M = M,
                    printout = False,
                    est_init = full_est)
    
    
    output_matrix = np.array([output['v1'], output['v2'], output['r']])
    
    return output_matrix

    

def jkse(model,
        full_est_params,
        blocksize = 1,
        printinfo = False,
        num_procs = 2):
    
    '''
    This runs the whole block jackknife 
    routine in parallel
    '''
    
    nblocks = int(np.ceil(model.z.shape[0]/blocksize))
    
    # construct blocks of indices
    indices = list(range(1, model.z.shape[0]))
    index_blocks = [indices[i * blocksize:(i + 1) * blocksize] for i in range((len(indices) + blocksize - 1) // blocksize )]
    
    # store full parameter estimate as array
    full_est = np.array([full_est_params['v1'], full_est_params['v2'], full_est_params['r']])
    
    jkse_toparallelize = partial(jkse_core, model = model, full_est = full_est)
    
    num_procs = num_procs
    pool = mp.Pool(num_procs)
    estimates_jk = pool.map(jkse_toparallelize, index_blocks)
    estimates_jk = np.array(estimates_jk)
    
    pseudovalues = nblocks * full_est - (nblocks - 1) * estimates_jk

    # calculate jackknife se
    jknife_cov = np.cov(pseudovalues.T, ddof=1) / nblocks
    jknife_var = np.diag(jknife_cov)
    jknife_se = np.sqrt(jknife_var)

    return jknife_se  
    

#  == Reading LD scores == #


def M(fh, N=2, common=False):
    '''Parses .l{N}.M files, split across num chromosomes. See docs/file_formats_ld.txt.'''
    parsefunc = lambda y: [float(z) for z in open(y, 'r').readline().split()]
    suffix = '.l' + str(N) + '.M'
    suffix += '_5_50'

    x = parsefunc(fh + suffix)

    return np.array(x).reshape((1, len(x)))


def read_ldscores(ldscore_path, ldcolnames):
    '''
    Reads in LD scores
    ldscore_path : string signifying where the ldscores are
    ldfiles : string - a glob identifier for all the ld score files
    ldcolnames : list - the columns present in the ldscore data
    Mfiles : string - a glob identifier for all the files which has M data
    Mcolnames : list - the columns present in the M file data
    '''
    files = glob.glob(f"{ldscore_path}")
    ldscores = pd.DataFrame(columns = ldcolnames)

    for file in files:
        snpi = pd.read_csv(file, compression='gzip', sep = "\t")
        ldscores = pd.concat([ldscores, snpi], sort = False)
        
    return ldscores

def read_mfiles(Mfilepath, Mcolnames):
    # Reading Number of Loci
    files_M = glob.glob(f"{Mfilepath}")
    nloci = pd.DataFrame(columns = Mcolnames)

    for file in files_M:
        
        chrom = int(file[-12:-10].strip('/'))
        
        nloci_i = pd.DataFrame({"M" : [M(file[:-10])[0, 0]],
                                "CHR" : [chrom]})
        
        nloci = pd.concat([nloci, nloci_i])

    nloci = nloci.reset_index(drop = True)

    return nloci


# == Some useful transformations == #

def transform_estimates(effect_estimated,
                        S, theta):
    '''
    Transforms theta (can also be z) and S data into
    the required format
    Format types can be:
    - population (1 dimensional)
    - direct_plus_averageparental (2 dimensional)
    - direct_plus_population (2 dimensional)
    - full (3 dimensional)

    Assumes input data is 3 dimensional
    '''
    if effect_estimated == "population":
        # == Keeping population effect == #
        Sdir = np.empty(len(S))
        for i in range(len(S)):
            Sdir[i] = np.array([[1.0, 0.5, 0.5]]) @ S[i] @ np.array([[1.0, 0.5, 0.5]]).T

        S = Sdir.reshape((len(S), 1, 1))
        theta = theta @ np.array([1.0, 0.5, 0.5])
        theta = theta.reshape((theta.shape[0], 1))
    elif effect_estimated == "direct_plus_averageparental":

        # == Combining indirect effects to make V a 2x2 matrix == #
        tmatrix = np.array([[1.0, 0.0],
                            [0.0, 0.5],
                            [0.0, 0.5]])
        Sdir = np.empty((len(S), 2, 2))
        for i in range(len(S)):
            Sdir[i] = tmatrix.T @ S[i] @ tmatrix
        S = Sdir.reshape((len(S), 2, 2))
        theta = theta @ tmatrix
        theta = theta.reshape((theta.shape[0], 2))
    elif effect_estimated == "direct_plus_population":

        # == keeping direct effect and population effect == #
        tmatrix = np.array([[1.0, 1.0],
                            [0.0, 0.5],
                            [0.0, 0.5]])
        Sdir = np.empty((len(S), 2, 2))
        for i in range(len(S)):
            Sdir[i] = tmatrix.T @ S[i] @ tmatrix

        S = Sdir.reshape((len(S), 2, 2))
        theta = theta @ tmatrix
        theta = theta.reshape((theta.shape[0], 2))
    elif effect_estimated == "full":
        pass


    return S, theta


@njit
def theta2z(theta, S, M):
    '''
    Transforms vector of theta values
    into z values based on S
    '''
    zval = np.empty_like(theta)
    for i in range(theta.shape[0]):
        Dmat = makeDmat(S[i], M)
        zval[i, :] = Dmat @ theta[i, :]

    return zval