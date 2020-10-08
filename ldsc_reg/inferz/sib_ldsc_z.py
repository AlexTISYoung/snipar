import numpy as np
from scipy.optimize import minimize
import scipy.optimize
from scipy.special import comb
from scipy.misc import derivative
import scipy.stats
from numba import jit, njit, prange

def extract_upper_triangle(x):
    
    # =============================== #
    # Extracts the upper triangular portion of 
    # a symmetric matrix
    # =============================== #
    
    n, m = x.shape
    assert n == m
    
    upper_triangle = x[np.triu_indices(n)]
    
    return upper_triangle

def return_to_symmetric(triangle_vec, final_size):
    
    # =============================== #
    # Given a vector of the upper triangular matrix,
    # get back the symmetric matrix
    # =============================== #
    
    X = np.zeros((final_size,final_size))
    X[np.triu_indices(X.shape[0], k = 0)] = triangle_vec
    X = X + X.T - np.diag(np.diag(X))
    
    return X

def extract_bounds(n):
    
    # =============================== #
    # From a number n, the function
    # outputs a list of bounds
    # for a var cov matrix of size
    # n x n
    # =============================== #
    
    # extract idx of flat array whcih are diagonals
    uptriangl_idx = np.array(np.triu_indices(n))
    diags = uptriangl_idx[0, :] == uptriangl_idx[1, :]
    
    # Construct list of bounds
    bounds_list = np.array([(None, None)] * len(diags))
    bounds_list[diags] = (1e-6, None)
    
    bounds_list_out = [tuple(i) for i in bounds_list]
    
    return bounds_list_out


def delete_obs_jk(var, start_idx, end_idx, end_cond):

    # ============================== #
    # var: numpy array
    # end_cond : boolean
    # Function helps take out observations
    # for a jackknife routine
    # ============================== #

    if end_cond:

        var_jk = np.delete(var, range(start_idx, end_idx), 
                             axis = 0)

    else:

        var_jk = np.delete(var, range(start_idx, var.shape[0]), 
                             axis = 0)
        var_jk = np.delete(var_jk, range(end_idx - var.shape[0]))
        
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
def _log_ll(V, z, S, r):

    """
    Returns the log likelihood matrix for a given SNP i as formulated by:

    .. math::
        l_i = -\frac{d}{2} log (2 \pi) - \frac{1}{2} log ( |I + r_i S_i^{-1/2} V S_i^{-1/2}| ) -
                \frac{1}{2} z_i^T (I + r_i S_i^{-1/2} V S_i^{-1/2}) ^{-1} z_i

    Inputs:
    V = dxd numpy matrix
    z = dx1 numpy matrix
    S = dxd numpy matrix
    r = scalar
    f = scalar

    Outputs:
    logll = scalar
    """

    S_inv_root = calc_inv_root(S)
    Sigma = np.identity(S.shape[0]) + r * S_inv_root @ V @ S_inv_root
    logdet = np.linalg.slogdet(Sigma)

    det = np.linalg.det(Sigma)
    if det > 1e-6 or det < -1e-6:
        Sigma_inv = np.linalg.inv(Sigma)
    else:
        Sigma_inv = np.linalg.pinv(Sigma)

    z = z.reshape(V.shape[0],1)
    d = V.shape[0]

    L = - (d/2) * np.log(2 * np.pi) \
        - (1/2) * logdet[0]*logdet[1] \
        - (1/2) * z.T @ Sigma_inv @ z

    return L[0, 0]


@njit
def _grad_ll_v(V, z, S, r):

    """
    Returns the gradient of the log likelihood wrt V for a given SNP i as formulated by:

    .. math::
        \frac{dl}{dV} = S^{-1/2} \Sigma_i^{-1} (\Sigma - z_i z_i^T) \Sigma_i^{-1} S^{-1/2}

    Inputs:
    V = dxd numpy matrix
    z = dx1 numpy matrix
    S = dxd numpy matrix
    u = 1 numpy matrix
    r = 1 numpy matrix
    f = 1 numpy matrix

    Outputs:
    grad_ll_v = dxd matrix 
    """

    S_inv_root = calc_inv_root(S)
    Sigma = np.identity(S.shape[0])+r * S_inv_root @ V @ S_inv_root

    det = np.linalg.det(Sigma)
    if det > 1e-6 or det < -1e-6:
        Sigma_inv = np.linalg.inv(Sigma)
    else:
        Sigma_inv = np.linalg.pinv(Sigma)

    z = z.reshape(z.shape[0],1)
    g = r * S_inv_root @ Sigma_inv @ (Sigma - z @ z.T) @ Sigma_inv @ S_inv_root
    return -(1/2) * g


@njit(parallel = True)
def neg_logll_grad(V, 
                   z, S, 
                   u, r):

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
    logllfunc = function which calculates logll
                (uses self._log_ll by default)
    gradfunc = function which calculated grad of logll
                (uses self._grad_ll_v by default)

    Outputs:
    -log_ll = 1x1 scalar
    -Gvec = dxd numpy matrix
    """

    # Unflatten V into a matrix
    d = S[0].shape[0]
    N = len(S)
    
    Gvec = np.zeros((N, d, d))
    log_ll = np.zeros(N)

    # Normalizg variables
    V_norm = V
    for i in prange(N):

        Si = S[i]
        zi = z[i, :].reshape((d, 1))
        ui = u[i]
        ri = r[i]

        Si = N * Si

        log_ll[i] = (1/ui) * _log_ll(V_norm , zi, Si, ri)
        Gvec[i, :, :] = (1/ui) * _grad_ll_v(V_norm, zi, Si, ri)

    return -log_ll.sum() , -Gvec.sum(axis = 0)


def neglike_wrapper(V, z, S, u, r, f):
    
    '''
    Wrapper for neg_logll_grad to convert V from an
    array of individual parameters to a symmetric
    matrix and solve for the negative log likelihood
    '''
    d = S[0].shape[0]
    N = S.shape[0]
    
    V = return_to_symmetric(V, d)
    
    normalizer = 2 * f  * (1 - f) if f is not None else np.ones(N)
    S = normalize_S(S, normalizer)
    
    logll, Gvec = neg_logll_grad(V, 
                               z, S, 
                               u, r)
    
    Gvec = extract_upper_triangle(Gvec)
    
    print(f"Logll : {logll}")
    print(f"V : {V}")
    
    return logll, Gvec


def _num_grad_V(V, z, S, r):
    """
    Returns numerical gradient vector of self._log_ll
    Mostly meant to check if self._grad_ll_v is working
    properly
        
    Inputs:
    V = dxd numpy matrix
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
        for j in range(0,V.shape[1]):
            dV = np.zeros((V.shape))
            dV[i,j] = 10 ** (-6)
            V_upper = V+dV
            V_lower = V-dV
            g[i,j] = (_log_ll(V_upper, z, S, r) - \
                        _log_ll(V_lower, z, S, r)) / (2 * 10 ** (-6))
    return g

class sibreg():
    
    def __init__(self, S, z = None, u = None, r = None, f = None):
        
        if S.ndim > 1:
            for s in S:
                n, m = s.shape
                assert n == m

        if z is None:
            print("Warning there is no value for z. Maybe consider simulating it")
        if u is None:
            print("No value for U given. Generating a vector of ones (all SNPs weighted equally)")
            u = np.ones(S.shape[0])
        if r is None:
            print("No value for r given. Generating a vector of ones for r")
            r = np.ones(S.shape[0])
        if f is None:
            print("Warning: No value given for allele frequencies. Some parameters won't be normalized.")
        
        self.z = None if z is None else z[~np.any(np.isnan(z), axis = 1)]
        self.S = S[~np.any(np.isnan(S), axis = (1, 2))]
        self.u = u[~np.isnan(u)]
        self.r = r[~np.isnan(r)]
        self.f = None if f is None else f[~np.isnan(f)]
    

    def simdata(self, V,  N, simr = False):
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
        
        if simr:
            self.r = np.random.uniform(low=1, high=5, size=N)
            print("Simulated LD scores!")
        
        r = self.r

        zhat_vec = np.empty((N, V.shape[1]))
        for i in range(N):
            
            Si = S[i]
            ri = r[i]
            
            V = np.array(V)
            Si = np.array(Si)
            S_inv = calc_inv_root(Si)
  
            # get shape of V
            d = V.shape[0]
            zeromat = np.zeros(d)

            # generate true effect vector
            if d > 1:
                sim = np.random.multivariate_normal(zeromat, np.eye(d) + ri * S_inv @ V @ S_inv)
            else:
                sim = np.random.normal(zeromat, np.eye(d) + ri * S_inv @ V @ S_inv)
            
            # Append to vector of effects
            zhat_vec[i, :] = sim
        

        print("Effect Vectors Simulated!")
        
        self.snp = np.arange(1, N+1, 1)
        self.pos = np.arange(1, N+1, 1)
        self.z = zhat_vec



    def solve(self,
              z = None, 
              S = None,
              u = None,
              r = None,
              f = None,
              neg_logll_grad = neglike_wrapper,
              est_init = None,
              printout = True):
        
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
        u = self.u if u is None else u
        r = self.r if r is None else r
        f = self.f if f is None else f

        # == Solves our MLE problem == #
        n, m = z.shape
        
        if est_init is not None:
            # Shape of initial varcov guess
            rowstrue = est_init.shape[0] == m
            colstrue = est_init.shape[1] == m

            if rowstrue & colstrue:
                pass
            else:
                if printout == True:
                    print("Warning: Initial Estimate given is not of the proper dimension")
                    print("Making 'optimal' matrix")
                    print("=================================================")
                
                est_init = np.zeros((m, m))
        else:
            if printout == True:
                print("No initial guess provided.")
                print("Making 'optimal' matrix")
                print("=================================================")
            
            est_init = np.zeros((m, m))
            
        
        # exporting for potential later reference
        self.est_init = est_init

        # extract array from est init
        est_init_array = extract_upper_triangle(est_init) 
        
        bounds = extract_bounds(m)     

        result = minimize(
            neg_logll_grad, 
            est_init_array,
            jac = True,
            args = (z, S, u, r, f),
            bounds = bounds,
            method = 'L-BFGS-B',
            options = {'ftol' : 1e-15}
            
        )

        output_matrix = return_to_symmetric(result.x, m)
        
        # re-normnalizing output matrix 
        self.output_matrix = output_matrix
        
        return output_matrix, result 

    def jackknife_se(self,
                  theta  = None, S = None,
                  r = None, u = None,
                  blocksize = 1):

        # Simple jackknife estimator for SE
        # Ref: https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/ldscore/jackknife.py#L231
        # Default value of blocksize = 1 is the normal jackknife

        theta = self.theta if (theta is None) else theta
        S = self.S if (S is None) else S
        r = self.r if (r is None) else r
        u = self.u if (u is None) else u

        
        assert theta.shape[0] == S.shape[0]

        nobs = theta.shape[0]
        
        estimates_jk = []
        
        start_idx = 0
        while True:
            
            end_idx = start_idx + blocksize
            end_idx_cond = end_idx <= theta.shape[0]
            
            # remove blocks of observations

            vars_jk = []

            for var in [theta, S, r, u]:

                var_jk = delete_obs_jk(var, start_idx, end_idx,
                                       end_idx_cond)
                vars_jk.append(var_jk)
            
            if start_idx < theta.shape[0]:
                # Get our estimate
                output_matrix, _ = self.solve(theta = vars_jk[0],
                                              S = vars_jk[1],
                                              r = vars_jk[2],
                                              u = vars_jk[3],
                                              printout = False,
                                              est_init = self.est_init)

                estimates_jk.append(output_matrix)

                start_idx += blocksize
            else:
                break
            
        estimates_jk = np.array(estimates_jk)
        full_est = self.output_matrix
        
        # calculate pseudo-values
        n_blocks = int(nobs/blocksize)
        pseudovalues = n_blocks * full_est - (n_blocks - 1) * estimates_jk
        
        # calculate jackknife se
        pseudovalues = pseudovalues.reshape((n_blocks, theta.shape[1] * theta.shape[1]))
        jknife_cov = np.cov(pseudovalues.T, ddof=1) / n_blocks
        jknife_var = np.diag(jknife_cov)
        jknife_se = np.sqrt(jknife_var)
    
        jknife_se  = jknife_se.reshape((theta.shape[1], theta.shape[1]))
        
        return jknife_se  
