import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float DTYPE_t

cpdef np.ndarray[np.double_t, ndim=3] normalize_S(np.ndarray[np.double_t, ndim=3] S, np.ndarray[np.double_t, ndim=1] norm):
    '''
    A function which normalizes a vector of S matrices
    '''
    
    cdef int N = S.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=3] S_norm = np.zeros_like(S)
    for idx in range(N):
        
        Si = S[idx]
        normi = norm[idx]
        S_norm[idx] = Si * normi
    
    return S_norm


cpdef np.ndarray[np.double_t, ndim=2] calc_inv_root(np.ndarray[np.double_t, ndim=2] S):
    '''
    A stable solver for S^{-1/2}
    '''
    
    cdef int ns = S.shape[0]
    cdef int ms = S.shape[1]
    
    if ~np.any(np.isnan(S)):
        S_eig = np.linalg.eig(S)
        l = np.zeros((ns, ms))
        np.fill_diagonal(l,np.power(S_eig[0],-0.5))
        S_inv_root = S_eig[1].dot(np.dot(l,S_eig[1].T))
    else:
        S_inv_root =  np.empty_like(S)
        S_inv_root[:] = np.nan
    return S_inv_root



cpdef np.ndarray[np.double_t, ndim=2] makeDmat(np.ndarray[np.double_t, ndim=2] S, int M):
    '''
    Makes D matrix = M [[1/sigma_1, 0], [0, 1/sigma_2]]
    '''

    cdef double sigma1 = np.sqrt(M * S[0, 0])
    cdef double sigma2 = np.sqrt(M * S[1, 1])
    
    cdef np.ndarray[np.double_t, ndim=2] Dmat = np.sqrt(M) * np.array([[1./sigma1, 0], 
                                                                        [0, 1./sigma2]])
    
    return Dmat


cpdef tuple standardize_mat(np.ndarray[np.double_t, ndim=2] V, 
                                                      np.ndarray[np.double_t, ndim=2] S, 
                                                      int M):
    '''
    Standardizes V and S matrices by constructing
    D = M [[1/sigma_1, 0], [0, 1/sigma_2]]. sigma_1 and sigma_2
    come from the S matrix provided.
    
    Then Vnew = Dmat @ V @ Dmat
    and Snew = Dmat @ S @ Dmat
    '''
    
    cdef np.ndarray[np.double_t, ndim=2] Dmat = makeDmat(S, M)
    cdef np.ndarray[np.double_t, ndim=2] Snew = Dmat @ S @ Dmat
    cdef np.ndarray[np.double_t, ndim=2] Vnew = Dmat @ V @ Dmat
    
    return Vnew, Snew


cpdef np.ndarray[np.double_t, ndim=2] V2Vmat(np.ndarray[np.double_t, ndim=1] V, int M):
    '''
    Transforms a 1 dimensional V array
    with 3 elements v1, v2 and r
    into a V matrix
    '''
    # getting important scalars
    cdef double v1 = V[0]
    cdef double v2 = V[1]
    cdef double r = V[2]

    cdef np.ndarray[np.double_t, ndim=2] Vmat = (1./M) * np.array([[v1, r * np.sqrt(v1 * v2)], 
                                                                  [r * np.sqrt(v1 * v2), v2]])

    return Vmat


cpdef np.ndarray[np.double_t, ndim=1] Vmat2V(np.ndarray[np.double_t, ndim=2] Vmat, int M):
    '''
    Makes a 2x2 V matrix
    into a 1 dimensional V
    array containing v1, v2 and
    r
    '''

    cdef double v1 = M * Vmat[0, 0]
    cdef double v2 = M * Vmat[1, 1]
    cdef double r = M * Vmat[0, 1]/np.sqrt(v1 * v2)
    cdef double r_check = M * Vmat[1, 0]/np.sqrt(v1 * v2)
    
    cdef np.ndarray[np.double_t, ndim=1] V = np.array([v1, v2, r])
    
    return V



cpdef double _log_ll(np.ndarray[np.double_t, ndim=1] V, 
                     np.ndarray[np.double_t, ndim=1] z, 
                     np.ndarray[np.double_t, ndim=2] S, 
                     double l, 
                     int N):

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
    
    cdef np.ndarray[np.double_t, ndim=2] Vmat = V2Vmat(V, N)


    Vnew, Snew = standardize_mat(Vmat, S, N)
    Sigma = Snew + l * Vnew
    logdet = np.linalg.slogdet(Sigma)
    
    print(Sigma)

    det = np.linalg.det(Sigma)
    if (det > 1e-6) or (det < -1e-6):
        Sigma_inv = np.linalg.inv(Sigma)
    else:
        Sigma_inv = np.linalg.pinv(Sigma)

    d = Vmat.shape[0]
    z2d = z.reshape(d,1)

    L = - (d/2.0) * np.log(2 * np.pi) \
        - (1.0/2.0) * logdet[0]*logdet[1] \
        - (1.0/2.0) * z2d.T @ Sigma_inv @ z2d

    logll = L[0, 0]

    return logll



cpdef _grad_ll_v(
    V, 
    z, 
    S, 
    l, 
    N
):

    """
    """

    cdef np.ndarray[np.double_t, ndim=2] Vmat = V2Vmat(V, N)
    d = S.shape[0]

    Vout = standardize_mat(Vmat, S, N)
    Vnew = Vout[0]
    Snew = Vout[1]
    
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
    
    grad = np.array([dl_dv1, dl_dv2, dl_dr])

    return grad



cpdef tuple neg_logll_grad(
    np.ndarray[np.double_t, ndim=1] V, 
    np.ndarray[np.double_t, ndim=2] z, 
    np.ndarray[np.double_t, ndim=3] S, 
    np.ndarray[np.double_t, ndim=1] l, 
    np.ndarray[np.double_t, ndim=1] u, 
    int M):

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
    cdef int N = len(S)
    
    cdef np.ndarray[np.double_t, ndim=1] G = np.zeros((3))
    cdef double log_ll = 0.0

    for i in range(N):

        Si = S[i]
        zi = z[i, :]
        ui = u[i]
        li = l[i]

        log_ll += (1/ui) * _log_ll(V, zi, Si, li, M)
        G += (1/ui) * _grad_ll_v(V, zi, Si, li, M) 

    return -log_ll , -G

# cpdef tuple neglike_wrapper(
#     np.ndarray[np.double_t, ndim=1] V,
#     np.ndarray[np.double_t, ndim=2] z,
#     np.ndarray[np.double_t, ndim=3] S,
#     np.ndarray[np.double_t, ndim=1] l,
#     np.ndarray[np.double_t, ndim=1] u,
#     np.ndarray[np.double_t, ndim=1] f,
#     double M
# ):
    
#     '''
#     Wrapper for neg_logll_grad to convert V from an
#     array of individual parameters to a symmetric
#     matrix and solve for the negative log likelihood
#     '''

#     cdef int N = S.shape[0]
    
#     cdef np.ndarray[np.double_t, ndim=1] normalizer = 2 * f * (1 - f) if f is not None else np.ones(N)
#     S = normalize_S(S, normalizer)
    
#     logll, Gvec = neg_logll_grad(V, 
#                                z, S, 
#                                l, u, M)
    
#     # print(f"Logll : {logll}")
#     # print(f"V : {V}")
    
#     return logll, Gvec



# cpdef np.ndarray[np.double_t, ndim=2] _num_grad_V(
#     np.ndarray[np.double_t, ndim=1] V,
#     np.ndarray[np.double_t, ndim=1] z,
#     np.ndarray[np.double_t, ndim=2] S,
#     np.ndarray[np.double_t, ndim=1] l,
#     double M
# ):
#     """
#     Returns numerical gradient vector of _log_ll
#     Mostly meant to check if _grad_ll_v is working
#     properly
        
#     Inputs:
#     V = dx1 numpy matrix
#     z = dx1 numpy matrix
#     S = dxd numpy matrix
#     u = 1 numpy matrix
#     r = 1 numpy matrix
#     f = 1 numpy matrix
        
#     Outputs:
#     g = dxd matrix 
#     """
    
#     cdef np.ndarray[np.double_t, ndim=2] g = np.zeros(V.shape)
#     cdef int vshape0 = V.shape[0]

#     for i in range(0,vshape0):
#         dV = np.zeros(V.shape)
#         dV[i] = 10 ** (-6)
#         V_upper = V+dV
#         V_lower = V-dV
#         g[i] = (_log_ll(V_upper, z, S, l, M) - \
#                     _log_ll(V_lower, z, S, l, M)) / (2 * 10 ** (-6))
#     return g


# def _num_grad2_V(V, z, S, l, M):
#     """
#     Calculates second derivative matrix (the Hessian) of 
#     the log likelihood at a particular observation

#     Used to calculate the standard errors of the estimates.
#     """

#     h = np.zeros((V.shape[0], V.shape[0]))

#     for i in range(V.shape[0]):
#         dV = np.zeros_like(V)
#         dV[i] = 10 ** (-6)
#         V_upper = V+dV
#         V_lower = V-dV
#         h[i, :] = (_grad_ll_v(V_upper, z, S, l, M) - \
#                     _grad_ll_v(V_lower, z, S, l, M)) / (2 * 10 ** (-6))
    
#     return h


# def _data_hessian(V, z, S, l, u, M):
#     """
#     Get hessian matrix at a particular value
#     of V across all data points
#     """

#     # Unflatten V into a matrix
#     N = len(S)
#     H = np.zeros((3, 3))

#     for i in range(N):

#         Si = S[i]
#         zi = z[i, :]
#         ui = u[i]
#         li = l[i]

#         H += (1/ui) * _num_grad2_V(V, zi, Si, li, M) 

#     return -H

# def get_hessian(V, z, S, l, u, f, M):
#     """
#     Get Hessian Matrix for dataset
#     """

#     N = S.shape[0]
    
#     normalizer = 2 * f * (1 - f) if f is not None else np.ones(N)
#     S = normalize_S(S, normalizer)
    
#     H = _data_hessian(V, z, S, l, u, M)

#     return H