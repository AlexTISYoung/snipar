import numpy as np

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