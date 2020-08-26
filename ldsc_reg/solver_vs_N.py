'''
Runs solver for multiple values of N and plots it
'''

import ldsc_72 as ld
import numpy as np
import pandas as pd
import h5py
import glob
import time

df = pd.DataFrame(columns = ["var1", "var2", "cov", "jkse_var1", "jkse_var2", "jkse_cov", "N"])
np.random.seed(123)
for N in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7]:


    np.random.seed(123)
    print(f"Simulating for N = {N} ...")
    N = int(N)
    S_size = int(N/2)
    S = np.array([np.array([[.5, 0], [0, .8]]),
        np.array([[0.5, 0], [0, 0.8]])] * S_size )
    V = np.identity(2) * 0.8


    model = ld.sibreg(S = S)
    model.simdata(V, N, simr = True)

    output_mat, result = model.solve(est_init = np.zeros((2, 2)))

    print(result)


    blocksize = 1 if N < 1e4 else 200
    jkse_mat = model.jackknife_se(blocksize = blocksize)


    var1 = output_mat[0, 0]
    var2 = output_mat[1, 1]
    cov = output_mat[0, 1]


    sevar1 = output_mat[0, 0]
    sevar2 = output_mat[1, 1]
    secov = output_mat[0, 1]

    df = df.append({"var1" : var1,
                   "var2" : var2, 
                   "cov" : cov, 
                   "jkse_var1" : sevar1, 
                   "jkse_var2" : sevar2, 
                   "jkse_cov" : secov, 
                   "N" : N},
                   ignore_index = True)



