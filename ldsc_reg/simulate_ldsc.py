'''
Runs a simple simulation from within ldsc_72
and solves it.

'''
import ldsc_72 as ld
import ldsc_71 as ld1
import numpy as np
import h5py
import glob
import time
startTime = time.time()

np.random.seed(123)

N = 100
S_size=  int(N/2)
S = np.array([np.array([[.5, 0], [0, .8]]),
    np.array([[0.5, 0], [0, 0.8]])] * S_size )
V = np.identity(2) * 0.5


model = ld.sibreg(S = S)
model.simdata(V, N, simr = True)

print("Solving Model...")

output_matrix, result = model.solve(est_init = np.zeros((2, 2)))

jkse = model.jackknife_se()

print("==============================================")
print("Output matrix: ",output_matrix)
print("Solver Output: ", result)
print("JK SE: ", jkse)
print("==============================================")

executionTime = (time.time() - startTime)
print('Execution time: ' + f'{executionTime:.2f}', " seconds")
