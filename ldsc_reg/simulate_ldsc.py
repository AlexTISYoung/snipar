
import ldsc_72 as ld
import ldsc_71 as ld1
import numpy as np
import h5py
import glob
import time
startTime = time.time()

np.random.seed(123)

N = 100
S = np.array([np.array([[5, 0], [0, 5]]),
      np.array([[2, 0], [0, 2]])] * 50 )# 50 = N/2
V = np.identity(2) * 10.0

model = ld.sibreg(S = S)
model.simdata(V, N, simr = True)

print("Solving Model...")

output_matrix, result = model.solve(est_init = np.zeros((2, 2)))

print("==============================================")
print("Output matrix: ",output_matrix)
print("Solver Output: ", result)
print("==============================================")

executionTime = (time.time() - startTime)
print('Execution time: ' + f'{executionTime:.2f}', " seconds")
