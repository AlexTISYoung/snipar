import ldsc_71 as ld
import numpy as np

N = 100
S = np.array([np.array([[5, 0], [0, 5]]),
    np.array([[2, 0], [0, 2]])] * 50 )# 50 = N/2
V = np.identity(2) * 10.0

model = ld.sibreg_71(S = S)
model.simdata(V, N)

output_matrix, result = model.solve()
