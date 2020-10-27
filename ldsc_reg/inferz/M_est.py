'''
Generates a plot of M vs Estimators to understnad
how changing the M used while estimating changes the
estimated parameters
'''

import numpy as np
import sib_ldsc_z as ld
from scipy.optimize import minimize
from scipy.special import comb
from scipy.misc import derivative
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(123)
N = int(1e4)
M = N
S = np.array([[[ 1. , -0.5],[-0.5,  1. ]]] * N)/M
V = np.array([[0.5, 0.25], [0.25, 0.5]])/M

v1, v2, r = ld.Vmat2V(V, M)


model = ld.sibreg(S = S)
model.simdata(V, N, M, simld = True)

# Solving for different M's
M_range = np.linspace(M-5000, M + 5000, 100)
output_mat = np.zeros((100, 3))

for idx, m in enumerate(M_range):
    
    print("Loop number: ", idx)
    output, result = model.solve(M = m)
    if result.success:
        output_mat[idx, 0] = output['v1']
        output_mat[idx, 1] = output['v2']
        output_mat[idx, 2] = output['r']
    else:
        output_mat[idx, 0] = np.nan
        output_mat[idx, 1] = np.nan
        output_mat[idx, 2] = np.nan


# Plotting
fig, ax = plt.subplots(3, figsize=(9,13))
fig.suptitle("M vs Estimates")
ax[0].plot(M_range/M, output_mat[:, 0])
ax[0].axhline(v1, color = "red", linestyle=':', label = "True V1")
ax[0].set_ylabel("$V_1$")
ax[1].plot(M_range/M, output_mat[:, 1])
ax[1].axhline(v2, color = "red", linestyle=':', label = "True V2")
ax[1].set_ylabel("$V_2$")
ax[2].plot(M_range/M, output_mat[:, 2])
ax[2].axhline(r, color = "red", linestyle=':', label = "True r")
ax[2].set_ylabel("$r$")
ax[2].set_xlabel("M used in Estimation/True M")
plt.savefig("ldsc_reg/inferz/est_vs_M.png")
