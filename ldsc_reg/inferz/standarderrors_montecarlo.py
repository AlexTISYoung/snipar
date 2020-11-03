import sib_ldsc_z as ld
import numpy as np
import matplotlib.pyplot as plt


nsim = 100
np.random.seed(123)
jkse_est = np.zeros((nsim, 3))
invhess_est = np.zeros((nsim, 3))

for i in range(nsim):

    # Simulating
    N = int(1e4)
    S = np.array([[[1e-4, -5 * 1e-5], [-5 * 1e-5, 1e-4]]] * N)
    V = np.array([[0.5, 0.25], [0.25, 0.5]])

    model = ld.sibreg(S = S)
    model.simdata(V/N, N, simld = True)

    # solving
    output, _ = model.solve()

    invhess_est[i, :] = np.diag(output['std_err_mat'])
    jkse_est[i, :] = model.jackknife_se(blocksize = int(N/100))
