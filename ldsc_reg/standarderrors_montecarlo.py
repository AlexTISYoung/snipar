import sib_ldsc_z as ld
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import seaborn as sns

nsim = 100
def simulate_se(N):

    # Simulating
    S = np.array([[[1, -5 * 1e-1], [-5 * 1e-1, 1]]] * N)/N
    V = np.array([[0.5, 0.25], [0.25, 0.5]])

    model = ld.sibreg(S = S)
    model.simdata(V/N, N, simld = True)

    # solving
    output, _ = model.solve()

    invhess_est = np.diag(output['std_err_mat'])
    jkse_est = ld.jkse(model, output, blocksize = int(N/100), num_procs = 4)
    estimates = np.array([output['v1'], output['v2'], output['r']])

    return [invhess_est, jkse_est, estimates]

def sim_data(nsim = 100):

    inv_hess = np.zeros((nsim, 3))
    jkse = np.zeros((nsim, 3))
    estimates = np.zeros((nsim, 3))
    for i in range(nsim):
        inv_hess[i, :], jkse[i, :], estimates[i, :] = simulate_se(N = int(1e3))

    return inv_hess, jkse, estimates


np.random.seed(123)
inv_hess, jkse, estimates = sim_data()
inv_hess = inv_hess/np.sqrt(2)

# plotting
fig, ax = plt.subplots(3, 1, figsize = (8, 16))
sns.distplot(inv_hess[:, 0], ax=ax[0])
sns.distplot(jkse[:, 0], ax=ax[0])
ax[0].vlines(np.std(estimates[:, 0]), 0, 250, linestyle = "--")
ax[0].set_xlabel("V1")
ax[0].legend(
    labels=["True Standard Errors", "Inverse Hessian","Block Jackknife"],
    loc='upper left'
)
sns.distplot(inv_hess[:, 1], ax=ax[1])
sns.distplot(jkse[:, 1], ax=ax[1])
ax[1].vlines(np.std(estimates[:, 1]), 0, 250, linestyle = "--")
ax[1].set_xlabel("V2")

sns.distplot(inv_hess[:, 2], ax=ax[2])
sns.distplot(jkse[:, 2], ax=ax[2])
ax[2].vlines(np.std(estimates[:, 2]), 0, 250, linestyle = "--")
ax[2].set_xlabel("r")

plt.savefig("ldsc_reg/inferz/standard_errors.png")