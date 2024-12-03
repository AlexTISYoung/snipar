import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from sklearn.model_selection import LeaveOneOut
import pandas as pd


def plot(alpha, alpha_cov, non_sampling=False):
    s = np.zeros(alpha.shape[0])
    loo = LeaveOneOut()
    loo.get_n_splits(alpha)
    for i, (index, _) in enumerate(loo.split(alpha)):
        if alpha_cov.ndim == 1:
            if not non_sampling:
                s[i] = np.var(alpha[index, 0] / np.sqrt(alpha_cov)[index]) #/ (np.var(alpha_max[index, 0] / np.sqrt(alpha_cov_max)[index,0,0]))
            else:
                s[i] = np.mean(alpha[index,0] ** 2 - alpha_cov[index]) #/ (np.mean(alpha_max[index,0] ** 2) - np.mean(alpha_cov_max[index,0,0]))
        else:
            if not non_sampling:
                s[i] = np.var(alpha[index, 0] / np.sqrt(np.diagonal(alpha_cov, axis1=1, axis2=2))[index,0]) #/ (np.var(alpha_max[index, 0] / np.sqrt(alpha_cov_max)[index,0,0]))
            else:
                s[i] = (np.mean(alpha[index,0] ** 2) - np.mean(alpha_cov[index,0,0])) #/(np.mean(alpha_max[index,0] ** 2) - np.mean(alpha_cov_max[index,0,0]))
    mean = np.mean(s)
    std = np.sqrt(np.sum(np.power(mean - s, 2)) * (s.shape[0] - 1) / s.shape[0])
    return mean, std

df = []
for F in [0, 0.001, 0.01, 0.1]:
    x = np.load(f'impute_cond_gau_results/3000SNPs_5000cousinpairs_0.5h2pop_{F}Fst.npz')
    m, std, = plot(x['alpha'], x['cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['cond. Gaussian'],
                            'non_sampling_var': m, 'non_sampling_var_std': std}))
    plot(x['alpha'], x['cov'], True)

pd.concat(df).to_csv('impute_cond_gau_results/sp_non_sampling_var.csv', index=False, sep='\t')
