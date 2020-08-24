import importlib
import ldsc_72 as ld2
import numpy as np
import pandas as pd
import h5py
from scipy.stats import norm


np.random.seed(123)

N = 100
S = np.array([np.array([[5, 0], [0, 5]]),
    np.array([[2, 0], [0, 2]])] * 50 )# 50 = N/2
V = np.identity(2) * 0.8


model = ld2.sibreg(S = S)
model.simdata(V, N, simr = True)

# === Outputing to LDSC Format - Direct Effects === #

#snp = np.arange(1, len(theta_dir) + 1, 1)
se = np.sqrt(S[:, 0, 0])/N
zval = theta[:, 0]/se
pval = 2*norm.sf(np.abs(zval))
pval[np.where(pval < 1e-32)] = 1e-32


simulated_data_out = pd.DataFrame({'chromosome' : np.array([22] * N),
				   					'snp' : model.snp,
                                    'pos' : model.pos,
			           				'A1' : np.array(["A"] * N),
			            			'A2' : np.array(["C"] * N),
				     				'N' : N,
                                    'b' : model.theta[:, 0],
				     				'se' : model.S[:, 0, 0],
				      				'p' : pval})


simulated_data_out.to_csv("ldsc_reg/ldcsores/simdata2ldsc_dir.csv", sep = ' ')

# main ldscore data
ldscores = pd.DataFrame({'CHR' : np.array([22] * N),
                        'SNP' : model.snp,
                        'BP' : model.pos,
                        'L2' : model.r})
ldscores.to_csv("ldsc_reg/ldscores/22.l2.ldscore.gz",
                compression = 'gzip',
                sep = " ")

# N data
with open("ldsc_reg/ldscores/22.l2.M", "w") as f:
    f.write(str(N))

with open("ldsc_reg/ldscores/22.l2.M_5_50", "w") as f:
    f.write(str(N))
