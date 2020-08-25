'''
This script is meant to do the "simple simulation"
from within the ldsc_72 package and convert it to a 
format which the ldsc package can read. Meant
for benchmarking
'''

import importlib
import ldsc_72 as ld2
import numpy as np
import pandas as pd
import h5py
from scipy.stats import norm
import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("N", help="The number of Loci we want", type=int)
args = parser.parse_args()

np.random.seed(123)

N = int(args.N)
print("Number of Loci:", N)
S_size = int(N/2)
# S = np.array([np.array([[.5, 0], [0, .8]]),
#     np.array([[0.9, 0], [0, 0.8]])] * S_size )
# V = np.identity(2) * 0.8
S = np.array([[[0.8]]] * N)
V = np.array([[[0.8]]])

model = ld2.sibreg(S = S)
model.simdata(V, N, simr = True)
print(model.solve(est_init = np.zeros_like(S[0])))

# === Outputing to LDSC Format - Direct Effects === #

#snp = np.arange(1, len(theta_dir) + 1, 1)
se = np.sqrt(S[:, 0, 0])/N
zval = model.theta[:, 0]/se
pval = 2*norm.sf(np.abs(zval))
pval[np.where(pval < 1e-32)] = 1e-32


simulated_data_out = pd.DataFrame({'chromosome' : np.array([22] * N),
				   					'snp' : np.array(list(map(lambda x: f"s{x}", model.snp))),
                                    'pos' : model.pos,
			           				'A1' : np.array(["A"] * N),
			            			'A2' : np.array(["C"] * N),
				     				'N' : N,
                                    'b' : model.theta[:, 0],
                                    # 'Z' : zval,
				     				'se' : model.S[:, 0, 0],
				      				'p' : pval})


simulated_data_out.to_csv("ldsc_reg/ldscores/simdata2ldsc_dir.csv", sep = ' ')


# main ldscore data
ldscores = pd.DataFrame({'CHR' : np.array([22] * N),
                        'SNP' : np.array(list(map(lambda x: f"s{x}", model.snp))),
                        'BP' : model.pos,
                        'L2' : model.r})
ldscores.to_csv("ldsc_reg/ldscores/22.l2.ldscore.gz",
                compression = 'gzip',
                sep = " ",
                index = False)

# N data
with open("ldsc_reg/ldscores/22.l2.M", "w") as f:
    f.write(str(N))

with open("ldsc_reg/ldscores/22.l2.M_5_50", "w") as f:
    f.write(str(N))


# writing empty chromosome files

for i in range(22):
    ldscores = pd.DataFrame({'CHR' : [i, i],
                            'SNP' : ["s99999", "s99999"],
                            'BP' : [99999, 99999],
                            'L2' : [99999, 99999]})
    
    ldscores.to_csv(f"ldsc_reg/ldscores/{i}.l2.ldscore.gz",
                compression = 'gzip',
                sep = " ",
                index = False)
    
    with open(f"ldsc_reg/ldscores/{i}.l2.M", "w") as f:
        f.write(str(0))

    with open(f"ldsc_reg/ldscores/{i}.l2.M_5_50", "w") as f:
        f.write(str(0))
