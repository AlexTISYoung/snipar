import importlib
import ldsc_72 as ld2
import numpy as np
import pandas as pd
import h5py
from helperfuncs import *


np.random.seed(123)

N = 100
S = np.array([np.array([[5, 0], [0, 5]]),
    np.array([[2, 0], [0, 2]])] * 50 )# 50 = N/2
V = np.identity(2) * 0.8


model = ld2.sibreg(S = S)
model.simdata(V, N, simr = True)


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
