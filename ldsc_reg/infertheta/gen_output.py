'''
Reads in Alex's simulations and generates
output which can be read by the ldsc package.

Doesn't output the LD scores. LD scores
are meant to be read from the baseline European
sample provided by LDSC or something similar
'''
import numpy as np
import pandas as pd
import h5py
import glob
from scipy.stats import norm


# == Direct Effect == #
print("=====================================")
print("Making CSV for Direct Effects")
print("=====================================")
# reading in  data
files = glob.glob("/disk/genetics/ukb/alextisyoung/vcinf/1/chr_*.hdf5")

file = files[0]
print("Reading in file: ", file)
hf = h5py.File(file, 'r')
metadata = hf.get('bim')[()]
chromosome = metadata[:, 0]
snp = metadata[:, 1]
pos = metadata[:, 3]
A1 = metadata[:, 4]
A2 = metadata[:, 5]
theta  = hf.get('estimate')[()][:, 0]
se  = hf.get('estimate_ses')[()][:, 0]
N = hf.get('N_L')[()]
S = hf.get('estimate_covariance')[()][:, 0, 0]
f = hf.get('freqs')[()]

for file in files[1:]:
    print("Reading in file: ", file)
    hf = h5py.File(file, 'r')
    metadata = hf.get('bim')[()]
    chromosome_file = metadata[:, 0]
    snp_file = metadata[:, 1]
    pos_file = metadata[:, 3]
    A1_file = metadata[:, 4]
    A2_file = metadata[:, 5]
    theta_file  = hf.get('estimate')[()][:, 0]
    se_file  = hf.get('estimate_ses')[()][:, 0]
    S_file = hf.get('estimate_covariance')[()][:, 0, 0]
    f_file = hf.get('freqs')[()]
    N_file = hf.get('N_L')[()]

    chromosome = np.append(chromosome, chromosome_file, axis = 0)
    snp = np.append(snp, snp_file, axis = 0)
    pos = np.append(pos, pos_file, axis = 0)
    A1 = np.append(A1, A1_file, axis = 0)
    A2 = np.append(A2, A2_file, axis = 0)
    theta = np.append(theta, theta_file, axis = 0)
    se = np.append(se, se_file, axis = 0)
    S = np.append(S, S_file, axis = 0)
    f = np.append(f, f_file, axis = 0)
    N = np.append(N, N_file, axis = 0)


# getting p value
zval = theta/se
pval = 2*norm.sf(np.abs(zval))


simulated_data_out = pd.DataFrame({'chromosome' : chromosome,
                                    'snp' : snp,
                                    'pos' : pos,
                                    'A1' : A1,
                                    'A2' : A2,
                                    'N' : N,
                                    'b' : theta,
                                    'se' : se,
                                    'p' : pval})

simulated_data_out['chromosome'] = simulated_data_out['chromosome'].astype(float)
simulated_data_out['snp'] = simulated_data_out['snp'].astype(str).str.replace("b'", "").str[:-1]
simulated_data_out['pos'] = simulated_data_out['pos'].astype(str).str.replace("b'", "").str[:-1]
simulated_data_out['A1'] = simulated_data_out['A1'].astype(str).str.replace("b'", "").str[:-1]
simulated_data_out['A2'] = simulated_data_out['A2'].astype(str).str.replace("b'", "").str[:-1]


simulated_data_out.to_csv("ldsc_reg/simulated_data_dir.csv", sep = ' ')
