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
print("Making CSV for Average Parental Effects")
print("=====================================")
# reading in  data
files = glob.glob("/disk/genetics/ukb/alextisyoung/vcinf/1/chr_*.hdf5")
for file in files:

    print("Reading in file: ", file)
    hf = h5py.File(file, 'r')
    metadata = hf.get('bim')[()]
    chromosome = metadata[:, 0]
    chr_num = str(chromosome[0]).replace("b'", "")[:-1]  
    snp = metadata[:, 1]
    pos = metadata[:, 3]
    A1 = metadata[:, 4]
    A2 = metadata[:, 5]
    N = hf.get('N_L')[()]

    simulated_data_out = pd.DataFrame({'chromosome' : chromosome,
                                    'snp' : snp,
                                    'pos' : pos,
                                    'A1' : A1,
                                    'A2' : A2,
                                    'N' : N})

    simulated_data_out['chromosome'] = simulated_data_out['chromosome'].astype(float)
    simulated_data_out['snp'] = simulated_data_out['snp'].astype(str).str.replace("b'", "").str[:-1]
    simulated_data_out['pos'] = simulated_data_out['pos'].astype(str).str.replace("b'", "").str[:-1]
    simulated_data_out['A1'] = simulated_data_out['A1'].astype(str).str.replace("b'", "").str[:-1]
    simulated_data_out['A2'] = simulated_data_out['A2'].astype(str).str.replace("b'", "").str[:-1]

    simulated_data_out.to_csv(f"ldsc_reg/inferz/ldsc_reg/{chr_num}.bim", sep = ' ')




# == Making CSV == #





