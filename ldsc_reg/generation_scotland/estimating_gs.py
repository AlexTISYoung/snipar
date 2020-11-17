'''
This script reads the HDF5 files for Generation Scotland
and estimates V

Trait: HDL
'''

import numpy as np
import h5py
import glob
import datetime
import matplotlib.pyplot as plt
import pandas as pd
import logging


import sys
import os
sys.path.append('/homes/nber/harij/gitrepos/SNIPar/ldsc_reg')
import sib_ldsc_z as ld

logging.basicConfig(filename= f"ldsc_reg/generation_scotland/estimating_gs_education.log", 
                         level = logging.INFO,
                         format = "%(message)s",
                         filemode = "w")

# getting all trait codes
traitnos = glob.glob("/disk/genetics/ukb/alextisyoung/GS20k_sumstats/traits/*/")
traitnos = [int(i[-3:-1].strip('/')) for i in traitnos]

traitcodes = pd.read_csv("/disk/genetics/ukb/alextisyoung/GS20k_sumstats/traits/traits.txt", 
            names = ['traitcode', 'trait'], 
            sep = ' ')

# == Reading in LD Scores == #

ldscore_path = "/disk/genetics/ukb/alextisyoung/GS20k_sumstats/ldscores/"
ldfiles = "*[0-9].l2.ldscore.gz"
ldcolnames = ["CHR", "SNP", "BP", "L2"]
Mfiles = "*[0-9].l2.M_5_50"
Mcolnames = ["M", "CHR"]

ldscores, mfile = ld.read_ldscores(ldscore_path, ldfiles, ldcolnames, Mfiles, Mcolnames)
ldscores['BP'] = ldscores['BP'].astype('int')

# trait 7 gives weird results
traitnos.remove(7)
for traitcode in traitnos:
    
    if traitcode == 7:
        traitname = "unknown"
    else:
        traitname = traitcodes.loc[traitcodes['traitcode'] == traitcode, 'trait'].values[0]

    print(f"Trait name: {traitname}")


    startTime = datetime.datetime.now()  
    logging.info(f"===============================")
    logging.info(f"Trait name:  {traitname}")
    logging.info(f"Start time:  {startTime}")

    # == Reading in Data == #
    print("=====================================")
    print("Reading in Data")
    print("=====================================")
    # reading in  data
    files = glob.glob(f"/disk/genetics/ukb/alextisyoung/GS20k_sumstats/traits/{traitcode}/chr_*.hdf5")

    file = files[0]
    print("Reading in file: ", file)
    hf = h5py.File(file, 'r')
    metadata = hf.get("bim")[()]
    chromosome = metadata[:, 0]
    snp = metadata[:, 3]
    theta  = hf.get('estimate')[()]
    se  = hf.get('estimate_ses')[()]
    N = hf.get('N_L')[()]
    S = hf.get('estimate_covariance')[()]
    f = hf.get('freqs')[()]

    # normalizing S
    sigma2 = hf.get('sigma2')[()]
    tau = hf.get('tau')[()]
    phvar = sigma2+sigma2/tau


    for file in files[1:]:
        print("Reading in file: ", file)
        hf = h5py.File(file, 'r')
        metadata = hf.get("bim")[()]
        chromosome_file = metadata[:, 0]  
        snp_file = metadata[:, 3]
        theta_file  = hf.get('estimate')[()]
        se_file  = hf.get('estimate_ses')[()]
        S_file = hf.get('estimate_covariance')[()]
        f_file = hf.get('freqs')[()]
        N_file = hf.get('N_L')[()]

        # normalizing S
        sigma2 = hf.get('sigma2')[()]
        tau = hf.get('tau')[()]

        chromosome = np.append(chromosome, chromosome_file, axis = 0)
        snp = np.append(snp, snp_file, axis = 0)
        theta = np.append(theta, theta_file, axis = 0)
        se = np.append(se, se_file, axis = 0)
        S = np.append(S, S_file, axis = 0)
        f = np.append(f, f_file, axis = 0)
        N = np.append(N, N_file, axis = 0)

    # Constructing dataframe of data
    zdata = pd.DataFrame({'CHR' : chromosome,
                        'BP' : snp,
                        'N' : N,
                        "f" : f,
                        'theta' : theta.tolist(),
                        'se' : se.tolist(),
                        "S" : S.tolist()})


    zdata['CHR'] = zdata['CHR'].astype(int)
    zdata['BP'] = zdata['BP'].astype(str).str.replace("b'", "").str[:-1]

    zdata['BP'] = zdata['BP'].astype('int')

    main_df = zdata.merge(ldscores, how = "inner", on = ["CHR", "BP"])

    # dropping NAs
    main_df = main_df.dropna()


    # transforming inputs

    S = np.array(list(main_df.S)) 
    theta = np.array(list(main_df.theta))
    f = np.array(list(main_df["f"]))
    l = np.array(list(main_df["L2"]))
    u = np.array(list(main_df["L2"]))

    # M = mfile['M'].sum()
    M = len(S)

    effect_estimated = "direct_plus_population"

    S, theta = ld.transform_estimates(effect_estimated, S, theta)


    # making z value
    zval = ld.theta2z(theta, S, M = M)

    # == Initializing model == #
    model = ld.sibreg(S = S, 
                    z = zval, 
                    l = l,
                    f = f,
                    u = u,
                    M = M) 


    output_matrix, result = model.solve() 

    logging.info(f"Output Matrix: {output_matrix}")
    logging.info(f"Result: {result}")

    estimationTime = (datetime.datetime.now() - startTime)
    logging.info(f"Estimation time (before calculating standard errors): {estimationTime}")
    

# # Get standard errors
# stderrors = model.jackknife_se(blocksize = 1000)

# logging.info(f"Jackknife Standard Errors: {stderrors}")

# executionTime = (datetime.datetime.now() - startTime)
# logging.info(f"Execution time: {executionTime}")