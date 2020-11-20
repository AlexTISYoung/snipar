import numpy as np
import h5py
import glob
import datetime
import matplotlib.pyplot as plt
import pandas as pd
import logging
import argparse

import sys
import os
sys.path.append('/homes/nber/harij/gitrepos/SNIPar/ldsc_reg')
import sib_ldsc_z as ld
    
    
    
logging.basicConfig(filename= f"ldsc_reg/ukb/estimating_ukb.log", 
                            level = logging.INFO,
                            format = "%(message)s",
                            filemode = "w")


traitnos = glob.glob("/disk/genetics/ukb/alextisyoung/haplotypes/relatives/traits/[0-9]*/")
traitnos = [int(i[-3:-1].strip('/')) for i in traitnos]

traitcodes = pd.read_csv("/disk/genetics/ukb/alextisyoung/haplotypes/relatives/traits/traits.txt", 
            names = ['traitcode', 'trait'], 
            sep = ' ')


# == Reading in LD Scores == #
ldscore_path = "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/ldscores/"
ldfiles = "*[0-9].l2.ldscore.gz"
ldcolnames = ["CHR", "SNP", "BP", "L2"]
Mfiles = "*[0-9].l2.M_5_50"
Mcolnames = ["M", "CHR"]

ldscores, mfile = ld.read_ldscores(ldscore_path, ldfiles, ldcolnames, Mfiles, Mcolnames)


for traitcode in traitnos:


    traitname = traitcodes.loc[traitcodes['traitcode'] == traitcode, 'trait'].values[0]
    print(f"Trait name: {traitname}")

    startTime = datetime.datetime.now()  
    logging.info(f"===============================")
    logging.info(f"Trait name:  {traitname}")
    logging.info(f"Start time:  {startTime}")

    # == Reading in data == #
    print("=====================================")
    print("Reading in Data")
    print("=====================================")
    # reading in  data
    files = glob.glob(f"/disk/genetics/ukb/alextisyoung/haplotypes/relatives/traits/{traitcode}/chr_*.hdf5")

    file = files[0]
    print("Reading in file: ", file)
    hf = h5py.File(file, 'r')
    metadata = hf.get("bim")[()]
    chromosome = metadata[:, 0]
    snp = metadata[:, 1]
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
        snp_file = metadata[:, 1]
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
                        'SNP' : snp,
                        'N' : N,
                        "f" : f,
                        'theta' : theta.tolist(),
                        'se' : se.tolist(),
                        "S" : S.tolist()})


    zdata['CHR'] = zdata['CHR'].astype(int)
    zdata['SNP'] = zdata['SNP'].astype(str).str.replace("b'", "").str[:-1]

    # Merging LD scores with main Data Frame
    main_df = zdata.merge(ldscores, how = "inner", on = ["CHR", "SNP"])

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

