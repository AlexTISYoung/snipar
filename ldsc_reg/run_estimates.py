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
import sib_ldsc_z as ld

if __name__ == '__main__':
    # command line arguments
    parser=argparse.ArgumentParser()
    parser.add_argument('path2file', type=str, 
                        help='''Path to hdf5 file with summary statistics. 
                                                Should include the file name.
                                                In case of multiple files include glob character.
                                                eg: /path/to/file/chr_*.hdf5''')
    parser.add_argument('-ldsc', '--ldsc_scores', type = str, required = True, help = '''Directory of LDSC scores.''')
    parser.add_argument('-m', '--mfiles', type = str, help = '''Directory of M files scores. If left blank M will be
                                                                the number of observations''')
    parser.add_argument('-l', '--logfile', type=str, help = '''The log file where the results will be stored.
                                                        If left blank no log file will be saved.''')
    parser.add_argument('-e', '--effect_transform', type = str, help ='''
                                                            How to convert the 3 dimensional data
                                                            into 2 dimensions. Options are direct_plus_population and
                                                            direct_plus_averageparental. Default is direct_plus_population.''')
    parser.add_argument('--rbound', dest = "rbounds", action='store_true')
    parser.add_argument('--no-rbound', dest = "rbounds", action = 'store_false')
    parser.set_defaults(rbounds=True)
    
    args=parser.parse_args()
    
    if args.logfile is not None:
        logging.basicConfig(filename= args.logfile, 
                                 level = logging.INFO,
                                 format = "%(message)s",
                                 filemode = "w")
    
    startTime = datetime.datetime.now()
    
    print(f"===============================")
    print(f"Start time:  {startTime}")
    
    if args.logfile is not None:
        logging.info(f"Start time:  {startTime}")

    # == Reading in data == #
    print("Reading in Data")
    # reading in  data
    filenames = args.path2file

    files = glob.glob(filenames)

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

    if len(files) > 1:
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
    
    zdata_n_message = f"Number of Observations before merging LD-Scores: {zdata.shape[0]}"
    
    print(zdata_n_message)
    if args.logfile is not None:
        logging.info(zdata_n_message)


    # == Reading in LD Scores == #
    ldscore_path = args.ldsc_scores
    ldcolnames = ["CHR", "SNP", "BP", "L2"]
    ldscores= ld.read_ldscores(ldscore_path, ldcolnames)

    # Merging LD scores with main Data Frame
    main_df = zdata.merge(ldscores, how = "inner", on = ["CHR", "SNP"])

    # dropping NAs
    main_df = main_df.dropna()
    
    maindata_n_message = f"Number of Observations after merging LD-Scores and dropping NAs: {main_df.shape[0]}"
    
    print(maindata_n_message)
    if args.logfile is not None:
        logging.info(maindata_n_message)

    # transforming inputs

    S = np.array(list(main_df.S)) 
    theta = np.array(list(main_df.theta))
    f = np.array(list(main_df["f"]))
    l = np.array(list(main_df["L2"]))
    u = np.array(list(main_df["L2"]))
    
    if args.mfiles is not None:
        Mfiles = args.mfiles
        Mcolnames = ["M", "CHR"]
        nloci = ld.read_mfiles(Mfilepath, Mcolnames)
        M = nloci['M'].sum()
    else:
        M = len(S)
    
    if args.effect_transform is not None:
        effect_estimated = args.effect_transform
    else:
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

    output_matrix, result = model.solve(rbounds = args.rbounds)
    
    estimates = {'v1' : output_matrix['v1'],
                'v2' : output_matrix['v2'],
                'r' : output_matrix['r']}
    
    std_errors = np.diag(output_matrix['std_err_mat'])
    
    estimationTime = (datetime.datetime.now() - startTime)
    
    print("----------------------------------")
    print(f"Estimates: {estimates}")
    print(f"Standard Errors: {std_errors}")
    print(f"Maximizer Output: {result}")
    print(f"Estimation time (before calculating standard errors): {estimationTime}")
    
    if args.logfile is not None:
        logging.info("----------------------------------")
        logging.info(f"Estimates: {estimates}")
        logging.info(f"Standard Errors: {std_errors}")
        logging.info(f"Maximizer Output: {result}")
        logging.info(f"Estimation time (before calculating standard errors): {estimationTime}")