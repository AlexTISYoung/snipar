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

def print_call(call):
    
    '''
    Gives the call given to python
    in a nice string
    '''
    
    message = ''
    for i in range(len(call)):
        if call[i][0] != "-":
            message += call[i]
            message += ' \\ \n'
        else:
            message += call[i] + " "
    
    return message[1:-1]
           
        

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
    
    parser.add_argument('--jkse', dest = "jkse", action = 'store_true', help = '''
    Specifies if one wants to estimate block Jack Knife Standard errors
    ''')
    parser.set_defaults(jkse=False)
    parser.add_argument('--jkse_nblocks', type = int,
                    help = "Number of blocks for block jack knife SE estimation.")
    parser.add_argument('--jkse_blocksize', type = int, help = "Block Size for Block Jackknife Standard Errors.")
    parser.add_argument('--jkse_cores', type = int, help = "Number of cores to use for block Jack Knife standard errors.")
    parser.set_defaults(jkse_cores = 2)
    
    parser.add_argument('-maf', '--maf-thresh', dest = "maf", type = float, help = """The threshold of minor allele frequency. All SNPs below this threshold
                                                                        are dropped. Default is 5 percent. Set number as the percentage i.e. 5 instead of 0.05""")
    parser.set_defaults(maf = 5.0)
    args=parser.parse_args()
    
    
    if args.jkse_blocksize is not None and args.jkse == False:
        print('''Option for Block Jack Knife block size was passed but wasn't told to actually estimate
                Block Jackknife Standard Errors. Script will run but will not estimate Block Jack Knife Standard
                Errors.''')
    
    
    if args.jkse_cores is not None and args.jkse == False:
        print('''Option for Block Jack Knife cores was passed but wasn't told to actually estimate
                Block Jackknife Standard Errors. Script will run but will not estimate Block Jack Knife Standard
                Errors.''')
    
    if args.logfile is not None:
        logging.basicConfig(filename= args.logfile, 
                                 level = logging.INFO,
                                 format = "%(message)s",
                                 filemode = "w")
        
    args_call = print_call(sys.argv)
    print(args_call)
    if args.logfile is not None:
        logging.info(f"Call: \n {args_call}")
    
    
    
    startTime = datetime.datetime.now()
    
    print(f"===============================")
    print(f"Start time:  {startTime}")
    
    if args.logfile is not None:
        logging.info(f"Start time:  {startTime}")

    # == Reading in data == #
    print("Reading in Data")
    # reading in  data
    filenames = args.path2file

    zdata = pd.DataFrame(
        columns = ['CHR', 'SNP', 'BP' ,'N',
                    "f" ,"A1", "A2" , 'theta', 
                    'se', "S" ]
    )
    for chromosome in range(1, 23):
        file = f"{filenames}/chr_{chromosome}.hdf5"
        print("Reading in file: ", file)
        hf = h5py.File(file, 'r')
        metadata = hf.get("bim")[()]
        chromosome = metadata[:, 0]
        snp = metadata[:, 1]
        bp = metadata[:, 3]
        A1 = metadata[:, 4]
        A2 = metadata[:, 5]
        theta  = hf.get('estimate')[()]
        se  = hf.get('estimate_ses')[()]
        N = hf.get('N_L')[()]
        S = hf.get('estimate_covariance')[()]
        f = hf.get('freqs')[()]

        # Constructing dataframe of data
        zdata_toappend = pd.DataFrame({'CHR' : chromosome,
                'SNP' : snp,
                'BP' : bp,
                'N' : N,
                "f" : f,
                "A1" : A1,
                "A2" : A2,
                'theta' : theta.tolist(),
                'se' : se.tolist(),
                "S" : S.tolist()})

        zdata = zdata.append(zdata_toappend)
    
    sigma2 = hf.get('sigma2')[()]
    tau = hf.get('tau')[()]
    phvar = sigma2+sigma2/tau

    
    # cleaning data
    zdata['CHR'] = zdata['CHR'].astype(int)
    zdata['SNP'] = zdata['SNP'].astype(str).str.replace("b'", "").str[:-1]
    zdata['BP'] = zdata['BP'].astype(str).str.replace("b'", "").str[:-1]
    zdata['BP'] = zdata['BP'].astype('int')
    
    # sorting by chromosome
    zdata = zdata.sort_values(by = ['CHR']).reset_index(drop = True)
    zdata['ordering'] = zdata.index
    
    
    zdata_n_message = f"Number of Observations before merging LD-Scores, before removing low MAF SNPs: {zdata.shape[0]}"
    
    print(zdata_n_message)
    if args.logfile is not None:
        logging.info(zdata_n_message)
    
    # dropping obs based on MAF
    zdata = zdata[zdata['f'] >= args.maf/100.0]
    
    zdata_n_message = f"Number of Observations before merging LD-Scores, after removing low MAF SNPs: {zdata.shape[0]}"
    
    print(zdata_n_message)
    if args.logfile is not None:
        logging.info(zdata_n_message)


    # == Reading in LD Scores == #
    ldscore_path = args.ldsc_scores
    ldcolnames = ["CHR", "SNP", "BP", "L2"]
    ldscores= ld.read_ldscores(ldscore_path, ldcolnames)
    ldscores['BP'] = ldscores['BP'].astype('int')

    # dropping NAs
    main_df = zdata.merge(ldscores, how = "inner", on = ["CHR", "SNP"])

    if main_df.shape[0] > 0:
        assert np.all(main_df.BP_x == main_df.BP_y)
        main_df = main_df.drop('BP_y', axis = 1)
        main_df = main_df.rename(columns = {'BP_x' : 'BP'})
        main_df = main_df.dropna()
    elif main_df.shape[0] == 0:
        print("No matches while matching LD-score data with main dataset using RSID.")
        print("Trying to match with BP.")
        main_df = zdata.merge(ldscores, how = "inner", on = ["CHR", "BP"])
        main_df = main_df.drop('SNP_y', axis = 1)
        main_df = main_df.rename(columns = {'SNP_x' : 'SNP'})
        main_df = main_df.dropna()

    maindata_n_message = f"Number of Observations after merging LD-Scores and dropping NAs: {main_df.shape[0]}"
    main_df = main_df.sort_values(by=['CHR', 'BP'])

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
    
    effect_message = f"Transforming effects into: {effect_estimated}"

    print(effect_message)
    if args.logfile is not None:
        logging.info(effect_message)

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

    output_matrix, result = model.solve(est_init = np.array([phvar/2, phvar/2, 0.0]),
                                        rbounds = args.rbounds)
    
    estimates = {'v1' : output_matrix['v1'],
                'v2' : output_matrix['v2'],
                'r' : output_matrix['r']}
    
    std_errors = np.diag(output_matrix['std_err_mat'])
    
    estimationTime = (datetime.datetime.now() - startTime)
    
    print("----------------------------------")
    print(f"Estimates: {estimates}")
    print(f"Standard Errors: {std_errors}")
    print(f"Maximizer Output: {result}")
    print(f"Estimation time: {estimationTime}")
    
    if args.logfile is not None:
        logging.info("----------------------------------")
        logging.info(f"Estimates: {estimates}")
        logging.info(f"Standard Errors: {std_errors}")
        logging.info(f"Maximizer Output: {result}")
        logging.info(f"Estimation time: {estimationTime}")
        
    if args.jkse:
        
        if args.jkse_nblocks is not None:
            nblocks_blocksize = np.ceil(zval.shape[0]/args.jkse_nblocks)
            blocksize = int(nblocks_blocksize)
        elif args.jkse_blocksize is not None and args.jkse_nblocks is None:
            # need nblocks to always override blocksize option
            blocksize = int(args.jkse_blocksize)
        else:
            nblocks = 200
            nblocks_blocksize = np.ceil(zval.shape[0]/nblocks)
            blocksize = int(nblocks_blocksize)

            

        print(f"Jack Knife Block Sizes = {blocksize}")
        print(f"Number of cores being used for Jack Knife: {args.jkse_cores}")
        print("Estimating Block Jackknife Standard Errors...")
        
        initguess = {'v1' : phvar/2, 'v2' : phvar/2, 'r' : 0.0} #output_matrix
        jkse, delvals = ld.jkse(model, initguess, blocksize = blocksize, num_procs=args.jkse_cores,
                        rbounds = args.rbounds)

        print(f"Block Jack Knife Standard Errors: {jkse}")
        
        estimationTime_jkse = (datetime.datetime.now() - startTime)
        
        print(f"Estimation time with Block Jack Knife Standard Error Estimation: {estimationTime_jkse}")
        

        if args.logfile is not None:
            logging.info(f"Jack Knife Block Sizes = {blocksize}")
            logging.info(f"Number of cores being used for Jack Knife: {args.jkse_cores}")
            logging.info(f"Block Jack Knife Standard Errors: {jkse}")
            logging.info(f"Estimation time with Block Jack Knife Standard Error Estimation: {estimationTime_jkse}")
            
            