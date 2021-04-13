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

sys.path.append("/homes/nber/harij/gitrepos/SNIPar/ldsc_reg")
import sib_ldsc_z as ld

import scipy.stats

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
    # parser.add_argument('-ldsc', '--ldsc_scores', type = str, required = True, help = '''Directory of LDSC scores.''')
    parser.add_argument('-m', '--mfiles', type = str, help = '''Directory of M files scores. If left blank M will be
                                                                the number of observations''')
    parser.add_argument('-l', '--logfile', type=str, help = '''The log file where the results will be stored.
                                                        If left blank no log file will be saved.''')
    parser.add_argument('-e', '--effect_transform', type = str, help ='''
                                                            How to convert the 3 dimensional data
                                                            into 2 dimensions. Options are direct_plus_population,
                                                            direct_plus_averageparental, full and population. 
                                                            Default is direct_plus_population.''')
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
    # parser.add_argument('--rbound-jkse', dest = "rbounds_jkse", action='store_true')
    # parser.add_argument('--no-rbound-jkse', dest = "rbounds_jkse", action = 'store_false')

    parser.add_argument('-maf', '--maf-thresh', dest = "maf", type = float, help = """The threshold of minor allele frequency. All SNPs below this threshold
                                                                        are dropped. Default is 5 percent. Set number as the percentage i.e. 5 instead of 0.05""")
    parser.set_defaults(maf = 5.0)
    parser.add_argument('--print-delete-values', dest = 'store_delvals', action = 'store_true',
    help = 'Option to save jack knife delete values. Saves to same location as log file.')

    # names of variable names
    parser.add_argument('-bim', default = "bim", type = str, help = "Name of bim column")
    parser.add_argument('-bim_chromosome', default = 0, type = int, help = "Column index of Chromosome in BIM variable")
    parser.add_argument('-bim_rsid', default = 1, type = int, help = "Column index of SNPID (RSID) in BIM variable")

    parser.add_argument('--rsid_readfrombim', type = str, 
                        help = '''Needs to be a comma seperated string of filename, BP-position, SNP-position, seperator.
                        If provided the variable bim_snp wont be used instead rsid's will be read
                        from the provided file set.''')
    parser.add_argument('-bim_bp', default = 3, type = int, help = "Column index of BP in BIM variable")
    parser.add_argument('-bim_a1', default = 4, type = int, help = "Column index of Chromosome in A1 variable")
    parser.add_argument('-bim_a2', default = 5, type = int, help = "Column index of Chromosome in A2 variable")

    parser.add_argument('-estimate', default = "estimate", type = str, help = "Name of estimate column")
    parser.add_argument('-estimate_ses', default = "estimate_ses", type = str, help = "Name of estimate_ses column")
    parser.add_argument('-Nname', default = "N_L", type = str, help = "Name of N column in hdf5 files.")
    parser.add_argument('-N', type = int, help = "GWAS N. One number used for all SNPs.")
    parser.add_argument('--Neff', action = 'store_true', help = "For N calculate effective N instead of actual GWAS N.")
    parser.add_argument('-estimate_covariance', default = "estimate_covariance", type = str, help = "Name of estimate_covariance column")
    parser.add_argument('-freqs', default = "freqs", type = str, help = "Name of freqs column")

    parser.add_argument('-sigma2', default = "sigma2", type = str, help = "Name of sigma2 column")
    parser.add_argument('-tau', default = "tau", type = str, help = "Name of tau column")


    parser.add_argument('-outsumstat', type = int, default = 1, help = "Column from estimates to output as the beta column")
    parser.add_argument('-merge-alleles', type = str, help = "Alleles to merge and keep")
    parser.add_argument('-out', type = str, help = "outputdir")

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


    # rsid files
    if args.rsid_readfrombim is not None:
        rsid_parts = args.rsid_readfrombim.split(",")
        rsidfiles = rsid_parts[0]
        bppos = int(rsid_parts[1])
        rsidpos = int(rsid_parts[2])
        file_sep = str(rsid_parts[3])
    
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
    metadata = hf.get(args.bim)[()]
    chromosome = metadata[:, args.bim_chromosome]
    bp = metadata[:, args.bim_bp]
    if args.rsid_readfrombim is not None:
        snp = np.zeros(bp.shape[0])
    else:
        snp = metadata[:, args.bim_rsid]
    A1 = metadata[:, args.bim_a1]
    A2 = metadata[:, args.bim_a2]
    theta  = hf.get(args.estimate)[()]
    se  = hf.get(args.estimate_ses)[()]
    S = hf.get(args.estimate_covariance)[()]
    f = hf.get(args.freqs)[()]
    if args.N is None and args.Neff is None:
        N = hf.get(args.Nname)[()]

    sigma2 = hf.get(args.sigma2)[()]
    tau = hf.get(args.tau)[()]
    phvar = sigma2+sigma2/tau
    
    hf.close()

    if len(files) > 1:
        for file in files[1:]:
            print("Reading in file: ", file)
            hf = h5py.File(file, 'r')
            metadata = hf.get(args.bim)[()]
            chromosome_file = metadata[:, args.bim_chromosome]
            bp_file = metadata[:, args.bim_bp]
            if args.rsid_readfrombim:
                snp_file = np.zeros(bp_file.shape[0])
            else:
                snp_file = metadata[:, args.bim_rsid]
            A1_file = metadata[:, args.bim_a1]
            A2_file = metadata[:, args.bim_a2]
            theta_file  = hf.get(args.estimate)[()]
            se_file  = hf.get(args.estimate_ses)[()]
            S_file = hf.get(args.estimate_covariance)[()]
            f_file = hf.get(args.freqs)[()]
            if args.N is None and args.Neff is None:
                N_file = hf.get(args.Nname)[()]

            chromosome = np.append(chromosome, chromosome_file, axis = 0)
            snp = np.append(snp, snp_file, axis = 0)
            bp = np.append(bp, bp_file, axis = 0)
            A1 = np.append(A1, A1_file, axis = 0)
            A2 = np.append(A2, A2_file, axis = 0)
            theta = np.append(theta, theta_file, axis = 0)
            se = np.append(se, se_file, axis = 0)
            S = np.append(S, S_file, axis = 0)
            f = np.append(f, f_file, axis = 0)
            if args.N is None and args.Neff is None:
                N = np.append(N, N_file, axis = 0)
            hf.close()

    if args.N is not None and args.Neff is None:
        N = np.repeat(args.N, theta.shape[0])
    elif args.Neff:
        # create placeholder N till we calculate population or average parental effects
        N = np.zeros(theta.shape[0])


    # Constructing dataframe of data
    zdata = pd.DataFrame({'CHR' : chromosome.astype(int),
                        'SNP' : snp.astype(str),
                        'BP' : bp.astype(int),
                        "f" : f,
                        "N" : N,
                        "A1" : A1.astype(str),
                        "A2" : A2.astype(str),
                        'theta' : theta.tolist(),
                        'se' : se.tolist(),
                        "S" : S.tolist()})
    

    if args.rsid_readfrombim is not None:
        rsidfiles = glob.glob(rsidfiles)
        snps = pd.DataFrame(columns = ["BP", "rsid"])
        for file in rsidfiles:
            snp_i = ld.return_rsid(file, bppos, rsidpos, file_sep)
            snps = snps.append(snp_i, ignore_index = True)
        
        snps = snps.drop_duplicates(subset=['BP'])
        zdata = zdata.merge(snps, how = "left", on = "BP")
        zdata = zdata.rename(columns = {"SNP" : "SNP_old"})
        zdata = zdata.rename(columns = {"rsid" : "SNP"})

    zdata_n_message = f"Number of Observations before merging LD-Scores, before removing low MAF SNPs: {zdata.shape[0]}"
    
    print(zdata_n_message)
    if args.logfile is not None:
        logging.info(zdata_n_message)
    
    # dropping obs based on MAF
    zdata = zdata[zdata['f'] >= args.maf/100.0]
    zdata = zdata[zdata['f'] <= 1-(args.maf/100.0)]

    zdata_n_message = f"Number of Observations before merging LD-Scores, after removing low MAF SNPs: {zdata.shape[0]}"
    print(zdata_n_message)

    zdata = zdata.sort_values(by=['CHR', 'BP'])


    S = np.array(list(zdata.S)) 
    theta = np.array(list(zdata.theta))
    f = np.array(list(zdata["f"]))
    chromosomes = np.array(list(zdata["CHR"]))
    snp = np.array(list(zdata["SNP"]))
    bp = np.array(list(zdata["BP"]))
    a1 = np.array(list(zdata["A1"]))
    a2 = np.array(list(zdata["A2"]))
    N = np.array(list(zdata["N"]))


    if args.effect_transform is not None:
        effect_estimated = args.effect_transform
    else:
        effect_estimated = "direct_plus_population"
    
    effect_message = f"Transforming effects into: {effect_estimated}"

    print(effect_message)
    if args.logfile is not None:
        logging.info(effect_message)

    S, theta = ld.transform_estimates(effect_estimated, S, theta)

    if args.Neff:
        N = np.round(phvar*(2*f*(1-f)*S[..., args.outsumstat, args.outsumstat])**(-1)).astype('int')

    if args.mfiles is not None:
        Mfiles = args.mfiles
        Mcolnames = ["M", "CHR"]
        nloci = ld.read_mfiles(Mfiles, Mcolnames)
        M = nloci['M'].sum()
    else:
        M = len(S)

    # making z value
    zval = ld.theta2z(theta, S, M = M)
    PVAL = lambda z: 2*scipy.stats.norm.sf(np.abs(z))
    pval = PVAL(zval)
    
    # making output
    dfout = pd.DataFrame(
        {
            'chr' : chromosomes,
            'SNP' : snp,
            'pos' : bp,
            'A1' : a1,
            'A2' : a2,
            'N' : N,
            'P' : pval[..., args.outsumstat],
            'beta' : theta[..., args.outsumstat],
            'Z' : zval[..., args.outsumstat],
            'SE' : np.sqrt(S[..., args.outsumstat, args.outsumstat]),
            'EAF' : f,
            'MAF' : 1.0 - f
        }
    )

    if args.merge_alleles is not None:
        print(f"N obs before merging alleles: {dfout.shape[0]}")
        ii = np.array([True for i in range(len(dfout))])
        old = ii.sum()
        merge_alleles = pd.read_csv(args.merge_alleles, sep = "\t")
        ii = dfout.SNP.isin(merge_alleles.SNP)
        drops = old - ii.sum()
        if ii.sum() == 0:
            pass
        dfout = dfout[ii].reset_index(drop=True)
        print(f"N obs after merging alleles: {dfout.shape[0]}")


    dfout = dfout.dropna()
    new_obs = dfout.shape[0]

    print(f"Number of observations dropped due to NAs: {zdata.shape[0] - new_obs}")

    dfout.to_csv(
        args.out,
        sep = " ",
        index = False
    )

    print("Done outputting sumstats file!")


