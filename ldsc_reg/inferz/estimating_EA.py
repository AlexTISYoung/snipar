'''
This script reads the HDF5 files for EA
and estimates V
'''
import sib_ldsc_z as ld
import numpy as np
import h5py
import glob
import time
import matplotlib.pyplot as plt
import pandas as pd

startTime = time.time()
print("Start time: ", time.localtime(startTime))

# == Direct Effect == #
print("=====================================")
print("Making CSV for Average Parental Effects")
print("=====================================")
# reading in  data
files = glob.glob("/disk/genetics/ukb/alextisyoung/genotypes/EA/chr_*.hdf5")

file = files[0]
print("Reading in file: ", file)
hf = h5py.File(file, 'r')
snp = hf.get("sid")[()]
nobs = len(snp)
chromosome = [int(i) for i in file[-7:-5].split("_") if i.isdigit()] * nobs
theta  = hf.get('estimate')[()]
se  = hf.get('estimate_ses')[()]
N = hf.get('N_L')[()]
S = hf.get('estimate_covariance')[()]
f = hf.get('freqs')[()]

for file in files[1:]:
    print("Reading in file: ", file)
    hf = h5py.File(file, 'r')
    snp_file = hf.get("sid")[()]
    nobs = len(snp_file)
    chromosome_file = [int(i) for i in file[-7:-5].split("_") if i.isdigit()] * nobs    
    theta_file  = hf.get('estimate')[()]
    se_file  = hf.get('estimate_ses')[()]
    S_file = hf.get('estimate_covariance')[()]
    f_file = hf.get('freqs')[()]
    N_file = hf.get('N_L')[()]

    chromosome = np.append(chromosome, chromosome_file, axis = 0)
    snp = np.append(snp, snp_file, axis = 0)
    theta = np.append(theta, theta_file, axis = 0)
    se = np.append(se, se_file, axis = 0)
    S = np.append(S, S_file, axis = 0)
    f = np.append(f, f_file, axis = 0)
    N = np.append(N, N_file, axis = 0)

# getting p value
zval = theta/se

# Constructing dataframe of data
zdata = pd.DataFrame({'CHR' : chromosome,
                    'SNP' : snp,
                    'N' : N,
                    'z' : zval.tolist(),
                    'se' : se.tolist(),
                    "S" : S.tolist()})


zdata['CHR'] = zdata['CHR'].astype(float)
zdata['SNP'] = zdata['SNP'].astype(str).str.replace("b'", "").str[:-1]


# == Reading in LD Scores == #
ldscore_path = "/var/genetics/pub/data/ld_ref_panel/eur_w_ld_chr/"

def M(fh, N=2, common=False):
    '''Parses .l{N}.M files, split across num chromosomes. See docs/file_formats_ld.txt.'''
    parsefunc = lambda y: [float(z) for z in open(y, 'r').readline().split()]
    suffix = '.l' + str(N) + '.M'
    suffix += '_5_50'

    x = parsefunc(fh + suffix)

    return np.array(x).reshape((1, len(x)))

# Reading LD Scores

files = glob.glob(f"{ldscore_path}/*[0-9].l2.ldscore.gz")
ldscores = pd.DataFrame(columns = ["CHR", "SNP", "BP", "CM", "MAF", "L2"])

for file in files:
    snpi = pd.read_csv(file, compression='gzip', sep = "\t")
    ldscores = pd.concat([ldscores, snpi])

# Reading Number of Loci
files_M = glob.glob(f"{ldscore_path}/*[0-9].l2.M_5_50")
nloci = pd.DataFrame(columns = ["M", "CHR"])

for file in files_M:
    
    if len(file) == 60:
        chrom = int(file[-11])
    else:
        chrom = int(file[-12:-10])
    
    nloci_i = pd.DataFrame({"M" : [M(file[:-10])[0, 0]],
                             "CHR" : [chrom]})
    
    nloci = pd.concat([nloci, nloci_i])

nloci = nloci.reset_index(drop = True)

ldscores = ldscores.merge(nloci, how = "left", on = "CHR")

# Merging LD scores with main Data Frame
main_df = zdata.merge(ldscores, how = "inner", on = ["CHR", "SNP"])

# dropping NAs
main_df = main_df.dropna()

# transforming inputs

S = np.array(list(main_df.S)) 
z = np.array(list(main_df.z))
f = np.array(list(main_df["MAF"]))
r = np.array(list(main_df["L2"]))
u = np.array(list(main_df["L2"]))

effect_estimated = "direct_plus_population"

if effect_estimated == "population":
    # == Keeping population effect == #
    Sdir = np.empty(len(S))
    for i in range(len(S)):
        Sdir[i] = np.array([[1.0, 0.5, 0.5]]) @ S[i] @ np.array([[1.0, 0.5, 0.5]]).T

    S = Sdir.reshape((len(S), 1, 1))
    z = z @ np.array([1.0, 0.5, 0.5])
    z = z.reshape((z.shape[0], 1))
elif effect_estimated == "direct_plus_averageparental":

    # == Combining indirect effects to make V a 2x2 matrix == #
    tmatrix = np.array([[1.0, 0.0],
                        [0.0, 0.5],
                        [0.0, 0.5]])
    Sdir = np.empty((len(S), 2, 2))
    for i in range(len(S)):
        Sdir[i] = tmatrix.T @ S[i] @ tmatrix
    S = Sdir.reshape((len(S), 2, 2))
    z = z @ tmatrix
    z = z.reshape((z.shape[0], 2))
elif effect_estimated == "direct_plus_population":

    # == keeping direct effect and population effect == #
    tmatrix = np.array([[1.0, 1.0],
                        [0.0, 0.5],
                        [0.0, 0.5]])
    Sdir = np.empty((len(S), 2, 2))
    for i in range(len(S)):
        Sdir[i] = tmatrix.T @ S[i] @ tmatrix

    S = Sdir.reshape((len(S), 2, 2))
    z = z @ tmatrix
    z = z.reshape((z.shape[0], 2))
elif effect_estimated == "full":
    pass

# == Initializing model == #
model = ld.sibreg(S = S, 
                z = z, 
                f = f,
                r = r,
                u = u) 

output_matrix, result = model.solve()

print(f"======================================")
print(f"Output Matrix: {output_matrix}")
print(f"Result: {result}")


executionTime = (time.time() - startTime)
print('Execution time: ' + f'{executionTime:.2f}', " seconds")
