'''
This script reads in the simulated data made by
Alex and solves for V.
'''
import sib_ldsc_z as ld
import numpy as np
import h5py
import glob
import time
import matplotlib.pyplot as plt
import pandas as pd

startTime = time.time()
print("Start time: ", startTime)

# == Direct Effect == #
print("=====================================")
print("Making CSV for Average Parental Effects")
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
theta  = hf.get('estimate')[()]
se  = hf.get('estimate_ses')[()]
N = hf.get('N_L')[()]
S = hf.get('estimate_covariance')[()]
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
    theta_file  = hf.get('estimate')[()]
    se_file  = hf.get('estimate_ses')[()]
    S_file = hf.get('estimate_covariance')[()]
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

# Constructing dataframe of data
zdata = pd.DataFrame({'CHR' : chromosome,
                    'SNP' : snp,
                    'pos' : pos,
                    'A1' : A1,
                    'A2' : A2,
                    'N' : N,
                    'z' : zval.tolist(),
                    'se' : se.tolist(),
                    "S" : S.tolist()})


zdata['CHR'] = zdata['CHR'].astype(float)
zdata['SNP'] = zdata['SNP'].astype(str).str.replace("b'", "").str[:-1]
zdata['pos'] = zdata['pos'].astype(str).str.replace("b'", "").str[:-1]
zdata['A1'] = zdata['A1'].astype(str).str.replace("b'", "").str[:-1]
zdata['A2'] = zdata['A2'].astype(str).str.replace("b'", "").str[:-1]

# == Reading in LD Scores == #
ldscore_path = "/var/genetics/pub/data/ld_ref_panel/eur_w_ld_chr/"

def M(fh, num=None, N=2, common=False):
    '''Parses .l{N}.M files, split across num chromosomes. See docs/file_formats_ld.txt.'''
    parsefunc = lambda y: [float(z) for z in open(y, 'r').readline().split()]
    suffix = '.l' + str(N) + '.M'
    suffix += '_5_50'

    if num is not None:
        x = np.sum([parsefunc(sub_chr(fh, i) + suffix) for i in get_present_chrs(fh, num+1)], axis=0)
    else:
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
