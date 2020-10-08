'''
This script reads in the simulated data made by
Alex and solves for V. Currently it is unable to
both with just the causal SNPs and the entire
set of SNPs.
'''
import sib_ldsc_z as ld
import numpy as np
import h5py
import glob
import time
import matplotlib.pyplot as plt


effect_estimated = "direct_plus_averageparental"
# effect estimated can be
# full (for a 3x3 V matrix)
# population
# direct_plus_averageparental
# direct_plus_population

print(f"Estimating {effect_estimated} effect")

startTime = time.time()
print("Start time: ", startTime)

# files = glob.glob('/disk/genetics/ukb/alextisyoung/vcinf/1/causal.hdf5')
# files = glob.glob("/disk/genetics/ukb/alextisyoung/vcinf/1/chr_*.hdf5")
files = glob.glob("C:/Users/Hariharan/Documents/genoecon_work/snipardata/causal.hdf5")
print("Reading files...")

# read in first file
file = files[0]
print("Reading file: ", file)
hf = h5py.File(file, 'r')
theta  = hf.get('estimate')[()]
S = hf.get('estimate_covariance')[()]
f = hf.get('freqs')[()]

if len(files) > 1:
    for file in files[1:]:
        print("Reading file: ", file)
        hf = h5py.File(file, 'r')
        theta_file  = hf.get('estimate')[()] 
        S_file = hf.get('estimate_covariance')[()]
        f_file = hf.get('freqs')[()]

        theta = np.append(theta, theta_file, axis = 0)
        S = np.append(S, S_file, axis = 0)
        f = np.append(f, f_file, axis = 0)
  

N = theta.shape[0]

print("S matrix:", S)
print("Theta Matrix: ", theta)
print("Initiating Model...")



# == Keeping only direct effects == #
# S = S[:,0 ,0].reshape((S.shape[0], 1, 1))
# theta = theta[:, 0].reshape((theta.shape[0], 1))

if effect_estimated == "population":
    # == Keeping population effect == #
    Sdir = np.empty(len(S))
    for i in range(len(S)):
        Sdir[i] = np.array([[1.0, 0.5, 0.5]]) @ S[i] @ np.array([[1.0, 0.5, 0.5]]).T

    S = Sdir.reshape((len(S), 1, 1))
    theta = theta @ np.array([1.0, 0.5, 0.5])
    theta = theta.reshape((theta.shape[0], 1))
elif effect_estimated == "direct_plus_averageparental":

    # == Combining indirect effects to make V a 2x2 matrix == #
    tmatrix = np.array([[1.0, 0.0],
                        [0.0, 0.5],
                        [0.0, 0.5]])
    Sdir = np.empty((len(S), 2, 2))
    for i in range(len(S)):
        Sdir[i] = tmatrix.T @ S[i] @ tmatrix
    S = Sdir.reshape((len(S), 2, 2))
    theta = theta @ tmatrix
    theta = theta.reshape((theta.shape[0], 2))
elif effect_estimated == "direct_plus_population":

    # == keeping direct effect and population effect == #
    tmatrix = np.array([[1.0, 1.0],
                        [0.0, 0.5],
                        [0.0, 0.5]])
    Sdir = np.empty((len(S), 2, 2))
    for i in range(len(S)):
        Sdir[i] = tmatrix.T @ S[i] @ tmatrix

    S = Sdir.reshape((len(S), 2, 2))
    theta = theta @ tmatrix
    theta = theta.reshape((theta.shape[0], 2))
elif effect_estimated == "full":
    pass

# == calcualting z == #
z = np.empty_like(theta)
z[:] = np.nan
for i in range(z.shape[0]):
    z[i, :] = ld.calc_inv_root(S[i]) @ theta[i, :].T

print("Z: ", z)

model = ld.sibreg(S = S, z = z, f = f) 

print("Solving Model...")


output_matrix, result = model.solve()

print("===================================")
print("Output matrix: ", output_matrix)
print("Solver Output: ", result)

executionTime = (time.time() - startTime)
print('Execution time: ' + f'{executionTime:.2f}', " seconds")



