'''
This script reads in the simulated data made by
Alex and solves for V. Currently it is unable to
both with just the causal SNPs and the entire
set of SNPs.
'''
import ldsc_72 as ld
import ldsc_71 as ld1
import numpy as np
import h5py
import glob
import time
startTime = time.time()

files = glob.glob('/disk/genetics/ukb/alextisyoung/vcinf/1/causal.hdf5')
# files = glob.glob("/disk/genetics/ukb/alextisyoung/vcinf/1/chr_*.hdf5")
#files = glob.glob("C:/Users/Hariharan/Documents/genoecon_work/snipardata/causal.hdf5")
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

# Keeping only direct effects
S = S[:,0 ,0].reshape((S.shape[0], 1, 1))
theta = theta[:, 0].reshape((theta.shape[0], 1))

# # amplifying direct effects
# Sdir = np.empty(len(S))
# for i in range(len(S)):
#   Sdir[i] = np.array([[1.0, 0.5, 0.5]]) @ S[i] @ np.array([1.0, 0.5, 0.5]).T

# S = Sdir.reshape((len(S), 1, 1))
# theta = theta @ np.array([1.0, 0.5, 0.5])
# theta = theta.reshape((theta.shape[0], 1))


model = ld.sibreg(S = S, theta = theta, f = f)

print("Solving Model...")


output_matrix, result = model.solve(est_init = np.atleast_2d(0.0))


print("Output matrix: ",output_matrix)
print("Solver Output: ", result)

executionTime = (time.time() - startTime)
print('Execution time: ' + f'{executionTime:.2f}', " seconds")
