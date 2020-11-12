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
import datetime
import matplotlib.pyplot as plt


startTime = datetime.datetime.now() 
print("Start time: ", startTime)

files = glob.glob('/disk/genetics/ukb/alextisyoung/vcinf/1/causal.hdf5')
# files = glob.glob("C:/Users/Hariharan/Documents/genoecon_work/snipardata/causal.hdf5")
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


effect_estimated = "direct_plus_averageparental"
print(f"Estimating {effect_estimated} effect")
S, theta = ld.transform_estimates(effect_estimated, S, theta)

# making z value
zval = ld.theta2z(theta, S, M = len(S))

model = ld.sibreg(S = S, z = z, f = f) 

print("Solving Model...")


output, result = model.solve()

print("===================================")
print("Output matrix: ", output)
print("Solver Output: ", result)

executionTime = (datetime.datetime.now() - startTime)
print('Execution time: ' + f'{executionTime:.2f}', " seconds")



