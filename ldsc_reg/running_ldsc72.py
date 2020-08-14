import ldsc_72 as ld
import ldsc_71 as ld1
import numpy as np
import h5py
import glob
import time
startTime = time.time()

files = glob.glob('/disk/genetics/ukb/alextisyoung/vcinf/1/causal.hdf5')

# read in first file
file = files[0]
hf = h5py.File(file, 'r')
theta  = hf.get('estimate')[()]
S = hf.get('estimate_covariance')[()]
f = hf.get('freqs')[()]

if len(files) > 1:
  for file in files[1:]:
      hf = h5py.File(file, 'r')
      theta_file  = hf.get('estimate')[()]
      S_file = hf.get('estimate_covariance')[()]
      f_file = hf.get('freqs')[()]
      
      theta = np.append(theta, theta_file, axis = 0)
      S = np.append(S, S_file, axis = 0)
      f = np.append(f, f_file, axis = 0)
  
N = theta.shape[0]

# N = 100
# S = np.array([np.array([[5, 0], [0, 5]]),
#       np.array([[2, 0], [0, 2]])] * 50 )# 50 = N/2
# V = np.identity(2) * 10.0

model = ld.sibreg(S = S, theta = theta, f = f)

output_matrix, result = model.solve()

print("Output matrix: ",output_matrix)
print("Solver Output: ", result)

executionTime = (time.time() - startTime)
print('Execution time: ' + str(executionTime), " seconds")