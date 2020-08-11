import ldsc_72 as ld
import numpy as np
import h5py


data = '/disk/genetics/ukb/alextisyoung/vcinf/1/chr_22.hdf5'
hf = h5py.File(data, 'r')
theta  = hf.get('estimate')[()]
S = hf.get('estimate_covariance')[()]
N = theta.shape[0]


#N = 100
#S = np.array([np.array([[5, 0], [0, 5]]),
#    np.array([[2, 0], [0, 2]])] * 50 )# 50 = N/2
#V = np.identity(2) * 10.0

model = ld.sibreg_72(S = S, theta = theta)

output_matrix, result = model.solve()

print("Output matrix: ",output_matrix)
print("Solver Output: ", result)