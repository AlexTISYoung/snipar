from msilib import make_id
import numpy as np
import h5py
from snipar.utilities import *

def pair_to_lower_tri(index1,index2):
    if index1>index2:
        return int(index1*(index1+1)/2)+index2
    else:
        return int(index2*(index2+1)/2)+index1

dir = '/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/grms/'

r = np.fromfile(dir+'/R.grm.bin',dtype=np.float32)
ids = np.loadtxt(dir+'/R.grm.id',dtype='<U20')
ped_ids = np.array([ids[x,0]+'_'+ids[x,1] for x in range(ids.shape[0])])
ids[:,1] = ped_ids
ids[:,0] = ped_ids

diag_indices = np.array([int(x*(x+1)/2)-1 for x in range(1,ids.shape[0]+1)])

r_diag = r[diag_indices]

# Make upper
r_upper = np.zeros(r.shape,dtype=r.dtype)
r_upper[:] = r 
r_upper[r_upper<0.05] = 0
r_upper.tofile(dir+'/R_upper.grm.bin')
np.savetxt(dir+'/R_upper.grm.id',ids,fmt='%s')

# Make lower
r_lower = np.zeros(r.shape,r.dtype)
r_lower[r<0.05] = r[r<0.05]
r_lower[diag_indices] = r_diag
r_lower.tofile(dir+'/R_lower.grm.bin')
np.savetxt(dir+'/R_lower.grm.id',ids,fmt='%s')

## Make sib
# read pedigree
imp = h5py.File(dir+'/../imputed/chr_1.hdf5','r')
ped = convert_str_array(imp['pedigree'])
ped = ped[1:ped.shape[0],:]
ped = ped[~np.array([x[0]=='_' for x in ped[:,0]]),:]
ped_dict = make_id_dict(ped,1)

# make binary grm
npair = int(ped.shape[0]*(ped.shape[0]+1)/2)
r_sib = np.zeros((npair),dtype=np.float32)

# fill diag
sib_diag = np.array([int(x*(x+1)/2)-1 for x in range(1,ped.shape[0]+1)])
r_sib[sib_diag] = 1

# fill sib
fams = np.unique(ped[:,0],return_counts=True)
f2 = fams[0][fams[1]>=2]
for fam in f2:
    fam_ids = ped[ped[:,0]==fam,1]
    fam_indices = np.array([ped_dict[x] for x in fam_ids])
    for i in range(fam_indices.shape[0]):
        for j in range(i):
            r_sib[pair_to_lower_tri(fam_indices[i],fam_indices[j])] = 1

# Write to file
r_sib_id = np.column_stack((ped[:,1],ped[:,1]))
r_sib.tofile(dir+'/R_sib.grm.bin')
np.savetxt(dir+'/R_sib.grm.id',r_sib_id ,fmt='%s')

# mgrm
mgrm = np.array([dir+'R_lower',dir+'R_upper',dir+'R_sib'])
np.savetxt(dir+'mgrm.txt',mgrm,fmt='%s')

