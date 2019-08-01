import numpy as np
import h5py, argparse
from pysnptools.snpreader import Bed

parser = argparse.ArgumentParser()
parser.add_argument('gts', type=str, help='Path to bed file with sibling genotypes')
parser.add_argument('ped', type=str, help='Path to pedigree file')
parser.add_argument('ncausal',type=int,help='Number of causal variants')
parser.add_argument('h2quad',type=float,help='heritability explained by combined direct, sib, and parental effects')
parser.add_argument('outprefix', type=str, help='prefix of output phenotyped file (.ped)')
parser.add_argument('--no_sib',action='store_true',default=False,help='No indirect genetic effects from sibs')
parser.add_argument('--dncor',type=float,help='Correlation between direct, sib, and parental effects (default 0.8)',default=0.8)
args=parser.parse_args()

print('reading genotypes')
gts_f = Bed(args.gts)
gts_ids = gts_f.iid
id_dict = {}
for i in range(0,gts_ids.shape[0]):
    id_dict[gts_ids[i,1]] =i

# Read causal genotypes
gts = gts_f.read().val
gts = np.array(gts)

# Find parents
ped = np.loadtxt(args.ped, dtype='S20')
# ped = np.loadtxt('relatedness/families.ped', dtype='S20', skiprows=1)

fams = np.unique(ped[:,0])

G = np.zeros((fams.shape[0]*2,4,gts.shape[1]),dtype=np.float32)
ped_out = np.zeros((fams.shape[0]*2,2),dtype='S20')
for i in range(0,fams.shape[0]):
    ped_fam = ped[ped[:,0]==fams[i],:]
    sib_1_index = id_dict[ped_fam[0,1]]
    G[2*i,0,:] = gts[sib_1_index,:]
    G[2 * i+1, 1, :] = gts[sib_1_index, :]
    sib_2_index = id_dict[ped_fam[1, 1]]
    G[2 * i+1, 0, :] = gts[sib_2_index, :]
    G[2*i,1,:] = gts[sib_2_index, :]
    G[np.array([2*i,2*i+1]),2,:] = gts[id_dict[ped_fam[2, 1]],:]
    G[np.array([2 * i, 2 * i + 1]), 3, :] = gts[id_dict[ped_fam[3, 1]], :]
    ped_out[2*i,:] = ped_fam[0,0:2]
    ped_out[2*i+1,:] = ped_fam[1,0:2]

print('simulating trait')
# Simulate genetic effects
if args.no_sib:
    effect_cov = np.zeros((3,3))
else:
    effect_cov = np.zeros((4, 4))

effect_cov[:] = args.dncor
np.fill_diagonal(effect_cov,1)
if args.no_sib:
    b = np.random.multivariate_normal(np.zeros((3)),effect_cov,(gts.shape[1]))
else:
    b = np.random.multivariate_normal(np.zeros((4)),effect_cov,(gts.shape[1]))
# # additive genetic component
if args.no_sib:
    A = G[:,0,:].dot(b[:,0])+G[:,2,:].dot(b[:,2])+G[:,3,:].dot(b[:,3])
else:
    A = G[:,0,:].dot(b[:,0])+G[:,1,:].dot(b[:,1])+G[:,2,:].dot(b[:,2])+G[:,3,:].dot(b[:,3])

# Simulate residual variance
e = np.random.randn(A.shape[0])

### Form phenotype
a_factor = np.sqrt(args.h2quad)*np.power(np.std(A,axis=0),-1)

# make phenotype
y = A*a_factor+np.sqrt(1-args.h2quad)*e

# Write phenotype
yout = np.hstack((ped_out,np.array(y,dtype=str).reshape((y.shape[0],1))))
np.savetxt(args.outprefix+'.ped',yout,fmt='%s')

# Write effects
if args.no_sib:
    b_out = np.vstack((np.array(b[:,0]*a_factor,dtype=str),np.array(b[:,1]*a_factor,dtype=str),
                       np.array(b[:,2]*a_factor,dtype=str),np.array(b[:,3]*a_factor,dtype=str)))
else:
    b_out = np.vstack((np.array(b[:,0]*a_factor,dtype=str),
                       np.array(b[:,2]*a_factor,dtype=str),np.array(b[:,3]*a_factor,dtype=str)))
b_out = b_out.T
np.savetxt(args.outprefix+'.effects.txt',b_out,fmt='%s')