import numpy as np
import h5py, argparse
from pysnptools.snpreader import Bed

parser = argparse.ArgumentParser()
parser.add_argument('gts', type=str, help='Path to bed file with sibling genotypes')
parser.add_argument('ped', type=str, help='Path to pedigree file')
parser.add_argument('ncausal',type=int,help='Number of causal variants')
parser.add_argument('h2sib',type=float,help='heritability of trait')
parser.add_argument('nrep',type=int,help='Number of phenotypes to simulate')
parser.add_argument('outprefix', type=str, help='Location to output csv file with association statistics')
parser.add_argument('--dncor',type=float,help='Correlation between direct and parental effects',default=0.8)
args=parser.parse_args()

print('reading genotypes')
gts_f = Bed(args.gts)
#gts_f = Bed('genotypes/causal_sim/causal_sim.bed')
gts_ids = gts_f.iid
id_dict = {}
for i in range(0,gts_ids.shape[0]):
    id_dict[gts_ids[i,1]] =i

# Sample causal variants
sid = gts_f.sid
causal_indices = np.sort(np.random.choice(np.array([x for x in range(0,sid.shape[0])]),size=args.ncausal,replace=False))
causal_sid = sid[causal_indices]

# Read causal genotypes
gts = gts_f[:,causal_indices].read().val
gts = np.array(gts)
gts_nans = np.isnan(gts)

# Mean impute NAs
for j in range(0,gts.shape[1]):
    gts[gts_nans[:,j],j] = np.mean(gts[np.logical_not(gts_nans[:,j]),j])

# Find parents
ped = np.loadtxt(args.ped, dtype='S20', skiprows=1)
# ped = np.loadtxt('relatedness/families.ped', dtype='S20', skiprows=1)

genotyped = np.array([x in id_dict for x in ped[:,1]])
ped = ped[genotyped,:]

sibcount = np.zeros((ped.shape[0]),dtype=int)
sibs_list = []
for i in xrange(0,ped.shape[0]):
    sibs_i = np.logical_and(ped[:,0]==ped[i,0],np.logical_and(ped[:,2]==ped[i,2],ped[:,3]==ped[i,3]))
    sibs_i[i] = False
    sibcount[i] = np.sum(sibs_i)
    if sibcount[i]>0:
        sibs_list.append(ped[sibs_i,1])

has_sibs = sibcount>0
print('Removing '+str(np.sum(sibcount==0))+' individuals without genotyped siblings')
ped = ped[has_sibs,:]

G = np.zeros((ped.shape[0],2,gts.shape[1]),dtype=np.float32)
G[:] = np.nan
G[:,0,:] = gts[np.array([id_dict[x] for x in ped[:,1]]),:]
gtcount = 0
for i in range(0,G.shape[0]):
    G[i, 1, :] = np.mean(gts[np.array([id_dict[x] for x in sibs_list[i]]), :], axis=0)

print('simulating trait')
# Simulate genetic effects
b = np.random.multivariate_normal(np.zeros((2)),np.array([[1,args.dncor],[args.dncor,1]]),(gts.shape[1],args.nrep))
# additive genetic component
a = G[:,0,:].dot(b[:,:,0])
a_par = G[:,1,:].dot(b[:,:,1])

# Simulate residual variance
e = np.random.randn(ped.shape[0],args.nrep)

### Form phenotype
## standardise
# genetic effects
A = a+a_par
a_factor = np.sqrt(args.h2sib)*np.power(np.std(A,axis=0),-1)

# make phenotype
y = A*a_factor+np.sqrt(1-args.h2sib)*e

# Write phenotype
yout = np.hstack((ped[:,0:2],np.array(y,dtype=str)))
np.savetxt(args.outprefix+'.ped',yout,fmt='%s')

# Write effects
b_out = np.hstack((causal_sid.reshape((causal_sid.shape[0],1)),np.array(b[:,:,0]*a_factor,dtype=str),np.array(b[:,:,1]*a_factor,dtype=str)))
np.savetxt(args.outprefix+'.effects.txt',b_out,fmt='%s')