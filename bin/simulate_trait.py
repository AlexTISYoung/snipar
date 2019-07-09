import numpy as np
import h5py, argparse
from pysnptools.snpreader import Bed

parser = argparse.ArgumentParser()
parser.add_argument('gts', type=str, help='Path to bed file with sibling genotypes')
parser.add_argument('ncausal',type=int,help='Number of causal variants')
parser.add_argument('h2',type=float,help='heritability of trait')
parser.add_argument('c2',type=float,help='proportion of variance explained by family environment')
parser.add_argument('nrep',type=int,help='Number of phenotypes to simulate')
parser.add_argument('outprefix', type=str, help='Location to output csv file with association statistics')
args=parser.parse_args()

print('reading genotypes')
gts_f = Bed(args.gts)
#gts_f = Bed('genotypes/causal_sim/causal_sim.bed')
gts_ids = gts_f.iid

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

print('simulating trait')
# Simulate genetic effects
b = np.random.randn(gts.shape[1],args.nrep)
# additive genetic component
a = gts.dot(b)

# Simulate family effects
fams = np.unique(gts_ids[:,0])
fam_dict = {}
for i in range(0,fams.shape[0]):
    fam_dict[fams[i]] = i

fam_means = np.random.randn(fams.shape[0],args.nrep)
fam_effects = np.zeros((gts.shape[0],args.nrep))
for i in range(0,gts.shape[0]):
    fam_effects[i] = fam_means[fam_dict[gts_ids[i,0]],:]

# Simulate residual variance
e = np.random.randn(gts.shape[0],args.nrep)

### Form phenotype
## standardise
# genetic effects
a_factor = np.sqrt(args.h2)*np.power(np.std(a,axis=0),-1)
b = b*a_factor
# make phenotype
y = a*a_factor+np.sqrt(args.c2)*(fam_effects/np.std(fam_effects,axis=0))+np.sqrt(1-args.h2-args.c2)*e

# Write phenotype
yout = np.hstack((gts_ids,np.array(y,dtype=str)))
np.savetxt(args.outprefix+'.ped',yout,fmt='%s')

# Write effects
b_out = np.hstack((causal_sid,np.array(b,dtype=str)))
np.savetxt(args.outprefix+'.effects.txt',yout,fmt='%s')