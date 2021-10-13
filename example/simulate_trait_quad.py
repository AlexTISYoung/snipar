import numpy as np
import argparse
from pysnptools.snpreader import Bed

def convert_str_array(x):
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.decode('UTF-8') for y in x])
    return x_out.reshape(x_shape)

parser = argparse.ArgumentParser()
parser.add_argument('gts', type=str, help='Path to bed file with sibling genotypes')
parser.add_argument('ped', type=str, help='Path to pedigree file')
parser.add_argument('h2quad',type=float,help='heritability explained by combined direct, sib, and parental effects')
parser.add_argument('outprefix', type=str, help='prefix of output phenotyped file (.ped)')
parser.add_argument('--no_sib',action='store_true',default=False,help='No indirect genetic effects from sibs')
parser.add_argument('--dncor',type=float,help='Correlation between direct, sib, and parental effects (default 0.8)',default=0.8)
parser.add_argument('--fsize',type=int,help='family size (number of children). Default is 2.',default=2)
args=parser.parse_args()

fsize = args.fsize

print('reading genotypes')
gts_f = Bed(args.gts,count_A1 = True)
gts_ids = gts_f.iid
#gts_ids = np.loadtxt('fsize4.fam', dtype='S20')
gts_ids = gts_ids[:,0:2]
id_dict = {}
for i in range(0,gts_ids.shape[0]):
    id_dict[gts_ids[i,1]] =i

# Read causal genotypes
gts = gts_f.read().val
gts = np.array(gts)

# Read bim
bim = np.loadtxt(args.gts.split('.bed')[0]+'.bim',dtype=str)

# Find parents
ped = convert_str_array(np.loadtxt(args.ped,dtype='S'))
ped = ped[1:ped.shape[0],:]

fams = np.unique(ped[:,0])

G = np.zeros((fams.shape[0]*fsize,4,gts.shape[1]),dtype=np.int8)
ped_out = np.zeros((fams.shape[0]*fsize,2),dtype='S20')
for i in range(0,fams.shape[0]):
    ped_fam = ped[ped[:,0]==fams[i],:]
    sib_indices = np.array([id_dict[ped_fam[j,1]] for j in range(0,fsize)])
    # Genotypes
    for j in range(0,fsize):
        G[fsize*i+j,0,:] = gts[sib_indices[j],:]
    # Average of other sib genotypes
    for j in range(0,fsize):
        sib_indices_j = np.delete(sib_indices,j)
        G[fsize*i+j,1,:] = np.sum(gts[sib_indices_j,:],axis=0)
    G[(fsize*i):(fsize*(i+1)),2,:] = gts[id_dict[ped_fam[fsize, 1]],:]
    G[(fsize*i):(fsize*(i+1)),3,:] = gts[id_dict[ped_fam[fsize+1, 1]],:]
    ped_out[(fsize*i):(fsize*(i+1)),:] = ped_fam[0:fsize,0:2]

print('simulating trait')
# Simulate genetic effects
if args.no_sib:
    effect_cov = np.zeros((3,3))
    G = G[:,np.array([0,2,3]),:]
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
    A = G[:,0,:].dot(b[:,0])+G[:,1,:].dot(b[:,1])+G[:,2,:].dot(b[:,2])
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
np.savetxt(args.outprefix+'.txt',yout,fmt='%s')

b = b*a_factor

# Write effects
if args.no_sib:
    b_out = np.vstack((np.array(b[:,0],dtype=str),
                       np.array(b[:,1],dtype=str),np.array(b[:,2],dtype=str)))
else:
    b_out = np.vstack((np.array(b[:,0],dtype=str),np.array(b[:,1],dtype=str),
                       np.array(b[:,2],dtype=str),np.array(b[:,3],dtype=str)))
b_out = b_out.T
np.savetxt(args.outprefix+'.effects.txt',b_out,fmt='%s')

# Write direct effect weights
bim_out = bim[:,[0,3,1,4,5]]
dout = open(args.outprefix+'.direct_weights.txt','w')
dout.write('chrom\tpos\tsid\tnt1\tnt2\traw_beta\tldpred_beta\n')
for i in range(b.shape[0]):
    doutline = 'chrom'+bim_out[i,0]+'\t'+bim_out[i,1]+'\t'+bim_out[i,2]+'\t'+bim_out[i,3]+\
              '\t'+bim_out[i,4]+'\t'+str(b[i,0])+'\t'+str(b[i,0])
    if i<(b.shape[0]-1):
        doutline += '\n'
    dout.write(doutline)
dout.close()