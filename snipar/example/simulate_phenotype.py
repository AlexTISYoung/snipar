from snipar.lmm import *
from pysnptools.snpreader import Bed

bedfile = 'whole_genome'
outdir = '/disk/genetics/ukb/alextisyoung/haplotypes/relatives/ibd/'

bed = Bed(bedfile+'.bed',count_A1=True)
bim = np.loadtxt(bedfile+'.bim',dtype=str)
ids = bed.iid

ncausal = 1500
h2 = 0.75

causal = np.sort(np.random.randint(0,bim.shape[0],ncausal))
bim = bim[causal,:]

gts = np.array(bed[:,causal].read().val)

a = np.random.randn(ncausal)

A = gts.dot(a)

c = np.sqrt(h2/np.var(A))

A = c*A

a = c*a

bim = np.hstack((bim,a.reshape(a.shape[0],1),np.zeros((a.shape[0],1))))

Y = A+np.sqrt(1-h2)*np.random.randn(A.shape[0])

Y_out = np.hstack((ids,Y.reshape(Y.shape[0],1)))

np.savetxt(outdir+'h2_0.75.fam',Y_out,fmt='%s')
np.savetxt(outdir+'h2_0.75_causal.txt',bim,fmt='%s')