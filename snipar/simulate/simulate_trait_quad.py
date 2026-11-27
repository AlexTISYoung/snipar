import numpy as np
import argparse
from pysnptools.snpreader import Bed

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

print("reading genotypes")
bed = Bed(args.gts, count_A1=True)
ids = {iid: index for index, iid in enumerate(bed.iid[:, 1])}
genotypes = bed.read().val
if not np.isfinite(genotypes).all():
    raise ValueError("Trait simulation requires complete genotypes")
bim = np.loadtxt(args.gts.removesuffix(".bed") + ".bim", dtype=str, ndmin=2)

ped = np.loadtxt(args.ped, dtype=str, ndmin=2)
if tuple(ped[0, :4]) == ("FID", "IID", "FATHER_ID", "MOTHER_ID"):
    ped = ped[1:]
families = {}
for fid, iid, father, mother, *_ in ped:
    if iid in ids and father in ids and mother in ids:
        families.setdefault((fid, father, mother), []).append(iid)
families = {key: children for key, children in families.items() if len(children) == fsize}
if not families:
    raise ValueError("The pedigree and genotypes do not contain a complete family of the requested size")
if not 0 <= args.h2quad <= 1 or fsize < 1:
    raise ValueError("Require heritability between zero and one and a positive family size")
n_effects = 3 if args.no_sib else 4
if not -1 / (n_effects - 1) <= args.dncor <= 1:
    raise ValueError("Direct and indirect effect correlation is outside the valid range")

print("simulating trait")
effect_cov = np.full((n_effects, n_effects), args.dncor)
np.fill_diagonal(effect_cov, 1)
b = np.random.multivariate_normal(np.zeros(n_effects), effect_cov, genotypes.shape[1])
scores = genotypes @ b
a_components, phenotype_ids = [], []
for (fid, father, mother), children in sorted(families.items()):
    rows = np.array([ids[iid] for iid in children])
    direct = scores[rows, 0]
    inherited = scores[ids[father], -2] + scores[ids[mother], -1]
    family_effects = direct + inherited
    if not args.no_sib:
        family_effects += scores[rows, 1].sum() - scores[rows, 1]
    a_components.extend(family_effects)
    phenotype_ids.extend((fid, iid) for iid in children)
A = np.asarray(a_components)
ped_out = np.asarray(phenotype_ids)
# Simulate residual variance
e = np.random.randn(A.shape[0])

### Form phenotype
genetic_sd = np.std(A)
if genetic_sd <= 0 and args.h2quad > 0:
    raise ValueError("Simulated genetic component has zero variance")
a_factor = np.sqrt(args.h2quad) / genetic_sd if args.h2quad > 0 else 0

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
