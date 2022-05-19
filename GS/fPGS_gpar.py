import numpy as np
from snipar.utilities import * 
from snipar.gtarray import gtarray
import snipar.read as read

dir = '/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/'
pgs_file = dir+'pgs/pgs_gpar_full.txt'
phenofile = dir+'processed_traits_noadj.txt'
phen_index = 16
hsq_file = dir+'grms/varcomps/16.hsq'
mgrm = dir+'grms/mgrm.txt'

### Read PGS ###
pgs_f = open(pgs_file, 'r')
pgs_header = pgs_f.readline().split(' ')
pgs_header[len(pgs_header) - 1] = pgs_header[len(pgs_header) - 1].split('\n')[0]

ncols = len(pgs_header)
pgs_cols = tuple([x for x in range(2, ncols)])
fams = np.loadtxt(pgs_file, usecols=0, dtype=str, skiprows=1)
X = np.ones((fams.shape[0],len(pgs_cols)+1))
X[:,1:X.shape[1]] = np.loadtxt(pgs_file, usecols=pgs_cols, skiprows=1)
pg = gtarray(X,
             np.loadtxt(pgs_file, usecols=1, dtype=str, skiprows=1),
             sid=np.array(['intercept']+pgs_header[2:ncols]),
             fams=fams)

## Read phenotype
y = read.phenotype.read_phenotype(phenofile, phen_index=phen_index)
print('Number of non-missing phenotype observations: ' + str(y.shape[0]))
# Remove individuals without phenotype observations from PGS
pg.filter_ids(y.ids)
y.filter_ids(pg.ids)
print('Final sample size of individuals with complete phenotype and PGS observations: ' + str(y.shape[0]))
# Scale phenotype
y.scale()

####### Estimate effects ###########
### Load variance components
hsq = open(hsq_file,'r')
varcomps = np.zeros((3))
hsq_line = hsq.readline()
while len(hsq_line)>0:
    hsq_line = hsq_line.split('\t')
    if hsq_line[0]=='V(G1)/Vp':
        varcomps[0] = float(hsq_line[1])
    if hsq_line[0]=='V(G2)/Vp':
        varcomps[1] = float(hsq_line[1])
    if hsq_line[0]=='V(G3)/Vp':
        varcomps[2] = float(hsq_line[1])
    hsq_line = hsq.readline()
hsq.close()

## Load GRMs ##
pheno_grm = np.zeros((y.shape[0],y.shape[0]),dtype=np.float32)
grms = np.loadtxt(mgrm,dtype=str)
for i in range(grms.shape[0]):
    print('Reading GRM: '+grms[i])
    grm_ids = np.loadtxt(grms[i]+'.grm.id',dtype=str,usecols=1)
    grm = np.zeros((grm_ids.shape[0],grm_ids.shape[0]),dtype=np.float32)
    grm[np.tril_indices(grm.shape[0])] = np.fromfile(grms[i]+'.grm.bin',dtype=np.float32)
    grm += grm.T
    np.fill_diagonal(grm,np.diag(grm)/2.0)
    grm_dict = make_id_dict(grm_ids)
    y_indices = np.array([grm_dict[x] for x in y.ids])
    pheno_grm += varcomps[i]*grm[np.ix_(y_indices,y_indices)]

# add residual
np.fill_diagonal(pheno_grm,np.diag(pheno_grm)+1-np.sum(varcomps))

# Estimate effects
pheno_grm_inv = np.linalg.inv(pheno_grm)
xtx = pg.gts.T.dot(pheno_grm_inv.dot(pg.gts))
xty = pg.gts.T.dot(pheno_grm_inv.dot(y.gts))
alpha = np.linalg.solve(xtx,xty)
alpha_cov = np.linalg.inv(xtx)
alpha_ses = np.sqrt(np.diag(alpha_cov))
# Save output
alpha_out = np.vstack((pg.sid,alpha.reshape((alpha.shape[0])),alpha_ses)).T
np.savetxt(dir+'pgs/fpgs_gpar.effects.txt',alpha_out,fmt='%s')
np.savetxt(dir+'pgs/fpgs_gpar.vcov.txt',alpha_cov,fmt='%s')

## Drop grandparents and estimate effects
pg_dim = pg.gts.shape[1]
pg.gts = pg.gts[:,0:(pg.gts.shape[1]-4)]
xtx = pg.gts.T.dot(pheno_grm_inv.dot(pg.gts))
xty = pg.gts.T.dot(pheno_grm_inv.dot(y.gts))
alpha = np.linalg.solve(xtx,xty)
alpha_cov = np.linalg.inv(xtx)
alpha_ses = np.sqrt(np.diag(alpha_cov))
# Save output
alpha_out = np.vstack((pg.sid[0:pg.gts.shape[1]],alpha.reshape((alpha.shape[0])),alpha_ses)).T
np.savetxt(dir+'pgs/fpgs_par.effects.txt',alpha_out,fmt='%s')
np.savetxt(dir+'pgs/fpgs_par.vcov.txt',alpha_cov,fmt='%s')