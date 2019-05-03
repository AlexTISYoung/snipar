import h5py, sibreg
import numpy as np
from pysnptools.snpreader import Bed, Pheno

dir = '/well/donnelly/ukbiobank_project_8874/ay/sib_imputation'

test_chr = h5py.File(dir+'/imputed_parental_hdf5/chr_22.hdf5','r')


iid = np.array(np.array(test_chr['ped'],dtype=int).T,dtype='S21')

sid = np.array(test_chr['vnames'])

quads = np.loadtxt('quads.fam',dtype='S21',skiprows=1)
quads = quads[np.logical_not(quads[:,0]=='NA'),:]

quad_rows = np.zeros((quads.shape[0]),dtype=int)

for i in xrange(0,quads.shape[0]):
    quad_rows[i] = np.where(iid[:,0]==quads[i,0])[0][0]

genotypes = np.array(test_chr['gts'])
genotypes = genotypes[quad_rows,:,:]

# Get parental genotypes
par = Bed(dir+'/quad_genotypes/chr22.bed')
# select subset to test
par_iid = par.iid
par_sid = par.sid
# Find ids from same families as those with phenotype data
par_gts = par.read()
par_gts = par_gts.val

# Match rsids
sid_match = np.zeros((par_sid.shape[0]),dtype=int)
for i in xrange(0,par_sid.shape[0]):
    sid_match[i] = np.where(sid==par_sid[i])[0][0]

genotypes = genotypes[:,:,sid_match]

# Find sib indices
sib_rows = np.zeros((quads.shape[0],2),dtype=int)
for i in xrange(0,quads.shape[0]):
    sib_rows[i,0] = np.where(par_iid[:, 0] == quads[i, 1])[0][0]
    sib_rows[i, 1] = np.where(par_iid[:, 0] == quads[i, 2])[0][0]

# Get sib genotypes
sib_par_gts = np.zeros((quads.shape[0],par_sid.shape[0],2),dtype=int)

for i in xrange(0,quads.shape[0]):
    sib_par_gts[i,:,0] = par_gts[sib_rows[i,0],:]
    sib_par_gts[i, :, 1] = par_gts[sib_rows[i, 1], :]

sib_par_gts[sib_par_gts<0] = -1

np.savetxt(dir+'/imputed_parental_hdf5/test/par_sib_gts.txt',
           np.hstack((sib_par_gts[:,:,0],sib_par_gts[:,:,1])),fmt='%d')

np.savetxt(dir+'/imputed_parental_hdf5/test/sib_gts.txt',
           np.hstack((genotypes[:,0,:],genotypes[:,1,:])),fmt='%d')

# Find par indices
par_rows = np.zeros((quads.shape[0],2),dtype=int)
for i in xrange(0,quads.shape[0]):
    par_rows[i,0] = np.where(par_iid[:, 0] == quads[i, 3])[0][0]
    par_rows[i, 1] = np.where(par_iid[:, 0] == quads[i, 4])[0][0]

quad_par_gts = np.zeros((quads.shape[0],par_sid.shape[0],2),dtype=int)

for i in xrange(0,quads.shape[0]):
    quad_par_gts[i,:,0] = par_gts[par_rows[i,0],:]
    quad_par_gts[i, :, 1] = par_gts[par_rows[i, 1], :]

par_sum = np.sum(quad_par_gts,axis=2)

np.savetxt(dir+'/imputed_parental_hdf5/test/par_gts.txt',par_sum,fmt='%d')

# Imputed par gts
np.savetxt(dir+'/imputed_parental_hdf5/test/imputed_par_gts.txt',genotypes[:,2,:])
