import numpy as np
from pysnptools.snpreader import Pheno
from sibreg import sibreg

def id_dict_make(ids):
## Make a dictionary mapping from IDs to indices ##
    if not type(ids)==np.ndarray:
        raise(ValueError('Unsupported ID type: should be numpy nd.array'))
    id_dict={}
    for id_index in xrange(0,len(ids)):
        id_dict[tuple(ids[id_index,:])]=id_index
    return id_dict


pheno = Pheno('phenotypes/height_resid.ped', missing='NA').read()

y = np.array(pheno.val)
pheno_ids = np.array(pheno.iid)
y=y[:,0]

y_not_nan = np.logical_not(np.isnan(y))
if np.sum(y_not_nan) < y.shape[0]:
    y = y[y_not_nan]
    pheno_ids = pheno_ids[y_not_nan, :]

# Read PCs
PCs = Pheno('pcs/UKB_PCs.ped', missing='NA').read()
pcs = np.array(PCs.val)
pc_ids = np.array(PCs.iid)
# remove nan
notnan = np.logical_not(np.isnan(pcs[:,0]))
pcs = pcs[notnan,:]
pc_ids = pc_ids[notnan,:]
# Make dict
pc_id_dict = id_dict_make(np.array(pc_ids))
# normalise
for j in xrange(0,pcs.shape[1]):
    pcs[:,j] = (pcs[:,j]-np.mean(pcs[:,j]))/np.std(pcs[:,j])

# Compute within family mean genotypes
pc_means = np.zeros(pcs.shape)
families = np.unique(pc_ids[:,0])
for family in families:
    family_rows = pc_ids[:,0]==family
    pc_mean = np.mean(pcs[family_rows,:],axis=0)
    pc_means[family_rows,:] = pc_mean

pc_diff = pcs - pc_means

# Write in GCTA format
np.savetxt('pcs/UKB_PC_BF_WF.gcta',np.hstack((pc_means,pc_diff)))
np.savetxt('pcs/UKB_PC_ids.txt',np.array(pc_ids,dtype=int),fmt='%d')

# match with pheno ids
id_match = np.zeros((y.shape[0]),dtype=int)
in_pc_ids = np.zeros((y.shape[0]),dtype=bool)
for i in xrange(0,y.shape[0]):
    pid = tuple(pheno_ids[i,:])
    if pid in pc_id_dict:
        id_match[i] = pc_id_dict[pid]
        in_pc_ids[i] = True

y = y[in_pc_ids]
pheno_ids = pheno_ids[in_pc_ids,:]

pc_means = pc_means[id_match,:]
pc_means = pc_means[in_pc_ids]
pc_diff = pc_diff[id_match,:]
pc_diff = pc_diff[in_pc_ids,:]

pc_all = np.hstack((np.ones((pc_means.shape[0],1)),pc_means,pc_diff))

# Write in GCTA format


model_l = sibreg.model(y,pc_all,pheno_ids[:,0])

sigma_2_init = np.var(y)*0.5

model_optim = model_l.optimize_model(np.array([sigma_2_init,1]))

print('Within family variance estimate: ' + str(round(model_optim['sigma2'] / model_optim['tau'], 4)))
print('Residual variance estimate: ' + str(round(model_optim['sigma2'], 4)))

alpha = model_l.alpha_mle(model_optim['tau'],model_optim['sigma2'],compute_cov = True)

alpha_est = alpha[0][1:alpha[0].shape[0]]
alpha_cov = alpha[1][1:alpha[0].shape[0], 1:alpha[0].shape[0]]
alpha_ses = np.sqrt(np.diag(alpha_cov))
alpha_out = np.zeros((alpha[0].shape[0]-1,3))
alpha_out[:,0] = alpha_est
alpha_out[:,1] = alpha_ses
alpha_out[:,2] = alpha_est/alpha_ses

np.savetxt('pcs/eduyears_resid_white_british_PC_ests_sibreg.txt',alpha_out)