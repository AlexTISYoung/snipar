#!/well/kong/users/wiw765/anaconda2/bin/python
import numpy as np
import numpy.ma as ma
from pysnptools.snpreader import Bed, Pheno
from sibreg import sibreg
import h5py, argparse, code

def make_id_dict(x,col=0):
    if len(x.shape)>1:
        x = x[:,col]
    id_dict = {}
    for i in range(0,x.shape[0]):
        id_dict[x[i]] = i
    return id_dict

def convert_str_array(x):
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.decode('UTF-8') for y in x])
    return x_out.reshape(x_shape)

def encode_str_array(x):
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.encode('ascii') for y in x])
    return x_out.reshape(x_shape)

def find_par_gts(pheno_ids,ped,fams,gts_id_dict):
    # Whether mother and father have observed/imputed genotypes
    par_status = np.zeros((pheno_ids.shape[0],2),dtype=int)
    par_status[:] = -1
    # Indices of obsered/imputed genotypes in relevant arrays
    gt_indices = np.zeros((pheno_ids.shape[0],3),dtype=int)
    gt_indices[:] = -1
    ## Build dictionaries
    # Where each individual is in the pedigree
    ped_dict = make_id_dict(ped,1)
    # Where the imputed data is for each family
    fam_dict = make_id_dict(fams)
    # Store family ID of each individual
    fam_labels = np.zeros((pheno_ids.shape[0]),dtype=fams.dtype)
    # Find status and find indices
    for i in range(0,pheno_ids.shape[0]):
        # Find index in genotypes
        if pheno_ids[i] in gts_id_dict:
            gt_indices[i,0] = gts_id_dict[pheno_ids[i]]
        # Find index in pedigree
        if pheno_ids[i] in ped_dict:
            ped_i = ped[ped_dict[pheno_ids[i]],:]
            fam_labels[i] = ped_i[0]
            # Check for observed father
            if ped_i[2] in gts_id_dict:
                gt_indices[i,1] = gts_id_dict[ped_i[2]]
                par_status[i,0] = 0
            # Check for observed mother
            if ped_i[3] in gts_id_dict:
                gt_indices[i, 2] = gts_id_dict[ped_i[3]]
                par_status[i,1] = 0
            # If parent not observed, look for imputation
            if ped_i[0] in fam_dict:
                imp_index = fam_dict[ped_i[0]]
                # Check if this is imputation of father, or mother, or both
                if ped_i[4] == 'False':
                    gt_indices[i,1] = imp_index
                    par_status[i,0] = 1
                if ped_i[5] == 'False':
                    gt_indices[i, 2] = imp_index
                    par_status[i, 1] = 1

    return par_status, gt_indices, fam_labels

def make_gts_matrix(gts,imp_gts,par_status,gt_indices,mean_normalise = True):
    if np.min(gt_indices)<0:
        raise(ValueError('Missing genotype index'))
    N = gt_indices.shape[0]
    G = np.zeros((N,3,gts.shape[1]),np.float32)
    for i in range(0,N):
        # Observed proband genotype
        G[i,0,:] = gts[gt_indices[i,0],:]
        # Paternal genotype
        if par_status[i,0] == 0:
            G[i,1,:] = gts[gt_indices[i,1],:]
        elif par_status[i,0] == 1:
            G[i,1,:] = imp_gts[gt_indices[i,1],:]
        else:
            ValueError('Paternal genotype neither imputed nor observed')
        # Maternal genotype
        if par_status[i, 1] == 0:
            G[i, 2, :] = gts[gt_indices[i, 2], :]
        elif par_status[i, 1] == 1:
            G[i, 2, :] = imp_gts[gt_indices[i, 2], :]
        else:
            ValueError('Maternal genotype neither imputed nor observed')
    if mean_normalise:
        for i in range(3):
            G[:,i,:] = G[:,i,:] - ma.mean(G[:,i,:],axis=0)
    G = G.transpose(2,0,1)
    return G

def fit_null_model(y,X,fam_labels,tau_init):
    # Optimize null model
    sigma_2_init = np.var(y)*tau_init/(1+tau_init)
    null_model = sibreg.model(y, X, fam_labels)
    null_optim = null_model.optimize_model(np.array([sigma_2_init,args.tau_init]))
    null_alpha = null_model.alpha_mle(null_optim['tau'],null_optim['sigma2'],compute_cov = True)
    # Return values
    return null_optim['sigma2'], null_optim['tau'], null_alpha


######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('gts',type=str,help='Path to bed file with sibling genotypes')
    parser.add_argument('pargts', type=str, help='Path to HDF5 file with imputed parental genotypes')
    parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
    parser.add_argument('outprefix',type=str,help='Location to output association statistic hdf5 file')
    parser.add_argument('--sibped',type=str,help='Path to pedigree file. By default uses pedigree in imputed parental genotype HDF5 file',default=None)
    parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                        default=1)
    parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                        default=1)
    parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)',default=0.01)
    parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)',default='NA')
    parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)',default=5)
    parser.add_argument('--append',action='store_true',default=False,help='Append results to existing output file with given outprefix (default overwrites existing')
    parser.add_argument('--no_covariate_estimates',action='store_true',default=False,help='Suppress output of covariate effect estimates')
    parser.add_argument('--fit_sib',action='store_true',default=False,help='Fit indirect effects from siblings')
    args=parser.parse_args()

    ######### Read Phenotype ########
    pheno = Pheno(args.phenofile, missing=args.missing_char).read()
    # pheno = Pheno('phenotypes/eduyears_resid.ped', missing='NA').read()
    y = np.array(pheno.val)
    pheno_ids = np.array(pheno.iid)[:,1]
    if y.ndim == 1:
        pass
    elif y.ndim == 2:
        y = y[:, args.phen_index - 1]
    else:
        raise (ValueError('Incorrect dimensions of phenotype array'))
    # Remove y NAs
    y_not_nan = np.logical_not(np.isnan(y))
    if np.sum(y_not_nan) < y.shape[0]:
        y = y[y_not_nan]
        pheno_ids = pheno_ids[y_not_nan]
    print('Number of non-missing phenotype observations: ' + str(y.shape[0]))

    ####### Find parental status #######
    ### Imputed parental file ###
    par_gts_f = h5py.File(args.pargts,'r')
    # Get families
    fams = convert_str_array(np.array(par_gts_f['families']))
    # Get pedigree
    ped = convert_str_array(np.array(par_gts_f['pedigree']))
    ped = ped[1:ped.shape[0],:]

    ### Genotype file ###
    gts_f = Bed(args.gts)
    # get ids of genotypes and make dict
    gts_ids = gts_f.iid[:,1]
    gts_id_dict = make_id_dict(gts_ids)

    ### Find parental status
    print('Checking for observed/imputed parental genotypes')
    par_status, gt_indices, fam_labels = find_par_gts(pheno_ids,ped,fams,gts_id_dict)
    # Find which individuals can be used
    none_missing = np.min(gt_indices, axis=1)
    none_missing = none_missing >= 0
    N = np.sum(none_missing)
    if N==0:
        raise (ValueError(
            'No individuals with phenotype observations and complete observed/imputed genotype observations'))
    print(str(N)+' individuals with phenotype observations and complete observed/imputed genotypes observations')
    # Take those that can be used
    pheno_ids = pheno_ids[none_missing]
    y = y[none_missing]
    gt_indices = gt_indices[none_missing,:]
    par_status = par_status[none_missing,:]
    ## Read genotypes
    observed_indices = np.sort(np.unique(np.hstack((gt_indices[:,0],
                                  gt_indices[par_status[:,0]==0,1],
                                  gt_indices[par_status[:,1]==0,2]))))
    # Get imputed indices
    imp_indices = np.sort(np.unique(np.hstack((gt_indices[par_status[:,0]==1,1],
                                  gt_indices[par_status[:,1]==1,2]))))
    print('Matching observed and imputed SNPs')
    # Match SNPs from imputed and observed
    imp_sid = convert_str_array(np.array(par_gts_f['sid']))
    obs_sid = gts_f.sid
    obs_sid_dict = make_id_dict(obs_sid)
    in_obs_sid = np.zeros((imp_sid.shape[0]),dtype=bool)
    obs_sid_index = np.zeros((imp_sid.shape[0]),dtype=int)
    for i in range(0,imp_sid.shape[0]):
        if imp_sid[i] in obs_sid_dict:
            in_obs_sid[i] = True
            obs_sid_index[i] = obs_sid_dict[imp_sid[i]]
    sid = imp_sid[in_obs_sid]
    if np.sum(in_obs_sid) == 0:
        ValueError('No SNPs in common between imputed and observed genotypes')
    # Read imputed parental genotypes
    print('Reading imputed parental genotypes')
    imp_gts = np.array(par_gts_f['imputed_par_gts'][imp_indices,:])
    imp_gts = imp_gts[:,in_obs_sid]
    fams = fams[imp_indices]
    # Read observed genotypes
    print('Reading observed genotypes')
    gts = gts_f[observed_indices, obs_sid_index].read().val
    gts = ma.array(gts,mask=np.isnan(gts))
    gts_id_dict = make_id_dict(gts_ids[observed_indices])
    # Filter genotypes based on frequency and missingness
    print('Filtering on frequency and MAF')
    freqs = ma.mean(gts,axis=0)/2.0
    missingness = ma.mean(gts.mask,axis=0)
    freqs_pass = np.logical_and(freqs > args.min_maf,freqs < (1-args.min_maf))
    print(str(freqs.shape[0]-np.sum(freqs_pass))+' SNPs with MAF<'+str(args.min_maf))
    missingness_pass = 100*missingness < args.max_missing
    print(str(freqs.shape[0] - np.sum(freqs_pass)) + ' SNPs with missingness >' + str(args.max_missing)+'%')
    filter_pass = np.logical_and(freqs_pass,missingness_pass)
    gts = gts[:,filter_pass]
    imp_gts = imp_gts[:,filter_pass]
    sid = sid[filter_pass]
    freqs = freqs[filter_pass]
    N_L = np.sum(np.logical_not(gts.mask),axis=0)
    print('After filtering, '+str(np.sum(filter_pass))+' SNPs remain')
    # Find indices in reduced data
    par_status, gt_indices, fam_labels = find_par_gts(pheno_ids, ped, fams, gts_id_dict)
    # Check for empty fam labels
    no_fam = np.array([len(x)==0 for x in fam_labels])
    if np.sum(no_fam)>0:
        ValueError('No family label from pedigree for some individuals')
    print('Constructing family based genotype matrix')
    ### Make genotype design matrix
    G = make_gts_matrix(gts,imp_gts,par_status,gt_indices)
    del gts
    del imp_gts
    # Fill NAs
    print('Imputing missing genotypes with population frequency')
    G[np.isnan(G)] = 0
    #### Fit null model ####
    print('Estimating variance components')
    null_model = sibreg.model(y,np.ones((y.shape[0],1)),fam_labels)
    sigma_2_init = np.var(y) * args.tau_init / (1 + args.tau_init)
    null_optim = null_model.optimize_model(np.array([sigma_2_init,args.tau_init]))
    sigma2 = null_optim['sigma2']
    tau = null_optim['tau']
    print('Family variance estimate: '+str(round(sigma2/tau,4)))
    print('Residual variance estimate: ' + str(round(sigma2,4)))
    ##### Transform ######
    #### Get inverse square root of Sigma ###
    print('Transforming genotypes and phenotypes')
    L = null_model.sigma_inv_root(tau,sigma2)
    #### Transform genotype and phenotype ###
    # Mean normalise
    y = y - np.mean(y)
    for label in L.keys():
        label_indices = null_model.label_indices[label]
        y[label_indices] = np.dot(L[label],y[label_indices])
        for i in range(3):
            G[:,label_indices,i] = np.dot(G[:,label_indices,i],L[label].T)
    ### Fit models for SNPs ###
    print('Estimating SNP effects')
    XTX = np.einsum('...ij,...ik', G, G)
    XTY = np.einsum('...ij,i',G,y)
    alpha = np.linalg.solve(XTX,XTY)
    alpha_cov = np.linalg.inv(XTX)
    alpha_ses = np.sqrt(np.diagonal(alpha_cov,axis1=1,axis2=2))
    ### Output file ###
    print('Writing output to '+args.outprefix+'.hdf5')
    outfile = h5py.File(args.outprefix+'.hdf5','w')
    outfile['sid'] = encode_str_array(sid)
    X_length = 3
    outfile.create_dataset('estimate_covariance',(sid.shape[0],X_length,X_length),dtype = 'f',chunks = True, compression = 'gzip', compression_opts=9)
    outfile.create_dataset('estimate', (sid.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile.create_dataset('estimate_ses', (sid.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile['estimate'][:] = alpha
    outfile['estimate_ses'][:] = alpha_ses
    outfile['estimate_covariance'][:] = alpha_cov
    outfile['sigma2'] = sigma2
    outfile['tau'] = tau
    outfile['N_L'] = N_L
    outfile['freqs'] = freqs
    outfile.close()