#!/apps/well/python/2.7.8/bin/python
import argparse
import numpy as np
from pysnptools.snpreader import Bed, Pheno
from scipy.stats import chi2, zscore
from sibreg import sibreg
import h5py
import code

####### Output functions ##########
def neglog10pval(x,df):
    return -np.log10(np.e)*chi2.logsf(x,df)

def vector_out(alpha_mle, no_sib, digits=6):
##Output parameter estimates along with standard errors ##
    ## Calculate test statistics
    if args.no_sib:
        X_length = 3
    else:
        X_length = 4
    alpha_est = alpha_l[0][1:X_length]
    alpha_cov = alpha_l[1][1:X_length,1:X_length]
    alpha_ses = np.sqrt(np.diag(alpha_cov))
    alpha_corr = np.dot(np.diag(1/alpha_ses).dot(alpha_cov),np.diag(1/alpha_ses))
    alpha_out = str(round(alpha_est[0],digits))+'\t'+str(round(alpha_ses[0],digits))+'\t'
    alpha_out += str(round(alpha_est[1],digits))+'\t'+str(round(alpha_ses[1],digits))+'\t'
    if args.no_sib:
        alpha_out += str(round(alpha_est[2],digits))+'\t'+str(round(alpha_ses[2],digits))+'\t'
        alpha_out += str(round(alpha_corr[0,1],digits))+'\t'+str(round(alpha_corr[0,2],digits))+'\t'+str(round(alpha_corr[2,1],digits))+'\n'
    else:
        alpha_out += str(round(alpha_corr[0, 1], digits)) + '\n'
    return alpha_out

def id_dict_make(ids):
## Make a dictionary mapping from IDs to indices ##
    if not type(ids)==np.ndarray:
        raise(ValueError('Unsupported ID type: should be numpy nd.array'))
    id_dict={}
    for id_index in xrange(0,len(ids)):
        id_dict[tuple(ids[id_index,:])]=id_index
    return id_dict

def read_covariates(covar_file,ids_to_match,missing):
## Read a covariate file and reorder to match ids_to_match ##
    # Read covariate file
    covar_f = Pheno(covar_file, missing=missing).read()
    ids = covar_f.iid
    # Get covariate values
    n_X=covar_f._col.shape[0]+1
    X=np.ones((covar_f.val.shape[0],n_X))
    X[:, 1:n_X] = covar_f.val
    # Get covariate names
    X_names = np.zeros((n_X), dtype='S10')
    X_names[0] = 'Intercept'
    X_names[1:n_X] = np.array(covar_f._col, dtype='S20')
    # Remove NAs
    NA_rows = np.isnan(X).any(axis=1)
    n_NA_row = np.sum(NA_rows)
    if n_NA_row>0:
        print('Number of rows removed from covariate file due to missing observations: '+str(np.sum(NA_rows)))
        X = X[~NA_rows]
        ids = ids[~NA_rows]
    id_dict = id_dict_make(ids)
    # Match with pheno_ids
    ids_to_match_tuples = [tuple(x) for x in ids_to_match]
    common_ids = id_dict.viewkeys() & set(ids_to_match_tuples)
    pheno_in = np.array([(tuple(x) in common_ids) for x in ids_to_match])
    match_ids = ids_to_match[pheno_in,:]
    X_id_match = np.array([id_dict[tuple(x)] for x in match_ids])
    X = X[X_id_match, :]
    return [X,X_names,pheno_in]

######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('genofile',type=str,help='Path to genotypes in HDF5')
    parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
    parser.add_argument('outprefix',type=str,help='Location to output csv file with association statistics')
    parser.add_argument('--mean_covar',type=str,help='Location of mean covariate file (default None)',
                        default=None)
    parser.add_argument('--fit_covariates',action='store_true',
                        help='Fit covariates for each locus. Default is to fit for null model and project out (mean) and rescale (variance)',
                        default=False)
    parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                        default=1)
    parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                        default=1)
    parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.05)',default=0.05)
    parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)',default='NA')
    parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)',default=5)
    parser.add_argument('--append',action='store_true',default=False,help='Append results to existing output file with given outprefix (default overwrites existing')
    parser.add_argument('--no_covariate_estimates',action='store_true',default=False,help='Suppress output of covariate effect estimates')
    parser.add_argument('--no_sib',action='store_true',default=False,help='Do not fit indirect genetic effects from sibs')
    args=parser.parse_args()

    ####################### Read in data #########################
    #### Read phenotype ###
    pheno = Pheno(args.phenofile, missing=args.missing_char).read()
    #pheno = Pheno('gwas/h2_0.4/h2_0.4.ped', missing='NA').read()
    y = np.array(pheno.val)
    pheno_ids = np.array(pheno.iid)
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
        pheno_ids = pheno_ids[y_not_nan,:]
    pheno_id_dict = id_dict_make(np.array(pheno_ids))
    pheno_fam_id_set = set(pheno_ids[:,0])
    print('Number of non-missing phenotype observations: ' + str(y.shape[0]))

    ### Get covariates
    ## Get mean covariates
    if not args.mean_covar == None:
        X, X_names, pheno_in = read_covariates(args.mean_covar,pheno_ids, args.missing_char)
        n_X = X.shape[1]
        # Remove rows with missing values
        if np.sum(pheno_in) < y.shape[0]:
            y = y[pheno_in]
            pheno_ids = pheno_ids[pheno_in,:]
        # Normalise non-constant cols
        X_stds = np.std(X[:, 1:n_X], axis=0)
        X[:, 1:n_X] = zscore(X[:, 1:n_X], axis=0)
    else:
        X = np.ones((int(y.shape[0]), 1))
        n_X = 1
        X_names = np.array(['Intercept'])

### Read genotypes ###
    print('Reading genotype file')
    test_chr = h5py.File(args.genofile,'r')
    #test_chr = h5py.File('imputed_parental_hdf5/chr_22.hdf5','r')
    # select subset to test
    iid = np.array(np.array(test_chr['ped'],dtype=int).T,dtype='S21')
    # Find ids from same families as those with phenotype data
    fid_with_phen = np.array([str(x) in pheno_fam_id_set for x in iid[:,0]])
    iid = iid[fid_with_phen,:]
    # Get frequencies
    freqs = np.array(test_chr['freqs'])
    for i in xrange(0,freqs.shape[0]):
        if freqs[i]>0.5:
            freqs[i]=1-freqs[i]
    print('Filtering on MAF')
    maf_pass = freqs>args.min_maf
    print(str(freqs.shape[0]-np.sum(maf_pass))+' SNPs below minimum MAF ('+str(args.min_maf)+') removed')
    genotypes = test_chr['gts'][:,:,maf_pass]
    genotypes = np.array(genotypes[fid_with_phen,:,:])
    sid = np.array(test_chr['vnames'])
    sid = sid[maf_pass]
    freqs = freqs[maf_pass]
    chr_length = genotypes.shape[2]
    print('Number of test loci: ' + str(genotypes.shape[2]))

### Intersect and match genotype data with covariate and phenotype data ###
    # Intersect with phenotype IDs
    geno_id_set = np.vstack((iid[:,np.array([0,1])],iid[:,np.array([0,2])]))
    geno_id_set = {tuple(geno_id_set[x,]) for x in xrange(0,geno_id_set.shape[0])}
    ids_in_common = {tuple(x) for x in pheno_ids} & geno_id_set
    pheno_ids_in_common = np.array([tuple(x) in ids_in_common for x in pheno_ids])
    y = y[pheno_ids_in_common]
    pheno_ids = pheno_ids[pheno_ids_in_common,:]
    pheno_id_dict = id_dict_make(pheno_ids)
    # Dictionary to look up rows for each family
    pheno_fam_dict = {}
    for fam in np.unique(pheno_ids[:,0]):
        pheno_fam_dict[fam] = np.where(pheno_ids[:,0]==fam)[0]
    # Get relevant covariate rows
    X = X[pheno_ids_in_common,:]

### Construct genetic covariate matrix
    print('Constructing genetic covariates')
    geno_fam_dict = {}
    for fam in np.unique(pheno_ids[:,0]):
        geno_fam_dict[fam] = np.where(iid[:,0]==fam)[0]
    if args.no_sib:
        G = np.zeros((y.shape[0], 2, genotypes.shape[2]))
        p_index = 1
    else:
        G = np.zeros((y.shape[0], 3, genotypes.shape[2]))
        p_index = 2
    for fam in np.unique(pheno_ids[:, 0]):
        gts_fam = genotypes[geno_fam_dict[fam],:,:]
        iid_fam = iid[geno_fam_dict[fam],1:3]
        pheno_id_fam = pheno_ids[pheno_fam_dict[fam],1]
        if args.no_sib:
            G_fam = np.zeros((pheno_id_fam.shape[0], 2, genotypes.shape[2]))

        else:
            G_fam = np.zeros((pheno_id_fam.shape[0], 3, genotypes.shape[2]))
        # Average imputed parental
        G_fam[:,p_index,:] = np.mean(gts_fam[:,2,:],axis=0)
        # Get proband genotypes
        for i in xrange(0,pheno_id_fam.shape[0]):
            sib = pheno_id_fam[i]
            sibmatch = np.where(sib==iid_fam[:,0])[0]
            if len(sibmatch)>0:
                sibmatch = sibmatch[0]
                G_fam[i,0,:] = gts_fam[sibmatch,0,:]
            else:
                sibmatch = np.where(sib == iid_fam[:, 1])[0][0]
                G_fam[i,0,:] = gts_fam[sibmatch,1,:]
        ## Compute sum of sib genotypes for each proband
        if not args.no_sib:
            G_fam[:,1,:] = np.sum(G_fam[:,0,:],axis=0)
            # remove proband genotype from sum
            G_fam[:,1,:] = G_fam[:,1,:] - G_fam[:,0,:]
            # add to sib column genotypes of any sibs with genotype but not phenotype data
            if iid_fam.shape[0]>pheno_id_fam.shape[0]:
                geno_sib_ids = np.unique(iid_fam[:,1:3])
                for sib in geno_sib_ids:
                    if sib not in pheno_id_fam:
                        sibmatch = np.where(sib == iid_fam[:, 0])[0]
                        if len(sibmatch) > 0:
                            sibmatch = sibmatch[0]
                            G_fam[:, 1, :] += gts_fam[sibmatch, 0, :]
                        else:
                            sibmatch = np.where(sib == iid_fam[:, 1])[0][0]
                            G_fam[:, 1, :] += gts_fam[sibmatch, 1, :]
        ## Set in full matrix
        G[pheno_fam_dict[fam],:,:] = G_fam
    del genotypes
    del G_fam

### Get sample size
    n = y.shape[0]
    if n == 0:
        raise (ValueError('No non-missing observations with both phenotype and genotype data'))
    print(str(n) + ' individuals in genotype file with no missing phenotype or covariate observations')
    n = float(n)

######### Initialise output files #######
    ## Output file
    if args.append:
        write_mode='ab'
    else:
        write_mode='wb'
    outfile=open(args.outprefix+'.models.gz',write_mode)
    if not args.append:
        if args.no_sib:
            header = 'SNP\tfrequency\tn\tdelta\tdelta_SE\tbeta\tbeta_SE\t'
            header += 'r_delta_beta\n'
        else:
            header = 'SNP\tfrequency\tn\tdelta\tdelta_SE\teta_s\teta_s_SE\tbeta\tbeta_SE\t'
            header += 'r_delta_eta_s\tr_delta_beta\tr_eta_s_beta\n'
        outfile.write(header)

######### Fit Null Model ##########
    ## Get initial guesses for null model
    print('Fitting Null Model')
    # Optimize null model
    sigma_2_init = np.var(y)*args.tau_init/(1+args.tau_init)
    #sigma_2_init = np.var(y) * 1 / (1 + 1)
    null_model = sibreg.model(y, X, pheno_ids[:,0])
    null_optim = null_model.optimize_model(np.array([sigma_2_init,args.tau_init]))
    print('Within family variance estimate: '+str(round(null_optim['sigma2']/null_optim['tau'],4)))
    print('Residual variance estimate: ' + str(round(null_optim['sigma2'],4)))
    #null_optim = null_model.optimize_model(np.array([sigma_2_init,1]))
    null_alpha = null_model.alpha_mle(null_optim['tau'],null_optim['sigma2'],compute_cov = True)
    #code.interact(local=locals())
    ## Record fitting of null model
    # Get print out for fixed mean effects
    alpha_out=np.zeros((n_X,2))
    alpha_out[:,0]=null_alpha[0]
    alpha_out[:,1]=np.sqrt(np.diag(null_alpha[1]))
    # Rescale
    if n_X>1:
        for i in xrange(0,2):
            alpha_out[1:n_X,i] = alpha_out[1:n_X,i]/X_stds
    if not args.append and not args.no_covariate_estimates and args.mean_covar is not None:
        np.savetxt(args.outprefix + '.null_covariate_effects.txt',
                   np.hstack((X_names.reshape((n_X, 1)), np.array(alpha_out, dtype='S20'))),
                   delimiter='\t', fmt='%s')

    # Fit SNP specific models
    ### Project out mean covariates
    if not args.fit_covariates:
        # Residual y
        y=y-X.dot(null_alpha[0])
        # Reformulate fixed_effects
        X=np.ones((int(n),1))
        n_X=1

    ############### Loop through loci and fit models ######################
    print('Fitting models for genome-wide SNPs')
    for loc in xrange(0,G.shape[2]):
        if args.no_sib:
            alpha_out = 'NA\tNA\tNA\tNA\tNA\n'
        else:
            alpha_out = 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n'
        # Optimize model for SNP
        if args.no_sib:
            X_length = 3
        else:
            X_length = 4
        X_l = np.ones((y.shape[0], X_length))
        X_l[:, 1:X_length] = G[:, :, loc]
        model_l = sibreg.model(y,X_l,pheno_ids[:,0])
        optim_l = model_l.optimize_model(np.array([null_optim['sigma2'],null_optim['tau']]))
        if optim_l['success']:
            alpha_l = model_l.alpha_mle(optim_l['tau'],optim_l['sigma2'],compute_cov = True)
            alpha_out = vector_out(alpha_l, args.no_sib)
        else:
            print('Maximisation of likelihood failed for for '+sid[loc])
        print(sid[loc]+' finished successfully')
        outfile.write(sid[loc] +'\t'+str(freqs[loc])+'\t'+str(int(n)) +'\t'+alpha_out+'\n')
    outfile.close()