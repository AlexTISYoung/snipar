#!/apps/well/python/2.7.8/bin/python
import argparse
import numpy as np
from pysnptools.snpreader import Bed, Pheno
from scipy.stats import chi2, zscore
from sibreg import sibreg
import code

####### Output functions ##########
def neglog10pval(x,df):
    return -np.log10(np.e)*chi2.logsf(x,df)

def vector_out(alpha_mle,digits=6):
##Output parameter estimates along with standard errors ##
    ## Calculate test statistics
    alpha_est = alpha_l[0][1:3]
    alpha_cov = alpha_l[1][1:3,1:3]
    alpha_ses = np.sqrt(np.diag(alpha_cov))
    alpha_out = str(round(alpha_est[0],digits))+'\t'+str(round(alpha_ses[0],digits))+'\t'
    alpha_out += str(round(alpha_est[1],digits))+'\t'+str(round(alpha_ses[1],digits))+'\t'
    alpha_out += str(round(alpha_cov[0,1]/(alpha_ses[0]*alpha_ses[1]),digits))+'\n'
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
    parser.add_argument('genofile',type=str,help='Path to genotypes in BED format')
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
    parser.add_argument('--fix_VC', action='store_true', default=False,
                        help='Fix the variance components to the values from the null model')
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
    test_chr = Bed(args.genofile)
    #test_chr = Bed('full_sib_genotypes/chr_22.bed')
    # select subset to test
    iid = test_chr.iid
    sid = test_chr.sid
    pos = test_chr.pos
    # Find ids from same families as those with phenotype data
    fid_with_phen = np.array([x in pheno_fam_id_set for x in iid[:,0]])
    test_chr = test_chr[fid_with_phen,:].read()
    genotypes = test_chr.val
    # Get genotype matrix
    if genotypes.ndim == 1:
        chr_length = 1
        genotypes = genotypes.reshape(genotypes.shape[0], 1)
    else:
        chr_length = genotypes.shape[1]
    print('Number of test loci: ' + str(genotypes.shape[1]))
    print('Genotypes for '+str(genotypes.shape[0])+' individuals read')
    # Get sample ids
    geno_id_dict = id_dict_make(np.array(test_chr.iid))

### Intersect and match genotype data with covariate and phenotype data ###
    # Intersect with phenotype IDs
    ids_in_common = {tuple(x) for x in pheno_ids} & geno_id_dict.viewkeys()
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
    # Vector to match genotype with phenotype ids
    geno_id_match = np.array([geno_id_dict[tuple(x)] for x in pheno_ids])
    # Dictionary to look up rows for each family
    geno_fam_dict = {}
    for fam in np.unique(pheno_ids[:,0]):
        geno_fam_dict[fam] = np.where(test_chr.iid[:,0]==fam)[0]

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
        header='SNP\tfrequency\tn\tWF\tWF_se\tBF\tBF_se\tr_WF_BF\n'
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
    for loc in xrange(0,chr_length):
        print(str(loc))
        alpha_out = 'NA\tNA\tNA\tNA\tNA\tNA\n'
        # Filler for output if locus doesn't pass thresholds
        allele_frq=np.nan
        # Get test genotypes
        test_gts=genotypes[:,loc]
        # Find missingness and allele freq
        test_gt_not_na=np.logical_not(np.isnan(test_gts))
        n_l=np.sum(test_gt_not_na)
        # Threshold on missingness
        missingness = 100.0 * (1 - float(n_l) / n)
        if missingness<args.max_missing:
            t=test_gts[test_gt_not_na]
            allele_frq=np.mean(t)/2
            if allele_frq>0.5:
                allele_frq=1-allele_frq
            if allele_frq>args.min_maf:
                del t
                # Compute within family mean genotypes
                g_mean = np.zeros((y.shape[0]))
                g_mean[:] = np.nan
                families = np.unique(pheno_ids[:,0])
                for fam in families:
                    g_fam = genotypes[geno_fam_dict[fam],loc]
                    g_fam_not_NA = np.logical_not(np.isnan(g_fam))
                    if np.sum(g_fam_not_NA) > 1:
                        g_mean[pheno_fam_dict[fam]] = np.mean(g_fam[g_fam_not_NA])
                # Get proband genotypes
                test_gts = test_gts[geno_id_match]
                # Remove those with missing proband genotypes and families with less than two observed sib genotypes
                not_na = np.logical_and(np.logical_not(np.isnan(test_gts)),np.logical_not(np.isnan(g_mean)))
                test_gts = test_gts[not_na]
                y_l = y[not_na]
                g_mean = g_mean[not_na]
                n_loc = g_mean.shape[0]
                fam_l = pheno_ids[not_na,0]
                # Optimize model for SNP
                X_l = np.ones((y_l.shape[0],3))
                X_l[:,1] = test_gts-g_mean
                X_l[:,2] = g_mean
                model_l = sibreg.model(y_l,X_l,fam_l)
                if not args.fix_VC:
                    optim_l = model_l.optimize_model(np.array([null_optim['sigma2'],null_optim['tau']]))
                    if optim_l['success']:
                        alpha_l = model_l.alpha_mle(optim_l['tau'],optim_l['sigma2'],compute_cov = True)
                    else:
                        print('Maximisation of likelihood failed for for ' + sid[loc])
                else:
                    alpha_l = model_l.alpha_mle(null_optim['tau'],null_optim['sigma2'],compute_cov = True)
                alpha_out = str(n_loc) + '\t' + vector_out(alpha_l)
            print('finished successfully')
        outfile.write(sid[loc] +'\t'+ str(allele_frq)+'\t'+alpha_out+'\n')
    outfile.close()