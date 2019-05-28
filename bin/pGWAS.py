#!/apps/well/python/2.7.8/bin/python
import numpy as np
import numpy.ma as ma
from pysnptools.snpreader import Bed, Pheno
from scipy.stats import chi2, zscore
from sibreg import sibreg
import h5py, argparse, code

####### Output functions ##########
def neglog10pval(x,df):
    return -np.log10(np.e)*chi2.logsf(x,df)

def vector_out(n,alpha_l, no_sib, n_X, digits=6):
    ##Output parameter estimates along with standard errors ##
    ## Calculate test statistics
    if no_sib:
        X_length = n_X+2
    else:
        X_length = n_X+3
    alpha_est = alpha_l[0][n_X:X_length]
    alpha_cov = alpha_l[1][n_X:X_length,n_X:X_length]
    alpha_ses = np.sqrt(np.diag(alpha_cov))
    alpha_corr = np.dot(np.diag(1/alpha_ses).dot(alpha_cov),np.diag(1/alpha_ses))
    alpha_out = str(n)+'\t'+str(round(alpha_est[0],digits))+'\t'+str(round(alpha_ses[0],digits))+'\t'
    alpha_out += str(round(alpha_est[1],digits))+'\t'+str(round(alpha_ses[1],digits))+'\t'
    if no_sib:
        alpha_out += str(round(alpha_corr[0, 1], digits)) + '\n'
    else:
        alpha_out += str(round(alpha_est[2],digits))+'\t'+str(round(alpha_ses[2],digits))+'\t'
        alpha_out += str(round(alpha_corr[0,1],digits))+'\t'+str(round(alpha_corr[0,2],digits))+'\t'+str(round(alpha_corr[2,1],digits))+'\n'
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
    parser.add_argument('sibgts',type=str,help='Path to bed file with sibling genotypes')
    parser.add_argument('pargts', type=str, help='Path to HDF5 file with imputed parental genotypes')
    parser.add_argument('sibped',type=str,help='Path to pedigree file with siblings sharing a family ID and non-siblings not')
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
    parser.add_argument('--fix_VC', action='store_true', default=False,
                        help='Fix the variance components to the values from the null model')
    args=parser.parse_args()

    ####################### Read in data #########################
    #### Read phenotype ###
    pheno = Pheno(args.phenofile, missing=args.missing_char).read()
    #pheno = Pheno('phenotypes/eduyears_resid.ped', missing='NA').read()
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
    pheno_id_dict = {}
    for i in xrange(0,y.shape[0]):
        pheno_id_dict[pheno_ids[i,1]] = i
    pheno_fams = set(pheno_ids[:,0])
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

### Read pedigree file ###
    ### Load pedigree
    ped = np.loadtxt(args.sibped, dtype='S20', skiprows=1)
    #ped = np.loadtxt('23andme/sibs.ped', dtype='S20', skiprows=1)
    sibs = np.unique(ped[:, 1])

    ### Create family dictionary
    fams = {}
    fam_ids = np.unique(ped[:, 0])
    for f in fam_ids:
        fams[f] = tuple(ped[ped[:, 0] == f, 1])
    # reverse lookup dict
    sib_fam_dict = {}
    for i in xrange(0, ped.shape[0]):
        sib_fam_dict[ped[i, 1]] = ped[i, 0]

### Read sibling genotypes ###
    #### Load genotypes
    gts_f = Bed(args.sibgts)
    #gts_f = Bed('genotypes/chr_22.bed')
    gts_ids = gts_f.iid
    # Build dict
    id_dict = {}
    for i in xrange(0, gts_ids.shape[0]):
        id_dict[gts_ids[i, 1]] = i

    # Find sibling indices
    sib_indices = []
    sibs_new = sibs
    for i in xrange(0, sibs.shape[0]):
        s = sibs[i]
        if s in id_dict:
            sib_indices.append(id_dict[s])
        else:
            print('No genotype data for ' + str(s))
            fam_i = sib_fam_dict[s]
            sfam = fams[fam_i]
            # Remove family if only two sibs and one missing
            if len(sfam) == 2:
                sindices = np.array([np.where(sibs == x)[0][0] for x in sfam])
                sibs_new = np.delete(sibs_new, sindices)
                fam_ids = np.delete(fam_ids, np.where(fam_ids == fam_i)[0][0])
                del fams[fam_i]
            # otherwise remove sib and keep fam
            else:
                sibs_new = np.delete(sibs, i)
                sfam = np.delete(sfam, np.where(sfam == s))
                fams[fam_i] = sfam

    sibs = sibs_new
    sib_indices = np.sort(sib_indices)

    # Read sibling genotypes
    gts = gts_f[sib_indices, :].read().val
    pos = gts_f.pos[:, 2]
    sid = gts_f.sid
    sid_dict = {}
    for i in range(0,sid.shape[0]):
        sid_dict[sid[i]] = i
    gts = ma.array(gts,mask=np.isnan(gts),dtype=int)

    # rebuild ID dictionary
    gts_ids = gts_ids[sib_indices]
    # Build dict
    id_dict = {}
    for i in xrange(0, gts_ids.shape[0]):
        id_dict[gts_ids[i, 1]] = i

    ### Read imputed parental genotypes ###
    print('Reading imputed parental genotype file')
    pargts_f = h5py.File(args.pargts,'r')
    #pargts_f = h5py.File('23andme/chr_22.hdf5','r')
    # get families
    par_fams = np.array(pargts_f['families'])
    # build family dictionary
    par_fam_dict = {}
    for i in range(0,par_fams.shape[0]):
        par_fam_dict[par_fams[i]] = i
    pargts = np.array(pargts_f['imputed_par_gts'])
    par_sid = np.array(pargts_f['sid'])
    par_sid_dict = {}
    for i in range(0,par_sid.shape[0]):
        par_sid_dict[par_sid[i]] = i

    ### Match SIDs of sibling and par gts ###
    in_sib_sid = np.zeros((par_sid.shape[0]),dtype=bool)
    sib_sid_indices = []
    for s in range(0,par_sid.shape[0]):
        sid_s = par_sid[s]
        if sid_s in sid_dict:
            in_sib_sid[s] = True
            sib_sid_indices.append(sid_dict[sid_s])
    if np.sum(in_sib_sid)>0:
        pargts = pargts[:,in_sib_sid]
        par_sid = par_sid[in_sib_sid]
        gts = gts[:,sib_sid_indices]
        print(str(gts.shape[1])+' variants in common between parental and sibling genotypes')
    else:
        raise(ValueError('No variants in common between sibling and parental genotypes'))

### Construct genetic covariate matrix
    # Find families with phenotype data, covariate data, at least two siblings genotyped, and
    fams_with_data = []
    sibs_with_pheno = []
    sibs_with_geno = []
    sample_size = 0
    for f in fam_ids:
        if f in par_fams and fams and pheno_fams:
            sibs = np.array(fams[f])
            sibs_in_gts = np.array([x in id_dict for x in sibs])
            if np.sum(sibs_in_gts) > 1:
                sibs_with_geno.append(sibs[sibs_in_gts])
                sibs_in_pheno = np.array([x in pheno_id_dict for x in sibs])
                sibs_in_pheno = np.logical_and(sibs_in_pheno,sibs_in_gts)
                nphen = np.sum(sibs_in_pheno)
                if nphen > 0:
                    fams_with_data.append(f)
                    sibs_with_pheno.append(sibs[sibs_in_pheno])
                    sample_size += nphen

    nfam = len(fams_with_data)
    print('Sample of '+str(sample_size)+' comprised of '+str(nfam)+' families with at least two siblings genotyped and at least one sibling phenotyped')

    print('Forming family-wise genotype matrix')
    if args.no_sib:
        gsize = 1
    else:
        gsize = 2
    G = ma.array(np.zeros((sample_size, gsize, gts.shape[1]),dtype=np.int8),
                 mask=np.zeros((sample_size, gsize, gts.shape[1]), dtype=bool))
    G_par = np.zeros((sample_size, gts.shape[1]), dtype=np.float32)
    start = 0
    y_new = np.zeros((sample_size))
    X_new = np.zeros((sample_size,X.shape[1]))
    fam_labels = np.zeros((sample_size),dtype='S20')
    for i in xrange(0,nfam):
        fam = fams_with_data[i]
        end = start + len(sibs_with_pheno[i])
        # Fill in family labels vector
        fam_labels[start:end] = fam
        # Fill in phenotype vector
        pheno_indices_i = np.array([pheno_id_dict[x] for x in sibs_with_pheno[i]])
        y_new[start:end] = y[pheno_indices_i]
        # Fill in phenotype vector
        X_new[start:end,:] = X[pheno_indices_i,:]
        # Fill in proband genotype column
        gindices = np.array([id_dict[x] for x in sibs_with_pheno[i]])
        G[start:end,0,:] = gts[gindices,:]
        G.mask[start:end,0,:] = gts.mask[gindices,:]
        # Fill in parental genotype column
        par_index_i = par_fam_dict[fams_with_data[i]]
        for j in xrange(start,end):
            G_par[j,:] = pargts[par_index_i,:]
        if not args.no_sib:
            # Fill in average sibling genotype column
            gindices = np.array([id_dict[x] for x in sibs_with_geno[i]])
            gsum = ma.sum(gts[gindices,:],axis=0)
            gmask = np.sum(gts[gindices, :].mask, axis=0) > 0
            for s in range(0,sibs_with_pheno[i].shape[0]):
                sname = sibs_with_pheno[i][s]
                G[start+s,1,:] = gsum-gts[id_dict[sname],:]
                G[start+s,1,:].mask = gmask
            G[start:end,1,:] = G[start:end,1,:]/float(sibs_with_geno[i].shape[0]-1)
        start = end

    y = y_new
    X = X_new
    families = np.array(fams_with_data)

######### Initialise output files #######
    ## Output file
    if args.append:
        write_mode='ab'
    else:
        write_mode='wb'
    outfile=open(args.outprefix+'.models.gz',write_mode)
    #outfile = open('23andme/test.models.gz', 'wb')
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
    null_model = sibreg.model(y, X, fam_labels)
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
        X=np.ones((sample_size,1))
        n_X=1

    ############### Loop through loci and fit models ######################
    print('Fitting models for genome-wide SNPs')
    # Optimize model for SNP
    if args.no_sib:
        X_length = n_X + 2
    else:
        X_length = n_X + 3
    freqs = ma.mean(G[:,0,:],axis=0)/2.0
    missingness = ma.mean(G.mask[:,0,:],axis=0)
    #for loc in xrange(0,G.shape[2]):
    for loc in xrange(41,42):
        if args.no_sib:
            alpha_out = 'NA\tNA\tNA\tNA\tNA\tNA\n'
        else:
            alpha_out = 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n'
        if freqs[loc] > args.min_maf and freqs[loc] < (1-args.min_maf) and missingness[loc] < args.max_missing:
            # Find NAs
            if args.no_sib:
                not_nans = np.logical_not(G[:,0,loc].mask)
            else:
                not_nans = np.logical_not(G[:,1,loc].mask)
            n_l = np.sum(not_nans)
            X_l = np.ones((n_l, X_length),dtype=np.float32)
            X_l[:, n_X:(X_length-1)] = G[not_nans, :, loc]
            X_l[:,X_length-1] = G_par[not_nans,loc]
            model_l = sibreg.model(y[not_nans],X_l,fam_labels[not_nans])
            if not args.fix_VC:
                optim_l = model_l.optimize_model(np.array([null_optim['sigma2'], null_optim['tau']]))
                if optim_l['success']:
                    alpha_l = model_l.alpha_mle(optim_l['tau'], optim_l['sigma2'], compute_cov=True)
                else:
                    print('Maximisation of likelihood failed for for ' + sid[loc])
            else:
                alpha_l = model_l.alpha_mle(null_optim['tau'], null_optim['sigma2'], compute_cov=True)
            alpha_out = vector_out(n_l,alpha_l, args.no_sib, n_X)
            code.interact(local=locals())
        outfile.write(sid[loc] +'\t'+str(freqs[loc])+'\t'+alpha_out)
    outfile.close()