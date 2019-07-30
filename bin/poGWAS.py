#!/well/kong/users/wiw765/anaconda2/bin/python
import numpy as np
import numpy.ma as ma
from pysnptools.snpreader import Bed, Pheno
from scipy.stats import chi2, zscore
from sibreg import sibreg
import h5py, argparse, code

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
    id_dict = {}
    for i in range(0,ids.shape[0]):
        id_dict[ids[i,1]] = i
    # Match with pheno_ids
    common_ids = id_dict.viewkeys() & set(ids_to_match[:,1])
    pheno_in = np.array([x in common_ids for x in ids_to_match[:,1]])
    match_ids = ids_to_match[pheno_in,1]
    X_id_match = np.array([id_dict[x] for x in match_ids])
    X = X[X_id_match, :]
    return [X,X_names,pheno_in]

######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('gts',type=str,help='Path to bed file with sibling genotypes')
    parser.add_argument('pargts', type=str, help='Path to HDF5 file with imputed parental genotypes')
    parser.add_argument('ped',type=str,help='Path to pedigree file')
    parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
    parser.add_argument('outprefix',type=str,help='Location to output csv file with association statistics')
    parser.add_argument('--covar',type=str,help='Location of covariate file (default None)',
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
    parser.add_argument('--no_sib',action='store_true',default=False,help='Remove individuals with genotyped siblings and do not fit indirect effects from sibs')
    parser.add_argument('--fit_sib',action='store_true',default=False,help='Only analyse individuals with genotyped siblings and fit indirect effects from sibs')
    parser.add_argument('--fit_VC', action='store_true', default=False,
                        help='Fit the variance components for each SNP (default is to use null model MLE)')
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
    if not args.covar == None:
        X, X_names, pheno_in = read_covariates(args.covar,pheno_ids, args.missing_char)
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
    ### Read pedigree file ###
    ### Load pedigree
    ped = np.loadtxt(args.ped, dtype='S20', skiprows=1)
    # ped = np.loadtxt('relatedness/families.ped', dtype='S20', skiprows=1)

    ### Read imputed parental genotypes ###
    print('Reading imputed parental genotype file')
    pargts_f = h5py.File(args.pargts,'r')
    #pargts_f = h5py.File('one_parent_genotyped/imputed/chr_22.hdf5','r')
    # get families
    par_ped = np.array(pargts_f['ped'])
    # build dictionary
    par_dict = {}
    for i in range(0,par_ped.shape[0]):
        par_dict[par_ped[i,1]] = i
    pargts = np.array(pargts_f['imputed_par_gts'])
    par_sid = np.array(pargts_f['sid'])
    par_sid_dict = {}
    for i in range(0,par_sid.shape[0]):
        par_sid_dict[par_sid[i]] = i

    father_genotyped = np.array(pargts_f['father_genotyped'])
    pargts_f.close()

    ### Read genotypes ###
    #### Load genotypes
    gts_f = Bed(args.gts)
    # gts_f = Bed('genotypes/chr_22.bed')
    gts_ids = gts_f.iid
    # Build dict
    id_dict = {}
    for i in xrange(0, gts_ids.shape[0]):
        id_dict[gts_ids[i, 1]] = i

    genotype_indices = []
    remove = []
    for i in range(0, par_ped.shape[0]):
        in_genotype = np.zeros((ped.shape[1]),dtype=bool)
        in_genotype[1:4] = np.array([x in id_dict for x in par_ped[i,1:4]])
        if in_genotype[1] and np.sum(in_genotype[2:4])==1:
            # Check for siblings
            genotype_indices = genotype_indices+[id_dict[x] for x in par_ped[i,in_genotype]]
        else:
            remove.append(i)

    # Remove rows without genotype data
    if len(remove)>0:
        remove = np.array(remove)
        par_ped = np.delete(par_ped,remove,axis=0)
        pargts = np.delete(pargts,remove,axis=0)
        father_genotyped = np.delete(father_genotyped,remove)

    # check for individuals with genotyped siblings
    genotype_indices = np.sort(np.unique(np.array(genotype_indices)))

    # Read genotypes
    gts = gts_f[genotype_indices, :].read().val
    pos = gts_f.pos[:, 2]
    sid = gts_f.sid
    gts = ma.array(gts, mask=np.isnan(gts), dtype=int)

    # rebuild ID dictionary
    gts_ids = gts_ids[genotype_indices,:]
    # Build dict
    id_dict = {}
    for i in xrange(0, gts_ids.shape[0]):
        id_dict[gts_ids[i, 1]] = i

    # remove/keep individuals on basis of genotyped siblings
    if args.no_sib or args.fit_sib:
        # Count number of siblings
        sibcount = np.zeros((par_ped.shape[0]),dtype=int)
        if args.fit_sib:
            sibs_list = []
        for i in xrange(0,par_ped.shape[0]):
            sibs_i = np.logical_and(par_ped[:,0]==par_ped[i,0],np.logical_and(par_ped[:,2]==par_ped[i,2],par_ped[:,3]==par_ped[i,3]))
            sibs_i[i] = False
            sibcount[i] = np.sum(sibs_i)
            if sibcount[i]>0 and args.fit_sib:
                sibs_list.append(par_ped[sibs_i,1])

        if args.no_sib:
            no_sibs = sibcount==0
            print('Removing ' + str(np.sum(sibcount > 0)) + ' individuals with genotyped siblings')
            par_ped = par_ped[no_sibs,:]
            pargts = pargts[no_sibs,:]
            father_genotyped = father_genotyped[no_sibs]

        if args.fit_sib:
            has_sibs = sibcount>0
            print('Removing '+str(np.sum(sibcount==0))+' individuals without genotyped siblings')
            par_ped = par_ped[has_sibs,:]
            pargts = pargts[has_sibs,:]
            father_genotyped = father_genotyped[has_sibs]


    ### Match SIDs of sibling and par gts ###
    in_gts_sid = np.zeros((par_sid.shape[0]),dtype=bool)
    gts_sid_indices = []
    for s in range(0,par_sid.shape[0]):
        sid_s = par_sid[s]
        if sid_s in par_sid_dict:
            in_gts_sid[s] = True
            gts_sid_indices.append(par_sid_dict[sid_s])
    if np.sum(in_gts_sid)>0:
        pargts = pargts[:,in_gts_sid]
        par_sid = par_sid[in_gts_sid]
        gts = gts[:,gts_sid_indices]
        print(str(gts.shape[1])+' variants in common between parental and sibling genotypes')
    else:
        raise(ValueError('No variants in common between sibling and parental genotypes'))
    sid = sid[gts_sid_indices]

### Construct genetic covariate matrix
    # Find parent-offspring pairs with genotype data and phenotyped offspring
    print('Forming family-wise genotype matrix')
    if args.fit_sib:
        G_plus = 1
    else:
        G_plus = 0
    G = ma.array(np.zeros((par_ped.shape[0], 3+G_plus, gts.shape[1]),dtype=np.float32),
                 mask=np.zeros((par_ped.shape[0], 3+G_plus,gts.shape[1]), dtype=bool))
    y_new = np.zeros((par_ped.shape[0]))
    y_new[:] = np.nan
    X_new = np.zeros((par_ped.shape[0],X.shape[1]))
    X_new[:] = np.nan
    for i in range(0,par_ped.shape[0]):
        # get phenotype
        if par_ped[i,1] in pheno_id_dict:
            yindex = pheno_id_dict[par_ped[i,1]]
            y_new[i] = y[yindex]
            X_new[i,:] = X[yindex,:]
        # get individual genotype
        G[i,0,:] = gts[id_dict[par_ped[i, 1]],:]
        # If fitting sib effects, get sib genotypes
        if args.fit_sib:
            G[i,1,:] = ma.mean(gts[np.array([id_dict[x] for x in sibs_list[i]]),:],axis=0)
        # get parental genotypes
        if father_genotyped[i]:
            G[i,1+G_plus,:] = gts[id_dict[par_ped[i, 2]],:]
            G[i,2+G_plus,:] = pargts[i,:]
        else:
            G[i, 1+G_plus, :] = pargts[i,:]
            G[i, 2+G_plus, :] = gts[id_dict[par_ped[i, 3]], :]

    del gts, pargts

    # Form phenotype
    y = y_new
    X = X_new

    y_not_nan = np.logical_not(np.isnan(y))
    G = G[y_not_nan,:]
    y = y[y_not_nan]
    X = X[y_not_nan,:]
    par_ped = par_ped[y_not_nan,:]

    print(str(G.shape[0])+' individuals with one parent genotyped and observed phenotype')

######### Fit Null Model ##########
    ## Get initial guesses for null model
    print('Fitting Null Model')
    # Optimize null model
    sigma_2_init = np.var(y)*args.tau_init/(1+args.tau_init)
    #sigma_2_init = np.var(y) * 1 / (1 + 1)
    null_model = sibreg.model(y, X, par_ped[:,0])
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
    if not args.append and not args.no_covariate_estimates and args.covar is not None:
        np.savetxt(args.outprefix + '.null_covariate_effects.txt',
                   np.hstack((X_names.reshape((n_X, 1)), np.array(alpha_out, dtype='S20'))),
                   delimiter='\t', fmt='%s')

    # Fit SNP specific models
    ### Project out mean covariates
    if not args.fit_covariates:
        # Residual y
        y=y-X.dot(null_alpha[0])
        # Reformulate fixed_effects
        X=np.ones((X.shape[0],1))
        n_X=1

    ######### Initialise output files #######
    ## Output file
    outfile = h5py.File(args.outprefix + '.hdf5', 'w')
    outfile['sid'] = sid
    X_length = n_X + 3 + G_plus
    outfile.create_dataset('xtx', (G.shape[2], X_length, X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile.create_dataset('xty', (G.shape[2], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)


    ############### Loop through loci and fit models ######################
    print('Fitting models for genome-wide SNPs')
    # Optimize model for SNP
    freqs = ma.mean(G[:,0,:],axis=0)/2.0
    missingness = ma.mean(G.mask[:,0,:],axis=0)
    N_L = np.zeros((G.shape[2]),dtype = int)
    for loc in xrange(0,G.shape[2]):
        if freqs[loc] > args.min_maf and freqs[loc] < (1-args.min_maf) and (100*missingness[loc]) < args.max_missing:
            # Find NAs
            not_nans = np.sum(G[:,:,loc].mask,axis = 1)==0
            n_l = np.sum(not_nans)
            N_L[loc] = n_l
            X_l = np.ones((n_l, X_length),dtype=np.float32)
            X_l[:,0:n_X] = X[not_nans,:]
            X_l[:, n_X:X_length] = G[not_nans, :, loc]
            model_l = sibreg.model(y[not_nans],X_l,par_ped[not_nans,0])
            if args.fit_VC:
                optim_l = model_l.optimize_model(np.array([null_optim['sigma2'], null_optim['tau']]))
                if optim_l['success']:
                    alpha_l = model_l.alpha_mle(optim_l['tau'], optim_l['sigma2'], compute_cov=True, xtx_out= True)
                else:
                    raise(ValueError('Maximisation of likelihood failed for for ' + sid[loc]))
            else:
                alpha_l = model_l.alpha_mle(null_optim['tau'], null_optim['sigma2'], compute_cov=True, xtx_out= True)
            outfile['xtx'][loc,:,:] = alpha_l[0]
            outfile['xty'][loc,:] = alpha_l[1]
        else:
            outfile['xtx'][loc, :, :] = np.nan
            outfile['xty'][loc, :] = np.nan
    outfile['sigma2'] = null_optim['sigma2']
    outfile['tau'] = null_optim['tau']
    outfile['N_L'] = N_L
    outfile['father_genotyped'] = father_genotyped
    outfile.close()