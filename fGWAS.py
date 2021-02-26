from sibreg.sibreg import *
import argparse
import numpy as np

def transform_phenotype(inv_root,y, fam_indices, null_mean = None):
    """
    Transform phenotype based on inverse square root of phenotypic covariance matrix.
    If the null model included covariates, the fitted mean is removed rather than the overall mean
    """
    # Mean normalise phenotype
    if null_mean is None:
        y = y - np.mean(y)
    else:
        y = y - null_mean
    # Transform by family
    for fam in fam_indices.keys():
        famsize = fam_indices[fam].shape[0]
        if famsize == 1:
            y[fam_indices[fam]] = inv_root[1] * y[fam_indices[fam]]
        else:
            y[fam_indices[fam]] = inv_root[famsize].dot(y[fam_indices[fam]])
    return y

def fit_models(y,G):
    """
    Perform repeated OLS to estimate SNP effects and sampling variance-covariance in transformed model
    """
    G.gts = G.gts.transpose(2,0,1)
    XTX = np.einsum('...ij,...ik', G.gts, G.gts)
    XTY = np.einsum('...ij,i',G.gts,y)
    alpha = np.linalg.solve(XTX,XTY)
    alpha_cov = np.linalg.inv(XTX)
    alpha_ses = np.sqrt(np.diagonal(alpha_cov,axis1=1,axis2=2))
    return alpha, alpha_cov, alpha_ses

def write_output(G, outprefix, parsum, sib, alpha, alpha_ses, alpha_cov, sigma2, tau, NAs, null_alpha = None):
    """
    Write fitted SNP effects and other parameters to output HDF5 file.
    """
    print('Writing output to ' + outprefix + '.hdf5')
    outfile = h5py.File(outprefix + '.hdf5', 'w')
    outbim = np.column_stack((G.chrom,G.sid,G.pos,G.alleles))
    outfile['bim'] = encode_str_array(outbim)
    X_length = 1
    outcols = ['direct']
    if sib:
        X_length += 1
        outcols.append('sib')
    if parsum:
        X_length += 1
        outcols.append('avg_parental')
    else:
        X_length += 2
        outcols = outcols + ['paternal','maternal']
    outfile.create_dataset('estimate_covariance', (G.sid.shape[0], X_length, X_length), dtype='f', chunks=True,
                           compression='gzip', compression_opts=9)
    outfile.create_dataset('estimate', (G.sid.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile.create_dataset('estimate_ses', (G.sid.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile['estimate'][:] = alpha
    outfile['estimate_cols'] = encode_str_array(np.array(outcols))
    outfile['estimate_ses'][:] = alpha_ses
    outfile['estimate_covariance'][:] = alpha_cov
    outfile['sigma2'] = sigma2
    outfile['tau'] = tau
    outfile['N'] = G.gts.shape[1]
    outfile['NAs'] = NAs
    outfile['freqs'] = G.freqs
    if null_alpha is not None:
        null_ses = np.sqrt(np.diag(null_alpha[1]))
        null_out = np.column_stack((null_alpha[0],null_ses))
        outfile['covariate_estimates'] = null_out
    outfile.close()

######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('pargts', type=str, help='HDF5 file with imputed parental genotypes (without .hdf5 suffix)')
    parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
    parser.add_argument('outprefix',type=str,help='Location to output association statistic hdf5 file')
    parser.add_argument('--bed',type=str,help='Bed file with observed genotypes (without .bed suffix).',default=None)
    parser.add_argument('--bgen',type=str,help='bgen file with observed genotypes (without .bgen suffix).',default=None)
    parser.add_argument('--parsum',action='store_true',help='Regress onto proband and sum of parental genotypes (useful when parental genotypes imputed from sibs only)',default = False)
    parser.add_argument('--fit_sib',action='store_true',help='Fit indirect effect from sibling ',default=False)
    parser.add_argument('--covar',type=str,help='Path to file with covariates: plain text file with columns FID, IID, covar1, covar2, ..', default=None)
    parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                        default=1)
    parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)', default=0.01)
    parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)', default=5)
    parser.add_argument('--min_info',type=float,help='Ignore SNPs with imputation INFO score below this threshold (default 0.99)', default=0.99)
    parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)', default='NA')
    parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                        default=1)
    parser.add_argument('--output_covar_ests',action='store_true',help='Output null model estimates of covariate fixed effects (default False)', default=False)
    args=parser.parse_args()

    ######### Read Phenotype ########
    y, pheno_ids = read_phenotype(args.phenofile, missing_char=args.missing_char, phen_index=args.phen_index)
    ######## Read covariates ########
    if args.covar is not None:
        print('Reading covariates')
        covariates = read_covariates(args.covar, missing_char=args.missing_char)
        # Match to pheno ids
        covariates.filter_ids(pheno_ids)
    ######## Check for bed/bgen #######
    if args.bed is None and args.bgen is None:
        raise(ValueError('Must supply either bed or bgen file with observed genotypes'))
    if args.bed is not None and args.bgen is not None:
        raise(ValueError('Both --bed and --bgen specified. Please specify one only'))
    if args.bed is not None:
        gts_f = args.bed+'.bed'
    elif args.bgen is not None:
        gts_f = args.bgen+'.bgen'
    ####### Construct family based genotype matrix #######
    G = get_gts_matrix(args.pargts+'.hdf5', gts_f, ids=pheno_ids, parsum=args.parsum, sib=args.fit_sib)
    # Check for empty fam labels
    no_fam = np.array([len(x) == 0 for x in G.fams])
    if np.sum(no_fam) > 0:
        ValueError('No family label from pedigree for some individuals')
    #### Filter SNPs ####
    print('Filtering based on MAF')
    G.filter_maf(args.min_maf)
    if args.bed is not None:
        print('Filtering based on missingness')
        G.filter_missingness(args.max_missing)
    if args.bgen is not None:
        print('Filtering based on imputation INFO')
        G.filter_info(args.min_info)
    print(str(G.shape[2])+' SNPs that pass filters')
    #### Fill NAs ####
    print('Imputing missing values with population frequencies')
    NAs = G.fill_NAs()
    #### Match phenotype ####
    y = match_phenotype(G,y,pheno_ids)
    #### Fit null model ####
    print('Estimating variance components')
    if args.covar is not None:
        # Match covariates #
        covariates.filter_ids(G.ids, verbose=False)
        null_model, sigma2, tau, null_alpha, null_alpha_cov = fit_sibreg_model(y, covariates.gts, G.fams, add_intercept=True,
                                                   tau_init=args.tau_init)

    else:
        null_model, sigma2, tau  = fit_sibreg_model(y, np.ones((y.shape[0], 1)), G.fams,
                                                                               tau_init = args.tau_init, return_fixed = False)
    print('Family variance estimate: '+str(round(sigma2/tau,4)))
    print('Residual variance estimate: ' + str(round(sigma2,4)))
    ##### Transform genotypes and phenotypes ######
    print('Transforming genotypes and phenotypes')
    L = null_model.sigma_inv_root(tau, sigma2)
    G.diagonalise(L)
    if args.covar is None:
        y = transform_phenotype(L, y, G.fam_indices)
    else:
        null_mean = null_alpha[0]+covariates.gts.dot(null_alpha[1:null_alpha.shape[0]])
        y = transform_phenotype(L, y, G.fam_indices, null_mean)
    ### Fit models for SNPs ###
    print('Estimating SNP effects')
    alpha, alpha_cov, alpha_ses = fit_models(y,G)
    ### Save output ###
    if args.output_covar_ests and args.covar is not None:
        write_output(G, args.outprefix, args.parsum, args.fit_sib, alpha, alpha_ses, alpha_cov, sigma2, tau, NAs, [null_alpha,null_alpha_cov])
    else:
        write_output(G, args.outprefix, args.parsum, args.fit_sib, alpha, alpha_ses, alpha_cov, sigma2, tau, NAs)