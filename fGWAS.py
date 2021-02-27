from sibreg.sibreg import *
import argparse, h5py
import numpy as np
from bgen_reader import open_bgen

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

def write_output(chrom, snp_ids, pos, alleles, outprefix, parsum, sib, alpha, alpha_ses, alpha_cov, sigma2, tau, freqs):
    """
    Write fitted SNP effects and other parameters to output HDF5 file.
    """
    print('Writing output to ' + outprefix + '.hdf5')
    outfile = h5py.File(outprefix + '.hdf5', 'w')
    outbim = np.column_stack((chrom,snp_ids,pos,alleles))
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
    outfile.create_dataset('estimate_covariance', (snp_ids.shape[0], X_length, X_length), dtype='f', chunks=True,
                           compression='gzip', compression_opts=9)
    outfile.create_dataset('estimate', (snp_ids.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile.create_dataset('estimate_ses', (snp_ids.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile['estimate'][:] = alpha
    outfile['estimate_cols'] = encode_str_array(np.array(outcols))
    outfile['estimate_ses'][:] = alpha_ses
    outfile['estimate_covariance'][:] = alpha_cov
    outfile['sigma2'] = sigma2
    outfile['tau'] = tau
    outfile['freqs'] = freqs
    outfile.close()

def compute_batch_boundaries(snp_ids,batch_size):
    nsnp = snp_ids.shape[0]
    n_blocks = np.int(np.ceil(float(nsnp)/float(batch_size)))
    block_bounds = np.zeros((n_blocks,2),dtype=int)
    start = 0
    for i in range(n_blocks-1):
        block_bounds[i,0] = start
        block_bounds[i,1] = start+batch_size
        start += batch_size
    block_bounds[n_blocks-1,:] = np.array([start,nsnp])
    return block_bounds

def process_batch(snp_ids, y, pheno_ids, pargts_f, gts_f, fit_null=False, tau=None, sigma2=None, null_alpha=None, covar=None, parsum=False,
                  fit_sib=False, max_missing=5, min_maf=0.01, min_info=0.9, tau_init=1, verbose = False):
    ####### Construct family based genotype matrix #######
    G = get_gts_matrix(pargts_f+'.hdf5', gts_f, snp_ids=snp_ids, ids=pheno_ids, parsum=parsum, sib=fit_sib)
    # Check for empty fam labels
    no_fam = np.array([len(x) == 0 for x in G.fams])
    if np.sum(no_fam) > 0:
        ValueError('No family label from pedigree for some individuals')
    G.compute_freqs()
    chrom, pos, alleles, freqs = G.chrom,  G.pos, G.alleles, G.freqs
    #### Filter SNPs ####
    if verbose:
        print('Filtering based on MAF')
    G.filter_maf(min_maf)
    gt_filetype = gts_f.split('.')[1]
    if gt_filetype=='bed':
        if verbose:
            print('Filtering based on missingness')
        G.filter_missingness(max_missing)
    if gt_filetype=='bgen':
        if verbose:
            print('Filtering based on imputation INFO')
        G.filter_info(min_info)
    if verbose:
        print(str(G.shape[2])+' SNPs that pass filters')
    #### Fill NAs ####
    if verbose:
        print('Imputing missing values with population frequencies')
    NAs = G.fill_NAs()
    #### Match phenotype ####
    y = match_phenotype(G,y,pheno_ids)
    #### Fit null model ####
    if fit_null:
        print('Estimating variance components')
    if covar is not None and fit_null:
        # Match covariates #
        covariates.filter_ids(G.ids, verbose=False)
        null_model, sigma2, tau, null_alpha, null_alpha_cov = fit_sibreg_model(y, covariates.gts, G.fams, add_intercept=True,
                                                   tau_init=tau_init)

    elif fit_null:
        null_model, sigma2, tau = fit_sibreg_model(y, np.ones((y.shape[0], 1)), G.fams,
                                                                               tau_init = tau_init, return_fixed = False)
    else:
        null_model = model(y, np.ones((y.shape[0], 1)), G.fams)
    if fit_null:
        print('Family variance estimate: '+str(round(sigma2/tau,4)))
        print('Residual variance estimate: ' + str(round(sigma2,4)))
    ##### Transform genotypes and phenotypes ######
    if verbose:
        print('Transforming genotypes and phenotypes')
    if tau is None or sigma2 is None:
        raise(ValueError('Must provide variance components if not fitting null model'))
    L = null_model.sigma_inv_root(tau, sigma2)
    G.diagonalise(L)
    if covar is None:
        y = transform_phenotype(L, y, G.fam_indices)
    else:
        null_mean = null_alpha[0]+covariates.gts.dot(null_alpha[1:null_alpha.shape[0]])
        y = transform_phenotype(L, y, G.fam_indices, null_mean)
    ### Fit models for SNPs ###
    if verbose:
        print('Estimating SNP effects')
    alpha, alpha_cov, alpha_ses = fit_models(y,G)
    return chrom, pos, alleles, freqs, G.sid, alpha, alpha_cov, alpha_ses, sigma2, tau, null_alpha

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
    parser.add_argument('--start',type=int,help='Start index of the SNPs to use in the observed genotype file, counting from zero',default = 0)
    parser.add_argument('--end',type=int,help='End index of SNPs in the observed genotype file. The script will use SNPs with indices in the range [start,end-1], indexing from zero.', default=None)
    parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)', default=0.01)
    parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)', default=5)
    parser.add_argument('--min_info',type=float,help='Ignore SNPs with imputation INFO score below this threshold (default 0.90)', default=0.90)
    parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)', default='NA')
    parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                        default=1)
    parser.add_argument('--output_covar_ests',action='store_true',help='Output null model estimates of covariate fixed effects (default False)', default=False)
    parser.add_argument('--batch_size',type=int,help='Batch size of SNPs in thousands to load at a time (reduce to reduce memory requirements)',default=100)
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
        snp_ids = np.loadtxt(args.bed+'.bim',dtype=str,usecols=1)
    elif args.bgen is not None:
        gts_f = args.bgen+'.bgen'
        bgen = open_bgen(gts_f, verbose=False)
        snp_ids = bgen.ids
    ####### Compute batches #######
    print('Found '+str(snp_ids.shape[0])+' variants in '+gts_f)
    if args.end is not None:
        snp_ids = snp_ids[args.start:args.end]
        print('Using SNPs with indices from '+str(args.start)+' to '+str(args.end))
    # Remove duplicates
    unique_snps, counts = np.unique(snp_ids, return_counts=True)
    non_duplicate = set(unique_snps[counts==1])
    if np.sum(counts>1)>0:
        print('Removing '+str(np.sum(counts>1))+' duplicate SNP ids')
    snp_ids = snp_ids[np.array([x in non_duplicate for x in snp_ids])]
    snp_dict = make_id_dict(snp_ids)
    # Compute batches
    batch_bounds = compute_batch_boundaries(snp_ids,args.batch_size*1000)
    print('Using '+str(batch_bounds.shape[0])+' batches')
    alpha_dim = 3
    if args.fit_sib:
        alpha_dim += 1
    elif args.parsum:
        alpha_dim = alpha_dim-1
    # Create output files
    alpha = np.zeros((snp_ids.shape[0],alpha_dim),dtype=np.float32)
    alpha[:] = np.nan
    alpha_cov = np.zeros((snp_ids.shape[0], alpha_dim, alpha_dim),dtype=np.float32)
    alpha_cov[:] = np.nan
    alpha_ses = np.zeros((snp_ids.shape[0],alpha_dim),dtype=np.float32)
    alpha_ses[:] = np.nan
    chrom = np.zeros((snp_ids.shape[0]),dtype=str)
    alleles = np.zeros((snp_ids.shape[0],2),dtype=str)
    pos = np.zeros((snp_ids.shape[0]),dtype=int)
    freqs = np.zeros((snp_ids.shape[0]),dtype=np.float32)
    ##############  Process batches of SNPs ##############
    ######## Fit null model in first batch ############
    batch_chrom, batch_pos, batch_alleles, batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses, sigma2, tau, null_alpha = process_batch(
                   snp_ids[batch_bounds[0, 0]:batch_bounds[0, 1]],
                   y, pheno_ids, args.pargts, gts_f,
                   fit_null=True, covar=args.covar, parsum=args.parsum, fit_sib=args.fit_sib,
                   max_missing=args.max_missing, min_maf=args.min_maf, min_info=args.min_info, tau_init=args.tau_init)
    # Fill in SNP info
    chrom[batch_bounds[0,0]:batch_bounds[0,1]] = np.array(batch_chrom,dtype=str)
    pos[batch_bounds[0,0]:batch_bounds[0,1]] = batch_pos
    alleles[batch_bounds[0,0]:batch_bounds[0,1],:] = batch_alleles
    freqs[batch_bounds[0,0]:batch_bounds[0,1]] = batch_freqs
    # Fill in fitted SNPs
    batch_indices = np.array([snp_dict[x] for x in batch_snps])
    alpha[batch_indices,:] = batch_alpha
    alpha_cov[batch_indices,:,:] = batch_alpha_cov
    alpha_ses[batch_indices,:] = batch_alpha_ses
    print('Done batch 1 out of '+str(batch_bounds.shape[0]))
    ########### Process remaining batches ###########
    for i in range(1,batch_bounds.shape[0]):
        batch_chrom, batch_pos, batch_alleles, batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses, sigma2, tau, null_alpha = process_batch(
            snp_ids[batch_bounds[i, 0]:batch_bounds[i, 1]], y, pheno_ids, args.pargts, gts_f, fit_null=False,
            tau=tau, sigma2=sigma2, null_alpha=null_alpha, covar=args.covar, parsum=args.parsum, fit_sib=args.fit_sib,
            max_missing=args.max_missing, min_maf=args.min_maf, min_info=args.min_info, tau_init=args.tau_init)
        # Fill in SNP info
        chrom[batch_bounds[i, 0]:batch_bounds[i, 1]] = batch_chrom
        pos[batch_bounds[i, 0]:batch_bounds[i, 1]] = batch_pos
        alleles[batch_bounds[i, 0]:batch_bounds[i, 1], :] = batch_alleles
        freqs[batch_bounds[i, 0]:batch_bounds[i, 1]] = batch_freqs
        # Fill in fitted SNPs
        batch_indices = np.array([snp_dict[x] for x in batch_snps])
        alpha[batch_indices, :] = batch_alpha
        alpha_cov[batch_indices, :, :] = batch_alpha_cov
        alpha_ses[batch_indices, :] = batch_alpha_ses
        print('Done batch '+str(i+1)+'out of '+str(batch_bounds.shape[0]))
    ######## Save output #########
    write_output(chrom, snp_ids, pos, alleles, args.outprefix, args.parsum, args.fit_sib, alpha, alpha_ses, alpha_cov,
                 sigma2, tau, freqs)
