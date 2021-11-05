from sibreg.sibreg import *
from sibreg.model import *
import argparse, h5py
import numpy as np
from bgen_reader import open_bgen
from scipy.stats import chi2
from math import log10
from pysnptools.snpreader import Bed


FORMAT = '%(asctime)-15s :: %(levelname)s :: %(filename)s :: %(funcName)s :: %(message)s'
# numeric_level = getattr(logging, loglevel.upper(), None)
# if not isinstance(numeric_level, int):
#     raise ValueError('Invalid log level: %s' % loglevel)
logging.basicConfig(format=FORMAT, level=logging.DEBUG if __debug__ else logging.INFO)
logger = logging.getLogger(__name__)


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

def write_output(chrom, snp_ids, pos, alleles, outprefix, parsum, sib, alpha, alpha_ses, alpha_cov, sigma_grm, sigma_sib, sigma_res, freqs):
    """
    Write fitted SNP effects and other parameters to output HDF5 file.
    """
    logger.info('Writing output to ' + outprefix + '.sumstats.hdf5')
    outfile = h5py.File(outprefix+'.sumstats.hdf5', 'w')
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
    outfile['sigma_grm'] = sigma_grm
    outfile['sigma_sib'] = sigma_sib
    outfile['sigma_res'] = sigma_res
    outfile['freqs'] = freqs
    outfile.close()

def outarray_effect(est, ses, freqs, vy):
    N_effective = vy/(2*freqs*(1-freqs)*np.power(ses,2))
    Z = est/ses
    P = -log10(np.exp(1))*chi2.logsf(np.power(Z,2),1)
    array_out = np.column_stack((N_effective,est,ses,Z,P))
    array_out = np.round(array_out, decimals=6)
    array_out[:,0] = np.round(array_out[:,0], 0)
    return array_out

def write_txt_output(chrom, snp_ids, pos, alleles, outprefix, parsum, sib, alpha, alpha_cov, sigma_grm, sigma_sib, sigma_res, freqs):
    outbim = np.column_stack((chrom, snp_ids, pos, alleles,np.round(freqs,3)))
    header = ['chromosome','SNP','pos','A1','A2','freq']
    # Which effects to estimate
    effects = ['direct']
    if sib:
        effects.append('sib')
    if not parsum:
        effects += ['paternal','maternal']
    effects += ['avg_parental','population']
    effects = np.array(effects)
    if not parsum:
        paternal_index = np.where(effects=='paternal')[0][0]
        maternal_index = np.where(effects=='maternal')[0][0]
    avg_par_index = np.where(effects=='avg_parental')[0][0]
    population_index = avg_par_index+1
    # Get transform matrix
    A = np.zeros((len(effects),alpha.shape[1]))
    A[0:alpha.shape[1],0:alpha.shape[1]] = np.identity(alpha.shape[1])
    if not parsum:
        A[alpha.shape[1]:(alpha.shape[1]+2), :] = 0.5
        A[alpha.shape[1], 0] = 0
        A[alpha.shape[1]+1, 0] = 1
    else:
        A[alpha.shape[1], :] = 1
    # Transform effects
    alpha = alpha.dot(A.T)
    alpha_ses_out = np.zeros((alpha.shape[0],A.shape[0]))
    corrs = ['r_direct_avg_parental','r_direct_population']
    if sib:
        corrs.append('r_direct_sib')
    if not parsum:
        corrs.append('r_paternal_maternal')
    ncor = len(corrs)
    alpha_corr_out = np.zeros((alpha.shape[0],ncor))
    for i in range(alpha_cov.shape[0]):
        alpha_cov_i = A.dot(alpha_cov[i,:,:].dot(A.T))
        alpha_ses_out[i,:] = np.sqrt(np.diag(alpha_cov_i))
        # Direct to average parental
        alpha_corr_out[i,0] = alpha_cov_i[0,avg_par_index]/(alpha_ses_out[i,0]*alpha_ses_out[i,avg_par_index])
        # Direct to population
        alpha_corr_out[i,1] = alpha_cov_i[0,population_index]/(alpha_ses_out[i,0]*alpha_ses_out[i,population_index])
        # Direct to sib
        if sib:
            alpha_corr_out[i,2] = alpha_cov_i[0,1]/(alpha_ses_out[i,0]*alpha_ses_out[i,1])
        # Paternal to maternal
        if not parsum:
            alpha_corr_out[i,ncor-1] = alpha_cov_i[paternal_index,maternal_index]/(alpha_ses_out[i,maternal_index]*alpha_ses_out[i,paternal_index])
    # Create output array
    vy = sigma_grm + sigma_sib + sigma_res
    outstack = [outbim]
    for i in range(len(effects)):
        outstack.append(outarray_effect(alpha[:,i],alpha_ses_out[:,i],freqs,vy))
        header += [effects[i]+'_N',effects[i]+'_Beta',effects[i]+'_SE',effects[i]+'_Z',effects[i]+'_log10_P']
    outstack.append(np.round(alpha_corr_out,6))
    header += corrs
    # Output array
    outarray = np.row_stack((np.array(header),np.column_stack(outstack)))
    logger.info('Writing text output to '+outprefix+'.sumstats.gz')
    np.savetxt(outprefix+'.sumstats.gz',outarray,fmt='%s')

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

def process_batch(snp_ids, pheno_ids, pargts_f, gts_f, parsum=False,
                  fit_sib=False, max_missing=5, min_maf=0.01, min_info=0.9, verbose=False, print_sample_info=False):
    ####### Construct family based genotype matrix #######
    G = get_gts_matrix(pargts_f+'.hdf5', gts_f, snp_ids=snp_ids, ids=pheno_ids, parsum=parsum, sib=fit_sib, print_sample_info=print_sample_info)
    # Check for empty fam labels
    no_fam = np.array([len(x) == 0 for x in G.fams])
    if np.sum(no_fam) > 0:
        ValueError('No family label from pedigree for some individuals')
    G.compute_freqs()
    #### Filter SNPs ####
    if verbose:
        logger.info('Filtering based on MAF')
    G.filter_maf(min_maf)
    gt_filetype = gts_f.split('.')[1]
    if gt_filetype=='bed':
        if verbose:
            logger.info('Filtering based on missingness')
        G.filter_missingness(max_missing)
    if gt_filetype=='bgen':
        if verbose:
            logger.info('Filtering based on imputation INFO')
        G.filter_info(min_info)
    if verbose:
        logger.info(str(G.shape[2])+' SNPs that pass filters')
    #### Fill NAs ####
    if verbose:
        logger.info('Imputing missing values with population frequencies')
    NAs = G.fill_NAs()
    
    alpha, alpha_cov, alpha_ses = model.fit_snps_eff(G.gts)
    return G.chrom, G.pos, G.alleles, G.freqs, G.sid, alpha, alpha_cov, alpha_ses 

######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('pargts', type=str, help='HDF5 file with imputed parental genotypes (without .hdf5 suffix)')
    parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
    parser.add_argument('--bed',type=str,help='Bed file with observed genotypes (without .bed suffix).',default=None)
    parser.add_argument('--bgen',type=str,help='bgen file with observed genotypes (without .bgen suffix).',default=None)
    parser.add_argument('--outprefix', type=str, help='Location to output association statistic hdf5 file. Outputs text output to outprefix.sumstats.gz and HDF5 output to outprefix.sumstats.hdf5', default='')
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
    parser.add_argument('--batch_size',type=int,help='Batch size of SNPs to load at a time (reduce to reduce memory requirements)',default=100000)
    parser.add_argument('--no_hdf5_out',action='store_true',help='Suppress HDF5 output of summary statistics',default=False)
    parser.add_argument('--no_txt_out',action='store_true',help='Suppress text output of summary statistics',default=False)

    parser.add_argument('--hapmap_bed', type=str, help='Bed file with observed hapmap3 snps (without suffix, chromosome number should be #).', default=None)
    parser.add_argument('--gcta_path', type=str, help='Path to gcta executable.', default=None)
    parser.add_argument('--grm_path', type=str, help='Path to gcta grm output (without prefix).', default=None)
    parser.add_argument('--plink_path', type=str, help='Path to plink2 executable.', default=None)

    parser.add_argument('--ibdseg_path', type=str, help='Path to KING ibdseg output (without .seg prefix).', default=None)

    parser.add_argument('--sparse_thres', type=float, help='Threshold of GRM/IBD sparsity', default=0.05)

    args=parser.parse_args()

    if args.ibdseg_path is not None and (args.grm_path is not None or args.gcta_path is not None):
        raise parser.error('Only one of ibdseg_path and grm_path/gcta_path should be supplied.')
    if args.ibdseg_path is not None and (args.grm_path is None or args.gcta_path is None):
        raise parser.error('Need one of ibdseg_path and grm_path/gcta_path should be supplied.') 
    if args.gcta_path is not None and (args.plink_path is None or args.hapmap_bed is None):
        raise parser.error('Should provide plink_path and hapmap_bed is gcta_path is supplied.')

    ######### Read Phenotype ########
    y, pheno_ids = read_phenotype(args.phenofile, missing_char=args.missing_char, phen_index=args.phen_index)
    ######## Read covariates ########
    if args.covar is not None:
        logger.info('Reading covariates')
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
        bed = Bed(gts_f,count_A1 = True)
        snp_ids = bed.sid
        pos = bed.pos[:,2]
        chrom = bed.pos[:,0]
        alleles = np.loadtxt(args.bed+'.bim',dtype=str,usecols=(4,5))
    elif args.bgen is not None:
        gts_f = args.bgen+'.bgen'
        bgen = open_bgen(gts_f, verbose=False)
        snp_ids = bgen.ids
        pos = np.array(bgen.positions)
        chrom = np.array(bgen.chromosomes)
        # If chromosomse unknown, set to zero
        chrom[[len(x)==0 for x in chrom]] = 0
        alleles = np.array([x.split(',') for x in bgen.allele_ids])
    ########## Construct GRM/IBD and sibship matrix ##########
    ids, fam_labels = get_ids_with_par(args.pargts, gts_f, pheno_ids)
    id_dict = make_id_dict(ids)
    logger.info('Building GRM...')
    if args.ibdseg_path is not None:
        grm = build_ibdseg_arr(args.ibdseg_path, id_dict=id_dict, keep=pheno_ids, thres=args.thres)
    else:
        if args.grm_path is None:
            run_gcta_grm(args.plink_path, args.gcta_path, args.hapmap_bed, args.outprefix, ids)
            grm_path = args.outprefix
        else:
            grm_path = args.grm_path
        grm = build_grm_arr(grm_path, id_dict=id_dict, thres=args.thres)
    logger.info('Building sibship matrix...')
    sib = build_sib_arr(fam_labels)
    ################## Optimize variance components ##################
    y = match_phenotype_(ids, y, pheno_ids)
    lmm = LinearMixedModel(y, (grm, sib), covar_X=covariates.gts)
    logger.info('Optimizing variance components...')
    lmm.scipy_optimize()
    sigma_grm, sigma_sib, sigma_res = lmm.varcomps
    ####### Compute batches #######
    logger.info('Found '+str(snp_ids.shape[0])+' variants in '+gts_f)
    if args.end is not None:
        snp_ids = snp_ids[args.start:args.end]
        pos = pos[args.start:args.end]
        chrom = chrom[args.start:args.end]
        alleles = alleles[args.start:args.end]
        logger.info('Using SNPs with indices from '+str(args.start)+' to '+str(args.end))
    # Remove duplicates
    unique_snps, counts = np.unique(snp_ids, return_counts=True)
    non_duplicate = set(unique_snps[counts==1])
    if np.sum(counts>1)>0:
        logger.info('Removing '+str(np.sum(counts>1))+' duplicate SNP ids')
        not_duplicated = np.array([x in non_duplicate for x in snp_ids])
        snp_ids = snp_ids[not_duplicated]
        pos = pos[not_duplicated]
        chrom = chrom[not_duplicated]
        alleles = alleles[not_duplicated,:]
    snp_dict = make_id_dict(snp_ids)
    # Compute batches
    batch_bounds = compute_batch_boundaries(snp_ids,args.batch_size)
    if batch_bounds.shape[0] == 1:
        logger.info('Using 1 batch')
    else:
        logger.info('Using '+str(batch_bounds.shape[0])+' batches')
    alpha_dim = 2
    if args.fit_sib:
        alpha_dim += 1
    if not args.parsum:
        alpha_dim += 1
    # Create output files
    alpha = np.zeros((snp_ids.shape[0],alpha_dim),dtype=np.float32)
    alpha[:] = np.nan
    alpha_cov = np.zeros((snp_ids.shape[0], alpha_dim, alpha_dim),dtype=np.float32)
    alpha_cov[:] = np.nan
    alpha_ses = np.zeros((snp_ids.shape[0],alpha_dim),dtype=np.float32)
    alpha_ses[:] = np.nan
    freqs = np.zeros((snp_ids.shape[0]),dtype=np.float32)
    freqs[:] = np.nan
    ##############  Process batches of SNPs ##############
    for i in range(0, batch_bounds.shape[0]):
        batch_chrom, batch_pos, batch_alleles, batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses = process_batch(
            snp_ids[batch_bounds[i, 0]:batch_bounds[i, 1]], pheno_ids, args.pargts, gts_f,
            parsum=args.parsum, fit_sib=args.fit_sib,
            max_missing=args.max_missing, min_maf=args.min_maf, min_info=args.min_info,)
        # Fill in fitted SNPs
        batch_indices = np.array([snp_dict[x] for x in batch_snps])
        alpha[batch_indices, :] = batch_alpha
        alpha_cov[batch_indices, :, :] = batch_alpha_cov
        alpha_ses[batch_indices, :] = batch_alpha_ses
        freqs[batch_indices] = batch_freqs
        logger.info('Done batch '+str(i+1)+' out of '+str(batch_bounds.shape[0]))
    ######## Save output #########
    if not args.no_hdf5_out:
        write_output(chrom, snp_ids, pos, alleles, args.outprefix, args.parsum, args.fit_sib, alpha, alpha_ses, alpha_cov,
                     sigma_grm, sigma_sib, sigma_res, freqs)
    if not args.no_txt_out:
        write_txt_output(chrom, snp_ids, pos, alleles, args.outprefix, args.parsum, args.fit_sib, alpha, alpha_cov,
                     sigma_grm, sigma_sib, sigma_res, freqs)
