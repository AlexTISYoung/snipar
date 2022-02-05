import argparse, h5py, code
import numpy as np
from bgen_reader import open_bgen
from pysnptools.snpreader import Bed
from scipy.stats import chi2
from math import log10
import snipar.read as read
import snipar.lmm as lmm
from snipar.utilities import *
import snipar.preprocess as preprocess
from numba import njit, prange

def transform_phenotype(inv_root, y, fam_indices, null_mean = None):
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

# @njit(parallel=True)
# def fit_models(y,G):
#     alpha = np.zeros((G.shape[2],G.shape[1]),dtype=np.float_)
#     alpha_cov = np.zeros((G.shape[2],G.shape[1],G.shape[1]),dtype=np.float_)
#     for i in prange(G.shape[2]):
#         not_na = np.sum(np.isnan(G[:,:,i]),axis=1)==0
#         xtx = G[not_na,:,i].T


def write_output(chrom, snp_ids, pos, alleles, outprefix, parsum, sib, alpha, alpha_ses, alpha_cov, sigma2, tau, freqs):
    """
    Write fitted SNP effects and other parameters to output HDF5 file.
    """
    print('Writing output to ' + outprefix + '.sumstats.hdf5')
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
    outfile['sigma2'] = sigma2
    outfile['tau'] = tau
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

def write_txt_output(chrom, snp_ids, pos, alleles, outprefix, parsum, sib, alpha, alpha_cov, sigma2, tau, freqs):
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
    vy = (1+1/tau)*sigma2
    outstack = [outbim]
    for i in range(len(effects)):
        outstack.append(outarray_effect(alpha[:,i],alpha_ses_out[:,i],freqs,vy))
        header += [effects[i]+'_N',effects[i]+'_Beta',effects[i]+'_SE',effects[i]+'_Z',effects[i]+'_log10_P']
    outstack.append(np.round(alpha_corr_out,6))
    header += corrs
    # Output array
    outarray = np.row_stack((np.array(header),np.column_stack(outstack)))
    print('Writing text output to '+outprefix+'.sumstats.gz')
    np.savetxt(outprefix+'.sumstats.gz',outarray,fmt='%s')

def compute_batch_boundaries(snp_ids,batch_size):
    nsnp = snp_ids.shape[0]
    n_blocks = int(np.ceil(float(nsnp)/float(batch_size)))
    block_bounds = np.zeros((n_blocks,2),dtype=int)
    start = 0
    for i in range(n_blocks-1):
        block_bounds[i,0] = start
        block_bounds[i,1] = start+batch_size
        start += batch_size
    block_bounds[n_blocks-1,:] = np.array([start,nsnp])
    return block_bounds

def process_batch(y, pedigree, tau, sigma2, snp_ids=None, bedfile=None, bgenfile=None, par_gts_f=None, parsum=False,
                  fit_sib=False, max_missing=5, min_maf=0.01, verbose=False, print_sample_info=False):
    ####### Construct family based genotype matrix #######
    G = read.get_gts_matrix(ped=pedigree, bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, snp_ids=snp_ids, 
                                ids=y.ids, parsum=parsum, sib=fit_sib, print_sample_info=print_sample_info)
    G.compute_freqs()
    #### Filter SNPs ####
    if verbose:
        print('Filtering based on MAF')
    G.filter_maf(min_maf)
    if verbose:
        print('Filtering based on missingness')
    G.filter_missingness(max_missing)
    if verbose:
        print(str(G.shape[2])+' SNPs that pass filters')
    #### Fill NAs ####
    if verbose:
        print('Imputing missing values with population frequencies')
    NAs = G.fill_NAs()
    #### Match phenotype ####
    y.filter_ids(G.ids)
    ##### Transform genotypes and phenotypes ######
    if verbose:
        print('Transforming genotypes and phenotypes')
    null_model = lmm.model(y.gts[:,0], np.ones((y.shape[0], 1)), y.fams)
    L = null_model.sigma_inv_root(tau, sigma2)
    G.diagonalise(L)
    transformed_y = transform_phenotype(L, y.gts[:,0], G.fam_indices)
    ### Fit models for SNPs ###
    if verbose:
        print('Estimating SNP effects')
    alpha, alpha_cov, alpha_ses = fit_models(transformed_y,G)
    return G.freqs, G.sid, alpha, alpha_cov, alpha_ses

def process_chromosome(chrom_out, y, pedigree, tau, sigma2, outprefix, bedfile=None, bgenfile=None, par_gts_f=None,
                        fit_sib=False, parsum=False, max_missing=5, min_maf=0.01, batch_size=10000, 
                        no_hdf5_out=False, no_txt_out=False):
    ######## Check for bed/bgen #######
    if bedfile is None and bgenfile is None:
        raise(ValueError('Must supply either bed or bgen file with observed genotypes'))
    if bedfile is not None and bgenfile is not None:
        raise(ValueError('Both --bed and --bgen specified. Please specify one only'))
    if bedfile is not None:
        bed = Bed(bedfile,count_A1 = True)
        snp_ids = bed.sid
        pos = bed.pos[:,2]
        alleles = np.loadtxt(bedfile.split('.bed')[0]+'.bim',dtype=str,usecols=(4,5))
        chrom = bed.pos[:,0]
    elif bgenfile is not None:
        bgen = open_bgen(bgenfile, verbose=False)
        snp_ids = bgen.ids
        # If SNP IDs are broken, try rsids
        if np.unique(snp_ids).shape[0] == 1:
            snp_ids = bgen.rsids
        pos = np.array(bgen.positions)
        alleles = np.array([x.split(',') for x in bgen.allele_ids])
        chrom = np.array(bgen.chromosomes)
        # If chromosomse unknown, set to chromosome inferred from filename
        chrom[[len(x)==0 for x in chrom]] = chrom_out
    ####### Compute batches #######
    print('Found '+str(snp_ids.shape[0])+' SNPs')
    # Remove duplicates
    unique_snps, counts = np.unique(snp_ids, return_counts=True)
    non_duplicate = set(unique_snps[counts==1])
    if np.sum(counts>1)>0:
        print('Removing '+str(np.sum(counts>1))+' duplicate SNP ids')
        not_duplicated = np.array([x in non_duplicate for x in snp_ids])
        snp_ids = snp_ids[not_duplicated]
        pos = pos[not_duplicated]
        chrom = chrom[not_duplicated]
        alleles = alleles[not_duplicated,:]
    snp_dict = make_id_dict(snp_ids)
    # Compute batches
    batch_bounds = compute_batch_boundaries(snp_ids,batch_size)
    if batch_bounds.shape[0] == 1:
        print('Using 1 batch')
    else:
        print('Using '+str(batch_bounds.shape[0])+' batches')
    alpha_dim = 2
    if fit_sib:
        alpha_dim += 1
    if not parsum:
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
    for i in range(0,batch_bounds.shape[0]):
        if i==0:
            print_sample_info = True
        else:
            print_sample_info = False
        batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses = process_batch(y, pedigree, 
                    tau, sigma2, snp_ids=snp_ids[batch_bounds[i, 0]:batch_bounds[i, 1]], bedfile=bedfile, bgenfile=bgenfile, 
                    par_gts_f = par_gts_f, parsum=parsum, fit_sib=fit_sib,max_missing=max_missing, min_maf=min_maf,
                    print_sample_info=print_sample_info)
        # Fill in fitted SNPs
        batch_indices = np.array([snp_dict[x] for x in batch_snps])
        alpha[batch_indices, :] = batch_alpha
        alpha_cov[batch_indices, :, :] = batch_alpha_cov
        alpha_ses[batch_indices, :] = batch_alpha_ses
        freqs[batch_indices] = batch_freqs
        print('Done batch '+str(i+1)+' out of '+str(batch_bounds.shape[0]))
    ######## Save output #########
    if not chrom_out==0:
        outprefix = outprefix+'chr_'+str(chrom_out)
    if not no_hdf5_out:
        write_output(chrom, snp_ids, pos, alleles, outprefix, parsum, fit_sib, alpha, alpha_ses, alpha_cov,
                     sigma2, tau, freqs)
    if not no_txt_out:
        write_txt_output(chrom, snp_ids, pos, alleles, outprefix, parsum, fit_sib, alpha, alpha_cov,
                     sigma2, tau, freqs)

######### Command line arguments #########
parser=argparse.ArgumentParser()
parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
parser.add_argument('outprefix', type=str, help='Location to output summary statistics files. Outputs text output to outprefix.sumstats.gz and HDF5 output to outprefix.sumstats.hdf5')
parser.add_argument('--bedfiles',type=str,help='Address of observed genotype files in .bed format (without .bed suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22.', default = None)
parser.add_argument('--bgenfiles',type=str,help='Address of observed genotype files in .bgen format (without .bgen suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22.', default = None)
parser.add_argument('--impfiles', type=str, help='Address of hdf5 files with imputed parental genotypes (without .hdf5 suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22.', default = None)
parser.add_argument('--pedigree',type=str,help='Address of pedigree file. Must be provided if not provided imputed parental genotypes.',default=None)
parser.add_argument('--parsum',action='store_true',help='Regress onto proband and sum of parental genotypes (useful when parental genotypes imputed from sibs only)',default = False)
parser.add_argument('--fit_sib',action='store_true',help='Fit indirect effect from sibling ',default=False)
parser.add_argument('--covar',type=str,help='Path to file with covariates: plain text file with columns FID, IID, covar1, covar2, ..', default=None)
parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                    default=1)
parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)', default=0.01)
parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)', default=5)
parser.add_argument('--output_covar_ests',action='store_true',help='Output null model estimates of covariate fixed effects (default False)', default=False)
parser.add_argument('--batch_size',type=int,help='Batch size of SNPs to load at a time (reduce to reduce memory requirements)',default=100000)
parser.add_argument('--no_hdf5_out',action='store_true',help='Suppress HDF5 output of summary statistics',default=False)
parser.add_argument('--no_txt_out',action='store_true',help='Suppress text output of summary statistics',default=False)
parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)', default='NA')
parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                    default=1)
args=parser.parse_args()

# Check arguments
if args.bedfiles is None and args.bgenfiles is None:
    raise(ValueError('Must provide one of --bedfiles and --bgenfiles'))
if args.bedfiles is not None and args.bgenfiles is not None:
    raise(ValueError('Both bedfiles and bgenfiles provided. Please provide only one'))
if args.impfiles is None and args.pedigree is None:
    raise(ValueError('Must provide pedigree if not providing imputed parental genotypes file(s)'))

# Find observed and imputed files
if args.impfiles is None:
    if args.bedfiles is not None:
        bedfiles, chroms = preprocess.parse_obsfiles(args.bedfiles, 'bed')
        bgenfiles = [None for x in range(chroms.shape[0])]
    elif args.bgenfiles is not None:
        bgenfiles, chroms = preprocess.parse_obsfiles(args.bgenfiles, 'bgen')
        bedfiles = [None for x in range(chroms.shape[0])]
    pargts_list = [None for x in range(chroms.shape[0])]
else:
    if args.bedfiles is not None:
        bedfiles, pargts_list, chroms = preprocess.parse_filelist(args.bedfiles, args.impfiles, 'bed')
        bgenfiles = [None for x in range(chroms.shape[0])]
    elif args.bgenfiles is not None:
        bgenfiles, pargts_list, chroms = preprocess.parse_filelist(args.bgenfiles, args.impfiles, 'bgen')
        bedfiles = [None for x in range(chroms.shape[0])]
if chroms.shape[0]==0:
    raise(ValueError('No input genotype files found'))

# Read phenotype and covariates
######### Read Phenotype ########
y = read.phenotype.read_phenotype(args.phenofile, missing_char=args.missing_char, phen_index=args.phen_index)
######## Read covariates ########
if args.covar is not None:
    print('Reading covariates')
    covariates = read.phenotype.read_covariates(args.covar, pheno_ids=y.ids, missing_char=args.missing_char)
    # Match to pheno ids
    covariates.filter_ids(y.ids)
else:
    covariates = None

# Read pedigree
if args.impfiles is None:
    print('Reading pedigree from '+str(args.pedigree))
    ped = np.loadtxt(args.pedigree,dtype=str)
    if ped.shape[1] < 4:
        raise(ValueError('Not enough columns in pedigree file'))
    elif ped.shape[1] > 4:
        print('Warning: pedigree file has more than 4 columns. The first four columns only will be used')
    # Remove rows with missing parents
    sibpairs, ped = preprocess.get_sibpairs_from_ped(ped)
    if sibpairs is not None:
        print('Found '+str(sibpairs.shape[0])+' sibling pairs in pedigree')
    else:
        print('Found 0 sibling pairs')
else:
    # Read pedigree
    par_gts_f = h5py.File(pargts_list[0],'r')
    ped = convert_str_array(par_gts_f['pedigree'])
    ped = ped[1:ped.shape[0]]
    # Remove control fams
    controls = np.array([x[0]=='_' for x in ped[:,0]])
    ped = ped[~controls,:]

####### Fit null model ######
# Match to pedigree
ped_dict = make_id_dict(ped,1)
y.filter_ids(ped[:,1])
print(str(y.shape[0])+' individuals with phenotype values found in pedigree')
ped_indices = np.array([ped_dict[x] for x in y.ids])
y.fams = ped[ped_indices,0]
print(' Fitting variance components')
if args.covar is not None:
    # Match covariates
    covariates.filter_ids(y.ids)
    # Fit null model
    null_model, sigma2, tau, null_alpha, null_alpha_cov = lmm.fit_model(y.gts[:,0], covariates.gts, y.fams, add_intercept=True,
                                                                        tau_init=args.tau_init)
    # Adjust for covariates
    y.gts = y.gts-(null_alpha[0]+covariates.gts.dot(null_alpha[1:null_alpha.shape[0]]))
else:
    # Fit null model
    null_model, sigma2, tau = lmm.fit_model(y.gts[:,0], np.ones((y.shape[0], 1)), y.fams,
                                            tau_init = args.tau_init, return_fixed = False)
    y.gts = y.gts-np.mean(y.gts)
print('Family variance estimate: '+str(round(sigma2/tau,4)))
print('Residual variance estimate: ' + str(round(sigma2,4)))

for i in range(chroms.shape[0]):
    if args.bedfiles is not None:
        print('Reading observed genotypes from '+bedfiles[i])
    if args.bgenfiles is not None:
        print('Reading observed genotypes from '+bgenfiles[i])
    if args.impfiles is not None:
        print('Reading imputed genotypes from '+pargts_list[i])
    if chroms.shape[0]>1:
        print('Estimating SNP effects for chromosome '+str(chroms[i]))
    else:
        print('Estimaing SNP effects')
    process_chromosome(chroms[i], y, ped, tau, sigma2, args.outprefix, bedfile=bedfiles[i], bgenfile=bgenfiles[i], 
                        par_gts_f=pargts_list[i], fit_sib=args.fit_sib, parsum=args.parsum, 
                        max_missing=args.max_missing, min_maf=args.min_maf, batch_size=args.batch_size, 
                        no_hdf5_out=args.no_hdf5_out, no_txt_out=args.no_txt_out)