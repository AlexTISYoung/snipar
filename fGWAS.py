from pysnptools.snpreader import Bed, Pheno
from sibreg.sibreg import *
import h5py, argparse, os, time, code
import numpy as np
import numpy.ma as ma

def read_phenotype(phenofile, missing_char = 'NA'):
    pheno = Pheno(args.phenofile, missing=args.missing_char).read()
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
    return y, pheno_ids

######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('gts',type=str,help='Path to bed file with sibling genotypes')
    parser.add_argument('pargts', type=str, help='Path to HDF5 file with imputed parental genotypes')
    parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
    parser.add_argument('outprefix',type=str,help='Location to output association statistic hdf5 file')
    parser.add_argument('--parsum',action='store_true',help='Regress onto proband and sum of parental genotypes (useful when parental genotypes imputed from sibs only)',default = False)
    parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                        default=1)
    parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                        default=1)
    parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)',default=0.01)
    parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)',default='NA')
    parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)',default=5)
    args=parser.parse_args()

    ######### Read Phenotype ########
    y , pheno_ids = read_phenotype(args.phenofile, missing_char=args.missing_char)

    ####### Construct family based genotype matrix #######
    G = get_gts_matrix(args.pargts, args.gts)

    code.interact(local=locals())
    # Filter genotypes based on frequency and missingness
    print('Filtering on frequency and MAF')
    freqs = ma.mean(gts,axis=0)/2.0
    missingness = ma.mean(gts.mask,axis=0)
    freqs_pass = np.logical_and(freqs > args.min_maf,freqs < (1-args.min_maf))
    print(str(freqs.shape[0]-np.sum(freqs_pass))+' SNPs with MAF<'+str(args.min_maf))
    missingness_pass = 100*missingness < args.max_missing
    print(str(freqs.shape[0] - np.sum(missingness_pass)) + ' SNPs with missingness >' + str(args.max_missing)+'%')
    filter_pass = np.logical_and(freqs_pass,missingness_pass)
    gts = gts[:,filter_pass]
    imp_gts = imp_gts[:,filter_pass]
    sid = sid[filter_pass,:]
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
    G = make_gts_matrix(gts,imp_gts,par_status,gt_indices, parsum = args.parsum)
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
        for i in range(G.shape[2]):
            G[:,label_indices,i] = np.dot(G[:,label_indices,i],L[label].T)
    ### Fit models for SNPs ###
    print('Estimating SNP effects')
    XTX = np.einsum('...ij,...ik', G, G)
    XTY = np.einsum('...ij,i',G,y)
    alpha = np.linalg.solve(XTX,XTY)
    alpha_cov = np.linalg.inv(XTX)
    alpha_ses = np.sqrt(np.diagonal(alpha_cov,axis1=1,axis2=2))
    t2 = time.time()
    print('Time: '+str(t2-t1)+' seconds')
    ### Output file ###
    print('Writing output to '+args.outprefix+'.hdf5')
    outfile = h5py.File(args.outprefix+'.hdf5','w')
    outfile['bim'] = encode_str_array(sid)
    if args.parsum:
        X_length = 2
        outcols = np.array(['direct','avg_parental'])
    else:
        X_length = 3
        outcols = np.array(['direct','paternal','maternal'])
    outfile.create_dataset('estimate_covariance',(sid.shape[0],X_length,X_length),dtype = 'f',chunks = True, compression = 'gzip', compression_opts=9)
    outfile.create_dataset('estimate', (sid.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile.create_dataset('estimate_ses', (sid.shape[0], X_length), dtype='f', chunks=True, compression='gzip',
                           compression_opts=9)
    outfile['estimate'][:] = alpha
    outfile['estimate_cols'] = encode_str_array(outcols)
    outfile['estimate_ses'][:] = alpha_ses
    outfile['estimate_covariance'][:] = alpha_cov
    outfile['sigma2'] = sigma2
    outfile['tau'] = tau
    outfile['N_L'] = N_L
    outfile['freqs'] = freqs
    outfile.close()
