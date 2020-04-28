#!/well/kong/users/wiw765/anaconda3/bin/python
import numpy as np
import numpy.ma as ma
from pysnptools.snpreader import Bed, Pheno
from sibreg import sibreg
import h5py, argparse, code
import pandas as pd

class gtarray(object):
    def __init__(self,garray,ids, sid = None, id_index = 0, alleles = None, fams= None):
        if type(garray) == np.ndarray or type(garray)==np.ma.core.MaskedArray:
            if type(garray) == np.ndarray:
                self.gts = ma.array(garray,mask=np.isnan(garray))
            else:
                self.gts = garray
            self.shape = garray.shape
            self.ndim = garray.ndim
            self.dtype = garray.dtype
        else:
            raise ValueError('Genotypes must be a numpy ndarray')
        if garray.shape[0] == ids.shape[0]:
            self.ids = ids
            self.id_dict = make_id_dict(ids, col = id_index)
        else:
            raise ValueError('Shape of genotypes and ids does not match')
        if sid is not None:
            if sid.shape[0] == garray.shape[1]:
                self.snp_index = 1
                self.sid = sid
                self.sid_dict = make_id_dict(sid)
            elif sid.shape[0] == garray.shape[2]:
                self.snp_index = 2
                self.sid = sid
                self.sid_dict = make_id_dict(sid)
            else:
                raise ValueError('Shape of SNP ids (sid) does not match shape of genotype array')
        if alleles is not None:
            if alleles.shape[0] == sid.shape[0]:
                self.alleles = alleles
            else:
                raise ValueError('Size of alleles does not match size of SNP ids')
        else:
            self.alleles = None

        if fams is not None:
            if fams.shape[0] == ids.shape[0] and fams.ndim==1:
                self.fams = fams
            else:
                raise ValueError('Fams not of same length as IDs')
        else:
            self.fams = None

        self.mean_normalised = False

    def mean_normalise(self):
        if not self.mean_normalised:
            if self.gts.ndim == 2:
                self.gts = self.gts - ma.mean(self.gts,axis=0)
            elif self.gts.ndim==3:
                for i in range(0, self.gts.shape[1]):
                    self.gts[:, i, :] = self.gts[:, i, :] - ma.mean(self.gts[:, i, :], axis=0)
            self.mean_normalised = True

    def scale(self):
        if self.gts.ndim == 2:
            self.gts = self.gts/ma.std(self.gts, axis=0)
        elif self.gts.ndim == 3:
            for i in range(0, self.gts.shape[1]):
                self.gts[:, i, :] = self.gts[:, i, :]/ma.std(self.gts[:, i, :], axis=0)

    def fill_NAs(self):
        if not self.mean_normalised:
            self.mean_normalise()
        self.gts[self.gts.mask] = 0

    def add(self,garray):
        if type(garray)==gtarray:
            pass
        else:
            raise ValueError('Must add to another gtarray')

        if not self.gts.ndim == garray.gts.ndim:
            raise ValueError('Arrays must have same number of dimensions')

        if self.gts.ndim == 2:
            if not self.gts.shape[1] == garray.gts.shape[1]:
                raise ValueError('Arrays must have same dimensions (apart from first)')

        if self.gts.ndim == 3:
            if not self.gts.shape[1:3] == garray.gts.shape[1:3]:
                raise ValueError('Arrays must have same dimensions (apart from first)')

        # Match IDs
        common_ids = list(self.id_dict.keys() & garray.id_dict.keys())
        if len(common_ids)==0:
            raise ValueError('No IDs in common')
        self_index = np.array([self.id_dict[x] for x in common_ids])
        other_index = np.array([garray.id_dict[x] for x in common_ids])

        # Out
        if self.ids.ndim == 1:
            ids_out = self.ids[self_index]
        else:
            ids_out = self.ids[self_index,:]

        if self.gts.ndim ==2:
            add_gts = self.gts[self_index,:]+garray.gts[other_index,:]
        else:
            add_gts = self.gts[self_index, :,:] + garray.gts[other_index, :,:]

        return gtarray(add_gts,ids_out,self.sid,alleles = self.alleles, fams = self.fams[self_index])


class pgs(object):
    def __init__(self,snp_ids,weights,alleles):
        if snp_ids.shape[0] == weights.shape[0] and alleles.shape[0] == weights.shape[0] and alleles.shape[1]==2:
            self.snp_ids = snp_ids
            self.snp_dict = make_id_dict(snp_ids)
            self.weights = weights
            self.alleles = alleles
        else:
            raise ValueError('All inputs must have the same dimension')

    def compute(self,garray, cols = None):
        if type(garray) == gtarray:
            garray.fill_NAs()
        else:
            raise ValueError('Must be of gtarray class')
        if garray.alleles is None:
            raise ValueError('Alleles of genotype matrix must be provided')
        # Match SNP IDs
        in_pgs_snps = np.array([x in self.snp_dict for x in garray.sid])
        nmatch = np.sum(in_pgs_snps)
        if nmatch==0:
            raise ValueError('No overlap between PGS SNPs and genotype SNPs')
        # Get weights
        matched_snps = garray.sid[in_pgs_snps]
        matched_alleles = garray.alleles[in_pgs_snps,:]
        snp_indices = np.zeros((nmatch),dtype=int)
        for i in range(0,nmatch):
            snp_indices[i] = self.snp_dict[matched_snps[i]]
        weights_compute = self.weights[snp_indices]
        alleles = self.alleles[snp_indices,:]

        # Match alleles and adjust weights
        a_match = np.logical_and(alleles[:,0] == matched_alleles[:, 0], alleles[:,1] == matched_alleles[:, 1])
        a_reverse = np.logical_and(alleles[:,0] == matched_alleles[:, 1], alleles[:,1] == matched_alleles[:, 0])
        a_nomatch = np.logical_and(np.logical_not(a_match), np.logical_not(a_reverse))
        n_nomatch = np.sum(a_nomatch)
        if n_nomatch > 0:
            print('Removing ' + str(n_nomatch) + ' SNPs due to allele mismatch between genotypes and PGS alleles')
            weights_compute[a_nomatch] = 0
        weights_compute[a_reverse] = -weights_compute[a_reverse]

        ### Compute PGS
        if garray.ndim == 2:
            pgs_val = garray.gts[:,in_pgs_snps].dot(weights_compute)
        elif garray.ndim ==3:
            pgs_val = np.zeros((garray.gts.shape[0],garray.gts.shape[1]),garray.dtype)
            for i in range(0,garray.gts.shape[1]):
                pgs_val[:,i] = garray.gts[:,i,in_pgs_snps].dot(weights_compute)

        return gtarray(pgs_val, garray.ids, sid = cols, fams = garray.fams)

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

def encode_str_array(x):
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.encode('ascii') for y in x])
    return x_out.reshape(x_shape)

def find_individuals_with_sibs(ids,ped,gts_ids, return_ids_only = False):
    # Find genotyped sibships of size > 1
    ped_dict = make_id_dict(ped, 1)
    ids_in_ped = np.array([x in ped_dict for x in gts_ids])
    gts_fams = np.array([ped[ped_dict[x], 0] for x in gts_ids[ids_in_ped]])
    fams, counts = np.unique(gts_fams, return_counts=True)
    sibships = set(fams[counts > 1])
    # Find individuals with genotyped siblings
    ids_in_ped = np.array([x in ped_dict for x in ids])
    ids = ids[ids_in_ped]
    ids_fams = np.array([ped[ped_dict[x], 0] for x in ids])
    ids_with_sibs = np.array([x in sibships for x in ids_fams])
    ids = ids[ids_with_sibs]
    ids_fams = ids_fams[ids_with_sibs]
    if return_ids_only:
        return ids
    else:
        return ids, ids_fams, gts_fams

def get_fam_means(ids,ped,gts,gts_ids,remove_proband = True):
    ids, ids_fams, gts_fams = find_individuals_with_sibs(ids,ped,gts_ids)
    fams = np.unique(ids_fams)
    fams_dict = make_id_dict(fams)
    # Compute sums of genotypes in each family
    fam_sums = np.zeros((fams.shape[0],gts.shape[1]),dtype=gts.dtype)
    fam_counts = np.zeros((fams.shape[0]),dtype=int)
    for i in range(0,fams.shape[0]):
        fam_indices = np.where(gts_fams==fams[i])[0]
        fam_sums[i,:] = np.sum(gts[fam_indices,:],axis=0)
        fam_counts[i] = fam_indices.shape[0]
    # Place in vector corresponding to IDs
    gts_id_dict = make_id_dict(gts_ids)
    G_sib = np.zeros((ids.shape[0],gts.shape[1]),dtype = np.float32)
    for i in range(0,ids.shape[0]):
        G_sib[i,:] = fam_sums[fams_dict[ids_fams[i]],:]
        n_i = fam_counts[fams_dict[ids_fams[i]]]
        if remove_proband:
            G_sib[i,:] = G_sib[i,:] - gts[gts_id_dict[ids[i]],:]
            n_i = n_i-1
        G_sib[i,:] = G_sib[i,:]/float(n_i)
    return gtarray(G_sib,ids)


def find_par_gts(pheno_ids,ped,fams,gts_id_dict):
    # Whether mother and father have observed/imputed genotypes
    par_status = np.zeros((pheno_ids.shape[0],2),dtype=int)
    par_status[:] = -1
    # Indices of obsered/imputed genotypes in relevant arrays
    gt_indices = np.zeros((pheno_ids.shape[0],3),dtype=int)
    gt_indices[:] = -1
    ## Build dictionaries
    ## Where each individual is in the pedigree
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
    #
    return par_status, gt_indices, fam_labels

def make_gts_matrix(gts,imp_gts,par_status,gt_indices):
    if np.min(gt_indices)<0:
        raise ValueError('Missing genotype index')
    N = gt_indices.shape[0]
    G = np.zeros((N,3,gts.shape[1]),np.float32)
    # Observed proband genotype
    G[:, 0, :] = gts[gt_indices[:, 0], :]
    # Observed fathers
    G[par_status[:,0]==0, 1, :] = gts[gt_indices[par_status[:,0]==0, 1], :]
    # Imputed fathers
    G[par_status[:,0]==1, 1, :] = imp_gts[gt_indices[par_status[:,0]==1, 1], :]
    # Observed mothers
    G[par_status[:,1]==0, 2, :] = gts[gt_indices[par_status[:,1]==0, 2], :]
    # Imputed mothers
    G[par_status[:,1]==1, 2, :] = imp_gts[gt_indices[par_status[:,1]==1, 2], :]
    return G


def get_gts_matrix(par_gts_f,gts_f,snp_ids,ids = None, sib = False):
    ####### Find parental status #######
    ### Imputed parental file ###
    par_gts_f = h5py.File(par_gts_f,'r')
    # Get families
    fams = convert_str_array(np.array(par_gts_f['families']))
    # Get pedigree
    ped = convert_str_array(np.array(par_gts_f['pedigree']))
    ped = ped[1:ped.shape[0],:]
    # Remove control families
    controls = np.array([x[0]=='_' for x in ped[:,0]])
    ped = ped[np.logical_not(controls),:]
    ### Genotype file ###
    bim = gts_f.split('.bed')[0] + '.bim'
    gts_f = Bed(gts_f)
    alleles = np.loadtxt(bim, dtype='U')[:,4:6]
    # get ids of genotypes and make dict
    gts_ids = gts_f.iid[:,1]
    gts_id_dict = make_id_dict(gts_ids)
    if ids is None:
        ids = gts_ids

    # Find mean of siblings
    if sib:
        ids = find_individuals_with_sibs(ids,ped,gts_ids, return_ids_only = True)
        print('Found '+str(ids.shape[0])+' individuals with genotyped siblings')

    ### Find parental status
    print('Checking for observed/imputed parental genotypes')
    par_status, gt_indices, fam_labels = find_par_gts(ids,ped,fams,gts_id_dict)
    # Find which individuals can be used
    none_missing = np.min(par_status, axis=1)
    none_missing = none_missing >= 0
    N = np.sum(none_missing)
    if N==0:
        raise ValueError(
            'No individuals with phenotype observations and complete observed/imputed genotype observations')
    print(str(N)+' individuals with phenotype observations and complete observed/imputed genotypes observations')
    # Take those that can be used
    gt_indices = gt_indices[none_missing,:]
    par_status = par_status[none_missing,:]
    ids = ids[none_missing]
    ## Read genotypes
    observed_indices = np.sort(np.unique(np.hstack((gt_indices[:,0],
                                  gt_indices[par_status[:,0]==0,1],
                                  gt_indices[par_status[:,1]==0,2]))))
    # Get imputed indices
    imp_indices = np.sort(np.unique(np.hstack((gt_indices[par_status[:,0]==1,1],
                                  gt_indices[par_status[:,1]==1,2]))))
    print('Matching observed and imputed SNPs')
    # Match SNPs from imputed and observed and restrict to those in list
    snp_set = set(snp_ids)
    imp_sid = convert_str_array(np.array(par_gts_f['sid']))
    obs_sid = gts_f.sid
    obs_sid_dict = make_id_dict(obs_sid)
    in_obs_sid = np.zeros((imp_sid.shape[0]),dtype=bool)
    obs_sid_index = np.zeros((imp_sid.shape[0]),dtype=int)
    for i in range(0,imp_sid.shape[0]):
        if imp_sid[i] in obs_sid_dict and imp_sid[i] in snp_set:
            in_obs_sid[i] = True
            obs_sid_index[i] = obs_sid_dict[imp_sid[i]]
    obs_sid_index = obs_sid_index[in_obs_sid]
    sid = imp_sid[in_obs_sid]
    alleles = alleles[obs_sid_index,:]
    if np.sum(in_obs_sid) == 0:
        raise ValueError('No SNPs in common between imputed, observed, and PGS')

    # Read imputed parental genotypes
    print('Reading imputed parental genotypes')
    imp_gts = np.array(par_gts_f['imputed_par_gts'][:,in_obs_sid])
    imp_gts = imp_gts[imp_indices,:]
    fams = fams[imp_indices]
    # Read observed genotypes
    print('Reading observed genotypes')
    gts = gts_f[observed_indices, obs_sid_index].read().val
    gts_ids = gts_f.iid[observed_indices,1]
    gts_id_dict = make_id_dict(gts_ids)
    # Find indices in reduced data
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, fams, gts_id_dict)
    print('Constructing family based genotype matrix')
    ### Make genotype design matrix
    if sib:
        G = np.zeros((ids.shape[0],4,gts.shape[1]),dtype = np.float32)
        G[:,np.array([0,2,3]),:] = make_gts_matrix(gts,imp_gts,par_status,gt_indices)
        G[:,1,:] = get_fam_means(ids, ped, gts, gts_ids, remove_proband=True).gts
    else:
        G = make_gts_matrix(gts, imp_gts, par_status, gt_indices)
    del gts
    del imp_gts
    return gtarray(G,ids,sid, alleles = alleles, fams = fam_labels)

def compute_pgs(par_gts_f,gts_f,pgs, sib = False):
    G = get_gts_matrix(par_gts_f,gts_f,pgs.snp_ids, sib = sib)
    if sib:
        cols = np.array(['proband','sibling','paternal','maternal'])
    else:
        cols = np.array(['proband','paternal','maternal'])
    return pgs.compute(G,cols)

def fit_null_model(y,X,fam_labels,tau_init):
    # Optimize null model
    sigma_2_init = np.var(y)*tau_init/(1+tau_init)
    null_model = sibreg.model(y, X, fam_labels)
    null_optim = null_model.optimize_model(np.array([sigma_2_init,args.tau_init]))
    null_alpha = null_model.alpha_mle(null_optim['tau'],null_optim['sigma2'],compute_cov = True)
    # Return values
    return null_optim['sigma2'], null_optim['tau'], null_alpha


def get_alpha_mle(y,pgs,fam_labels, add_intercept = False):
    # Initialise var pars
    sigma_2_init = np.var(y) / 2.0
    model = sibreg.model(y, pgs, fam_labels, add_intercept = add_intercept)
    optim = model.optimize_model(np.array([sigma_2_init,1]))
    print('Family variance estimate: '+str(round(optim['sigma2']/optim['tau'],4)))
    print('Residual variance estimate: ' + str(round(optim['sigma2'],4)))
    alpha = model.alpha_mle(optim['tau'],optim['sigma2'],compute_cov = True)
    return alpha

######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('outprefix',type=str,help='Location to output association statistic hdf5 file')
    parser.add_argument('--gts_list',type=str,help='File with list of bed files of observed genotypes', default = None)
    parser.add_argument('--pargts_list', type=str, help='File with list of imputed parental genotype HDF5 files', default = None)
    parser.add_argument('--weights',type=str,help='Location of the PGS allele weights', default = None)
    parser.add_argument('--phenofile',type=str,help='Location of the phenotype file',default = None)
    parser.add_argument('--pgs', type=str, help='Location of the pre-computed PGS', default=None)
    parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                        default=1)
    parser.add_argument('--fit_sib',action='store_true',default=False,help='Fit indirect effects from siblings')
    parser.add_argument('--sibdiff',action='store_true',default = False,help='Fit sibling difference in PGS model')
    parser.add_argument('--sibped',type=str,help='Path to pedigree file. By default uses pedigree in imputed parental genotype HDF5 file',default=None)
    parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                        default=1)
    parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)',default=0.01)
    parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)',default='NA')
    parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)',default=5)
    parser.add_argument('--append',action='store_true',default=False,help='Append results to existing output file with given outprefix (default overwrites existing')
    parser.add_argument('--no_covariate_estimates',action='store_true',default=False,help='Suppress output of covariate effect estimates')
    args=parser.parse_args()

    if args.weights is not None:
        if args.gts_list is None:
            raise ValueError('Weights provided but no observed genotypes provided')
        if args.pargts_list is None:
            raise ValueError('Weights provided but no imputed parental genotypes provided')
        print('Computing PGS from weights file')
        ####### Read PGS #######
        weights = pd.read_csv(args.weights, delimiter='\t')
        w_sid = np.array(weights.loc[:, 'sid'], dtype='U')
        print('Read weights for '+str(w_sid.shape[0])+' SNPs')
        beta = np.array(weights.loc[:, 'ldpred_beta'])
        a1 = np.array(weights.loc[:, 'nt1'], dtype='U')
        a2 = np.array(weights.loc[:, 'nt2'], dtype='U')
        p = pgs(w_sid,beta,np.vstack((a1,a2)).T)

        ###### Compute PGS ########
        G_list = []
        gts_list = np.loadtxt(args.gts_list,dtype='U')
        pargts_list = np.loadtxt(args.pargts_list,dtype='U')
        if not gts_list.shape[0] == pargts_list.shape[0]:
            raise ValueError('Lists of imputed and observed genotype files not of same length')
        print('Computing PGS')
        print('Using '+str(pargts_list[0])+' and '+str(gts_list[0]))
        pg = compute_pgs(pargts_list[0],gts_list[0],p, sib = args.fit_sib)
        for i in range(1,gts_list.shape[0]):
            print('Using ' + str(pargts_list[i]) + ' and ' + str(gts_list[i]))
            pg = pg.add(compute_pgs(pargts_list[i],gts_list[i],p, sib = args.fit_sib))
        # Normalise PGS
        pg.mean_normalise()
        pg.scale()
        print('PGS computed')

        ####### Write PGS to file ########
        print('Writing PGS to '+args.outprefix+'.pgs.hdf5')
        pgs_out = h5py.File(args.outprefix+'.pgs.hdf5','w')
        pgs_out['pgs'] = pg.gts
        pgs_out['ids'] = encode_str_array(pg.ids)
        pgs_out['cols'] = encode_str_array(pg.sid)
        pgs_out['fams'] = encode_str_array(pg.fams)
        pgs_out.close()
    elif args.pgs is not None:
        if args.phenofile is None:
            raise ValueError('Pre-computed PGS provided but no phenotype provided')
        print('Reading PGS from '+args.pgs)
        pgs_f = h5py.File(args.pgs, 'r')
        pg = gtarray(np.array(pgs_f['pgs']),
                     convert_str_array(np.array(pgs_f['ids'])),
                     sid = convert_str_array(np.array(pgs_f['cols'])),
                     fams = convert_str_array(np.array(pgs_f['fams'])))
        print('Normalising PGS')
        pg.mean_normalise()
        pg.scale()
        pgs_f.close()
    else:
        raise ValueError('Weights or PGS must be provided')

    if args.phenofile is not None:
        print('Fitting PGS for '+str(args.phenofile))
        pheno = Pheno(args.phenofile, missing=args.missing_char).read()
        # pheno = Pheno('phenotypes/eduyears_resid.ped', missing='NA').read()
        y = np.array(pheno.val)
        pheno_ids = np.array(pheno.iid)
        if y.ndim == 1:
            pass
        elif y.ndim == 2:
            y = y[:, args.phen_index - 1]
        else:
            raise ValueError('Incorrect dimensions of phenotype array')
        # Remove y NAs
        y_not_nan = np.logical_not(np.isnan(y))
        if np.sum(y_not_nan) < y.shape[0]:
            y = y[y_not_nan]
            pheno_ids = pheno_ids[y_not_nan, 1]
        print('Number of non-missing phenotype observations: ' + str(y.shape[0]))
        in_pgs = np.array([x in pg.id_dict for x in pheno_ids])
        y = y[in_pgs]
        pheno_ids = pheno_ids[in_pgs]
        gt_indices = np.array([pg.id_dict[x] for x in pheno_ids])
        print('Final sample size: '+str(gt_indices.shape[0]))
        alpha_imp = get_alpha_mle(y, pg.gts[gt_indices,:], pg.fams[gt_indices], add_intercept = True)
        # Estimate proband only model
        alpha_proband = get_alpha_mle(y, pg.gts[gt_indices, 0], pg.fams[gt_indices], add_intercept=True)
        # Get print out for fixed mean effects
        alpha_out = np.zeros((pg.sid.shape[0]+1, 2))
        alpha_out[0:pg.sid.shape[0], 0] = alpha_imp[0][1:(1+pg.sid.shape[0])]
        alpha_out[0:pg.sid.shape[0], 1] = np.sqrt(np.diag(alpha_imp[1])[1:(1+pg.sid.shape[0])])
        alpha_out[pg.sid.shape[0],0] = alpha_proband[0][1]
        alpha_out[pg.sid.shape[0],1] = np.sqrt(np.diag(alpha_proband[1])[1])
        print('Saving estimates to '+args.outprefix+ '.pgs_effects.txt')
        outcols = np.hstack((pg.sid,np.array(['associative']))).reshape((pg.sid.shape[0]+1,1))
        np.savetxt(args.outprefix + '.pgs_effects.txt',
                   np.hstack((outcols, np.array(alpha_out, dtype='S20'))),
                   delimiter='\t', fmt='%s')
        print('Saving sampling covariance matrix of estimates to ' + args.outprefix + '.pgs_vcov.txt')
        np.savetxt(args.outprefix + '.pgs_vcov.txt', alpha_imp[1][1:(1+pg.sid.shape[0]),1:(1+pg.sid.shape[0])])

        if args.sibdiff:
            ped = np.zeros((pg.fams.shape[0],2),dtype=pg.fams.dtype)
            ped[:,0] = pg.fams
            ped[:,1] = pg.ids
            G_sdiff = np.zeros((pg.gts.shape[0],2))
            G_sdiff[:,0] = pg.gts[:,0]
            G_sdiff[:,1] = get_fam_means(pg.ids, ped, pg.gts[:,0], pg.ids, remove_proband=False)
            alpha_sdiff = get_alpha_mle(y, G_sdiff[gt_indices,:], pg.fams[gt_indices], add_intercept=True)
            alpha_sdiff_out  = np.zeros((2,2))
            alpha_sdiff_out[:, 0] = alpha_sdiff[0][1:3]
            alpha_sdiff_out[:, 1] = np.sqrt(np.diag(alpha_sdiff[1])[1:3])
            np.savetxt(args.outprefix + '.pgs_effects.txt',
                       np.hstack((np.array(['direct','between-family']), np.array(alpha_sdiff_out, dtype='S20'))),
                       delimiter='\t', fmt='%s')