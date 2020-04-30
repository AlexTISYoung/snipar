import numpy as np
import numpy.ma as ma
from pysnptools.snpreader import Bed, Pheno
from scipy.optimize import fmin_l_bfgs_b
import h5py

class model(object):
    """Define a linear model with within-class correlations.

    Args:
        y : :class:`~numpy:numpy.array`
            1D array of phenotype observations
        X : :class:`~numpy:numpy.array`
            Design matrix for the fixed mean effects.
        labels : :class:`~numpy:numpy.array`
            1D array of sample labels

    Returns:
        model : :class:`sibreg.model`

    """
    def __init__(self,y,X,labels, add_intercept = False):
        if y.shape[0] == X.shape[0] and X.shape[0] == labels.shape[0]:
            pass
        else:
            raise(ValueError('inconsistent sample sizes of response, covariates, and labels'))
        # Get sample size
        self.n = X.shape[0]
        if X.ndim == 1:
            X = X.reshape((self.n,1))
        if add_intercept:
            X = np.hstack((np.ones((self.n,1),dtype=X.dtype),X))
        self.X = X
        # Label mapping
        self.label_counts = dict()
        self.label_indices = dict()
        for l in range(0,labels.shape[0]):
            if labels[l] not in self.label_counts:
                self.label_counts[labels[l]]=1
                self.label_indices[labels[l]] = [l]
            else:
                self.label_counts[labels[l]]+=1
                self.label_indices[labels[l]].append(l)
        self.y_lab = dict()
        self.X_lab = dict()
        for label in self.label_indices.keys():
            self.y_lab[label]=y[self.label_indices[label]]
            self.X_lab[label]=X[self.label_indices[label],:]
        self.n_labels = len(self.y_lab.keys())
        # response
        self.y=y
        self.labels=labels

    def alpha_mle(self, tau, sigma2, compute_cov = False, xtx_out = False):
        """
        Compute the MLE of alpha given variance parameters

        Args:
            sigma2 : :class:`float`
                variance of model residuals
            tau : :class:`float`
                ratio of variance of model residuals to variance explained by mean differences between classes

        Returns:
            alpha : :class:`~numpy:numpy.array`
                MLE of alpha

        """
        X_T_X = np.zeros((self.X.shape[1],self.X.shape[1]),dtype = np.float64)
        X_T_y = np.zeros((self.X.shape[1]), dtype = np.float64)

        for label in self.y_lab.keys():
            sigma_u = sigma2/tau
            Sigma_lab = sigma_u*np.ones((self.label_counts[label],self.label_counts[label]))
            np.fill_diagonal(Sigma_lab,sigma_u+sigma2)
            Sigma_lab_inv = np.linalg.inv(Sigma_lab)
            X_T_X = X_T_X + np.dot(self.X_lab[label].T,Sigma_lab_inv.dot(self.X_lab[label]))
            X_T_y = X_T_y + np.dot(self.X_lab[label].T, Sigma_lab_inv.dot(self.y_lab[label]))

        if xtx_out:
            return [X_T_X,X_T_y.reshape((self.X.shape[1]))]
        else:
            alpha = np.linalg.solve(X_T_X,X_T_y)
            alpha = alpha.reshape((alpha.shape[0],))

            if compute_cov:
                alpha_cov = np.linalg.inv(X_T_X)
                return [alpha,alpha_cov]
            else:
                return alpha


    # Compute likelihood of data given beta, alpha
    def likelihood_and_gradient(self, sigma2, tau):
        """
        Compute the loss function, which is -2 times the likelihood along with its gradient

        Args:
            sigma2 : :class:`float`
                variance of model residuals
            tau : :class:`float`
                ratio of variance of model residuals to variance explained by mean differences between classes

        Returns:
            L, grad : :class:`float`
                loss function and gradient, divided by sample size

        """
        ## Likelihood
        alpha = self.alpha_mle(tau, sigma2)
        resid = self.y - self.X.dot(alpha)
        RSS = np.sum(np.square(resid))

        L = self.n * np.log(sigma2)+RSS/sigma2

        ## Gradient with respect to sigma2
        grad_sigma2 = self.n/sigma2-RSS/np.square(sigma2)

        ## Gradient with respect to tau
        grad_tau = 0

        for label in self.y_lab.keys():
            resid_label=resid[self.label_indices[label]]
            resid_sum = np.sum(resid_label)
            resid_square_sum = np.square(resid_sum)
            # Add to likelihood
            L = L - resid_square_sum/(sigma2*(tau+self.label_counts[label]))+np.log(1+self.label_counts[label]/tau)
            # Add to grad sigma2
            grad_sigma2+=resid_square_sum/(np.square(sigma2)*(tau+self.label_counts[label]))
            # Add to grad tau
            grad_tau+=(resid_square_sum/sigma2-self.label_counts[label]*(1+self.label_counts[label]/tau))/np.square(tau+self.label_counts[label])

        # Overall gradient vector
        grad = np.hstack((grad_sigma2,grad_tau))

        return L/self.n, grad/self.n

    def optimize_model(self,init_params):
        """
        Find the parameters that minimise the loss function for a given regularisation parameter

        Args:
            init_param : :class:`array`
                initial values for residual variance (sigma^2_epsilon) followed by ratio
                of residual variance to within-class variance (tau)

        Returns:
            optim : :class:`dict`
                dictionary with keys: 'success', whether optimisation was successful (bool);
                'warnflag', output of L-BFGS-B algorithm giving warnings; 'sigma2', MLE of
                residual variance; 'tau', MLE of ratio of residual variance to within-class variance;
                'likelihood', maximum of likelihood.

        """
        # Paramtere boundaries
        parbounds=[(0.00001, None),(0.00001, None)]
        # Optimize
        optimized = fmin_l_bfgs_b(func=lik_and_grad,x0=init_params,
                                args=(self.y, self.X, self.labels),
                                  bounds = parbounds)

        # Get MLE
        optim = {}
        optim['success'] = True
        optim['warnflag'] = optimized[2]['warnflag']
        if optim['warnflag'] != 0:
            print('Optimization unsuccessful.')
            optim['success'] = False
        optim['sigma2'] = optimized[0][0]
        optim['tau'] = optimized[0][1]
        # Get parameter covariance
        optim['likelihood'] = -0.5 * np.float64(self.n) * (optimized[1] + np.log(2 * np.pi))

        return optim

    def predict(self,X):
        """
        Predict new observations based on model regression coefficients

        Args:
            X : :class:`array`
                matrix of covariates to predict from

        Returns:
            y : :class:`array`
                predicted values
                
        """
        if hasattr(self,'alpha'):
            return X.dot(self.alpha)
        else:
            raise(AttributeError('Model does not have known regression coefficients. Try optimizing model first'))

    def set_alpha(self,alpha):
        self.alpha = alpha

def lik_and_grad(pars,*args):
    # Wrapper for function to pass to L-BFGS-B
    y, X, labels = args
    mod = model(y,X,labels)
    return mod.likelihood_and_gradient(pars[0],pars[1])

def simulate(n,alpha,sigma2,tau):
    """Simulate from a linear model with correlated observations within-class. The mean for each class
     is drawn from a normal distribution.

    Args:
        n : :class:`int`
            sample size
        alpha : :class:`~numpy:numpy.array`
            value of regression coefficeints
        sigma2 : :class:`float`
            variance of residuals
        tau : :class:`float`
            ratio of variance of residuals to variance of distribution of between individual means

    Returns:
        model : :class:`regrnd.model`
            linear model with repeated observations
            
    """
    c = alpha.shape[0]
    #X = np.random.randn((n * c)).reshape((n, c))
    X_cov = np.ones((c,c))
    np.fill_diagonal(X_cov,1.2)
    X = np.random.multivariate_normal(np.zeros((c)),X_cov,n).reshape((n, c))
    labels = np.random.choice(n//10,n)
    random_effects = np.sqrt(sigma2//tau)*np.random.randn(n)
    y = X.dot(alpha)+random_effects[labels-1]+np.random.randn(n)*np.sqrt(sigma2)
    return model(y,X,labels)

class gtarray(object):
    def __init__(self,garray,ids, sid = None, id_index = 0, alleles = None, fams= None, par_status = None):
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

        if par_status is not None:
            if par_status.shape[0] == ids.shape[0] and par_status.shape[1] == 2:
                self.par_status = par_status
            else:
                raise ValueError('Incompatible par status array')
        else:
            self.par_status = None

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
    gts_fams = np.zeros((gts_ids.shape[0]),dtype=gts_ids.dtype)
    gts_fams[ids_in_ped] = np.array([ped[ped_dict[x], 0] for x in gts_ids[ids_in_ped]])
    fams, counts = np.unique(gts_fams[ids_in_ped], return_counts=True)
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
    if remove_proband:
        gts_id_dict = make_id_dict(gts_ids)
    G_sib = np.zeros((ids.shape[0],gts.shape[1]),dtype = np.float32)
    for i in range(0,ids.shape[0]):
        fam_index = fams_dict[ids_fams[i]]
        G_sib[i,:] = fam_sums[fam_index,:]
        n_i = fam_counts[fam_index]
        if remove_proband:
            G_sib[i,:] = G_sib[i,:] - gts[gts_id_dict[ids[i]],:]
            n_i = n_i-1
        G_sib[i,:] = G_sib[i,:]/float(n_i)
    return gtarray(G_sib,ids)


def find_par_gts(pheno_ids,ped,gts_id_dict,fams = None):
    # Whether mother and father have observed/imputed genotypes
    par_status = np.zeros((pheno_ids.shape[0],2),dtype=int)
    par_status[:] = -1
    # Indices of obsered/imputed genotypes in relevant arrays
    gt_indices = np.zeros((pheno_ids.shape[0],3),dtype=int)
    gt_indices[:] = -1
    ## Build dictionaries
    ## Where each individual is in the pedigree
    ped_dict = make_id_dict(ped,1)
    if fams is not None:
        # Where the imputed data is for each family
        fam_dict = make_id_dict(fams)
        # Store family ID of each individual
    fam_labels = np.zeros((pheno_ids.shape[0]),dtype='U')
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
            if fams is not None:
                if ped_i[0] in fam_dict:
                    imp_index = fam_dict[ped_i[0]]
                    # Check if this is imputation of father, or mother, or both
                    if ped_i[4] == 'False':
                        gt_indices[i,1] = imp_index
                        par_status[i,0] = 1
                    if ped_i[5] == 'False':
                        gt_indices[i, 2] = imp_index
                        par_status[i, 1] = 1
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
    par_status, gt_indices, fam_labels = find_par_gts(ids,ped,gts_id_dict, fams = fams)
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
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, gts_id_dict, fams = fams)
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
    return gtarray(G,ids,sid, alleles = alleles, fams = fam_labels, par_status = par_status)

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
    null_model = model(y, X, fam_labels)
    null_optim = null_model.optimize_model(np.array([sigma_2_init,tau_init]))
    null_alpha = null_model.alpha_mle(null_optim['tau'],null_optim['sigma2'],compute_cov = True)
    # Return values
    return null_optim['sigma2'], null_optim['tau'], null_alpha


def get_alpha_mle(y,pgs,fam_labels, add_intercept = False):
    # Initialise var pars
    sigma_2_init = np.var(y) / 2.0
    rmodel = model(y, pgs, fam_labels, add_intercept = add_intercept)
    optim = rmodel.optimize_model(np.array([sigma_2_init,1]))
    print('Family variance estimate: '+str(round(optim['sigma2']/optim['tau'],4)))
    print('Residual variance estimate: ' + str(round(optim['sigma2'],4)))
    alpha = rmodel.alpha_mle(optim['tau'],optim['sigma2'],compute_cov = True)
    return alpha

