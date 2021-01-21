import numpy as np
import numpy.ma as ma
from pysnptools.snpreader import Bed, Pheno
from scipy.optimize import fmin_l_bfgs_b
import h5py, code

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

    def sigma_inv_root(self,tau,sigma2):
        sigma_u = sigma2 / tau
        sigma2_nsqrt = dict()
        famsizes = np.unique(list(self.label_counts.values()))
        sigma2_nsqrt[1] = np.power(sigma_u+sigma2,-0.5)
        famsizes = famsizes[famsizes>1]
        for famsize in famsizes:
            Sigma_lab = sigma_u*np.ones((famsize,famsize))
            np.fill_diagonal(Sigma_lab,sigma_u+sigma2)
            vals, vectors = np.linalg.eigh(Sigma_lab)
            vals = np.power(vals,0.25)
            vectors = vectors/vals
            sigma2_nsqrt[famsize] = vectors.dot(vectors.T)
        return sigma2_nsqrt

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
    """Define a genotype or PGS array that stores individual IDs, family IDs, and SNP information.

    Args:
        garray : :class:`~numpy:numpy.array`
            2 or 3 dimensional numpy array of genotypes/PGS values. First dimension is individuals. For a 2 dimensional array, the second dimension is SNPs or PGS values.
            For a 3 dimensional array, the second dimension indexes the individual and his/her relatives' genotypes (for example: proband, paternal, and maternal); and
            the third dimension is the SNPs.
        ids : :class:`~numpy:numpy.array`
            vector of individual IDs
        sid : :class:`~numpy:numpy.array`
            vector of SNP ids, equal in length size of last dimension of array
        alleles : :class:`~numpy:numpy.array`
            [L x 2] matrix of ref and alt alleles for the SNPs. L must match size of sid
        pos : :class:`~numpy:numpy.array`
            vector of SNP positions; must match size of sid
        chrom : :class:`~numpy:numpy.array`
            vector of SNP chromosomes; must match size of sid
        fams : :class:`~numpy:numpy.array`
            vector of family IDs; must match size of ids
        par_status : :class:`~numpy:numpy.array'
             [N x 2] numpy matrix that records whether parents have observed or imputed genotypes/PGS, where N matches size of ids.
             The first column is for the father of that individual; the second column is for the mother of that individual.
             If the parent is neither observed nor imputed, the value is -1; if observed, 0; and if imputed, 1.

    Returns:
        G : :class:`sibreg.gtarray`

    """
    def __init__(self, garray, ids, sid=None, alleles=None, pos=None, chrom=None, fams=None, par_status=None):
        if type(garray) == np.ndarray or type(garray)==np.ma.core.MaskedArray:
            if type(garray) == np.ndarray:
                self.gts = ma.array(garray,mask=np.isnan(garray))
            else:
                self.gts = garray
            self.shape = garray.shape
            self.ndim = garray.ndim
            self.dtype = garray.dtype
            self.freqs = None
        else:
            raise ValueError('Genotypes must be a numpy ndarray')
        if garray.shape[0] == ids.shape[0]:
            self.ids = ids
            self.id_dict = make_id_dict(ids)
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
            if self.sid is not None:
                if alleles.shape[0] == self.sid.shape[0]:
                    self.alleles = alleles
                else:
                    raise ValueError('Size of alleles does not match size of genotypes')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.alleles = None
        if pos is not None:
            if self.sid is not None:
                if pos.shape[0] == self.sid.shape[0]:
                    self.pos = pos
                else:
                    raise ValueError('Size of position vector does not match size of genotypes')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.pos = None
        if chrom is not None:
            if self.sid is not None:
                if chrom.shape[0] == self.sid.shape[0]:
                    self.chrom = chrom
                else:
                    raise ValueError('Size of position vector does not match size of genotypes')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.chrom = None

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

        if np.sum(self.gts.mask)>0:
            self.has_NAs = True
        else:
            self.has_NAs = False

    def compute_freqs(self):
        """
        Computes the frequencies of the SNPs. Stored in self.freqs.
        """
        if self.ndim == 2:
            self.freqs = ma.mean(self.gts,axis=0)/2.0
        elif self.ndim == 3:
            self.freqs = ma.mean(self.gts[:,0,:], axis=0) / 2.0

    def filter_snps(self, min_maf = 0, max_missing = 0):
        """
        Filter SNPs based on having minor allele frequency (MAF) greater than min_maf, and have % missing observations less than max_missing.
        """
        if self.freqs is None:
            self.compute_freqs()
        if self.ndim == 2:
            missingness = ma.mean(self.gts.mask,axis=0)
        elif self.ndim == 3:
            missingness = ma.mean(self.gts.mask,axis = (0,1))
        freqs_pass = np.logical_and(self.freqs > min_maf, self.freqs < (1 - min_maf))
        print(str(self.freqs.shape[0] - np.sum(freqs_pass)) + ' SNPs with MAF<' + str(min_maf))
        missingness_pass = 100 * missingness < max_missing
        print(str(self.freqs.shape[0] - np.sum(missingness_pass)) + ' SNPs with missingness >' + str(max_missing) + '%')
        filter_pass = np.logical_and(freqs_pass, missingness_pass)
        self.freqs = self.freqs[filter_pass]
        if self.ndim == 2:
            self.gts = self.gts[:,filter_pass]
        elif self.ndim == 3:
            self.gts = self.gts[:,:,filter_pass]
        self.shape = self.gts.shape
        if self.sid is not None:
            self.sid = self.sid[filter_pass]
            self.sid_dict = make_id_dict(self.sid)
        if self.pos is not None:
            self.pos = self.pos[filter_pass]
        if self.alleles is not None:
            self.alleles = self.alleles[filter_pass]
        if self.chrom is not None:
            self.chrom = self.chrom[filter_pass]

    def filter_ids(self,keep_ids):
        """
        Keep only individuals with ids given by keep_ids
        """
        in_ids = np.array([x in self.id_dict for x in keep_ids])
        n_filtered = np.sum(in_ids)
        if n_filtered==0:
            raise(ValueError('No individuals would be left after filtering'))
        else:
            print('After filtering, '+str(n_filtered)+' individuals remain')
            indices = np.array([self.id_dict[x] for x in keep_ids[in_ids]])
            if self.ndim == 2:
                self.gts = self.gts[indices, :]
            elif self.ndim == 3:
                self.gts = self.gts[indices, :, :]
            self.ids = self.ids[indices]
            self.id_dict = make_id_dict(self.ids)
            self.shape = self.gts.shape


    def mean_normalise(self):
        """
        This normalises the SNPs/PGS columns to have mean-zero.
        """
        if not self.mean_normalised:
            if self.gts.ndim == 2:
                self.gts = self.gts - ma.mean(self.gts,axis=0)
            elif self.gts.ndim==3:
                for i in range(0, self.gts.shape[1]):
                    self.gts[:, i, :] = self.gts[:, i, :] - ma.mean(self.gts[:, i, :], axis=0)
            self.mean_normalised = True

    def scale(self):
        """
        This normalises the SNPs/PGS columns to have variance 1.
        """
        if self.gts.ndim == 2:
            self.gts = self.gts/ma.std(self.gts, axis=0)
        elif self.gts.ndim == 3:
            for i in range(0, self.gts.shape[1]):
                self.gts[:, i, :] = self.gts[:, i, :]/ma.std(self.gts[:, i, :], axis=0)

    def fill_NAs(self):
        """
        This normalises the SNP columns to have mean-zero, then fills in NA values with zero.
        """
        if not self.mean_normalised:
            self.mean_normalise()
        NAs = np.sum(self.gts.mask, axis=0)
        self.gts[self.gts.mask] = 0
        self.gts.mask = False
        self.has_NAs = False
        return NAs


    def add(self,garray):
        """
        Adds another gtarray of the same dimension to this array and returns the sum. It matches IDs before summing.
        """
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

    def diagonalise(self,inv_root):
        """
        This will transform the genotype array based on the inverse square root of the phenotypic covariance matrix
        from the family based linear mixed model.
        """
        if self.fams is None:
            raise(ValueError('Family labels needed for diagonalization'))
        if not self.mean_normalised:
            self.mean_normalise()
        if self.has_NAs:
            self.fill_NAs()
        unique_fams, famsizes = np.unique(self.fams, return_counts = True)
        fam_indices = dict()
        # Transform
        for fam in unique_fams:
            fam_indices[fam] = np.where(self.fams == fam)[0]
            famsize = fam_indices[fam].shape[0]
            if self.ndim == 2:
                if famsize == 1:
                    self.gts[fam_indices[fam], :] = inv_root[1]*self.gts[fam_indices[fam],:]
                else:
                    self.gts[fam_indices[fam],:] = inv_root[famsize].dot(self.gts[fam_indices[fam],:])
            elif self.ndim == 3:
                if famsize == 1:
                    self.gts[fam_indices[fam], : , :] = inv_root[1]*self.gts[fam_indices[fam], : , :]
                else:
                    for j in range(self.shape[1]):
                        self.gts[fam_indices[fam],j, :] = inv_root[famsize].dot(self.gts[fam_indices[fam],j, :])
        self.fam_indices = fam_indices

class pgs(object):
    """Define a polygenic score based on a set of SNPs with weights and ref/alt allele pairs.

    Args:
        snp_ids : :class:`~numpy:numpy.array`
            vector of SNP ids
        snp_ids : :class:`~numpy:numpy.array`
            vector of weights of equal length to snp_ids
        alleles : :class:`~numpy:numpy.array`
            [L x 2] matrix of ref and alt alleles for the SNPs. L must match size of snp_ids

    Returns:
        pgs : :class:`sibreg.pgs`

    """
    def __init__(self,snp_ids,weights,alleles):
        if snp_ids.shape[0] == weights.shape[0] and alleles.shape[0] == weights.shape[0] and alleles.shape[1]==2:
            self.snp_ids = snp_ids
            self.snp_dict = make_id_dict(snp_ids)
            self.weights = weights
            self.alleles = alleles
        else:
            raise ValueError('All inputs must have the same dimension')

    def compute(self,garray, cols = None):
        """Compute polygenic score values from a given genotype array. Finds the SNPs in the genotype array
        that have weights in the pgs and matching alleles, and computes the PGS based on these SNPs and the
        weights after allele-matching.


        Args:
            garray : :class:`sbreg.gtarray`
                genotype array to compute PGS values for
            cols : :class:`numpy:numpy.array`
                names to give the columns in the output gtarray

        Returns:
            pg : :class:`sibreg.gtarray`
                2d gtarray with PGS values. If a 3d gtarray is input, then each column corresponds to
                the second dimension on the input gtarray (for example, individual, paternal, maternal PGS).
                If a 2d gtarray is input, then there will be only one column in the output gtarray. The
                names given in 'cols' are stored in 'sid' attribute of the output.

        """
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
    """
    Make a dictionary that maps from the values in the given column (col) to their row-index in the input array
    """
    if len(x.shape)>1:
        x = x[:,col]
    id_dict = {}
    for i in range(0,x.shape[0]):
        id_dict[x[i]] = i
    return id_dict

def convert_str_array(x):
    """
    Convert an ascii array to unicode array (UTF-8)
    """
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.decode('UTF-8') for y in x])
    return x_out.reshape(x_shape)


def encode_str_array(x):
    """
    Encode a unicode array as an ascii array
    """
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.encode('ascii') for y in x])
    return x_out.reshape(x_shape)

def find_individuals_with_sibs(ids,ped,gts_ids, return_ids_only = False):
    """
    Used in get_gts_matrix and get_fam_means to find the individuals in ids that have genotyped siblings.
    """
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

def get_fam_means(ids,ped,gts,gts_ids,remove_proband = True, return_famsizes = False):
    """
    Used in get_gts_matrix to find the mean genotype in each sibship (family) for each SNP or for a PGS.
    The gtarray that is returned is indexed based on the subset of ids provided from sibships of size 2 or greater.
    If remove_proband=True, then the genotype/PGS of the index individual is removed from the fam_mean given for that individual.
    """
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
    if return_famsizes:
        return [gtarray(G_sib, ids),fam_counts,fam_sums]
    else:
        return gtarray(G_sib,ids)


def find_par_gts(pheno_ids, ped, fams, gts_id_dict):
    """
    Used in get_gts_matrix to find whether individuals have imputed or observed parental genotypes, and to
    find the indices of the observed/imputed parents in the observed/imputed genotype arrays.
    'par_status' codes whether an individual has parents that are observed or imputed or neither.
    'gt_indices' records the relevant index of the parent in the observed/imputed genotype arrays
    'fam_labels' records the family of the individual based on the pedigree
    """
    # Whether mother and father have observed/imputed genotypes
    par_status = np.zeros((pheno_ids.shape[0],2),dtype=int)
    par_status[:] = -1
    # Indices of obsered/imputed genotypes in relevant arrays
    gt_indices = np.zeros((pheno_ids.shape[0],3),dtype=int)
    gt_indices[:] = -1
    ## Build dictionaries
    # Where each individual is in the pedigree
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
            ped_i = ped[ped_dict[pheno_ids[i]], :]
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
                if ped_i[4] == 'False' and not par_status[i,0] == 0:
                    gt_indices[i, 1] = imp_index
                    par_status[i, 0] = 1
                if ped_i[5] == 'False' and not par_status[i,0] == 0:
                    gt_indices[i, 2] = imp_index
                    par_status[i, 1] = 1
    return par_status, gt_indices, fam_labels

def make_gts_matrix(gts,imp_gts,par_status,gt_indices, parsum = False):
    """
    Used in get_gts_matrix to construct the family based genotype matrix given
    observed/imputed genotypes. 'gt_indices' has the indices in the observed/imputed genotype arrays;
    and par_status codes whether the parents are observed (0) or imputed (1).
    """
    if np.min(gt_indices)<0:
        raise(ValueError('Missing genotype index'))
    N = gt_indices.shape[0]
    if parsum:
        gdim = 2
    else:
        gdim = 3
    G = np.zeros((N,gdim,gts.shape[1]),np.float32)
    # Proband genotypes
    G[:,0,:] = gts[gt_indices[:,0],:]
    # Paternal genotypes
    G[par_status[:, 0] == 0, 1 ,:] = gts[gt_indices[par_status[:, 0] == 0, 1], :]
    G[par_status[:, 0] == 1, 1, :] = imp_gts[gt_indices[par_status[:, 0] == 1, 1], :]
    # Maternal genotypes
    if parsum:
        G[par_status[:, 1] == 0, 1, :] += gts[gt_indices[par_status[:, 1] == 0, 2], :]
        G[par_status[:, 1] == 1, 1, :] += imp_gts[gt_indices[par_status[:, 1] == 1, 2], :]
    else:
        G[par_status[:, 1] == 0, 2, :] = gts[gt_indices[par_status[:, 1] == 0, 2], :]
        G[par_status[:, 1] == 1, 2, :] = imp_gts[gt_indices[par_status[:, 1] == 1, 2], :]
    return G


def get_gts_matrix(par_gts_f,gts_f,snp_ids = None,ids = None, sib = False, compute_controls = False, parsum = False):
    """Reads observed and imputed genotypes and constructs a family based genotype matrix for the individuals with
    observed/imputed parental genotypes, and if sib=True, at least one genotyped sibling.

    Args:
        par_gts_f : :class:`str`
            path to HDF5 file with imputed parental genotypes
        gts_f : :class:`str`
            path to bed file with observed genotypes
        snp_ids : :class:`numpy.ndarray`
            If provided, only obtains the subset of SNPs specificed that are present in both imputed and observed genotypes
        ids : :class:`numpy.ndarray`
            If provided, only obtains the ids with observed genotypes and imputed/observed parental genotypes (and observed sibling genotypes if sib=True)
        sib : :class:`bool`
            Retrieve genotypes for individuals with at least one genotyped sibling along with the average of their siblings' genotypes and observed/imputed parental genotypes. Default False.
        compute_controls : :class:`bool`
            Compute polygenic scores for control families (families with observed parental genotypes set to missing). Default False.
        parsum : :class:`bool`
            Return the sum of maternal and paternal observed/imputed genotypes rather than separate maternal/paternal genotypes. Default False.

    Returns:
        G : :class:`sibreg.gtarray`
            Genotype array for the subset of genotyped individuals with complete imputed/obsereved parental genotypes. The array is [N x k x L], where
            N is the number of individuals; k depends on whether sib=True and whether parsum=True; and  L is the number of SNPs. If sib=False and parsum=False,
            then k=3 and this axis indexes individual's genotypes, individual's father's imputed/observed genotypes, individual's mother's imputed/observed genotypes.
            If sib=True and parsum=False, then k=4, and this axis indexes the individual, the sibling, the paternal, and maternal genotypes in that order. If parsum=True and sib=False,
            then k=2, and this axis indexes the individual and sum of paternal and maternal genotypes; etc.
            If compute_controls=True, then a list is returned, where the first element is as above, and the following elements give equivalent genotyping arrays for control families where the mother has been set
            to missing, the father has been set to missing, and both parents have been set to missing.

    """
    ####### Find parental status #######
    ### Imputed parental file ###
    par_gts_f = h5py.File(par_gts_f,'r')
    # Get pedigree
    ped = convert_str_array(np.array(par_gts_f['pedigree']))
    ped = ped[1:ped.shape[0],:]
    # Remove control families
    controls = np.array([x[0]=='_' for x in ped[:,0]])
    G = [get_gts_matrix_given_ped(ped[np.logical_not(controls),:],par_gts_f,gts_f,
                                  snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum)]
    if compute_controls:
        G.append(get_gts_matrix_given_ped(ped[np.array([x[0:3]=='_p_' for x in ped[:,0]]),],par_gts_f,gts_f,
                                  snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum))
        G.append(
            get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_m_' for x in ped[:, 0]]),], par_gts_f, gts_f,
                                  snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum))
        G.append(
            get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_o_' for x in ped[:, 0]]),], par_gts_f, gts_f,
                                  snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum))
        return G
    else:
        return G[0]


def get_indices_given_ped(ped, fams, gts_ids, ids=None, sib=False):
    """
    Used in get_gts_matrix_given_ped to get the ids of individuals with observed/imputed parental genotypes and, if sib=True, at least one genotyped sibling.
    It returns those ids along with the indices of the relevant individuals and their first degree relatives in the observed genotypes (observed indices),
    and the indices of the imputed parental genotypes for those individuals.
    """
    # Made dictionary for observed genotypes
    gts_id_dict = make_id_dict(gts_ids)
    # If IDs not provided, use all individuals with observed genotypes
    if ids is None:
        ids = gts_ids
    # Find individuals with siblings
    if sib:
        ids = find_individuals_with_sibs(ids, ped, gts_ids, return_ids_only=True)
        print('Found ' + str(ids.shape[0]) + ' individuals with genotyped siblings')
    ### Find parental status
    print('Checking for observed/imputed parental genotypes')
    par_status, gt_indices, fam_labels = find_par_gts(ids, ped, fams, gts_id_dict)
    # Find which individuals can be used
    none_missing = np.min(par_status, axis=1)
    none_missing = none_missing >= 0
    N = np.sum(none_missing)
    if N == 0:
        raise ValueError(
            'No individuals with phenotype observations and complete observed/imputed genotype observations')
    print(str(N) + ' individuals with phenotype observations and complete observed/imputed genotypes observations')
    # Take those that can be used
    gt_indices = gt_indices[none_missing, :]
    par_status = par_status[none_missing, :]
    ids = ids[none_missing]
    # Find indices of individuals and their parents in observed genotypes
    observed_indices = np.sort(np.unique(np.hstack((gt_indices[:, 0],
                                                    gt_indices[par_status[:, 0] == 0, 1],
                                                    gt_indices[par_status[:, 1] == 0, 2]))))
    # Get indices of imputed parents
    imp_indices = np.sort(np.unique(np.hstack((gt_indices[par_status[:, 0] == 1, 1],
                                               gt_indices[par_status[:, 1] == 1, 2]))))
    # Return ids with imputed/observed parents
    return ids, observed_indices, imp_indices


def match_observed_and_imputed_snps(gts_f, par_gts_f, bim, snp_ids=None):
    """
    Used in get_gts_matrix_given_ped to match observed and imputed SNPs and return SNP information on shared SNPs.
    Removes SNPs that have duplicated SNP ids.
    in_obs_sid contains the SNPs in the imputed genotypes that are present in the observed SNPs
    obs_sid_index contains the index in the observed SNPs of the common SNPs
    """
    # Match SNPs from imputed and observed and restrict to those in list
    if snp_ids is None:
        snp_ids = gts_f.sid
    # Get bim info
    alleles = np.loadtxt(bim, dtype='U', usecols=(4, 5))
    pos = np.loadtxt(bim, dtype=int, usecols=3)
    chromosome = np.loadtxt(bim, dtype=int, usecols=0)
    # Remove duplicate ids
    unique_snps, snp_indices, snp_counts = np.unique(snp_ids, return_index=True, return_counts=True)
    snp_set = set(snp_ids[snp_indices[snp_counts == 1]])
    if len(snp_set) < snp_ids.shape[0]:
        print(str(snp_ids.shape[0]-len(snp_set))+' SNPs with duplicate IDs removed')
    # Read and match SNP ids
    imp_bim = convert_str_array(np.array(par_gts_f['bim_values']))
    imp_sid = imp_bim[:, 1]
    obs_sid = gts_f.sid
    obs_sid_dict = make_id_dict(obs_sid)
    in_obs_sid = np.zeros((imp_sid.shape[0]), dtype=bool)
    obs_sid_index = np.zeros((imp_sid.shape[0]), dtype=int)
    for i in range(0, imp_sid.shape[0]):
        if imp_sid[i] in obs_sid_dict and imp_sid[i] in snp_set:
            in_obs_sid[i] = True
            obs_sid_index[i] = obs_sid_dict[imp_sid[i]]
    if np.sum(in_obs_sid) == 0:
        raise ValueError('No SNPs in common between imputed and observed data')
    obs_sid_index = obs_sid_index[in_obs_sid]
    sid = imp_sid[in_obs_sid]
    alleles = alleles[obs_sid_index, :]
    chromosome = chromosome[obs_sid_index]
    pos = pos[obs_sid_index]
    return chromosome, sid, pos, alleles, in_obs_sid, obs_sid_index

def get_gts_matrix_given_ped(ped, par_gts_f, gts_f, snp_ids=None, ids=None, sib=False, parsum=False):
    """
    Used in get_gts_matrix: see get_gts_matrix for documentation
    """
    ### Genotype file ###
    bim = gts_f.split('.bed')[0] + '.bim'
    gts_f = Bed(gts_f,count_A1=True)
    # get ids of genotypes and make dict
    gts_ids = gts_f.iid[:, 1]
    # Get families with imputed parental genotypes
    fams = convert_str_array(np.array(par_gts_f['families']))
    ### Find ids with observed/imputed parents and indices of those in observed/imputed data
    ids, observed_indices, imp_indices = get_indices_given_ped(ped, fams, gts_ids, ids=ids, sib=sib)
    ### Match observed and imputed SNPs ###
    print('Matching observed and imputed SNPs')
    chromosome, sid, pos, alleles, in_obs_sid, obs_sid_index = match_observed_and_imputed_snps(gts_f, par_gts_f, bim, snp_ids=snp_ids)
    # Read imputed parental genotypes
    print('Reading imputed parental genotypes')
    if (imp_indices.shape[0]*in_obs_sid.shape[0]) < (np.sum(in_obs_sid)*fams.shape[0]):
        imp_gts = np.array(par_gts_f['imputed_par_gts'][imp_indices, :])
        imp_gts = imp_gts[:,np.arange(in_obs_sid.shape[0])[in_obs_sid]]
    else:
        imp_gts = np.array(par_gts_f['imputed_par_gts'][:,np.arange(in_obs_sid.shape[0])[in_obs_sid]])
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
        if parsum:
            G = np.zeros((ids.shape[0], 3, gts.shape[1]), dtype=np.float32)
            G[:, np.array([0, 2]), :] = make_gts_matrix(gts, imp_gts, par_status, gt_indices, parsum=parsum)
        else:
            G = np.zeros((ids.shape[0],4,gts.shape[1]), dtype=np.float32)
            G[:,np.array([0,2,3]),:] = make_gts_matrix(gts, imp_gts, par_status, gt_indices, parsum=parsum)
        G[:,1,:] = get_fam_means(ids, ped, gts, gts_ids, remove_proband=True).gts
    else:
        G = make_gts_matrix(gts, imp_gts, par_status, gt_indices, parsum=parsum)
    del gts
    del imp_gts
    return gtarray(G, ids, sid, alleles=alleles, pos=pos, chrom=chromosome, fams=fam_labels, par_status=par_status)

def compute_pgs(par_gts_f, gts_f, pgs, sib=False, compute_controls=False):
    """Compute a polygenic score (PGS) for the individuals with observed genotypes and observed/imputed parental genotypes.

    Args:
        par_gts_f : :class:`str`
            path to HDF5 file with imputed parental genotypes
        gts_f : :class:`str`
            path to bed file with observed genotypes
        pgs : :class:`sibreg.pgs`
            the PGS, defined by the weights for a set of SNPs and the alleles of those SNPs
        sib : :class:`bool`
            Compute the PGS for genotyped individuals with at least one genotyped sibling and observed/imputed parental genotypes. Default False.
        compute_controls : :class:`bool`
            Compute polygenic scores for control families (families with observed parental genotypes set to missing). Default False.

    Returns:
        pg : :class:`sibreg.gtarray`
            Return the polygenic score as a genotype array with columns: individual's PGS, mean of their siblings' PGS, observed/imputed paternal PGS,
            observed/imputed maternal PGS

    """
    G = get_gts_matrix(par_gts_f, gts_f, snp_ids=pgs.snp_ids, sib=sib, compute_controls=compute_controls)
    if sib:
        cols = np.array(['proband', 'sibling', 'paternal', 'maternal'])
    else:
        cols = np.array(['proband', 'paternal', 'maternal'])
    if compute_controls:
        return [pgs.compute(x,cols) for x in G]
    else:
        return pgs.compute(G,cols)


def fit_sibreg_model(y, X, fam_labels, add_intercept=False, tau_init=1, return_model=True, return_vcomps=True, return_fixed=True):
    """Compute the MLE for the fixed effects in a family-based linear mixed model.

    Args:
        y : :class:`~numpy:numpy.array`
            vector of phenotype values
        X: :class:`~numpy:numpy.array`
            regression design matrix for fixed effects
        fam_labels : :class:`~numpy:numpy.array`
            vector of family labels: residual correlations in y are modelled between family members (that share a fam_label)
        add_intercept : :class:`bool`
            whether to add an intercept to the fixed effect design matrix

    Returns:
       model : :class:`sibreg.model`
            the sibreg model object, if return_model=True
       vcomps: :class:`float`
            the MLEs for the variance parameters: sigma2 (residual variance) and tau (ratio between sigma2 and family variance), if return_vcomps=True
       alpha : :class:`~numpy:numpy.array`
            MLE of fixed effects, if return_fixed=True
       alpha_cov : :class:`~numpy:numpy.array`
            sampling variance-covariance matrix for MLE of fixed effects, if return_fixed=True

    """
    # Optimize  model
    sigma_2_init = np.var(y)*tau_init/(1+tau_init)
    null_model = model(y, X, fam_labels, add_intercept=add_intercept)
    null_optim = null_model.optimize_model(np.array([sigma_2_init,tau_init]))
    # Create return list
    return_list = []
    if return_model:
        return_list.append(null_model)
    if return_vcomps:
        return_list += [null_optim['sigma2'], null_optim['tau']]
    if return_fixed:
        return_list += null_model.alpha_mle(null_optim['tau'], null_optim['sigma2'], compute_cov=True)
    return return_list

def read_phenotype(phenofile, missing_char = 'NA', phen_index = 1):
    """Read a phenotype file and remove missing values.

    Args:
        phenofile : :class:`str`
            path to plain text phenotype file with columns FID, IID, phenotype1, phenotype2, ...
        missing_char : :class:`str`
            The character that denotes a missing phenotype value; 'NA' by default.
        phen_index : :class:`int`
           The index of the phenotype (counting from 1) if multiple phenotype columns present in phenofile

    Returns:
        y : :class:`~numpy:numpy.array`
            vector of non-missing phenotype values from specified column of phenofile
        pheno_ids: :class:`~numpy:numpy.array`
            corresponding vector of individual IDs (IID)
    """
    pheno = Pheno(phenofile, missing=missing_char).read()
    y = np.array(pheno.val)
    pheno_ids = np.array(pheno.iid)[:,1]
    if y.ndim == 1:
        pass
    elif y.ndim == 2:
        y = y[:, phen_index - 1]
    else:
        raise (ValueError('Incorrect dimensions of phenotype array'))
    # Remove y NAs
    y_not_nan = np.logical_not(np.isnan(y))
    if np.sum(y_not_nan) < y.shape[0]:
        y = y[y_not_nan]
        pheno_ids = pheno_ids[y_not_nan]
    print('Number of non-missing phenotype observations: ' + str(y.shape[0]))
    return y, pheno_ids

def match_phenotype(G,y,pheno_ids):
    """Match a phenotype to a genotype array by individual IDs.

    Args:
        G : :class:`gtarray`
            genotype array to match phenotype to
        y : :class:`~numpy:numpy.array`
            vector of phenotype values
        pheno_ids: :class:`~numpy:numpy.array`
            vector of individual IDs corresponding to phenotype vector, y

    Returns:
       y : :class:`~numpy:numpy.array`
            vector of phenotype values matched by individual IDs to the genotype array

    """
    in_G_dict = np.array([x in G.id_dict for x in pheno_ids])
    y = y[in_G_dict]
    pheno_ids = pheno_ids[in_G_dict]
    pheno_id_dict = make_id_dict(pheno_ids)
    y = y[[pheno_id_dict[x] for x in G.ids]]
    return y
