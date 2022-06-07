import numpy as np
import numpy.ma as ma
from snipar.utilities import make_id_dict
from scipy.optimize import fmin_l_bfgs_b
from numba import njit, prange

class gtarray(object):
    """Define a genotype or PGS array that stores individual IDs, family IDs, and SNP information.

    Args:
        garray : :class:`~numpy:numpy.array`
            2 or 3 dimensional numpy array of genotypes/PGS values. First dimension is individuals. For a 2 dimensional array, the second dimension is SNPs or PGS values.
            For a 3 dimensional array, the second dimension indexes the individual and his/her relatives' genotypes (for example: proband, paternal, and maternal); and
            the third dimension is the SNPs.
        ids : :class:`~numpy:numpy.array`
            vector of individual IDs of same length as first dimension of garray
        sid : :class:`~numpy:numpy.array`
            vector of SNP ids, equal in length, L, to last dimension of array
        alleles : :class:`~numpy:numpy.array`
            [L x 2] matrix of ref and alt alleles for the SNPs. L must match size of sid
        pos : :class:`~numpy:numpy.array`
            vector of SNP positions; must match size of sid
        chrom : :class:`~numpy:numpy.array`
            vector of SNP chromosomes; must match size of sid
        map : :class:`~numpy:numpy.array`
            vector of SNP chromosomes; must match size of sid
        fams : :class:`~numpy:numpy.array`
            vector of family IDs; must match size of ids
        par_status : :class:`~numpy:numpy.array'
             [N x 2] numpy matrix that records whether parents have observed or imputed genotypes/PGS, where N matches size of ids.
             The first column is for the father of that individual; the second column is for the mother of that individual.
             If the parent is neither observed nor imputed, the value is -1; if observed, 0; and if imputed, 1.

    Returns:
        G : :class:`snipar.gtarray`

    """
    def __init__(self, garray, ids, sid=None, alleles=None, pos=None, chrom=None, map=None, error_probs=None, fams=None, par_status=None):
        if type(garray) == np.ndarray or type(garray) == np.ma.core.MaskedArray:
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
                    raise ValueError('Size of map does not match number of SNPs')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.chrom = None
        if map is not None:
            if self.sid is not None:
                if map.shape[0] == self.sid.shape[0]:
                    self.map = map
                else:
                    raise ValueError('Size of map does not match number of SNPs')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.map = None
        if error_probs is not None:
            if self.sid is not None:
                if error_probs.shape[0] == self.sid.shape[0]:
                    self.error_probs = error_probs
                else:
                    raise ValueError('Size of map does not match number of SNPs')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.error_probs = None
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

        self.info = None

    def compute_freqs(self):
        """
        Computes the frequencies of the SNPs. Stored in self.freqs.
        """
        if self.ndim == 2:
            self.freqs = ma.mean(self.gts,axis=0)/2.0
        elif self.ndim == 3:
            self.freqs = ma.mean(self.gts[:,0,:], axis=0) / 2.0

    def filter(self, filter_pass):
        if self.freqs is not None:
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
        if self.map is not None:
            self.map = self.map[filter_pass]
        if self.error_probs is not None:
            self.error_probs = self.error_probs[filter_pass]

    def filter_maf(self, min_maf = 0.01, verbose=False):
        """
        Filter SNPs based on having minor allele frequency (MAF) greater than min_maf, and have % missing observations less than max_missing.
        """
        if self.freqs is None:
            self.compute_freqs()
        freqs_pass = np.logical_and(self.freqs > min_maf, self.freqs < (1 - min_maf))
        if verbose:
            print(str(self.freqs.shape[0] - np.sum(freqs_pass)) + ' SNPs with MAF<' + str(min_maf))
        self.filter(freqs_pass)

    def filter_missingness(self, max_missing = 5, verbose=False):
        if self.ndim == 2:
            missingness = ma.mean(self.gts.mask,axis=0)
        elif self.ndim == 3:
            missingness = ma.mean(self.gts.mask,axis = (0,1))
        missingness_pass = 100 * missingness < max_missing
        if verbose:
            print(str(self.freqs.shape[0] - np.sum(missingness_pass)) + ' SNPs with missingness >' + str(max_missing) + '%')
        self.filter(missingness_pass)

    def compute_info(self):
        if self.freqs is None:
            self.compute_freqs()
        if self.ndim == 2:
            self.variances = np.var(self.gts, axis=0)
        elif self.ndim==3:
            self.variances = np.var(self.gts[:,0,:], axis=0)
        self.info = self.variances/(2.0*self.freqs*(1-self.freqs))

    def filter_info(self, min_info = 0.99, verbose=False):
        if self.info is None:
            self.compute_info()
        info_pass = self.info > min_info
        if verbose:
            print(str(self.info.shape[0] - np.sum(info_pass)) + ' SNPs with INFO <' + str(min_info))
        self.filter(info_pass)

    def filter_ids(self,keep_ids, verbose=False):
        """
        Keep only individuals with ids given by keep_ids
        """
        in_ids = np.array([x in self.id_dict for x in keep_ids])
        n_filtered = np.sum(in_ids)
        if n_filtered==0:
            raise(ValueError('No individuals would be left after filtering'))
        else:
            if verbose:
                print('After filtering, '+str(n_filtered)+' individuals remain')
            indices = np.array([self.id_dict[x] for x in keep_ids[in_ids]])
            if self.ndim == 2:
                self.gts = self.gts[indices, :]
            elif self.ndim == 3:
                self.gts = self.gts[indices, :, :]
            self.ids = self.ids[indices]
            self.id_dict = make_id_dict(self.ids)
            self.shape = self.gts.shape
            if self.fams is not None:
                self.fams = self.fams[indices]

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
        if len(common_ids) == 0:
            raise ValueError('No IDs in common')
        self_index = np.array([self.id_dict[x] for x in common_ids])
        other_index = np.array([garray.id_dict[x] for x in common_ids])

        # Out
        if self.ids.ndim == 1:
            ids_out = self.ids[self_index]
        else:
            ids_out = self.ids[self_index, :]

        if self.gts.ndim ==2:
            add_gts = self.gts[self_index, :]+garray.gts[other_index, :]
        else:
            add_gts = self.gts[self_index, :, :] + garray.gts[other_index, :, :]

        return gtarray(add_gts, ids_out, self.sid, alleles=self.alleles, fams=self.fams[self_index], par_status=self.par_status[self_index,:])

    def diagonalise(self,inv_root):
        """
        This will transform the genotype array based on the inverse square root of the phenotypic covariance matrix
        from the family based linear mixed model.
        """
        if self.fams is None:
            raise(ValueError('Family labels needed for diagonalization'))
        if not self.mean_normalised:
            self.mean_normalise()
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

    def estimate_r(self):
        # Check pgs columns
        if 'paternal' in self.sid:
            paternal_index = np.where(self.sid=='paternal')[0][0]
        else:
            raise(ValueError('No paternal PGS column found'))
        if 'maternal' in self.sid:
            maternal_index = np.where(self.sid=='maternal')[0][0]
        else:
            raise(ValueError('No maternal PGS column found'))
        # count fams
        fams = np.unique(self.fams, return_counts=True)
        # Mapping of sibships to pgs rows
        fam_dict = {}
        for fam in fams[0]:
            fam_dict[fam] = np.where(self.fams==fam)[0]
        # record genotype status of parents in each fam
        par_status_fams = np.zeros((fams[0].shape[0],2),dtype=int)
        for i in range(fams[0].shape[0]):
            par_status_fams[i,:] = self.par_status[fam_dict[fams[0][i]][0]]
        # get family sizes including parents
        fsizes = fams[1]+np.sum(par_status_fams==0,axis=1)
        ## pgs matrix with observed sibs and parents for each fam
        pgs_fam = np.zeros((fams[0].shape[0],np.max(fsizes)))
        pgs_fam[:] = np.nan
        # record whether sib or parent
        is_sib = np.zeros((fams[0].shape[0],np.max(fsizes)),dtype=bool)
        is_father = np.zeros((fams[0].shape[0],np.max(fsizes)),dtype=bool)
        is_mother = np.zeros((fams[0].shape[0],np.max(fsizes)),dtype=bool)
        # populate
        for i in range(fams[0].shape[0]):
            # Fill in sibs
            sib_indices = fam_dict[fams[0][i]]
            pgs_fam[i,0:sib_indices.shape[0]] = self.gts[sib_indices,0]
            is_sib[i,0:sib_indices.shape[0]] = True
            npar = 0
            if par_status_fams[i,0] == 0:
                pgs_fam[i,sib_indices.shape[0]] = self.gts[sib_indices[0],paternal_index]
                is_father[i,sib_indices.shape[0]] = True
                npar = 1
            if par_status_fams[i,1] == 0:
                pgs_fam[i,sib_indices.shape[0]+npar] = self.gts[sib_indices[0],maternal_index]
                is_mother[i,sib_indices.shape[0]+npar] = True
        is_parent = np.logical_or(is_father,is_mother)
        # normalize
        pgs_fam[is_sib] = (pgs_fam[is_sib]-np.mean(pgs_fam[is_sib]))/np.std(pgs_fam[is_sib])
        pgs_fam[is_mother] = (pgs_fam[is_mother]-np.mean(pgs_fam[is_mother]))/np.std(pgs_fam[is_mother])
        pgs_fam[is_father] = (pgs_fam[is_father]-np.mean(pgs_fam[is_father]))/np.std(pgs_fam[is_father])
        ### find correlation between maternal and paternal pgis
        # grid search over r
        print('Finding MLE for correlation between parents scores')
        ## Initialize with correlation from sibs and between parents
        # correlation between parents
        bpg = np.sum(self.par_status==0,axis=1)==2
        n_bpg = np.sum(bpg)
        if n_bpg>0:
            r_bpg = np.corrcoef(self.gts[bpg,1],self.gts[bpg,2])[0,1]
        else:
            r_bpg = 0
        # Correlation from sibs
        sib_2 = np.sum(is_sib,axis=1)>1
        n_sib_2 = np.sum(sib_2)
        if n_sib_2>0:
            r_sib = 2*np.corrcoef(pgs_fam[sib_2,0],pgs_fam[sib_2,1])[0,1]-1
        else:
            r_sib = 0
        # Initialize at weighted average
        r_init = n_bpg*r_bpg+n_sib_2*r_sib
        r_init = r_init/(n_bpg+n_sib_2)
        # Find MLE
        optimized = fmin_l_bfgs_b(func=pgs_cor_lik, x0=r_init,
                                args=(pgs_fam, is_sib, is_parent),
                                approx_grad=True,
                                bounds=[(-0.999,0.999)])
        #print('r='+str(round(r_mle,3)))
        if optimized[2]['warnflag']==0:
            r=optimized[0][0]
        else:
            print('Could not find MLE for correlation. Returning weighted average from siblings and parents.')
            r=r_init
        # fam size dict
        fsizes = dict(zip(list(fams[0]),list(fams[1])))
        return r, fsizes
    
    def am_adj(self):
        print('Estimating correlation between maternal and paternal PGSs assuming equilibrium')
        r, fsizes = self.estimate_r()
        print('Estimated correlation between maternal and paternal PGSs: '+str(round(r,4)))
        if r>0:    
            # Check pgs columns
            if 'paternal' in self.sid:
                paternal_index = np.where(self.sid=='paternal')[0][0]
            else:
                raise(ValueError('No paternal PGS column found'))
            if 'maternal' in self.sid:
                maternal_index = np.where(self.sid=='maternal')[0][0]
            else:
                raise(ValueError('No maternal PGS column found'))
            # Adjust imputed parental PGSs
            npar = np.sum(self.par_status==0,axis=1)
            print('Adjuting imputed PGSs for assortative mating')
            for i in range(self.gts.shape[0]):
                # No parents genotyped
                if npar==0:
                    self.gts[i,[paternal_index, maternal_index]] = npg_am_adj(r,fsizes[self.fams[i]])*self.gts[i,[paternal_index, maternal_index]]
                # One parent genotyped
                if npar==1:
                    # Father imputed
                    if self.par_status[i,paternal_index] == 1:
                        self.gts[i,paternal_index] = opg_am_adj(self.gts[i,paternal_index],self.gts[i,maternal_index],r,fsizes[self.fams[i]])
                    # Mother imputed
                    if self.par_status[i,maternal_index] == 1:
                        self.gts[i,maternal_index] = opg_am_adj(self.gts[i,maternal_index],self.gts[i,paternal_index],r,fsizes[self.fams[i]])
        else:
            print('Estimated correlation is negative, so not performing assortative mating adjustment')
        return r
        
def opg_am_adj(pgi_imp, pgi_obs, r, n):
    rcoef = r/((2**n)*(1+r)-r)
    return rcoef*pgi_obs+(1+rcoef)*pgi_imp

def npg_am_adj(r,n):
  rcoef = (1+r)/(1+(1-(1/2)**(n-1))*r)
  return rcoef

@njit
def pgs_corr_matrix(r,is_sib_fam,is_parent_fam):
    n_sib = np.sum(is_sib_fam)
    n_par = np.sum(is_parent_fam)
    # Full matrix
    R = np.zeros((n_sib+n_par,n_sib+n_par), dtype=np.float_)
    np.fill_diagonal(R, 1.0)
    # sib submatrix
    r_sib = (1+r)/2
    R_sib = r_sib*np.ones((n_sib,n_sib),dtype=np.float_)
    np.fill_diagonal(R_sib,1.0)
    R[0:n_sib,0:n_sib] = R_sib
    # sib to parent
    if n_par>0:
        R[0:n_sib,n_sib:(n_sib+n_par)] = r_sib
        R[n_sib:(n_sib+n_par),0:n_sib] = r_sib
    # parent-to-parent
    if n_par>1:
        R[n_sib+n_par-1,n_sib+n_par-2] = r
        R[n_sib+n_par-2,n_sib+n_par-1] = r
    # return 
    return R

@njit
def pgs_corr_likelihood_fam(r,pg_fam,is_sib_fam,is_parent_fam):
    sib_or_parent = np.logical_or(is_sib_fam,is_parent_fam)
    R = pgs_corr_matrix(r,is_sib_fam,is_parent_fam)
    slogdet_R = np.linalg.slogdet(R)
    pg_vec = pg_fam[sib_or_parent].reshape((np.sum(sib_or_parent),1))
    L_fam = slogdet_R[0]*slogdet_R[1]+pg_vec.T @ np.linalg.inv(R) @ pg_vec
    return L_fam[0,0]

@njit(parallel=True)
def pgs_corr_likelihood(r,pgs_fam,is_sib,is_parent):
    L = 0.0
    for i in prange(pgs_fam.shape[0]):
        L += pgs_corr_likelihood_fam(r, pgs_fam[i,:], is_sib[i,:], is_parent[i,:])
    return L

def pgs_cor_lik(r, *args):
    pgs_fam, is_sib, is_parent = args
    return pgs_corr_likelihood(r[0], pgs_fam, is_sib, is_parent)

def simulate_r_inf(r, nfam, nsib, npar):
    is_sib = np.zeros((nfam,nsib+npar),dtype=np.bool_)
    is_sib[:,0:nsib] = True
    is_parent = np.zeros((nfam,nsib+npar),dtype=np.bool_)
    is_parent[:,nsib:(nsib+npar)] = True
    R = pgs_corr_matrix(r,is_sib[0,:],is_parent[0,:])
    pgs_fam = np.random.multivariate_normal(np.zeros((nsib+npar)), R, size=nfam)
    # correlation between parents
    if npar==2:
        bpg = np.ones((nfam),dtype=bool)
    else:
        bpg = np.zeros((nfam),dtype=bool)
    n_bpg = np.sum(bpg)
    if n_bpg>0:
        r_bpg = np.corrcoef(pgs_fam[:,nsib],pgs_fam[:,nsib+1])[0,1]
    else:
        r_bpg = 0
    # Correlation from sibs
    if nsib>1:
        sib_2 = np.ones((nfam),dtype=bool)
    else:
        sib_2 = np.zeros((nfam),dtype=bool)
    n_sib_2 = np.sum(sib_2)
    if n_sib_2>0:
        r_sib = 2*np.corrcoef(pgs_fam[sib_2,0],pgs_fam[sib_2,1])[0,1]-1
    else:
        r_sib = 0
    # Initialize at weighted average
    r_init = n_bpg*r_bpg+n_sib_2*r_sib
    r_init = r_init/(n_bpg+n_sib_2)
    # Find MLE
    optimized = fmin_l_bfgs_b(func=pgs_cor_lik, x0=r_init,
                            args=(pgs_fam, is_sib, is_parent),
                            approx_grad=True,
                            bounds=[(-0.999,0.999)])
    return np.array([r_init,optimized[0][0]])