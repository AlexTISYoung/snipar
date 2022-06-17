from snipar.gtarray import gtarray
import numpy as np
from snipar.read import get_gts_matrix
from snipar.utilities import *
from scipy.optimize import fmin_l_bfgs_b
from numba import njit, prange
import numpy.ma as ma

class pgs(object):
    """Define a polygenic score based on a set of SNPs with weights and ref/alt allele pairs.

    Args:
        snp_ids : :class:`~numpy:numpy.array`
            [L] vector of SNP ids
        weights : :class:`~numpy:numpy.array`
            [L] vector of weights of equal length to snp_ids
        alleles : :class:`~numpy:numpy.array`
            [L x 2] matrix of ref and alt alleles for the SNPs. L must match size of snp_ids

    Returns:
        pgs : :class:`snipar.pgs`

    """
    def __init__(self,snp_ids,weights,alleles):
        if snp_ids.shape[0] == weights.shape[0] and alleles.shape[0] == weights.shape[0] and alleles.shape[1]==2:
            self.snp_ids = snp_ids
            self.snp_dict = make_id_dict(snp_ids)
            self.weights = weights
            self.alleles = alleles
        else:
            raise ValueError('All inputs must have the same dimension')
    
    def remove_zeros(self):
        zero_weight = self.weights==0
        n_zero = np.sum(zero_weight)
        if n_zero>0:
            print('Filtering '+str(n_zero)+' SNPs with zero weight')
            self.snp_ids = self.snp_ids[~zero_weight]
            self.snp_dict = make_id_dict(self.snp_ids)
            self.weights = self.weights[~zero_weight]
            self.alleles = self.alleles[~zero_weight]
        else:
            print('No zero weight SNPs found')

    def compute(self, garray, cols=None):
        """Compute polygenic score values from a given genotype array. Finds the SNPs in the genotype array
        that have weights in the pgs and matching alleles, and computes the PGS based on these SNPs and the
        weights after allele-matching.


        Args:
            garray : :class:`sbreg.gtarray`
                genotype array to compute PGS values for
            cols : :class:`numpy:numpy.array`
                names to give the columns in the output gtarray

        Returns:
            pg : :class:`snipar.gtarray`
                2d gtarray with PGS values. If a 3d gtarray is input, then each column corresponds to
                the second dimension on the input gtarray (for example, individual, paternal, maternal PGS).
                If a 2d gtarray is input, then there will be only one column in the output gtarray. The
                names given in 'cols' are stored in 'sid' attribute of the output.

        """
        if not type(garray) == gtarray:
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
            pgs_val = ma.dot(garray.gts[:, in_pgs_snps],weights_compute)
        elif garray.ndim == 3:
            pgs_val = np.zeros((garray.gts.shape[0], garray.gts.shape[1]), garray.dtype)
            for i in range(0, garray.gts.shape[1]):
                pgs_val[:, i] = ma.dot(garray.gts[:, i, in_pgs_snps], weights_compute)

        return pgarray(pgs_val, garray.ids, sid=cols, fams=garray.fams, par_status=garray.par_status, ped=garray.ped)
        
def read_weights(weights, SNP='SNP', beta_col='b', A1='A1', A2='A2', sep=None):
    if sep is None:
        weights = np.loadtxt(weights,dtype=str)
    else:
        weights = np.loadtxt(weights,dtype=str, delimiter=sep)
    colnames = weights[0,:]
    weights = weights[1:weights.shape[0],:]
    print('Read weights for '+str(weights.shape[0])+' variants')
    beta = np.array(weights[:,np.where(colnames == beta_col)[0][0]],dtype=np.float64)
    allele_indices = np.array([np.where(colnames==A1)[0][0],np.where(colnames==A2)[0][0]])
    return pgs(weights[:,np.where(colnames==SNP)[0][0]],
            beta,
            weights[:,allele_indices])

def read_pgs(pgs_file):
    pgs_f = open(pgs_file,'r')
    pgs_header = pgs_f.readline().split(' ')
    pgs_header[len(pgs_header)-1] = pgs_header[len(pgs_header)-1].split('\n')[0]
    ncols = len(pgs_header)
    if pgs_header[0]=='FID' and pgs_header[1]=='IID':
        precols = 2
    else:
        raise(ValueError('First two columns of PGS file must be FID and IID'))
    if pgs_header[2]=='FATHER_ID' and pgs_header[3]=='MOTHER_ID':
        precols = 4
    pgs_cols = tuple([x for x in range(precols,ncols)])
    # Read pedigree
    if precols==4:
        ped = np.loadtxt(pgs_file,usecols=tuple([x for x in range(4)]),dtype=str, skiprows=1)
    else:
        ped = None
    if ped is not None:
        par_status = np.array(ped=='NA',dtype=int)
    else:
        par_status = None
    pg = pgarray(np.loadtxt(pgs_file,usecols = pgs_cols, skiprows=1),
                    np.loadtxt(pgs_file,usecols = 1, dtype=str, skiprows=1),
                    sid=np.array(pgs_header[precols:ncols]),
                    fams=np.loadtxt(pgs_file,usecols = 0, dtype=str, skiprows=1),
                    par_status = par_status,
                    ped=ped)
    return pg


def compute(pgs, bedfile=None, bgenfile=None, par_gts_f=None, ped=None, sib=False, compute_controls=False, verbose=True):
    """Compute a polygenic score (PGS) for the individuals with observed genotypes and observed/imputed parental genotypes.

    Args:
        par_gts_f : :class:`str`
            path to HDF5 file with imputed parental genotypes
        gts_f : :class:`str`
            path to bed file with observed genotypes
        pgs : :class:`snipar.pgs`
            the PGS, defined by the weights for a set of SNPs and the alleles of those SNPs
        sib : :class:`bool`
            Compute the PGS for genotyped individuals with at least one genotyped sibling and observed/imputed parental genotypes. Default False.
        compute_controls : :class:`bool`
            Compute polygenic scores for control families (families with observed parental genotypes set to missing). Default False.

    Returns:
        pg : :class:`snipar.gtarray`
            Return the polygenic score as a genotype array with columns: individual's PGS, mean of their siblings' PGS, observed/imputed paternal PGS,
            observed/imputed maternal PGS

    """
    G = get_gts_matrix(bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, ped=ped, snp_ids=pgs.snp_ids, sib=sib, compute_controls=compute_controls, verbose=verbose)
    if sib:
        cols = np.array(['proband', 'sibling', 'paternal', 'maternal'])
    else:
        cols = np.array(['proband', 'paternal', 'maternal'])
    if compute_controls:
        pgs_out = [pgs.compute(x,cols) for x in G[0:3]]
        if sib:
            o_cols = np.array(['proband', 'sibling', 'parental'])
        else:
            o_cols = np.array(['proband','parental'])
        pgs_out.append(pgs.compute(G[3], o_cols))
        return pgs_out
    else:
        return pgs.compute(G,cols)

#def write(pg,filename,scale_PGS = False):
#    if scale_PGS:
#        # Rescale by observed proband PGS
#        pg.gts = pg.gts / np.std(pg.gts[:, 0])
#    ####### Write PGS to file ########
#    pg_out = np.column_stack((pg.fams,pg.ids,pg.gts))
#    pg_header = np.column_stack((np.array(['FID','IID']).reshape(1,2),pg.sid.reshape(1,pg.sid.shape[0])))
#    pg_out = np.row_stack((pg_header,pg_out))
#    print('Writing PGS to ' + filename)
#    np.savetxt(filename, pg_out, fmt='%s')
#    return None

class pgarray(gtarray):

    def add(self,garray):
        """
        Adds another gtarray of the same dimension to this array and returns the sum. It matches IDs before summing.
        """
        if type(garray)==pgarray:
            pass
        else:
            raise ValueError('Must add to another pgarray')

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
        if self.ped is None:
            ped = garray.ped
        else:
            ped = self.ped

        return pgarray(add_gts, ids_out, self.sid, alleles=self.alleles, fams=self.fams[self_index], par_status=self.par_status[self_index,:], ped=ped)

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
                if npar[i]==0:
                    self.gts[i,[paternal_index, maternal_index]] = npg_am_adj(r,fsizes[self.fams[i]])*self.gts[i,[paternal_index, maternal_index]]
                # One parent genotyped
                if npar[i]==1:
                    # Father imputed
                    if self.par_status[i,0] == 1:
                        self.gts[i,paternal_index] = opg_am_adj(self.gts[i,paternal_index],self.gts[i,maternal_index],r,fsizes[self.fams[i]])
                    # Mother imputed
                    if self.par_status[i,1] == 1:
                        self.gts[i,maternal_index] = opg_am_adj(self.gts[i,maternal_index],self.gts[i,paternal_index],r,fsizes[self.fams[i]])
        else:
            print('Estimated correlation is negative, so not performing assortative mating adjustment')
        return r

    def scale(self):
        # Rescale by observed proband PGS
        proband_sd = np.std(self.gts[:, 0])
        self.gts = self.gts / proband_sd
        return proband_sd

    def write(self, filename, scale=False):
        if scale:
            self.scale()
        ####### Write PGS to file ########  
        parent_genotyped = self.par_status == 0
        ped_dict = make_id_dict(self.ped,1)
        ped_indices = np.array([ped_dict[x] for x in self.ids])
        parent_ids = self.ped[ped_indices,2:4]
        parent_ids[~parent_genotyped] = 'NA'
        pg_out = np.column_stack((self.fams,self.ids,parent_ids,self.gts))
        pg_header = np.column_stack((np.array(['FID','IID','FATHER_ID','MOTHER_ID']).reshape(1,4),self.sid.reshape(1,self.sid.shape[0])))
        pg_out = np.row_stack((pg_header,pg_out))
        print('Writing PGS to ' + filename)
        np.savetxt(filename, pg_out, fmt='%s')
        return pg_out
    
    def compute_grandpar(self, r):
        # Check we have required information
        if self.ped is None:
            raise(ValueError('Pedigree needed for grandparental model'))
        if self.par_status is None:
            raise(ValueError('Parental genotype status needed to compute grandparents'))
        # Check pgs columns
        if 'paternal' in self.sid:
            paternal_index = np.where(self.sid=='paternal')[0][0]
        else:
            raise(ValueError('No paternal PGS column found'))
        if 'maternal' in self.sid:
            maternal_index = np.where(self.sid=='maternal')[0][0]
        else:
            raise(ValueError('No maternal PGS column found'))
        # Find individuals with both parents genotyped
        bpg = np.sum(self.par_status==0,axis=1)==2
        n_bpg = np.sum(bpg)
        print('Found '+str(n_bpg)+' individuals with both parents genotyped for which grandparental scores will be computed')
        ## Create grandparental PGS matrix
        gpar = np.zeros((self.gts.shape[0],4))
        gpar[:] = np.nan
        if n_bpg>0:
            # Find their parents' IDs
            ped_dict = make_id_dict(self.ped,1)
            bpg_ped = self.ped[[ped_dict[x] for x in self.ids[bpg]],:]
            ## Fill
            for i in range(bpg_ped.shape[0]):
                i_index = self.id_dict[bpg_ped[i,1]]
                ## Paternal grandparental scores
                # Find imputed grandparents
                if bpg_ped[i,2] in self.id_dict:
                    gpar[i_index, 0:2] = self.gts[self.id_dict[bpg_ped[i,2]],[paternal_index,maternal_index]]
                else:
                    # linear imputation from father if no imputed
                    gpar[i_index, 0:2] = (1+r)*self.gts[i_index, paternal_index]/2.0
                ## Maternal grandparental scores
                # Find imputed grandparents
                if bpg_ped[i,3] in self.id_dict:
                    gpar[i_index, 2:4] = self.gts[self.id_dict[bpg_ped[i,3]],[paternal_index,maternal_index]]
                else:
                    # linear imputation from mother if no imputed
                    gpar[i_index, 2:4] = (1+r)*self.gts[i_index, maternal_index]/2.0
        ## Append to PGS matrix
        self.gts = np.hstack((self.gts,gpar))
        self.sid = np.hstack((self.sid,np.array(['gpp','gpm','gmp','gmm'])))
        return bpg_ped

        
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