from snipar.gtarray import gtarray
import numpy as np
from snipar.read import get_gts_matrix
from snipar.utilities import *

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
            pgs_val = garray.gts[:, in_pgs_snps].dot(weights_compute)
        elif garray.ndim == 3:
            pgs_val = np.zeros((garray.gts.shape[0], garray.gts.shape[1]), garray.dtype)
            for i in range(0, garray.gts.shape[1]):
                pgs_val[:, i] = garray.gts[:, i, in_pgs_snps].dot(weights_compute)

        return gtarray(pgs_val, garray.ids, sid=cols, fams=garray.fams, par_status=garray.par_status)
    
    def am_adj(self):
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
        is_parent = np.zeros((fams[0].shape[0],np.max(fsizes)),dtype=bool)
        # populate
        for i in range(fams[0].shape[0]):
            # Fill in sibs
            sib_indices = fam_dict[fams[0][i]]
            pgs_fam[i,0:sib_indices.shape[0]] = self.gts[sib_indices,0]
            is_sib[i,0:sib_indices.shape[0]] = True
            npar = 0
            if par_status_fams[i,0] == 0:
                pgs_fam[i,sib_indices.shape[0]] = self.gts[sib_indices[0],paternal_index]
                is_parent[i,sib_indices.shape[0]] = True
                npar = 1
            if par_status_fams[i,1] == 0:
                pgs_fam[i,sib_indices.shape[0]+npar] = self.gts[sib_indices[0],maternal_index]
                is_parent[i,sib_indices.shape[0]+npar] = True
        # normalize
        pgs_fam = (pgs_fam-np.mean(pgs_fam[~np.isnan(pgs_fam)]))/np.std(pgs_fam[~np.isnan(pgs_fam)])
        ### find correlation between maternal and paternal pgis
        # grid search over r
        print('Finding MLE for correlation between parents scores')
        rvals = -0.999+np.arange(1000)*(2*0.999/999)
        L_vals = np.zeros(1000)
        for i in range(1000):
            L_vals[i] = pgs_corr_likelihood(rvals[i],pgs_fam,is_sib,is_parent)
        r_mle = rvals[np.argmin(L_vals)]
        print('r='+str(round(r_mle,3)))

def pgs_corr_matrix(r,is_sib_fam,is_parent_fam):
    n_sib = np.sum(is_sib_fam)
    n_par = np.sum(is_parent_fam)
    # Full matrix
    R = np.identity(n_sib+n_par,dtype=np.float_)
    # sib submatrix
    r_sib = (1+r)/2
    R_sib = (1-r_sib)*np.identity(n_sib,dtype=np.float_)+r_sib*np.ones((n_sib,n_sib),dtype=np.float_)
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

def pgs_corr_likelihood_fam(r,pg_fam,is_sib_fam,is_parent_fam):
    sib_or_parent = np.logical_or(is_sib_fam,is_parent_fam)
    R = pgs_corr_matrix(r,is_sib_fam,is_parent_fam)
    slogdet_R = np.linalg.slogdet(R)
    pg_vec = pg_fam[sib_or_parent].reshape((np.sum(sib_or_parent),1))
    return slogdet_R[0]*slogdet_R[1]+pg_vec.T @ np.linalg.inv(R) @ pg_vec

def pgs_corr_likelihood(r,pgs_fam,is_sib,is_parent):
    L = 0
    for i in range(pgs_fam.shape[0]):
        L += pgs_corr_likelihood_fam(r,pgs_fam[i,:],is_sib[i,:],is_parent[i,:])
    return L


pgs_corr_matrix(0.5,np.array([True,True,False]),np.array([False,False,True]))


        

        

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

def write(pg,filename,scale_PGS = False):
    if scale_PGS:
        # Rescale by observed proband PGS
        pg.gts = pg.gts / np.std(pg.gts[:, 0])
    ####### Write PGS to file ########
    pg_out = np.column_stack((pg.fams,pg.ids,pg.gts))
    pg_header = np.column_stack((np.array(['FID','IID']).reshape(1,2),pg.sid.reshape(1,pg.sid.shape[0])))
    pg_out = np.row_stack((pg_header,pg_out))
    print('Writing PGS to ' + filename)
    np.savetxt(filename, pg_out, fmt='%s')
    return None