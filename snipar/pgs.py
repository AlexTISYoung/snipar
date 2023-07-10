from snipar.gtarray import gtarray
import numpy as np
from snipar.read import get_gts_matrix
from snipar.utilities import *
from scipy.optimize import fmin_l_bfgs_b
from numba import njit, prange
import numpy.ma as ma
import snipar.slmm as slmm
from pysnptools.snpreader import Bed
import numdifftools as nd
from snipar.pgs_am import *

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
            print('No overlap between PGS SNPs and genotype SNPs')
            return None
        else:
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
            geno_means = np.mean(garray.gts[:, in_pgs_snps],axis=0)
            garray.gts = garray.gts[:, in_pgs_snps]-geno_means
            pgs_val = ma.dot(garray.gts[:, in_pgs_snps],weights_compute)
        elif garray.ndim == 3:
            geno_means = np.mean(garray.gts[:, 0, in_pgs_snps], axis=0)
            pgs_val = np.zeros((garray.gts.shape[0], garray.gts.shape[1]), garray.dtype)
            for i in range(0, garray.gts.shape[1]):
                if cols[i]=='parental':
                    sf = 2.0                
                else:
                    sf = 1.0
                garray.gts[:, i, in_pgs_snps] = garray.gts[:, i, in_pgs_snps]-sf*geno_means
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
        par_status = np.array(ped[:,2:4]=='NA',dtype=int)
    else:
        par_status = None
    pg = pgarray(np.genfromtxt(pgs_file,usecols = pgs_cols, skip_header=1, missing_values=['nan','NA','N/A']),
                    np.loadtxt(pgs_file,usecols = 1, dtype=str, skiprows=1),
                    sid=np.array(pgs_header[precols:ncols]),
                    fams=np.loadtxt(pgs_file,usecols = 0, dtype=str, skiprows=1),
                    par_status = par_status,
                    ped=ped)
    return pg


def compute(pgs, bedfile=None, bgenfile=None, par_gts_f=None, ped=None, sib=False, compute_controls=False, verbose=True, batch_size=None):
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
    # Check for SNP overlap
    if bedfile is not None:
        bed = Bed(bedfile, count_A1=True)
        snp_ids = bed.sid
    if bgenfile is not None:
        bgen = open_bgen(bgenfile)
        snp_ids = bgen.ids
        if np.unique(snp_ids).shape[0] == 1:
            snp_ids = bgen.rsids
    pgs_snp_set = set(pgs.snp_ids)
    in_snp_set = np.array([x in pgs_snp_set for x in snp_ids])
    if np.sum(in_snp_set)==0:
        print('No overlap between variants in weights file and observed genotypes')
        return None
    else:
        ## Get genotype matrix
        snps_in_pgs = snp_ids[in_snp_set]
        # Get batch size
        if batch_size is not None:
            batch_size = min(batch_size, snps_in_pgs.shape[0])
        else: 
            batch_size = snps_in_pgs.shape[0]
        # Set batch boundaries 
        n_batches = int(np.ceil(snps_in_pgs.shape[0]/batch_size))
        batch_boundaries = np.zeros((n_batches,2),dtype=int)
        for i in range(n_batches-1):
            batch_boundaries[i,:] = [i*batch_size,(i+1)*batch_size]
        batch_boundaries[n_batches-1,:] = [(n_batches-1)*batch_size,snps_in_pgs.shape[0]]
        ## Compute PGS, reading snps in batches
        # First batch
        G = get_gts_matrix(bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, ped=ped, 
                           snp_ids=snps_in_pgs[batch_boundaries[0,0]:batch_boundaries[0,1]], 
                           sib=sib, compute_controls=compute_controls, verbose=verbose)
        if sib:
            if G.shape[1]==4:
                cols = np.array(['proband', 'sibling', 'paternal', 'maternal'])
            else:
                cols = np.array(['proband','sibling','parental'])
        else:
            if G.shape[1]==3:
                cols = np.array(['proband', 'paternal', 'maternal'])
            else:
                cols = np.array(['proband','parental'])
        if compute_controls:
            pgs_out = [pgs.compute(x,cols) for x in G[0:3]]
            if sib:
                o_cols = np.array(['proband', 'sibling', 'parental'])
            else:
                o_cols = np.array(['proband','parental'])
            pgs_out.append(pgs.compute(G[3], o_cols))
        else:
            pgs_out = pgs.compute(G,cols)
        # Remaining batches
        for i in range(1,n_batches):
            del G
            G = get_gts_matrix(bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, ped=ped, 
                               snp_ids=snps_in_pgs[batch_boundaries[i,0]:batch_boundaries[i,1]], 
                               sib=sib, compute_controls=compute_controls, verbose=False)
            if compute_controls:
                pgs_out_i = [pgs.compute(x,cols) for x in G[0:3]]
                if sib:
                    o_cols = np.array(['proband', 'sibling', 'parental'])
                else:
                    o_cols = np.array(['proband','parental'])
                pgs_out_i.append(pgs.compute(G[3], o_cols))
                pgs_out = [pgs_out[x].add(pgs_out_i[x]) for x in range(0, len(pgs_out))]
            else:
                pgs_out = pgs_out.add(pgs.compute(G,cols))
        return pgs_out
        

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

    def filter_bpg(self):
        if self.par_status is None:
            raise(ValueError('Parental genotype status unknown so cannot restrict to both parents genotyped sample'))
        else:
            bpg_ids = self.ids[np.sum(self.par_status==0, axis=1)==2]
            self.filter_ids(bpg_ids)

    def estimate_r(self, return_se=True, parents_only=False):
        # Check pgs columns
        if 'paternal' in self.sid:
            paternal_index = np.where(self.sid=='paternal')[0][0]
        if 'maternal' in self.sid:
            maternal_index = np.where(self.sid=='maternal')[0][0]
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
        if np.sum(is_mother)>0:
            pgs_fam[is_mother] = (pgs_fam[is_mother]-np.mean(pgs_fam[is_mother]))/np.std(pgs_fam[is_mother])
        if np.sum(is_father)>0:
            pgs_fam[is_father] = (pgs_fam[is_father]-np.mean(pgs_fam[is_father]))/np.std(pgs_fam[is_father])
        # fam size dict
        fsizes = dict(zip(list(fams[0]),list(fams[1])))
        ### find correlation between maternal and paternal pgis
        print('Finding MLE for correlation between parents scores')
        ## Initialize with correlation from sibs and between parents
        # correlation between parents
        bpg = np.sum(self.par_status==0,axis=1)==2
        n_bpg = np.sum(bpg)
        if n_bpg>0:
            r_bpg = np.corrcoef(self.gts[bpg,1],self.gts[bpg,2])[0,1]
        else:
            r_bpg = 0
        if parents_only:
            if n_bpg>0:
                r_se = (1-r_bpg**2)/np.sqrt(n_bpg-1)
                return r_bpg, r_se, fsizes
            else:
                return np.nan, np.nan, fsizes
        else:
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
            if optimized[2]['warnflag']==0:
                r=optimized[0][0]
                hess = nd.Hessian(pgs_cor_lik)
                hess = float(hess([r], pgs_fam, is_sib, is_parent))
                r_se = np.sqrt(2/hess)
            else:
                print('Could not find MLE for correlation. Returning weighted average from siblings and parents.')
                r=r_init
                r_se = (1-r**2)/np.sqrt(n_bpg+n_sib_2-1)
            return r, r_se, fsizes
    
    def am_adj(self):
        print('Estimating correlation between maternal and paternal PGSs assuming equilibrium')
        r, r_se, fsizes = self.estimate_r()
        r_z = r/r_se
        print('Estimated correlation between maternal and paternal PGSs: '+str(round(r,4))+' S.E.='+str(round(r_se,4)))
        if r_z>1.5:    
            # Check pgs columns
            if 'paternal' in self.sid and 'maternal' in self.sid:
                paternal_index = np.where(self.sid=='paternal')[0][0]
                maternal_index = np.where(self.sid=='maternal')[0][0]
                parental_index = None
            elif 'parental' in self.sid:
                parental_index = np.where(self.sid=='parental')[0][0]
            else:
                raise(ValueError('No parental PGS values to adjust'))
            # Adjust imputed parental PGSs
            npar = np.sum(self.par_status==0,axis=1)
            print('Adjuting imputed PGSs for assortative mating')
            for i in range(self.gts.shape[0]):
                # No parents genotyped
                if npar[i]==0:
                    if parental_index is None:
                        self.gts[i,[paternal_index, maternal_index]] = npg_am_adj(r,fsizes[self.fams[i]])*self.gts[i,[paternal_index, maternal_index]]
                    else:
                        self.gts[i,parental_index] = npg_am_adj(r,fsizes[self.fams[i]])*self.gts[i,parental_index] 
                # One parent genotyped
                if npar[i]==1:
                    # Father imputed
                    if self.par_status[i,0] == 1:
                        self.gts[i,paternal_index] = opg_am_adj(self.gts[i,paternal_index],self.gts[i,maternal_index],r,fsizes[self.fams[i]])
                    # Mother imputed
                    if self.par_status[i,1] == 1:
                        self.gts[i,maternal_index] = opg_am_adj(self.gts[i,maternal_index],self.gts[i,paternal_index],r,fsizes[self.fams[i]])
        else:
            print('Signal to noise ratio for correlation estimate too small, so not performing assortative mating adjustment. Use --force_am_adj to override')
        return r, r_se

    def scale(self, sf=None):
        if sf is None:
            if 'proband' in self.sid:
                proband_index = np.where(self.sid=='proband')[0][0]
            else:
                raise(ValueError('Cannot scale as no proband column found and no scale factor provided'))
            sf = np.std(self.gts[:, proband_index])
        # Rescale by observed proband PGS
        self.gts = self.gts / sf
        return sf

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
            # Find the mean of the mothers and fathers
            parent_means = np.mean(self.gts[bpg,:],axis=0)[[paternal_index, maternal_index]]
            ## Fill
            for i in range(bpg_ped.shape[0]):
                i_index = self.id_dict[bpg_ped[i,1]]
                ## Paternal grandparental scores
                # Find imputed grandparents
                if bpg_ped[i,2] in self.id_dict:
                    gpar[i_index, 0:2] = self.gts[self.id_dict[bpg_ped[i,2]],[paternal_index,maternal_index]]
                else:
                    # linear imputation from father if no imputed
                    gpar[i_index, 0:2] = parent_means[0]+(1+r)*(self.gts[i_index, paternal_index]-parent_means[0])/2.0
                ## Maternal grandparental scores
                # Find imputed grandparents
                if bpg_ped[i,3] in self.id_dict:
                    gpar[i_index, 2:4] = self.gts[self.id_dict[bpg_ped[i,3]],[paternal_index,maternal_index]]
                else:
                    # linear imputation from mother if no imputed
                    gpar[i_index, 2:4] = parent_means[1]+(1+r)*(self.gts[i_index, maternal_index]-parent_means[1])/2.0
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

def fit_pgs_model(y, pg, ngen, ibdrel_path=None, covariates=None, fit_sib=False, parsum=False, gparsum=False, outprefix=None, sparse_thresh=0.025):
    pg.gts = ma.array(pg.gts,fill_value=np.nan)
    pg.gts = pg.gts.filled()
    if ngen in [1,2,3]:
        pass
    else:
        raise(ValueError('ngen must be 1, 2, or 3'))
    # Check if IDs are aligned
    if not np.array_equal(y.ids,pg.ids):
        pg.filter_ids(y.ids)
        y.filter_ids(pg.ids)
    # Check covariate IDs are aligned
    if covariates is not None:
        if not np.array_equal(pg.ids,covariates.ids):
            covariates.filter_ids(pg.ids)
    ## Fit model
    if ngen==1:
        print('Fitting 1 generation model (proband only)')
        alpha, alpha_cols = make_and_fit_model(y, pg, ['proband'], ibdrel_path=ibdrel_path, covariates=covariates, sparse_thresh=sparse_thresh)
    elif ngen==2 or ngen==3:
        if fit_sib:
            if 'sib' in pg.sid:
                pg_cols = ['proband','sibling']
            else:
                raise(ValueError('Sibling PGS not found (use --fit_sib when calculating PGS)'))
        else:
            pg_cols = ['proband']
        if parsum:
            if 'maternal' in pg.sid and 'paternal' in pg.sid:
                parcols = np.sort(np.array([np.where(pg.sid=='maternal')[0][0],np.where(pg.sid=='paternal')[0][0]]))
                trans_matrix = np.identity(pg.gts.shape[1])
                trans_matrix[:,parcols[0]] += trans_matrix[:,parcols[1]]
                trans_matrix = np.delete(trans_matrix,parcols[1],1)
                pg.gts = pg.gts.dot(trans_matrix)
                pg.sid = np.delete(pg.sid,parcols[1])
                pg.sid[parcols[0]] = 'parental'
            elif 'parental' in pg.sid:
                pass
            else:
                raise(ValueError('Maternal and paternal PGS not found so cannot sum (--parsum option given)'))
            pg_cols.append('parental')
        else:
            pg_cols += ['paternal','maternal']
        if ngen==2:
            print('Fitting 2 generation model (proband and observed/imputed parents)')
            alpha, alpha_cols = make_and_fit_model(y, pg, pg_cols, ibdrel_path=ibdrel_path, covariates=covariates, sparse_thresh=sparse_thresh)
        elif ngen==3:
            print('Fitting 3 generation model: observed proband and observed parents, and observed/imputed grandparents')
            if gparsum:
                if 'gpp' in pg.sid and 'gpm' in pg.sid and 'gmp' in pg.sid and 'gmm' in pg.sid:
                    # Sum of paternal grandparents
                    gparcols = np.sort(np.array([np.where(pg.sid==x)[0][0] for x in ['gpp','gpm']]))
                    trans_matrix = np.identity(pg.gts.shape[1])
                    trans_matrix[:,gparcols[0]] += trans_matrix[:,gparcols[1]]
                    trans_matrix = np.delete(trans_matrix,gparcols[1],1)
                    pg.gts = pg.gts.dot(trans_matrix)
                    pg.sid = np.delete(pg.sid,gparcols[1])
                    pg.sid[gparcols[0]] = 'gp'
                    # Sum of maternal grandparents
                    gparcols = np.sort(np.array([np.where(pg.sid==x)[0][0] for x in ['gmp','gmm']]))
                    trans_matrix = np.identity(pg.gts.shape[1])
                    trans_matrix[:,gparcols[0]] += trans_matrix[:,gparcols[1]]
                    trans_matrix = np.delete(trans_matrix,gparcols[1],1)
                    pg.gts = pg.gts.dot(trans_matrix)
                    pg.sid = np.delete(pg.sid,gparcols[1])
                    pg.sid[gparcols[0]] = 'gm'
                elif 'gp' in pg.sid and 'gm' in pg.sid:
                    pass
                else:
                    raise(ValueError('Grandparental PGSs not found so cannot sum (--gparsum option given)'))
                pg_cols += ['gp','gm']
            else:
                pg_cols += ['gpp','gpm','gmp','gmm']
            alpha, alpha_cols = make_and_fit_model(y, pg, pg_cols, ibdrel_path=ibdrel_path, covariates=covariates, sparse_thresh=sparse_thresh)
    # Save to file
    if outprefix is not None:
        write_estimates(outprefix+'.'+str(ngen), alpha, alpha_cols)
    return alpha, alpha_cols

def make_and_fit_model(y, pg, pg_cols, ibdrel_path=None, covariates=None, sparse_thresh=0.025):
    pg_col_indices = [np.where(pg.sid==x)[0][0] for x in pg_cols]
    if covariates is not None:
        X = np.hstack((covariates.gts,np.array(pg.gts[:,pg_col_indices])))
        X_cols = np.hstack((np.array(['intercept']),covariates.sid,pg_cols))
    else:
        X = np.array(pg.gts[:, pg_col_indices])
        X_cols = np.hstack((np.array(['intercept']),pg_cols))
    # Check for NAs
    no_NA = np.sum(np.isnan(X),axis=1)==0
    print('Sample size: '+str(np.sum(no_NA)))
    ## Read GRM
    id_dict = make_id_dict(pg.ids[no_NA])
    varcomp_lst = make_grms(pg.fams[no_NA], ibdrel_path=ibdrel_path, id_dict=id_dict, keep=pg.ids[no_NA], sparse_thresh=sparse_thresh)
    ## Define LMM
    slmm_model = slmm.LinearMixedModel(np.array(y.gts[no_NA,0],dtype=float), 
                    varcomp_arr_lst=varcomp_lst, covar_X=np.array(X[no_NA,:],dtype=float), add_intercept=True)
    # Optimize Lmm
    slmm_model.scipy_optimize()
    # Print variance components
    print(f'Variance components: {list(i for i in slmm_model.varcomps)}')
    ZT_Vinv_Z_imp = slmm_model.Z.T @ slmm_model.Vinv_Z
    alpha = [np.linalg.solve(ZT_Vinv_Z_imp, slmm_model.Z.T @ slmm_model.Vinv_y), np.linalg.inv(ZT_Vinv_Z_imp)]
    return alpha, X_cols

def write_estimates(outprefix, alpha, cols):
    cols = cols.reshape((cols.shape[0],1))
    # Find cols
    alpha_out = np.zeros((alpha[0].shape[0], 2))
    alpha_out[:, 0] = alpha[0]
    alpha_out[:, 1] = np.sqrt(np.diag(alpha[1]))
    print('Saving effect estimates to '+outprefix+ '.effects.txt')
    np.savetxt(outprefix + '.effects.txt',
                np.hstack((cols, np.array(alpha_out, dtype='S'))),
                delimiter='\t', fmt='%s')
    vcols = np.zeros((1,cols.shape[0]+1),dtype=cols.dtype)
    vcols[0,1:vcols.shape[1]] = cols[:,0]
    alpha_cov_out = np.vstack((vcols, np.hstack((cols,alpha[1]))))
    print('Saving sampling variance-covariance matrix to '+outprefix+ '.vcov.txt')
    np.savetxt(outprefix+ '.vcov.txt',alpha_cov_out, fmt='%s')

def make_grms(fams, ibdrel_path=None, id_dict=None, keep=None, sparse_thresh=0.05):
    if ibdrel_path is not None:
        grm_data, grm_row_ind, grm_col_ind = slmm.build_ibdrel_arr(
            ibdrel_path, id_dict=id_dict, keep=keep, thres=sparse_thresh)
    ## Build sparse sib GRM
    sib_data, sib_row_ind, sib_col_ind = slmm.build_sib_arr(fams)
    ## GRM list
    if 'grm_data' in locals():
        varcomp_lst = (
            (grm_data, grm_row_ind, grm_col_ind),
            (sib_data, sib_row_ind, sib_col_ind),
        )
    else:
        varcomp_lst = (
            (sib_data, sib_row_ind, sib_col_ind),
        )
    return varcomp_lst

######## Assortative mating adjustment in two-generation model ###########

def h2f_parse(h2f_str):
    h2f_split = h2f_str.split(',')
    if not len(h2f_split)==2:
        raise(ValueError('Invalid h2f input. Please use estimate,se format, such as 0.5,0.01.'))
    h2f = float(h2f_split[0])
    h2f_se = float(h2f_split[1].split('\n')[0])
    return h2f, h2f_se

def am_adj_2gen_calc(delta, delta_se, ab, ab_se, r_delta_ab, h2f, h2f_se, rk, rk_se, is_beta=False, verbose=True):
    # Get estimates
    if is_beta:
        estimates = v_eta_delta(delta, delta_se, h2f, h2f_se, rk, beta=ab)
        estimates['beta'] = ab
    else:
        estimates = v_eta_delta(delta, delta_se, h2f, h2f_se, rk, alpha=ab)
        estimates['alpha'] = ab
    estimates['rk'] = rk
    estimates['delta'] = delta
    
    # Get ses
    ses = {'rk':rk_se,
           'k':k_se(delta, delta_se, h2f, h2f_se, rk, rk_se),
           'r':r_se(delta, delta_se, h2f, h2f_se, rk, rk_se),
           'h2_eq':h2eq_se(delta, delta_se, h2f, h2f_se, rk, rk_se),
           'rho':rho_se(delta, delta_se, h2f, h2f_se, rk, rk_se),
           'delta':delta_se,
           'v_eta_delta':se_v_eta_delta(delta, delta_se, ab, ab_se, r_delta_ab, h2f, h2f_se, rk, rk_se, is_beta=is_beta)}
    if is_beta:
        ses['alpha_delta'] = se_alpha_from_beta(delta, delta_se, ab, ab_se, r_delta_ab, h2f, h2f_se, rk, rk_se)
        ses['beta'] = ab_se
        ses['r_delta_beta'] = r_delta_ab
    else:
        ses['alpha_delta'] = se_alpha_from_alpha(delta, delta_se, ab, ab_se, r_delta_ab, h2f, h2f_se, rk, rk_se)
        ses['alpha'] = ab_se
        ses['r_delta_alpha'] = r_delta_ab
    # Print
    if verbose:
        for par in ['rk','k','r','h2_eq','rho','alpha_delta','v_eta_delta']:
            print(par+': '+str(round(estimates[par],4))+' ('+str(round(ses[par],4))+')')
    # Return
    return estimates, ses

def write_2gen_adj_ests(estimates,ses, outprefix=''):
    pars = np.array(['delta','alpha','rk','k','r','h2_eq','rho','alpha_delta','v_eta_delta'])
    outarray = np.zeros((pars.shape[0],2))
    parcount = 0
    for par in pars:
        outarray[parcount,:] = [estimates[par],ses[par]]
        parcount += 1
    outarray = np.hstack((pars.reshape((pars.shape[0],1)),outarray))
    outarray = np.vstack((np.array(['parameter','estimate','SE']).reshape((1,3)),outarray))
    np.savetxt(outprefix+'.am_adj_pars.txt',outarray,fmt='%s')
    return outarray

def am_adj_2gen(estimates, estimate_cols, h2f, h2f_se, rk=None, rk_se=None, pg=None, y_std=1, pg_std=1):
    """
    Adjust 2-generation model results for assortative mating (assuming equilibrium)
    """
    ##### Check we have r, or a pg object to estimate r with #####
    if pg is None:
        if rk is None:
            raise(ValueError('Need estimate of correlation between parents PGIs if not providing pg object to use to estimate'))
        if rk_se is None:
            raise(ValueError('Need standard error of correlation between parents PGIs if not providing pg object to use to estimate')) 
    ##### Find estimates and their sampling variance-covariance #####
    if 'proband' in estimate_cols:
        proband_index = np.where(estimate_cols=='proband')[0][0]
    else:
        raise(ValueError('Proband coefficient (direct effect) estimate needed for adjustment'))
    # If separate maternal and paternal NTC estimates, take average
    if 'paternal' in estimate_cols and 'maternal' in estimate_cols:
        paternal_index = np.where(estimate_cols=='paternal')[0][0]
        maternal_index = np.where(estimate_cols=='maternal')[0][0]
        est_cols = [proband_index, paternal_index, maternal_index]
        estimate = estimates[0][est_cols].reshape((3,1))
        estimate_cov = estimates[1][np.ix_(est_cols,est_cols)]
        A =  np.array([[1,0,0],[0,0.5,0.5]])
        estimate = A @ estimate
        estimate_cov = A @ estimate_cov @ A.T
    elif 'parental' in estimate_cols:
        est_cols = [proband_index, np.where(estimate_cols=='parental')[0][0]]
        estimate = estimates[0][est_cols].reshape((2,1))
        estimate_cov = estimates[1][np.ix_(est_cols,est_cols)]
    else:
        raise(ValueError('Need parental NTC estimate(s) for adjustment')) 
    # Adjust for non-normalized pgs/phenotype
    estimate = estimate*pg_std/y_std
    estimate_cov = estimate_cov*((pg_std/y_std)**2)
    ##### Estimate correlation between parents' PGIs ######
    if rk is None:
        rk, rk_se, fam_sizes = pg.estimate_r(parents_only=True)
    ##### Estimate equilibrium quantities ######
    print('Computing equilibrium adjusted quantities')
    adj_estimates, adj_ses = am_adj_2gen_calc(estimate[0,0], np.sqrt(estimate_cov[0,0]), 
                                              estimate[1,0], np.sqrt(estimate_cov[0,0]), 
                                              estimate_cov[0,1]/np.sqrt(estimate_cov[0,0]*estimate_cov[1,1]), 
                                              h2f, h2f_se, rk, rk_se, verbose=True)
    return adj_estimates, adj_ses