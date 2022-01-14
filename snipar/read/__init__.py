import snipar.read.bed as bed
import snipar.read.bgen as bgen
import snipar.read.phenotype as phenotype
import h5py
import numpy as np
from snipar.utilities import convert_str_array

def get_gts_matrix(par_gts_f, gts_f, snp_ids = None,ids = None, sib = False, compute_controls = False, parsum = False, start=0, end=None, print_sample_info=False):
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
        G : :class:`snipar.gtarray`
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
    # Compute genotype matrices
    if gts_f[(len(gts_f)-4):len(gts_f)] == '.bed':
        G = [bed.get_gts_matrix_given_ped(ped[np.logical_not(controls),:],par_gts_f,gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum, start=start, end=end,
                                      print_sample_info = print_sample_info)]
        if compute_controls:
            G.append(bed.get_gts_matrix_given_ped(ped[np.array([x[0:3]=='_p_' for x in ped[:,0]]),],par_gts_f,gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum, start=start, end=end,
                                              print_sample_info = print_sample_info))
            G.append(
                bed.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_m_' for x in ped[:, 0]]),], par_gts_f, gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum, start=start, end=end,
                                         print_sample_info = print_sample_info))
            G.append(
                bed.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_o_' for x in ped[:, 0]]),], par_gts_f, gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum, start=start, end=end,
                                         print_sample_info = print_sample_info))
            return G
        else:
            return G[0]
    elif gts_f[(len(gts_f)-5):len(gts_f)]  == '.bgen':
        G = [bgen.get_gts_matrix_given_ped(ped[np.logical_not(controls),:],par_gts_f,gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum, start=start, end=end,
                                           print_sample_info = print_sample_info)]
        if compute_controls:
            G.append(bgen.get_gts_matrix_given_ped(ped[np.array([x[0:3]=='_p_' for x in ped[:,0]]),],par_gts_f,gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum, start=start, end=end,
                                                   print_sample_info = print_sample_info))
            G.append(
                bgen.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_m_' for x in ped[:, 0]]),], par_gts_f, gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum, start=start, end=end,
                                              print_sample_info = print_sample_info))
            G.append(
                bgen.get_gts_matrix_given_ped(ped[np.array([x[0:3] == '_o_' for x in ped[:, 0]]),], par_gts_f, gts_f,
                                      snp_ids=snp_ids, ids=ids, sib=sib, parsum=parsum, start=start, end=end,
                                              print_sample_info = print_sample_info))
            return G
        else:
            return G[0]
    else:
        raise(ValueError('Unknown filetype for observed genotypes file: '+str(gts_f)))
