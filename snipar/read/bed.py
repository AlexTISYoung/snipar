import snipar.preprocess as preprocess
import numpy as np
from pysnptools.snpreader import Bed
from snipar.gtarray import gtarray
from snipar.utilities import *

def match_observed_and_imputed_snps(gts_f, par_gts_f, bim, snp_ids=None, start=0, end=None):
    """
    Used in get_gts_matrix_given_ped to match observed and imputed SNPs and return SNP information on shared SNPs.
    Removes SNPs that have duplicated SNP ids.
    in_obs_sid contains the SNPs in the imputed genotypes that are present in the observed SNPs
    obs_sid_index contains the index in the observed SNPs of the common SNPs
    """
    # Match SNPs from imputed and observed and restrict to those in list
    if snp_ids is None:
        snp_ids = gts_f.sid
        if end is None:
            end = snp_ids.shape[0]
        snp_ids = snp_ids[start:end]
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

def get_snps(gts_f,bim,snp_ids=None):
    # Match SNPs in bed and get SNP info
    if snp_ids is None:
        snp_ids = gts_f.sid
    # Get bim info
    alleles = np.loadtxt(bim, dtype='U', usecols=(4, 5))
    pos = np.loadtxt(bim, dtype=int, usecols=3)
    chromosome = np.loadtxt(bim, dtype=int, usecols=0)
    # Remove duplicate ids
    unique_snps, snp_indices, snp_counts = np.unique(snp_ids, return_index=True, return_counts=True)
    snp_ids = snp_ids[snp_indices[snp_counts == 1]]
    # Read and match SNP ids
    obs_sid = gts_f.sid
    obs_sid_dict = make_id_dict(obs_sid)
    in_obs_sid = np.array([x in obs_sid_dict for x in snp_ids])
    if np.sum(in_obs_sid) == 0:
        raise ValueError('No SNPs found in bed file')
    obs_sid_index = np.array([obs_sid_dict[x] for x in snp_ids[in_obs_sid]])
    sid = obs_sid[obs_sid_index]
    alleles = alleles[obs_sid_index, :]
    chromosome = chromosome[obs_sid_index]
    pos = pos[obs_sid_index]
    return chromosome, sid, pos, alleles, obs_sid_index


def get_gts_matrix_given_ped(ped, bedfile, par_gts_f=None, snp_ids=None, ids=None, sib=False, parsum=False, verbose=False, print_sample_info = False):
    """
    Used in get_gts_matrix: see get_gts_matrix for documentation
    """
    ### Genotype file ###
    bim = bedfile.split('.bed')[0] + '.bim'
    gts_f = Bed(bedfile,count_A1=True)
    # get ids of genotypes and make dict
    gts_ids = gts_f.iid[:, 1]
    # Get families with imputed parental genotypes
    if par_gts_f is not None:
        imp_fams = convert_str_array(par_gts_f['families'])
    else:
        imp_fams = None
    ### Find ids with observed/imputed parents and indices of those in observed/imputed data
    ids, observed_indices, imp_indices = preprocess.get_indices_given_ped(ped, gts_ids, imp_fams=imp_fams, ids=ids, 
                                                                                sib=sib, verbose=print_sample_info)
    ### Match observed and imputed SNPs ###
    if par_gts_f is not None:
        if verbose:
            print('Matching observed and imputed SNPs')
        chromosome, sid, pos, alleles, in_obs_sid, obs_sid_index = match_observed_and_imputed_snps(gts_f, par_gts_f, bim, snp_ids=snp_ids)
        # Read imputed parental genotypes
        if verbose:
            print('Reading imputed parental genotypes')
        if (imp_indices.shape[0]*in_obs_sid.shape[0]) < (np.sum(in_obs_sid)*imp_fams.shape[0]):
            imp_gts = np.array(par_gts_f['imputed_par_gts'][imp_indices, :])
            imp_gts = imp_gts[:,np.arange(in_obs_sid.shape[0])[in_obs_sid]]
        else:
            imp_gts = np.array(par_gts_f['imputed_par_gts'][:,np.arange(in_obs_sid.shape[0])[in_obs_sid]])
            imp_gts = imp_gts[imp_indices,:]
        imp_fams = imp_fams[imp_indices]
    else:
        chromosome, sid, pos, alleles, obs_sid_index = get_snps(gts_f, bim, snp_ids=snp_ids)
        imp_gts = None
    # Read observed genotypes
    if verbose:
        print('Reading observed genotypes')
    gts = gts_f[observed_indices, obs_sid_index].read().val
    gts_ids = gts_f.iid[observed_indices,1]
    gts_id_dict = make_id_dict(gts_ids)
    # Find indices in reduced data
    par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict, imp_fams=imp_fams)
    if verbose:
        print('Constructing family based genotype matrix')
    ### Make genotype design matrix
    if sib:
        if parsum:
            G = np.zeros((ids.shape[0], 3, gts.shape[1]), dtype=np.float32)
            G[:, np.array([0, 2]), :] = preprocess.make_gts_matrix(gts, par_status, gt_indices, imp_gts=imp_gts, parsum=parsum)
        else:
            G = np.zeros((ids.shape[0],4,gts.shape[1]), dtype=np.float32)
            G[:,np.array([0,2,3]),:] = preprocess.make_gts_matrix(gts, par_status, gt_indices, imp_gts=imp_gts, parsum=parsum)
        G[:,1,:] = preprocess.get_fam_means(ids, ped, gts, gts_ids, remove_proband=True).gts
    else:
        G = preprocess.make_gts_matrix(gts, par_status, gt_indices, parsum=parsum, imp_gts=imp_gts)
    del gts
    if imp_gts is not None:
        del imp_gts
    return gtarray(G, ids, sid, alleles=alleles, pos=pos, chrom=chromosome, fams=fam_labels, par_status=par_status)

def read_sibs_from_bed(bedfile,sibpairs):
    bed = Bed(bedfile, count_A1=True)
    ids = bed.iid
    id_dict = make_id_dict(ids, 1)
    # Find sibpairs in bed
    in_bed = np.vstack((np.array([x in id_dict for x in sibpairs[:,0]]),
                        np.array([x in id_dict for x in sibpairs[:, 1]]))).T
    both_in_bed = np.sum(in_bed,axis=1)==2
    # Remove pairs without both in bedfile
    if np.sum(both_in_bed)<sibpairs.shape[0]:
        print(str(sibpairs.shape[0]-np.sum(both_in_bed))+' sibpairs do not both have genotypes')
        sibpairs = sibpairs[both_in_bed,:]
    # Find indices of sibpairs
    sibindices = np.sort(np.array([id_dict[x] for x in sibpairs.flatten()]))
    gts = np.zeros((sibindices.shape[0],bed.sid.shape[0]),dtype=np.float32)
    gts[:] = bed[sibindices,:].read().val
    return gtarray(garray = gts, ids = ids[sibindices, 1], sid = bed.sid, pos = np.array(bed.pos[:,2],dtype=int))

