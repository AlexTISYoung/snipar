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
    imp_bim_cols = convert_str_array(np.array(par_gts_f['bim_columns']))
    imp_bim = convert_str_array(np.array(par_gts_f['bim_values']))
    # Get imputed SNP ids
    found_snp_ids = False
    if 'rsid' in imp_bim_cols:
        imp_sid = imp_bim[:,np.where(imp_bim_cols=='rsid')[0][0]]
        if np.unique(imp_sid).shape[0] == 0:
            found_snp_ids = False
        else:
            found_snp_ids = True
    if not found_snp_ids:
        if 'id' in imp_bim_cols:
            imp_sid = imp_bim[:,np.where(imp_bim_cols=='id')[0][0]]
        else:
            raise(ValueError('Cannot find imputed SNP ids'))
    # Get imputed allele ids
    if 'allele_ids' in imp_bim_cols:
        imp_alleles = np.array([x.split(',') for x in imp_bim[:,np.where(imp_bim_cols=='allele_ids')[0][0]]])
    elif 'allele1' in imp_bim_cols and 'allele2' in imp_bim_cols:
        imp_alleles = imp_bim[:,[np.where(imp_bim_cols=='allele1')[0][0],np.where(imp_bim_cols=='allele2')[0][0]]]
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
    imp_alleles = imp_alleles[in_obs_sid,:]
    chromosome = chromosome[obs_sid_index]
    pos = pos[obs_sid_index]
    # Check for allele flip/mismatch
    allele_match = np.logical_and(alleles[:,0]==imp_alleles[:,0],alleles[:,1]==imp_alleles[:,1])
    if np.sum(allele_match) < alleles.shape[0]:
        allele_flip = np.logical_and(alleles[:,0]==imp_alleles[:,1],alleles[:,1]==imp_alleles[:,0])
        allele_mismatch = np.logical_not(np.logical_or(allele_match,allele_flip))
        n_mismatch = np.sum(allele_mismatch)
        if n_mismatch == alleles.shape[0]:
            raise(ValueError('Np alleles match between observed and imputed genotypes'))
        elif n_mismatch > 0:
            print('Removing '+str(n_mismatch)+' SNPs with mismatched alleles between imputed and observed')
            chromosome = chromosome[~allele_mismatch]
            sid = sid[~allele_mismatch]
            pos = pos[~allele_mismatch]
            alleles = alleles[~allele_mismatch]
            imp_alleles = imp_alleles[~allele_mismatch]
            in_obs_sid[np.where(in_obs_sid)[0][allele_mismatch]] = False
            obs_sid_index =  obs_sid_index[~allele_mismatch]
    allele_flip = np.logical_and(alleles[:,0]==imp_alleles[:,1],alleles[:,1]==imp_alleles[:,0])
    return chromosome, sid, pos, alleles, allele_flip, in_obs_sid, obs_sid_index

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
    if ids is None:
        ids = gts_ids
    # Get families with imputed parental genotypes
    if par_gts_f is not None:
        imp_fams = convert_str_array(par_gts_f['families'])
    else:
        imp_fams = None
    ### Find ids with observed/imputed parents and indices of those in observed/imputed data
    ids, observed_indices, imp_indices, parcount = preprocess.get_indices_given_ped(ped, gts_ids, imp_fams=imp_fams, ids=ids, 
                                                                                sib=sib, verbose=print_sample_info)
    if np.sum(parcount>0)==0 and not parsum:
        if verbose:
            print('No individuals with genotyped parents found. Using sum of imputed maternal and paternal genotypes to prevent collinearity.')
        parsum = True
    elif 100 > np.sum(parcount>0) > 0 and not parsum:
        if verbose:
            print('Warning: low number of individuals with observed parental genotypes. Consider using the --parsum argument to prevent issues due to collinearity.')
    ### Match observed and imputed SNPs ###
    if par_gts_f is not None:
        if verbose:
            print('Matching observed and imputed SNPs')
        chromosome, sid, pos, alleles, allele_flip, in_obs_sid, obs_sid_index = match_observed_and_imputed_snps(gts_f, par_gts_f, bim, snp_ids=snp_ids)
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
        # Check for allele flip
        nflip = np.sum(allele_flip)
        if nflip>0:
            print('Flipping alleles of '+str(nflip)+' SNPs to match observed genotypes')
            imp_gts[:,allele_flip] = 2-imp_gts[:,allele_flip]
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

def read_PO_pairs_from_bed(ped,bedfile):
    # Read bed
    bed = Bed(bedfile, count_A1=True)
    ids = bed.iid
    id_dict = make_id_dict(ids, 1)
    ## Find parent-offspring pairs
    # genotyped individuals
    genotyped = np.array([x in id_dict for x in ped[:, 1]])
    ped = ped[genotyped, :]
    # with genotyped father
    father_genotyped = np.array([x in id_dict for x in ped[:, 2]])
    # with genotyped mother
    mother_genotyped = np.array([x in id_dict for x in ped[:, 3]])
    # either
    opg = np.logical_or(father_genotyped, mother_genotyped)
    opg_ped = ped[opg, :]
    # number of pairs
    npair = np.sum(father_genotyped) + np.sum(mother_genotyped)
    if npair == 0:
        raise(ValueError('No parent-offspring pairs in  '+str(bedfile)+' for genotype error probability estimation'))
    print(str(npair)+' parent-offspring pairs found in '+bedfile)
    if npair*bed.sid.shape[0] < 10**5:
        print('Warning: limited information for estimation of genotyping error probability.')
    ## Read genotypes
    all_ids = np.unique(np.hstack((opg_ped[:, 1],
                                   ped[father_genotyped, 2],
                                   ped[mother_genotyped, 3])))
    all_ids_indices = np.sort(np.array([id_dict[x] for x in all_ids]))
    gts = bed[all_ids_indices, :].read().val
    return gtarray(gts,ids = ids[all_ids_indices, 1], sid=bed.sid), opg_ped, npair
