import snipar.preprocess as preprocess
import numpy as np
from pysnptools.snpreader import Bed
from snipar.gtarray import gtarray
from snipar.utilities import *
from snipar.read.variants import match_variants


def _match_variants(gts_f, par_gts_f, bim, snp_ids=None, start=0, end=None):
    alleles = np.loadtxt(bim, dtype=str, usecols=(4, 5), ndmin=2)
    return match_variants(gts_f.sid, gts_f.pos[:, 0].astype(int), gts_f.pos[:, 2].astype(int),
                          alleles, par_gts_f, snp_ids, start, end)


def match_observed_and_imputed_snps(gts_f, par_gts_f, bim, snp_ids=None, start=0, end=None):
    matched = _match_variants(gts_f, par_gts_f, bim, snp_ids, start, end)
    return matched.with_imputed_mask(par_gts_f['bim_values'].shape[0])


def get_snps(gts_f,bim,snp_ids=None, match=False):
    # Match SNPs in bed and get SNP info
    if snp_ids is None:
        snp_ids = gts_f.sid
    # Get bim info
    alleles = np.loadtxt(bim, dtype='U', usecols=(4, 5), ndmin=2)
    pos = np.loadtxt(bim, dtype=int, usecols=3, ndmin=1)
    chromosome = np.loadtxt(bim, dtype=int, usecols=0, ndmin=1)
    # Remove duplicate ids
    if not match: # this will change the ordering of snps, which affects the robust estimator
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


def get_gts_matrix_given_ped(ped, imp_fams, bedfile, par_gts_f=None ,snp_ids=None, ids=None, sib=False, sib_diff=False, parsum=False, start=0, end=None, include_unrel=False, robust=False, trios_sibs=False, match_snp_ids=False, verbose=False, print_sample_info = False, matched_variants=None):
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
    # if par_gts_f is not None:
    #     imp_fams = convert_str_array(np.array(par_gts_f['families']))
    # else:
    #     imp_fams = None
    if sib_diff:
        # Find ids with sibs
        ids, observed_indices = preprocess.get_indices_given_ped_sibs(ped, gts_ids, ids=ids, verbose=False)
    elif trios_sibs:
        ids, observed_indices, trios_indices, sibs_indices = preprocess.get_indices_given_ped_trios_sibs(ped, gts_ids, ids=ids, verbose=False)
    else:
        # Find ids with observed/imputed parents and indices of those in observed/imputed data
        ids, observed_indices, imp_indices, parcount = preprocess.get_indices_given_ped(ped, gts_ids, imp_fams=imp_fams, ids=ids, 
                                                                                    sib=sib or sib_diff, include_unrel=include_unrel, verbose=print_sample_info) 
        if np.sum(parcount>0)==0 and not parsum:
            if verbose:
                print('No individuals with genotyped parents found. Using sum of imputed maternal and paternal genotypes to prevent collinearity.')
            parsum = True
        elif 100 > np.sum(parcount>0) > 0 and not parsum:
            if verbose:
                print('Warning: low number of individuals with observed parental genotypes. Consider using the --parsum argument to prevent issues due to collinearity.')
    ### Match observed and imputed SNPs ###
    if par_gts_f is not None and not sib_diff:
        if verbose:
            print('Matching observed and imputed SNPs')
        if matched_variants is None:
            matched_variants = _match_variants(gts_f, par_gts_f, bim, snp_ids, start, end)
        chromosome, sid, pos, alleles, allele_flip, imp_columns, obs_sid_index = matched_variants
        # Read imputed parental genotypes
        if verbose:
            print('Reading imputed parental genotypes')
        if imp_indices.size * par_gts_f['imputed_par_gts'].shape[1] < imp_columns.size * imp_fams.size:
            imp_gts = np.array(par_gts_f['imputed_par_gts'][imp_indices, :])
            imp_gts = imp_gts[:,imp_columns]
            if robust:
                num_obs_par_al = np.array(par_gts_f['num_observed_parental_alleles'][imp_indices, :])
                num_obs_par_al = num_obs_par_al[:,imp_columns]
        else:
            imp_gts = np.array(par_gts_f['imputed_par_gts'][:,imp_columns])
            imp_gts = imp_gts[imp_indices,:]
            if robust:
                num_obs_par_al = np.array(par_gts_f['num_observed_parental_alleles'][:,imp_columns])
                num_obs_par_al = num_obs_par_al[imp_indices,:]
        imp_fams = imp_fams[imp_indices]
        # Check for allele flip
        nflip = np.sum(allele_flip)
        if nflip>0:
            if verbose:
                print('Flipping alleles of '+str(nflip)+' SNPs to match observed genotypes')
            imp_gts[:,allele_flip] = 2-imp_gts[:,allele_flip]
    else:
        chromosome, sid, pos, alleles, obs_sid_index = get_snps(gts_f, bim, snp_ids=snp_ids, match=match_snp_ids)
        imp_gts = None
    # Read observed genotypes
    if verbose:
        print('Reading observed genotypes')
    # gts = gts_f[observed_indices, obs_sid_index].read().val
    gts = gts_f[:, obs_sid_index].read().val
    gts_ids_reduced = gts_f.iid[observed_indices,1]
    gts_id_dict = make_id_dict(gts_ids_reduced)
    # Find indices in reduced data
    if sib_diff:
        gt_indices, fam_labels = preprocess.find_gts(ids, ped, gts_id_dict)
        par_status = None
    elif trios_sibs:
        # par_status, gt_indices = preprocess.find_trio_gts(ids, ped, gts_id_dict)
        par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict)
    else:
        par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict, imp_fams=imp_fams)
    if verbose:
        print('Constructing family based genotype matrix')
    ### Make genotype design matrix
    if sib:
        if parsum:
            G = np.zeros((ids.shape[0], 3, gts.shape[1]), dtype=np.float32)
            G[:, np.array([0, 2]), :] = preprocess.make_gts_matrix(gts[observed_indices], par_status, gt_indices, imp_gts=imp_gts, parsum=parsum)
        else:
            G = np.zeros((ids.shape[0],4,gts.shape[1]), dtype=np.float32)
            G[:,np.array([0,2,3]),:] = preprocess.make_gts_matrix(gts[observed_indices], par_status, gt_indices, imp_gts=imp_gts, parsum=parsum)
        # need to use original gts and gts_ids, since some of the sibs might not have phenotypes
        # gts_all = np.sum(gts_f.read((None,obs_sid_index), np.float32)[:,:,np.array([0,2])],axis=2)
        G[:,1,:] = preprocess.get_fam_means(ids, ped, gts, gts_ids, remove_proband=True).gts
    elif sib_diff:
        G = np.full((ids.shape[0], 2, gts.shape[1]), fill_value=-1, dtype=np.float32)
        G[:,0,:] = gts[observed_indices][gt_indices,:]
        # need to use original gts and gts_ids, since some of the sibs might not have phenotypes
        # gts_all = np.sum(gts_f.read((None,obs_sid_index), np.float32)[:,:,np.array([0,2])],axis=2)
        G[:,1,:] = preprocess.get_fam_means(ids, ped, gts, gts_ids, remove_proband=False).gts
    elif trios_sibs:
        # gts_all = np.sum(gts_f.read((None,obs_sid_index), np.float32)[:,:,np.array([0,2])],axis=2)
        # if parsum:
        #     G = np.zeros((ids.shape[0], 2, gts.shape[1]), dtype=np.float32)
        #     if sibs_indices.shape[0] > 0:
        #         G[sibs_indices,1,:] = preprocess.get_fam_means(ids[sibs_indices], ped, gts[observed_indices], gts_ids, remove_proband=False).gts
        # else:
        trios_indices = np.arange(0, trios_indices.shape[0]) # trios are in the upper rows of observed_indices based on construction 
        sibs_indices = trios_indices.shape[0] + np.arange(0, sibs_indices.shape[0])
        G = np.zeros((ids.shape[0],3,gts.shape[1]), dtype=np.float32)
        G[:, :, :] = preprocess.make_gts_matrix(gts[observed_indices], par_status, gt_indices, imp_gts=None, parsum=False)
        if sibs_indices.shape[0] > 0:
            G[sibs_indices,0,:] = gts[observed_indices][sibs_indices,:]
            G[sibs_indices,1,:] = preprocess.get_fam_means(ids[sibs_indices], ped, gts, gts_ids, remove_proband=False).gts
            # if not parsum:
            G[sibs_indices,2,:] = G[sibs_indices,1,:]
    else:
        # G = preprocess.make_gts_matrix(gts, par_status, gt_indices, imp_gts=imp_gts, parsum=parsum)
        G = preprocess.make_gts_matrix(gts[observed_indices], par_status, gt_indices, imp_gts=imp_gts, parsum=parsum)
    del gts
    if imp_gts is not None:
        del imp_gts
    if robust:
        num_obs_par_al = preprocess.make_num_obs_par_al_matrix(num_obs_par_al, par_status, gt_indices, G.shape[0])
        return gtarray(G, ids, sid, alleles=alleles, pos=pos, chrom=chromosome, fams=fam_labels, par_status=par_status, num_obs_par_al=num_obs_par_al) #, ped=ped, standard_f=standard_f)
    if trios_sibs:
        G = gtarray(G, ids, sid, alleles=alleles, pos=pos, chrom=chromosome, fams=fam_labels, par_status=par_status)
        G.complete_trios_inds = trios_indices
        G.sibs_inds = sibs_indices
        return G
    return gtarray(G, ids, sid, alleles=alleles, pos=pos, chrom=chromosome, fams=fam_labels, par_status=par_status, ped=ped)


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
