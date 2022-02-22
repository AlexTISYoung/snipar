from snipar.read.bed import read_PO_pairs_from_bed
from snipar.read.bgen import read_PO_pairs_from_bgen
import numpy as np
from snipar.utilities import *
from numba import njit, prange

class g_error(object):
    def __init__(self, error_ests, ME, sum_het, sid):
        if error_ests.shape[0] == ME.shape[0] and ME.shape[0] == sum_het.shape[0] and sum_het.shape[0] ==sid.shape[0]:
            self.error_ests = error_ests
            self.ME = ME
            self.sum_het = sum_het
            self.sid = sid
            self.sid_dict = make_id_dict(self.sid)
    def bayes_shrink(self, alpha, beta):
        self.error_ests = (self.ME+alpha)/(self.sum_het+beta)

def estimate_genotyping_error_rate(ped, bedfiles=None, bgenfiles=None, min_maf=0.01):
    genome_errors = []
    # Estimate per-SNP errors for each chromosome
    if bedfiles is not None:
        nsnp = np.zeros((bedfiles.shape[0]), dtype=int)
        for i in range(bedfiles.shape[0]):
            ME_chr = mendelian_errors(ped, bedfile=bedfiles[i], min_maf=min_maf)
            genome_errors.append(ME_chr)
            nsnp[i] = ME_chr.sid.shape[0]
    elif bgenfiles is not None:
        nsnp = np.zeros((bgenfiles.shape[0]), dtype=int)
        for i in range(bgenfiles.shape[0]):
            ME_chr = mendelian_errors(ped, bgenfile=bgenfiles[i], min_maf=min_maf)
            genome_errors.append(ME_chr)
            nsnp[i] = ME_chr.sid.shape[0]
    else:
        raise(ValueError('Must provide bed files or bgen files'))
    ## Estimate empirical bayes prior parameters
    # Collect MLE error rates across genome
    genome_error_rates = np.zeros((np.sum(nsnp)))
    sum_het = np.zeros((np.sum(nsnp)))
    snp_start = 0
    for i in range(nsnp.shape[0]):
        snp_end = snp_start+nsnp[i]
        genome_error_rates[snp_start:snp_end] = genome_errors[i].error_ests
        sum_het[snp_start:snp_end] = genome_errors[i].sum_het
        snp_start = snp_end
    # Estimate mean and variance of MLE error estimates
    mean_error = np.mean(genome_error_rates)
    var_error = np.var(genome_error_rates)
    if var_error > 0:
        mean_inv_het = np.mean(1/sum_het)
        beta = 1/(var_error/mean_error-mean_inv_het)
        if beta > 0:
            alpha = mean_error*beta
            for i in range(nsnp.shape[0]):
                genome_errors[i].bayes_shrink(alpha, beta)
            return mean_error, genome_errors
        else:
            return mean_error, None
    else:
        return mean_error, None

def mendelian_errors(ped, bedfile=None, bgenfile=None, min_maf=0.01):
    if bedfile is not None:
        gts, opg_ped, npair = read_PO_pairs_from_bed(ped, bedfile=bedfile)
    elif bgenfile is not None:
        gts, opg_ped, npair = read_PO_pairs_from_bgen(ped, bgenfile=bgenfile)
    #print('Finding indices of parent-offspring pairs')
    ## Get indices
    pair_indices = np.zeros((npair,2),dtype=int)
    pair_count = 0
    for i in range(opg_ped.shape[0]):
        o_index = gts.id_dict[opg_ped[i,1]]
        if opg_ped[i, 2] in gts.id_dict:
            pair_indices[pair_count,:] = np.array([o_index,gts.id_dict[opg_ped[i,2]]])
            pair_count += 1
        if opg_ped[i, 3] in gts.id_dict:
            pair_indices[pair_count,:] = np.array([o_index,gts.id_dict[opg_ped[i,3]]])
            pair_count += 1
    # Filter on MAF
    #print('Filtering on MAF')
    gts.filter_maf(min_maf)
    ## Count Mendelian errors
    #print('Counting mendelain errors')
    ME = count_ME(np.array(gts.gts,dtype=np.float_), pair_indices)
    #print('Counted mendelain errors')
    # Estimate error probability
    N_pair = np.sum(np.logical_not(gts.gts[pair_indices[:, 0], :].mask) * np.logical_not(gts.gts[pair_indices[:, 1], :].mask),
                    axis=0)
    sum_het = N_pair * gts.freqs * (1 - gts.freqs)
    error_mle = ME/sum_het
    return g_error(error_mle, ME, sum_het, gts.sid)

@njit(parallel=True)
def count_ME(gts,pair_indices):
    ME = np.zeros((gts.shape[1]), dtype=np.int_)
    # Count Mendelian errors
    for j in prange(gts.shape[1]):
        for i in range(pair_indices.shape[0]):
            if np.abs(gts[pair_indices[i, 0], j] - gts[pair_indices[i, 1], j]) > 1:
                ME[j] += 1
    return ME