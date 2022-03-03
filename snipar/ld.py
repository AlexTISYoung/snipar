from numba import njit, prange
import numpy as np
from snipar.utilities import make_id_dict
from snipar.map import map_from_bed
from pysnptools.snpreader import Bed

#### Compute LD-scores ####
@njit(parallel=True)
def compute_ld_scores(gts,map,max_dist = 1):
    ldscores = np.zeros((gts.shape[1]),dtype=np.float64)
    for i in prange(gts.shape[1]):
        ldscore_i = 1
        if i>0:
            j = i-1
            dist = map[i]-map[j]
            while dist < max_dist and j>=0:
                ldscore_i += r2_est(gts[...,i],gts[...,j])
                j -= 1
                if j>=0:
                    dist = map[i]-map[j]
        if i<(gts.shape[1]-1):
            j = i + 1
            dist = map[j] - map[i]
            while dist < max_dist and j < gts.shape[1]:
                ldscore_i += r2_est(gts[..., i], gts[..., j])
                j += 1
                if j < gts.shape[1]:
                    dist = map[j] - map[i]
        ldscores[i] = ldscore_i
    return ldscores

## Unbiased estimator of R^2 between SNPs
@njit
def r2_est(g1,g2):
    not_nan = np.logical_not(np.logical_or(np.isnan(g1), np.isnan(g2)))
    r2 = np.power(np.corrcoef(g1[not_nan],g2[not_nan])[0,1],2)
    return r2-(1-r2)/(np.sum(not_nan)-2)

def ldscores_from_bed(bedfile, chrom, ld_wind, ld_out = None):
    # Get map
    map_snps, map = map_from_bed(bedfile, chrom)
    not_na = ~np.isnan(map)
    map = map[not_na]
    map_snps = map_snps[not_na]
    map_snp_dict = make_id_dict(map_snps)
    # Read genotypes
    print('Reading genotypes')
    bed = Bed(bedfile, count_A1 = True)
    bed_in_map = np.array([x in map_snp_dict for x in bed.sid])
    gts = bed[:,bed_in_map].read().val
    sid = bed.sid[bed_in_map]
    bim = bed.pos[bed_in_map,:]
    chrom = np.array(bim[:,0],dtype=int).reshape((sid.shape[0],1))
    pos = np.array(bim[:,2],dtype=int).reshape((sid.shape[0],1))
    map = map[np.array([map_snp_dict[x] for x in sid])]
    print('Computing LD scores')
    ldscores = compute_ld_scores(gts,map,ld_wind)
    if ld_out is not None:
        ld_outfile = ld_out+str(chrom[0])+'.l2.ldscore.gz'
        print('Writing LD-scores to '+ld_outfile)
        ld_outarray = np.hstack((chrom, sid.reshape((sid.shape[0],1)),pos,ldscores.reshape((sid.shape[0],1))))
        ld_outarray = np.vstack((np.array(['CHR', 'SNP', 'BP', 'L2']).reshape((1,4)),
                        ld_outarray))
        np.savetxt(ld_outfile, ld_outarray, fmt='%s')
    return ldscores, sid

