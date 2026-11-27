from numba import njit, prange
import numpy as np
from snipar.utilities import make_id_dict
from snipar.map import map_from_bed
from pysnptools.snpreader import Bed

@njit(cache=True, parallel=True)
def compute_ld_scores(gts, map, max_dist=1):
    # Consecutive samples for one SNP should occupy consecutive memory.
    gts = np.asfortranarray(gts)
    ldscores = np.ones(gts.shape[1], dtype=np.float64)
    for i in prange(gts.shape[1]):
        for j in range(i - 1, -1, -1):
            if map[i] - map[j] >= max_dist:
                break
            ldscores[i] += r2_est(gts[:, i], gts[:, j])
        for j in range(i + 1, gts.shape[1]):
            if map[j] - map[i] >= max_dist:
                break
            ldscores[i] += r2_est(gts[:, i], gts[:, j])
    return ldscores


@njit(cache=True)
def r2_est(g1, g2):
    """Pairwise complete, sample-size corrected r-squared for callable genotypes."""
    n = 0
    x = y = xx = yy = xy = 0.0
    for i in range(g1.size):
        a, b = np.float64(g1[i]), np.float64(g2[i])
        if np.isfinite(a) and np.isfinite(b):
            n += 1
            x += a
            y += b
            xx += a * a
            yy += b * b
            xy += a * b
    if n < 3:
        return 0.0
    # Centered moments without temporary subsets. Tolerate cancellation when
    # one dosage is constant or indistinguishable from constant.
    eps = 4 * np.finfo(np.float64).eps * n
    vx = xx - x*x/n
    vy = yy - y*y/n
    if vx <= eps*xx or vy <= eps*yy:
        return 0.0
    cov = xy - x*y/n
    r2 = min(1.0, cov*cov / (vx*vy))
    return r2 - (1.0 - r2) / (n - 2)

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
    if not bed_in_map.any():
        raise ValueError("No BED SNPs have genetic map positions")
    gts = bed[:, bed_in_map].read(order="F").val
    sid = bed.sid[bed_in_map]
    bim = bed.pos[bed_in_map,:]
    chromosomes = bim[:, :1].astype(int)
    pos = bim[:, 2:3].astype(int)
    map = map[np.array([map_snp_dict[x] for x in sid])]
    if np.any(np.diff(map) < 0):
        raise ValueError("SNPs must be in increasing genetic map order")
    print('Computing LD scores')
    ldscores = compute_ld_scores(gts,map,ld_wind)
    if ld_out is not None:
        ld_outfile = f"{ld_out}{chromosomes[0, 0]}.l2.ldscore.gz"
        print('Writing LD-scores to '+ld_outfile)
        ld_outarray = np.hstack((chromosomes, sid.reshape((sid.shape[0],1)),pos,ldscores.reshape((sid.shape[0],1))))
        ld_outarray = np.vstack((np.array(['CHR', 'SNP', 'BP', 'L2']).reshape((1,4)),
                        ld_outarray))
        np.savetxt(ld_outfile, ld_outarray, fmt='%s')
    return ldscores, sid