from numba import get_num_threads, njit, prange
import numpy as np
from snipar.utilities import make_id_dict
from snipar.map import map_from_bed
from snipar.observed import open_observed
from snipar.resources import available_cpus, plan_batches

@njit(cache=True, parallel=True)
def compute_ld_scores(gts, map, max_dist=1, columns=None):
    # Consecutive samples for one SNP should occupy consecutive memory.
    gts = np.asfortranarray(gts)
    if columns is None:
        columns = np.arange(gts.shape[1])
    ldscores = np.ones(len(columns), dtype=np.float64)
    for i in prange(len(columns)):
        for j in range(i - 1, -1, -1):
            if map[i] - map[j] >= max_dist:
                break
            ldscores[i] += r2_est(gts[:, columns[i]], gts[:, columns[j]])
        for j in range(i + 1, len(columns)):
            if map[j] - map[i] >= max_dist:
                break
            ldscores[i] += r2_est(gts[:, columns[i]], gts[:, columns[j]])
    return ldscores


@njit(cache=True)
def r2_est(g1, g2):
    """Pairwise complete, sample-size corrected r-squared for callable genotypes."""
    n = 0
    x = y = xx = yy = xy = 0.0
    for i in range(g1.size):
        a, b = np.float64(g1[i]), np.float64(g2[i])
        if np.isfinite(a) and np.isfinite(b) and a != -127 and b != -127:
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

@njit(cache=True, parallel=True)
def _ld_partner_tile(center, partner, center_map, partner_map, center_start, partner_start, max_dist, scores):
    for i in prange(center.shape[1]):
        for j in range(partner.shape[1]):
            if center_map[i] - partner_map[j] >= max_dist:
                continue
            if partner_map[j] - center_map[i] >= max_dist:
                break
            if center_start + i != partner_start + j:
                scores[i] += r2_est(center[:, i], partner[:, j])


def _bounded_bed_ld(source, markers, positions, max_dist, batch_size=None):
    _, width = plan_batches(len(markers), 2 * len(source.ids), batch_size=batch_size)
    if width == len(markers):
        return compute_ld_scores(source.read(slice(None), markers), positions, max_dist)
    scores = np.ones(len(markers))
    for start in range(0, len(markers), width):
        end = min(start + width, len(markers))
        center = source.read(slice(None), markers[start:end])
        # The cM search bounds avoid reading tiles with no relevant partners.
        first = np.searchsorted(positions, positions[start] - max_dist, side='left')
        last = np.searchsorted(positions, positions[end - 1] + max_dist, side='right')
        for offset in range(first, last, width):
            stop = min(offset + width, last)
            partner = center if offset == start and stop == end else source.read(slice(None), markers[offset:stop])
            _ld_partner_tile(center, partner, positions[start:end], positions[offset:stop],
                             start, offset, max_dist, scores[start:end])
            del partner
        del center
    return scores


def ldscores_from_bed(bedfile, chrom, ld_wind, ld_out=None):
    map_snps, genetic_map = map_from_bed(bedfile, chrom)
    valid = ~np.isnan(genetic_map)
    map_snp_dict = make_id_dict(map_snps[valid])
    known_map = genetic_map[valid]
    with open_observed(bedfile=bedfile, num_threads=min(get_num_threads(), available_cpus())) as source:
        matches = np.fromiter((map_snp_dict.get(sid, -1) for sid in source.sid), dtype=int, count=len(source.sid))
        variants = np.flatnonzero(matches >= 0)
        if not variants.size:
            raise ValueError('No BED SNPs have genetic map positions')
        genetic_map = known_map[matches[variants]]
        if np.any(np.diff(genetic_map) < 0):
            raise ValueError('SNPs must be in increasing genetic map order')
        print('Computing LD scores')
        scores = _bounded_bed_ld(source, variants, genetic_map, ld_wind)
        sid = source.sid[variants]
        if ld_out is not None:
            import gzip

            filename = f'{ld_out}{chrom}.l2.ldscore.gz'
            print('Writing LD-scores to ' + filename)
            with gzip.open(filename, 'wt', compresslevel=3) as output:
                output.write('CHR SNP BP L2\n')
                output.writelines(f'{chrom} {name} {pos} {score}\n' for name, pos, score in
                                  zip(sid, source.reader.bp_position[variants], scores))
    return scores, sid
