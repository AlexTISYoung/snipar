from numba import get_num_threads, njit, prange, types
from numba.extending import intrinsic
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
    return _r2_from_moments(n, x, y, xx, yy, xy)


@njit(cache=True, inline='always')
def _r2_from_moments(n, x, y, xx, yy, xy):
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


@intrinsic
def _popcount(typing_context, word):
    signature = types.int64(types.uint64)

    def codegen(context, builder, signature, args):
        function = builder.module.declare_intrinsic('llvm.ctpop', [args[0].type])
        return builder.call(function, args)

    return signature, codegen


@njit(cache=True, parallel=True)
def pack_ld_hardcalls(values, output):
    """Pack sample × marker 0/1/2 calls and negative missing calls into uint64 planes.

    The output has shape (markers, 3, ceil(samples / 64)), with called,
    at-least-one and two copies as the three planes. Blocks can be written into
    disjoint marker slices of a preallocated array or memory map.
    """
    if output.shape != (values.shape[1], 3, (values.shape[0] + 63) // 64):
        raise ValueError('Packed LD output shape differs from the genotype shape')
    for marker in prange(values.shape[1]):
        for word in range(output.shape[2]):
            start = word * 64
            stop = min(start + 64, values.shape[0])
            called = once = twice = np.uint64(0)
            for sample in range(start, stop):
                call = values[sample, marker]
                if call >= 0:
                    bit = np.uint64(1) << np.uint64(sample - start)
                    called |= bit
                    if call > 0:
                        once |= bit
                        if call == 2:
                            twice |= bit
            output[marker, 0, word] = called
            output[marker, 1, word] = once
            output[marker, 2, word] = twice


@njit(cache=True)
def _packed_moments(packed, columns):
    moments = np.zeros((len(columns), 3), dtype=np.int64)
    ever_called = np.zeros(packed.shape[2], dtype=np.uint64)
    for i, column in enumerate(columns):
        n = x1 = x2 = 0
        for word in range(packed.shape[2]):
            called = packed[column, 0, word]
            ever_called[word] |= called
            n += _popcount(called)
            x1 += _popcount(packed[column, 1, word])
            x2 += _popcount(packed[column, 2, word])
        moments[i, 0] = n
        moments[i, 1] = x1 + x2
        moments[i, 2] = x1 + 3*x2
    n_ever_called = 0
    for word in ever_called:
        n_ever_called += _popcount(word)
    return moments, n_ever_called


@njit(cache=True, inline='always')
def _packed_product(x1, x2, y1, y2):
    # g = I[g >= 1] + I[g == 2]. Missing genotypes have both bits clear.
    return _popcount(x1 & y1) + _popcount(x1 & y2) + _popcount(x2 & y1) + _popcount(x2 & y2)


@njit(cache=True, inline='always')
def _packed_r2(first, second, first_moments, second_moments, ever_called):
    xy = 0
    # When both sets have the size of their union across all markers, both
    # markers have identical called sets. Reuse their centered input moments.
    if first_moments[0] == ever_called and second_moments[0] == ever_called:
        for word in range(first.shape[1]):
            xy += _packed_product(first[1, word], first[2, word], second[1, word], second[2, word])
        return _r2_from_moments(ever_called, np.float64(first_moments[1]), np.float64(second_moments[1]),
                                np.float64(first_moments[2]), np.float64(second_moments[2]), np.float64(xy))

    n = x1 = x2 = y1 = y2 = 0
    for word in range(first.shape[1]):
        called_x, called_y = first[0, word], second[0, word]
        ax1, ax2 = first[1, word], first[2, word]
        ay1, ay2 = second[1, word], second[2, word]
        n += _popcount(called_x & called_y)
        x1 += _popcount(ax1 & called_y)
        x2 += _popcount(ax2 & called_y)
        y1 += _popcount(ay1 & called_x)
        y2 += _popcount(ay2 & called_x)
        xy += _packed_product(ax1, ax2, ay1, ay2)
    return _r2_from_moments(n, np.float64(x1 + x2), np.float64(y1 + y2),
                            np.float64(x1 + 3*x2), np.float64(y1 + 3*y2), np.float64(xy))


@njit(cache=True, parallel=True)
def compute_ld_scores_packed(packed, map, max_dist=1, columns=None):
    """Compute LD scores for packed hard calls in increasing genetic-map order."""
    if columns is None:
        columns = np.arange(packed.shape[0])
    moments, ever_called = _packed_moments(packed, columns)
    scores = np.ones(len(columns), dtype=np.float64)
    for i in prange(len(columns)):
        first, first_moments = packed[columns[i]], moments[i]
        score = 1.0
        for j in range(i - 1, -1, -1):
            if map[i] - map[j] >= max_dist:
                break
            score += _packed_r2(first, packed[columns[j]], first_moments, moments[j], ever_called)
        for j in range(i + 1, len(columns)):
            if map[j] - map[i] >= max_dist:
                break
            score += _packed_r2(first, packed[columns[j]], first_moments, moments[j], ever_called)
        scores[i] = score
    return scores


@njit(cache=True, parallel=True)
def _ld_partner_tile(center, partner, center_moments, partner_moments, center_map, partner_map,
                     center_start, partner_start, max_dist, samples, scores):
    for i in prange(center.shape[0]):
        first, first_moments = center[i], center_moments[i]
        score = scores[i]
        for j in range(partner.shape[0]):
            if center_map[i] - partner_map[j] >= max_dist:
                continue
            if partner_map[j] - center_map[i] >= max_dist:
                break
            if center_start + i != partner_start + j:
                score += _packed_r2(first, partner[j], first_moments, partner_moments[j], samples)
        scores[i] = score


def _read_packed_ld(source, markers, words):
    packed = np.empty((len(markers), 3, words), dtype=np.uint64)
    pack_ld_hardcalls(source.read(slice(None), markers), packed)
    return packed


def _bounded_bed_ld(source, markers, positions, max_dist, batch_size=None):
    samples = len(source.ids)
    words = (samples + 63) // 64
    # At most two packed tiles, their moments and one raw int8 input tile.
    _, width = plan_batches(len(markers), samples + 48 * (words + 1), batch_size=batch_size)
    if width == len(markers):
        return compute_ld_scores_packed(_read_packed_ld(source, markers, words), positions, max_dist)
    scores = np.ones(len(markers))
    for start in range(0, len(markers), width):
        end = min(start + width, len(markers))
        center = _read_packed_ld(source, markers[start:end], words)
        center_moments, _ = _packed_moments(center, np.arange(len(center)))
        # The cM search bounds avoid reading tiles with no relevant partners.
        first = np.searchsorted(positions, positions[start] - max_dist, side='left')
        last = np.searchsorted(positions, positions[end - 1] + max_dist, side='right')
        for offset in range(first, last, width):
            stop = min(offset + width, last)
            if offset == start and stop == end:
                partner, partner_moments = center, center_moments
            else:
                partner = _read_packed_ld(source, markers[offset:stop], words)
                partner_moments, _ = _packed_moments(partner, np.arange(len(partner)))
            _ld_partner_tile(center, partner, center_moments, partner_moments,
                             positions[start:end], positions[offset:stop], start, offset,
                             max_dist, samples, scores[start:end])
            del partner, partner_moments
        del center, center_moments
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
