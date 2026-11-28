from contextlib import ExitStack
import gzip
from pathlib import Path
import re
from shutil import disk_usage
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from typing import NamedTuple
from numba import get_num_threads, njit, prange
from snipar.map import decode_map_from_pos, get_map_positions
from snipar.ld import compute_ld_scores_packed, pack_ld_hardcalls
import numpy as np
from snipar.observed import allele_stats, open_observed
from snipar.read.family import IdIndex
from snipar.resources import plan_batches, get_memory_budget
from snipar.utilities import outfile_name

####### Transition Matrix ######
@njit
def transition_matrix(cM):
    """Compute probabilities of transitioning between IBD states as a function of genetic distance.
    Args:
        cM : :class:`float`
            genetic distance in centiMorgans (cM)
    Returns:
        log(P) : :class:`~numpy:numpy.array`
            the natural logarithm (element-wise) of the matrix of transition probabilities
    """
    P = np.identity(3)
    r_array = np.array([[-1.0,1.0,0.0],[0.5,-1.0,0.5],[0.0,1.0,-1.0]],dtype=np.float64)
    r = 1.0 - np.power((1.0 + np.exp(-cM / 50)) / 2.0, 4)
    P += r_array*r
    return np.log(P)

#
@njit
def p_ibd_0(f):
    """Compute Joint-PMF for sibling pair genotypes given IBD0.
    Args:
        f : :class:`float`
            allele frequency
    Returns:
        P : :class:`~numpy:numpy.array`
            matrix of probabilities
    """
    P_vec = np.array([(1.0-f)**2.0,2.0*f*(1-f),f**2.0]).reshape((3,1))
    return P_vec @ P_vec.T

@njit
def p_ibd_1(f):
    """Compute Joint-PMF for sibling pair genotypes given IBD1.
    Args:
        f : :class:`float`
            allele frequency
    Returns:
        P : :class:`~numpy:numpy.array`
            matrix of probabilities
    """
    q = 1 - f
    return np.array([[q**3, f*q*q, 0.0],
                     [f*q*q, f*q, f*f*q],
                     [0.0, f*f*q, f**3]])

@njit
def p_ibd_2(f):
    """Compute Joint-PMF for sibling pair genotypes given IBD2.
    Args:
        f : :class:`float`
            allele frequency
    Returns:
        P : :class:`~numpy:numpy.array`
            matrix of probabilities
    """
    P = np.zeros((3,3),dtype=np.float64)
    fdiag = np.array([(1.0-f)**2.0,2.0*f*(1.0-f),f**2.0])
    np.fill_diagonal(P,fdiag.reshape((3,1)))
    return P

@njit
def p_obs_given_IBD(g1_obs,g2_obs,f,p):
    """Compute Joint-PMF for sibling pair genotypes given IBD0.
    Args:
        g1_obs : :class:`integer`
            observed genotype for sibling 1
        g2_obs : :class:`integer`
            observed genotype for sibling 2
        f : :class:`float`
            allele frequency
        p : :class:'float'
            genotyping error probability
    Returns:
        log(P) : :class:`~numpy:numpy.array`
            vector giving log-probabilities of observing g1_obs,g2_obs give IBD 0,1,2
    """
    P_out = np.zeros((3))
    # Get probabilities of observed genotypes given true
    err_vectors = np.array([[1.0 - p, p / 2.0, 0.0],
                            [p, 1.0 - p, p],
                            [0.0,   p / 2.0, 1.0 - p]])
    # Get probabilities of true genotypes given IBD
    probs = np.zeros((3, 3, 3))
    probs[0, ...] = p_ibd_0(f)
    probs[1, ...] = p_ibd_1(f)
    probs[2, ...] = p_ibd_2(f)
    # Compute probabilities conditional on IBD
    for i in range(3):
        P_out[i] = err_vectors[g1_obs, ...].T @ probs[i,...] @ err_vectors[g2_obs, ...]
    # Return
    return np.log(P_out)

@njit
def make_dynamic(g1, g2, freqs, map, weights, error_probs):
    """Make state-matrix and pointer matrix for a sibling pair by dynamic programming
    Args:
        g1 : :class:`~numpy:numpy.array`
            integer vector of first sibling's genotypes
        g2 : :class:`~numpy:numpy.array`
            integer vector of first sibling's genotypes
        freqs : :class:`~numpy:numpy.array`
            floating point vector of allele frequencies
        map : :class:`~numpy:numpy.array`
            floating point vector of genetic positions in cM
        weights : :class:`~numpy:numpy.array`
            floating point vector of SNP weights (usually inverse LD-scores)
        p : :class:'float'
            genotyping error probability
    Returns:
        state_matrix : :class:`~numpy:numpy.array`
            matrix where each column gives the prob of max prob path to that state, where each row is IBD 0,1,2
        pointers : :class:`~numpy:numpy.array`
            integer vector giving the pointer to which state from previous position lead to max at this position
    """
    state_matrix = np.zeros((3, g1.shape[0]), dtype=np.float64)
    pointers = np.zeros((3, g1.shape[0]), dtype=np.int8)
    if not g1.size:
        return state_matrix, pointers
    # Check for nans
    not_nan = np.logical_not(np.logical_or(np.isnan(g1), np.isnan(g2)))
    # Initialise
    state_matrix[:,0] = np.log(np.array([0.25, 0.5, 0.25],dtype=np.float64))
    if not_nan[0]:
        state_matrix[:, 0] += weights[0]*p_obs_given_IBD(np.int8(g1[0]), np.int8(g2[0]), freqs[0], error_probs[0])
    # Compute
    for l in range(1, g1.shape[0]):
        if not_nan[l]:
            probs = weights[l]*p_obs_given_IBD(np.int8(g1[l]), np.int8(g2[l]), freqs[l], error_probs[l])
        else:
            probs = np.zeros((3))
        tmatrix = transition_matrix(map[l]-map[l-1])
        for i in range(3):
            tprobs = tmatrix[:, i]+state_matrix[:, l-1]
            state_matrix[i, l] = np.max(tprobs)+probs[i]
            pointers[i, l] = np.argmax(tprobs)
    return state_matrix, pointers

@njit
def viterbi(state_matrix,pointers):
    """Get viterbi path from state_matrix and pointers output of make_dynamic
    Args:
        state_matrix : :class:`~numpy:numpy.array`
            matrix where each column gives the prob of max prob path to that state, where each row is IBD 0,1,2
        pointers : :class:`~numpy:numpy.array`
            integer vector giving the pointer to which state from previous position lead to max at this position
    Returns:
            path : :class:`~numpy:numpy.array`
            integer vector giving the Viterbi path through the IBD states
    """
    path = np.zeros(state_matrix.shape[1],dtype=np.int8)
    if not path.size:
        return path
    path[path.shape[0]-1] = np.argmax(state_matrix[:, state_matrix.shape[1]-1])
    for i in range(1,path.shape[0]):
        path[path.shape[0]-(i+1)] = pointers[path[path.shape[0]-i], path.shape[0]-i]
    return path

@njit(cache=True)
def _chromosome_hmm(freqs, genetic_map, weights, error_probs):
    n = len(freqs)
    transitions = np.empty((n, 3, 3))
    emissions = np.empty((n, 3, 3, 3))
    for marker in range(n):
        if marker:
            transitions[marker] = transition_matrix(genetic_map[marker] - genetic_map[marker - 1])
        f = freqs[marker]
        q = 1 - f
        e = error_probs[marker]
        # Conditional on the shared haplotype's allele, the unshared alleles
        # of the two siblings are independent draws from the population.
        error = np.array([[1-e, e/2, 0], [e, 1-e, e], [0, e/2, 1-e]])
        for x in range(3):
            x0 = q * error[x, 0] + f * error[x, 1]
            x1 = q * error[x, 1] + f * error[x, 2]
            for y in range(3):
                y0 = q * error[y, 0] + f * error[y, 1]
                y1 = q * error[y, 1] + f * error[y, 2]
                if weights[marker] == 0:
                    emissions[marker, x, y] = 0
                    continue
                probabilities = (
                    (q*x0 + f*x1) * (q*y0 + f*y1),
                    q*x0*y0 + f*x1*y1,
                    q*q*error[x, 0]*error[y, 0] + 2*q*f*error[x, 1]*error[y, 1]
                    + f*f*error[x, 2]*error[y, 2],
                )
                for state in range(3):
                    emissions[marker, x, y, state] = weights[marker] * np.log(probabilities[state])
    return transitions, emissions


@njit(cache=True)
def infer_ibd(sibpairs, gts, freqs, map, weights, error_probs):
    """Return unsmoothed paths for an in-memory dosage matrix."""
    if not gts.shape[1]:
        return np.empty((len(sibpairs), 0), dtype=np.int8)
    transitions, emissions = _chromosome_hmm(freqs, map, weights, error_probs)
    return _infer_prepared(sibpairs, gts, np.arange(gts.shape[1]), transitions, emissions, 0)


@njit(inline='always')
def _call_at(genotypes, row, column, layout):
    if layout == 0:
        value = genotypes[row, column]
        return int(value) if np.isfinite(value) and value >= 0 else -1
    if layout == 1:
        call = (int(genotypes[row, column // 4]) >> (2 * (column % 4))) & 3
        return -1 if call == 3 else int(call)
    call = (int(genotypes[column, row // 4]) >> (2 * (row % 4))) & 3
    return (2, -1, 1, 0)[call]


@njit(inline='always')
def _emission_index(call1, call2):
    high, low = max(call1, call2), min(call1, call2)
    return high * (high + 1) // 2 + low


@njit(cache=True)
def _posterior_tables(transitions, emissions):
    """Make scaled likelihoods with only one entry per unordered genotype pair."""
    n = len(emissions)
    rates = np.zeros((n, 2))
    likelihoods = np.empty((n, 6, 3))
    for marker in range(n):
        if marker:
            rates[marker, 0] = np.exp(transitions[marker, 0, 0])
            rates[marker, 1] = np.exp(transitions[marker, 0, 1])
        for high in range(3):
            for low in range(high + 1):
                values = emissions[marker, high, low]
                peak = max(values[0], values[1], values[2])
                code = _emission_index(high, low)
                for state in range(3):
                    likelihoods[marker, code, state] = np.exp(values[state] - peak) if np.isfinite(peak) else 0
    return rates, likelihoods


@njit(inline='always')
def _scale(p0, p1, p2):
    total = p0 + p1 + p2
    # The forward workspace is float32. If forward and backward support barely
    # overlap, use the log-space fallback before rounding can affect the result.
    if not total > 1e-25:
        return 0., 0., 0.
    inverse = 1. / total
    return p0 * inverse, p1 * inverse, p2 * inverse


@njit(inline='always')
def _condition(a0, a1, a2, marker, call1, call2, likelihoods):
    if min(call1, call2) >= 0:
        index = _emission_index(call1, call2)
        a0 *= likelihoods[marker, index, 0]
        a1 *= likelihoods[marker, index, 1]
        a2 *= likelihoods[marker, index, 2]
    return _scale(a0, a1, a2)


@njit(cache=True)
def _backward_posterior(output, gts, row1, row2, markers, rates, likelihoods, layout):
    b0, b1, b2 = 1., 1., 1.
    for marker in range(len(markers) - 1, -1, -1):
        p0, p1, p2 = _scale(output[marker, 0] * b0, output[marker, 1] * b1, output[marker, 2] * b2)
        if p0 + p1 + p2 == 0:
            return False
        output[marker, 0], output[marker, 1], output[marker, 2] = p0, p1, p2
        if marker:
            call1 = _call_at(gts, row1, markers[marker], layout)
            call2 = _call_at(gts, row2, markers[marker], layout)
            b0, b1, b2 = _condition(b0, b1, b2, marker, call1, call2, likelihoods)
            stay, move = rates[marker, 0], rates[marker, 1]
            b0, b1, b2 = (
                stay * b0 + move * b1,
                move * .5 * (b0 + b2) + stay * b1,
                move * b1 + stay * b2,
            )
    return True


@njit(inline='always')
def _log_add(a, b):
    hi, lo = max(a, b), min(a, b)
    if lo == -np.inf:
        return hi
    return hi + np.log1p(np.exp(lo - hi))


@njit(inline='always')
def _log_step(a0, a1, a2, marker, call1, call2, transitions, emissions, backward):
    if backward and min(call1, call2) >= 0:
        a0 += emissions[marker, call1, call2, 0]
        a1 += emissions[marker, call1, call2, 1]
        a2 += emissions[marker, call1, call2, 2]
    if marker:
        stay = transitions[marker, 0, 0]
        move = transitions[marker, 0, 1]
        half = transitions[marker, 1, 0]
        out, inward = (move, half) if backward else (half, move)
        a0, a1, a2 = (
            _log_add(a0 + stay, a1 + out),
            _log_add(a1 + stay, _log_add(a0, a2) + inward),
            _log_add(a1 + out, a2 + stay),
        )
    if not backward and min(call1, call2) >= 0:
        a0 += emissions[marker, call1, call2, 0]
        a1 += emissions[marker, call1, call2, 1]
        a2 += emissions[marker, call1, call2, 2]
    peak = max(a0, a1, a2)
    if peak == -np.inf:
        return a0, a1, a2
    return a0 - peak, a1 - peak, a2 - peak


_LOG_BLOCK = 256


@njit(cache=True)
def _log_posterior(output, gts, row1, row2, markers, transitions, emissions, layout):
    """Resolve extremely small overlaps in log space with checkpointed alpha vectors."""
    n = len(markers)
    checkpoints = np.empty(((n + _LOG_BLOCK - 1) // _LOG_BLOCK, 3))
    a0, a1, a2 = np.log(.25), np.log(.5), np.log(.25)
    for marker in range(n):
        if marker % _LOG_BLOCK == 0:
            checkpoints[marker // _LOG_BLOCK, 0] = a0
            checkpoints[marker // _LOG_BLOCK, 1] = a1
            checkpoints[marker // _LOG_BLOCK, 2] = a2
        call1 = _call_at(gts, row1, markers[marker], layout)
        call2 = _call_at(gts, row2, markers[marker], layout)
        a0, a1, a2 = _log_step(a0, a1, a2, marker, call1, call2, transitions, emissions, False)
    if max(a0, a1, a2) == -np.inf:
        # The observations have zero likelihood under every possible path.
        output.fill(np.nan)
        return

    alpha = np.empty((_LOG_BLOCK, 3))
    b0, b1, b2 = 0., 0., 0.
    for block in range(len(checkpoints) - 1, -1, -1):
        start = block * _LOG_BLOCK
        end = min(start + _LOG_BLOCK, n)
        a0, a1, a2 = checkpoints[block, 0], checkpoints[block, 1], checkpoints[block, 2]
        for marker in range(start, end):
            call1 = _call_at(gts, row1, markers[marker], layout)
            call2 = _call_at(gts, row2, markers[marker], layout)
            a0, a1, a2 = _log_step(a0, a1, a2, marker, call1, call2, transitions, emissions, False)
            alpha[marker - start, 0], alpha[marker - start, 1], alpha[marker - start, 2] = a0, a1, a2
        for marker in range(end - 1, start - 1, -1):
            p0 = alpha[marker - start, 0] + b0
            p1 = alpha[marker - start, 1] + b1
            p2 = alpha[marker - start, 2] + b2
            peak = max(p0, p1, p2)
            p0, p1, p2 = _scale(np.exp(p0 - peak), np.exp(p1 - peak), np.exp(p2 - peak))
            output[marker, 0], output[marker, 1], output[marker, 2] = p0, p1, p2
            if marker:
                call1 = _call_at(gts, row1, markers[marker], layout)
                call2 = _call_at(gts, row2, markers[marker], layout)
                b0, b1, b2 = _log_step(b0, b1, b2, marker, call1, call2, transitions, emissions, True)


@njit(cache=True, parallel=True)
def _infer_prepared(sibpairs, gts, markers, transitions, emissions, layout,
                    posterior=None, rates=None, likelihoods=None):
    n = markers.size
    ibd = np.empty((len(sibpairs), n), dtype=np.int8)
    for i in prange(len(sibpairs)):
        row1, row2 = sibpairs[i, 0], sibpairs[i, 1]
        previous = np.log(np.array([.25, .5, .25]))
        current = np.empty(3)
        pointers = np.empty((n, 3), dtype=np.int8)
        call1 = _call_at(gts, row1, markers[0], layout)
        call2 = _call_at(gts, row2, markers[0], layout)
        if min(call1, call2) >= 0:
            previous += emissions[0, call1, call2]
        a0, a1, a2 = .25, .5, .25
        if posterior is not None:
            a0, a1, a2 = _condition(a0, a1, a2, 0, call1, call2, likelihoods)
            posterior[i, 0, 0], posterior[i, 0, 1], posterior[i, 0, 2] = a0, a1, a2
        for marker in range(1, n):
            call1 = _call_at(gts, row1, markers[marker], layout)
            call2 = _call_at(gts, row2, markers[marker], layout)
            observed = min(call1, call2) >= 0
            if posterior is not None:
                stay, move = rates[marker, 0], rates[marker, 1]
                half = .5 * move * a1
                a0, a1, a2 = stay * a0 + half, move * (a0 + a2) + stay * a1, half + stay * a2
                a0, a1, a2 = _condition(a0, a1, a2, marker, call1, call2, likelihoods)
                posterior[i, marker, 0], posterior[i, marker, 1], posterior[i, marker, 2] = a0, a1, a2
            for dest in range(3):
                best = 0
                maximum = previous[0] + transitions[marker, 0, dest]
                for source in range(1, 3):
                    probability = previous[source] + transitions[marker, source, dest]
                    if probability > maximum:
                        best, maximum = source, probability
                pointers[marker, dest] = best
                current[dest] = maximum
                if observed:
                    current[dest] += emissions[marker, call1, call2, dest]
            previous, current = current, previous
        state = np.argmax(previous)
        ibd[i, n - 1] = state
        for marker in range(n - 1, 0, -1):
            state = pointers[marker, state]
            ibd[i, marker - 1] = state
        if posterior is not None:
            output = posterior[i]
            if not _backward_posterior(output, gts, row1, row2, markers, rates, likelihoods, layout):
                _log_posterior(output, gts, row1, row2, markers, transitions, emissions, layout)
    return ibd


def infer_ibd_posteriors(sibpairs, gts, freqs, genetic_map, weights, error_probs):
    """Return marginal P(IBD = 0, 1, 2) for each pair and in-memory marker."""
    posterior = np.empty((len(sibpairs), gts.shape[1], 3), dtype=np.float32)
    if gts.shape[1]:
        transitions, emissions = _chromosome_hmm(freqs, genetic_map, weights, error_probs)
        rates, likelihoods = _posterior_tables(transitions, emissions)
        _infer_prepared(sibpairs, gts, np.arange(gts.shape[1]), transitions, emissions, 0,
                        posterior, rates, likelihoods)
    return posterior


@njit(cache=True, parallel=True)
def _pack_calls(values, out):
    for row in prange(values.shape[0]):
        for byte in range(out.shape[1]):
            word = 0
            for shift in range(min(4, values.shape[1] - 4 * byte)):
                call = values[row, 4 * byte + shift]
                code = int(call) if np.isfinite(call) and call >= 0 else 3
                word |= code << (2 * shift)
            out[row, byte] = word


class segment(object):
    def __init__(self,start_index,end_index,start_bp,end_bp,start_snp,end_snp,length,state):
        self.start = start_index
        self.end = end_index
        self.start_bp = start_bp
        self.end_bp = end_bp
        self.start_snp = start_snp
        self.end_snp = end_snp
        self.length = length
        self.state = state
    def to_text(self,id1,id2,chr,end=False):
        seg_txt = str(id1)+'\t'+str(id2)+'\t'+str(self.state)+'\t'+str(chr)+'\t'+str(self.start_bp)+'\t'+str(self.end_bp)+'\t'+str(self.start_snp)+'\t'+str(self.end_snp)+'\t'+str(self.length)
        if end:
            return seg_txt
        else:
            return seg_txt+'\n'

def find_segments(path, map, snps, pos):
    if not len(path):
        return []
    changes = np.flatnonzero(path[1:] != path[:-1]) + 1
    boundaries = np.concatenate(([0], changes, [len(path)]))
    segments = []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        last = end - 1
        segments.append(segment(start, end, pos[start], pos[last], snps[start], snps[last],
                                map[last] - map[start], path[start]))
    return segments

def smooth_segments(path, map, snps, pos, min_length):
    segments = find_segments(path, map, snps, pos)
    if len(segments) <= 1:
        return path, segments
    changed = False
    for i, current in enumerate(segments):
        if current.length >= min_length:
            continue
        if i == 0:
            replacement = segments[1].state
        elif i == len(segments) - 1:
            replacement = segments[i - 1].state
        elif segments[i - 1].state == segments[i + 1].state:
            replacement = segments[i - 1].state
        else:
            continue
        # Segment indices are half-open, physical coordinates remain inclusive.
        path[current.start:current.end] = replacement
        changed = True
    return path, find_segments(path, map, snps, pos) if changed else segments

def smooth_ibd(ibd,map,snps,pos,min_length):
    allsegs = []
    for i in range(ibd.shape[0]):
        ibd_path, segments = smooth_segments(ibd[i,:],map,snps,pos,min_length)
        ibd[i,:] = ibd_path
        allsegs.append(segments)
    return ibd, allsegs

_SEGMENT_HEADER = "ID1\tID2\tIBDType\tChr\tstart_coordinate\tstop_coordinate\tstartSNP\tstopSNP\tlength\n"


def write_segs(sibpairs, allsegs, chr, outfile):
    with gzip.open(outfile, "wt", compresslevel=3) as handle:
        handle.write(_SEGMENT_HEADER)
        for i, ((id1, id2), segments) in enumerate(zip(sibpairs, allsegs)):
            lines = [s.to_text(id1, id2, chr, end=True) for s in segments]
            if lines:
                handle.write("\n".join(lines))
                if i != len(sibpairs) - 1:
                    handle.write("\n")


def write_segs_from_matrix(ibd, sibpairs, snps, pos, map, chrom, outfile, *, return_segments=True):
    if return_segments:
        allsegs = [find_segments(path, map, snps, pos) for path in ibd]
        write_segs(sibpairs, allsegs, chrom, outfile)
        return allsegs
    # Simulated chromosomes can contain millions of short segments. Stream
    # one sibling pair at a time for callers that only need the file.
    snps, positions = list(snps), [str(p) for p in pos]
    with gzip.open(outfile, "wt", compresslevel=3) as handle:
        handle.write(_SEGMENT_HEADER)
        for i, ((id1, id2), path) in enumerate(zip(sibpairs, ibd)):
            if not len(path):
                continue
            changes = np.flatnonzero(path[1:] != path[:-1]) + 1
            starts = np.concatenate(([0], changes))
            ends = np.concatenate((changes - 1, [len(path) - 1]))
            states = path[starts].tolist()
            lengths = (map[ends] - map[starts]).tolist()
            prefix = f"{id1}\t{id2}\t"
            lines = [
                f"{prefix}{state}\t{chrom}\t{positions[a]}\t{positions[b]}\t{snps[a]}\t{snps[b]}\t{length}"
                for a, b, state, length in zip(starts.tolist(), ends.tolist(), states, lengths)
            ]
            handle.write("\n".join(lines))
            if i != len(sibpairs) - 1:
                handle.write("\n")

def _workspace(stack, scratch, name, shape, dtype, memory, order='C'):
    if np.prod(shape, dtype=np.int64) * np.dtype(dtype).itemsize <= memory:
        return np.empty(shape, dtype=dtype, order=order)
    result = np.memmap(scratch / name, mode='w+', dtype=dtype, shape=shape, order=order)
    stack.callback(result._mmap.close)
    return result


def _ibd_map(source, chrom, positions, mapfile, *, ordered_bim=True):
    if mapfile is not None:
        return get_map_positions(mapfile, SimpleNamespace(pos=positions), min_map_prop=None, extrapolate=True)
    if source.bed:
        # Keep the map in float64. The BED reader exposes a float32 cM
        # field, which can move a pair across the exact LD window boundary.
        genetic_map = np.loadtxt(Path(source.reader.location).with_suffix('.bim'),
                                 dtype=np.float64, usecols=2, ndmin=1)
        if np.any(genetic_map != 0) and (len(genetic_map) == 1 or np.any(genetic_map != genetic_map[0])):
            if not np.isfinite(genetic_map).all() or np.min(genetic_map) < 0 or np.max(genetic_map) > 5000:
                raise ValueError('BIM genetic map must be finite and within [0, 5000] cM')
            if ordered_bim and np.any(np.diff(genetic_map) < 0):
                raise ValueError('BIM genetic map must be nondecreasing')
            return genetic_map
    return decode_map_from_pos(chrom, positions, extrapolate=True)


def _joint_families(row_pairs, eligible=None):
    """Identify complete connected components of three or four siblings."""
    from snipar.ibd_joint import joint_pair_groups

    labels = np.full(len(row_pairs), -1, dtype=np.int64)
    groups = []
    rows = np.flatnonzero(eligible) if eligible is not None else None
    if rows is not None and not len(rows):
        return labels, groups
    for edges in joint_pair_groups(row_pairs if rows is None else row_pairs[rows]):
        if rows is not None:
            edges = rows[edges]
        labels[edges] = len(groups)
        groups.append(np.unique(row_pairs[edges]))
    return labels, groups


@njit(cache=True)
def _family_calls(genotypes, members, markers, layout):
    calls = np.empty((len(members), len(markers)), dtype=np.int8)
    for member in range(len(members)):
        for marker in range(len(markers)):
            calls[member, marker] = _call_at(genotypes, members[member], markers[marker], layout)
    return calls


def _infer_block(row_pairs, genotypes, markers, transitions, emissions, layout, tables):
    if tables is None:
        return _infer_prepared(row_pairs, genotypes, markers, transitions, emissions, layout), None
    posterior = np.empty((len(row_pairs), len(markers), 3), dtype=np.float32)
    paths = _infer_prepared(row_pairs, genotypes, markers, transitions, emissions, layout, posterior, *tables)
    return paths, posterior


def _infer_joint_block(row_pairs, labels, groups, genotypes, markers, transitions, emissions,
                       layout, tables, freqs, genetic_map, weights, error_probs):
    from snipar.ibd_joint import infer_joint_family

    paths = np.empty((len(row_pairs), len(markers)), dtype=np.int8)
    posterior = None if tables is None else np.empty((*paths.shape, 3), dtype=np.float32)
    fallback = np.flatnonzero(labels < 0)
    if fallback.size:
        pair_paths, pair_posterior = _infer_block(row_pairs[fallback], genotypes, markers,
                                                  transitions, emissions, layout, tables)
        paths[fallback] = pair_paths
        if posterior is not None:
            posterior[fallback] = pair_posterior
        del pair_paths, pair_posterior
    joint_rows = np.flatnonzero(labels >= 0)
    joint_rows = joint_rows[np.argsort(labels[joint_rows], kind='stable')]
    boundaries = np.flatnonzero(labels[joint_rows[:-1]] != labels[joint_rows[1:]]) + 1
    for in_group in np.split(joint_rows, boundaries):
        label = labels[in_group[0]]
        members = groups[label]
        calls = _family_calls(genotypes, members, markers, layout)
        local_pairs = np.searchsorted(members, row_pairs[in_group])
        group_paths, group_posterior = infer_joint_family(
            calls, local_pairs, freqs, genetic_map, weights, error_probs, posterior=posterior is not None,
        )
        paths[in_group] = group_paths
        if posterior is not None:
            posterior[in_group] = group_posterior
        del calls, group_paths, group_posterior
    return paths, posterior


def _write_ibd_batches(pairs, row_pairs, genotypes, layout, markers, transitions, emissions,
                       snps, positions, genetic_map, min_length, chrom, outprefix, ibdmatrix,
                       pair_batch, pair_indices, posterior_callback, tables, joint_families,
                       freqs, weights, error_probs, phased_inference=None):
    labels, groups = joint_families
    if groups:
        print(f'Using joint IBD inference for {len(groups)} complete three- or four-sibling families')
    with ExitStack() as stack:
        output = stack.enter_context(gzip.open(outfile_name(outprefix, '.ibd.segments.gz', chrom),
                                               'wt', compresslevel=3))
        output.write(_SEGMENT_HEADER)
        matrix = None
        if ibdmatrix:
            matrix = stack.enter_context(gzip.open(outfile_name(outprefix, '.ibdmatrix.gz', chrom),
                                                   'wt', compresslevel=3))
            matrix.write('sib1 sib2 ' + ' '.join(snps) + '\n')
        for start in range(0, len(pairs), pair_batch):
            end = min(start + pair_batch, len(pairs))
            joint_labels = labels[start:end] if groups else None
            if phased_inference is not None:
                paths, posterior = phased_inference(row_pairs[start:end])
            elif joint_labels is not None and np.any(joint_labels >= 0):
                paths, posterior = _infer_joint_block(
                    row_pairs[start:end], joint_labels, groups, genotypes, markers, transitions, emissions,
                    layout, tables, freqs, genetic_map, weights, error_probs,
                )
            else:
                paths, posterior = _infer_block(row_pairs[start:end], genotypes, markers,
                                                transitions, emissions, layout, tables)
            if posterior_callback is not None:
                posterior_callback(chrom, pairs[start:end], snps, positions, genetic_map, posterior,
                                   pair_indices[start:end])
            for i, ((id1, id2), path) in enumerate(zip(pairs[start:end], paths)):
                if joint_labels is not None and joint_labels[i] >= 0:
                    segments = find_segments(path, genetic_map, snps, positions)
                else:
                    path, segments = smooth_segments(path, genetic_map, snps, positions, min_length)
                output.writelines(seg.to_text(id1, id2, chrom) for seg in segments)
                if matrix is not None:
                    matrix.write(f'{id1} {id2} ' + ' '.join(path.astype(str)) + '\n')
            del paths, posterior


class IBDChromosomeFiles(NamedTuple):
    chromosome: int | str
    stem: str | None
    origin: str | None


def _input_chromosome(source, chrom):
    if chrom is not None:
        return chrom
    chromosomes = np.unique(source.reader.chromosome if source.bed else source.reader.chromosomes)
    if len(chromosomes) != 1:
        raise ValueError('Provide observed genotypes for exactly one chromosome')
    chrom = chromosomes[0] or 0
    label = re.fullmatch(r'(?:chr)?([0-9]+)(?:\.0+)?', str(chrom), flags=re.IGNORECASE)
    return int(label[1]) if label is not None else chrom


def _save_origin_sidecar(source, chrom, mapfile, outprefix, duos, positions, packed=None, read_rows=None):
    from snipar.ibd_prephase import origin_grid, pack_phased, read_phased_calls, write_parent_origins

    markers, ids, alleles = origin_grid(source.reader, positions)
    origin_pos = positions[markers]
    origin_map = _ibd_map(source, chrom, origin_pos, mapfile)
    path = outfile_name(outprefix, '.ibd.origin.hdf5', chrom)
    print(f'Inferring parent of origin for {len(duos)} parent-offspring pairs')
    with ExitStack() as stack:
        packed_markers = markers
        if packed is None:
            read_rows = np.unique(duos)
            memory = get_memory_budget()
            m = len(markers)
            shape = (len(read_rows), (m + 1) // 2)
            packed_bytes = int(np.prod(shape, dtype=np.int64))
            in_memory = packed_bytes <= memory // 5
            scratch = Path(stack.enter_context(TemporaryDirectory(prefix='snipar-origin-')))
            if not in_memory and packed_bytes > disk_usage(scratch).free:
                raise OSError(f'Origin inference scratch needs {packed_bytes:,} bytes')
            packed = _workspace(stack, scratch, 'genotypes', shape, np.uint8, memory // 5)
            fixed = source.fixed_bytes + (packed_bytes if in_memory else 0)
            _, width = plan_batches(m, 40 * len(read_rows), fixed_bytes=fixed, memory_limit=memory)
            if width < m:
                if width < 2:
                    raise MemoryError('Phased packing needs enough memory to decode two variants')
                width -= width % 2
            for start in range(0, m, width):
                stop = min(m, start + width)
                calls = read_phased_calls(source, read_rows, markers[start:stop])
                pack_phased(calls, packed[:, start // 2:(stop + 1) // 2])
                del calls
            packed_markers = np.arange(m, dtype=np.intp)
        local_duos = np.searchsorted(read_rows, duos)
        write_parent_origins(path, chrom, local_duos, source.ids[read_rows], packed, packed_markers, origin_pos,
                             ids, alleles, origin_map, source.threads)
    return path, markers, local_duos


def infer_ibd_chr(sibpairs, error_prob, error_probs, outprefix, bedfile=None, bgenfile=None, chrom=None,
                  min_length=0.01, mapfile=None, ibdmatrix=False, ld_out=False, min_maf=0.01,
                  max_missing=5, max_error=0.01, num_threads=None, *, posterior_callback=None, joint=False,
                  joint_pair_mask=None, phased_ibd=None, phase_switch_rate=.002, pedigree=None,
                  all_parental_origin=False, allow_no_siblings=False, origin_only=False):
    """Decode one chromosome in bounded pair batches and optionally report marginal probabilities.

    ``posterior_callback`` receives chromosome, pair IDs, marker IDs, physical
    positions, genetic positions, float32 P(IBD0, IBD1, IBD2), and the pair row
    indices in the original ``sibpairs`` array. Consume the probabilities one
    batch at a time to maintain bounded memory use.

    An optional ``joint_pair_mask`` identifies pairs eligible for the joint
    model before restricting the input to a batch of the full sibling graph.
    By default, phased inputs use the phased pair model. A pedigree enables
    duoHMM for children with one or two genotyped parents. The original-slot
    calls are written for imputation or parent-of-origin haplotype export.
    """
    threads = get_num_threads() if num_threads is None else num_threads
    if phased_ibd and joint:
        raise ValueError('Phased and joint unphased IBD are separate models')
    # Marker decoding and Numba computation execute in separate phases.
    with ExitStack() as stack:
        source = stack.enter_context(open_observed(bedfile, bgenfile, num_threads=threads))
        sids = source.sid
        if not len(sids):
            raise ValueError('No observed variants on chromosome')
        has_phase = not source.bed and source.phased
        if phased_ibd and not has_phase:
            raise ValueError('Phased IBD requires phased BGEN or phased VCF genotypes')
        phased = has_phase and (phased_ibd is None or phased_ibd) and not joint and not origin_only
        make_origins = has_phase and pedigree is not None and (phased or all_parental_origin or origin_only)
        if phased and not 0 <= phase_switch_rate <= .5:
            raise ValueError('Phase switch probability must be between 0 and 0.5')
        reader = source.reader
        positions = reader.bp_position if source.bed else np.asarray(reader.positions)
        chrom = _input_chromosome(source, chrom)
        if mapfile is not None:
            mapfile = str(mapfile).replace('@', str(chrom))
        sample_index = IdIndex(source.ids)
        row_pairs = sample_index.lookup(sibpairs)
        present = (row_pairs >= 0).all(axis=1)
        has_pairs = np.any(present) and not origin_only
        if not has_pairs and not (allow_no_siblings or origin_only):
            raise ValueError('No genotyped sibling pairs found')
        origin_path = origin_markers = None
        duos = np.empty((0, 2), dtype=np.intp)
        duo_roles = np.empty(0, dtype=np.intp)
        if make_origins:
            from snipar.ibd_prephase import observed_parent_duos

            children = None if all_parental_origin or origin_only else np.unique(row_pairs[present])
            duos, duo_roles = observed_parent_duos(pedigree, source.ids, children, sample_index=sample_index)
        del sample_index
        if len(duos) and (not has_pairs or not phased):
            origin_path, _, _ = _save_origin_sidecar(source, chrom, mapfile, outprefix, duos, positions)
        if not has_pairs:
            return IBDChromosomeFiles(chrom, None, origin_path)
        pair_indices = np.flatnonzero(present) if posterior_callback is not None else None
        if joint and joint_pair_mask is not None:
            joint_pair_mask = np.asarray(joint_pair_mask, dtype=bool)
            if joint_pair_mask.shape != present.shape:
                raise ValueError('Joint IBD eligibility must have one entry per input sibling pair')
            joint_pair_mask = joint_pair_mask[present]
        pairs, row_pairs = sibpairs[present], row_pairs[present]
        rows = np.unique(row_pairs)
        n, m = len(rows), len(sids)
        read_rows = rows
        packing = 4
        if phased:
            from snipar.ibd_prephase import PhasedPairInference, pack_phased, read_phased_calls

            read_rows = np.union1d(rows, duos)
            sibling_rows = np.searchsorted(read_rows, rows)
            packing = 2
            print('Using phased sibling IBD')
        print(f'Inferring IBD for {len(pairs)} sibling pairs on chromosome {chrom}')
        ld_rows = np.random.randint(0, n, size=min(n, 1000))
        memory = get_memory_budget()
        current_dir = Path(stack.enter_context(TemporaryDirectory(prefix='snipar-ibd-')))
        packed_shape = (len(read_rows), (m + packing - 1) // packing)
        ld_shape = (m, 3, (len(ld_rows) + 63) // 64)
        packed_bytes = np.prod(packed_shape, dtype=np.int64)
        ld_bytes = np.prod(ld_shape, dtype=np.int64) * np.dtype(np.uint64).itemsize
        in_memory_packed = packed_bytes <= memory // 5
        in_memory_ld = ld_bytes <= memory // 5
        needed = (0 if in_memory_packed else packed_bytes) + (0 if in_memory_ld else ld_bytes)
        layout = 1
        available_disk = disk_usage(current_dir).free
        if needed > available_disk:
            if source.bed and (in_memory_ld or ld_bytes < available_disk):
                # BED is already packed. This fallback rereads its variant-major
                # pages during tracing but does not require extra genotype disk.
                with open(reader.location, 'rb') as stream:
                    if stream.read(3) != b'\x6c\x1b\x01':
                        raise ValueError('IBD requires a variant-major PLINK BED file')
                genotypes = np.memmap(reader.location, mode='r', dtype=np.uint8, offset=3,
                                      shape=(m, (len(source.ids) + 3) // 4))
                stack.callback(genotypes._mmap.close)
                layout = 2
            else:
                raise OSError(f'IBD scratch needs {needed:,} bytes, with {available_disk:,} bytes available')
        if layout == 1:
            genotypes = _workspace(stack, current_dir, 'genotypes', packed_shape, np.uint8, memory // 5)
            row_pairs = np.searchsorted(read_rows, row_pairs)
        ld_genotypes = _workspace(stack, current_dir, 'ld', ld_shape, np.uint64, memory // 5)
        fixed = source.fixed_bytes + (packed_bytes if layout == 1 and in_memory_packed else 0)
        fixed += ld_bytes if in_memory_ld else 0
        read_bytes = 40 * len(read_rows) if phased else n * (2 if source.bed else 24)
        _, variant_batch = plan_batches(m, read_bytes, fixed_bytes=fixed,
                                        memory_limit=memory)
        if layout == 1 and variant_batch < m:
            if variant_batch < packing:
                raise MemoryError(f'IBD packing needs enough memory to decode {packing} variants')
            variant_batch -= variant_batch % packing
        freqs = np.empty(m)
        missing = np.empty(m)
        for start in range(0, m, variant_batch):
            end = min(start + variant_batch, m)
            columns = np.arange(start, end)
            if phased:
                calls = read_phased_calls(source, read_rows, columns)
                pack_phased(calls, genotypes[:, start // 2:(end + 1) // 2])
                # Each called allele is 0 or 1. Any sum involving a missing
                # allele stays negative and is excluded by the dosage kernels.
                values = np.add(calls[:, :, 0], calls[:, :, 1])[sibling_rows]
                del calls
            else:
                values = source.read(rows, columns, require_calls=True)
            total, called = allele_stats(values)
            freqs[start:end] = np.where(called > 0, total / np.maximum(2 * called, 1), np.nan)
            missing[start:end] = 100 * (n - called) / n
            pack_ld_hardcalls(values[ld_rows, :], ld_genotypes[start:end])
            if layout == 1 and not phased:
                _pack_calls(values, genotypes[:, start // 4:(end + 3) // 4])
            del values
        p_error = np.full(m, error_prob, dtype=float)
        if error_probs is not None:
            index = IdIndex(error_probs.sid).lookup(sids)
            valid = index >= 0
            p_error[valid] = error_probs.error_ests[index[valid]]
        markers = np.flatnonzero((freqs > min_maf) & (freqs < 1 - min_maf)
                                 & (missing <= max_missing) & (p_error <= max_error))
        if not len(markers):
            raise ValueError('No SNPs remain after MAF, missingness and genotyping error filtering')
        marker_positions = positions[markers]
        order_by_position = posterior_callback is not None or joint or phased
        if order_by_position and np.any(marker_positions[1:] < marker_positions[:-1]):
            markers = markers[np.argsort(marker_positions, kind='stable')]
        if mapfile is None and source.bed:
            genetic_map = _ibd_map(source, chrom, positions, None, ordered_bim=not order_by_position)[markers]
        else:
            genetic_map = _ibd_map(source, chrom, positions[markers], mapfile)
        keep = np.isfinite(genetic_map)
        markers, genetic_map = markers[keep], genetic_map[keep]
        if not len(markers) or np.any(np.diff(genetic_map) < 0):
            raise ValueError('SNPs must have increasing genetic map positions')
        print('Computing LD weights')
        # Keep the original variant-major LD cache on disk and pass the retained
        # source column indices to the compiled kernel.
        ld = compute_ld_scores_packed(ld_genotypes, genetic_map, 1., markers)
        weights = 1 / np.where(np.isfinite(ld) & (ld > 0), ld, 1)
        snps, pos = sids[markers], positions[markers]
        marker_freqs, marker_errors = freqs[markers], p_error[markers]
        if phased:
            if len(duos):
                origin_path, origin_markers, duos = _save_origin_sidecar(
                    source, chrom, mapfile, outprefix, duos, positions, genotypes, read_rows,
                )
            phased_inference = PhasedPairInference(
                stack, genotypes, markers, marker_freqs, genetic_map, weights, marker_errors, phase_switch_rate,
                len(pairs), threads, duos=duos, duo_roles=duo_roles, origin_path=origin_path,
                origin_markers=origin_markers,
                posterior=posterior_callback is not None,
            )
            pair_batch = phased_inference.batch_size
            transitions = emissions = tables = None
        else:
            phased_inference = None
            transitions, emissions = _chromosome_hmm(marker_freqs, genetic_map, weights, marker_errors)
            tables = _posterior_tables(transitions, emissions) if posterior_callback is not None else None
        joint_families = _joint_families(row_pairs, joint_pair_mask) if joint else (None, [])
        if not phased:
            nmarkers = len(markers)
            output_bytes = 1 if tables is None else 13
            pair_bytes = nmarkers * output_bytes + 1
            worker_bytes = 3 * nmarkers + 128
            if tables is not None:
                worker_bytes += 24 * ((nmarkers + _LOG_BLOCK - 1) // _LOG_BLOCK + _LOG_BLOCK)
            fixed_bytes = min(get_num_threads(), len(pairs)) * worker_bytes
            # Projection needs an extra batch when inference mixes the two models.
            if joint_families[1]:
                pair_bytes *= 2
                fixed_bytes += (4 + 6 * output_bytes) * nmarkers
            batch_cap = 10000 if tables is None else max(1, min(10000, 64 * 1024**2 // pair_bytes))
            _, pair_batch = plan_batches(len(pairs), pair_bytes, batch_size=batch_cap, fixed_bytes=fixed_bytes)
        print(f'Tracing IBD in blocks of up to {pair_batch} sibling pairs')
        _write_ibd_batches(pairs, row_pairs, genotypes, layout, markers, transitions, emissions,
                           snps, pos, genetic_map, min_length, chrom, outprefix, ibdmatrix,
                           pair_batch, pair_indices, posterior_callback, tables, joint_families,
                           marker_freqs, weights, marker_errors, phased_inference)
        if ld_out:
            filename = outfile_name(outprefix, '.l2.ldscore.gz', chrom)
            with gzip.open(filename, 'wt', compresslevel=3) as output:
                output.write('CHR SNP BP L2\n')
                output.writelines(f'{chrom} {sid} {bp} {score}\n' for sid, bp, score in zip(snps, pos, ld))
        return IBDChromosomeFiles(chrom, outfile_name(outprefix, '.ibd', chrom), origin_path)
