"""Phase-aware IBD inference for a full-sibling pair.

State (C, U0, U1) records whether the siblings' input-haplotype labels align
directly (C=0) or are crossed (C=1), and whether each of sibling one's two
haplotypes is shared. The eight states remove arbitrary absolute parental
labels under a common recombination map. The two IBD0 states have identical
emissions but retain their relative alignments across an IBD0 interval.

This model requires user-specified phase switch rates and allele error rates.
It assumes unrelated founders and single-marker allele-frequency emissions;
LD weights are optional. In IBD2, equal switch rates make the source of a
relative orientation change unidentifiable from the pair alone.
"""

from typing import NamedTuple

from numba import njit
import numpy as np


# State code is 4*C + 2*U0 + U1.
_IBD_COUNT = np.array([0, 1, 1, 2, 0, 1, 1, 2], dtype=np.int8)


class PhasedPairResult(NamedTuple):
    state_path: np.ndarray
    hard_ibd: np.ndarray
    state_posterior: np.ndarray
    ibd_posterior: np.ndarray
    phase_switch_posterior: np.ndarray


@njit(cache=True)
def _transition(distance, s1, s2, matrix, first_switch, second_switch):
    """Sum independent meiosis and switch parities for every source-target pair."""
    matrix[:, :] = 0
    first_switch[:, :] = 0
    second_switch[:, :] = 0
    # Two meioses contribute to the sharing parity for each parent.
    q = -.5 * np.expm1(-distance / 25)
    for source in range(8):
        c, u0, u1 = source >> 2, (source >> 1) & 1, source & 1
        for flip_first in range(2):
            a = u1 if flip_first else u0
            b = u0 if flip_first else u1
            p_first = s1 if flip_first else 1 - s1
            for flip_second in range(2):
                target_c = c ^ flip_first ^ flip_second
                p_second = s2 if flip_second else 1 - s2
                for change0 in range(2):
                    p_change0 = q if change0 else 1 - q
                    for change1 in range(2):
                        p_change1 = q if change1 else 1 - q
                        target = 4 * target_c + 2 * (a ^ change0) + (b ^ change1)
                        probability = p_first * p_second * p_change0 * p_change1
                        matrix[source, target] += probability
                        if flip_first:
                            first_switch[source, target] += probability
                        if flip_second:
                            second_switch[source, target] += probability


@njit(cache=True, inline="always")
def _allele_likelihood(observed, allele, error):
    if observed < 0:
        return 1.
    return 1 - error if observed == allele else error


@njit(cache=True)
def _emissions(phase, f, errors, weight, likelihood, log_likelihood):
    if weight == 0:
        likelihood[:] = 1
        log_likelihood[:] = 0
        return True
    l1, l2 = np.empty((2, 2)), np.empty((2, 2))
    for hap in range(2):
        for allele in range(2):
            l1[hap, allele] = _allele_likelihood(phase[0, hap], allele, errors[0])
            l2[hap, allele] = _allele_likelihood(phase[1, hap], allele, errors[1])
    q = 1 - f
    maximum = 0.
    for c in range(2):
        shared, independent = np.empty(2), np.empty(2)
        for hap in range(2):
            other = hap ^ c
            shared[hap] = q * l1[hap, 0] * l2[other, 0] + f * l1[hap, 1] * l2[other, 1]
            independent[hap] = ((q * l1[hap, 0] + f * l1[hap, 1]) *
                                (q * l2[other, 0] + f * l2[other, 1]))
        for u0 in range(2):
            for u1 in range(2):
                state = 4 * c + 2 * u0 + u1
                p0 = shared[0] if u0 else independent[0]
                p1 = shared[1] if u1 else independent[1]
                likelihood[state] = p0 * p1
                maximum = max(maximum, likelihood[state])
    if maximum == 0:
        return False
    normalizer = np.log(maximum)
    for state in range(8):
        probability = likelihood[state]
        log_p = weight * (np.log(probability) - normalizer) if probability > 0 else -np.inf
        log_likelihood[state] = log_p
        likelihood[state] = np.exp(log_p)
    return True


@njit(cache=True)
def _infer_haplotypes(haplotypes, freqs, genetic_map, switches, errors, weights,
                      state_path, posterior, switch_posterior):
    markers = len(freqs)
    alpha = np.empty((markers, 8))
    pointer = np.empty((markers, 8), np.int8)
    transitions = np.empty((3, 8, 8))
    matrix, first_switch, second_switch = transitions
    log_previous = np.zeros(8)
    work = np.empty((4, 8))
    log_current, likelihood, log_likelihood, current = work
    for t in range(markers):
        if not _emissions(haplotypes[:, t, :], freqs[t], errors[:, t], weights[t], likelihood, log_likelihood):
            return t
        if t == 0:
            current[:] = 1. / 8
            log_current[:] = log_previous
        else:
            _transition(genetic_map[t] - genetic_map[t - 1], switches[0, t - 1], switches[1, t - 1],
                        matrix, first_switch, second_switch)
            for dest in range(8):
                maximum, best, total = -np.inf, 0, 0.
                for source in range(8):
                    probability = matrix[source, dest]
                    if probability > 0:
                        candidate = log_previous[source] + np.log(probability)
                        if candidate > maximum:
                            best, maximum = source, candidate
                        total += alpha[t - 1, source] * probability
                pointer[t, dest] = best
                log_current[dest] = maximum
                current[dest] = total
        log_current += log_likelihood
        maximum = np.max(log_current)
        current *= likelihood
        scale = current.sum()
        if not np.isfinite(maximum) or scale <= 0:
            return t
        log_previous[:] = log_current - maximum
        alpha[t] = current / scale

    state = np.argmax(log_previous)
    beta = np.ones(8)
    for t in range(markers - 1, -1, -1):
        state_path[t] = state
        current[:] = alpha[t] * beta
        scale = current.sum()
        if scale <= 0:
            return t
        posterior[t] = current / scale
        if not t:
            break

        state = pointer[t, state]
        _emissions(haplotypes[:, t, :], freqs[t], errors[:, t], weights[t], likelihood, log_likelihood)
        likelihood *= beta
        _transition(genetic_map[t] - genetic_map[t - 1], switches[0, t - 1], switches[1, t - 1],
                    matrix, first_switch, second_switch)
        both_weight = first_weight = second_weight = 0.
        for source in range(8):
            total = first = second = 0.
            for dest in range(8):
                p = likelihood[dest]
                total += matrix[source, dest] * p
                first += first_switch[source, dest] * p
                second += second_switch[source, dest] * p
            beta[source] = total
            both_weight += alpha[t - 1, source] * total
            first_weight += alpha[t - 1, source] * first
            second_weight += alpha[t - 1, source] * second
        if both_weight <= 0:
            return t
        switch_posterior[0, t - 1] = first_weight / both_weight
        switch_posterior[1, t - 1] = second_weight / both_weight
        beta /= beta.sum()
    return -1


def infer_phased_pair(haplotypes, freqs, genetic_map, phase_switch, allele_error, *, weights=None):
    """Decode a pair from biallelic phased hard calls of shape (2, markers, 2).

    Calls are 0 or 1, with -1 and NaN representing missing alleles. Genetic
    coordinates are nondecreasing cM positions. ``phase_switch`` broadcasts
    to (2, markers - 1) and ``allele_error`` broadcasts to (2, markers); the
    first axis selects the sibling. Supply 0.5 for a switch rate at a boundary
    where the next input phase block has arbitrary orientation.

    The joint state path uses code 4*C + 2*U0 + U1. Posterior array shapes are
    (markers, 8), (markers, 3) and (2, markers - 1), respectively. Switch
    posteriors describe toggles of input phase labels, including uncertain
    phase-block boundaries. Parental identity is arbitrary.
    """
    haps = np.asarray(haplotypes)
    if haps.ndim != 3 or haps.shape[0] != 2 or haps.shape[2] != 2 or not haps.shape[1]:
        raise ValueError("Provide phased alleles of shape (2, positive marker count, 2)")
    if haps.dtype.kind not in 'biuf' or not np.all(np.isnan(haps) | np.isin(haps, (-1, 0, 1))):
        raise ValueError("Phased alleles must be 0, 1, -1 or NaN")
    if haps.dtype.kind == 'f':
        haps = np.where(np.isnan(haps), -1, haps)
    haps = np.ascontiguousarray(haps, dtype=np.int8)
    markers = haps.shape[1]
    freqs = np.asarray(freqs, dtype=float)
    genetic_map = np.asarray(genetic_map, dtype=float)
    weights = np.ones(markers) if weights is None else np.asarray(weights, dtype=float)
    if any(x.shape != (markers,) or not np.all(np.isfinite(x)) for x in (freqs, genetic_map, weights)):
        raise ValueError("Marker vectors must be finite and match the phased alleles")
    if np.any(np.diff(genetic_map) < 0):
        raise ValueError("Genetic coordinates must be nondecreasing")
    if np.any(weights < 0):
        raise ValueError("Marker weights must be nonnegative")
    try:
        phase_switch = np.broadcast_to(np.asarray(phase_switch, dtype=float), (2, markers - 1))
        allele_error = np.broadcast_to(np.asarray(allele_error, dtype=float), (2, markers))
    except ValueError as exc:
        raise ValueError("Switch and allele error rates must match the interval and marker axes") from exc
    if any(not np.all(np.isfinite(x) & (x >= 0) & (x <= 1)) for x in (freqs, phase_switch, allele_error)):
        raise ValueError("Allele frequencies, phase switch and allele error probabilities must be in [0, 1]")

    path = np.empty(markers, np.int8)
    posterior = np.empty((markers, 8), np.float32)
    switches = np.empty((2, markers - 1), np.float32)
    failed = _infer_haplotypes(haps, freqs, genetic_map, phase_switch, allele_error, weights, path, posterior, switches)
    if failed >= 0:
        raise ValueError(f"No phased IBD path through marker {failed + 1}; check error rates and allele frequencies")
    ibd_posterior = np.column_stack([posterior[:, _IBD_COUNT == k].sum(axis=1) for k in range(3)])
    return PhasedPairResult(path, _IBD_COUNT[path], posterior, ibd_posterior, switches)
