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

Optional duoHMM labels identify the homologs transmitted from a common
observed parent. When labels are informative, the model tracks both absolute
child orientations and projects its 16 states back to the eight pair states.
"""

from typing import NamedTuple

from numba import njit
import numpy as np


# State code is 4*C + 2*U0 + U1.
_IBD_COUNT = np.array([0, 1, 1, 2, 0, 1, 1, 2], dtype=np.int8)
# State code is 8*O1 + 4*O2 + 2*K + U. O1 and O2 index each child's
# known-parent haplotype; K and U indicate sharing from known and other parent.
_PAIR_STATES = np.array([0, 1, 2, 3, 4, 5, 6, 7, 4, 6, 5, 7, 0, 2, 1, 3], dtype=np.int8)


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
    if first_switch is not None:
        first_switch[:, :] = 0
    if second_switch is not None:
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
                        if flip_first and first_switch is not None:
                            first_switch[source, target] += probability
                        if flip_second and second_switch is not None:
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
def _labeled_emissions(phase, f, errors, weight, origins, copies, origin_errors, copy_errors,
                       log_likelihood, pair_likelihood, pair_log_likelihood):
    if not _emissions(phase, f, errors, weight, pair_likelihood, pair_log_likelihood):
        return False
    first, second = origins
    a, b = origin_errors
    first_match, first_mismatch = np.log1p(-a), np.log(a)
    second_match, second_mismatch = np.log1p(-b), np.log(b)
    both_copies_known = copies[0] >= 0 and copies[1] >= 0
    same_copy = copies[0] == copies[1]
    # The reported equality is wrong exactly when one of the copy labels is wrong.
    a, b = copy_errors
    parity_error = a * (1 - b) + b * (1 - a)
    match, mismatch = np.log1p(-parity_error), np.log(parity_error)
    for state, projected in enumerate(_PAIR_STATES):
        log_p = pair_log_likelihood[projected]
        if first >= 0:
            log_p += first_match if (state >> 3) == first else first_mismatch
        if second >= 0:
            log_p += second_match if ((state >> 2) & 1) == second else second_mismatch
        if both_copies_known:
            log_p += match if bool(state & 2) == same_copy else mismatch
        log_likelihood[state] = log_p
    log_likelihood -= np.max(log_likelihood)
    return True


@njit(cache=True, inline="always")
def _store_posterior(probability, scale, output):
    if len(output) == 8:
        output[:] = probability / scale
    else:
        output[:] = 0
        for state in range(8):
            output[_IBD_COUNT[state]] += np.float32(probability[state] / scale)


@njit(cache=True, nogil=True)
def _infer_haplotypes(haplotypes, freqs, genetic_map, switches, errors, weights,
                      state_path, posterior, switch_posterior):
    markers = len(freqs)
    alpha = np.empty((markers if posterior is not None else 0, 8))
    pointer = np.empty((markers, 8), np.int8)
    matrix = np.empty((8, 8))
    first_switch = np.empty((8, 8)) if switch_posterior is not None else None
    second_switch = np.empty((8, 8)) if switch_posterior is not None else None
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
                        if posterior is not None:
                            total += alpha[t - 1, source] * probability
                pointer[t, dest] = best
                log_current[dest] = maximum
                current[dest] = total
        log_current += log_likelihood
        maximum = np.max(log_current)
        if not np.isfinite(maximum):
            return t
        if posterior is not None:
            current *= likelihood
            scale = current.sum()
            if scale <= 0:
                return t
            alpha[t] = current / scale
        log_previous[:] = log_current - maximum

    state = np.argmax(log_previous)
    if posterior is None:
        for t in range(markers - 1, 0, -1):
            state_path[t] = state
            state = pointer[t, state]
        state_path[0] = state
        return -1

    beta = np.ones(8)
    for t in range(markers - 1, -1, -1):
        state_path[t] = state
        current[:] = alpha[t] * beta
        scale = current.sum()
        if scale <= 0:
            return t
        _store_posterior(current, scale, posterior[t])
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
                if first_switch is not None and second_switch is not None:
                    first += first_switch[source, dest] * p
                    second += second_switch[source, dest] * p
            beta[source] = total
            both_weight += alpha[t - 1, source] * total
            first_weight += alpha[t - 1, source] * first
            second_weight += alpha[t - 1, source] * second
        if both_weight <= 0:
            return t
        if switch_posterior is not None:
            switch_posterior[0, t - 1] = first_weight / both_weight
            switch_posterior[1, t - 1] = second_weight / both_weight
        beta /= beta.sum()
    return -1


@njit(cache=True, inline="always")
def _log_add(first, second):
    if second > first:
        first, second = second, first
    if second == -np.inf:
        return first
    return first + np.log1p(np.exp(second - first))


@njit(cache=True, inline="always")
def _mix_sum(values, bit, stay, toggle):
    for lower in range(0, 16, 2 * bit):
        for state in range(lower, lower + bit):
            other = state + bit
            first, second = values[state], values[other]
            values[state] = _log_add(first + stay, second + toggle)
            values[other] = _log_add(second + stay, first + toggle)


@njit(cache=True, inline="always")
def _mix_max(values, sources, bit, stay, toggle):
    for lower in range(0, 16, 2 * bit):
        for state in range(lower, lower + bit):
            other = state + bit
            first, second = values[state], values[other]
            first_source, second_source = sources[state], sources[other]
            direct, flipped = first + stay, second + toggle
            if flipped > direct or (flipped == direct and second_source < first_source):
                values[state], sources[state] = flipped, second_source
            else:
                values[state], sources[state] = direct, first_source
            direct, flipped = second + stay, first + toggle
            if flipped > direct or (flipped == direct and first_source < second_source):
                values[other], sources[other] = flipped, first_source
            else:
                values[other], sources[other] = direct, second_source


@njit(cache=True, inline="always")
def _interval_logs(distance, first_switch, second_switch, stay, toggle):
    q = -.5 * np.expm1(-distance / 25)
    for bit, change in enumerate((q, q, second_switch, first_switch)):
        stay[bit] = np.log1p(-change)
        toggle[bit] = np.log(change)


@njit(cache=True, inline="always")
def _log_dot(left, right, mask):
    value = -np.inf
    for state in range(16):
        value = _log_add(value, left[state] + right[state ^ mask])
    return value


@njit(cache=True, nogil=True)
def _infer_labeled_haplotypes(haplotypes, freqs, genetic_map, switches, errors, weights,
                              state_path, posterior, switch_posterior, origins, copies, origin_errors, copy_errors):
    markers = len(freqs)
    alpha = np.empty((markers if posterior is not None else 0, 16))
    pointer = np.empty((markers, 16), dtype=np.int8)
    log_path = np.zeros(16)
    log_emission = np.empty(16)
    pair_likelihood, pair_log_likelihood = np.empty(8), np.empty(8)
    stay, toggle = np.empty(4), np.empty(4)
    for t in range(markers):
        if not _labeled_emissions(haplotypes[:, t, :], freqs[t], errors[:, t], weights[t], origins[:, t], copies[:, t],
                                  origin_errors[:, t], copy_errors[:, t], log_emission, pair_likelihood,
                                  pair_log_likelihood):
            return t
        if t:
            _interval_logs(genetic_map[t] - genetic_map[t - 1], switches[0, t - 1], switches[1, t - 1], stay, toggle)
            for state in range(16):
                pointer[t, state] = state
            for bit in range(4):
                _mix_max(log_path, pointer[t], 1 << bit, stay[bit], toggle[bit])
        log_path += log_emission
        maximum = np.max(log_path)
        if not np.isfinite(maximum):
            return t
        log_path -= maximum
        if posterior is not None:
            if t:
                alpha[t] = alpha[t - 1]
                for bit in range(4):
                    _mix_sum(alpha[t], 1 << bit, stay[bit], toggle[bit])
                alpha[t] += log_emission
            else:
                alpha[t] = log_emission
            alpha[t] -= np.max(alpha[t])

    state = np.argmax(log_path)
    for t in range(markers - 1, 0, -1):
        state_path[t] = _PAIR_STATES[state]
        state = pointer[t, state]
    state_path[0] = _PAIR_STATES[state]
    if posterior is None:
        return -1

    beta = np.zeros(16)
    work, without_second = np.empty(16), np.empty(16)
    for t in range(markers - 1, -1, -1):
        work[:] = alpha[t] + beta
        maximum = np.max(work)
        pair_likelihood[:] = 0
        for state, projected in enumerate(_PAIR_STATES):
            pair_likelihood[projected] += np.exp(work[state] - maximum)
        _store_posterior(pair_likelihood, pair_likelihood.sum(), posterior[t])
        if not t:
            break

        _labeled_emissions(haplotypes[:, t, :], freqs[t], errors[:, t], weights[t], origins[:, t], copies[:, t],
                           origin_errors[:, t], copy_errors[:, t], log_emission, pair_likelihood, pair_log_likelihood)
        _interval_logs(genetic_map[t] - genetic_map[t - 1], switches[0, t - 1], switches[1, t - 1], stay, toggle)
        beta += log_emission
        _mix_sum(beta, 1, stay[0], toggle[0])
        _mix_sum(beta, 2, stay[1], toggle[1])
        if switch_posterior is not None:
            without_second[:] = beta
            _mix_sum(without_second, 8, stay[3], toggle[3])
            second = _log_dot(alpha[t - 1], without_second, 4) + toggle[2]
        _mix_sum(beta, 4, stay[2], toggle[2])
        if switch_posterior is not None:
            first = _log_dot(alpha[t - 1], beta, 8) + toggle[3]
        _mix_sum(beta, 8, stay[3], toggle[3])
        if switch_posterior is not None:
            total = _log_dot(alpha[t - 1], beta, 0)
            switch_posterior[0, t - 1] = np.exp(first - total)
            switch_posterior[1, t - 1] = np.exp(second - total)
        beta -= np.max(beta)
    return -1


def _parent_labels(labels, markers, name):
    if labels is None:
        return None
    labels = np.asarray(labels)
    if labels.shape != (2, markers) or labels.dtype.kind not in 'biu' or not np.isin(labels, (-1, 0, 1)).all():
        raise ValueError(f"{name} must contain -1, 0 or 1 as integers of shape (2, markers)")
    return np.ascontiguousarray(labels, dtype=np.int8)


def _parent_label_errors(errors, markers, name):
    try:
        errors = np.broadcast_to(np.asarray(errors, dtype=float), (2, markers))
    except ValueError as exc:
        raise ValueError(f"{name} must broadcast to (2, markers)") from exc
    if not np.all((errors > 0) & (errors <= .5)):
        raise ValueError(f"{name} must contain probabilities in (0, 0.5]")
    return errors


def infer_phased_pair(haplotypes, freqs, genetic_map, phase_switch, allele_error, *, weights=None,
                      parent_origin=None, parent_haplotype=None, parent_origin_error=.01,
                      parent_haplotype_error=.01):
    """Decode a pair from biallelic phased hard calls of shape (2, markers, 2).

    Calls are 0 or 1, with -1 and NaN representing missing alleles. Genetic
    coordinates are nondecreasing cM positions. ``phase_switch`` broadcasts
    to (2, markers - 1) and ``allele_error`` broadcasts to (2, markers); the
    first axis selects the sibling. Supply 0.5 for a switch rate at a boundary
    where the next input phase block has arbitrary orientation.

    Optional ``parent_origin`` and ``parent_haplotype`` arrays contain the
    integer labels 0 or 1 or -1 for unresolved, with shape (2, markers).
    Both rows describe transmission from the same observed parent.
    A parent origin is the original input haplotype index of the transmitted
    copy in that child. A parental haplotype is the index in the observed
    parent's input phase. Where both parental copies are called, their
    equality informs IBD sharing from that parent. Origin information from
    either child is propagated through the chromosome with the original
    phase switch rates. Label error inputs broadcast to (2, markers) and
    apply per call as soft observation probabilities in (0, 0.5]. These
    rates represent assumed reliability of the duo labels.

    The joint state path uses code 4*C + 2*U0 + U1. Posterior array shapes are
    (markers, 8), (markers, 3) and (2, markers - 1), respectively. Switch
    posteriors describe toggles of input phase labels, including uncertain
    phase-block boundaries. When family evidence is informative, the state
    path projects one coherent path of both children's parental orientations.
    """
    path, posterior, switches = _decode_phased_pair(
        haplotypes, freqs, genetic_map, phase_switch, allele_error, weights, parent_origin, parent_haplotype,
        parent_origin_error, parent_haplotype_error, posterior=True, full=True,
    )
    ibd_posterior = np.column_stack([posterior[:, _IBD_COUNT == k].sum(axis=1) for k in range(3)])
    return PhasedPairResult(path, _IBD_COUNT[path], posterior, ibd_posterior, switches)


def infer_phased_pair_ibd(haplotypes, freqs, genetic_map, phase_switch, allele_error, *, posterior=True,
                          weights=None, parent_origin=None, parent_haplotype=None, parent_origin_error=.01,
                          parent_haplotype_error=.01):
    """Return the IBD count path and, optionally, P(IBD0, IBD1, IBD2).

    Inputs follow :func:`infer_phased_pair`. Set ``posterior=False`` to decode
    only the hard path and obtain ``None`` for the second result.
    """
    path, probabilities, _ = _decode_phased_pair(
        haplotypes, freqs, genetic_map, phase_switch, allele_error, weights, parent_origin, parent_haplotype,
        parent_origin_error, parent_haplotype_error, posterior=posterior, full=False,
    )
    return _IBD_COUNT[path], probabilities


def _decode_phased_pair(haplotypes, freqs, genetic_map, phase_switch, allele_error, weights, parent_origin,
                        parent_haplotype, parent_origin_error, parent_haplotype_error, *, posterior, full):
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

    origins = _parent_labels(parent_origin, markers, 'parent_origin')
    copies = _parent_labels(parent_haplotype, markers, 'parent_haplotype')
    origin_errors = copies_errors = None
    if origins is not None or copies is not None:
        origin_errors = _parent_label_errors(parent_origin_error, markers, 'parent_origin_error')
        copies_errors = _parent_label_errors(parent_haplotype_error, markers, 'parent_haplotype_error')
        has_origins = origins is not None and np.any((origins >= 0) & (origin_errors < .5))
        has_copies = copies is not None and np.any(np.all((copies >= 0) & (copies_errors < .5), axis=0))
        if has_origins or has_copies:
            if origins is None:
                origins = np.full((2, markers), -1, dtype=np.int8)
            if copies is None:
                copies = np.full((2, markers), -1, dtype=np.int8)
        else:
            origins = copies = origin_errors = copies_errors = None
    path = np.empty(markers, np.int8)
    probabilities = np.empty((markers, 8 if full else 3), np.float32) if posterior else None
    switches = np.empty((2, markers - 1), np.float32) if full else None
    if origins is None:
        failed = _infer_haplotypes(haps, freqs, genetic_map, phase_switch, allele_error, weights, path,
                                  probabilities, switches)
    else:
        failed = _infer_labeled_haplotypes(haps, freqs, genetic_map, phase_switch, allele_error, weights, path,
                                          probabilities, switches, origins, copies, origin_errors, copies_errors)
    if failed >= 0:
        raise ValueError(f"No phased IBD path through marker {failed + 1}; check error rates and allele frequencies")
    return path, probabilities, switches
