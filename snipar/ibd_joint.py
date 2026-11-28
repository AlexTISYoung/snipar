"""Exact inheritance-vector inference for three- and four-sibling families.

The two founder homologs of each parent have arbitrary labels. By defining
that parent's first homolog as the one transmitted to sibling zero at each
marker, the state records two partitions rather than the full transmission
vector. The transition integrates over the anchor sibling's recombination.
This is the Lander–Green inheritance-vector model for unrelated founders,
with independent loci conditional on transmissions (or LD-weighted
likelihoods when weights differ from one).

Only unphased diploid dosage calls enter the emissions. The ordering of
statistically phased alleles does not establish parental origin. Projected
pair marginals do not preserve the joint distribution needed for posterior
imputation that combines inheritance information from more than two siblings.
"""

import math

from numba import njit
import numpy as np


@njit(cache=True, inline="always")
def _homolog(partition, child):
    if child == 0:
        return 0
    return (partition >> (child - 1)) & 1


@njit(cache=True)
def _design(n_siblings, pair_rows):
    groups = 1 << (n_siblings - 1)
    states = groups * groups
    distances = np.zeros(groups, np.int8)
    allele_counts = np.zeros(16, np.int8)
    for i in range(1, groups):
        distances[i] = distances[i >> 1] + (i & 1)
    for i in range(1, 16):
        allele_counts[i] = allele_counts[i >> 1] + (i & 1)
    true_dosage = np.empty((states, 16, n_siblings), np.int8)
    sharing = np.empty((states, len(pair_rows)), np.int8)
    for state in range(states):
        paternal, maternal = state // groups, state % groups
        for child in range(n_siblings):
            f = _homolog(paternal, child)
            m = _homolog(maternal, child) + 2
            for pattern in range(16):
                true_dosage[state, pattern, child] = ((pattern >> f) & 1) + ((pattern >> m) & 1)
        for i in range(len(pair_rows)):
            a, b = pair_rows[i]
            sharing[state, i] = (int(_homolog(paternal, a) == _homolog(paternal, b)) +
                                 int(_homolog(maternal, a) == _homolog(maternal, b)))
    return distances, allele_counts, true_dosage, sharing


@njit(cache=True)
def _emission(calls, frequency, weight, error, allele_counts, true_dosage, evidence, log_likelihood):
    """Integrate true Mendelian genotypes over the four founder alleles."""
    if weight == 0:
        log_likelihood[:] = 0
        return True
    n_sibs = len(calls)
    for i in range(n_sibs):
        call = calls[i]
        for true in range(3):
            if call < 0:
                evidence[i, true] = 1.
            elif call == true:
                evidence[i, true] = 1. - error
            elif abs(call - true) == 1:
                evidence[i, true] = error / 2 if true == 1 else error
            else:
                evidence[i, true] = 0.
    prior = np.empty(5)
    for count in range(5):
        prior[count] = frequency ** count * (1. - frequency) ** (4 - count)
    # Interchanging parents gives the same unphased emission likelihood.
    groups = 1 << (n_sibs - 1)
    maximum = 0.
    for f in range(groups):
        for m in range(f, groups):
            state = f * groups + m
            probability = 0.
            for pattern in range(16):
                term = prior[allele_counts[pattern]]
                for i in range(n_sibs):
                    term *= evidence[i, true_dosage[state, pattern, i]]
                probability += term
            log_likelihood[state] = probability
            log_likelihood[m * groups + f] = probability
            maximum = max(maximum, probability)
    if maximum == 0:
        return False
    normalizer = np.log(maximum)
    for state in range(len(log_likelihood)):
        probability = log_likelihood[state]
        log_likelihood[state] = weight * (np.log(probability) - normalizer) if probability > 0 else -np.inf
    return True


@njit(cache=True, inline="always")
def _transition(distance_cm, n_sibs, distances, kernel, log_kernel):
    # Haldane recombination in each meiosis, on the sex-averaged map.
    r = -.5 * np.expm1(-distance_cm / 50.)
    q = 1. - r
    probabilities = np.empty(n_sibs)
    for changed in range(n_sibs):
        probabilities[changed] = r ** changed * q ** (n_sibs - changed) + r ** (n_sibs - changed) * q ** changed
    for flip_mask in range(len(distances)):
        probability = probabilities[distances[flip_mask]]
        kernel[flip_mask] = probability
        log_kernel[flip_mask] = np.log(probability) if probability > 0 else -np.inf


@njit(cache=True)
def _max_transition(previous, log_kernel, tmp, first_choice, dest, pointer):
    """Apply the paternal and maternal max-product transfers successively."""
    groups = len(log_kernel)
    for f in range(groups):
        for prior_m in range(groups):
            best = 0
            maximum = previous[prior_m] + log_kernel[f]
            for prior_f in range(1, groups):
                candidate = previous[prior_f * groups + prior_m] + log_kernel[f ^ prior_f]
                if candidate > maximum:
                    best = prior_f
                    maximum = candidate
            temp_index = f * groups + prior_m
            tmp[temp_index] = maximum
            first_choice[temp_index] = best * groups + prior_m
    for f in range(groups):
        for m in range(groups):
            first = f * groups
            best = first_choice[first]
            maximum = tmp[first] + log_kernel[m]
            for prior_m in range(1, groups):
                index = first + prior_m
                predecessor = first_choice[index]
                candidate = tmp[index] + log_kernel[m ^ prior_m]
                if candidate > maximum or (candidate == maximum and predecessor < best):
                    best = predecessor
                    maximum = candidate
            dest[first + m] = maximum
            pointer[first + m] = best


@njit(cache=True)
def _sum_transition(previous, kernel, tmp, dest):
    """Multiply by the symmetric transition without a dense joint-state matrix."""
    groups = len(kernel)
    for f in range(groups):
        for prior_m in range(groups):
            total = 0.
            for prior_f in range(groups):
                total += previous[prior_f * groups + prior_m] * kernel[f ^ prior_f]
            tmp[f * groups + prior_m] = total
    for f in range(groups):
        for m in range(groups):
            total = 0.
            first = f * groups
            for prior_m in range(groups):
                total += tmp[first + prior_m] * kernel[m ^ prior_m]
            dest[first + m] = total


@njit(cache=True, inline="always")
def _advance(marker, log_previous, sum_previous, calls, freqs, genetic_map, weights, errors,
             distances, counts, dosage, evidence, log_emissions, likelihood, kernel, log_kernel,
             log_work, sum_work, temp, first_choice, pointer, posterior):
    if not _emission(calls[:, marker], freqs[marker], weights[marker], errors[marker],
                     counts, dosage, evidence, log_emissions):
        return False
    if marker:
        _transition(genetic_map[marker] - genetic_map[marker - 1], calls.shape[0], distances, kernel, log_kernel)
        _max_transition(log_previous, log_kernel, temp, first_choice, log_work, pointer)
        if posterior:
            _sum_transition(sum_previous, kernel, temp, sum_work)
    else:
        log_work[:] = log_previous
        if posterior:
            sum_work[:] = sum_previous
    log_work += log_emissions
    maximum = np.max(log_work)
    if not np.isfinite(maximum):
        return False
    log_previous[:] = log_work - maximum
    if posterior:
        likelihood[:] = np.exp(log_emissions)
        sum_work *= likelihood
        total = sum_work.sum()
        if total <= 0:
            return False
        sum_previous[:] = sum_work / total
    return True


@njit(cache=True)
def _infer(calls, pair_rows, freqs, genetic_map, weights, errors, hard, marginal, block_size):
    n_sibs, n_markers = calls.shape
    n_blocks = (n_markers + block_size - 1) // block_size
    distances, counts, dosage, sharing = _design(n_sibs, pair_rows)
    states = len(dosage)
    posterior = marginal.size > 0
    log_checkpoints = np.empty((n_blocks, states))
    sum_checkpoints = np.empty((n_blocks if posterior else 0, states))
    alpha_block = np.empty((block_size if posterior else 0, states))
    emission_block = np.empty_like(alpha_block)
    pointers = np.empty((block_size, states), np.int8)
    log_previous = np.zeros(states)
    sum_previous = np.full(states, 1. / states)
    work = np.empty((5, states))
    log_work, sum_work, temp, likelihood, log_emissions = work
    evidence = np.empty((n_sibs, 3))
    first_choice = np.empty(states, np.int8)
    kernel = np.empty(len(distances))
    log_kernel = np.empty_like(kernel)

    for block in range(n_blocks):
        start = block * block_size
        end = min(start + block_size, n_markers)
        log_checkpoints[block] = log_previous
        if posterior:
            sum_checkpoints[block] = sum_previous
        for marker in range(start, end):
            if not _advance(marker, log_previous, sum_previous, calls, freqs, genetic_map, weights, errors,
                            distances, counts, dosage, evidence, log_emissions, likelihood, kernel, log_kernel,
                            log_work, sum_work, temp, first_choice, pointers[0], posterior):
                return marker

    state = np.argmax(log_previous)
    beta = np.ones(states)
    joint = np.empty(states)
    for block in range(n_blocks - 1, -1, -1):
        start = block * block_size
        end = min(start + block_size, n_markers)
        log_previous[:] = log_checkpoints[block]
        if posterior:
            sum_previous[:] = sum_checkpoints[block]
        # Recompute only this block to recover its local traceback and forward
        # probabilities. Loci outside the block use boundary checkpoints.
        for marker in range(start, end):
            local = marker - start
            _advance(marker, log_previous, sum_previous, calls, freqs, genetic_map, weights, errors,
                     distances, counts, dosage, evidence, log_emissions, likelihood, kernel, log_kernel,
                     log_work, sum_work, temp, first_choice, pointers[local], posterior)
            if posterior:
                alpha_block[local] = sum_previous
                emission_block[local] = likelihood
        for marker in range(end - 1, start - 1, -1):
            local = marker - start
            for i in range(len(pair_rows)):
                hard[i, marker] = sharing[state, i]
            if posterior:
                joint[:] = alpha_block[local] * beta
                total = joint.sum()
                if total <= 0:
                    return marker
                joint /= total
                for i in range(len(pair_rows)):
                    p0 = p1 = p2 = 0.
                    for j in range(states):
                        share = sharing[j, i]
                        if share == 0:
                            p0 += joint[j]
                        elif share == 1:
                            p1 += joint[j]
                        else:
                            p2 += joint[j]
                    marginal[i, marker, 0] = p0
                    marginal[i, marker, 1] = p1
                    marginal[i, marker, 2] = p2
                if marker:
                    likelihood[:] = emission_block[local] * beta
                    _transition(genetic_map[marker] - genetic_map[marker - 1], n_sibs,
                                distances, kernel, log_kernel)
                    _sum_transition(likelihood, kernel, temp, beta)
                    scale = beta.sum()
                    if scale <= 0:
                        return marker
                    beta /= scale
            if marker:
                state = pointers[local, state]
    return -1


def infer_joint_family(calls, pair_rows, freqs, genetic_map, weights, error_probs, *, posterior=False):
    """Return a joint Viterbi path projected onto the requested sibling pairs.

    ``calls`` has shape (siblings, markers) with unphased 0/1/2 dosages and
    missing calls encoded as -1 or NaN. ``pair_rows`` indexes rows of ``calls``.
    All marker vectors share that column order. When ``posterior`` is true,
    also return exact pair marginals of shape (pairs, markers, 3).
    """
    calls = np.asarray(calls)
    if calls.ndim != 2 or calls.shape[0] not in (3, 4):
        raise ValueError("Joint IBD requires three or four full siblings")
    n_sibs, n_markers = calls.shape
    if calls.dtype.kind not in 'biuf' or not np.all(np.isnan(calls) | np.isin(calls, (-1, 0, 1, 2))):
        raise ValueError("Joint IBD calls must be 0, 1, 2, -1 or NaN")
    if calls.dtype.kind == 'f':
        calls = np.where(np.isnan(calls), -1, calls)
    calls = np.asfortranarray(calls, dtype=np.int8)
    pair_rows = np.asarray(pair_rows)
    if (pair_rows.ndim != 2 or pair_rows.shape[1] != 2 or not len(pair_rows)
            or pair_rows.dtype.kind not in 'iu' or np.any(pair_rows < 0)
            or np.any(pair_rows >= n_sibs) or np.any(pair_rows[:, 0] == pair_rows[:, 1])):
        raise ValueError("Each IBD pair must contain two distinct sibling row indices")
    pair_rows = np.ascontiguousarray(pair_rows, dtype=np.intp)
    freqs, genetic_map, weights, error_probs = [np.asarray(x, dtype=float)
                                                for x in (freqs, genetic_map, weights, error_probs)]
    if any(x.shape != (n_markers,) or not np.all(np.isfinite(x))
           for x in (freqs, genetic_map, weights, error_probs)):
        raise ValueError("Joint IBD marker vectors must be finite and match the calls")
    if np.any(np.diff(genetic_map) < 0):
        raise ValueError("Joint IBD genetic positions must be nondecreasing")
    if np.any((freqs < 0) | (freqs > 1)) or np.any((error_probs < 0) | (error_probs > 1)):
        raise ValueError("Joint IBD allele frequencies and error probabilities must be in [0, 1]")
    if np.any(weights < 0):
        raise ValueError("Joint IBD marker weights must be nonnegative")
    hard = np.empty((len(pair_rows), n_markers), np.int8)
    marginal = np.empty((len(pair_rows), n_markers, 3) if posterior else (0, 0, 3), np.float32)
    if n_markers:
        block_size = max(1, math.isqrt(n_markers))
        failed = _infer(calls, pair_rows, freqs, genetic_map, weights, error_probs, hard, marginal, block_size)
        if failed >= 0:
            raise ValueError(f"No joint Mendelian path through marker {failed + 1}; check genotype error rates and allele frequencies")
    return hard, marginal if posterior else None


@njit(cache=True)
def _pair_components(pairs, n_members):
    parents = np.arange(n_members)
    sizes = np.ones(n_members, dtype=np.intp)
    for left, right in pairs:
        while parents[left] != left:
            parents[left] = parents[parents[left]]
            left = parents[left]
        while parents[right] != right:
            parents[right] = parents[parents[right]]
            right = parents[right]
        if left != right:
            if sizes[left] < sizes[right]:
                left, right = right, left
            parents[right] = left
            sizes[left] += sizes[right]
    labels = np.empty(len(pairs), dtype=np.intp)
    for i in range(len(pairs)):
        first = pairs[i, 0]
        while parents[first] != first:
            first = parents[first]
        labels[i] = first
    return labels


def joint_pair_groups(sibpairs):
    """Return original pair-row indices of complete three- or four-sibling cliques.

    Incomplete, repeated, overlapping or larger pair components are left to
    pairwise inference. The input may use IIDs or row numbers.
    """
    pairs = np.asarray(sibpairs)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError("Sibling pairs must have shape (pairs, 2)")
    if not len(pairs):
        return []
    # Chromosome readers often already use compact, nonnegative row numbers.
    if pairs.dtype.kind in 'iu' and pairs.min() >= 0 and pairs.max() < pairs.size:
        n_members = int(pairs.max()) + 1
        pairs = pairs.astype(np.intp, copy=False)
    else:
        members, indices = np.unique(pairs.ravel(), return_inverse=True)
        n_members = len(members)
        pairs = indices.reshape(-1, 2)
    labels = _pair_components(pairs, n_members)
    counts = np.bincount(labels)[labels]
    eligible = np.flatnonzero((counts == 3) | (counts == 6))
    if not eligible.size:
        return []
    eligible = eligible[np.argsort(labels[eligible], kind='stable')]
    split = np.flatnonzero(labels[eligible[1:]] != labels[eligible[:-1]]) + 1
    groups = []
    for edges in np.split(eligible, split):
        group_pairs = np.sort(pairs[edges], axis=1)
        n_sibs = len(np.unique(group_pairs))
        if (n_sibs in (3, 4) and len(edges) == math.comb(n_sibs, 2)
                and np.all(group_pairs[:, 0] != group_pairs[:, 1])
                and len(np.unique(group_pairs, axis=0)) == len(edges)):
            groups.append(edges)
    return sorted(groups, key=lambda edges: edges[0])
