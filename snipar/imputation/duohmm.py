"""Infer chromosome-local parent of origin from statistically phased duos.

The four states identify the IBD haplotype in the observed parent and in the
child, as in O'Connell et al. (2014), doi:10.1371/journal.pgen.1004234. A change
in parental state can be due to a crossover or a phase switch. A change in
child state reflects a phase switch. This implementation estimates the
transmitted haplotype labels without changing the input genotypes.
"""

import numpy as np
from numba import config, get_num_threads, njit, prange, set_num_threads


# Phase switches are modeled between successive heterozygous markers, where
# a change of phase is observable. The probabilities are fixed per duo.
_SWITCH_PROBABILITY = 0.001
_MISMATCH_PROBABILITY = 0.001


def recombination_probabilities(positions, genetic_positions=None):
    """Return per-marker crossover probabilities for an ordered chromosome.

    The entry at each marker corresponds to the interval ending there. Use
    genetic positions in cM when available, or approximate at 1 cM per Mb.
    """
    positions = np.asarray(positions, dtype=float)
    if (positions.ndim != 1 or not np.isfinite(positions).all()
            or (positions < 0).any() or (np.diff(positions) < 0).any()):
        raise ValueError("Physical positions must be finite and nondecreasing")
    if genetic_positions is None:
        distances = np.diff(positions) * 1e-8
    else:
        genetic_positions = np.asarray(genetic_positions, dtype=float)
        if (genetic_positions.shape != positions.shape or not np.isfinite(genetic_positions).all()
                or (genetic_positions < 0).any() or (np.diff(genetic_positions) < 0).any()):
            raise ValueError("Genetic positions must be finite, nondecreasing and aligned with markers")
        distances = np.diff(genetic_positions) / 100
    probabilities = np.zeros(positions.size)
    probabilities[1:] = -0.5 * np.expm1(-2 * distances)
    return probabilities


@njit(cache=True, inline="always")
def _called(allele):
    return allele == 0 or allele == 1


@njit(cache=True, inline="always")
def _marker_weights(parent, child, marker, crossover, transitions, emissions):
    parent_switch = 0.0
    child_switch = 0.0
    if _called(parent[marker, 0]) and _called(parent[marker, 1]) and parent[marker, 0] != parent[marker, 1]:
        parent_switch = _SWITCH_PROBABILITY
    if _called(child[marker, 0]) and _called(child[marker, 1]) and child[marker, 0] != child[marker, 1]:
        child_switch = _SWITCH_PROBABILITY
    parent_change = crossover + parent_switch - 2 * crossover * parent_switch
    transitions[0] = (1 - parent_change) * (1 - child_switch)
    transitions[1] = (1 - parent_change) * child_switch
    transitions[2] = parent_change * (1 - child_switch)
    transitions[3] = parent_change * child_switch

    # With an unobserved allele, the likelihood for the inherited copy
    # averages over a uniform background. Missing calls are never treated
    # as matches or used to assign transmitted alleles.
    for state in range(4):
        p = parent[marker, state // 2]
        c = child[marker, state % 2]
        if not _called(p) or not _called(c):
            emissions[state] = 0.5
        elif p == c:
            emissions[state] = 1 - _MISMATCH_PROBABILITY
        else:
            emissions[state] = _MISMATCH_PROBABILITY


@njit(cache=True, inline="always")
def _mix_states(state, transitions, i):
    # A bit in (previous_state ^ current_state) identifies each haplotype
    # label that changed. Transition probabilities are symmetric.
    return (state[0] * transitions[i] + state[1] * transitions[i ^ 1]
            + state[2] * transitions[i ^ 2] + state[3] * transitions[i ^ 3])


@njit(cache=True, inline="always")
def _normalize(source, target):
    total = source[0] + source[1] + source[2] + source[3]
    for i in range(4):
        target[i] = source[i] / total


@njit(cache=True)
def _forward(parent, child, recombination, initial, final, history):
    final[:] = initial
    transitions = np.empty(4)
    emissions = np.empty(4)
    updated = np.empty(4)
    for marker in range(recombination.size):
        _marker_weights(parent, child, marker, recombination[marker], transitions, emissions)
        for state in range(4):
            updated[state] = _mix_states(final, transitions, state) * emissions[state]
        _normalize(updated, final)
        if history is not None:
            history[marker, :] = final


@njit(cache=True, inline="always")
def _call_labels(parent, child, marker, minimum, alpha, beta):
    w0 = alpha[0] * beta[0]
    w1 = alpha[1] * beta[1]
    w2 = alpha[2] * beta[2]
    w3 = alpha[3] * beta[3]
    total = w0 + w1 + w2 + w3
    if not (_called(child[marker, 0]) and _called(child[marker, 1])):
        return -1, -1
    if w0 + w2 >= minimum * total:
        child_haplotype, p0, p1 = 0, w0, w2
    elif w1 + w3 >= minimum * total:
        child_haplotype, p0, p1 = 1, w1, w3
    else:
        return -1, -1

    inherited = child[marker, child_haplotype]
    if parent[marker, 0] != inherited and parent[marker, 1] != inherited:
        return -1, -1
    # Call the parental copy only when its joint posterior with the child
    # copy is high. A child origin can be clear with an ambiguous parent copy.
    if p0 >= minimum * total and parent[marker, 0] == inherited:
        return child_haplotype, 0
    if p1 >= minimum * total and parent[marker, 1] == inherited:
        return child_haplotype, 1
    return child_haplotype, -1


@njit(cache=True)
def _backward(parent, child, recombination, history, initial, final, child_origin, parent_copy, minimum):
    final[:] = initial
    transitions = np.empty(4)
    emissions = np.empty(4)
    updated = np.empty(4)
    for marker in range(recombination.size - 1, -1, -1):
        child_origin[marker], parent_copy[marker] = _call_labels(parent, child, marker, minimum,
                                                               history[marker], final)
        _marker_weights(parent, child, marker, recombination[marker], transitions, emissions)
        for state in range(4):
            emissions[state] *= final[state]
        for state in range(4):
            updated[state] = _mix_states(emissions, transitions, state)
        _normalize(updated, final)


@njit(cache=True, parallel=True)
def _advance_pairs(phased, pairs, recombination, previous_alpha):
    next_alpha = np.empty_like(previous_alpha)
    for i in prange(len(pairs)):
        parent, child = pairs[i]
        _forward(phased[parent], phased[child], recombination, previous_alpha[i], next_alpha[i], None)
    return next_alpha


@njit(cache=True, parallel=True)
def _call_pairs(phased, pairs, recombination, previous_alpha, next_beta, minimum, workers):
    child_origin = np.empty((len(pairs), recombination.size), dtype=np.int8)
    parent_copy = np.empty_like(child_origin)
    previous_beta = np.empty_like(next_beta)
    for worker in prange(workers):
        history = np.empty((recombination.size, 4))
        final_alpha = np.empty(4)
        for i in range(worker, len(pairs), workers):
            parent, child = pairs[i]
            _forward(phased[parent], phased[child], recombination, previous_alpha[i], final_alpha, history)
            _backward(phased[parent], phased[child], recombination, history, next_beta[i], previous_beta[i],
                      child_origin[i], parent_copy[i], minimum)
    return child_origin, parent_copy, previous_beta


def _block_inputs(phased, pairs, recombination, threads):
    phased = np.asarray(phased)
    pairs = np.asarray(pairs, dtype=np.intp)
    recombination = np.asarray(recombination, dtype=float)
    if phased.ndim != 3 or phased.shape[2] != 2:
        raise ValueError("Phased genotypes must have shape (samples, markers, 2)")
    if (pairs.ndim != 2 or pairs.shape[1] != 2
            or (pairs < 0).any() or (pairs >= phased.shape[0]).any()):
        raise ValueError("Duo rows must contain valid parent and child sample indices")
    if (recombination.shape != (phased.shape[1],) or not np.isfinite(recombination).all()
            or (recombination < 0).any() or (recombination > 0.5).any()):
        raise ValueError("Recombination probabilities must be in [0, 0.5] and aligned with markers")
    if threads < 1:
        raise ValueError("threads must be positive")
    workers = max(1, min(len(pairs), threads, config.NUMBA_NUM_THREADS))
    return phased, pairs, recombination, workers


def _boundary(boundary, n_pairs):
    if boundary is None:
        return np.full((n_pairs, 4), 0.25)
    boundary = np.asarray(boundary, dtype=float)
    if boundary.shape != (n_pairs, 4):
        raise ValueError("HMM boundary state must have four probabilities per duo")
    return boundary


def advance_parental_origin(phased, pairs, recombination, previous_alpha=None, *, threads=1):
    """Return normalized forward probabilities at the last marker of a block.

    ``pairs`` contains (parent, child) sample row indices. ``recombination``
    is a slice from the chromosome-wide array, including the interval from
    the previous block's last marker. ``previous_alpha`` is the checkpoint
    at that preceding marker or ``None`` at the start of the chromosome.
    """
    phased, pairs, recombination, workers = _block_inputs(phased, pairs, recombination, threads)
    previous_alpha = _boundary(previous_alpha, len(pairs))
    old_threads = get_num_threads()
    set_num_threads(workers)
    try:
        return _advance_pairs(phased, pairs, recombination, previous_alpha)
    finally:
        set_num_threads(old_threads)


def call_parental_origin(phased, pairs, recombination, previous_alpha=None, next_beta=None, *,
                         threads=1, minimum_posterior=0.99):
    """Return child origins, parent copies and the left backward checkpoint.

    Both label matrices are int8 of shape (pairs, markers), with 0 or 1 for
    a confident input haplotype index and -1 when unresolved. ``next_beta``
    holds the likelihood of future markers given the state at this block's
    last marker. Use ``None`` at the chromosome end and pass the returned
    checkpoint when decoding the adjacent block to the left.
    """
    phased, pairs, recombination, workers = _block_inputs(phased, pairs, recombination, threads)
    if not 0.5 < minimum_posterior <= 1:
        raise ValueError("minimum_posterior must be greater than 0.5 and at most 1")
    previous_alpha = _boundary(previous_alpha, len(pairs))
    next_beta = _boundary(next_beta, len(pairs))
    old_threads = get_num_threads()
    set_num_threads(workers)
    try:
        return _call_pairs(phased, pairs, recombination, previous_alpha, next_beta, minimum_posterior, workers)
    finally:
        set_num_threads(old_threads)


def infer_parental_origin(parent_haplotypes, child_haplotypes, positions, genetic_positions=None, *,
                          minimum_posterior=0.99):
    """Return one child's observed-parent haplotype indices or -1 if unresolved.

    Haplotypes have shape (markers, 2). A call of 0 or 1 identifies the
    transmitted copy in the child, so the other parent's transmitted copy
    is ``1 - origin``. Alleles other than 0 and 1 are treated as missing.
    """
    phased = np.stack((parent_haplotypes, child_haplotypes))
    recombination = recombination_probabilities(positions, genetic_positions)
    origins, _, _ = call_parental_origin(phased, [[0, 1]], recombination,
                                         minimum_posterior=minimum_posterior)
    return origins[0]
