from itertools import combinations, product

import numpy as np
import pytest

from snipar.ibd_joint import _design, _emission, infer_joint_family, joint_pair_groups


def _full_inheritance_reference(calls, pairs, freqs, genetic_map, weights, errors):
    """Dense reference over all 4**n labeled transmissions."""
    n_sibs, n_markers = calls.shape
    vectors = np.array(list(product((0, 1), repeat=2 * n_sibs))).reshape(-1, 2, n_sibs)
    n_states = len(vectors)
    alleles = np.array(list(product((0, 1), repeat=4))).reshape(-1, 2, 2)
    true = (alleles[:, 0, vectors[:, 0, :]] + alleles[:, 1, vectors[:, 1, :]]).transpose(1, 0, 2)
    sharing = (vectors[:, :, pairs[:, 0]] == vectors[:, :, pairs[:, 1]]).sum(axis=1)
    founder_counts = alleles.sum(axis=(1, 2))
    emissions = []
    for site, (f, w, e) in enumerate(zip(freqs, weights, errors)):
        genotyping_error = np.array([[1 - e, e / 2, 0],
                                     [e, 1 - e, e],
                                     [0, e / 2, 1 - e],
                                     [1, 1, 1]])
        likelihood = genotyping_error[np.where(calls[:, site] < 0, 3, calls[:, site]), true]
        founder_prob = f ** founder_counts * (1 - f) ** (4 - founder_counts)
        emissions.append((likelihood.prod(axis=2) @ founder_prob) ** w)
    emissions = np.array(emissions)
    distances = (vectors[:, None] != vectors[None, :]).sum(axis=(2, 3))
    switches = (1 - np.exp(-np.diff(genetic_map) / 50)) / 2
    transitions = np.array([r ** distances * (1 - r) ** (2 * n_sibs - distances) for r in switches])

    alpha = np.empty_like(emissions)
    alpha[0] = emissions[0] / n_states
    for t in range(1, n_markers):
        alpha[t] = (alpha[t - 1] @ transitions[t - 1]) * emissions[t]
    beta = np.ones_like(emissions)
    for t in range(n_markers - 2, -1, -1):
        beta[t] = transitions[t] @ (emissions[t + 1] * beta[t + 1])
    gamma = alpha * beta
    gamma /= gamma.sum(axis=1, keepdims=True)
    marginals = np.stack([gamma @ (sharing == share) for share in range(3)], axis=2).transpose(1, 0, 2)

    # Collapse the dense reference to the same relative transmission partitions
    # by summing over all labels in each target class.
    relative = vectors[:, :, 1:] ^ vectors[:, :, :1]
    code = relative @ (2 ** np.arange(n_sibs - 1))
    _, representatives, state_groups = np.unique(code, axis=0, return_index=True, return_inverse=True)
    membership = state_groups.ravel()[:, None] == np.arange(len(representatives))
    collapsed_transitions = transitions[:, representatives, :] @ membership
    return marginals, emissions[:, representatives], collapsed_transitions, sharing[representatives]


def _viterbi_score(emissions, transitions, states_allowed):
    with np.errstate(divide='ignore'):
        score = np.log(emissions[0]) - np.log(emissions.shape[1])
        score[~states_allowed[0]] = -np.inf
        for t, transition in enumerate(transitions, 1):
            score = np.max(score[:, None] + np.log(transition), axis=0) + np.log(emissions[t])
            score[~states_allowed[t]] = -np.inf
    return np.max(score)


@pytest.mark.parametrize('n_sibs', [3, 4])
def test_exact_joint_posteriors_and_valid_hard_traceback(n_sibs):
    rng = np.random.default_rng(182 + n_sibs)
    loci = 13  # Spans multiple checkpoint boundaries.
    calls = rng.integers(-1, 3, size=(n_sibs, loci), dtype=np.int8)
    freqs = rng.uniform(.15, .45, loci)
    increments = np.array([0, .5, 0, .1, 5, 20, .2, .02, .5, 2, 0, .3, .5])
    genetic_map = np.cumsum(increments)
    weights = rng.uniform(.25, 1., loci)
    weights[3] = 0
    errors = rng.uniform(.01, .04, loci)
    pairs = np.array(list(combinations(range(n_sibs), 2)))[::-1, ::-1].copy()
    hard, pair_posterior = infer_joint_family(calls, pairs, freqs, genetic_map, weights, errors, posterior=True)
    hard_only, absent_posterior = infer_joint_family(calls, pairs, freqs, genetic_map, weights, errors)
    assert absent_posterior is None
    np.testing.assert_array_equal(hard_only, hard)
    expected, emissions, transitions, state_sharing = _full_inheritance_reference(
        calls, pairs, freqs, genetic_map, weights, errors)
    np.testing.assert_allclose(pair_posterior, expected, rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(pair_posterior.sum(axis=2), 1., atol=1e-7)

    # The projected hard decisions must fit one entire path with maximum joint
    # probability, even when parental role symmetry admits several traces.
    allowed = np.all(state_sharing[None] == hard.T[:, None, :], axis=2)
    unconstrained = np.ones_like(allowed)
    assert np.isfinite(_viterbi_score(emissions, transitions, allowed))
    np.testing.assert_allclose(_viterbi_score(emissions, transitions, allowed),
                               _viterbi_score(emissions, transitions, unconstrained), rtol=0, atol=1e-10)


def test_emission_enforces_higher_order_mendelian_constraint():
    pairs = np.array([[0, 1], [0, 2], [1, 2]])
    _, alt_allele_counts, dosage, sharing = _design(3, pairs)
    # Anchor inherits F0/M0, second sibling F0/M1, third sibling F1/M0.
    # A homozygous-reference call in both other siblings forces both of the
    # anchor's alleles to reference. Each pair in isolation allows (1, 0, 0).
    state = 2 * 4 + 1
    np.testing.assert_array_equal(sharing[state], [1, 1, 0])
    log_p = np.empty(16)
    evidence = np.empty((3, 3))
    calls = np.array([1, 0, 0], np.int8)
    assert _emission(calls, .3, 1., 0., alt_allele_counts, dosage, evidence, log_p)
    assert np.isneginf(log_p[state])
    assert _emission(calls, .3, 1., .002, alt_allele_counts, dosage, evidence, log_p)
    assert np.isfinite(log_p[state])


def test_all_missing_calls_retain_the_mendelian_pair_prior():
    calls = np.full((4, 7), np.nan)
    pairs = np.array([[0, 2], [3, 1]])
    n_loci = calls.shape[1]
    hard, marginals = infer_joint_family(calls, pairs, np.full(n_loci, .3), np.arange(n_loci),
                                         np.ones(n_loci), np.zeros(n_loci), posterior=True)
    np.testing.assert_allclose(marginals, np.broadcast_to([.25, .5, .25], marginals.shape), rtol=0, atol=1e-7)
    np.testing.assert_array_equal(hard, 2)


def test_joint_path_rejects_inconsistent_data_under_zero_error():
    with pytest.raises(ValueError, match='No joint Mendelian path through marker 2'):
        infer_joint_family(np.array([[0, 2], [0, 2], [0, 2]]), [[0, 2]], [0., 0.], [0., 0.], [1., 1.], [0., 0.])


@pytest.mark.parametrize('numeric', [False, True])
def test_only_complete_small_cliques_are_grouped(numeric):
    edges = [
        ('a', 'c'), ('p', 'q'), ('i', 'j'), ('u', 'v'), ('k', 'j'),
        ('q', 'r'), ('a', 'b'), ('u', 'w'), ('i', 'k'), ('b', 'c'),
        ('p', 'q'), ('u', 'x'), ('v', 'w'), ('v', 'x'), ('w', 'x'),
        *combinations('12345', 2),
        ('x0', 'x1'), ('x1', 'x2'), ('x2', 'x0'), ('x2', 'x3'),
        ('self', 'self'),
    ]
    if numeric:
        ids = {key: i * 2000 for i, key in enumerate(sorted({iid for edge in edges for iid in edge}))}
        edges = [[ids[left], ids[right]] for left, right in edges]
    observed = joint_pair_groups(edges)
    expected = [[0, 6, 9], [2, 4, 8], [3, 7, 11, 12, 13, 14]]
    assert [group.tolist() for group in observed] == expected
