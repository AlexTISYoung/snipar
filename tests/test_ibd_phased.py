from itertools import product

import numpy as np
import pytest

from snipar.ibd_phased import _emissions, _transition, infer_phased_pair, infer_phased_pair_ibd


def _labeled_reference(haplotypes, freqs, genetic_map, switches, errors, weights, *,
                       origins=None, copies=None, origin_error=None, copy_error=None):
    states = np.array(list(product((0, 1), repeat=4)))  # O1, O2, Bp, Bq
    markers = len(freqs)
    founder = np.array(list(product((0, 1), repeat=4)))
    emission = np.empty((markers, len(states)))
    state_codes = np.empty(len(states), int)
    for i, (o1, o2, bp, bq) in enumerate(states):
        t = np.empty((len(founder), 2, 2), int)
        t[:, 0, o1] = founder[:, 0]
        t[:, 0, 1 - o1] = founder[:, 2]
        t[:, 1, o2] = founder[:, 0 if bp else 1]
        t[:, 1, 1 - o2] = founder[:, 2 if bq else 3]
        counts = founder.sum(axis=1)
        for marker in range(markers):
            observed = haplotypes[:, marker, :]
            error = errors[:, marker, None]
            factors = np.where(observed < 0, 1., np.where(observed == t, 1 - error, error))
            f = freqs[marker]
            emission[marker, i] = (factors.prod(axis=(1, 2)) @ (f ** counts * (1 - f) ** (4 - counts))) ** weights[marker]
            if origins is not None:
                for sibling, copy in enumerate((o1, o2)):
                    call = origins[sibling, marker]
                    if call >= 0:
                        e = origin_error[sibling, marker]
                        emission[marker, i] *= (1 - e) if call == copy else e
            if copies is not None:
                probability = 0.
                for first in range(2):
                    term = .5
                    for sibling, copy in enumerate((first, first if bp else 1 - first)):
                        call = copies[sibling, marker]
                        if call >= 0:
                            e = copy_error[sibling, marker]
                            term *= (1 - e) if call == copy else e
                    probability += term
                emission[marker, i] *= probability
        state_codes[i] = 4 * (o1 ^ o2) + 2 * (bp if o1 == 0 else bq) + (bq if o1 == 0 else bp)

    changes = states[:, None] != states
    rate = (1 - np.exp(-np.diff(genetic_map) / 25)) / 2
    all_rates = np.column_stack((switches.T, rate, rate))
    transition = np.prod(np.where(changes, all_rates[:, None, None, :], 1 - all_rates[:, None, None, :]), axis=3)
    alpha = np.empty_like(emission)
    alpha[0] = emission[0] / 16
    trace = np.empty_like(emission, dtype=np.int8)
    best = emission[0] / 16
    for marker in range(1, markers):
        alpha[marker] = (alpha[marker - 1] @ transition[marker - 1]) * emission[marker]
        scores = best[:, None] * transition[marker - 1]
        trace[marker] = scores.argmax(axis=0)
        best = scores.max(axis=0) * emission[marker]
    path = np.empty(markers, dtype=np.int8)
    state = best.argmax()
    for marker in range(markers - 1, 0, -1):
        path[marker] = state_codes[state]
        state = trace[marker, state]
    path[0] = state_codes[state]
    beta = np.ones_like(emission)
    for marker in range(markers - 2, -1, -1):
        beta[marker] = transition[marker] @ (emission[marker + 1] * beta[marker + 1])
    gamma = alpha * beta
    gamma /= gamma.sum(axis=1, keepdims=True)
    orientation_posteriors = gamma @ (state_codes[:, None] == np.arange(8))
    switch_posteriors = np.empty((2, markers - 1))
    for marker in range(markers - 1):
        joint = alpha[marker, :, None] * transition[marker] * (emission[marker + 1] * beta[marker + 1])
        for sibling in range(2):
            switch_posteriors[sibling, marker] = joint[changes[:, :, sibling]].sum() / joint.sum()
    return orientation_posteriors, switch_posteriors, path


def test_eight_state_model_matches_labeled_reference_and_switch_posteriors():
    rng = np.random.default_rng(92317)
    markers = 11
    haplotypes = rng.integers(-1, 2, (2, markers, 2))
    freqs = rng.uniform(.15, .49, markers)
    genetic_map = np.cumsum(rng.uniform(0, 2, markers))
    genetic_map[3] = genetic_map[2]  # No meiotic switch at this interval.
    switches = rng.uniform(.001, .3, (2, markers - 1))
    switches[0, 6] = .5  # The next input phase block has no known orientation.
    errors = rng.uniform(.001, .02, (2, markers))
    weights = rng.uniform(.25, 1, markers)
    weights[7] = 0
    expected, expected_switches, _ = _labeled_reference(haplotypes, freqs, genetic_map, switches, errors, weights)
    result = infer_phased_pair(haplotypes, freqs, genetic_map, switches, errors, weights=weights)
    np.testing.assert_allclose(result.state_posterior, expected, atol=1e-7)
    np.testing.assert_allclose(result.phase_switch_posterior, expected_switches, atol=1e-7)
    np.testing.assert_allclose(result.ibd_posterior.sum(axis=1), 1, atol=1e-7)
    np.testing.assert_array_equal(result.hard_ibd, (result.state_path & 1) + ((result.state_path >> 1) & 1))


def test_ibd_zero_alignments_have_same_emissions_and_different_transitions():
    likelihood, logs = np.empty(8), np.empty(8)
    assert _emissions(np.array([[0, 1], [1, 0]]), .21, np.array([.01, .02]), 1., likelihood, logs)
    np.testing.assert_allclose(likelihood[0], likelihood[4], rtol=0, atol=1e-15)
    transitions = [np.empty((8, 8)) for _ in range(3)]
    _transition(.1, 0., 0., *transitions)
    matrix = transitions[0]
    assert matrix[0, 1] > 0 and matrix[4, 1] == 0
    assert matrix[4, 5] > 0 and matrix[0, 5] == 0
    np.testing.assert_allclose(matrix.sum(axis=1), 1.)


def test_ibd_two_relative_flip_does_not_identify_sibling_with_equal_switch_rates():
    markers, cut = 26, 13
    haplotypes = np.repeat(np.array([[[0, 1]], [[0, 1]]]), markers, axis=1)
    haplotypes[0, cut:] = haplotypes[0, cut:, ::-1]
    freqs = np.full(markers, .5)
    genetic_map = np.zeros(markers)  # The relative switch must be in input phase.
    symmetric = infer_phased_pair(haplotypes, freqs, genetic_map, .001, 0.)
    np.testing.assert_array_equal(symmetric.hard_ibd, np.full(markers, 2))
    np.testing.assert_array_equal(symmetric.state_path >> 2, np.arange(markers) >= cut)
    np.testing.assert_allclose(symmetric.phase_switch_posterior[:, cut - 1], [.5, .5], atol=1e-6)
    asymmetric = infer_phased_pair(haplotypes, freqs, genetic_map, [[.03], [.0001]], 0.)
    assert asymmetric.phase_switch_posterior[0, cut - 1] > .99
    assert asymmetric.phase_switch_posterior[1, cut - 1] < .01


def test_ibd_one_can_identify_a_flip_of_the_sharing_haplotype():
    markers, cut = 34, 17
    haplotypes = np.repeat(np.array([[[0, 1]], [[0, 0]]]), markers, axis=1)
    haplotypes[0, cut:] = haplotypes[0, cut:, ::-1]
    result = infer_phased_pair(haplotypes, np.full(markers, .5), np.zeros(markers), .001, 0.)
    np.testing.assert_array_equal(result.hard_ibd, np.ones(markers))
    assert result.phase_switch_posterior[0, cut - 1] > .99
    np.testing.assert_allclose(result.phase_switch_posterior[1, cut - 1], .001, atol=1e-7)


def test_missing_haplotypes_retain_initial_ibd_prior():
    calls = np.full((2, 1, 2), np.nan)
    result = infer_phased_pair(calls, [.31], [0.], .01, .02)
    np.testing.assert_allclose(result.state_posterior[0], np.full(8, 1/8))
    np.testing.assert_allclose(result.ibd_posterior[0], [.25, .5, .25])
    assert result.phase_switch_posterior.shape == (2, 0)


def test_duo_evidence_matches_labeled_reference_with_uneven_calls():
    rng = np.random.default_rng(21268)
    markers = 15
    haplotypes = rng.integers(-1, 2, (2, markers, 2))
    freqs = rng.uniform(.1, .49, markers)
    genetic_map = np.cumsum(rng.uniform(0, 2, markers))
    switches = rng.uniform(.01, .25, (2, markers - 1))
    errors = rng.uniform(.001, .02, (2, markers))
    weights = rng.uniform(.25, 1, markers)
    weights[3] = 0
    origins = rng.integers(-1, 2, (2, markers))
    copies = rng.integers(-1, 2, (2, markers))
    origin_error = rng.uniform(.01, .15, (2, markers))
    copy_error = rng.uniform(.01, .25, (2, markers))
    expected, expected_switches, expected_path = _labeled_reference(
        haplotypes, freqs, genetic_map, switches, errors, weights,
        origins=origins, copies=copies, origin_error=origin_error, copy_error=copy_error,
    )
    result = infer_phased_pair(haplotypes, freqs, genetic_map, switches, errors, weights=weights,
                               parent_origin=origins, parent_haplotype=copies,
                               parent_origin_error=origin_error, parent_haplotype_error=copy_error)
    np.testing.assert_allclose(result.state_posterior, expected, atol=1e-7)
    np.testing.assert_allclose(result.phase_switch_posterior, expected_switches, atol=1e-7)
    np.testing.assert_array_equal(result.state_path, expected_path)
    np.testing.assert_array_equal(result.hard_ibd, (result.state_path & 1) + ((result.state_path >> 1) & 1))


def test_origin_calls_on_different_markers_anchor_relative_phase():
    haplotypes = np.full((2, 3, 2), -1)
    origin = np.array([[0, -1, -1], [-1, -1, 1]], dtype=np.int8)
    result = infer_phased_pair(haplotypes, [.3] * 3, [0.] * 3, 0., .01, parent_origin=origin)
    np.testing.assert_allclose(result.state_posterior[:, 4:].sum(axis=1), (1 - .01) ** 2 + .01 ** 2)
    np.testing.assert_array_equal(result.state_path >> 2, 1)
    np.testing.assert_array_equal(result.phase_switch_posterior, 0)
    np.testing.assert_allclose(result.ibd_posterior, np.tile([.25, .5, .25], (3, 1)))


def test_inconsistent_origin_labels_do_not_force_an_impossible_switch():
    haplotypes = np.full((2, 2, 2), -1)
    origins = [[0, 1], [0, 0]]
    result = infer_phased_pair(haplotypes, [.31, .31], [0., 0.], 0., .01, parent_origin=origins)
    np.testing.assert_array_equal(result.phase_switch_posterior, 0)
    np.testing.assert_array_equal(result.state_posterior[0], result.state_posterior[1])
    np.testing.assert_allclose(result.state_posterior[:, :4].sum(axis=1), .5)
    assert result.state_path[0] == result.state_path[1]


def test_long_conflicting_origin_runs_preserve_rare_paths():
    markers = 700
    haplotypes = np.full((2, markers, 2), -1, dtype=np.int8)
    origins = np.zeros((2, markers), dtype=np.int8)
    origins[0, markers // 2:] = 1
    result = infer_phased_pair(haplotypes, np.full(markers, .3), np.zeros(markers), 0., .01, parent_origin=origins)
    np.testing.assert_array_equal(result.phase_switch_posterior, 0)
    np.testing.assert_allclose(result.state_posterior[:, :4].sum(axis=1), .5, atol=1e-7)
    np.testing.assert_allclose(result.ibd_posterior, np.tile([.25, .5, .25], (markers, 1)))


def test_parental_copy_equality_informs_sharing_and_survives_parent_phase_flips():
    calls = np.full((2, 3, 2), -1)
    origins = np.zeros((2, 3), dtype=np.int8)
    copies = np.array([[0, 0, 1], [0, 1, 1]], dtype=np.int8)
    copy_error = np.array([[.02], [.03]])
    args = (calls, [.31] * 3, [0., 2000., 4000.], .001, .01)
    result = infer_phased_pair(*args, parent_origin=origins, parent_haplotype=copies,
                               parent_haplotype_error=copy_error)
    parity_error = .02 * .97 + .03 * .98
    same = np.array([1 - parity_error, parity_error, 1 - parity_error])
    np.testing.assert_allclose(result.ibd_posterior, np.column_stack((.5 * (1 - same), [.5] * 3, .5 * same)),
                               atol=1e-7)
    np.testing.assert_array_equal((result.state_path >> 1) & 1, [1, 0, 1])

    flipped_copies = copies ^ np.array([0, 1, 1], dtype=np.int8)
    flipped = infer_phased_pair(*args, parent_origin=origins, parent_haplotype=flipped_copies,
                                parent_haplotype_error=copy_error)
    for original, changed in zip(result, flipped):
        np.testing.assert_array_equal(original, changed)


def test_absent_and_uninformative_duo_labels_reproduce_the_raw_result():
    haplotypes = np.array([[[0, 1], [1, 1]], [[1, 0], [0, 1]]])
    args = (haplotypes, [.5, .2], [0., .04], [[.001], [.002]], .01)
    untreated = infer_phased_pair(*args)
    missing = infer_phased_pair(*args, parent_origin=np.full((2, 2), -1), parent_haplotype=[[0, 1], [-1, -1]])
    uninformative = infer_phased_pair(*args, parent_origin=[[0, 1], [1, 0]], parent_origin_error=.5,
                                      parent_haplotype=[[0, 1], [0, 1]], parent_haplotype_error=[[.01], [.5]])
    for ref, no_calls, no_information in zip(untreated, missing, uninformative):
        np.testing.assert_array_equal(ref, no_calls)
        np.testing.assert_array_equal(ref, no_information)


@pytest.mark.parametrize('use_duo', [False, True])
def test_ibd_only_decoder_matches_full_outputs(use_duo):
    rng = np.random.default_rng(109)
    markers = 101
    haplotypes = rng.integers(-1, 2, (2, markers, 2), dtype=np.int8)
    args = (haplotypes, rng.uniform(.1, .5, markers), np.arange(markers) * .07, .002, .01)
    options = dict(weights=rng.uniform(.1, 1., markers))
    if use_duo:
        options.update(parent_origin=rng.integers(-1, 2, (2, markers)),
                       parent_haplotype=rng.integers(-1, 2, (2, markers)))
    full = infer_phased_pair(*args, **options)
    hard, probabilities = infer_phased_pair_ibd(*args, **options)
    hard_only, absent = infer_phased_pair_ibd(*args, posterior=False, **options)
    np.testing.assert_array_equal(hard, full.hard_ibd)
    np.testing.assert_array_equal(hard_only, full.hard_ibd)
    np.testing.assert_array_equal(probabilities, full.ibd_posterior)
    assert absent is None


@pytest.mark.parametrize(('option', 'invalid'), [
    ('parent_origin', [[0], [1]]),
    ('parent_origin', [[0., 0.], [1., .2]]),
    ('parent_haplotype', [[-1, 3], [0, 1]]),
    ('parent_origin_error', 0.),
    ('parent_haplotype_error', .6),
])
def test_invalid_duo_label_inputs(option, invalid):
    options = {'parent_origin': [[0, 1], [1, 0]], option: invalid}
    with pytest.raises(ValueError, match=option):
        infer_phased_pair(np.full((2, 2, 2), -1), [.3, .3], [0., 1.], .01, .01, **options)
