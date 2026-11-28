import numpy as np
import pytest

from snipar.imputation.duohmm import (
    advance_parental_origin,
    call_parental_origin,
    infer_parental_origin,
    recombination_probabilities,
)


def _duos_with_switches():
    rng = np.random.default_rng(25)
    marker = np.arange(1000)
    parent = rng.integers(0, 2, (marker.size, 2), dtype=np.int8)
    child = rng.integers(0, 2, (2, marker.size, 2), dtype=np.int8)
    parent_copy = np.stack((marker >= 400, marker < 280)).astype(np.int8)
    child_origin = np.stack((marker >= 650, marker < 540)).astype(np.int8)
    for i in range(2):
        child[i, marker, child_origin[i]] = parent[marker, parent_copy[i]]

    parent[800:] = parent[800:, ::-1]
    parent_copy[:, 800:] ^= 1
    parent[20] = -11
    child[:, 50] = -11
    parent[120] = 1
    child[:, 120] = 0
    return np.concatenate((parent[None], child)), marker * 2000, child_origin, parent_copy


def test_duos_infer_both_haplotype_indices_and_handle_missing_calls():
    phased, positions, child_origin, parent_copy = _duos_with_switches()
    called_origin, called_copy, _ = call_parental_origin(
        phased, [[0, 1], [0, 2]], recombination_probabilities(positions), threads=2,
    )
    for called, true, haplotypes in (
        (called_origin, child_origin, phased[1:]),
        (called_copy, parent_copy, np.broadcast_to(phased[0], phased[1:].shape)),
    ):
        heterozygous = haplotypes[:, :, 0] != haplotypes[:, :, 1]
        typed = (called >= 0) & heterozygous
        assert (typed.sum(axis=1) / heterozygous.sum(axis=1) > 0.98).all()
        np.testing.assert_array_equal(called[typed], true[typed])
        np.testing.assert_array_equal(called[:, [20, 50, 120]], -np.ones((2, 3), dtype=np.int8))

    np.testing.assert_array_equal(infer_parental_origin(*phased[:2], positions), called_origin[0])


def test_duo_with_only_shared_heterozygotes_has_unresolved_orientation():
    haplotypes = np.tile(np.array([0, 1], dtype=np.int8), (100, 1))
    recombination = recombination_probabilities(np.arange(len(haplotypes)))
    origin, parent_copy, _ = call_parental_origin(np.stack((haplotypes, haplotypes)), [[0, 1]], recombination)
    assert np.all(origin == -1)
    assert np.all(parent_copy == -1)


def test_forward_and_backward_checkpoints_preserve_calls_across_chunk_boundaries():
    phased, positions, _, _ = _duos_with_switches()
    recombination = recombination_probabilities(positions)
    pairs = np.array([[0, 1], [0, 2]])
    whole_origin, whole_copy, whole_beta = call_parental_origin(phased, pairs, recombination)
    boundaries = [0, 1, 137, 400, 542, 650, 800, phased.shape[1]]
    blocks = [slice(left, right) for left, right in zip(boundaries, boundaries[1:])]
    checkpoints = [None]
    for block in blocks[:-1]:
        checkpoints.append(advance_parental_origin(phased[:, block], pairs, recombination[block], checkpoints[-1]))

    beta = None
    for block, alpha in reversed(list(zip(blocks, checkpoints))):
        origin, copy, beta = call_parental_origin(phased[:, block], pairs, recombination[block], alpha, beta)
        np.testing.assert_array_equal(origin, whole_origin[:, block])
        np.testing.assert_array_equal(copy, whole_copy[:, block])
    np.testing.assert_array_equal(beta, whole_beta)


def test_genetic_map_controls_interval_crossover_probabilities():
    positions = [500, 500_500, 1_000_500]
    approximate = recombination_probabilities(positions)
    mapped = recombination_probabilities(positions, [1, 1, 2])
    assert approximate[0] == mapped[0] == mapped[1] == 0
    assert approximate[1] == approximate[2]
    assert approximate[2] < mapped[2]
    assert mapped[2] == pytest.approx((1 - np.exp(-0.02)) / 2)
    with pytest.raises(ValueError, match="nondecreasing"):
        recombination_probabilities(positions, [0, 2, 1])
