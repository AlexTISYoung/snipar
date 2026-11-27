import numpy as np
import pytest

from snipar.ld import compute_ld_scores, compute_ld_scores_packed, pack_ld_hardcalls


@pytest.mark.parametrize('samples', [2, 63, 64, 65, 155])
def test_packed_ld_matches_pairwise_complete_corrected_ld(samples):
    rng = np.random.default_rng(122)
    calls = rng.integers(0, 3, (samples, 37), dtype=np.int8)
    calls[:, 1] = 1
    calls[:, 2] = -127
    calls[rng.random(calls.shape) < .2] = -127
    calls[:, 4] = rng.integers(0, 3, samples, dtype=np.int8)
    calls[:, 5] = rng.integers(0, 3, samples, dtype=np.int8)
    packed = np.empty((calls.shape[1], 3, (samples + 63) // 64), dtype=np.uint64)

    # Independently decoded blocks write directly into their workspace slices.
    pack_ld_hardcalls(calls[:, :13], packed[:13])
    pack_ld_hardcalls(np.asfortranarray(calls[:, 13:]), packed[13:])
    columns = np.array([4, 2, 1, 5, 33, 14, 17, 13, 0, 22, 35])
    positions = np.array([0., 0., .5, .6, 1., 1., 1.5, 1.6, 2., 3., 5.])
    as_float = calls.astype(np.float32)
    as_float[calls < 0] = np.nan
    np.testing.assert_array_equal(compute_ld_scores_packed(packed, positions, 1., columns),
                                  compute_ld_scores(as_float, positions, 1., columns))


def test_equal_called_counts_can_have_different_missingness_patterns():
    calls = np.array([
        [0, 2, -127],
        [1, 0, 1],
        [1, 0, 2],
        [2, 1, 1],
        [2, -127, 0],
        [-127, 2, 0],
        [-127, -127, -127],
    ], dtype=np.int8)
    packed = np.empty((3, 3, 1), dtype=np.uint64)
    pack_ld_hardcalls(calls, packed)
    positions = np.array([.1, .2, .3])
    np.testing.assert_array_equal(compute_ld_scores_packed(packed, positions), compute_ld_scores(calls, positions))
