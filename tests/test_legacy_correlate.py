import numpy as np
import pytest

from snipar import correlate


def _summary(n, chrom=1, ldscores=None):
    effects = np.linspace(-1., 1., n)
    return correlate.sumstats(
        chrom, np.array([f'{chrom}:{i}' for i in range(n)]), np.arange(n),
        np.repeat('A', n), np.repeat('C', n), np.full(n, .3),
        effects, np.full(n, .01), effects / 2, np.full(n, .01),
        effects * 1.5, np.full(n, .01), np.full(n, .2), np.full(n, .4),
        ldscores=ldscores,
    )


def test_legacy_sumstats_filter_and_concatenate_keep_columns_aligned():
    data = _summary(6, ldscores=np.ones(6))
    data.concatenate(_summary(5, chrom=2))
    assert len(data.ldscores) == len(data.sid) == 11
    assert data.ldscores.mask[6:].all()
    data.direct_SE[1] = np.ma.masked
    data.filter_NAs()
    data.filter_corrs(6)
    np.testing.assert_array_equal(data.sid, ['1:0', '1:2', '1:3', '1:4', '1:5'])


def test_legacy_jackknife_accepts_masked_summary_statistics():
    rng = np.random.default_rng(74)
    z1, other = rng.normal(size=(2, 60))
    z2 = .5*z1 + other
    v = np.full(len(z1), .01)
    zeros = np.zeros(len(z1))
    ld = rng.uniform(1, 3, size=len(z1))
    mask = np.arange(len(z1)) % 11 == 0
    actual = correlate.jacknife_est(np.ma.array(z1, mask=mask), z2, v, v, zeros, ld, 4)
    expected = correlate.jacknife_est(z1[~mask], z2[~mask], v[~mask], v[~mask], zeros[~mask], ld[~mask], 4)
    for result, reference in zip(actual, expected):
        np.testing.assert_allclose(result, reference)
    with pytest.raises(ValueError, match='n_blocks'):
        correlate.jacknife_est(z1, z2, v, v, zeros, ld, 1)
