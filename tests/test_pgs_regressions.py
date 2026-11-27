import numpy as np
import pytest

from snipar.gtarray import gtarray
from snipar import pgs, pgs_am


def test_compute_two_dimensional_partial_overlap_without_mutating_input():
    scores = pgs.pgs(np.array(['a', 'c']), np.array([2., 3.]),
                     np.array([['A', 'C'], ['T', 'G']]))
    genotypes = np.array([[0., 2., 2.], [2., 0., 0.], [np.nan, 1., 1.]])
    array = gtarray(genotypes, np.array(['x', 'y', 'z']), sid=np.array(['a', 'b', 'c']),
                    alleles=np.array([['A', 'C'], ['A', 'C'], ['G', 'T']]))
    score = scores.compute(array)
    np.testing.assert_allclose(score.gts, [[-5.], [5.], [0.]])
    np.testing.assert_allclose(array.gts.filled(np.nan), genotypes)
    np.testing.assert_array_equal(array.sid, ['a', 'b', 'c'])


def test_parental_score_uses_double_mean_and_preserves_missing(tmp_path):
    scores = pgs.pgs(np.array(['a']), np.array([2.]), np.array([['A', 'C']]))
    genotypes = np.array([[[0.], [3.]], [[2.], [np.nan]]])
    array = gtarray(genotypes, np.array(['x', 'y']), sid=np.array(['a']),
                    alleles=np.array([['A', 'C']]), fams=np.array(['1', '2']))
    output = scores.compute(array, cols=np.array(['proband', 'parental']))
    np.testing.assert_allclose(output.gts.filled(np.nan), [[-2., 2.], [2., np.nan]])
    np.testing.assert_allclose(array.gts.filled(np.nan), genotypes)
    path = tmp_path / 'scores.txt'
    output.write(str(path))
    np.testing.assert_allclose(pgs.read_pgs(path).gts.filled(np.nan), [[-2., 2.], [2., np.nan]])


@pytest.mark.parametrize('header, row', [('FID IID proband', '000 010 1.25'),
                                        ('FID\tIID\tFATHER_ID\tMOTHER_ID\tproband', '000\t010\tNA\t002\tNA')])
def test_read_single_observation_and_single_score_column(tmp_path, header, row):
    path = tmp_path / 'scores.txt'
    path.write_text(header + '\n' + row + '\n')
    scores = pgs.read_pgs(path)
    assert scores.gts.shape == (1, 1)
    np.testing.assert_array_equal(scores.ids, ['010'])
    np.testing.assert_array_equal(scores.fams, ['000'])
    if 'FATHER' in header:
        assert np.ma.is_masked(scores.gts[0, 0])
        np.testing.assert_array_equal(scores.par_status, [[1, 0]])
    else:
        assert scores.gts[0, 0] == 1.25


def test_score_batch_add_matches_by_id_preserves_order_and_uses_observed_batch():
    first = pgs.pgarray(np.array([[1.], [np.nan], [np.nan]]), np.array(['z', 'x', 'y']),
                        sid=np.array(['proband']))
    second = pgs.pgarray(np.array([[2.], [np.nan], [3.]]), np.array(['x', 'y', 'z']),
                         sid=np.array(['proband']))
    total = first.add(second)
    np.testing.assert_array_equal(total.ids, ['z', 'x', 'y'])
    np.testing.assert_allclose(total.gts.filled(np.nan), [[4.], [2.], [np.nan]])


def test_parsum_regression_does_not_drop_rows_with_missing_unused_scores():
    rng = np.random.default_rng(12)
    fams = np.repeat(np.arange(30).astype(str), 2)
    ids = np.arange(len(fams)).astype(str)
    scores = rng.normal(size=(len(ids), 4))
    scores[:, 3] = np.nan  # model 2 doesn't use grandparent scores
    y = (scores[:, 0] + .2*(scores[:, 1] + scores[:, 2]) + rng.normal(size=len(ids))
         + np.repeat(rng.normal(size=30), 2))
    score_array = pgs.pgarray(scores, ids, sid=np.array(['proband', 'paternal', 'maternal', 'gpp']), fams=fams)
    phenotype = gtarray(y[:, None], ids)
    fit = pgs.fit_pgs_model(phenotype, score_array, 2, parsum=True, return_fit=True)
    np.testing.assert_array_equal(fit['ids'], ids)
    np.testing.assert_array_equal(fit['cols'], ['intercept', 'proband', 'parental'])
    np.testing.assert_array_equal(score_array.sid, ['proband', 'paternal', 'maternal', 'gpp'])
    np.testing.assert_allclose(score_array.gts.filled(np.nan), scores)


def test_simulate_am_summary_uses_sampling_covariance(monkeypatch):
    estimates_and_se = tuple(value for i in range(7) for value in
                             (np.array([i, i + 2.]), np.array([.2, .2])))
    observed = []
    def simulate(*args):
        observed.extend(args)
        return estimates_and_se
    monkeypatch.setattr(pgs_am, 'simulate_ests', simulate)
    covariance = np.eye(3)
    output = pgs_am.check_se_calc(2, .1, .2, .3, covariance, .4, .01, .2, .01)
    assert observed[4] is covariance
    np.testing.assert_allclose(output[::2], np.arange(1, 8.))
    np.testing.assert_allclose(output[1::2], .2)


def test_grouped_pgs_likelihood_matches_dense_family_covariances():
    rng = np.random.default_rng(71)
    counts = np.array([[0, 2], [1, 1], [1, 2], [2, 0], [3, 1], [4, 2]])
    group_sizes = np.arange(2, len(counts) + 2)
    moments = []
    dense = []
    for (ns, n_parents), nx in zip(counts, group_sizes):
        siblings, parents = rng.normal(size=(nx, ns)), rng.normal(size=(nx, n_parents))
        sr = siblings.sum(axis=1) / np.sqrt(max(ns, 1))
        pr = parents.sum(axis=1) / np.sqrt(max(n_parents, 1))
        moments.append([sr @ sr, pr @ pr, sr @ pr,
                        np.sum(siblings ** 2) - sr @ sr, np.sum(parents ** 2) - pr @ pr])
        dense.append(np.column_stack((siblings, parents)))
    for r in [-.99, -.6, 0., .6, .99]:
        expected = 0.
        for (ns, n_parents), sample in zip(counts, dense):
            sibling = np.arange(ns + n_parents) < ns
            cov = pgs.pgs_corr_matrix(r, sibling, ~sibling)
            expected += len(sample) * np.linalg.slogdet(cov)[1]
            expected += np.sum(sample.T * np.linalg.solve(cov, sample.T))
        np.testing.assert_allclose(pgs._pgs_profile_loglik(r, counts, group_sizes, np.array(moments)),
                                   expected, rtol=1e-12)


def test_parental_only_correlation_counts_each_family_once_with_sibling_column():
    # One family has many siblings and shouldn't receive proportionally more weight.
    parents = np.array([[-1.4, 2.], [.2, -.3], [1.2, .6], [1.7, 1.]])
    sizes = [6, 1, 3, 1]
    family = np.repeat(np.arange(len(sizes)).astype(str), sizes)
    n = len(family)
    columns = ['proband', 'sibling', 'paternal', 'maternal']
    genotypes = np.column_stack((np.arange(n), np.ones(n), np.repeat(parents, sizes, axis=0)))
    scores = pgs.pgarray(genotypes, np.arange(n).astype(str), sid=np.array(columns),
                         fams=family, par_status=np.zeros((n, 2), dtype=int))
    observed, se, family_sizes = scores.estimate_r(parents_only=True)
    expected = np.corrcoef(parents.T)[0, 1]
    np.testing.assert_allclose(observed, expected)
    np.testing.assert_allclose(se, (1 - expected**2) / np.sqrt(len(sizes) - 1))
    assert list(family_sizes.values()) == sizes


def test_parental_sum_adjustment_does_not_require_separate_parent_score_columns(monkeypatch):
    scores = pgs.pgarray(np.array([[1., 2.], [2., 3.]]), np.array(['x', 'y']),
                         sid=np.array(['proband', 'parental']), fams=np.array(['f1', 'f2']),
                         par_status=np.ones((2, 2), dtype=int))
    monkeypatch.setattr(scores, 'estimate_r', lambda: (.4, .1, {'f1': 1, 'f2': 1}))
    scores.am_adj()
    np.testing.assert_allclose(scores.gts[:, 1], np.array([2., 3.]) * 1.4)


def test_assortative_mating_profile_handles_zero_population_effect():
    uncontested = pgs_am.estimate_rho(.25, .01, .45, .02, .21)
    actual = pgs_am.alpha_from_beta(.25, .01, 0., .45, .02, .21, return_all=False)
    np.testing.assert_allclose(actual, -1 / (1 + uncontested['r']))
    assert np.isnan(pgs_am.estimate_r(0., .01, .45, .02, 0., return_all=False))


def test_batched_pgs_skips_weighted_snp_absent_from_imputed_genotypes(tmp_path):
    import h5py
    from pysnptools.snpreader import Bed, SnpData

    ids = np.array(['child1', 'child2', 'father', 'mother'])
    snps = np.array(['a', 'b', 'c'])
    bedfile = tmp_path / 'observed.bed'
    Bed.write(str(bedfile), SnpData(iid=np.column_stack((ids, ids)), sid=snps,
                                    pos=np.column_stack((np.ones(3), np.zeros(3), np.arange(1, 4))),
                                    val=np.array([[0., 2., 1.], [2., 0., 0.],
                                                  [2., 1., 0.], [1., 2., 1.]])), count_A1=True)
    bim = np.loadtxt(bedfile.with_suffix('.bim'), dtype=str, usecols=(1, 4, 5))
    imputed = tmp_path / 'imputed.hdf5'
    with h5py.File(imputed, 'w') as handle:
        handle['pedigree'] = np.array([['FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'has_father', 'has_mother'],
                                        ['family1', 'child1', 'father', 'mother', 'True', 'True'],
                                        ['family1', 'child2', 'father', 'mother', 'True', 'True']], dtype='S')
        handle['families'] = np.array(['family1'], dtype='S')
        handle['bim_columns'] = np.array(['rsid', 'allele1', 'allele2'], dtype='S')
        handle['bim_values'] = bim[[2, 0]].astype('S')
        handle['imputed_par_gts'] = [[0., 0.]]
    definition = pgs.pgs(snps, np.array([2., 50., 4.]), bim[:, 1:])
    args = dict(bedfile=str(bedfile), par_gts_f=imputed, verbose=False)
    all_scores = pgs.compute(definition, **args)
    batches = pgs.compute(definition, batch_size=1, **args)
    np.testing.assert_array_equal(batches.ids, all_scores.ids)
    np.testing.assert_allclose(batches.gts, all_scores.gts)
    np.testing.assert_allclose(all_scores.gts,
                               [[0., 0., 2.], [0., 0., 2.]])
