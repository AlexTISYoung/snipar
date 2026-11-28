from itertools import product
from types import SimpleNamespace

import numpy as np
import pytest
from pysnptools.snpreader import Bed, SnpData

from snipar import ibd
from snipar.map import decode_map_from_pos, get_map_positions


def test_ibd_marginals_match_enumerated_hidden_paths():
    genotypes = np.array([[0., 2., np.nan, 0., 1.], [1., 0., 2., 0., 2.]])
    freqs = np.array([.1, .31, .4, .18, .3])
    positions = np.array([.9, .9, 1.5, 2.6, 8.2])
    weights = np.array([1., 0., .6, 1., .7])
    errors = np.array([0., 0., .02, .01, 0.])
    pairs = np.array([[0, 1], [1, 0]])
    emission = np.ones((len(freqs), 3))
    for marker, (a, b) in enumerate(genotypes.T):
        if weights[marker] and np.isfinite(a + b):
            emission[marker] = np.exp(weights[marker] * ibd.p_obs_given_IBD(
                int(a), int(b), freqs[marker], errors[marker],
            ))
    transitions = [np.exp(ibd.transition_matrix(d)) for d in np.diff(positions)]
    expected = np.zeros_like(emission)
    for path in product(range(3), repeat=len(freqs)):
        probability = (.25, .5, .25)[path[0]] * emission[0, path[0]]
        for marker, (a, b) in enumerate(zip(path, path[1:])):
            probability *= transitions[marker][a, b] * emission[marker + 1, b]
        for marker, state in enumerate(path):
            expected[marker, state] += probability
    expected /= expected.sum(axis=1, keepdims=True)

    actual = ibd.infer_ibd_posteriors(pairs, genotypes, freqs, positions, weights, errors)
    assert actual.dtype == np.float32
    np.testing.assert_allclose(actual[0], expected, rtol=2e-7)
    np.testing.assert_array_equal(actual[0], actual[1])


def test_ibd_marginals_remain_finite_when_forward_and_backward_evidence_underflows():
    size = 4096
    genotypes = np.zeros((2, size), dtype=np.int8)
    genotypes[1, 0] = 2
    # Only IBD0 can explain the first observation with no genotyping error.
    # The remaining markers strongly favor IBD2 but share its first map position.
    probabilities = ibd.infer_ibd_posteriors(
        np.array([[0, 1]]), genotypes, np.full(size, .8), np.zeros(size), np.ones(size), np.zeros(size),
    )
    np.testing.assert_array_equal(probabilities, np.tile([1., 0., 0.], (1, size, 1)))


def test_missing_genotypes_leave_ibd_prior_stationary():
    size = 18
    args = (np.full(size, .13), np.arange(size, dtype=float) * 10, np.ones(size), np.zeros(size))
    actual = ibd.infer_ibd_posteriors(np.array([[0, 1]]), np.full((2, size), np.nan), *args)
    np.testing.assert_array_equal(actual, np.tile([.25, .5, .25], (1, size, 1)))


def test_ibd_map_extends_sparse_anchors_with_endpoint_slopes(tmp_path):
    filename = tmp_path / 'chr.map'
    filename.write_text('pposition gposition\n100 1.0\n200 2.5\n400 3.5\n')
    positions = np.array([0, 130, 300, 500])
    np.testing.assert_allclose(ibd._ibd_map(None, 1, positions, filename), [-.5, 1.45, 3., 4.])
    np.testing.assert_allclose(get_map_positions(filename, SimpleNamespace(pos=positions), min_map_prop=None),
                               [1., 1.45, 3., 3.5])

    positions = np.array([17_076_000, 17_076_255, 17_082_000, 51_000_000])
    default = decode_map_from_pos(22, positions)
    extrapolated = decode_map_from_pos(22, positions, extrapolate=True)
    assert np.isnan(default[[0, 3]]).all()
    np.testing.assert_array_equal(default[1:3], extrapolated[1:3])
    assert np.isfinite(extrapolated).all() and np.all(np.diff(extrapolated) >= 0)


def test_ibd_posterior_batches_keep_pair_indices_and_sort_physical_positions(tmp_path, monkeypatch):
    ids = np.array(['one', 'two', 'three', 'four'])
    snps = np.array(['s1', 's2', 's3', 's4'])
    bp = np.array([7, 3, 3, 6])
    gts = np.array([[1., 0., 2., 1.], [0., 1., 1., 0.], [2., 1., 0., 1.], [1., 2., 0., 2.]])
    prefix = tmp_path / 'observed'
    Bed.write(str(prefix), SnpData(iid=np.column_stack((ids, ids)), sid=snps, val=gts,
                                  pos=np.column_stack((np.ones(len(bp)), bp * .17, bp))), count_A1=True)
    pairs = np.array([['one', 'two'], ['absent', 'four'], ['three', 'four']])
    batches = []
    original_plan = ibd.plan_batches

    def cap_pairs(items, item_bytes, **kwargs):
        if items == 2:
            kwargs['batch_size'] = 1
        return original_plan(items, item_bytes, **kwargs)

    monkeypatch.setattr(ibd, 'plan_batches', cap_pairs)

    def capture(chrom, batch_pairs, markers, positions, cM, posterior, pair_indices):
        assert chrom == 1
        np.testing.assert_array_equal(batch_pairs, pairs[pair_indices])
        np.testing.assert_array_equal(markers, ['s2', 's3', 's4', 's1'])
        np.testing.assert_array_equal(positions, [3, 3, 6, 7])
        np.testing.assert_allclose(cM, [.51, .51, 1.02, 1.19])
        np.testing.assert_allclose(posterior.sum(axis=2), 1., rtol=2e-7)
        assert posterior.shape == (1, 4, 3)
        batches.extend(pair_indices.tolist())

    ibd.infer_ibd_chr(pairs, .001, None, str(tmp_path / 'ibd'), bedfile=str(prefix), chrom=1,
                      posterior_callback=capture, num_threads=2)
    assert batches == [0, 2]


@pytest.mark.parametrize('label, expected', [('1.0', 1), ('chr1', 1), ('scaffold', 'scaffold')])
def test_ibd_infers_chromosome_labels_for_map_and_outputs(tmp_path, label, expected):
    ids = np.array([['one', 'one'], ['two', 'two']])
    prefix = tmp_path / 'observed'
    values = np.array([[0., 2.], [1., 0.]])
    Bed.write(str(prefix), SnpData(iid=ids, val=values, sid=np.array(['rs1', 'rs2']),
                                  pos=np.array([[1, 0, 1000], [1, 0, 2000]])), count_A1=True)
    bim = prefix.with_suffix('.bim')
    lines = bim.read_text().splitlines()
    bim.write_text(''.join(f'{label}\t' + '\t'.join(line.split()[1:]) + '\n' for line in lines))
    mapfile = None
    if label == 'scaffold':
        mapfile = tmp_path / 'scaffold.map'
        mapfile.write_text('pposition gposition\n1000 1\n2000 2\n')

    ibd.infer_ibd_chr(ids[:, 0].reshape(1, 2), .001, None, str(tmp_path / 'ibd'),
                      bedfile=str(prefix), mapfile=mapfile, num_threads=1)
    assert (tmp_path / f'ibd_chr_{expected}.ibd.segments.gz').is_file()


def test_mixed_pair_and_joint_ibd_preserves_original_edge_order():
    from snipar.ibd_joint import infer_joint_family

    gts = np.random.default_rng(102).integers(3, size=(5, 15), dtype=np.int8)
    rows = np.array([[0, 1], [4, 3], [2, 0], [1, 2]])
    marker = np.arange(gts.shape[1])
    freqs = np.full(len(marker), .3)
    cM = marker * .9
    weights = np.ones(len(marker))
    error = np.full(len(marker), .01)
    transitions, emissions = ibd._chromosome_hmm(freqs, cM, weights, error)
    tables = ibd._posterior_tables(transitions, emissions)
    labels, groups = ibd._joint_families(rows)
    np.testing.assert_array_equal(labels, [0, -1, 0, 0])
    filtered_labels, _ = ibd._joint_families(rows, [True, False, True, True])
    np.testing.assert_array_equal(filtered_labels, labels)
    excluded_labels, excluded_groups = ibd._joint_families(rows, [False, False, False, False])
    np.testing.assert_array_equal(excluded_labels, [-1, -1, -1, -1])
    assert not excluded_groups
    # The batch ends before the final two edges of the joint family.
    paths, posterior = ibd._infer_joint_block(
        rows[:2], labels[:2], groups, gts, marker, transitions, emissions, 0, tables, freqs, cM, weights, error,
    )
    joint_paths, joint_probs = infer_joint_family(gts[:3], rows[[0]], freqs, cM, weights, error, posterior=True)
    pair_paths, pair_probs = ibd._infer_block(rows[[1]], gts, marker, transitions, emissions, 0, tables)
    np.testing.assert_array_equal(paths, np.concatenate((joint_paths, pair_paths)))
    np.testing.assert_array_equal(posterior, np.concatenate((joint_probs, pair_probs)))
