import h5py
import numpy as np
import pandas as pd
import pytest

from snipar.ibd_output import Relatedness
from snipar.ibd_joint import joint_pair_groups
from snipar.scripts.ibd import _pair_batches


def test_posterior_output_and_length_weighted_genome_relatedness(tmp_path):
    pairs = np.array([['s1', 's2'], ['s3', 's4'], ['s5', 's6']])
    output = Relatedness(pairs)
    prefix = str(tmp_path / 'chr_@')
    snps = np.array(['a', 'b', 'c'])
    positions = np.array([10, 1_000_010, 2_000_010])
    first = np.tile(np.array([0, 1, 0], dtype=np.float32), (1, 3, 1))
    second = np.tile(np.array([1, 0, 0], dtype=np.float32), (1, 3, 1))
    third = np.tile(np.array([0, 0, 1], dtype=np.float32), (1, 3, 1))
    with output.chromosome_batch(np.array([0, 2, 1]), prefix) as store:
        store(1, pairs[:1], snps, positions, np.array([0, 1, 3]), first, np.array([0]))
        store(1, pairs[1:2], snps, positions, np.array([0, 1, 3]), second, np.array([2]))
    with output.chromosome_batch(np.array([0]), prefix) as store:
        store(2, pairs[:1], snps, np.array([10, 2_000_010, 3_000_010]), np.array([0, .7, 1]),
              third, np.array([0]))
    path = tmp_path / 'chr_1.ibd.posterior.hdf5'
    with h5py.File(path) as f:
        assert f.attrs['format_version'] == 1
        assert f.attrs['chromosome'] == '1'
        np.testing.assert_array_equal(f['variant_ids'].asstr()[:], snps)
        np.testing.assert_array_equal(f['pairs'].asstr()[:], pairs[:2])
        np.testing.assert_array_equal(f['positions'][:], positions)
        np.testing.assert_array_equal(f['probabilities'][:], np.concatenate([first, second]))
    stats = pd.read_csv(output.write(prefix), sep='\t')
    np.testing.assert_array_equal(stats['N_CHR'], [2, 1])
    np.testing.assert_allclose(stats['PropIBD_cM'], [.625, 0])
    np.testing.assert_allclose(stats['PropIBD_Mb'], [.8, 0])
    np.testing.assert_array_equal(stats['ID1'], ['s1', 's3'])


def test_failed_posterior_write_keeps_previous_output(tmp_path):
    target = tmp_path / 'example_chr_1.ibd.posterior.hdf5'
    target.write_bytes(b'previous')
    pairs = np.array([['s1', 's2']])
    with pytest.raises(RuntimeError, match='halt'):
        with Relatedness(pairs).chromosome_batch(np.array([0]), str(tmp_path / 'example')) as store:
            store(1, pairs, np.array(['a']), np.array([1]), np.array([0.]),
                  np.array([[[.25, .5, .25]]]), np.array([0]))
            raise RuntimeError('halt')
    assert target.read_bytes() == b'previous'
    assert list(tmp_path.iterdir()) == [target]


def test_joint_clique_stays_together_across_argument_batches():
    pairs = np.array([('a', 'b'), ('d', 'e'), ('a', 'c'), ('f', 'g'), ('b', 'c')])
    batches = _pair_batches(pairs, 3, joint_groups=joint_pair_groups(pairs))
    assert any({0, 2, 4} <= set(rows) for rows in batches)
    assert sorted(np.concatenate(batches)) == list(range(len(pairs)))
