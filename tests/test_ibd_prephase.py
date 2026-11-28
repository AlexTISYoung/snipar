"""Exercise phased IBD and reuse its duoHMM calls through the normal entry points."""

import h5py
import numpy as np
import pytest

from snipar.config import nan_integer
from snipar.ibd_prephase import _unpack, observed_parent_duos, pack_phased
from snipar.scripts import ibd, impute


def test_packed_haplotypes_preserve_phase_and_partial_missing_calls():
    calls = np.array([[[0, 1], [nan_integer, 1], [0, nan_integer]],
                      [[1, 0], [1, 1], [nan_integer, nan_integer]]], dtype=np.int8)
    packed = np.empty((2, 2), dtype=np.uint8)
    pack_phased(calls, packed)
    np.testing.assert_array_equal(_unpack(packed, np.array([1, 0]), np.array([2, 0, 1]), nan_integer),
                                  calls[[1, 0]][:, [2, 0, 1]])


def test_observed_parent_duos_require_exactly_one_genotyped_parent():
    ped = [['one', 'c1', 'p1', 'm1'], ['two', 'c2', 'p1', 'p2'],
           ['zero', 'c3', 'p3', 'm3'], ['singleton', 'c4', '0', 'p2']]
    ids = np.array(['p2', 'outside', 'c2', 'c1', 'p1', 'c4', 'c3', '0'])
    children = np.array([1, 2, 3, 5, 6])
    np.testing.assert_array_equal(observed_parent_duos(ped, ids, children), [[4, 3], [0, 5]])


def test_phased_ibd_reuses_origin_in_imputation(tmp_path):
    pedigree = tmp_path / 'pedigree.txt'
    pedigree.write_text('FID IID FATHER_ID MOTHER_ID\nf c1 father mother\nf c2 father mother\n'
                        'g c3 f2 m2\ng c4 f2 m2\n')
    positions = np.arange(101, 181)
    vcf = tmp_path / 'genotypes.vcf'
    records = [
        '##fileformat=VCFv4.3\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfather\tc1\tc2\tc3\tc4\n',
    ]
    for marker, bp in enumerate(positions):
        calls = ['0|0' if marker % 4 == 0 else '0|1',
                 '0|1' if marker < 40 else '1|0', '0|1', '1|0', '0|1']
        if marker == 0:
            calls = ['0|0'] * len(calls)
        records.append(f'1\t{bp}\ts{marker}\tA\tG\t.\tPASS\t.\tGT\t' + '\t'.join(calls) + '\n')
    vcf.write_text(''.join(records))
    genetic_map = tmp_path / 'map.txt'
    genetic_map.write_text('pposition gposition\n110 0\n170 0.001\n')
    common = ['--vcf', str(vcf), '--pedigree', str(pedigree), '--threads', '2']
    for name, extra in [('phase', []), ('raw', ['--no-duohmm'])]:
        np.random.seed(1407)
        ibd.main(ibd.parser.parse_args([
            *common, '--phased-ibd', '--p_error', '0.001', '--map', str(genetic_map), '--posterior',
            '--out', str(tmp_path / f'{name}_@'), *extra,
        ]))

    origin_path = tmp_path / 'phase_1.ibd.origin.hdf5'
    with h5py.File(origin_path) as origin:
        np.testing.assert_array_equal(origin['pairs'][:], [[b'father', b'c1'], [b'father', b'c2']])
        np.testing.assert_array_equal(origin['positions'][:], positions)
        np.testing.assert_array_equal(origin['origin'][0, [30, 50]], [0, 1])
        np.testing.assert_array_equal(origin['parent_copy'][0, [30, 50]], [0, 0])
    assert not (tmp_path / 'raw_1.ibd.origin.hdf5').exists()
    with (h5py.File(tmp_path / 'phase_1.ibd.posterior.hdf5') as prephased,
          h5py.File(tmp_path / 'raw_1.ibd.posterior.hdf5') as raw):
        assert len(prephased['positions']) == len(positions) - 1
        np.testing.assert_array_equal(prephased['pairs'][:], raw['pairs'][:])
        unanchored = np.flatnonzero(prephased['pairs'][:, 0] == b'c3')[0]
        np.testing.assert_array_equal(prephased['probabilities'][unanchored], raw['probabilities'][unanchored])

    for suffix in ['segments.gz', 'posterior.hdf5']:
        (tmp_path / f'local_1.ibd.{suffix}').symlink_to(tmp_path / f'phase_1.ibd.{suffix}')
    for name in ['phase', 'local']:
        impute.main(impute.parser.parse_args([
            *common, '--ibd', str(tmp_path / f'{name}_@.ibd'), '--parent_origin_map', str(genetic_map),
            '--chunks', '3', '-silent_progress', '--out', str(tmp_path / f'imputed_{name}_@'),
        ]))
    with (h5py.File(tmp_path / 'imputed_phase_1.hdf5') as reused,
          h5py.File(tmp_path / 'imputed_local_1.hdf5') as local):
        assert set(reused) == set(local)
        for key in reused:
            np.testing.assert_array_equal(reused[key][:], local[key][:], err_msg=key)
        assert np.isfinite(reused['robust_imputed_par_gts'][:]).any()


def test_unphased_vcf_rejected_in_phased_ibd_mode(tmp_path):
    vcf = tmp_path / 'genotypes.vcf'
    vcf.write_text('##fileformat=VCFv4.3\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n'
                   '1\t100\ts1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n')
    pedigree = tmp_path / 'pedigree.txt'
    pedigree.write_text('f c1 father mother\nf c2 father mother\n')
    with pytest.raises(ValueError, match='phas'):
        ibd.main(ibd.parser.parse_args([
            '--phased-ibd', '--p_error', '0.001', '--vcf', str(vcf), '--pedigree', str(pedigree),
            '--out', str(tmp_path / 'result'), '--threads', '1',
        ]))
