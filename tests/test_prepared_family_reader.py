import h5py
import numpy as np
import pytest
from bed_reader import to_bed

import snipar.read as read
from snipar.read.family import IdIndex, _read_imputed
from snipar.read.variants import NoMatchingVariantsError


def _bed(path, ids, variants, values):
    to_bed(path, np.asarray(values, dtype=np.float32), count_A1=True, properties=dict(
        iid=ids, sid=variants, chromosome=['1'] * len(variants), bp_position=np.arange(len(variants)) + 1,
        allele_1=['A'] * len(variants), allele_2=['G'] * len(variants),
    ))


def test_prepared_reader_reuses_plan_and_reads_only_relevant_relatives(tmp_path, monkeypatch):
    path = tmp_path / 'family.bed'
    ids = ['sib', 'unrelated', 'mother', 'child', 'father', 'not_a_sibling']
    values = [[0, 2, np.nan], [2, 2, 2], [1, 2, 1], [2, 0, 1], [0, 1, np.nan], [1, 1, 1]]
    _bed(path, ids, ['z', 'a', 'b'], values)
    pedigree = np.array([['family', 'child', 'father', 'mother'], ['family', 'sib', 'father', 'mother']])
    with read.prepare_gts_matrix(bedfile=path, ped=pedigree, ids=['child'], sib=True) as reader:
        plan = reader.sample_plans
        assert reader.group_shapes == ((1, 4),)
        assert reader.n_observed == 4
        variants = reader.match_variants(match_snp_ids=True)
        np.testing.assert_array_equal(variants.sid, ['z', 'a', 'b'])
        chunks = [reader.read(variants.subset(slice(i, i + 1))).gts.data for i in range(3)]
        order_f = reader.read(variants, order='F').gts.data
        np.testing.assert_array_equal(np.concatenate(chunks, axis=2), order_f)
        np.testing.assert_array_equal(order_f[0], np.asarray(values)[[3, 0, 4, 2]])
    monkeypatch.setattr(read, '_sample_plan', lambda *args, **kwargs: pytest.fail('Rebuilt sample plan'))
    with read.prepare_gts_matrix(bedfile=path, sample_plans=plan) as reader:
        assert reader.sample_plans is plan
        np.testing.assert_array_equal(reader.read(variants).gts.data, order_f)


def test_variant_metadata_can_be_read_without_building_family_plan(tmp_path):
    path = tmp_path / 'duplicates.bed'
    _bed(path, ['a', 'b'], ['z', 'dup', 'a', 'dup'], [[0, 1, 2, 1], [2, 2, 2, 0]])
    with read.prepare_gts_matrix(bedfile=path) as reader:
        np.testing.assert_array_equal(reader.observed_variants.sid, ['z', 'dup', 'a', 'dup'])
        np.testing.assert_array_equal(reader.match_variants().sid, ['a', 'z'])
        np.testing.assert_array_equal(reader.match_variants(match_snp_ids=True).sid, ['z', 'a'])
        with pytest.raises(NoMatchingVariantsError):
            reader.match_variants(['dup'])


def test_sample_id_index_retains_last_duplicate_and_uses_int32():
    lookup = IdIndex(np.array(['12', '13', '02', '12', '14', '02']))
    rows = lookup.lookup(['02', '00', '14', '13', '12', '19'])
    assert rows.dtype == np.int32
    np.testing.assert_array_equal(rows, [5, -1, 4, 1, 3, -1])


def test_king_duplicate_filter_uses_chunks_and_preserves_leading_zeros(tmp_path, monkeypatch):
    path = tmp_path / 'king.seg'
    path.write_text('ID1\tID2\tInfType\n003\t002\tFS\n001\t002\tDup/MZ\n004\t003\tDup/MZ\n')
    read_csv = read.pd.read_csv

    def small_chunks(*args, **kwargs):
        assert kwargs.get('chunksize')
        return read_csv(*args, **dict(kwargs, chunksize=2))

    monkeypatch.setattr(read.pd, 'read_csv', small_chunks)
    ids = np.array(['004', '003', '001', '1', '002'])
    np.testing.assert_array_equal(read._exclude_duplicate_twins(ids, path), ['003', '1', '002'])


def test_h5_gather_handles_sparse_rows_and_reordered_duplicate_variants(tmp_path):
    values = np.arange(120, dtype=np.float64).reshape(12, 10)
    path = tmp_path / 'imputed.h5'
    with h5py.File(path, 'w') as handle:
        dataset = handle.create_dataset('calls', data=values, chunks=(4, 3))
        for rows in ([1, 2, 3], [0, 2, 6, 9, 11], []):
            for columns in ([1, 4, 8], [3, 4, 5], [9, 1, 1], []):
                actual = _read_imputed(dataset, np.array(rows, dtype=int), np.array(columns, dtype=int))
                assert actual.dtype == np.float32
                np.testing.assert_array_equal(actual, values[np.ix_(rows, columns)])


def test_prepared_reader_threads_preserve_family_design(tmp_path):
    rng = np.random.default_rng(37)
    families, loci = 1400, 120
    labels = np.array([f'{i:05}' for i in range(families)])
    children_a = np.char.add('A', labels)
    children_b = np.char.add('B', labels)
    dads = np.char.add('D', labels)
    moms = np.char.add('M', labels)
    ids = np.concatenate((children_a, children_b, dads, moms))
    values = rng.integers(0, 3, size=(len(ids), loci), dtype=np.int8).astype(np.float32)
    values[:4, :2] = np.nan
    path = tmp_path / 'parallel.bed'
    _bed(path, ids, [f'rs{i:04}' for i in range(loci)], values)
    ped = np.concatenate((np.column_stack((labels, children_a, dads, moms)),
                          np.column_stack((labels, children_b, dads, moms))))
    with read.prepare_gts_matrix(bedfile=path, ped=ped, sib=True, num_threads=1) as serial:
        variants = serial.match_variants(match_snp_ids=True)
        plan = serial.sample_plans
        expected = serial.read(variants).gts
    with read.prepare_gts_matrix(bedfile=path, sample_plans=plan, num_threads=2) as parallel:
        actual = parallel.read(variants).gts
    assert expected.size >= 1_000_000  # Exercise the parallel assembly threshold.
    np.testing.assert_array_equal(actual.data, expected.data)
    np.testing.assert_array_equal(actual.mask, expected.mask)


def test_prepared_reader_imputes_flipped_alleles_and_keeps_robust_counts(tmp_path):
    bed = tmp_path / 'observed.bed'
    ids = ['child2', 'mother', 'dad', 'child1', 'child3']
    vals = np.array([[1, 0, 2], [2, 0, 1], [0, 2, 1], [2, 1, 0], [0, 2, 1]])
    _bed(bed, ids, ['b', 'c', 'a'], vals)
    imputation = tmp_path / 'imputed.h5'
    with h5py.File(imputation, 'w') as handle:
        handle['pedigree'] = np.array([
            ['FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'has_father', 'has_mother'],
            ['fam2', 'child2', 'missing-father', 'mother', 'False', 'True'],
            ['fam3', 'child3', 'other-father', 'other-mother', 'False', 'False'],
            ['fam1', 'child1', 'dad', 'mother', 'True', 'True'],
        ], dtype='S')
        handle['families'] = np.array(['fam3', 'unused', 'fam2'], dtype='S')
        handle['bim_columns'] = np.array(['id', 'allele1', 'allele2'], dtype='S')
        handle['bim_values'] = np.array([['a', 'A', 'G'], ['b', 'G', 'A']], dtype='S')
        handle['imputed_par_gts'] = [[.25, .75], [0., 0.], [1.25, .50]]
        handle['num_observed_parental_alleles'] = [[1, 2], [0, 0], [3, 3]]
    with read.prepare_gts_matrix(bedfile=bed, par_gts_f=imputation, robust=True,
                                 ids=['child1', 'child2', 'child3']) as reader:
        variants = reader.match_variants()
        assert reader.n_imputed == 2
        assert reader.n_samples == 3
        genotypes = reader.read(variants)
    np.testing.assert_array_equal(variants.sid, ['a', 'b'])
    np.testing.assert_array_equal(variants.flipped, [False, True])
    np.testing.assert_array_equal(genotypes.gts.data, [
        [[0, 2], [1, 0], [1, 2]],
        [[2, 1], [1.25, 1.50], [1, 2]],
        [[1, 0], [.25, 1.25], [.25, 1.25]],
    ])
    np.testing.assert_array_equal(genotypes.num_obs_par_al, [[4, 4], [3, 3], [1, 2]])


def test_native_bgen_dosage_matches_probability_calculation(tmp_path):
    pytest.importorskip('bgen')
    from pathlib import Path
    from shutil import copyfile
    from snipar.bgen_io import BgenFile

    source = Path(__file__).resolve().parents[1] / 'snipar/example/example_data/chr_1.bgen'
    path = tmp_path / 'native.bgen'
    copyfile(source, path)
    with BgenFile(path) as reader:
        inds = [7, -1, 7], [2, 0, 0, 3]
        probabilities = reader.read(inds, np.float32)
        dosage = reader.read_dosage(*inds, order='C')
    assert dosage.flags.c_contiguous
    np.testing.assert_array_equal(dosage, probabilities[..., 0] + probabilities[..., 2])
