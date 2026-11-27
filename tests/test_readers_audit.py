from pathlib import Path

import h5py
import numpy as np
import pytest
from pysnptools.snpreader import Bed, SnpData

from snipar.errors import mendelian_errors
from snipar.read import build_ped_from_par_gts, get_gts_matrix, get_ids_with_par
from snipar.read.bed import read_sibs_from_bed
from snipar.read.bgen import _first_allele_dosage
from snipar.read.phenotype import read_phenotype
from snipar.read.variants import observed_variant_indices
from snipar.scripts import make_rdr_grms
from snipar.utilities import coord2linear, encode_str_array, linear2coord, parse_filelist


def test_phenotypes_preserve_ids_and_accept_mixed_whitespace(tmp_path):
    path = tmp_path / 'phenotypes.txt'
    path.write_text('FID   IID\tstatus other\nfam 0001 11 0\nfam NA 12 1\nfam 003 -9 2\nfam 004 nan 3\n')
    result = read_phenotype(path, column_index=1, na_values='-9')
    np.testing.assert_array_equal(result.ids, ['0001', 'NA'])
    np.testing.assert_array_equal(result.gts[:, 0], [11, 12])
    np.testing.assert_array_equal(read_phenotype(path, column='other').gts[:, 0], [0, 1, 2, 3])
    with pytest.raises(ValueError, match='out of range'):
        read_phenotype(path, column_index=0)


def test_variant_indices_reject_duplicates_in_observed_data_even_with_requests():
    observed = np.array(['z', 'a', 'dup', 'b', 'dup', '.'])
    np.testing.assert_array_equal(observed_variant_indices(observed, None, False), [1, 3, 0])
    np.testing.assert_array_equal(observed_variant_indices(observed, ['b', 'dup', 'z'], True), [3, 0])


def test_bgen_first_allele_dosage_handles_both_probability_layouts():
    class Bgen:
        def __init__(self, probabilities):
            self.probabilities = probabilities

        def read(self, index, dtype):
            return np.array(self.probabilities, dtype=dtype)[None, :, :]

    np.testing.assert_allclose(_first_allele_dosage(Bgen([[.2, .3, .5]]), None, None), [[.7]])
    np.testing.assert_allclose(_first_allele_dosage(Bgen([[.2, .8, .5, .5]]), None, None), [[.7]])


def test_packed_coordinates_remain_exact_beyond_float_precision():
    row = 100_000_001
    for col in (0, 3, row):
        assert linear2coord(coord2linear(row, col)) == (row, col)
    with pytest.raises(ValueError):
        linear2coord(-1)


def test_paired_file_patterns_accept_one_shared_imputation(tmp_path):
    for filename in ('gt_2.bed', 'gt_4.bed', 'imputed.hdf5'):
        (tmp_path / filename).touch()
    obs, imp, chromosomes = parse_filelist(str(tmp_path / 'gt_@'), str(tmp_path / 'imputed'), 'bed')
    np.testing.assert_array_equal(chromosomes, [2, 4])
    assert [Path(f).name for f in obs] == ['gt_2.bed', 'gt_4.bed']
    assert len(imp) == 2 and imp[0] == imp[1]


@pytest.fixture
def genotypes(tmp_path):
    bedfile = tmp_path / 'observed.bed'
    ids = np.array(['0001', '0002', 'father', 'mother'])
    ped = np.array([['very_long_family_identifier', iid, 'father', 'mother'] for iid in ids[:2]])
    snps = np.array(['rs1', 'rs2', 'rs3'])
    positions = np.array([[1, 0, 1], [1, 0, 2], [1, 0, 3]])
    values = np.array([[0., 1., 1.], [2., 0., np.nan], [1., 2., 1.], [2., 1., np.nan]])
    Bed.write(str(bedfile), SnpData(iid=np.column_stack((ids, ids)), sid=snps,
                                    pos=positions, val=values), count_A1=True)
    return bedfile, ped, values


def test_sibling_reader_counts_each_individual_once(genotypes):
    bed, _, values = genotypes
    pairs = [['0001', '0002'], ['0002', 'father'], ['father', 'absent']]
    result = read_sibs_from_bed(str(bed), np.array(pairs))
    np.testing.assert_array_equal(result.ids, ['0001', '0002', 'father'])
    np.testing.assert_array_equal(result.gts.filled(np.nan), values[:3])
    with pytest.raises(ValueError, match='No sibling pairs'):
        read_sibs_from_bed(str(bed), np.array([['0001', 'unknown']]))


def test_family_reader_preserves_order_controls_and_long_family_ids(tmp_path, genotypes):
    bed, ped, values = genotypes
    path = tmp_path / 'imputed.hdf5'
    family = ped[0, 0]
    families = [family, f'_p_{family}', f'_m_{family}', f'_o_{family}']
    rows = [['FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'has_father', 'has_mother']]
    for fid, father, mother in zip(families, ['father', 'father', '', ''], ['mother', '', 'mother', '']):
        rows.extend([fid, iid, father, mother, str(bool(father)), str(bool(mother))] for iid in ped[:, 1])
    alleles = np.loadtxt(bed.with_suffix('.bim'), dtype=str, usecols=(1, 4, 5))
    with h5py.File(path, 'w') as file:
        file['pedigree'] = encode_str_array(rows)
        file['families'] = encode_str_array(families)
        file['bim_columns'] = encode_str_array(['rsid', 'allele1', 'allele2'])
        file['bim_values'] = encode_str_array(alleles)
        file['imputed_par_gts'] = np.array([[0., 0., 0.], [1., 1., 1.], [0., 1., 2.], [.5, .5, .5]])
    observed_ped, _ = build_ped_from_par_gts(path)
    assert observed_ped.shape == (2, 6)
    ids, labels = get_ids_with_par(bed, ped, None, np.array(['0002', '0001']))
    np.testing.assert_array_equal(labels, [family, family])
    all_genotypes = get_gts_matrix(bedfile=bed, par_gts_f=path, ids=ids, parsum=True,
                                   compute_controls=True)
    for i, matrix in enumerate(all_genotypes):
        np.testing.assert_array_equal(matrix.ids, ['0002', '0001'])
        np.testing.assert_array_equal(matrix.fams, [families[i]] * 2)
        np.testing.assert_array_equal(matrix.gts[:, 0, :].filled(np.nan), values[:2][::-1])
    np.testing.assert_array_equal(all_genotypes[3].gts[:, 1, :], np.ones((2, 3)))
    with pytest.raises(ValueError, match='require an imputation file'):
        get_gts_matrix(bedfile=bed, ped=ped, compute_controls=True)


def test_mendelian_errors_exclude_missing_pairs(genotypes):
    bed, ped, _ = genotypes
    result = mendelian_errors(ped, bedfile=str(bed))
    np.testing.assert_array_equal(result.ME, [1, 1, 0])
    np.testing.assert_allclose(result.sum_het, [4 * .625 * .375, 1, .25])


def test_rdr_grm_hdf5_and_gcta_agree_with_matrix_reference(tmp_path, monkeypatch, genotypes):
    bed, ped, _ = genotypes
    monkeypatch.chdir(tmp_path)
    # Use a second trio with two different observed parents so both generations
    # have nonconstant dosages in the same two variants.
    ped[1, 2:] = ['0001', 'mother']
    pedigree = tmp_path / 'pedigree.txt'
    np.savetxt(pedigree, ped, fmt='%s', header='FID IID FATHER_ID MOTHER_ID', comments='')
    args = [str(bed), 'out', '--ped', str(pedigree)]
    make_rdr_grms.main([*args, '--hdf5_out'])
    make_rdr_grms.main(args)
    with h5py.File(tmp_path / 'out_grms.hdf5') as source:
        assert source['N'][()] == 2
        np.testing.assert_array_equal(source['ped'].asstr()[:], ped[:, :2])
        expected = np.array([[1., -1.], [-1., 1.]])
        for key, reference in [('R', expected), ('R_par', expected), ('R_o_par', np.zeros((2, 2)))]:
            np.testing.assert_allclose(source[key][:], reference, atol=1e-7)
            np.testing.assert_allclose(np.fromfile(f'out_{key}.grm.bin', dtype='<f4'),
                                       reference[np.tril_indices(2)], atol=1e-7)
            np.testing.assert_array_equal(np.fromfile(f'out_{key}.grm.N.bin', dtype='<f4'), [2, 2, 2])


def test_tutorial_bias_uses_standard_ols_with_intercept():
    from scipy.stats import linregress
    from snipar.example.example_data.estimate_sim_effects import slope_and_se

    true = np.array([0., 1., 2., 3., 4., 9.])
    fitted = np.array([1.1, 3., 5.5, 7., 8.9, np.nan])
    oracle = linregress(true[:-1], fitted[:-1])
    np.testing.assert_allclose(slope_and_se(true, fitted), [oracle.slope, oracle.stderr])


def test_example_phenotype_simulation_respects_seed_and_genetic_variance(tmp_path, genotypes):
    from snipar.example.simulate_phenotype import main

    bed, _, values = genotypes
    main(['--bed', str(bed), '--outdir', str(tmp_path), '--ncausal', '3', '--h2', '1', '--seed', '17'])
    output = tmp_path / 'h2_1.0.fam'
    effects = np.loadtxt(tmp_path / 'h2_1.0_causal.txt', usecols=6)
    expected = np.where(np.isnan(values), np.nanmean(values, axis=0), values) @ effects
    np.testing.assert_allclose(np.loadtxt(output, usecols=2), expected)
    np.testing.assert_allclose(expected.var(), 1.)


def test_native_bgen_validates_indices_and_preserves_duplicates(tmp_path):
    import shutil
    from snipar import example

    pytest.importorskip('bgen')
    from snipar.bgen_io import BgenFile

    path = tmp_path / 'chr_1.bgen'
    shutil.copyfile(Path(example.__file__).parent / 'example_data' / path.name, path)
    with BgenFile(path) as reader:
        ref = reader.read(([0, -1], [-1, -1, 0]))
        mask = np.isin(np.arange(len(reader.ids)), [0, len(reader.ids) - 1])
        np.testing.assert_array_equal(ref[:, [2, 0]], reader.read(([0, -1], mask)))
        np.testing.assert_array_equal(ref[:, 0], ref[:, 1])
        assert reader.read(([], [])).shape == (0, 0, 4)
        for invalid in ([True], [[0]], [0.3], [len(reader.samples)]):
            with pytest.raises(IndexError):
                reader.read((invalid, [0]))
