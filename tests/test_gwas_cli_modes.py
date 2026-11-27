"""Exercise GWAS analysis modes on the same trios, duos, sibships and singletons."""

import gzip
from itertools import combinations

import h5py
import numpy as np
import pandas as pd
import pytest

from snipar.scripts import gwas, impute
from snipar.slmm import _identified_information
from snipar.vcf import convert_vcf


@pytest.fixture(scope='module')
def family_inputs(tmp_path_factory):
    root = tmp_path_factory.mktemp('gwas_cli')
    rng = np.random.default_rng(4298)
    positions = np.arange(101, 165)
    frequencies = rng.uniform(.2, .5, len(positions))
    phases, pedigree, phenotypes, segments, relationships = {}, [], [], [], []
    for family in range(24):
        fid = f'f{family}'
        father, mother = f'{fid}_p', f'{fid}_m'
        parents = rng.random((2, 2, len(positions))) < frequencies
        if family % 4 in (0, 1):
            phases[father] = parents[0]
        if family % 4 in (0, 2):
            phases[mother] = parents[1]
        transmissions = rng.integers(0, 2, (3, 2))
        if family % 4 == 1:
            transmissions[:, 1] = 0
        if family % 4 == 2:
            transmissions[:, 0] = 0
        children = [f'{fid}_{child}' for child in range(3)]
        family_effect = rng.normal(scale=.4)
        for iid, transmission in zip(children, transmissions):
            calls = parents[np.arange(2), transmission]
            phases[iid] = calls
            pedigree.append((fid, iid, father, mother))
            phenotypes.append((fid, iid, calls[:, 1].sum() + family_effect + rng.normal(), rng.normal()))
        for i, j in combinations(range(3), 2):
            state = int(np.count_nonzero(transmissions[i] == transmissions[j]))
            segments.append((children[i], children[j], 22, 101, 164, state))
            relationships.append((children[i], children[j], .4))
    for person in range(12):
        iid = f'unrelated{person}'
        phases[iid] = rng.random((2, len(positions))) < frequencies
        pedigree.append((iid, iid, '0', '0'))
        phenotypes.append((iid, iid, phases[iid][:, 1].sum() + rng.normal(), rng.normal()))
    for parent in ('father', 'mother'):
        fid = f'duo_{parent}'
        father, mother = f'{fid}_p', f'{fid}_m'
        iid = f'{fid}_child'
        parents = rng.random((2, 2, len(positions))) < frequencies
        phases[father if parent == 'father' else mother] = parents[0 if parent == 'father' else 1]
        phases[iid] = parents[:, 0]
        pedigree.append((fid, iid, father, mother))
        phenotypes.append((fid, iid, phases[iid][:, 1].sum() + rng.normal(), rng.normal()))
    vcf = root / 'input.vcf.gz'
    with gzip.open(vcf, 'wt') as output:
        output.write('##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
        output.write('\t'.join(phases) + '\n')
        for index, position in enumerate(positions):
            calls = '\t'.join(f'{int(sample[0, index])}|{int(sample[1, index])}' for sample in phases.values())
            output.write(f'22\t{position}\ts{index}\tA\tG\t.\tPASS\t.\tGT\t{calls}\n')
    for name, header, rows in (
        ('pedigree.txt', 'FID IID FATHER_ID MOTHER_ID', pedigree),
        ('phenotype.txt', 'FID IID trait alternative', phenotypes),
        ('grm.seg', 'ID1 ID2 PropIBD', relationships),
    ):
        with (root / name).open('w') as output:
            output.write(header + '\n')
            output.writelines(' '.join(map(str, row)) + '\n' for row in rows)
    with gzip.open(root / 'truth.segments.gz', 'wt') as output:
        output.write('ID1 ID2 Chr start_coordinate stop_coordinate IBDType\n')
        output.writelines(' '.join(map(str, row)) + '\n' for row in segments)
    bed = convert_vcf(vcf, root / 'bed').bed
    impute.main(impute.parser.parse_args(['--vcf', str(vcf), '--pedigree', str(root / 'pedigree.txt'),
                                        '--ibd', str(root / 'truth'), '--out', str(root / 'parents'),
                                        '--threads', '1', '-silent_progress']))
    assert (root / 'parents.hdf5').is_file()
    return root, vcf, bed


@pytest.mark.parametrize('mode, with_imputation, options, columns', [
    ('young', True, [], {'paternal_Beta', 'maternal_Beta'}),
    ('meta', False, [], {'direct_Beta'}),
    ('sibling', False, ['--sib_diff'], {'direct_Beta'}),
    ('unified', True, ['--impute_unrel'], {'paternal_Beta', 'maternal_Beta'}),
    ('unified_observed', False, ['--impute_unrel'], {'paternal_Beta', 'maternal_Beta'}),
    ('robust', True, ['--robust'], {'direct_Beta'}),
    ('parental_sum', True, ['--parsum'], {'avg_NTC_Beta'}),
    ('sibling_effect', True, ['--fit_sib'], {'sib_Beta'}),
    ('grm_only', True, ['--no_sib_var', '--vc_list', '.4', '.6'], {'direct_Beta'}),
])
def test_gwas_modes_agree_between_vcf_and_bed(family_inputs, tmp_path, capsys,
                                             mode, with_imputation, options, columns):
    root, vcf, bed = family_inputs
    pedigree = ['--imp', str(root / 'parents')] if with_imputation else ['--pedigree', str(root / 'pedigree.txt')]
    if mode == 'grm_only':
        options = [*options, '--grm', str(root / 'grm.seg'), '--sparse_thresh', '.1', '--cpus', '2']
    else:
        options = [*options, '--vc_list', '.4', '.6']
    common = [str(root / 'phenotype.txt'), *pedigree, *options, '--threads', '1', '--batch_size', '20']
    outputs = []
    for genotype in ('bed', 'vcf'):
        out = tmp_path / genotype
        gwas.main(gwas.parser.parse_args([*common, f'--{genotype}', str(bed if genotype == 'bed' else vcf),
                                         '--out', str(out)]))
        log = capsys.readouterr().out
        if mode in ('unified', 'unified_observed'):
            expected = 86 if with_imputation else 30
            assert f'Using {expected} individuals eligible' in log
        stats = pd.read_csv(f'{out}.sumstats.gz', sep=r'\s+')
        assert columns <= set(stats)
        # The robust estimate requires four identifiable design-specific fits.
        assert len(stats) == 64 and stats.direct_Beta.notna().sum() >= (20 if mode == 'robust' else 50)
        outputs.append(out)
    with h5py.File(f'{outputs[0]}.sumstats.hdf5') as bed_result, h5py.File(f'{outputs[1]}.sumstats.hdf5') as vcf_result:
        for field in ('estimate', 'estimate_ses', 'estimate_covariance', 'freqs'):
            np.testing.assert_allclose(vcf_result[field][:], bed_result[field][:], equal_nan=True)


def test_gwas_phenotype_selection_and_output_options(family_inputs, tmp_path):
    root, vcf, _ = family_inputs
    phenotypes = pd.read_csv(root / 'phenotype.txt', sep=r'\s+', dtype={'IID': str})
    kept = phenotypes.iloc[:48].copy()
    covariates = kept[['FID', 'IID']].assign(age=np.arange(len(kept)) % 7)
    covariates.to_csv(tmp_path / 'covar', sep=' ', index=False)
    kept['IID'].to_csv(tmp_path / 'keep', index=False, header=False)
    kept.loc[0, 'alternative'] = np.nan
    kept.to_csv(tmp_path / 'phenotypes', sep=' ', index=False, na_rep='MISSING')
    common = [str(tmp_path / 'phenotypes'), '--vcf', str(vcf), '--imp', str(root / 'parents'),
              '--vc_list', '.4', '.6', '--covar', str(tmp_path / 'covar'), '--keep', str(tmp_path / 'keep'),
              '--min_maf', '.05', '--max_missing', '0', '--missing_char', 'MISSING', '--no_hdf5_out']
    results = []
    for selector in (('--phen', 'alternative'), ('--phen_index', '2')):
        prefix = tmp_path / selector[0][2:]
        gwas.main(gwas.parser.parse_args([*common, *selector, '--out', str(prefix), '--vc_out', str(prefix)]))
        np.testing.assert_array_equal(np.loadtxt(prefix), [.4, .6])
        assert not prefix.with_suffix('.sumstats.hdf5').exists()
        results.append(pd.read_csv(f'{prefix}.sumstats.gz', sep=r'\s+'))
    pd.testing.assert_frame_equal(*results)


def test_gls_skips_roundoff_scale_collinearity_within_a_batch():
    # A small positive eigenvalue caused by covariance projection roundoff
    # must not turn collinear genotype columns into an identifiable design.
    eps = np.finfo(float).eps
    stable = np.array([[3., .25], [.25, 4.]])
    nearly_collinear = np.array([[1., 1.], [1. - 8 * eps, 1.]])
    information = np.array([stable, nearly_collinear])
    valid, covariance = _identified_information(information, np.diagonal(information, axis1=1, axis2=2), 24)
    np.testing.assert_array_equal(valid, [0])
    np.testing.assert_allclose(covariance[0], np.linalg.inv(stable), rtol=1e-14)


def test_gwas_quality_filters_match_input_subset(family_inputs, tmp_path):
    root, vcf, _ = family_inputs
    annotated = tmp_path / 'annotated.vcf.gz'
    scores = tmp_path / 'quality.tsv'
    selected = set()
    with gzip.open(vcf, 'rt') as original, gzip.open(annotated, 'wt') as target, scores.open('w') as table:
        target.write(next(original))
        target.write('##INFO=<ID=R2,Number=1,Type=Float,Description="Imputation quality">\n')
        target.write(next(original))
        table.write('ID\tINFO\n')
        for i, line in enumerate(original):
            fields = line.split('\t', 8)
            keep = i % 3 != 1
            score = '.99' if keep else '.9899'
            if keep:
                selected.add(fields[2])
            table.write(f'{fields[2]}\t{score}\n')
            fields[7] = f'R2={score}'
            target.write('\t'.join(fields))

    bed = convert_vcf(vcf, tmp_path / 'bed', variant_ids=selected).bed
    bgen = convert_vcf(vcf, tmp_path / 'bgen', output='phased').bgen
    common = [str(root / 'phenotype.txt'), '--imp', str(root / 'parents'), '--vc_list', '.4', '.6',
              '--threads', '1', '--cpus', '2', '--batch_size', '8']
    for name, options in {
        'bed': ['--bed', bed],
        'vcf': ['--vcf', str(annotated)],
        'bgen': ['--bgen', bgen, '--bgen-quality', str(scores)],
    }.items():
        gwas.main(gwas.parser.parse_args([*common, *options, '--out', str(tmp_path / name)]))

    with h5py.File(tmp_path / 'bed.sumstats.hdf5') as expected:
        assert set(expected['bim'].asstr()[:, 1]) == selected
        for name in ('vcf', 'bgen'):
            with h5py.File(tmp_path / f'{name}.sumstats.hdf5') as actual:
                np.testing.assert_array_equal(actual['bim'][:], expected['bim'][:])
                for field in ('estimate', 'estimate_ses', 'estimate_covariance', 'freqs'):
                    np.testing.assert_allclose(actual[field][:], expected[field][:], rtol=1e-7)
