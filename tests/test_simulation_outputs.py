"""Check the simulator's exported multi-chromosome weights and unphased dosages."""

import csv
import gzip
import json

import h5py
import numpy as np
import pytest
from pysnptools.snpreader import Bed

import snipar.simulate as simulation
from snipar.genotype_cli import genotype_inputs
from snipar.scripts import simulate


def test_vcf_founder_selection_and_batched_phased_input(tmp_path, monkeypatch, capsys):
    source = tmp_path / 'founders.vcf'
    rng = np.random.default_rng(336)
    haplotypes = rng.integers(0, 2, (8, 5, 2))
    with source.open('w') as output:
        output.write('##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
        output.write(''.join(f'\ts{i}' for i in range(len(haplotypes))) + '\n')
        for j in range(haplotypes.shape[1]):
            calls = '\t'.join(f'{a}|{b}' for a, b in haplotypes[:, j, :])
            output.write(f'22\t{j + 1}\tv{j}\tA\tG\t.\tPASS\t.\tGT\t{calls}\n')
    args = simulate.parser.parse_args(['2', '.5', str(tmp_path / 'sim'), '--vcf', str(source), '--nfam', '2'])
    monkeypatch.setattr(simulation, 'decode_map_from_pos', lambda chrom, pos: np.asarray(pos, dtype=float))
    monkeypatch.setattr(simulation, 'plan_batches', lambda *args, **kwargs: (1, 2))
    with genotype_inputs(args, output='phased', vcf_samples=simulate._vcf_founders) as resolved:
        observed, *_ = simulation.haps_from_bgen(resolved.bgen, None, nfam=2)
    assert '4 samples as phased BGEN' in capsys.readouterr().out
    expected = (haplotypes[:4] == 0).reshape(2, 2, 5, 2).transpose(1, 0, 2, 3)
    np.testing.assert_array_equal(observed[0], expected)


def test_two_chromosome_selection_outputs(tmp_path, monkeypatch):
    rng = np.random.default_rng(1028)
    chroms = [1, 2]
    maps = [np.arange(6, dtype=float), np.arange(5, dtype=float)]
    positions = [np.arange(1, len(m) + 1) for m in maps]
    ids = [np.array(['variant_{}_{}'.format(c, i) for i in range(len(m))])
           for c, m in zip(chroms, maps)]
    haps = [rng.integers(0, 2, (32, 2, len(m), 2)).astype(bool) for m in maps]
    alleles = [np.tile(['A', 'G'], (len(m), 1)) for m in maps]
    monkeypatch.setattr(simulate, 'haps_from_bgen',
                        lambda *args: (haps, maps, ids, alleles, positions, chroms))

    args = simulate.parser.parse_args([
        '8', '.5', str(tmp_path) + '/', '--bgen', 'chr_@', '--nfam', '16',
        '--n_stabilizing', '1', '--selection_variance', 'inf', '--unphased_impute',
    ])
    simulate.main(args)

    with (tmp_path / 'causal_effects.txt').open() as handle:
        effects = list(csv.DictReader(handle, delimiter=' '))
    assert [row['SNP'] for row in effects] == list(np.concatenate(ids))
    assert all(row['A1'] == 'A' and row['A2'] == 'G' for row in effects)
    exported_betas = np.array([float(row['direct']) for row in effects])

    dosages = []
    for chrom, snps in zip(chroms, ids):
        bed = Bed(str(tmp_path / ('chr_%d' % chrom)), count_A1=True).read()
        np.testing.assert_array_equal(bed.sid, snps)
        dosages.append(bed.val)
        with h5py.File(tmp_path / ('unphased_impute_chr_%d.hdf5' % chrom)) as h5:
            imputed = h5['imputed_par_gts'][:]
            assert imputed.shape == (16, len(snps))
            assert np.isfinite(imputed).all() and np.all((imputed >= 0) & (imputed <= 4))
            np.testing.assert_array_equal(h5['families'][:].astype(str), bed.iid[::2, 0])
            np.testing.assert_array_equal(h5['bim_values'][:, 1].astype(str), snps)

        expected = np.full(imputed.shape, np.nan)
        frequencies = np.mean(bed.val, axis=0) / 2
        families = {iid: i for i, iid in enumerate(bed.iid[::2, 1])}
        with gzip.open(tmp_path / ('chr_%d.segments.gz' % chrom), 'rt') as handle:
            for segment in csv.DictReader(handle, delimiter='\t'):
                i = families[segment['ID1']]
                for j in range(int(segment['start_coordinate']) - 1, int(segment['stop_coordinate'])):
                    expected[i, j] = simulate.impute_from_sibs(
                        bed.val[2*i, j], bed.val[2*i+1, j], int(segment['IBDType']), frequencies[j])
        np.testing.assert_allclose(imputed, expected)

    with (tmp_path / 'pedigree.txt').open() as handle:
        children = [row for row in csv.DictReader(handle, delimiter=' ')
                    if row['FID'].startswith('1_')]
    metadata = json.loads((tmp_path / 'selection_metadata.json').read_text())
    scores = np.column_stack(dosages) @ exported_betas - metadata['genetic_origin']
    np.testing.assert_allclose(scores, [float(row['DIRECT']) for row in children],
                               rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize('option', ['--maf', '--min_maf'])
def test_ordinary_simulation_validates_requested_allele_frequencies(tmp_path, option):
    args = simulate.parser.parse_args(['4', '.5', str(tmp_path) + '/', '--nfam', '8',
                                      option, '0.7'])
    with pytest.raises(ValueError, match='MAF must be between'):
        simulate.main(args)
