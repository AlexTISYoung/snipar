from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys

import h5py
import numpy as np
import pytest

from snipar import example, gwas
from snipar.gtarray import gtarray
from snipar.read.variants import bgen_variant_ids, match_variants


DATA = Path(example.__file__).parent / 'example_data'


@pytest.mark.parametrize('bgen_metadata', [True, False])
def test_match_unique_variants_and_alleles(tmp_path, bgen_metadata):
    sid = np.array(['a', 'copy', 'b', 'c', 'copy', 'wrong', 'd'])
    pos = np.arange(sid.size)
    chrom = np.full(sid.size, '1')
    alleles = np.tile(['A', 'G'], (sid.size, 1))
    ids = ['b', 'd', 'copy', 'missing', 'a', 'wrong', 'd']
    pairs = ['G,A', 'A,G', 'A,G', 'A,G', 'A,G', 'G,C', 'A,G']
    with h5py.File(tmp_path / 'imputed.hdf5', 'w') as f:
        if bgen_metadata:
            f['bim_columns'] = np.array(['id', 'rsid', 'allele_ids'], dtype='S')
            f['bim_values'] = np.column_stack((['.'] * len(ids), ids, pairs)).astype('S')
        else:
            f['bim_columns'] = np.array(['rsid', 'id', 'allele1', 'allele2'], dtype='S')
            f['bim_values'] = np.column_stack((['.'] * len(ids), ids, [p.split(',') for p in pairs])).astype('S')
        matched = match_variants(sid, chrom, pos, alleles, f)
        np.testing.assert_array_equal(matched.sid, ['b', 'a'])
        np.testing.assert_array_equal(matched.observed, [2, 0])
        np.testing.assert_array_equal(matched.imputed, [0, 4])
        np.testing.assert_array_equal(matched.flipped, [True, False])
        with pytest.raises(ValueError, match='No alleles match'):
            match_variants(sid, chrom, pos, alleles, f, snp_ids=['wrong'])
        with pytest.raises(ValueError, match='No SNPs'):
            match_variants(sid, chrom, pos, alleles, f, snp_ids=['copy', 'd'])


def test_bgen_aliases_require_unambiguous_matches(tmp_path):
    # A primary ID wins over another variant's alias. Duplicates and missing
    # identifiers may not match, even when their alleles agree.
    sid = np.array(['v0', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8', '.'])
    aliases = np.array(['rs0', 'rs1', 'v0', 'repeat', 'repeat', 'rs5', 'rs6', '.', '', '.'])
    imputed_ids = ['rs1', 'v0', 'rs0', 'repeat', 'rs5', 'rs5', 'rs6', '.', '']
    pos = np.arange(sid.size)
    with h5py.File(tmp_path / 'imputed.hdf5', 'w') as f:
        f['bim_columns'] = np.array(['id', 'allele_ids'], dtype='S')
        f['bim_values'] = np.column_stack((imputed_ids, ['A,G'] * len(imputed_ids))).astype('S')
        matched = match_variants(sid, np.full(sid.size, '1'), pos, np.tile(['A', 'G'], (sid.size, 1)),
                                 f, aliases=aliases)
        np.testing.assert_array_equal(matched.sid, ['v1', 'v0', 'v6'])
        np.testing.assert_array_equal(matched.imputed, [0, 1, 6])
    # Files containing just one SNP can have a usable SNPID and no rsID.
    result, _ = bgen_variant_ids(SimpleNamespace(ids=['only'], rsids=['.']))
    np.testing.assert_array_equal(result, ['only'])


def test_filtered_parental_allele_counts_stay_aligned():
    counts = np.arange(12).reshape(3, 4)
    g = gtarray(np.zeros((3, 3, 4)), ids=np.array(['a', 'b', 'c']),
                sid=np.array(['s1', 's2', 's3', 's4']), num_obs_par_al=counts)
    g.filter(np.array([False, True, False, True]))
    g.filter_ids(np.array(['c', 'a']))
    np.testing.assert_array_equal(g.num_obs_par_al, [[9, 11], [1, 3]])


def test_all_filtered_batch_skips_robust_sibling_read(monkeypatch):
    # The robust estimator normally requests sibling genotypes in a second
    # pass. A monomorphic batch should produce no fitted SNPs before that pass.
    ids = np.array(['a', 'b'])
    sid = np.array(['monomorphic'])
    g = gtarray(np.zeros((2, 3, 1)), ids=ids, sid=sid, fams=np.array(['f', 'f']))

    def read_genotypes(**kwargs):
        assert not kwargs['sib_diff']
        return g

    monkeypatch.setattr(gwas, '_var_dict', {'y_': np.ones(2), 'varcomps_': np.ones(1)})
    monkeypatch.setattr(gwas.read, 'get_gts_matrix', read_genotypes)
    freqs, tested, *_ = gwas.process_batch(sid, None, 0, None, pheno_ids=ids, robust=True)
    assert freqs.size == tested.size == 0


def test_empty_chromosome_names_inputs(tmp_path):
    imputation = tmp_path / 'parents.hdf5'
    with h5py.File(imputation, 'w') as f:
        f['bim_columns'] = np.array(['id', 'allele_ids'], dtype='S')
        f['bim_values'] = np.array([['not_in_bed', 'A,G']], dtype='S')
    with pytest.raises(ValueError, match=r'No SNPs.*chromosome 1.*chr_1.bed.*parents.hdf5'):
        gwas.process_chromosome(1, None, [], None, None, [], str(tmp_path / 'out'),
                                bedfile=str(DATA / 'chr_1.bed'), par_gts_f=imputation, parsum=True)


@pytest.mark.parametrize('fmt', ['bed', 'bgen'])
def test_gwas_batches_use_sparse_variant_overlap(tmp_path, fmt):
    genotype = DATA / 'chr_1'
    if fmt == 'bgen':
        # Some backends create a metadata cache next to the genotypes.
        genotype = tmp_path / 'chr_1'
        # Give the matched variants distinct SNPIDs while retaining their rsIDs
        # in the packaged phased BGEN. The one-byte substitutions leave all
        # BGEN lengths, offsets and genotype probabilities intact.
        bgen = (DATA / 'chr_1.bgen').read_bytes()
        for number, snpid in [(700, b'A'), (40, b'B'), (15, b'C'), (905, b'D')]:
            rsid = f'1_rs{number}'.encode()
            identifiers = len(rsid).to_bytes(2, 'little') + rsid + b'\x00\x00'
            old = b'\x01\x001' + identifiers
            assert bgen.count(old) == 1
            bgen = bgen.replace(old, b'\x01\x00' + snpid + identifiers)
        genotype.with_suffix('.bgen').write_bytes(bgen)
    # A small sample of sibling pairs with imputed parents, taken from the
    # packaged example. Most of the 1,000 observed SNPs have no imputed data.
    pedigree = np.loadtxt(DATA / 'pedigree.txt', dtype=str)
    families = np.arange(40).astype(str)
    children = np.char.endswith(pedigree[:, 1], '_0') | np.char.endswith(pedigree[:, 1], '_1')
    rows = np.isin(pedigree[:, 0], families) & children
    flags = np.full((rows.sum(), 2), 'False')
    pedigree = np.vstack(([*pedigree[0], 'has_father', 'has_mother'],
                          np.column_stack((pedigree[rows], flags))))
    variants = [700, 40, 15, 905]
    bim = np.loadtxt(DATA / 'chr_1.bim', dtype=str)[variants]
    bim[1, 4:6] = ['G', 'C']  # incompatible alleles
    bim[2, 4:6] = bim[2, 4:6][::-1]  # opposite coding
    rng = np.random.default_rng(78)
    imputation = tmp_path / 'parents.hdf5'
    with h5py.File(imputation, 'w') as f:
        f['families'] = families.astype('S')
        f['pedigree'] = pedigree.astype('S')
        f['bim_columns'] = np.array(['Chr', 'id', 'morgans', 'coordinate', 'allele1', 'allele2'], dtype='S')
        f['bim_values'] = bim.astype('S')
        f['imputed_par_gts'] = rng.uniform(0, 2, size=(families.size, len(variants)))
    command = [sys.executable, '-m', 'snipar.scripts.gwas', str(DATA / 'phenotype.txt'),
               '--' + fmt, str(genotype), '--imp', str(imputation.with_suffix('')),
               '--out', str(tmp_path / 'sparse'), '--batch_size', '1', '--cpus', '2',
               '--vc_list', '0.5', '0.5', '--parsum']
    run = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert run.returncode == 0, run.stdout
    assert 'Using 3 batches' in run.stdout
    with h5py.File(tmp_path / 'sparse.sumstats.hdf5') as output:
        estimates = output['estimate'][:]
        fitted = np.flatnonzero(np.isfinite(estimates[:, 0]))
        expected = ['C', 'A', 'D'] if fmt == 'bgen' else ['1_rs15', '1_rs700', '1_rs905']
        np.testing.assert_array_equal(output['bim'].asstr()[fitted, 1], expected)
        assert np.all(np.isfinite(estimates[fitted]))
        assert np.all(output['estimate_ses'][:][fitted] > 0)
        assert np.isfinite(output['estimate_covariance'][:][fitted]).all()
