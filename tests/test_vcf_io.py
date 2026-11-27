import gzip
from pathlib import Path
import struct
import zlib

from bed_reader import open_bed
import numpy as np
import pytest

from snipar.example.snipar_example_data import load_example_data
from snipar.utilities import open_bgen
from snipar.vcf import convert_vcf


def _vcf(records, samples=('S1', 'S2', 'S3', 'S4', 'S5'), version='4.2'):
    header = (f'##fileformat=VCFv{version}\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
    return (header + '\t'.join(samples) + '\n' + '\n'.join('\t'.join(row) for row in records) + '\n').encode()


def _bgzf(payload):
    compressor = zlib.compressobj(wbits=-15)
    data = compressor.compress(payload) + compressor.flush()
    prefix = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00'
    return prefix + struct.pack('<H', 25 + len(data)) + data + struct.pack('<II', zlib.crc32(payload), len(payload))


def test_bundled_example_vcf_matches_the_bgen_haplotypes(tmp_path):
    data = tmp_path / 'example_data'
    load_example_data(data)
    converted = convert_vcf(data / 'chr_1.vcf.gz', tmp_path / 'phased', output='phased')
    assert converted.chromosomes == (1,)
    with open_bgen(data / 'chr_1.bgen') as original, open_bgen(converted.bgen + '.bgen') as phased:
        np.testing.assert_array_equal(phased.samples, original.samples)
        np.testing.assert_array_equal(phased.rsids, original.rsids)
        np.testing.assert_array_equal(phased.positions, original.positions)
        np.testing.assert_array_equal(phased.allele_ids, original.allele_ids)
        for start in range(0, original.nvariants, 64):
            tile = (slice(None), slice(start, start + 64))
            np.testing.assert_array_equal(phased.read(tile, dtype=np.float32),
                                          original.read(tile, dtype=np.float32))


@pytest.mark.parametrize('compression', ['none', 'gzip', 'bgzf'])
def test_vcf_dosages_and_phases_roundtrip(tmp_path, compression):
    source = tmp_path / 'calls.vcf'
    payload = _vcf([
        ('chr1', '12', 'var_12', 'A', 'G', '.', 'PASS', 'AF=0.5', 'GT', '0|0', '0|1', '1|0', '1|1', '.|.'),
        ('chr1', '24', '.', 'T', 'TC', '.', '.', '.', 'DP:GT', '3:|0|1', '4:1|0', '5:1|1', '9:0|0', '7:1|.'),
    ], version='4.5')
    if compression == 'gzip':
        source = source.with_suffix('.vcf.gz')
        payload = gzip.compress(payload)
    elif compression == 'bgzf':
        source = source.with_suffix('.vcf.bgz')
        payload = _bgzf(payload[:128]) + _bgzf(payload[128:]) + _bgzf(b'')
    source.write_bytes(payload)
    bed = convert_vcf(source, tmp_path / 'doses')
    bgen = convert_vcf(source, tmp_path / 'phases', output='auto')
    assert bed.n_samples == bgen.n_samples == 5
    assert bed.n_records == bed.n_variants == bgen.n_variants == 2
    assert bgen.phase_fallback is None
    with open_bed(bed.bed + '.bed') as observed, open_bgen(bgen.bgen + '.bgen') as phased:
        np.testing.assert_array_equal(observed.iid, phased.samples)
        np.testing.assert_array_equal(observed.sid, ['var_12', '1:24:T:TC'])
        np.testing.assert_array_equal(observed.sid, phased.ids)
        np.testing.assert_array_equal(observed.bp_position, phased.positions)
        np.testing.assert_array_equal(phased.allele_ids, ['A,G', 'T,TC'])
        haplotypes = phased.read()[..., [0, 2]]
        np.testing.assert_array_equal(haplotypes[:, 0], [[1, 1], [1, 0], [0, 1], [0, 0], [np.nan, np.nan]])
        np.testing.assert_array_equal(haplotypes[:, 1], [[1, 0], [0, 1], [0, 0], [1, 1], [np.nan, np.nan]])
        np.testing.assert_array_equal(haplotypes.sum(axis=2), observed.read())


def test_vcf_filters_sites_before_decoding_and_maf_after_sample_selection(tmp_path):
    source = tmp_path / 'calls.vcf'
    source.write_bytes(_vcf([
        ('1', '1', 'rare', 'A', 'G', '.', 'PASS', 'AF=0.5;AC=1;AN=10', 'GT', 'bad', 'bad', 'bad', 'bad', 'bad'),
        ('1', '2', 'excluded', 'A', 'G', '.', 'q10', 'AF=0.2', 'GT', 'bad', 'bad', 'bad', 'bad', 'bad'),
        ('1', '3', 'exclude_id', 'A', 'G', '.', 'PASS', 'AF=0.5', 'GT', 'bad', 'bad', 'bad', 'bad', 'bad'),
        ('1', '4', 'multi', 'A', 'G,T', '.', 'PASS', 'AF=0.3,0.2', 'GT', 'bad', 'bad', 'bad', 'bad', 'bad'),
        ('1', '5', 'sv', 'A', '<DEL>', '.', 'PASS', '.', 'GT', 'bad', 'bad', 'bad', 'bad', 'bad'),
        ('1', '6', 'subset_common', 'A', 'G', '.', 'PASS', 'AF=bad', 'GT', '0/1', '0/0', '0/0', '0/0', '0/0'),
        ('1', '7', 'subset_fixed', 'A', 'G', '.', 'PASS', 'AF=0.3', 'GT', '1/1', '0/0', '0/0', '0/0', '1/1'),
        ('1', '8', 'ref_rare', 'A', 'G', '.', 'PASS', 'AF=0.999', 'GT', 'bad', 'bad', 'bad', 'bad', 'bad'),
        ('1', '9', 'no_gt', 'A', 'G', '.', 'PASS', '.', 'DP', '0', '0', '0', '0', '0'),
    ]))
    result = convert_vcf(source, tmp_path / 'filtered', min_info_maf=.2, min_maf=.2, pass_only=True,
                         sample_ids=['S1', 'S5'], variant_ids=['rare', 'excluded', 'multi', 'sv', 'subset_common',
                                                              'subset_fixed', 'ref_rare', 'no_gt'])
    assert (result.n_records, result.n_variants, result.n_samples) == (9, 1, 2)
    with open_bed(result.bed + '.bed') as calls:
        np.testing.assert_array_equal(calls.iid, ['S1', 'S5'])
        np.testing.assert_array_equal(calls.sid, ['subset_common'])
        np.testing.assert_array_equal(calls.read()[:, 0], [1, 2])


@pytest.mark.parametrize('calls,reason', [
    (('0/1:4', '1|0:4'), 'unphased heterozygous calls'),
    (('0|1:4', '1|0:8'), 'heterozygous calls in multiple phase sets'),
])
def test_vcf_cannot_use_local_phase_as_chromosome_wide_phase(tmp_path, calls, reason):
    source = tmp_path / 'mixed.vcf.gz'
    source.write_bytes(gzip.compress(_vcf([
        ('1', '5', 's1', 'A', 'G', '.', 'PASS', '.', 'GT:PS', calls[0], '0|0:.', '1|0:7'),
        ('1', '15', 's2', 'A', 'G', '.', 'PASS', '.', 'GT:PS', calls[1], '0|0:.', '0|1:7'),
    ], samples=('S1', 'S2', 'S3'))))
    result = convert_vcf(source, tmp_path / 'auto', output='auto')
    assert result.bed and not result.bgen
    assert reason in result.phase_fallback
    assert not list(Path(result.bed).parent.glob('*.bgen'))
    with open_bed(result.bed + '.bed') as reader:
        np.testing.assert_array_equal(reader.read(), [[1, 1], [2, 2], [1, 1]])
    with pytest.raises(ValueError, match=reason):
        convert_vcf(source, tmp_path / 'phased', output='phased')
    assert not list((tmp_path / 'phased').iterdir())


def test_vcf_multichromosome_sources_align_sample_ids(tmp_path):
    for chrom, samples, calls in [(1, ('001', 'NA'), ('0/0', '1/1')), (2, ('NA', '001'), ('1/1', '0/0'))]:
        source = tmp_path / f'input_{chrom}.vcf'
        source.write_bytes(_vcf([(str(chrom), '5', f's{chrom}', 'A', 'G', '.', 'PASS', '.', 'GT', *calls)], samples))
    result = convert_vcf(tmp_path / 'input_@', tmp_path / 'output', chromosomes=['2', '1'])
    assert result.chromosomes == (1, 2)
    assert result.bed.endswith('chr_@')
    for chrom in result.chromosomes:
        with open_bed(result.bed.replace('@', str(chrom)) + '.bed') as reader:
            np.testing.assert_array_equal(reader.iid, ['NA', '001'])
            np.testing.assert_array_equal(reader.read()[:, 0], [0, 2])


def test_vcf_accepts_diploid_phasing_flags_and_downgrades_allele_specific_sets(tmp_path):
    source = tmp_path / 'latest.vcf'
    source.write_bytes(_vcf([
        ('1', '7', 's1', 'A', 'G', '.', 'PASS', '.', 'GT:PSL', '|1|0:a,b', '|0|1:a,b'),
    ], samples=('S1', 'S2'), version='4.5'))
    result = convert_vcf(source, tmp_path / 'bed', output='auto')
    assert 'FORMAT/PSL' in result.phase_fallback
    with open_bed(result.bed + '.bed') as reader:
        np.testing.assert_array_equal(reader.read()[:, 0], [1, 1])


def test_vcf_phased_output_is_readable_by_native_bgen_backend(tmp_path):
    pytest.importorskip('bgen')
    from snipar.bgen_io import BgenFile

    source = tmp_path / 'calls.vcf'
    source.write_bytes(_vcf([
        ('22', '7', 's1', 'G', 'A', '.', 'PASS', '.', 'GT', '0|1', '1|0', '.|.'),
    ], samples=('sample_0', '010', 'NA')))
    result = convert_vcf(source, tmp_path / 'bgen', output='phased')
    with BgenFile(result.bgen + '.bgen') as reader:
        np.testing.assert_array_equal(reader.samples, ['sample_0', '010', 'NA'])
        np.testing.assert_array_equal(reader.rsids, ['s1'])
        np.testing.assert_array_equal(reader.read()[:, 0, :], [[1, 0, 0, 1], [0, 1, 1, 0], [np.nan] * 4])


def test_vcf_requires_diploid_calls_and_leaves_destinations_intact(tmp_path):
    source = tmp_path / 'bad.vcf'
    source.write_bytes(_vcf([
        ('1', '7', 's1', 'A', 'G', '.', 'PASS', '.', 'GT', '0', '0/0', '1/1', '0/1', '0/1'),
    ]))
    target = tmp_path / 'output'
    target.mkdir()
    sentinel = target / 'genotypes.bed'
    sentinel.write_bytes(b'original')
    with pytest.raises(ValueError, match=r':1:7: expected diploid'):
        convert_vcf(source, target)
    assert list(target.iterdir()) == [sentinel]
    assert sentinel.read_bytes() == b'original'


def test_vcf_rejects_empty_filtered_input(tmp_path):
    source = tmp_path / 'calls.vcf'
    source.write_bytes(_vcf([
        ('X', '7', 's1', 'A', 'G', '.', 'PASS', '.', 'GT', '0', '0/0', '1/1', '0/1', '0/1'),
    ]))
    with pytest.raises(ValueError, match='No biallelic diploid VCF variants'):
        convert_vcf(source, tmp_path / 'output')
