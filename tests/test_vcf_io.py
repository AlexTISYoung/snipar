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


def _vcf(records, samples=('S1', 'S2', 'S3', 'S4', 'S5'), version='4.2', info_headers=''):
    header = (f'##fileformat=VCFv{version}\n{info_headers}##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
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


def test_vcf_info_quality_filters_before_genotype_decoding(tmp_path):
    source = tmp_path / 'calls.vcf'
    kept = ('boundary', 'last_field')
    scores = [
        ('boundary', 'R2=.8;DR2=6e-1;AF=.3'),
        ('last_field', 'AF=.5;DR2=0.999;R2=0.9'),
        ('low', 'R2=.7;DR2=.9'),
        ('one_low', 'R2=.9;DR2=.2'),
        ('prefix_match', 'OTHER_R2=.9;DR2=.9'),
        ('wrong_case', 'r2=.9;DR2=.9'),
        ('missing', 'DR2=.9'),
        ('flag', 'R2;DR2=.9'),
        ('missing_value', 'R2=.;DR2=.9'),
        ('nan', 'R2=NaN;DR2=.9'),
        ('infinite', 'R2=inf;DR2=.9'),
        ('nonscalar', 'R2=.9,1;DR2=.9'),
        ('malformed', 'R2=0_.9;DR2=.9'),
        ('duplicate', 'R2=.9;DR2=.9;R2=.2'),
        ('maf', 'R2=.9;DR2=.9;AF=.05'),
    ]
    records = [('1', str(pos), sid, 'A', 'G', '.', 'PASS', info, 'GT',
                *(['0|1', '0|0'] if sid in kept else ['bad', 'bad']))
               for pos, (sid, info) in enumerate(scores, 1)]
    source.write_bytes(_vcf(records, samples=('S1', 'S2')))
    result = convert_vcf(source, tmp_path / 'filtered', min_quality={'R2': .8, 'DR2': .6}, min_info_maf=.2)
    assert result.n_records == len(records)
    assert result.n_quality_removed == len(records) - len(kept) - 1
    assert dict(result.quality_filters) == {'R2': .8, 'DR2': .6}
    with open_bed(result.bed + '.bed') as reader:
        np.testing.assert_array_equal(reader.sid, kept)
        np.testing.assert_array_equal(reader.read(), [[1, 1], [2, 2]])


def test_vcf_auto_quality_uses_declared_metric_and_can_be_overridden(tmp_path):
    source = tmp_path / 'calls.vcf'
    source.write_bytes(_vcf([
        ('1', '1', 'equivalent', 'A', 'G', '.', 'PASS', 'INFO=.99;DR2=.3', 'GT', '0|0', '0|1'),
        ('1', '2', 'low', 'A', 'G', '.', 'PASS', 'DR2=1;INFO=.98', 'GT', '0|1', '0|1'),
        ('1', '3', 'missing', 'A', 'G', '.', 'PASS', 'DR2=1', 'GT', '0|1', '0|0'),
    ], samples=('S1', 'S2'), info_headers=(
        '##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R2">\n'
        '##INFO=<ID=INFO,Description="A quality score, with comma",Number=1,Type=Float>\n'
    )))
    for name, options, expected, expected_filters in [
        ('default', {}, ['equivalent'], {'INFO': .99}),
        ('override', {'min_quality': {'DR2': .9}}, ['low', 'missing'], {'DR2': .9}),
        ('disabled', {'min_quality': {}}, ['equivalent', 'low', 'missing'], {}),
    ]:
        result = convert_vcf(source, tmp_path / name, **options)
        assert dict(result.quality_filters) == expected_filters
        with open_bed(result.bed + '.bed') as reader:
            np.testing.assert_array_equal(reader.sid, expected)


def test_vcf_auto_quality_is_selected_per_file_from_scalar_numeric_fields(tmp_path):
    (tmp_path / 'chr_1.vcf').write_bytes(_vcf([
        ('1', '1', 'low', 'A', 'G', '.', 'PASS', 'Rsq=.5', 'GT', 'bad', 'bad'),
        ('1', '2', 'high', 'A', 'G', '.', 'PASS', 'Rsq=1', 'GT', '0|1', '0|0'),
    ], samples=('S1', 'S2'), info_headers='##INFO=<ID=Rsq,Number=A,Type=Float,Description="Dosage R2">\n'))
    (tmp_path / 'chr_2.vcf').write_bytes(_vcf([
        ('2', '1', 'no_metric', 'A', 'G', '.', 'PASS', 'INFO=label', 'GT', '0|0', '0|1'),
    ], samples=('S1', 'S2'), info_headers=(
        '##INFO=<ID=INFO,Number=1,Type=String,Description="A label">\n'
        '##INFO=<ID=DR2,Number=R,Type=Float,Description="Values per allele">\n'
    )))
    result = convert_vcf(tmp_path / 'chr_@.vcf', tmp_path / 'output')
    assert result.quality_filters == (('Rsq', .99),)
    assert result.n_quality_removed == 1
    for chrom, expected in [(1, 'high'), (2, 'no_metric')]:
        with open_bed(result.bed.replace('@', str(chrom)) + '.bed') as reader:
            np.testing.assert_array_equal(reader.sid, [expected])


@pytest.mark.parametrize('filter_, message', [
    ({'INFO/DR2': .9}, 'Invalid VCF INFO quality tag'),
    ({'R2': np.inf}, 'must be a finite number'),
    ({'R2': np.nan}, 'must be a finite number'),
    ({'R2': 'NA'}, 'must be a finite number'),
    ([('R2', .9)], 'must map INFO tags'),
])
def test_vcf_rejects_invalid_info_quality_options_before_reading(tmp_path, filter_, message):
    with pytest.raises(ValueError, match=message):
        convert_vcf(tmp_path / 'absent.vcf', tmp_path / 'output', min_quality=filter_)


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


def test_vcf_header_selector_aligns_the_subset_and_checks_panel_membership(tmp_path):
    for chrom, samples, calls in [(1, ('001', 'NA', 'extra'), ('0/0', '1/1', '0/0')),
                                 (2, ('extra', 'NA', '001'), ('1/1', '1/1', '0/0'))]:
        (tmp_path / f'chr{chrom}.vcf').write_bytes(_vcf([
            (str(chrom), '5', f's{chrom}', 'A', 'G', '.', 'PASS', '.', 'GT', *calls),
        ], samples))
    headers = []

    def select(samples):
        headers.append(samples)
        return {'NA'}

    source = tmp_path / 'chr@.vcf'
    result = convert_vcf(source, tmp_path / 'output', select_samples=select)
    assert headers == [('001', 'NA', 'extra')]
    assert result.n_samples == 1
    for chrom in (1, 2):
        with open_bed(result.bed.replace('@', str(chrom)) + '.bed') as reader:
            np.testing.assert_array_equal(reader.iid, ['NA'])
            np.testing.assert_array_equal(reader.read()[:, 0], [0])

    chromosome_two = tmp_path / 'chr2.vcf'
    chromosome_two.write_bytes(chromosome_two.read_bytes().replace(b'extra', b'other'))
    with pytest.raises(ValueError, match='same set of samples'):
        convert_vcf(source, tmp_path / 'mismatch', select_samples=select)


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
