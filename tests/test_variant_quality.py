import gzip

import pytest

from snipar.variant_quality import DEFAULT_MIN_QUALITY, read_bgen_quality


def test_bgen_quality_default_is_inclusive_and_uses_info_before_r2(tmp_path):
    source = tmp_path / 'quality.txt'
    source.write_text('## per-variant scores\n#ID r2 Info\n'
                      '001 0.5 0.99\n'
                      'rsBelow 1 0.989999999\n'
                      'rsMissing 1 .\n'
                      'rsNonfinite 1 inf\n'
                      'rsInvalid 1 bad\n'
                      'NA 1 1\n')
    result = read_bgen_quality(source)
    assert result.filters == {'Info': DEFAULT_MIN_QUALITY}
    assert result.passing_ids == {'001', 'NA'}
    assert result.known_ids == {'001', 'rsBelow', 'rsMissing', 'rsNonfinite', 'rsInvalid', 'NA'}


def test_bgen_quality_named_rules_and_identifier_column(tmp_path):
    source = tmp_path / 'quality.tsv.gz'
    with gzip.open(source, 'wt') as stream:
        stream.write('Marker\tR2\tScore\tPos\n'
                     'rsPass\t0.85\t10\t20\n'
                     'rsLow\t0.85\t9.9\t30\n'
                     'rsEmpty\t\t11\t40\n'
                     'rsNan\tNaN\t11\t50\n')
    filters = {'R2': .85, 'Score': 10.}
    result = read_bgen_quality(source, filters, id_column='Marker')
    assert result.filters == filters
    assert result.passing_ids == {'rsPass'}
    assert result.known_ids == {'rsPass', 'rsLow', 'rsEmpty', 'rsNan'}


def test_bgen_biobank_mfi_uses_unique_primary_ids(tmp_path):
    source = tmp_path / 'ukb_mfi_chr1_v3.txt'
    source.write_text('id_01 rs1 101 C G 0.1 C 0.99\n'
                      'id_02 rs1 101 C A 0.01 A 0.98\n'
                      'id_03 rs3 302 C T 0.03 T 1.00\n')
    result = read_bgen_quality(source)
    assert result.filters == {'INFO': .99}
    assert result.passing_ids == {'id_01', 'id_03'}
    assert result.known_ids == {'id_01', 'id_02', 'id_03'}

    with pytest.raises(ValueError, match=r"Duplicate variant ID 'rs1'.*line 2"):
        read_bgen_quality(source, id_column='RSID')


@pytest.mark.parametrize(('contents', 'rules', 'id_column', 'message'), [
    ('', None, None, 'Empty BGEN quality table'),
    ('ID INFO\n', None, None, 'no variant IDs'),
    ('ID INFO INFO\nrs1 1 1\n', None, None, 'Duplicate column names'),
    ('name INFO\nrs1 1\n', None, None, 'No recognized variant ID column'),
    ('ID score\nrs1 1\n', None, None, 'No recognized imputation quality column'),
    ('ID INFO\nrs1 1\n', {'info': .9}, None, "no quality column 'info'"),
    ('ID INFO\nrs1 1\n', {'INFO': float('inf')}, None, 'must be finite'),
    ('ID INFO\nrs1 1\n', None, 'id', "no ID column 'id'"),
    ('ID INFO\nrs1\n', None, None, 'Expected 2 columns.*line 2'),
    ('ID INFO\nrs1 0.1\nrs1 1.0\n', None, None, "Duplicate variant ID 'rs1'.*line 3"),
])
def test_bgen_quality_rejects_ambiguous_or_malformed_tables(tmp_path, contents, rules, id_column, message):
    source = tmp_path / 'quality.txt'
    source.write_text(contents)
    with pytest.raises(ValueError, match=message):
        read_bgen_quality(source, rules, id_column=id_column)
