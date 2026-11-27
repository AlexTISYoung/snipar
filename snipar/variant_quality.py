"""Select variants using per-variant imputation quality scores."""

import gzip
import math
from itertools import chain
from pathlib import Path
from typing import NamedTuple


DEFAULT_MIN_QUALITY = .99
QUALITY_FIELDS = ('INFO', 'R2', 'RSQ', 'DR2', 'IMPINFO')
_ID_FIELDS = ('ID', 'SNPID', 'SNP_ID', 'VARIANT_ID', 'ALTERNATE_IDS', 'SNP', 'RSID', 'RS_ID')
_MFI_COLUMNS = ('ID', 'RSID', 'POS', 'A1', 'A2', 'MAF', 'MINOR_ALLELE', 'INFO')


class BgenQuality(NamedTuple):
    passing_ids: set[str]
    known_ids: set[str]
    filters: dict[str, float]


def _default_column(columns, candidates):
    for candidate in candidates:
        if candidate in columns:
            return candidate
        matches = [name for name in columns if name.upper() == candidate]
        if len(matches) > 1:
            raise ValueError(f'Ambiguous BGEN quality table column for {candidate!r}: {matches}')
        if matches:
            return matches[0]
    return None


def _table_rows(stream):
    for line_number, line in enumerate(stream, 1):
        first_text = line.lstrip()
        if first_text and not first_text.startswith(('##', '# ', '#\t')):
            yield line_number, line.rstrip('\r\n')


def _numeric_passes(text, minimum):
    try:
        value = float(text)
    except ValueError:
        return False
    return math.isfinite(value) and value >= minimum


def read_bgen_quality(path, min_quality=None, *, id_column=None):
    """Return ID sets for a named quality table or a headerless UK Biobank MFI file.

    The quality table may be plain text or gzip-compressed. Its variant IDs
    must be unique. Threshold names in ``min_quality`` are exact column names.
    By default the first recognized quality metric uses a minimum of 0.99.
    """
    path = Path(path)
    opener = gzip.open if path.suffix.lower() in ('.gz', '.bgz', '.bgzf') else open
    with opener(path, 'rt', encoding='utf-8-sig') as stream:
        rows = _table_rows(stream)
        first_row = next(rows, None)
        if first_row is None:
            raise ValueError(f'Empty BGEN quality table: {path}')
        first = first_row[1].split()
        if len(first) == len(_MFI_COLUMNS) and first[2].isdecimal():
            columns = _MFI_COLUMNS
            rows = chain((first_row,), rows)
        else:
            first[0] = first[0].removeprefix('#')
            columns = first

        index = {name: i for i, name in enumerate(columns)}
        if len(index) != len(columns):
            raise ValueError(f'Duplicate column names in BGEN quality table: {path}')
        if id_column is None:
            id_column = _default_column(columns, _ID_FIELDS)
            if id_column is None:
                raise ValueError(f'No recognized variant ID column in BGEN quality table: {path}. '
                                 'Set --bgen-quality-id to its column name.')
        if id_column not in index:
            raise ValueError(f'BGEN quality table {path} has no ID column {id_column!r}')
        id_index = index[id_column]
        if min_quality is None:
            score_column = _default_column(columns, QUALITY_FIELDS)
            if score_column is None:
                raise ValueError(f'No recognized imputation quality column in BGEN quality table: {path}. '
                                 'Set --min-quality to its column name.')
            min_quality = {score_column: DEFAULT_MIN_QUALITY}

        filters = {}
        for name, minimum in min_quality.items():
            if name not in index:
                raise ValueError(f'BGEN quality table {path} has no quality column {name!r}')
            minimum = float(minimum)
            if not math.isfinite(minimum):
                raise ValueError(f'Minimum quality for {name!r} must be finite')
            filters[name] = minimum
        checks = tuple((index[name], minimum) for name, minimum in filters.items())

        known, passing = set(), set()
        width = len(columns)
        for line_number, line in rows:
            fields = line.split('\t') if '\t' in line else line.split()
            if len(fields) != width:
                fields = line.split()
            if len(fields) != width:
                raise ValueError(f'Expected {width} columns in BGEN quality table {path} at line {line_number}')
            identifier = fields[id_index].strip()
            if identifier in ('', '.'):
                continue
            if identifier in known:
                raise ValueError(f'Duplicate variant ID {identifier!r} in BGEN quality table {path} at line {line_number}')
            known.add(identifier)
            if all(_numeric_passes(fields[i], minimum) for i, minimum in checks):
                passing.add(identifier)
    if not known:
        raise ValueError(f'BGEN quality table has no variant IDs: {path}')
    return BgenQuality(passing, known, filters)
