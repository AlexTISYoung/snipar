"""Stream VCF hard calls into bounded, random-access genotype files."""

from collections.abc import Mapping
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
import gzip
import math
from pathlib import Path
import re
import struct
from tempfile import TemporaryDirectory
import zlib

import numpy as np
from numba import njit

from snipar.variant_quality import DEFAULT_MIN_QUALITY as DEFAULT_VARIANT_QUALITY, QUALITY_FIELDS


_INFO_ATTRIBUTES = re.compile(rb'(?:^|,)\s*([A-Za-z0-9_.-]+)=\s*("(?:\\.|[^"\\])*"|[^,]+)')
_QUALITY_RANK = {field.encode(): index for index, field in enumerate(QUALITY_FIELDS)}


@dataclass(frozen=True)
class VCFConversion:
    bed: str | None
    bgen: str | None
    chromosomes: tuple[int, ...]
    n_samples: int
    n_variants: int
    n_records: int
    phase_fallback: str | None = None
    quality_filters: tuple[tuple[str, float], ...] = ()
    n_quality_removed: int = 0


def _vcf_files(source, chromosomes):
    source = str(source)
    suffixes = ('', '.vcf', '.vcf.gz', '.vcf.bgz', '.vcf.bgzf')
    if '@' not in source:
        for suffix in suffixes:
            candidate = Path(source + suffix)
            if candidate.is_file():
                yield None, candidate
                return
        raise FileNotFoundError(f'VCF not found: {source}')
    found = False
    for chrom in chromosomes:
        prefix = source.replace('@', str(chrom))
        for suffix in suffixes:
            candidate = Path(prefix + suffix)
            if candidate.is_file():
                found = True
                yield chrom, candidate
                break
    if not found:
        raise FileNotFoundError(f'No VCF files match {source}')


@contextmanager
def _vcf_lines(source):
    with source.open('rb', buffering=1024 * 1024) as file:
        if file.peek(2)[:2] == b'\x1f\x8b':
            with gzip.GzipFile(fileobj=file) as compressed:
                yield compressed
        else:
            yield file


def _sample_header(lines, source, *, detect_quality=False):
    quality_tag = None
    quality_rank = len(_QUALITY_RANK)
    for line in lines:
        if line.startswith(b'##'):
            if detect_quality and quality_rank and line.startswith(b'##INFO=<'):
                attributes = {key: value.strip() for key, value in
                              _INFO_ATTRIBUTES.findall(line[8:].rstrip().removesuffix(b'>'))}
                tag = attributes.get(b'ID', b'')
                rank = _QUALITY_RANK.get(tag.upper(), len(_QUALITY_RANK))
                # Number=A is also a scalar because we keep biallelic sites.
                if (rank < quality_rank and attributes.get(b'Number') in (b'1', b'A') and
                        attributes.get(b'Type') in (b'Float', b'Integer')):
                    quality_tag = tag.decode()
                    quality_rank = rank
            continue
        fields = line.rstrip(b'\r\n').split(b'\t')
        if fields[:9] != [b'#CHROM', b'POS', b'ID', b'REF', b'ALT', b'QUAL', b'FILTER', b'INFO', b'FORMAT']:
            break
        samples = tuple(fields[9:])
        if not samples or len(set(samples)) != len(samples) or any(
            not sample or any(char <= 32 for char in sample) for sample in samples
        ):
            raise ValueError(f'VCF sample IDs must be nonempty, unique and contain no whitespace: {source}')
        return samples, quality_tag
    raise ValueError(f'VCF has no valid genotype header: {source}')


def _variant_fields(line):
    """Read the site columns without allocating a copy of the sample fields."""
    fields = []
    start = 0
    for _ in range(9):
        end = line.find(b'\t', start)
        if end < 0:
            raise ValueError('VCF record has fewer columns than its genotype header')
        fields.append(line[start:end])
        start = end + 1
    return fields, start


def _info_maf(info):
    # Prefer allele counts because AF fields are often rounded. These values
    # describe the input INFO population, which may differ from a later subset.
    ac, an, freq = None, None, None
    for field in info.split(b';'):
        try:
            if field.startswith(b'AC='):
                ac = int(field[3:])
            elif field.startswith(b'AN='):
                an = int(field[3:])
            elif field.startswith(b'AF='):
                freq = float(field[3:])
        except ValueError:
            continue
    if ac is not None and an is not None and 0 <= ac <= an and an > 0:
        return min(ac, an - ac) / an
    if freq is not None and 0 <= freq <= 1:
        return min(freq, 1 - freq)
    return None


def _info_quality_searches(min_quality):
    """Prepare delimiter-aware searches once for the requested INFO tags."""
    if min_quality is None:
        return ()
    if not isinstance(min_quality, Mapping):
        raise ValueError('VCF quality filters must map INFO tags to minimum scores')
    searches = []
    for tag, minimum in min_quality.items():
        if not isinstance(tag, str) or not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', tag):
            raise ValueError(f'Invalid VCF INFO quality tag: {tag!r}')
        try:
            minimum = float(minimum)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError(f'VCF INFO minimum for {tag} must be a finite number') from error
        if not math.isfinite(minimum):
            raise ValueError(f'VCF INFO minimum for {tag} must be a finite number')
        prefix = tag.encode() + b'='
        searches.append((prefix, b';' + prefix, minimum))
    return tuple(searches)


def _passes_info_quality(info, searches):
    # A handful of requested scores should not require parsing or allocating
    # every annotation in the record. Check delimiters rather than substrings.
    for prefix, delimited, minimum in searches:
        start = len(prefix)
        if not info.startswith(prefix):
            position = info.find(delimited)
            if position < 0:
                return False
            start += position + 1
        end = info.find(b';', start)
        if end < 0:
            value = info[start:]
        else:
            if info.find(delimited, end) >= 0:
                return False
            value = info[start:end]
        if b'_' in value or value != value.strip():
            return False
        try:
            score = float(value)
        except ValueError:
            return False
        if not math.isfinite(score) or score < minimum:
            return False
    return True


@njit(cache=True, inline='always')
def _decode_gt(line, start, end, a, b, phased, sample):
    a[sample] = -1
    b[sample] = -1
    phased[sample] = False
    # VCF 4.5 allows an explicit phasing separator before the first allele.
    first_phased = True
    if end > start and (line[start] == 47 or line[start] == 124):
        first_phased = line[start] == 124
        start += 1
    if end - start == 1 and line[start] == 46:  # .
        return True
    if end - start != 3:
        return False
    first, separator, second = line[start], line[start + 1], line[start + 2]
    if separator != 47 and separator != 124:  # / and |
        return False
    if (first != 46 and first != 48 and first != 49) or (second != 46 and second != 48 and second != 49):
        return False
    if first != 46 and second != 46:
        a[sample] = first - 48
        b[sample] = second - 48
        phased[sample] = first_phased and separator == 124
    return True


@njit(cache=True)
def _decode_samples(line, start, end, gt_field, ps_field, psl_field, source_indices, sample_order,
                    a, b, phased, phase_sets):
    # Most statistical-phasing outputs store GT by itself at a constant
    # stride. Decode only selected samples on this fast path.
    last_field = max(gt_field, ps_field, psl_field)
    if last_field == 0 and end - start == 4 * len(sample_order) - 1:
        fixed_width = True
        for pos in range(start + 3, end, 4):
            if line[pos] != 9:
                fixed_width = False
                break
        if fixed_width:
            for i, source in enumerate(source_indices):
                pos = start + 4 * source
                if not _decode_gt(line, pos, pos + 3, a, b, phased, i):
                    return 1
                phase_sets[i] = -1
            return 0

    source = 0
    while start < end:
        if source == len(sample_order):
            return 2
        sample = sample_order[source]
        stop = start
        while stop < end and line[stop] != 9:
            stop += 1
        if sample >= 0:
            part = 0
            field_start = start
            got_gt = False
            phase_sets[sample] = -1
            for pos in range(start, stop + 1):
                if pos != stop and line[pos] != 58:  # :
                    continue
                if part == gt_field:
                    if not _decode_gt(line, field_start, pos, a, b, phased, sample):
                        return 1
                    got_gt = True
                elif part == ps_field and phase_sets[sample] != -3 and pos > field_start and line[field_start] != 46:
                    value = 0
                    for digit in range(field_start, pos):
                        if line[digit] < 48 or line[digit] > 57 or value > 214748364:
                            return 3
                        value = 10 * value + line[digit] - 48
                    if value > 2147483647:
                        return 3
                    phase_sets[sample] = value
                elif part == psl_field:
                    for char in range(field_start, pos):
                        if line[char] != 46 and line[char] != 44:  # . and ,
                            phase_sets[sample] = -3
                            break
                if part >= last_field:
                    break
                field_start = pos + 1
                part += 1
            if not got_gt:
                a[sample] = -1
                b[sample] = -1
                phased[sample] = False
        source += 1
        start = stop + 1
    return 0 if source == len(sample_order) else 2


@njit(cache=True)
def _allele_counts(a, b):
    alt, called = 0, 0
    for i in range(len(a)):
        if a[i] >= 0:
            alt += a[i] + b[i]
            called += 2
    return alt, called


@njit(cache=True)
def _track_phasing(a, b, phased, phase_sets, previous_sets):
    has_phase = False
    for sample in range(len(a)):
        has_phase |= phased[sample]
        if a[sample] < 0 or a[sample] == b[sample]:
            continue
        if not phased[sample]:
            return False, 1
        if phase_sets[sample] == -3:
            return False, 3
        previous = previous_sets[sample]
        if previous == -2:
            previous_sets[sample] = phase_sets[sample]
        elif previous != phase_sets[sample]:
            return False, 2
    return has_phase, 0


@njit(cache=True)
def _pack_bed(a, b, genotypes):
    genotypes[:] = 0
    for i in range(len(a)):
        # BED stores 00=A1/A1, 01=missing, 10=A1/A2, and 11=A2/A2.
        count = a[i] + b[i]
        encoded = 1 if count < 0 else count + (count > 0)
        genotypes[i // 4] |= encoded << (2 * (i % 4))


@njit(cache=True)
def _pack_phased(a, b, ploidy, probabilities):
    probabilities[:] = 0
    for i in range(len(a)):
        ploidy[i] = 130 if a[i] < 0 else 2
        # One bit for the first-allele probability on each haplotype.
        # BGEN implicitly stores the complementary probabilities.
        encoded = (a[i] == 0) + 2 * (b[i] == 0)
        probabilities[i // 4] |= encoded << (2 * (i % 4))


class _BedWriter:
    def __init__(self, stack, prefix, samples):
        self.bed = stack.enter_context(prefix.with_suffix('.bed').open('xb'))
        self.bim = stack.enter_context(prefix.with_suffix('.bim').open('xb'))
        with prefix.with_suffix('.fam').open('xb') as fam:
            fam.writelines(sample + b' ' + sample + b' 0 0 0 -9\n' for sample in samples)
        self.bed.write(b'\x6c\x1b\x01')
        self.packed = np.empty((len(samples) + 3) // 4, dtype=np.uint8)

    def write(self, chrom, pos, sid, ref, alt, a, b):
        self.bim.write(b'\t'.join((chrom, sid, b'0', str(pos).encode(), ref, alt)) + b'\n')
        _pack_bed(a, b, self.packed)
        self.bed.write(self.packed)


def _string16(value):
    if len(value) > 65535:
        raise ValueError('VCF identifier is too long for BGEN')
    return struct.pack('<H', len(value)) + value


class _PhasedWriter:
    """Write BGEN layout 2 with zlib and lossless one-bit hard calls."""

    def __init__(self, stack, prefix, samples):
        self.file = stack.enter_context(prefix.with_suffix('.bgen').open('xb'))
        if samples[0] == b'sample_0':
            # The existing BGEN readers interpret this ID as a placeholder.
            with prefix.with_suffix('.sample').open('xb') as sample_file:
                sample_file.write(b'ID_1 ID_2 missing\n0 0 0\n')
                sample_file.writelines(sample + b' ' + sample + b' 0\n' for sample in samples)
        samples_size = 8 + sum(2 + len(sample) for sample in samples)
        self.file.write(struct.pack('<4I4s3I', 20 + samples_size, 20, 0, len(samples), b'bgen',
                                    (1 << 31) | (2 << 2) | 1, samples_size, len(samples)))
        self.file.writelines(_string16(sample) for sample in samples)
        self.raw = np.empty(10 + len(samples) + (len(samples) + 3) // 4, dtype=np.uint8)
        self.raw[:8] = np.frombuffer(struct.pack('<IHBB', len(samples), 2, 2, 2), dtype=np.uint8)
        self.ploidy = self.raw[8:8 + len(samples)]
        self.raw[8 + len(samples):10 + len(samples)] = 1  # phased, bits per probability
        self.packed = self.raw[10 + len(samples):]
        self.count = 0

    def write(self, chrom, pos, sid, ref, alt, a, b):
        _pack_phased(a, b, self.ploidy, self.packed)
        compressed = zlib.compress(self.raw, level=1)
        self.file.write(b''.join((_string16(sid), _string16(sid), _string16(chrom),
                                  struct.pack('<IHI', pos, 2, len(ref)), ref,
                                  struct.pack('<I', len(alt)), alt,
                                  struct.pack('<II', len(compressed) + 4, len(self.raw)), compressed)))
        self.count += 1

    def finish(self):
        self.file.seek(8)
        self.file.write(struct.pack('<I', self.count))


def convert_vcf(source, destination, *, chromosomes=None, output='bed', sample_ids=None,
                variant_ids=None, min_maf=0, min_info_maf=0, min_quality=None, pass_only=False,
                split_chromosomes=False):
    """Convert plain, gzip or BGZF VCF into per-autosome analysis files.

    ``bed`` retains diploid dosages. ``phased`` requires a globally consistent
    phase frame per sample and chromosome. ``auto`` uses phased BGEN when
    possible and otherwise preserves the calls as BED. Output variants use VCF
    IDs and REF-first allele pairs. Unidentified variants receive a
    ``CHROM:POS:REF:ALT`` ID.

    ``min_maf`` uses called alleles of the selected samples. The independent,
    opt-in ``min_info_maf`` filter uses the input's INFO/AC and INFO/AN or
    INFO/AF population. Sites with missing or invalid metadata stay eligible.
    ``min_quality`` maps case-sensitive INFO tags to inclusive minimum scores.
    Variants must have a finite scalar score for each requested tag and meet all
    thresholds. By default, apply a 0.99 cutoff when a recognized numeric
    quality tag is declared in the file header. An empty mapping disables it.
    Site and ID filters run before even scanning the sample fields.
    """
    if output not in {'bed', 'auto', 'phased'}:
        raise ValueError('VCF output must be bed, auto or phased')
    if not (0 <= min_maf <= .5 and 0 <= min_info_maf <= .5):
        raise ValueError('VCF minimum allele frequencies must be in [0, 0.5]')
    requested_quality = _info_quality_searches(min_quality)
    requested_chroms = tuple(dict.fromkeys(range(1, 23) if chromosomes is None else map(int, chromosomes)))
    if not requested_chroms or any(chrom not in range(1, 23) for chrom in requested_chroms):
        raise ValueError('VCF input supports autosomes 1 through 22')
    allowed_chroms = {str(chrom).encode(): chrom for chrom in requested_chroms}
    requested_samples = None if sample_ids is None else {str(sample).encode() for sample in sample_ids}
    requested_variants = None if variant_ids is None else {str(sid).encode() for sid in variant_ids}
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    samples = None
    bed, bgen, phase_sets_by_chrom = {}, {}, {}
    phase_possible, saw_phase = output != 'bed', False
    phase_fallback = None
    applied_quality = {}
    n_records = n_variants = n_quality_removed = 0

    with TemporaryDirectory(dir=destination) as work, ExitStack() as stack:
        work = Path(work)
        phased_handles = stack.enter_context(ExitStack())
        for expected_chrom, source_file in _vcf_files(source, requested_chroms):
            with _vcf_lines(source_file) as lines:
                header, auto_quality_tag = _sample_header(lines, source_file, detect_quality=min_quality is None)
                quality_searches = (requested_quality if auto_quality_tag is None else
                                    _info_quality_searches({auto_quality_tag: DEFAULT_VARIANT_QUALITY}))
                for prefix, _, minimum in quality_searches:
                    applied_quality[prefix[:-1].decode()] = minimum
                selected = (header if requested_samples is None else
                            tuple(sample for sample in header if sample in requested_samples))
                if samples is None:
                    samples = selected
                    if not samples or requested_samples is not None and requested_samples != set(samples):
                        raise ValueError(f'Requested VCF samples are missing from {source_file}')
                    a, b = (np.empty(len(samples), dtype=np.int8) for _ in range(2))
                    phased = np.empty(len(samples), dtype=np.bool_)
                    phase_sets = np.empty(len(samples), dtype=np.int32)
                elif selected != samples and set(selected) != set(samples):
                    raise ValueError('VCF files must contain the same set of selected samples')
                if header == samples:
                    source_indices = sample_order = np.arange(len(samples), dtype=np.intp)
                else:
                    output_index = {sample: i for i, sample in enumerate(samples)}
                    sample_order = np.array([output_index.get(sample, -1) for sample in header], dtype=np.intp)
                    source_indices = np.empty(len(samples), dtype=np.intp)
                    included = np.flatnonzero(sample_order >= 0)
                    source_indices[sample_order[included]] = included

                for line in lines:
                    if not line or line[:1] in (b'\n', b'\r', b'#'):
                        continue
                    n_records += 1
                    fields, start = _variant_fields(line)
                    contig, position, sid, ref, alt, _, status, info, format_ = fields
                    contig = contig.removeprefix(b'chr')
                    chromosome = allowed_chroms.get(contig)
                    if chromosome is None or (expected_chrom is not None and chromosome != expected_chrom):
                        continue
                    if not ref.isalpha() or not alt.isalpha() or ref == alt or pass_only and status != b'PASS':
                        continue
                    if sid == b'.' or not sid:
                        sid = b':'.join((contig, position, ref, alt))
                    if requested_variants is not None and sid not in requested_variants:
                        continue
                    if quality_searches and not _passes_info_quality(info, quality_searches):
                        n_quality_removed += 1
                        continue
                    if min_info_maf:
                        maf = _info_maf(info)
                        if maf is not None and maf < min_info_maf and not math.isclose(maf, min_info_maf):
                            continue
                    sample_fields = format_.split(b':')
                    if b'GT' not in sample_fields:
                        continue
                    ps_field = sample_fields.index(b'PS') if phase_possible and b'PS' in sample_fields else -1
                    psl_field = sample_fields.index(b'PSL') if phase_possible and b'PSL' in sample_fields else -1
                    end = len(line)
                    if line[-1:] == b'\n':
                        end -= 1
                    if line[end - 1:end] == b'\r':
                        end -= 1
                    decoded = _decode_samples(np.frombuffer(line, dtype=np.uint8), start, end,
                                              sample_fields.index(b'GT'), ps_field, psl_field, source_indices,
                                              sample_order, a, b, phased, phase_sets)
                    if decoded:
                        reasons = ('expected diploid biallelic GT fields', 'sample count differs from header',
                                   'invalid phase-set value')
                        raise ValueError(f'VCF {source_file}:{contig.decode()}:{position.decode()}: {reasons[decoded - 1]}')
                    if min_maf:
                        alt_count, called = _allele_counts(a, b)
                        if not called or min(alt_count, called - alt_count) / called < min_maf:
                            continue
                    pos = int(position)
                    if not 0 < pos <= 2**31 - 1:
                        raise ValueError(f'VCF position out of range: {source_file}:{contig.decode()}:{pos}')
                    if phase_possible:
                        if chromosome not in phase_sets_by_chrom:
                            phase_sets_by_chrom[chromosome] = np.full(len(samples), -2, dtype=np.int32)
                        has_phase, conflict = _track_phasing(a, b, phased, phase_sets,
                                                            phase_sets_by_chrom[chromosome])
                        saw_phase |= has_phase
                        if conflict:
                            reasons = ('unphased heterozygous calls',
                                       'heterozygous calls in multiple phase sets for one sample',
                                       'allele-specific phase sets (FORMAT/PSL)')
                            phase_fallback = reasons[conflict - 1]
                            if output == 'phased':
                                raise ValueError(f'Cannot write phased BGEN for {source_file}:{chromosome}:{pos}: '
                                                 f'{phase_fallback}')
                            phase_possible = False
                            phase_sets_by_chrom.clear()
                            phased_handles.close()
                            bgen.clear()
                        else:
                            if chromosome not in bgen:
                                bgen[chromosome] = _PhasedWriter(phased_handles, work / f'chr_{chromosome}', samples)
                            bgen[chromosome].write(contig, pos, sid, ref, alt, a, b)
                    if output != 'phased':
                        if chromosome not in bed:
                            bed[chromosome] = _BedWriter(stack, work / f'chr_{chromosome}', samples)
                        bed[chromosome].write(contig, pos, sid, ref, alt, a, b)
                    n_variants += 1

        if not n_variants:
            if n_quality_removed:
                raise ValueError('No biallelic diploid VCF variants remain after filtering '
                                 f'({n_quality_removed:,} failed INFO quality thresholds)')
            raise ValueError('No biallelic diploid VCF variants remain after filtering')
        use_bgen = phase_possible and (output == 'phased' or saw_phase)
        for writer in bgen.values():
            writer.finish()
        stack.close()
        chroms = tuple(sorted(bgen if use_bgen else bed))
        split_chromosomes = split_chromosomes or len(chroms) > 1 or '@' in str(source)
        suffixes = ('.bgen', '.sample') if use_bgen else ('.bed', '.bim', '.fam')
        output_files = []
        for chrom in chroms:
            stem = f'chr_{chrom}' if split_chromosomes else 'genotypes'
            for suffix in suffixes:
                original = work / f'chr_{chrom}{suffix}'
                if original.exists():
                    target = destination / (stem + suffix)
                    if target.exists():
                        raise FileExistsError(f'VCF conversion output already exists: {target}')
                    output_files.append((original, target))
        for original, target in output_files:
            original.rename(target)
    prefix = str(destination / ('chr_@' if split_chromosomes else 'genotypes'))
    return VCFConversion(None if use_bgen else prefix, prefix if use_bgen else None,
                         chroms, len(samples), n_variants, n_records, phase_fallback,
                         tuple(applied_quality.items()), n_quality_removed)
