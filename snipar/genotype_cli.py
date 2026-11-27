from argparse import Namespace
from contextlib import contextmanager
from functools import wraps
from math import isfinite
from pathlib import Path
from tempfile import TemporaryDirectory


def add_vcf_options(parser, *, phased_only=False, bgen_quality=False):
    calls = 'Phased' if phased_only else 'Phased or unphased'
    parser.add_argument('--vcf', help=f'{calls} diploid VCF (.vcf, .vcf.gz or .vcf.bgz). '
                                      'Use @ in a filename to substitute chromosome numbers.')
    add_vcf_filters(parser)
    add_quality_options(parser, bgen=bgen_quality)


def add_vcf_filters(parser):
    parser.add_argument('--vcf-variants', type=Path,
                        help='Keep VCF variant IDs listed one per line before reading genotypes.')
    parser.add_argument('--vcf-pass-only', action='store_true',
                        help='Keep only VCF records whose FILTER is PASS.')
    parser.add_argument('--vcf-min-info-maf', type=float, default=0.,
                        help='Early filter using VCF INFO/AC and INFO/AN or INFO/AF. '
                             'Records without usable frequencies are retained. This uses the VCF annotation cohort.')
    parser.add_argument('--vcf-min-maf', type=float, default=0.,
                        help='Minimum MAF among called VCF genotypes before preparing files. '
                             'This uses all selected VCF samples and is independent of the analysis MAF filter.')


def add_quality_options(parser, *, bgen=False):
    parser.add_argument('--min-quality', action='append', metavar='FIELD[=MIN]',
                        help='Minimum variant quality for an INFO tag or BGEN quality column. '
                             'Repeat to require every threshold. Default minimum: 0.99.')
    parser.add_argument('--no-quality-filter', action='store_true',
                        help='Disable automatic variant quality filtering.')
    if bgen:
        parser.add_argument('--bgen-quality', metavar='FILE',
                            help='BGEN variant quality table or UK Biobank MFI file. '
                                 'Use @ for per-chromosome tables.')
        parser.add_argument('--bgen-quality-id', metavar='COLUMN',
                            help='Variant ID column in a headered BGEN quality table.')


def quality_thresholds(args):
    specified = getattr(args, 'min_quality', None)
    if getattr(args, 'no_quality_filter', False):
        if specified:
            raise ValueError('Use either --min-quality or --no-quality-filter')
        return {}
    if not specified:
        return None
    from snipar.variant_quality import DEFAULT_MIN_QUALITY

    filters = {}
    for rule in specified:
        field, separator, minimum = rule.partition('=')
        try:
            cutoff = float(minimum) if separator else DEFAULT_MIN_QUALITY
        except ValueError as error:
            raise ValueError(f'Expected FIELD or FIELD=MIN for --min-quality, got {rule!r}') from error
        if not field or any(character.isspace() for character in field) or not isfinite(cutoff):
            raise ValueError(f'Expected a field and finite minimum for --min-quality, got {rule!r}')
        if field in filters:
            raise ValueError(f'Duplicate --min-quality field {field!r}')
        filters[field] = cutoff
    return filters


def read_id_list(path):
    with Path(path).open() as source:
        return {identifier for line in source if (identifier := line.strip())}


def prepare_vcf(args, destination, *, variant_ids=None, chromosomes=None, **options):
    from snipar.vcf import convert_vcf

    ids_path = getattr(args, 'vcf_variants', None)
    if ids_path is not None:
        requested = read_id_list(ids_path)
        variant_ids = requested if variant_ids is None else requested.intersection(variant_ids)
    result = convert_vcf(
        args.vcf, destination, variant_ids=variant_ids,
        chromosomes=None if chromosomes is None else (tuple(map(int, chromosomes)) or None),
        min_maf=getattr(args, 'vcf_min_maf', 0.),
        min_info_maf=getattr(args, 'vcf_min_info_maf', 0.),
        pass_only=getattr(args, 'vcf_pass_only', False),
        min_quality=quality_thresholds(args),
        **options,
    )
    output = 'phased BGEN' if result.bgen is not None else 'BED'
    print(f'Prepared {result.n_variants:,} of {result.n_records:,} VCF records and '
          f'{result.n_samples:,} samples as {output}', flush=True)
    if result.phase_fallback:
        print(f'Using unphased genotypes because the VCF contains {result.phase_fallback}', flush=True)
    if result.quality_filters:
        print(f'VCF quality filtering excluded {result.n_quality_removed:,} records '
              f'({_quality_description(result.quality_filters)})', flush=True)
    return result


def _quality_description(rules):
    return ', '.join(f'{field} >= {minimum:g}' for field, minimum in rules)


@contextmanager
def _bgen_inputs(args, table, thresholds, chromosomes):
    from snipar.bgen_filter import BGEN_SUBSET_SUFFIX, filter_bgen
    from snipar.utilities import parse_obsfiles
    from snipar.variant_quality import read_bgen_quality

    if getattr(args, 'bed', None) is not None or getattr(args, 'bim', None) is not None:
        raise ValueError('BGEN quality filtering requires BGEN as the only variant input')
    prefix = str(args.bgen).removesuffix('.bgen')
    table = str(table)
    split = '@' in prefix
    by_chromosome = '@' in table
    if by_chromosome and not split:
        raise ValueError('A chromosome pattern for --bgen-quality requires @ in --bgen')
    files, chroms = parse_obsfiles(prefix, 'bgen', chromosomes=chromosomes)
    id_column = getattr(args, 'bgen_quality_id', None)
    shared_quality = None if by_chromosome else read_bgen_quality(table, thresholds, id_column=id_column)
    with TemporaryDirectory(prefix='snipar-bgen-quality-') as directory:
        output = Path(directory) / (('chr_@' if split else 'genotypes') + BGEN_SUBSET_SUFFIX.removesuffix('.bgen'))
        for source, chromosome in zip(files, chroms):
            quality = (read_bgen_quality(table.replace('@', str(chromosome)), thresholds, id_column=id_column)
                       if by_chromosome else shared_quality)
            target = Path(str(output).replace('@', str(chromosome)) + '.bgen')
            filtered = filter_bgen(source, target, quality.passing_ids, known_ids=quality.known_ids)
            print(f'BGEN quality retained {filtered.n_kept:,} of {filtered.n_variants:,} variants in {source} '
                  f'({_quality_description(quality.filters.items())})', flush=True)
        del quality, shared_quality
        yield Namespace(**(vars(args) | {'bgen': str(output)}))


@contextmanager
def genotype_inputs(args, *, output='bed', variant_ids=None, chromosomes=None):
    """Provide seekable genotype inputs for the lifetime of a CLI invocation."""
    table = getattr(args, 'bgen_quality', None)
    if table is not None and (getattr(args, 'bgen', None) is None or getattr(args, 'vcf', None) is not None):
        raise ValueError('--bgen-quality requires --bgen as the observed genotype input')
    if getattr(args, 'bgen_quality_id', None) is not None and table is None:
        raise ValueError('--bgen-quality-id requires --bgen-quality')
    if chromosomes is None:
        chromosomes = getattr(args, 'chr_range', None)
    chromosomes = None if chromosomes is None else tuple(chromosomes) or None
    if getattr(args, 'vcf', None) is None:
        if (getattr(args, 'vcf_variants', None) is not None or getattr(args, 'vcf_pass_only', False)
                or getattr(args, 'vcf_min_info_maf', 0.) or getattr(args, 'vcf_min_maf', 0.)):
            raise ValueError('VCF filters require --vcf')
        thresholds = quality_thresholds(args)
        if thresholds and table is None:
            raise ValueError('--min-quality requires --vcf or a --bgen-quality table')
        if table is not None and thresholds != {}:
            with _bgen_inputs(args, table, thresholds, chromosomes) as resolved:
                yield resolved
        else:
            yield args
        return
    if getattr(args, 'bed', None) is not None or getattr(args, 'bgen', None) is not None:
        raise ValueError('Provide --vcf by itself as the observed genotype input')
    if getattr(args, 'fam', None) is not None or getattr(args, 'bim', None) is not None:
        raise ValueError('--fam and --bim require BED or BGEN input')
    split = bool(chromosomes) or any('@' in str(getattr(args, name, None))
                                     for name in ('vcf', 'imp', 'ibd', 'out'))
    with TemporaryDirectory(prefix='snipar-vcf-') as directory:
        converted = prepare_vcf(args, directory, variant_ids=variant_ids, output=output,
                                chromosomes=chromosomes, split_chromosomes=split)
        yield Namespace(**(vars(args) | {'bed': converted.bed, 'bgen': converted.bgen}))


def with_genotypes(*, output='bed'):
    def decorate(main):
        @wraps(main)
        def run(args):
            with genotype_inputs(args, output=output) as resolved:
                return main(resolved)
        return run
    return decorate
