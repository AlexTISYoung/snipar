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


def family_vcf_samples(args, available_ids, *, ids=None, sib=False, compute_controls=False,
                       expand_ids=None, regroup_pedigree=False, pedigree=None):
    """Select observed probands and the relatives needed by a family design."""
    import numpy as np

    def pedigrees():
        if args.imp is None:
            if pedigree is not None:
                yield pedigree
            elif args.pedigree is not None:
                ped = np.loadtxt(args.pedigree, dtype=str, ndmin=2)
                if ped.shape[1] < 4:
                    raise ValueError('Not enough columns in pedigree file')
                if regroup_pedigree:
                    from snipar.pedigree import get_sibpairs_from_ped

                    if ped.size and tuple(ped[0, :2]) == ('FID', 'IID'):
                        ped = ped[1:]
                    _, ped = get_sibpairs_from_ped(ped[:, :4], return_pairs=False)
                yield ped
            return

        import h5py

        prefix = str(args.imp)
        if '@' in prefix:
            chromosomes = getattr(args, 'chr_range', None)
            chromosomes = range(1, 23) if chromosomes is None else dict.fromkeys(map(int, chromosomes)) or range(1, 23)
            if '@' in str(args.vcf):
                from snipar.vcf import _vcf_files

                chromosomes = (chrom for chrom, _ in _vcf_files(args.vcf, chromosomes))
            files = (Path(prefix.replace('@', str(chrom)) + '.hdf5') for chrom in chromosomes)
        else:
            files = (Path(prefix + '.hdf5'),)
        for file in files:
            if file.is_file():
                try:
                    with h5py.File(file, 'r') as source:
                        ped = source['pedigree'].asstr(encoding='utf-8')[1:].astype(str)
                except (OSError, KeyError, TypeError, ValueError):
                    # A candidate chromosome can disappear after site filtering.
                    # Let the command validate only the chromosomes it analyzes.
                    yield None
                    return
                if ped.ndim != 2 or ped.shape[1] < 4:
                    yield None
                    return
                if not compute_controls:
                    controls = np.char.startswith(ped[:, 0], '_')
                    if controls.any():
                        ped = ped[~controls]
                yield ped

    requested = None if ids is None else set(ids)
    required = set() if requested is None else requested.copy()
    have_pedigree = False
    for ped in pedigrees():
        if ped is None:
            return None
        have_pedigree = True
        if requested is None:
            rows = ped
            required.update(ped[:, 1])
        else:
            probands = requested if expand_ids is None else set(expand_ids(ped, requested))
            required.update(probands)
            rows = ped[np.fromiter(map(probands.__contains__, ped[:, 1]), dtype=bool, count=len(ped))]
        required.update(rows[:, 2])
        required.update(rows[:, 3])
        if sib and requested is not None:
            families = set(rows[:, 0])
            members = np.fromiter(map(families.__contains__, ped[:, 0]), dtype=bool, count=len(ped))
            required.update(ped[members, 1])
    if not have_pedigree:
        return None
    required.intersection_update(available_ids)
    return required


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
def genotype_inputs(args, *, output='bed', variant_ids=None, chromosomes=None, vcf_samples=None):
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
    # Observed VCF MAF is independent of analysis sample selection.
    resolved = Namespace(**vars(args))
    select_samples = None
    if vcf_samples is not None and not getattr(args, 'vcf_min_maf', 0.):
        def select_samples(available_ids):
            return vcf_samples(resolved, available_ids)
    with TemporaryDirectory(prefix='snipar-vcf-') as directory:
        converted = prepare_vcf(resolved, directory, variant_ids=variant_ids, output=output,
                                chromosomes=chromosomes, split_chromosomes=split, select_samples=select_samples)
        resolved.bed, resolved.bgen = converted.bed, converted.bgen
        yield resolved


def with_genotypes(*, output='bed', vcf_samples=None):
    def decorate(main):
        @wraps(main)
        def run(args):
            with genotype_inputs(args, output=output, vcf_samples=vcf_samples) as resolved:
                return main(resolved)
        return run
    return decorate
