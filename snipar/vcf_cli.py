from argparse import Namespace
from contextlib import contextmanager
from functools import wraps
from pathlib import Path
from tempfile import TemporaryDirectory


def add_vcf_options(parser, *, phased_only=False):
    calls = 'Phased' if phased_only else 'Phased or unphased'
    parser.add_argument('--vcf', help=f'{calls} diploid VCF (.vcf, .vcf.gz or .vcf.bgz). '
                                      'Use @ in a filename to substitute chromosome numbers.')
    add_vcf_filters(parser)


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
        **options,
    )
    output = 'phased BGEN' if result.bgen is not None else 'BED'
    print(f'Prepared {result.n_variants:,} of {result.n_records:,} VCF records and '
          f'{result.n_samples:,} samples as {output}', flush=True)
    if result.phase_fallback:
        print(f'Using unphased genotypes because the VCF contains {result.phase_fallback}', flush=True)
    return result


@contextmanager
def vcf_inputs(args, *, output='bed', variant_ids=None, chromosomes=None):
    """Provide seekable genotype inputs for the lifetime of a CLI invocation."""
    if getattr(args, 'vcf', None) is None:
        if (getattr(args, 'vcf_variants', None) is not None or getattr(args, 'vcf_pass_only', False)
                or getattr(args, 'vcf_min_info_maf', 0.) or getattr(args, 'vcf_min_maf', 0.)):
            raise ValueError('VCF filters require --vcf')
        yield args
        return
    if getattr(args, 'bed', None) is not None or getattr(args, 'bgen', None) is not None:
        raise ValueError('Provide --vcf by itself as the observed genotype input')
    if getattr(args, 'fam', None) is not None or getattr(args, 'bim', None) is not None:
        raise ValueError('--fam and --bim require BED or BGEN input')
    if chromosomes is None:
        chromosomes = getattr(args, 'chr_range', None)
    chromosomes = None if chromosomes is None else tuple(chromosomes)
    split = bool(chromosomes) or any('@' in str(getattr(args, name, None))
                                     for name in ('vcf', 'imp', 'ibd', 'out'))
    with TemporaryDirectory(prefix='snipar-vcf-') as directory:
        converted = prepare_vcf(args, directory, variant_ids=variant_ids, output=output,
                                chromosomes=chromosomes, split_chromosomes=split)
        yield Namespace(**(vars(args) | {'bed': converted.bed, 'bgen': converted.bgen}))


def with_vcf(*, output='bed'):
    def decorate(main):
        @wraps(main)
        def run(args):
            with vcf_inputs(args, output=output) as resolved:
                return main(resolved)
        return run
    return decorate
