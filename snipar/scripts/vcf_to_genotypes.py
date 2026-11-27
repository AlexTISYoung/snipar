#!/usr/bin/env python
"""Prepare reusable BED or phased BGEN genotypes from VCF files."""

import argparse
from pathlib import Path

from snipar.docgen import get_parser_doc
from snipar.numrange import NumRangeAction, parseNumRange
from snipar.genotype_cli import add_quality_options, add_vcf_filters, prepare_vcf, read_id_list


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('vcf', help='Input VCF (.vcf, .vcf.gz or .vcf.bgz). Use @ to substitute chromosome numbers.')
parser.add_argument('out', type=Path, help='Directory for the prepared files.')
parser.add_argument('--format', choices=('bed', 'auto', 'phased'), default='bed',
                    help='Write BED (default), preserve consistent phase when possible, or require phased BGEN.')
parser.add_argument('--samples', type=Path, help='Keep sample IDs listed one per line, in source order.')
parser.add_argument('--chr_range', type=parseNumRange, nargs='*', action=NumRangeAction,
                    help='Chromosomes to keep. Supply a range or individual chromosome numbers.')
add_vcf_filters(parser)
add_quality_options(parser)
__doc__ += '\n\nArgs:\n' + get_parser_doc(parser)


def main(argv=None):
    args = parser.parse_args(argv)
    try:
        samples = None if args.samples is None else read_id_list(args.samples)
        prepared = prepare_vcf(args, args.out, chromosomes=args.chr_range, output=args.format,
                               sample_ids=samples, split_chromosomes=bool(args.chr_range))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f'Genotype prefix {prepared.bed or prepared.bgen}')


if __name__ == '__main__':
    main()
