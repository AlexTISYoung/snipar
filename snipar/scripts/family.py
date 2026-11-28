#!/usr/bin/env python
"""Infer sibling IBD, impute parental genotypes, or export parent-of-origin haplotypes.

Use ``--impute`` to infer IBD automatically for genotyped sibling pairs and
impute missing parents. Use ``--ibd`` to infer IBD on its own. Phased inputs
use the phased sibling model and duoHMM by default.

Args:
@parser@
"""

import argparse
from itertools import combinations
import logging
from pathlib import Path
import re

import numpy as np
import pandas as pd

from snipar.genotype_cli import add_quality_options, add_vcf_filters, genotype_inputs
from snipar.imputation.preprocess_data import prepare_imputation_pedigree
from snipar.numrange import NumRangeAction, parseNumRange
from snipar.observed import open_observed
from snipar.pedigree import create_pedigree, get_sibpairs_from_king, get_sibpairs_from_ped
from snipar.scripts import ibd, impute
from snipar.utilities import get_parser_doc, outfile_name

_PEDIGREE_COLUMNS = ('FID', 'IID', 'FATHER_ID', 'MOTHER_ID')
_EXPORT_SUFFIXES = ('.vcf.gz', '.vcf', '.bgen')


def _relationships(args):
    if args.pedigree:
        values = np.loadtxt(args.pedigree, dtype=str, ndmin=2)
        if not values.size or values.shape[1] < 4:
            raise ValueError('A pedigree requires FID, IID, FATHER_ID, and MOTHER_ID columns')
        values = values[:, :4]
        if tuple(values[0]) == _PEDIGREE_COLUMNS:
            values = values[1:]
    elif args.agesex:
        values = create_pedigree(args.king, args.agesex).loc[:, list(_PEDIGREE_COLUMNS)].to_numpy(dtype=str)
    else:
        values = None
    if values is not None and args.pedigree_nan != '0':
        values[:, 2:4][values[:, 2:4] == args.pedigree_nan] = '0'
    pairs = get_sibpairs_from_king(args.king) if args.king else get_sibpairs_from_ped(values)[0]
    if args.impute and args.pedigree:
        # Imputation also recognizes a known parent shared within a family
        # when the other parent's ID is absent.
        candidates = values[(values[:, 2] == '0') ^ (values[:, 3] == '0')]
        groups = {}
        for fid, child, father, mother in candidates:
            if not fid.startswith('_'):
                groups.setdefault((fid, father, mother), []).append(child)
        extra = [pair for children in groups.values() for pair in combinations(children, 2)]
        if extra:
            pairs = np.asarray(extra) if pairs is None else np.concatenate((pairs, extra))
    pairs = np.empty((0, 2), dtype=str) if pairs is None else pairs
    return pairs, values


def _vcf_samples(args, available, relationships, pedigree):
    pairs, values = relationships
    available = set(available)
    infer_ibd = args.ibd or (args.impute and args.ibd_input is None)
    calibrate = infer_ibd and args.p_error is None
    phase_siblings = infer_ibd and ibd.phase_mode(args) is not False and not args.no_duohmm
    selected = set()
    if infer_ibd:
        for first, second in pairs:
            if first in available and second in available:
                selected.update((first, second))
    if values is not None:
        offspring = selected.copy()
        for _, child, father, mother in values:
            parents = available.intersection((father, mother)) - {'0'}
            if parents and ((phase_siblings and child in offspring)
                            or (child in available and (args.parent_origin_out or calibrate))):
                selected.add(child)
                selected.update(parents)
    if args.impute:
        prepared = prepare_imputation_pedigree(pedigree, available, control=args.c)
        args._vcf_pedigree = prepared
        selected.update(iid.decode() for iid in prepared.sample_ids)
    if not selected:
        raise ValueError('No genotyped relatives for the requested family analysis')
    return selected


def _detected_chromosomes(args):
    bedfiles, bgenfiles, chromosomes = ibd.input_files(args)
    if chromosomes[0]:
        return chromosomes.tolist()
    with open_observed(bedfiles[0], bgenfiles[0]) as source:
        labels = np.unique(source.reader.chromosome if source.bed else source.reader.chromosomes)
    match = re.fullmatch(r'(?:chr)?([0-9]+)(?:\.0+)?', str(labels[0]), re.IGNORECASE) if len(labels) == 1 else None
    if match is None:
        raise ValueError('Supply one chromosome per genotype file, or use --chrom to label a single file')
    return [int(match[1])]


def _export_targets(pattern, suffix, chromosomes):
    if '@' in pattern:
        return [Path(pattern.replace('@', str(chrom))) for chrom in chromosomes]
    if len(chromosomes) == 1:
        return [Path(pattern)]
    return [Path(outfile_name(pattern.removesuffix(suffix), suffix, chrom)) for chrom in chromosomes]


def main(args):
    """Run the selected operations with a single genotype preparation."""
    if not (args.ibd or args.impute or args.parent_origin_out):
        raise ValueError('Select --ibd, --impute, or --parent-origin-out')
    if args.ibd_input and (args.ibd or not args.impute):
        raise ValueError('--ibd-input is for --impute with precomputed IBD')
    if args.ibd_is_king and not args.ibd_input:
        raise ValueError('--ibd-is-king requires --ibd-input')
    if args.parent_origin and not args.ibd_input:
        raise ValueError('--parent-origin requires --ibd-input')
    if args.parent_origin_out and (args.bed or args.no_duohmm):
        raise ValueError('--parent-origin-out requires phased genotypes and duoHMM')
    if not args.pedigree and not args.king:
        raise ValueError('Supply --pedigree or --king')
    if (args.impute or args.parent_origin_out) and not (args.pedigree or args.agesex):
        raise ValueError('Imputation and origin export require --pedigree or both --king and --agesex')
    infer_ibd = args.ibd or (args.impute and args.ibd_input is None)
    suffix = next((ending for ending in _EXPORT_SUFFIXES if str(args.parent_origin_out).endswith(ending)), None)
    if args.parent_origin_out and suffix is None:
        raise ValueError('--parent-origin-out must end in .vcf, .vcf.gz, or .bgen')

    pairs, values = _relationships(args)
    relationships = pairs, values
    pedigree = None if values is None else pd.DataFrame(values, columns=_PEDIGREE_COLUMNS)
    args = argparse.Namespace(**vars(args))
    for kind in ('bed', 'bgen'):
        path = getattr(args, kind)
        if path is not None:
            setattr(args, kind, str(path).removesuffix(f'.{kind}'))
    original = argparse.Namespace(**vars(args))
    genotype_format = 'phased' if args.phased_ibd else 'auto'
    select_samples = lambda resolved, available: _vcf_samples(resolved, available, relationships, pedigree)
    with genotype_inputs(args, output=genotype_format, vcf_samples=select_samples,
                         split_chromosomes=True) as prepared:
        results = ()
        if infer_ibd or args.parent_origin_out:
            ibd_args = argparse.Namespace(**(vars(prepared) | {
                '_vcf_relationships': relationships, 'batches': 1,
                'posterior': args.posterior or (args.impute and args.ibd_posterior is None),
            }))
            results = ibd.run_prepared(ibd_args, all_parental_origin=bool(args.parent_origin_out or args.impute),
                                       allow_no_siblings=True, origin_only=not infer_ibd)
            if args.ibd and not any(result.stem for result in results):
                raise ValueError('No genotyped full-sibling pairs for --ibd')
        chromosomes = [result.chromosome for result in results] if results else _detected_chromosomes(prepared)
        if args.parent_origin_out:
            from snipar.imputation.origin_export import export_parental_haplotypes

            targets = _export_targets(args.parent_origin_out, suffix, chromosomes)
            if not any(result.origin for result in results):
                raise ValueError('No phased genotyped parent–offspring pairs for origin export')
            inputs = set()
            if original.vcf:
                from snipar.vcf import _vcf_files

                inputs.update(path.resolve() for _, path in _vcf_files(original.vcf, chromosomes))
            if original.bgen:
                inputs.update(Path(original.bgen.replace('@', str(chrom)) + '.bgen').resolve()
                              for chrom in chromosomes)
            if any(target.resolve() in inputs for target in targets):
                raise ValueError('Parent-of-origin output must use a different path from the input genotypes')
            for result, target in zip(results, targets):
                if result.origin:
                    export_parental_haplotypes(target, result.chromosome, values, [result.origin])
                    print(f'Wrote parent-of-origin haplotypes to {target}')
        if args.impute:
            impute_args = argparse.Namespace(**(vars(prepared) | {
                'ibd': args.ibd_input, 'no_parent_origin': args.no_duohmm, 'pedigree_nan': '0',
                'parent_origin_map': args.parent_origin_map or args.map,
            }))
            count = impute.run_prepared(
                impute_args, pedigree=pedigree, chromosomes=chromosomes,
                ibd_prefixes={result.chromosome: result.stem for result in results} if infer_ibd else None,
                origin_paths={result.chromosome: result.origin for result in results}
                if (infer_ibd or args.parent_origin_out) and args.parent_origin is None else None,
            )
            if not count and not (args.ibd or args.parent_origin_out):
                raise ValueError('No families have the required sibling or parental genotypes for imputation')


parser = argparse.ArgumentParser(description='Sibling IBD, parental imputation, and parent-of-origin haplotypes.')
parser.add_argument('--ibd', action='store_true', help='Infer sibling IBD and write segments and relatedness.')
parser.add_argument('--impute', action='store_true',
                    help='Impute missing parental genotypes. Automatically infer IBD for sibling pairs.')
parser.add_argument('--parent-origin-out', metavar='FILE',
                    help='Write paternal|maternal offspring haplotypes to .vcf, .vcf.gz or .bgen. Use @ for chromosome numbers.')
parser.add_argument('--ibd-input', metavar='PREFIX', help='Use precomputed IBD segments for imputation.')
parser.add_argument('--ibd-posterior', metavar='FILE', help='Use precomputed marginal IBD posteriors for imputation.')
parser.add_argument('--ibd-is-king', '--ibd_is_king', action='store_true', help='Read --ibd-input segments in KING format.')
genotypes = parser.add_mutually_exclusive_group(required=True)
genotypes.add_argument('--bgen', help='Phased or unphased BGEN prefix. Use @ for chromosome numbers.')
genotypes.add_argument('--bed', help='PLINK BED prefix. Use @ for chromosome numbers.')
genotypes.add_argument('--vcf', help='Phased or unphased diploid VCF. Use @ for chromosome numbers.')
add_vcf_filters(parser)
add_quality_options(parser, bgen=True)
relationships = parser.add_mutually_exclusive_group(required=True)
relationships.add_argument('--pedigree', help='Pedigree with FID, IID, FATHER_ID and MOTHER_ID columns.')
relationships.add_argument('--king', help='KING kinship file, optionally with --agesex for pedigree reconstruction.')
parser.add_argument('--agesex')
parser.add_argument('--pedigree_nan', default='0', help='Missing parent ID in the pedigree (default 0).')
parser.add_argument('--out', default='family', help='Output prefix. Replace @ with the chromosome number.')
parser.add_argument('--chr_range', nargs='*', type=parseNumRange, action=NumRangeAction)
parser.add_argument('--chrom', type=int, help='Chromosome number for a single genotype file.')
parser.add_argument('--map', help='Genetic map in SHAPEIT format. Use @ for chromosome numbers.')
parser.add_argument('--p_error', type=float, help='Genotyping error probability. Estimated when informative duos are available.')
model = parser.add_mutually_exclusive_group()
model.add_argument('--phased-ibd', action='store_true', help='Require input phase for sibling IBD.')
model.add_argument('--unphased-ibd', action='store_true', help='Use dosages for sibling IBD.')
model.add_argument('--joint', action='store_true', help='Use the joint dosage model for complete sibships of size three or four.')
origin = parser.add_mutually_exclusive_group()
origin.add_argument('--no-duohmm', '--no_parent_origin', action='store_true', help='Disable duoHMM parent-of-origin inference.')
origin.add_argument('--parent-origin', metavar='FILE', help='Use precomputed duoHMM labels for imputation with --ibd-input.')
parser.add_argument('--phase-switch-rate', type=float, default=.002, help='Switch probability per retained heterozygote.')
parser.add_argument('--min_maf', type=float, default=.01, help='IBD minimum minor allele frequency.')
parser.add_argument('--max_missing', type=float, default=5, help='IBD maximum missing genotype percentage.')
parser.add_argument('--max_error', type=float, default=.01, help='IBD maximum per-variant genotyping error probability.')
parser.add_argument('--min_length', type=float, default=.01, help='IBD minimum segment length for smoothing in cM.')
parser.add_argument('--posterior', action='store_true', help='Write marginal IBD probabilities. Enabled when inferring IBD for imputation.')
parser.add_argument('--ibdmatrix', action='store_true', help='Write the hard-call IBD matrix.')
parser.add_argument('--ld_out', action='store_true', help='Write variant LD scores used for IBD weighting.')
parser.add_argument('--threads', type=int, help='Maximum threads per chromosome.')
parser.add_argument('--processes', type=int, default=1, help='Maximum concurrent chromosomes during imputation.')
parser.add_argument('--chunks', type=int, default=1, help='Minimum imputation marker blocks per chromosome.')
parser.add_argument('--parent_origin_map', help='Genetic map for locally fitted duoHMM calls.')
parser.add_argument('--pcs', help='Principal components used for imputation frequency estimates.')
parser.add_argument('--pc_num', type=int)
parser.add_argument('-find_optimal_pc', action='store_true')
parser.add_argument('-use_backup', action='store_true')
parser.add_argument('-silent_progress', action='store_true')
parser.add_argument('-c', action='store_true', help='Include masked controls for imputation.')
parser.add_argument('--start', type=int, help='Start imputation at this marker index.')
parser.add_argument('--end', type=int, help='Stop imputation before this marker index.')
parser.add_argument('--output_compression', choices=['gzip', 'lzf'])
parser.add_argument('--output_compression_opts', type=int)
parser.set_defaults(bim=None, fam=None)
__doc__ = __doc__.replace('@parser@', get_parser_doc(parser))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s: %(message)s')
    main(parser.parse_args())
