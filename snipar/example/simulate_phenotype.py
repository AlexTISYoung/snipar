"""Simulate an additive phenotype from observed genotypes in a PLINK BED file."""

import argparse
from pathlib import Path

import numpy as np
from pysnptools.snpreader import Bed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bed', default='whole_genome', help='PLINK BED filename or prefix')
    parser.add_argument('--outdir', type=Path, default=Path('.'), help='Output directory')
    parser.add_argument('--ncausal', type=int, default=1500, help='Number of causal variants')
    parser.add_argument('--h2', type=float, default=.75, help='Target additive heritability')
    parser.add_argument('--seed', type=int, help='Random seed')
    args = parser.parse_args(argv)
    bed = Bed(args.bed, count_A1=True)
    if not 0 <= args.h2 <= 1 or not 1 <= args.ncausal <= bed.sid_count:
        parser.error('h2 must be within [0, 1] and ncausal within the available variant count')
    rng = np.random.default_rng(args.seed)
    causal = np.sort(rng.choice(bed.sid_count, size=args.ncausal, replace=False))
    gts = bed[:, causal].read(dtype=np.float32).val
    observed = np.isfinite(gts)
    means = np.nansum(gts, axis=0) / np.maximum(observed.sum(axis=0), 1)
    np.copyto(gts, means, where=~observed)
    effects = rng.standard_normal(args.ncausal)
    genetic = gts @ effects
    sd = genetic.std()
    if args.h2 > 0 and sd == 0:
        raise ValueError('The selected variants have no additive genetic variance')
    scale = np.sqrt(args.h2) / sd if args.h2 else 0
    effects *= scale
    phenotype = genetic * scale + rng.normal(scale=np.sqrt(1 - args.h2), size=len(genetic))
    bim = np.loadtxt(Path(args.bed).with_suffix('.bim') if args.bed.endswith('.bed')
                     else f'{args.bed}.bim', dtype=str, ndmin=2)
    prefix = args.outdir / f'h2_{args.h2}'
    np.savetxt(f'{prefix}.fam', np.column_stack((bed.iid, phenotype)), fmt='%s')
    np.savetxt(f'{prefix}_causal.txt', np.column_stack((bim[causal], effects, np.zeros_like(effects))), fmt='%s')


if __name__ == '__main__':
    main()
