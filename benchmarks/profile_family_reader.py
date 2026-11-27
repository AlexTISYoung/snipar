"""Profile bounded family-design reads on a BED fixture from scale_workloads.py."""

import argparse
from contextlib import nullcontext
import os
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd

from snipar import read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture', type=Path, help='Directory from scale_workloads.py bed')
    parser.add_argument('--width', type=int, default=8, help='Number of variants per block')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--unprepared', action='store_true', help='Measure the public one-block API on the baseline')
    args = parser.parse_args()
    if min(args.width, args.threads) < 1:
        parser.error('--width and --threads must be positive')
    os.environ.setdefault('PST_NUM_THREADS', str(args.threads))
    path = args.fixture / 'chr_1.bed'
    ped = pd.read_csv(args.fixture / 'pedigree.txt', sep=r'\s+', dtype=str).to_numpy(dtype=str)
    options = dict(ped=ped, bedfile=path, ids=ped[:, 1], sib=True)
    if args.unprepared:
        variants = pd.read_csv(path.with_suffix('.bim'), sep=r'\s+', header=None,
                               dtype=str, usecols=[1])[1].to_numpy(dtype=str)
        context = nullcontext()
    else:
        context = read.prepare_gts_matrix(**options, num_threads=args.threads)
    started = perf_counter()
    with context as genotype_reader:
        if genotype_reader is not None:
            alignment = genotype_reader.match_variants(match_snp_ids=True)
            variants = alignment.sid
            print(f'Plan for {genotype_reader.n_samples:,} probands in {perf_counter() - started:.3f}s')
        for index in range(0, len(variants), args.width):
            block = slice(index, index + args.width)
            if genotype_reader is None:
                matrix = read.get_gts_matrix(**options, snp_ids=variants[block], match_snp_ids=True)
            else:
                matrix = genotype_reader.read(alignment.subset(block))
            values = matrix.gts.data
            print(f'Block {index // args.width + 1}, missing {np.isnan(values).sum()}, '
                  f'sum {np.nansum(values, dtype=np.float64):.12g}', flush=True)
            del matrix, values
    print(f'Total {perf_counter() - started:.3f}s')
    try:
        import resource
    except ImportError:
        return
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f'Peak RSS {rss * (1 if sys.platform == "darwin" else 1024) / 2**30:.3f} GiB')


if __name__ == '__main__':
    main()
