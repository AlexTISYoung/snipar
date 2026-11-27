"""Encode A1-counted BED hard calls as unphased diploid BGEN for the benchmark."""

import argparse
from pathlib import Path
from time import perf_counter

from bed_reader import open_bed
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bed', type=Path)
    parser.add_argument('bgen', type=Path)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--batch', type=int, default=1024)
    parser.add_argument('--a2-first', action='store_true', help='Encode A2 as BGEN allele 1 (REF in the EUR fixture).')
    args = parser.parse_args()
    if args.threads < 1 or args.batch < 1:
        parser.error('--threads and --batch must be positive')
    for path in (args.bgen, Path(str(args.bgen) + '.bgi')):
        if path.exists():
            parser.error(f'{path} already exists')
    try:
        from bgen import BgenReader, BgenWriter
    except ModuleNotFoundError as error:
        if error.name != 'bgen':
            raise
        parser.error('Install the optional bgen package to write this fixture')

    started = perf_counter()
    args.bgen.parent.mkdir(parents=True, exist_ok=True)
    bed = open_bed(args.bed, count_A1=True, num_threads=args.threads)
    samples, variants = bed.shape
    first, second = ((bed.allele_2, bed.allele_1) if args.a2_first else (bed.allele_1, bed.allele_2))
    # The probabilities are P(first/first), P(first/second) and P(second/second).
    states = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0], [np.nan, np.nan, np.nan]], dtype=np.float32)
    with BgenWriter(args.bgen, samples, list(bed.iid), compression='zstd', layout=2) as writer:
        for left in range(0, variants, args.batch):
            right = min(variants, left + args.batch)
            hard_calls = bed.read((slice(None), slice(left, right)), dtype='int8', order='F')
            if args.a2_first:
                np.subtract(2, hard_calls, out=hard_calls, where=hard_calls >= 0)
            hard_calls[hard_calls < 0] = 3
            for j in range(left, right):
                writer.add_variant(str(bed.sid[j]), str(bed.sid[j]), str(bed.chromosome[j]),
                                   int(bed.bp_position[j]), [str(first[j]), str(second[j])],
                                   states[hard_calls[:, j - left]], ploidy=2, phased=False, bit_depth=8)
    print(f'Wrote {variants:,} variants and {samples:,} samples in {perf_counter() - started:.3f} s', flush=True)

    with BgenReader(args.bgen) as reader:
        np.testing.assert_array_equal(reader.samples, bed.iid)
        assert len(reader) == variants
        checked = np.linspace(0, variants - 1, min(100, variants), dtype=int)
        original = bed.read((slice(None), checked), dtype='float32')
        if args.a2_first:
            original = 2 - original
        for col, index in enumerate(checked):
            variant = reader[int(index)]
            assert not variant.is_phased
            assert variant.alleles == [first[index], second[index]]
            probabilities = variant.probabilities
            actual = 2 * probabilities[:, 0] + probabilities[:, 1]
            np.testing.assert_array_equal(actual, original[:, col])
    print(f'Checked sample order, allele order and every call for {len(checked)} evenly spaced variants')


if __name__ == '__main__':
    main()
