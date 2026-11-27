"""Generate repeatable large BED inputs for scaling benchmarks."""

import argparse
from pathlib import Path

import numpy as np


def bed_fixture(destination, families, variants):
    destination.mkdir(parents=True, exist_ok=True)
    stem = destination / 'chr_1'
    rng = np.random.default_rng(293)
    # Four people per family, with two probands followed by both parents.
    phenotype_rng = np.random.default_rng(40657)
    with (
        stem.with_suffix('.fam').open('w') as sample_file,
        (destination / 'pedigree.txt').open('w') as pedigree,
        (destination / 'phenotype.txt').open('w') as phenotype,
    ):
        pedigree.write('FID IID FATHER_ID MOTHER_ID\n')
        phenotype.write('FID IID phenotype\n')
        for start in range(0, families, 25_000):
            stop = min(families, start + 25_000)
            group = phenotype_rng.normal(size=(stop - start, 3))
            group[:, :2] += .5 * group[:, 2, None]
            for f, (first_y, second_y) in zip(range(start, stop), group[:, :2]):
                family = f'F{f:07}'
                first, second, father, mother = [f'{family}_{role}' for role in ('a', 'b', 'p', 'm')]
                sample_file.writelines(f'{family}\t{iid}\t0\t0\t0\t-9\n' for iid in (first, second, father, mother))
                pedigree.write(f'{family} {first} {father} {mother}\n{family} {second} {father} {mother}\n')
                phenotype.write(f'{family} {first} {first_y:.10f}\n{family} {second} {second_y:.10f}\n')
    with stem.with_suffix('.bim').open('w') as metadata, (destination / 'weights.txt').open('w') as weights:
        weights.write('SNP A1 A2 b\n')
        for j in range(variants):
            metadata.write(f'1\trs{j:08}\t{j * .01}\t{1000 + j}\tA\tG\n')
            weights.write(f'rs{j:08} A G {(-1.) ** j / (j + 1):.10g}\n')
    # Variant-major PLINK BED. Generate hard calls in fixed-size tiles rather
    # than allocating a sample-by-variant matrix to create the benchmark.
    codes = np.repeat(np.array([1, 0, 2, 3], dtype=np.uint8), [2, 24, 50, 24])
    with stem.with_suffix('.bed').open('wb') as output:
        output.write(b'\x6c\x1b\x01')
        for start in range(0, variants, 16):
            values = codes[rng.integers(0, 100, size=(min(16, variants - start), families, 4), dtype=np.uint8)]
            packed = values[..., 0] | (values[..., 1] << 2) | (values[..., 2] << 4) | (values[..., 3] << 6)
            packed.tofile(output)


def imputation_fixture(destination, families, variants):
    destination.mkdir(parents=True, exist_ok=True)
    stem = destination / 'chr_1'
    rng = np.random.default_rng(122)
    people = 2 * families
    with stem.with_suffix('.fam').open('w') as fam, (destination / 'pedigree.txt').open('w') as ped:
        ped.write('FID IID FATHER_ID MOTHER_ID\n')
        for first in range(0, families, 25_000):
            indices = range(first, min(families, first + 25_000))
            fam.write(''.join(f'F{i}\tC{i}\tP{i}\t0\t0\t-9\nF{i}\tP{i}\t0\t0\t1\t-9\n' for i in indices))
            ped.write(''.join(f'F{i} C{i} P{i} 0\n' for i in indices))
    with stem.with_suffix('.bim').open('w') as bim:
        bim.writelines(f'1\tv{i}\t0\t{i+1}\tA\tG\n' for i in range(variants))
    encoding = np.array([3, 2, 0], dtype=np.uint8)
    padded = (people + 3) // 4 * 4
    with stem.with_suffix('.bed').open('wb') as bed:
        bed.write(b'\x6c\x1b\x01')
        for start in range(0, variants, 16):
            width = min(16, variants - start)
            inherited, untransmitted, other = (rng.random((width, families)) < .32 for _ in range(3))
            transmitted = inherited.astype(np.uint8)
            calls = np.zeros((width, padded), dtype=np.uint8)
            calls[:, :people:2] = encoding[transmitted + other]
            calls[:, 1:people:2] = encoding[transmitted + untransmitted]
            g = calls.reshape(width, padded // 4, 4)
            packed = g[..., 0] | (g[..., 1] << 2) | (g[..., 2] << 4) | (g[..., 3] << 6)
            packed.tofile(bed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    fixture = commands.add_parser('bed')
    fixture.add_argument('output', type=Path)
    fixture.add_argument('--families', type=int, default=20000)
    fixture.add_argument('--variants', type=int, default=1000)
    impute = commands.add_parser('impute-bed')
    impute.add_argument('output', type=Path)
    impute.add_argument('--families', type=int, default=50_000)
    impute.add_argument('--variants', type=int, default=1_024)
    args = parser.parse_args()
    if args.families < 1 or args.variants < 1:
        parser.error('--families and --variants must be positive')
    generate = bed_fixture if args.command == 'bed' else imputation_fixture
    generate(args.output, args.families, args.variants)


if __name__ == '__main__':
    main()
