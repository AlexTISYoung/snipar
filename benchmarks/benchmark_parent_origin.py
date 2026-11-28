"""Compare phased SNIPar imputation with and without inferred parent of origin."""

import argparse
from contextlib import ExitStack
import gzip
import json
import os
from pathlib import Path
import platform
import subprocess

from bed_reader import open_bed
import h5py
import numpy as np

from prepare_1000g_eur import (
    BASES, UNPHASED_GT, read_diploid_gt, vcf_header, write_bed_genotypes,
    write_marker_map, write_vcf_genotypes, write_vcf_header,
)
from profile_1000g_eur import step_command, timed_command
from validate_1000g_eur import check_imputation


def prepare_duos(source, pedigree, source_map, out, start, end):
    """Hold out one parent in each of the disjoint EUR pedigree groups."""
    (out / 'truth').mkdir(parents=True)
    with ExitStack() as stack:
        panel = stack.enter_context(gzip.open(source, 'rb'))
        sample_ids = vcf_header(panel)
        sample_index = {iid: index for index, iid in enumerate(sample_ids)}
        with pedigree.open() as file:
            header = next(file).split()
            if header != 'FamilyID SampleID FatherID MotherID Sex Population Superpopulation'.split():
                raise ValueError('Expected the IGSR pedigree with population labels')
            eur = [row for line in file if (row := line.split())[-1] == 'EUR' and row[1] in sample_index]
        eur_ids = {row[1] for row in eur}
        trios = {}
        for fid, child, father, mother, *_ in sorted(eur):
            if father in eur_ids and mother in eur_ids:
                trios.setdefault(fid, (child, father, mother))
        selected = [person for trio in trios.values() for person in trio]
        if not selected or len(set(selected)) != len(selected):
            raise ValueError('Expected one complete disjoint trio per selected pedigree')

        indices = np.array([sample_index[person] for person in selected])
        rows = np.arange(len(trios))
        observed_parent = 1 + rows % 2
        observed = np.column_stack((3 * rows, 3 * rows + observed_parent)).ravel()
        parent_truth = np.column_stack((3 * rows + 1, 3 * rows + 2)).ravel()

        ped = stack.enter_context((out / 'pedigree.txt').open('w'))
        fam = stack.enter_context((out / 'truth/chr_22.fam').open('w'))
        ped.write('FID IID FATHER_ID MOTHER_ID\n')
        for (fid, (child, father, mother)), parent in zip(trios.items(), observed_parent):
            ped.write(f'{fid} {child} {father} {mother}\n')
            ped.write(f'{fid} {(father, mother)[parent - 1]} 0 0\n')
            fam.write(f'{fid} {father} 0 0 1 -9\n{fid} {mother} 0 0 2 -9\n')

        unphased = stack.enter_context(gzip.open(out / 'unphased.vcf.gz', 'wb', compresslevel=1))
        write_vcf_header(unphased, [selected[index] for index in observed])
        bed = stack.enter_context((out / 'truth/chr_22.bed').open('wb'))
        bed.write(b'\x6c\x1b\x01')
        bim = stack.enter_context((out / 'truth/chr_22.bim').open('w'))
        cm_map = np.loadtxt(source_map, skiprows=1, usecols=(0, 2))
        codes = UNPHASED_GT.copy()
        codes[2] = codes[1]
        n_variants = 0
        for line in panel:
            fields = line.split(b'\t', 9)
            position = int(fields[1])
            if position < start:
                continue
            if position > end:
                break
            if fields[3] not in BASES or fields[4] not in BASES:
                continue
            haplotypes = read_diploid_gt(fields, len(sample_ids), indices)
            observed_calls = haplotypes[observed]
            if np.any(observed_calls > 1):
                continue
            count = int(observed_calls.sum())
            if min(count, observed_calls.size - count) / observed_calls.size < .01:
                continue
            truth_calls = haplotypes[parent_truth]
            if np.any(truth_calls > 1):
                continue
            sid = f"22:{position}:{fields[3].decode()}:{fields[4].decode()}".encode()
            fields[2] = sid
            write_vcf_genotypes(unphased, fields, observed_calls, codes, count)
            write_bed_genotypes(bed, truth_calls)
            cm = np.interp(position, cm_map[:, 0], cm_map[:, 1])
            bim.write(f'22\t{sid.decode()}\t{cm:.9f}\t{position}\t{fields[4].decode()}\t{fields[3].decode()}\n')
            n_variants += 1
    if not n_variants:
        raise ValueError('No eligible markers in the selected interval')
    write_marker_map(out / 'map.txt', cm_map, out / 'truth/chr_22.bim')
    print(json.dumps(dict(chrom=22, start=start, end=end, families=len(trios),
                          observed_fathers=int(np.count_nonzero(observed_parent == 1)),
                          observed_mothers=int(np.count_nonzero(observed_parent == 2)), variants=n_variants)))


def phase_fixture(fixture, source_map, out, java, jar, threads):
    """Normalize genotype order and phase only the visible individuals with Beagle."""
    unphased = fixture / 'unphased.vcf.gz'
    if not unphased.exists():
        unphased = out / 'unphased.vcf.gz'
        with (gzip.open(fixture / 'unphased/chr_22.vcf.gz', 'rb') as source,
              gzip.open(unphased, 'wb', compresslevel=1) as destination):
            for line in source:
                if line.startswith(b'#'):
                    destination.write(line)
                    continue
                site = line.split(b'\t', 9)
                site[9] = site[9].replace(b'1/0', b'0/1')
                destination.write(b'\t'.join(site))
    beagle_map = out / 'beagle.map'
    map_rows = np.loadtxt(source_map, skiprows=1, usecols=(0, 2))
    with beagle_map.open('w') as file:
        for position, cm in map_rows:
            file.write(f'chr22 . {cm:.9f} {int(position)}\n')
    phased = out / 'beagle'
    subprocess.run([str(java), '-Xmx4g', '-jar', str(jar), f'gt={unphased}', f'map={beagle_map}',
                    f'out={phased}', 'impute=false', 'seed=20261127', f'nthreads={threads}'], check=True)


def read_phase_panel(path):
    with gzip.open(path, 'rb') as panel:
        samples = vcf_header(panel)
        indices = np.arange(len(samples))
        variants, calls = [], []
        for line in panel:
            fields = line.split(b'\t', 9)
            variants.append((fields[2].decode(), int(fields[1]), fields[3], fields[4]))
            calls.append(read_diploid_gt(fields, len(samples), indices).astype(np.int8))
    return samples, variants, np.ascontiguousarray(np.stack(calls).transpose(1, 0, 2))


def score_origin(calls, labels, expected, eligible):
    called = eligible & (labels >= 0)
    allele = np.take_along_axis(calls, labels[:, :, None], axis=2)[:, :, 0]
    correct = called & (allele == expected)
    return dict(eligible=int(eligible.sum()), called=int(called.sum()), correct=int(correct.sum()))


def check_origin(fixture, phased_vcf, threads, truth_kind):
    """Check simulated transmissions or informative genotypes from the masked trio parent."""
    from snipar.imputation.duohmm import call_parental_origin, recombination_probabilities

    samples, variants, phased = read_phase_panel(phased_vcf)
    sample_rows = {iid: index for index, iid in enumerate(samples)}
    pedigree = np.loadtxt(fixture / 'pedigree.txt', skiprows=1, dtype=str)
    pair_ids, parent_sex = [], []
    for _, child, father, mother in pedigree:
        observed = [(parent, sex) for sex, parent in enumerate((father, mother)) if parent in sample_rows]
        if child in sample_rows and len(observed) == 1:
            (parent, sex), = observed
            pair_ids.append((parent, child))
            parent_sex.append(sex)
    pairs = np.array([(sample_rows[parent], sample_rows[child]) for parent, child in pair_ids])
    ids, positions, _, _ = zip(*variants)
    bp, cm = np.loadtxt(fixture / 'map.txt', skiprows=1, unpack=True)
    recombination = recombination_probabilities(positions, np.interp(positions, bp, cm))
    child_hap, parent_hap, _ = call_parental_origin(phased, pairs, recombination, threads=threads)
    parent, child = phased[pairs[:, 0]], phased[pairs[:, 1]]
    parent_het, child_het = parent.sum(axis=2) == 1, child.sum(axis=2) == 1
    if truth_kind == 'simulation':
        truth_samples, truth_variants, truth_calls = read_phase_panel(fixture / 'chr_22.vcf.gz')
        if truth_variants != variants:
            raise ValueError('Statistical and simulated phase panels must cover the same markers in order')
        truth_rows = {iid: index for index, iid in enumerate(truth_samples)}
        truth = truth_calls[[truth_rows[iid] for iid in samples]]
        np.testing.assert_array_equal(truth.sum(axis=2), phased.sum(axis=2))
        expected = truth[pairs[:, 1], :, parent_sex]
        return dict(child_het=score_origin(child, child_hap, expected, child_het),
                    double_het=score_origin(child, child_hap, expected, parent_het & child_het),
                    parent_het=score_origin(parent, parent_hap, expected, parent_het))

    bed = open_bed(fixture / 'truth/chr_22.bed')
    truth_rows = {iid: index for index, iid in enumerate(bed.iid)}
    parents_by_child = {child_id: (father, mother) for _, child_id, father, mother in pedigree}
    withheld_rows = [truth_rows[parents_by_child[child_id][1 - sex]]
                     for (_, child_id), sex in zip(pair_ids, parent_sex)]
    truth_index = {sid: index for index, sid in enumerate(bed.sid)}
    cols = np.array([truth_index[sid] for sid in ids])
    withheld = bed.read(index=(withheld_rows, cols), dtype='float32')
    informative = parent_het & child_het & ((withheld == 0) | (withheld == 2))
    return dict(double_het=score_origin(child, child_hap, 1 - withheld / 2, informative))


def check_paired_support(out):
    with (h5py.File(out / 'statistical_phase_0.hdf5', 'r') as baseline,
          h5py.File(out / 'parent_origin_0.hdf5', 'r') as updated):
        for key in ('families', 'parental_status', 'bim_values'):
            np.testing.assert_array_equal(baseline[key][:], updated[key][:])
        old, new = baseline['imputed_par_gts'], updated['imputed_par_gts']
        both_missing = ~baseline['parental_status'][:, :2].astype(bool).any(axis=1)
        only_old, only_new = 0, 0
        for start in range(0, old.shape[1], 2048):
            a, b = old[:, start:start + 2048], new[:, start:start + 2048]
            old_finite, new_finite = np.isfinite(a), np.isfinite(b)
            only_old += int((old_finite & ~new_finite).sum())
            only_new += int((new_finite & ~old_finite).sum())
            np.testing.assert_array_equal(a[both_missing], b[both_missing])
        return dict(statistical_phase_only=only_old, parent_origin_only=only_new)


def compare(fixture, phased_vcf, out, threads, repeat, ibd=None, origin_truth=None):
    """Time both imputation modes on the same inputs and score masked genotypes."""
    environment = os.environ | dict.fromkeys(('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                              'VECLIB_MAXIMUM_THREADS'), '1')
    environment['NUMBA_NUM_THREADS'] = str(threads)
    common = ['--vcf', str(phased_vcf), '--chr_range', '22', '--pedigree', str(fixture / 'pedigree.txt'),
              '--threads', str(threads), '--no-quality-filter']
    if ibd is not None:
        common += ['--ibd', str(ibd)]
    result = dict(platform=platform.platform(), python=platform.python_version(), threads=threads, modes={})
    modes = {'statistical_phase': ['--no_parent_origin'],
             'parent_origin': ['--parent_origin_map', str(fixture / 'map.txt')]}
    measurements = {mode: [] for mode in modes}
    # Warm up once, then alternate execution order within each measured pair.
    for run in range(repeat + 1):
        for mode in modes if run % 2 else reversed(modes):
            stem = out / f'{mode}_{run}'
            options = [*common, *modes[mode], '--out', str(stem), '-silent_progress']
            measurement = timed_command(step_command('impute', options), stem.with_suffix('.log'), environment)
            if run:
                measurements[mode].append(measurement)
                print(f'{mode} run {run}  {measurement["seconds"]:.3f} s  '
                      f'{measurement["peak_rss_mib"]:.1f} MiB', flush=True)
            else:
                accuracy = check_imputation(fixture, stem.with_suffix('.hdf5'))
                result['modes'][mode] = dict(accuracy=accuracy)
    for mode, times in measurements.items():
        result['modes'][mode].update(seconds=float(np.median([run['seconds'] for run in times])),
                                     peak_rss_mib=float(np.median([run['peak_rss_mib'] for run in times])),
                                     runs=times)
    result['unshared_finite_dosages'] = check_paired_support(out)
    if any(result['unshared_finite_dosages'].values()):
        for mode, paired in zip(modes, reversed(modes)):
            with h5py.File(out / f'{paired}_0.hdf5', 'r') as other:
                result['modes'][mode]['common_support_accuracy'] = check_imputation(
                    fixture, out / f'{mode}_0.hdf5', paired_dosages=other['imputed_par_gts'])
    if origin_truth:
        result['orientation'] = check_origin(fixture, phased_vcf, threads, origin_truth)
    (out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prepare = commands.add_parser('prepare-duos', help='Create a masked real-trio fixture')
    prepare.add_argument('source', type=Path, help='IGSR 30x GRCh38 phased chr22 VCF')
    prepare.add_argument('pedigree', type=Path, help='IGSR pedigree with population labels')
    prepare.add_argument('map', type=Path, help='SHAPEIT GRCh38 chr22 map')
    prepare.add_argument('output', type=Path)
    prepare.add_argument('--start', type=int, default=17_000_000)
    prepare.add_argument('--end', type=int, default=25_000_000)
    phase = commands.add_parser('phase', help='Obtain statistical phase with no pedigree input')
    phase.add_argument('fixture', type=Path)
    phase.add_argument('map', type=Path, help='SHAPEIT GRCh38 chr22 map')
    phase.add_argument('output', type=Path)
    phase.add_argument('--java', default='java')
    phase.add_argument('--beagle-jar', type=Path, required=True)
    phase.add_argument('--threads', type=int, default=4)
    run = commands.add_parser('compare', help='Run both modes and measure accuracy and runtime')
    run.add_argument('fixture', type=Path)
    run.add_argument('vcf', type=Path, help='Statistically phased Beagle VCF')
    run.add_argument('output', type=Path)
    run.add_argument('--ibd', type=Path, help='Sibship IBD file prefix, e.g., ibd_22.ibd')
    run.add_argument('--threads', type=int, default=4)
    run.add_argument('--repeat', type=int, default=3)
    run.add_argument('--origin-truth', choices=('trios', 'simulation'),
                     help='Validate inferred orientation against Mendelian or simulated truth')
    args = parser.parse_args()
    if args.command == 'prepare-duos':
        if not 1 <= args.start < args.end:
            parser.error('Provide a valid chromosome 22 interval')
        prepare_duos(args.source, args.pedigree, args.map, args.output, args.start, args.end)
    else:
        if args.threads < 1 or (args.command == 'compare' and args.repeat < 1):
            parser.error('Threads and repetitions must be positive')
        args.output.mkdir(parents=True, exist_ok=True)
        if args.command == 'phase':
            phase_fixture(args.fixture, args.map, args.output, args.java, args.beagle_jar, args.threads)
        else:
            compare(args.fixture, args.vcf, args.output, args.threads, args.repeat, args.ibd, args.origin_truth)


if __name__ == '__main__':
    main()
