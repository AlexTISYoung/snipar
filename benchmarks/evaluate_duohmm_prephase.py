"""Benchmark duo-guided IBD and parental imputation on unchanged Beagle calls.

The previous evaluation supplies Beagle genotypes, allele-aligned parental
truth, and dosage/joint/standalone phased references. The output directory
holds new inference products and paired scores. Run ``--stage timing``
after inference for the alternating warm, four-block k5 reuse comparison.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from dataclasses import dataclass
import json
import os
from pathlib import Path
from statistics import median
from types import SimpleNamespace

import h5py
import numpy as np

from benchmark_ibd_imputation import METRICS, _family_score, _mask_results
from evaluate_empirical_ibd_quads import score_ibd as score_empirical
from evaluate_sibling_ibd import _canonical_pair, _segments, common_markers, evaluate_cli
from profile_1000g_eur import timed_command


NEW_IBD = ('raw-phased', 'duo-phased')
OLD_IBD = ('baseline-pair', 'current-pair', 'current-joint', 'phased-statistical')


@dataclass(frozen=True)
class Cohort:
    name: str
    size: int
    fixture: Path
    old_ibd: Path
    old_imputation: Path
    bgen: Path
    output: Path
    empirical: bool

    @property
    def ibd(self):
        return self.output / 'ibd'

    @property
    def imputation(self):
        return self.output / 'imputation'


def cohort_inputs(previous, output, sizes, empirical):
    for size in sizes:
        name = f'k{size}'
        yield Cohort(name, size, previous / 'data' / name, previous / 'ibd' / name,
                     previous / 'imputation' / name / 'statistical',
                     previous / 'imputation' / 'genotypes' / name / 'statistical' / 'genotypes',
                     output / name, False)
    if empirical:
        for mask in ('one_parent', 'zero_parents'):
            name = f'empirical-{mask}'
            prior = previous / 'empirical-results' / mask
            yield Cohort(name, 2, previous / 'empirical' / mask, prior / 'ibd', prior / 'beagle',
                         output / name / 'genotypes' / 'genotypes', output / name, True)


def link_input(source, target):
    source = source.resolve(strict=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_symlink() and target.resolve() == source:
        return
    target.symlink_to(source, target_is_directory=source.is_dir())


def prepare_cohort(cohort):
    for method in OLD_IBD:
        path = cohort.old_ibd / method
        if path.is_dir():
            link_input(path, cohort.ibd / method)
    for method in ('pair_hard', 'pair_hard_origin', 'pair_posterior', 'joint_hard', 'joint_fallback_hard',
                   'phased_hard', 'phased_posterior'):
        suffix = '_22.hdf5' if cohort.empirical else '.hdf5'
        source = cohort.old_imputation / (method + suffix)
        if source.is_file():
            link_input(source, cohort.imputation / f'{method}.hdf5')
            timing = source.with_suffix('.run.json')
            if timing.exists():
                link_input(timing, cohort.imputation / f'{method}.run.json')

    bgen = cohort.bgen.with_suffix('.bgen')
    if cohort.empirical and not bgen.is_file():
        from snipar.vcf import convert_vcf

        source = cohort.fixture / 'beagle' / 'beagle.vcf.gz'
        result = convert_vcf(source, cohort.bgen.parent, chromosomes=[22], output='phased', min_quality={})
        if Path(result.bgen) != cohort.bgen:
            raise ValueError(f'Unexpected phased genotype prefix: {result.bgen}')
    if not bgen.is_file():
        raise FileNotFoundError(bgen)
    details = dict(fixture=str(cohort.fixture), genotype=str(bgen), previous_ibd=str(cohort.old_ibd),
                   previous_imputation=str(cohort.old_imputation), switch_rate_per_retained_heterozygote=.002,
                   genotyping_error=.001, phase_allele_error=.0005)
    (cohort.output / 'inputs.json').write_text(json.dumps(details, indent=2) + '\n')


def run_ibd(cohort, method, python, threads, environment):
    out = cohort.ibd / method
    out.mkdir(exist_ok=True)
    stem = out / 'chr_22.ibd'
    posterior, origin = (Path(str(stem) + suffix) for suffix in ('.posterior.hdf5', '.origin.hdf5'))
    timing = out / 'timing.json'
    needs_origin = method == 'duo-phased' and cohort.name != 'empirical-zero_parents'
    if timing.is_file() and posterior.is_file() and (not needs_origin or origin.is_file()):
        return
    seed = 20261127 if cohort.empirical else 20261128
    launcher = (f'import numpy as np, runpy; np.random.seed({seed}); '
                'runpy.run_module("snipar.scripts.ibd", run_name="__main__")')
    command = [str(python), '-c', launcher, '--bgen', str(cohort.bgen), '--phased-ibd', '--posterior',
               '--pedigree', str(cohort.fixture / 'pedigree.txt'), '--map', str(cohort.fixture / 'map.txt'),
               '--chrom', '22', '--threads', str(threads), '--min_maf', '.01', '--max_missing', '5',
               '--p_error', '.001', '--max_error', '.01', '--min_length', '0', '--phase-switch-rate', '.002',
               '--ld_out', '--out', str(out / 'chr_@')]
    if method == 'raw-phased':
        command.append('--no-duohmm')
    measurement = timed_command(command, out / 'run.log', environment)
    if not posterior.is_file():
        raise FileNotFoundError(posterior)
    if needs_origin and not origin.is_file():
        raise FileNotFoundError(origin)
    timing.write_text(json.dumps({'wall_seconds': measurement.pop('seconds'), **measurement, 'command': command},
                                 indent=2) + '\n')
    print(cohort.name, method, 'IBD complete', flush=True)


def imputation_arms(cohort, stage):
    if stage == 'reference':
        if not (cohort.imputation / 'phased_hard.hdf5').exists():
            yield 'phased_hard', 'phased-statistical', False, 'recompute'
        return
    yield 'raw_cli_hard', 'raw-phased', False, 'recompute'
    yield 'duo_hard', 'duo-phased', False, 'sidecar'
    yield 'duo_hard_recomputed', 'duo-phased', False, 'recompute'
    if cohort.size == 2:
        yield 'raw_cli_posterior', 'raw-phased', True, 'recompute'
        yield 'duo_posterior', 'duo-phased', True, 'sidecar'
        yield 'duo_posterior_recomputed', 'duo-phased', True, 'recompute'


def run_imputation(cohort, arm, python, threads, environment):
    name, method, soft, mode = arm
    target = cohort.imputation / name
    output = target.with_suffix('.hdf5')
    if output.exists():
        return
    source = cohort.ibd / method / 'chr_22.ibd'
    suffixes = ('.segments.gz', '.posterior.hdf5', '.origin.hdf5')
    segments, posterior, origin = (Path(f'{source}{suffix}') for suffix in suffixes)
    # Keep both automatic sidecars out of the prefix, then explicitly select each.
    isolated = cohort.output / 'hard-input' / method / 'chr_22.ibd'
    link_input(segments, Path(f'{isolated}.segments.gz'))
    command = [str(python), '-m', 'snipar.scripts.impute', '--bgen', str(cohort.bgen), '--ibd', str(isolated),
               '--pedigree', str(cohort.fixture / 'pedigree.txt'), '--out', str(target), '--threads', str(threads),
               '--parent_origin_map', str(cohort.fixture / 'map.txt'), '-silent_progress']
    if soft:
        command += ['--ibd-posterior', str(posterior.resolve(strict=True))]
    if mode == 'sidecar' and origin.is_file():
        command += ['--parent-origin', str(origin)]
    measurement = timed_command(command, target.with_suffix('.log'), environment)
    if not output.is_file():
        raise FileNotFoundError(output)
    target.with_suffix('.run.json').write_text(json.dumps({**measurement, 'command': command}, indent=2) + '\n')
    print(cohort.name, name, 'imputation complete', flush=True)


def score_dosages(cohort, bootstrap, seed):
    paths = {path.stem: path for path in sorted(cohort.imputation.glob('*.hdf5'))}
    scores = _family_score(paths, cohort.fixture)
    references = ('pair_hard_origin', 'pair_posterior', 'joint_hard', 'phased_hard', 'phased_posterior',
                  'raw_cli_hard', 'raw_cli_posterior', 'duo_hard_recomputed', 'duo_posterior_recomputed')
    results = {mask: _mask_results(mask, *moments, bootstrap, seed, references=references)
               for mask, moments in scores.items()}
    for methods in results.values():
        for name, result in methods.items():
            timing = paths[name].with_suffix('.run.json')
            if timing.is_file():
                result['timing'] = json.loads(timing.read_text())
    return results


def score_sharing(cohort):
    methods = tuple(method for method in (*OLD_IBD, *NEW_IBD) if (cohort.ibd / method).is_dir())
    if cohort.empirical:
        return score_empirical(cohort.fixture.parent, cohort.fixture, cohort.ibd, methods=methods)
    markers = common_markers(cohort.ibd)
    for method in NEW_IBD:
        if method in methods:
            args = SimpleNamespace(fixture=cohort.fixture, output=cohort.ibd)
            evaluate_cli(args, method, markers)
    return {method: json.loads((cohort.ibd / method / 'metrics.json').read_text())
            for method in methods if (cohort.ibd / method / 'metrics.json').is_file()}


def verify_ineligible_pairs(cohort):
    """The phased comparator and duo arm must agree without a single observed parent."""
    paths = [cohort.ibd / method / 'chr_22.ibd.posterior.hdf5' for method in NEW_IBD]
    if not all(path.is_file() for path in paths):
        return None
    visible = set(np.loadtxt(cohort.fixture / 'chr_22.fam', usecols=1, dtype=str, ndmin=1))
    ped = np.loadtxt(cohort.fixture / 'pedigree.txt', dtype=str, skiprows=1, ndmin=2)
    parents = {child: int(father in visible) + int(mother in visible) for _, child, father, mother in ped}
    segment_tables = [_segments(path.parent / 'chr_22.ibd.segments.gz') for path in paths]
    tested = dict(zero_parent_pairs=0, two_parent_pairs=0, one_parent_pairs_changed=0)
    with h5py.File(paths[0]) as raw, h5py.File(paths[1]) as duo:
        for key in ('pairs', 'positions', 'variant_ids', 'genetic_map_cm'):
            np.testing.assert_array_equal(raw[key][:], duo[key][:])
        for index, (left, right) in enumerate(raw['pairs'].asstr()[:]):
            pair = _canonical_pair(left, right)
            parent_count = parents[left]
            hard_equal = segment_tables[0][pair] == segment_tables[1][pair]
            probabilities_equal = np.array_equal(raw['probabilities'][index], duo['probabilities'][index])
            if parent_count != 1:
                if not hard_equal or not probabilities_equal:
                    raise ValueError(f'Duo changed an IBD pair with {parent_count} observed parents: {pair}')
                tested['zero_parent_pairs' if parent_count == 0 else 'two_parent_pairs'] += 1
            elif not (hard_equal and probabilities_equal):
                tested['one_parent_pairs_changed'] += 1
    return tested


def compare_origin_reuse(cohort):
    """Compare both dosage support and estimates on the same IBD predictions."""
    result = {}
    for method in ('duo_hard', 'duo_posterior'):
        reused = cohort.imputation / f'{method}.hdf5'
        recomputed = cohort.imputation / f'{method}_recomputed.hdf5'
        if not reused.is_file() or not recomputed.is_file():
            continue
        with h5py.File(reused) as a, h5py.File(recomputed) as b:
            for name in ('families', 'pos', 'bim_values'):
                np.testing.assert_array_equal(a[name][:], b[name][:])
            for channel in ('imputed_par_gts', 'robust_imputed_par_gts', 'num_observed_parental_alleles',
                            'num_expected_observed_parental_alleles'):
                if channel not in a:
                    continue
                left, right = a[channel], b[channel]
                assert left.shape == right.shape
                same_support, equal, evaluated = 0, 0, 0
                maximum = 0.
                for start in range(0, left.shape[1], 2048):
                    block = slice(start, start + 2048)
                    x, y = left[:, block], right[:, block]
                    finite = np.isfinite(x) & np.isfinite(y)
                    same_support += np.count_nonzero(np.isfinite(x) == np.isfinite(y))
                    equal += np.count_nonzero(finite & (x == y))
                    evaluated += np.count_nonzero(finite)
                    maximum = max(maximum, np.max(np.abs(x[finite] - y[finite]), initial=0.))
                result.setdefault(method, {})[channel] = dict(cells=int(np.prod(left.shape)),
                    common_finite=int(evaluated), equal_common=int(equal), same_support=int(same_support),
                    max_absolute_difference=float(maximum))
    origin_path = cohort.ibd / 'duo-phased' / 'chr_22.ibd.origin.hdf5'
    if origin_path.is_file():
        with h5py.File(origin_path) as sidecar:
            labels = sidecar['origin']
            calls = copies = 0
            for start in range(0, labels.shape[1], 4096):
                stop = start + 4096
                calls += np.count_nonzero(labels[:, start:stop] >= 0)
                copies += np.count_nonzero(sidecar['parent_copy'][:, start:stop] >= 0)
            result['sidecar'] = dict(duos=int(labels.shape[0]), markers=int(labels.shape[1]),
                                     origin_call_fraction=float(calls / labels.size),
                                     parent_copy_call_fraction=float(copies / labels.size))
    return result


def score(cohorts, out, bootstrap, seed):
    results = {}
    for cohort in cohorts:
        results[cohort.name] = dict(ibd=score_sharing(cohort), imputation=score_dosages(cohort, bootstrap, seed),
                                    paired_controls=verify_ineligible_pairs(cohort), origin_reuse=compare_origin_reuse(cohort))
        print(cohort.name, 'scored', flush=True)
    (out / 'metrics.json').write_text(json.dumps(results, indent=2) + '\n')
    with (out / 'imputation.tsv').open('w') as stream:
        names = ('cohort', 'mask', 'method', 'families', 'n_truth', 'n_finite', 'n_common', 'coverage', *METRICS)
        writer = csv.DictWriter(stream, names, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for name, result in results.items():
            for mask, methods in result['imputation'].items():
                for method, entry in methods.items():
                    writer.writerow(dict(cohort=name, method=method, **entry) | entry['estimates'])


def time_origin_reuse(cohort, environment):
    destination = cohort.output / 'origin-timing-chunks4'
    destination.mkdir(exist_ok=True)
    names = {'reused': 'duo_hard', 'recomputed': 'duo_hard_recomputed'}
    default = {method: json.loads((cohort.imputation / f'{name}.run.json').read_text())
               for method, name in names.items()}
    measurements = {method: [] for method in names}
    for run in range(4):
        order = tuple(names) if run % 2 == 0 else tuple(reversed(names))
        outputs = {}
        for method in order:
            target = destination / f'{method}-{run}'
            command = default[method]['command'].copy()
            command[command.index('--out') + 1] = str(target)
            command += ['--chunks', '4']
            result = timed_command(command, target.with_suffix('.log'), environment)
            measurements[method].append(dict(run=run, warmup=run == 0, **result, command=command))
            outputs[method] = target.with_suffix('.hdf5')
            print(method, run, result['seconds'], flush=True)
        with h5py.File(outputs['reused']) as reused, h5py.File(outputs['recomputed']) as recomputed:
            for key in ('imputed_par_gts', 'num_observed_parental_alleles'):
                dataset = reused[key]
                for start in range(0, dataset.shape[1], 2048):
                    block = slice(start, start + 2048)
                    np.testing.assert_array_equal(dataset[:, block], recomputed[key][:, block])
    elapsed = {method: median(row['seconds'] for row in rows if not row['warmup'])
               for method, rows in measurements.items()}
    sampled_rss = {method: median(row['peak_rss_mib'] for row in rows if not row['warmup'])
                   for method, rows in measurements.items()}
    result = dict(chunks=4, median_seconds=elapsed, median_peak_rss_mib=sampled_rss,
                  speedup=elapsed['recomputed'] / elapsed['reused'], exact_output_parity=True,
                  default_chunk_measurements=default, runs=measurements)
    (destination / 'timing.json').write_text(json.dumps(result, indent=2) + '\n')
    print('warm median seconds', elapsed, 'peak MiB', sampled_rss, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('previous', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--python', type=Path, required=True)
    parser.add_argument('--stage', choices=('prepare', 'reference', 'duo', 'score', 'timing', 'all'), default='all')
    parser.add_argument('--sizes', nargs='+', type=int, choices=(2, 3, 4, 5), default=(2, 3, 4, 5))
    parser.add_argument('--empirical', action='store_true')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--bootstrap', type=int, default=999)
    parser.add_argument('--seed', type=int, default=94821)
    args = parser.parse_args()
    if min(args.threads, args.jobs, args.bootstrap) < 1:
        parser.error('threads, jobs and bootstrap count must be positive')
    out = args.output.resolve()
    previous = args.previous.resolve()
    if out == previous:
        parser.error('Use a separate output directory from the prior benchmark')
    cohorts = list(cohort_inputs(previous, out, args.sizes, args.empirical))
    environment = os.environ | dict.fromkeys(('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                             'VECLIB_MAXIMUM_THREADS'), '1')
    environment['NUMBA_NUM_THREADS'] = str(args.threads)
    for cohort in cohorts:
        prepare_cohort(cohort)
    if args.stage in ('reference', 'all'):
        for cohort in cohorts:
            for arm in imputation_arms(cohort, 'reference'):
                run_imputation(cohort, arm, args.python, args.threads, environment)
    if args.stage in ('duo', 'all'):
        def evaluate(cohort):
            for method in NEW_IBD:
                run_ibd(cohort, method, args.python, args.threads, environment)
            for arm in imputation_arms(cohort, 'duo'):
                run_imputation(cohort, arm, args.python, args.threads, environment)

        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            list(pool.map(evaluate, cohorts))
    if args.stage in ('score', 'all'):
        score(cohorts, out, args.bootstrap, args.seed)
    if args.stage == 'timing':
        k5 = next((cohort for cohort in cohorts if cohort.size == 5), None)
        if k5 is None:
            parser.error('Include size 5 to time origin reuse')
        time_origin_reuse(k5, environment)


if __name__ == '__main__':
    main()
