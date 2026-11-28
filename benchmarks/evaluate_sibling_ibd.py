"""Benchmark sibling IBD methods on identical genotype and marker inputs.

Use a chromosome 22 fixture from prepare_ibd_evaluation.py for simulation
truth, or pass --skip-scoring to write inference outputs for a compatible
fixture. Each method directory contains the segments, timings, and available
posterior sidecars and accuracy metrics.
"""

import argparse
from collections import defaultdict, deque
from contextlib import nullcontext
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

from bed_reader import open_bed
import h5py
import numpy as np

from snipar.ibd import _SEGMENT_HEADER, find_segments
from snipar.ibd_output import Relatedness
from snipar.ibd_phased import infer_phased_pair


METHODS = {
    'baseline-pair': (),
    'current-pair': ('--posterior',),
    'current-joint': ('--posterior', '--joint'),
}
PHASED_INPUTS = {
    'phased-true': 'chr_22.vcf.gz',
    'phased-switch': 'switch/chr_22.vcf.gz',
    'phased-statistical': 'beagle/beagle.vcf.gz',
}


def _write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def run_cli(args, method):
    out = args.output / method
    out.mkdir(parents=True, exist_ok=True)
    if not args.rerun and (out / 'timing.json').exists():
        return
    python = args.legacy_python if method == 'baseline-pair' else args.current_python
    launcher = (f"import numpy as np, runpy; np.random.seed({args.seed}); "
                "runpy.run_module('snipar.scripts.ibd', run_name='__main__')")
    command = [str(python), '-c', launcher, '--bed', str(args.fixture / 'chr_22'),
               '--pedigree', str(args.fixture / 'pedigree.txt'), '--map', str(args.fixture / 'map.txt'),
               '--chrom', '22', '--threads', str(args.threads), '--p_error', str(args.p_error),
               '--min_length', '0', '--min_maf', '.01', '--max_missing', '5', '--max_error', '.01',
               '--ld_out', '--out', str(out / 'chr_@'), *METHODS[method]]
    start = perf_counter()
    with (out / 'run.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    elapsed = perf_counter() - start
    _write_json(out / 'timing.json', {'wall_seconds': elapsed, 'command': command})
    print(f'{args.fixture.name} {method}: {elapsed:.2f}s', flush=True)


def _ld_table(path):
    with gzip.open(path, 'rt') as source:
        header = next(source).split()
        if header != ['CHR', 'SNP', 'BP', 'L2']:
            raise ValueError(f'Unexpected LD columns in {path}')
        table = np.loadtxt(source, dtype=str, ndmin=2)
    return table[:, 1], table[:, 2].astype(np.int64), table[:, 3].astype(float)


def common_markers(output):
    file = output / 'current-pair' / 'chr_22.ibd.posterior.hdf5'
    with h5py.File(file) as source:
        pairs = source['pairs'].asstr()[:]
        snps = source['variant_ids'].asstr()[:]
        positions = source['positions'][:]
        cm = source['genetic_map_cm'][:]
    ld = None
    for method in METHODS:
        path = output / method / 'chr_22.l2.ldscore.gz'
        if not path.is_file():
            continue
        ids, bp, values = _ld_table(path)
        np.testing.assert_array_equal(ids, snps, err_msg=method)
        np.testing.assert_array_equal(bp, positions, err_msg=method)
        if ld is None:
            ld = values
        else:
            np.testing.assert_allclose(values, ld, rtol=1e-12, atol=1e-12, err_msg=f'{method} LD')
    if ld is None:
        raise ValueError('Missing matched LD scores')
    return pairs, snps, positions, cm, 1 / np.where(np.isfinite(ld) & (ld > 0), ld, 1)


def _canonical_pair(left, right):
    return (left, right) if left < right else (right, left)


def _segments(path):
    rows = defaultdict(list)
    with gzip.open(path, 'rt') as file:
        for row in csv.DictReader(file, delimiter='\t'):
            rows[_canonical_pair(row['ID1'], row['ID2'])].append((
                int(row['start_coordinate']), int(row['stop_coordinate']),
                row['startSNP'], row['stopSNP'], int(row['IBDType']),
            ))
    return rows


def _hard_from_segments(segments, snps, positions):
    hard = np.full(len(positions), -1, dtype=np.int8)
    for start, stop, first_snp, last_snp, state in segments:
        a, b = np.searchsorted(positions, [start, stop], side='left')
        a_end, b_end = np.searchsorted(positions, [start, stop], side='right')
        possible_start = np.flatnonzero(snps[a:a_end] == first_snp)
        possible_stop = np.flatnonzero(snps[b:b_end] == last_snp)
        if not possible_start.size or not possible_stop.size:
            raise ValueError('IBD segment boundary does not match the common SNP set')
        first, last = a + possible_start[0], b + possible_stop[-1] + 1
        if np.any((hard[first:last] != -1) & (hard[first:last] != state)):
            raise ValueError('Conflicting hard IBD segments overlap')
        hard[first:last] = state
    if np.any(hard < 0):
        raise ValueError('A method did not cover the common marker set')
    return hard


def _read_truth(fixture, cm):
    path = fixture / 'truth' / 'transmissions.tsv'
    if not path.is_file():
        return None, None
    fam = np.loadtxt(fixture / 'chr_22.fam', dtype=str, usecols=1, ndmin=1)
    observed = set(fam)
    transmissions, visible = {}, {}
    with path.open() as file:
        for row in csv.DictReader(file, delimiter='\t'):
            iid = row['IID']
            if iid not in observed:
                continue
            paternal_and_maternal = []
            for parent in ('paternal', 'maternal'):
                raw = row[f'{parent}_crossover_cM']
                crossovers = np.fromstring(raw, sep=',') if raw != '.' else np.empty(0)
                start = int(row[f'{parent}_start'])
                paternal_and_maternal.append((start + np.searchsorted(crossovers, cm, side='right')) % 2)
            transmissions[iid] = np.asarray(paternal_and_maternal, np.int8)
            visible[iid] = int(row['FATHER_ID'] in observed) + int(row['MOTHER_ID'] in observed)
    return transmissions, visible


class Scores:
    def __init__(self):
        self.confusion = np.zeros((3, 3), np.int64)
        self.posterior_confusion = np.zeros((3, 3), np.int64)
        self.brier = 0.
        self.n_coherent = self.n_family_markers = self.n_pairs = 0
        self.calibration = np.zeros((3, 10, 3))

    def update(self, truth, path, posterior):
        self.n_pairs += 1
        self.confusion += np.bincount(3 * truth + path, minlength=9).reshape(3, 3)
        if posterior is None:
            return
        if posterior.shape != (len(truth), 3) or not np.isfinite(posterior).all():
            raise ValueError('Malformed IBD probabilities')
        np.testing.assert_allclose(posterior.sum(axis=1), 1., atol=1e-6)
        posterior_map = np.argmax(posterior, axis=1)
        self.posterior_confusion += np.bincount(3 * truth + posterior_map, minlength=9).reshape(3, 3)
        for state in range(3):
            observed = truth == state
            p = posterior[:, state]
            bins = np.minimum((p * 10).astype(np.int8), 9)
            self.calibration[state, :, 0] += np.bincount(bins, minlength=10)
            self.calibration[state, :, 1] += np.bincount(bins, weights=p, minlength=10)
            self.calibration[state, :, 2] += np.bincount(bins, weights=observed, minlength=10)
            self.brier += float(np.square(p - observed).sum(dtype=np.float64))

    def coherence(self, pairs, paths):
        names, rows = np.unique(pairs, return_inverse=True)
        n = len(names)
        if n < 3:
            return
        if n > 5 or len(pairs) != n * (n - 1) // 2:
            raise ValueError('Coherence metric requires complete sibling-family pairs of size up to five')
        relative = ((np.arange(2 ** (n - 1))[:, None] >> np.arange(n - 1)) & 1)
        inheritance = np.column_stack((np.zeros(len(relative), dtype=int), relative))
        rows = rows.reshape(-1, 2)
        share = (inheritance[:, rows[:, 0]] == inheritance[:, rows[:, 1]]).astype(np.int64)
        radix = 3 ** np.arange(len(pairs))
        known = (share @ radix)
        possible = np.unique((known[:, None] + known[None, :]).ravel())
        encoded_paths = radix @ paths
        self.n_coherent += int(np.isin(encoded_paths, possible).sum())
        self.n_family_markers += paths.shape[1]

    @staticmethod
    def _accuracy(confusion):
        n = confusion.sum()
        if not n:
            return {}
        rows = confusion.sum(axis=1)
        accuracy = np.divide(np.diag(confusion), rows, out=np.zeros(3), where=rows > 0)
        return dict(accuracy=float(np.trace(confusion) / n),
                    balanced_accuracy=float(accuracy[rows > 0].mean()),
                    accuracy_by_true_state=[float(a) if count else None for a, count in zip(accuracy, rows)],
                    confusion=confusion.tolist())

    def summary(self):
        result = dict(n_pairs=self.n_pairs, pair_markers=int(self.confusion.sum()),
                      hard=self._accuracy(self.confusion))
        if self.n_family_markers:
            result.update(n_family_markers=self.n_family_markers,
                          coherent_fraction=self.n_coherent / self.n_family_markers)
        n = self.posterior_confusion.sum()
        if n:
            ec = np.abs(self.calibration[:, :, 1] - self.calibration[:, :, 2]).sum(axis=1) / n
            result['posterior'] = dict(**self._accuracy(self.posterior_confusion),
                                       brier_score=float(self.brier / n),
                                       calibration_error_by_state=ec.tolist(),
                                       calibration_bins=[
                                           [{'count': int(count), 'expected': float(predicted / count),
                                             'observed': float(observed / count), 'lower': i / 10}
                                            for i, (count, predicted, observed) in enumerate(table) if count]
                                           for table in self.calibration
                                       ])
        return result


class Evaluation:
    def __init__(self, fixture, cm, pairs, *, skip_scoring=False):
        self.transmissions, self.visible = (None, None) if skip_scoring else _read_truth(fixture, cm)
        self.groups = defaultdict(list)
        self.records = {k: Scores() for k in ('all', '0', '1', '2')}
        ped = np.loadtxt(fixture / 'pedigree.txt', dtype=str, skiprows=1, ndmin=2)
        label = dict(zip(ped[:, 1], ped[:, 0]))
        for index, (left, right) in enumerate(pairs):
            if label[left] != label[right]:
                raise ValueError('Sibling pairs must belong to the same family')
            self.groups[label[left]].append(index)

    def pair(self, left, right, path, posterior):
        if self.transmissions is None:
            return
        truth = np.count_nonzero(self.transmissions[left] == self.transmissions[right], axis=0).astype(np.int8)
        for key in ('all', str(self.visible[left])):
            self.records[key].update(truth, path, posterior)

    def family(self, pairs, paths):
        if self.transmissions is None:
            return
        for key in ('all', str(self.visible[pairs[0, 0]])):
            self.records[key].coherence(pairs, paths)

    def summary(self, path, **options):
        if self.transmissions is not None:
            result = {f'{name}_visible_parents' if name != 'all' else name: score.summary()
                      for name, score in self.records.items() if score.n_pairs}
            _write_json(path, dict(**options, groups=result))


def evaluate_cli(args, method, markers):
    dest = args.output / method
    pairs, snps, positions, cm, _ = markers
    evaluation = Evaluation(args.fixture, cm, pairs)
    if evaluation.transmissions is None:
        return
    segments = _segments(dest / 'chr_22.ibd.segments.gz')
    if set(segments) != {_canonical_pair(*p) for p in pairs}:
        raise ValueError(f'{method} pair coverage is mismatched')
    posterior_path = dest / 'chr_22.ibd.posterior.hdf5'
    context = h5py.File(posterior_path, 'r') if posterior_path.exists() else nullcontext()
    with context as sidecar:
        if sidecar is not None:
            np.testing.assert_array_equal(sidecar['positions'][:], positions)
            np.testing.assert_array_equal(sidecar['variant_ids'].asstr()[:], snps)
            np.testing.assert_array_equal(sidecar['genetic_map_cm'][:], cm)
            posterior_index = {_canonical_pair(*p): i for i, p in enumerate(sidecar['pairs'].asstr()[:])}
            if set(posterior_index) != set(segments):
                raise ValueError('Posterior and hard IBD cover different pairs')
        for indices in evaluation.groups.values():
            family_pairs = pairs[indices]
            paths = []
            for left, right in family_pairs:
                key = _canonical_pair(left, right)
                hard = _hard_from_segments(segments[key], snps, positions)
                posterior = None if sidecar is None else sidecar['probabilities'][posterior_index[key]]
                evaluation.pair(left, right, hard, posterior)
                paths.append(hard)
            evaluation.family(family_pairs, np.asarray(paths))
    settings = json.loads((dest / 'timing.json').read_text())
    evaluation.summary(dest / 'metrics.json', method=method, **settings)


def verify_joint_fallback(output):
    pair, joint = output / 'current-pair', output / 'current-joint'
    segments = 'chr_22.ibd.segments.gz'
    with gzip.open(pair / segments, 'rb') as a, gzip.open(joint / segments, 'rb') as b:
        if a.read() != b.read():
            raise ValueError('Joint fallback changed the pairwise hard IBD paths')
    posterior = 'chr_22.ibd.posterior.hdf5'
    with h5py.File(pair / posterior) as a, h5py.File(joint / posterior) as b:
        for name in ('pairs', 'variant_ids', 'positions', 'genetic_map_cm'):
            np.testing.assert_array_equal(a[name][:], b[name][:])
        np.testing.assert_equal(a['probabilities'].shape, b['probabilities'].shape)
        for left, right in zip(a['probabilities'], b['probabilities']):
            np.testing.assert_array_equal(left, right)
    _write_json(joint / 'fallback.json', dict(identical_hard_segments=True, identical_posterior=True))


def _source_columns(fixture, snps, positions):
    bim = np.loadtxt(fixture / 'chr_22.bim', dtype=str, ndmin=2)
    occurrence = defaultdict(deque)
    for i, row in enumerate(bim):
        occurrence[int(row[3]), row[1]].append(i)
    indices = np.array([occurrence[pos, sid].popleft() for pos, sid in zip(positions, snps)])
    return bim[indices], indices


def _observed(fixture, pairs, indices):
    with open_bed(fixture / 'chr_22.bed', count_A1=True) as bed:
        order = dict(zip(bed.iid, range(len(bed.iid))))
        samples = sorted(set(pairs.ravel()), key=order.__getitem__)
        calls = bed.read((np.array([order[iid] for iid in samples]), indices), dtype='int8')
    mask = calls >= 0
    freqs = np.sum(calls, axis=0, where=mask) / (2 * mask.sum(axis=0))
    return samples, calls, freqs


def _phased_genotypes(file, bim, samples, bed_calls):
    ids = defaultdict(list)
    for i, row in enumerate(bim):
        _, sid, _, position, a1, a2 = row
        ids[(int(position), *sorted((a1, a2)))].append((i, sid, a1))
    haps = np.empty((*bed_calls.shape, 2), dtype=np.int8)
    filled = np.zeros(len(bim), dtype=bool)
    with gzip.open(file, 'rb') as lines:
        for line in lines:
            if line.startswith(b'#CHROM\t'):
                header = line.split()[9:]
                rows = dict(zip(header, range(len(header))))
                select = np.array([rows[iid.encode()] for iid in samples])
                break
        else:
            raise ValueError(f'No sample header in {file}')
        fast_offsets = select[:, None] * 4 + np.array([0, 2])
        for raw in lines:
            fields = raw.rstrip(b'\r\n').split(b'\t', 9)
            position = int(fields[1])
            sid, ref, alt = [value.decode() for value in fields[2:5]]
            candidates = [v for v in ids.get((position, *sorted((ref, alt))), ()) if not filled[v[0]]]
            if not candidates:
                continue
            if len(candidates) > 1:
                candidates = [v for v in candidates if v[1] == sid]
                if len(candidates) != 1:
                    raise ValueError(f'Ambiguous genotype match at {position} in {file}')
            (target, _, a1), = candidates
            values = fields[9]
            if fields[8] == b'GT' and len(values) == 4 * len(header) - 1:
                data = np.frombuffer(values, dtype=np.int8)
                alleles = data[fast_offsets] - 48
                ordered = data[fast_offsets[:, 0] + 1] == ord('|')
            else:
                cells = values.split(b'\t')
                column = fields[8].split(b':').index(b'GT')
                calls = np.array([list(cells[index].split(b':')[column]) for index in select], dtype=np.int8)
                if calls.shape != (len(samples), 3):
                    raise ValueError(f'Invalid diploid genotype in {file} at {position}')
                alleles = calls[:, ::2] - 48
                ordered = calls[:, 1] == ord('|')
            supplied = bed_calls[:, target] >= 0
            if np.any(supplied & ~ordered & (alleles[:, 0] != alleles[:, 1])):
                raise ValueError(f'Unphased heterozygous genotype in {file} at {position}')
            if np.any(supplied[:, None] & ((alleles < 0) | (alleles > 1))):
                raise ValueError(f'Missing or invalid phased allele in {file} at {position}')
            if a1 == ref:
                alleles = 1 - alleles
            np.testing.assert_array_equal(alleles[supplied].sum(axis=1), bed_calls[supplied, target],
                                          err_msg=f'{file} at {position} differs from BED dosages')
            alleles[~supplied] = -1
            haps[:, target] = alleles
            filled[target] = True
    if not filled.all():
        raise ValueError(f'{file} is missing {np.count_nonzero(~filled)} markers in the common set')
    return haps


def run_phased(args, method, source, markers, observed_data):
    dest = args.output / method
    setting = dict(phase_ld_power=args.phase_ld_power, allele_error=args.p_error / 2, input=str(source))
    if method == 'phased-switch':
        setting['switch_rate_per_original_interval'] = args.corruption_rate
    elif method == 'phased-statistical':
        setting['switch_rate_per_retained_heterozygote'] = args.statistical_switch_rate
    else:
        setting['switch_rate_per_interval'] = 0.
    dest.mkdir(parents=True, exist_ok=True)
    pairs, snps, positions, cm, weights = markers
    evaluation = Evaluation(args.fixture, cm, pairs, skip_scoring=args.skip_scoring)
    write_posterior = args.write_phased_posterior or max(map(len, evaluation.groups.values())) == 1
    existing = dest / 'timing.json'
    if (not args.rerun and existing.exists() and (dest / 'chr_22.ibd.segments.gz').is_file() and
            (args.skip_scoring or (dest / 'metrics.json').is_file()) and
            (not write_posterior or (dest / 'chr_22.ibd.posterior.hdf5').is_file())):
        saved = {'phase_ld_power': 1.} | json.loads(existing.read_text())
        if all(saved.get(name) == value for name, value in setting.items()):
            return
    weights = weights ** args.phase_ld_power
    bim, indices, samples, bed_calls, freqs = observed_data
    start = perf_counter()
    haps = _phased_genotypes(source, bim, samples, bed_calls)
    loading = perf_counter() - start
    lookup = {iid: i for i, iid in enumerate(samples)}
    if method == 'phased-switch':
        # Compose independent switch parities over skipped original markers.
        switches = (1 - (1 - 2 * args.corruption_rate) ** np.diff(indices)) * .5
    elif method == 'phased-statistical':
        h = haps[:, 1:]
        switches = args.statistical_switch_rate * ((h[:, :, 0] >= 0) & (h[:, :, 1] >= 0) &
                                                    (h[:, :, 0] != h[:, :, 1]))
    else:
        switches = 0.
    prefix = str(dest / 'chr_@')
    relatedness = Relatedness(pairs)
    start = perf_counter()
    scoring_time = 0.
    target = dest / 'chr_22.ibd.segments.gz'
    context = relatedness.chromosome_batch(np.arange(len(pairs)), prefix) if write_posterior else nullcontext()
    with TemporaryDirectory(dir=dest, prefix='.segments-') as stage:
        pending = Path(stage) / target.name
        with context as store, gzip.open(pending, 'wt', compresslevel=3) as segments:
            segments.write(_SEGMENT_HEADER)
            for family_rows in evaluation.groups.values():
                family_pairs = pairs[family_rows]
                paths = []
                for row in family_rows:
                    left, right = pairs[row]
                    call_rows = [lookup[left], lookup[right]]
                    phase_rate = switches[call_rows] if np.ndim(switches) == 2 else switches
                    posterior = infer_phased_pair(haps[call_rows], freqs, cm, phase_rate,
                                                  args.p_error / 2, weights=weights)
                    path = posterior.hard_ibd
                    for segment in find_segments(path, cm, snps, positions):
                        segments.write(segment.to_text(left, right, 22))
                    if store is not None:
                        store(22, pairs[row:row+1], snps, positions, cm,
                              posterior.ibd_posterior[None, :], np.array([row]))
                    score_start = perf_counter()
                    evaluation.pair(left, right, path, posterior.ibd_posterior)
                    scoring_time += perf_counter() - score_start
                    paths.append(path)
                score_start = perf_counter()
                evaluation.family(family_pairs, np.asarray(paths))
                scoring_time += perf_counter() - score_start
        pending.replace(target)
    if write_posterior:
        relatedness.write(prefix)
    elapsed = perf_counter() - start - scoring_time
    timing = dict(wall_seconds=elapsed, vcf_loading_seconds=loading, **setting)
    _write_json(dest / 'timing.json', timing)
    evaluation.summary(dest / 'metrics.json', method=method, **timing)
    print(f'{args.fixture.name} {method}: {elapsed:.2f}s plus {loading:.2f}s VCF loading', flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--legacy-python', type=Path, help='Executable loading the baseline SNIPar package')
    parser.add_argument('--current-python', type=Path, default=Path(sys.executable),
                        help='Executable loading the current SNIPar package')
    parser.add_argument('--stage', choices=('all', 'unphased', 'phased'), default='all')
    parser.add_argument('--phase-arms', nargs='+', choices=PHASED_INPUTS, default=list(PHASED_INPUTS))
    parser.add_argument('--seed', type=int, default=20261128)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--p-error', type=float, default=.001)
    parser.add_argument('--corruption-rate', type=float, default=.001)
    parser.add_argument('--statistical-switch-rate', type=float, default=.002)
    parser.add_argument('--phase-ld-power', type=float, default=1.,
                        help='Power applied to the inverse LD weights for phased emissions')
    parser.add_argument('--write-phased-posterior', action='store_true', help='Write phased sidecars in larger sibships too')
    parser.add_argument('--verify-joint-fallback', action='store_true',
                        help='Also run --joint on unsupported sibships and compare to current pairwise outputs')
    parser.add_argument('--skip-scoring', action='store_true', help='Create inference outputs without simulation truth')
    parser.add_argument('--rerun', action='store_true')
    args = parser.parse_args(argv)
    if not np.isfinite(args.phase_ld_power) or args.phase_ld_power < 0:
        parser.error('--phase-ld-power must be finite and nonnegative')
    args.fixture, args.output = args.fixture.resolve(), args.output.resolve()
    if args.stage != 'phased':
        if args.legacy_python is None:
            parser.error('--legacy-python is required for unphased comparisons')
        run_cli(args, 'baseline-pair')
        run_cli(args, 'current-pair')
        pair_file = args.output / 'current-pair' / 'chr_22.ibd.posterior.hdf5'
        from snipar.ibd_joint import joint_pair_groups
        with h5py.File(pair_file) as source:
            has_joint = bool(joint_pair_groups(source['pairs'].asstr()[:]))
        if has_joint or args.verify_joint_fallback:
            run_cli(args, 'current-joint')
        if not has_joint and args.verify_joint_fallback:
            verify_joint_fallback(args.output)
    markers = common_markers(args.output)
    if args.stage != 'phased' and not args.skip_scoring:
        for method in METHODS:
            if (args.output / method / 'timing.json').is_file():
                evaluate_cli(args, method, markers)
    if args.stage != 'unphased':
        pairs, snps, positions, _, _ = markers
        bim, indices = _source_columns(args.fixture, snps, positions)
        observed = (bim, indices, *_observed(args.fixture, pairs, indices))
        for method in args.phase_arms:
            source = args.fixture / PHASED_INPUTS[method]
            if source.is_file():
                log = source.with_suffix('').with_suffix('.log')
                if method == 'phased-statistical' and log.is_file() and 'End time:' not in log.read_text():
                    print(f'Waiting for Beagle to finish before decoding {source}', flush=True)
                    continue
                run_phased(args, method, source, markers, observed)


if __name__ == '__main__':
    main()
