"""Run matched hard and posterior imputation on repeatedly masked real sibling quads."""

import argparse
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

from bed_reader import open_bed
import h5py
import numpy as np

from benchmark_ibd_imputation import _family_score, _mask_results
from calibrate_ibd_posteriors import calibration
from evaluate_sibling_ibd import _canonical_pair, _hard_from_segments, _segments, common_markers, run_cli
from evaluate_sibling_ibd import main as run_benchmark_ibd
from snipar.ibd import _chromosome_hmm, _infer_block, _posterior_tables
from snipar.ld import compute_ld_scores_packed, pack_ld_hardcalls


def reference_sensitivity(fixture, pairs, markers, cm, weights, truth, valid, baseline):
    """Separate frequency estimation and LD weighting in the small real cohort."""
    with open_bed(fixture / 'chr_22.bed', count_A1=True) as bed:
        observed = {iid: i for i, iid in enumerate(bed.iid)}
        variants = {sid: i for i, sid in enumerate(bed.sid)}
        columns = np.array([variants[sid] for sid in markers])
        sibling_rows = np.array(sorted({observed[iid] for iid in pairs.flat}))
        row_index = {row: i for i, row in enumerate(sibling_rows)}
        pair_rows = np.array([[row_index[observed[iid]] for iid in pair] for pair in pairs])
        reference_rows = np.flatnonzero(np.char.startswith(np.asarray(bed.fid, dtype=str), 'freq'))
        siblings = bed.read((sibling_rows, columns), dtype='int8')
        reference = bed.read((reference_rows, columns), dtype='int8')
    if np.any(siblings < 0) or np.any(reference < 0):
        raise ValueError('The empirical fixture must have called sibling and reference genotypes')
    child_frequency, reference_frequency = (data.mean(axis=0) / 2 for data in (siblings, reference))
    # The fixture filters MAF with pooled reference samples, so all frequencies are interior.
    packed = np.empty((len(markers), 3, (len(reference_rows) + 63) // 64), dtype=np.uint64)
    pack_ld_hardcalls(reference, packed)
    reference_weights = 1 / compute_ld_scores_packed(packed, cm, 1.)
    errors = np.full(len(markers), .001)
    results = {}
    for method, frequencies, ld_weights in (
        ('default', child_frequency, weights),
        ('reference_freq', reference_frequency, weights),
        ('reference_ld', child_frequency, reference_weights),
        ('reference_both', reference_frequency, reference_weights),
    ):
        transitions, emissions = _chromosome_hmm(frequencies, cm, ld_weights, errors)
        hard, posterior = _infer_block(pair_rows, siblings, np.arange(len(markers)), transitions,
                                       emissions, 0, _posterior_tables(transitions, emissions))
        if method == 'default':
            np.testing.assert_array_equal(baseline, hard)
        results[method] = calibration(posterior[valid], truth[valid], hard[valid])
    return results


def score_ibd(root, fixture, output):
    pairs, markers, positions, cm, weights = common_markers(output)
    with np.load(root / 'heldout_truth.npz') as truth:
        all_variants = {sid: j for j, sid in enumerate(truth['variant_ids'])}
        all_pairs = {_canonical_pair(*pair): j for j, pair in enumerate(truth['pairs'])}
        lookup = [all_variants[sid] for sid in markers]
        pair_rows = [all_pairs[_canonical_pair(*pair)] for pair in pairs]
        states = truth['states'][np.ix_(pair_rows, lookup)]
    valid = states >= 0
    if not valid.any():
        raise ValueError('No held-out Mendelian labels among the common IBD markers')
    results = {}
    decoded = {}
    methods = ('baseline-pair', 'current-pair', 'phased-statistical')
    for method in methods:
        stem = output / method / 'chr_22'
        segment_file = stem.with_suffix('.ibd.segments.gz')
        if not segment_file.is_file():
            continue
        segments = _segments(segment_file)
        hard = np.stack([_hard_from_segments(segments[_canonical_pair(*pair)], markers, positions)
                         for pair in pairs])
        decoded[method] = hard
        results[method] = {'accuracy': float((hard[valid] == states[valid]).mean()), 'n': int(valid.sum())}
        posterior_file = stem.with_suffix('.ibd.posterior.hdf5')
        if posterior_file.is_file():
            with h5py.File(posterior_file) as source:
                pair_lookup = {_canonical_pair(*pair): row for row, pair in enumerate(source['pairs'].asstr()[:])}
                np.testing.assert_array_equal(source['variant_ids'].asstr()[:], markers)
                posterior_rows = [pair_lookup[_canonical_pair(*pair)] for pair in pairs]
                results[method + '-posterior'] = calibration(source['probabilities'][:][posterior_rows][valid],
                                                              states[valid], hard[valid])
    np.testing.assert_array_equal(decoded['baseline-pair'], decoded['current-pair'])
    results['reference_sensitivity'] = reference_sensitivity(
        fixture, pairs, markers, cm, weights, states, valid, decoded['baseline-pair'],
    )
    return results


def impute(fixture, output, ibd_outputs, *, current, legacy, phase, threads):
    genotype = ['--bed', str(fixture / 'chr_@')] if phase == 'unphased' else [
        '--vcf', str(fixture / 'beagle/beagle.vcf.gz'), '--no-quality-filter']
    methods = {'legacy_pair': ('baseline-pair', legacy, False, ()),
               'pair_hard': ('current-pair', current, False, ('--no_parent_origin',)),
               'pair_posterior': ('current-pair', current, True, ())}
    if phase != 'unphased':
        methods['pair_hard_origin'] = ('current-pair', current, False, ())
        methods['phased_hard'] = ('phased-statistical', current, False, ())
        methods['phased_posterior'] = ('phased-statistical', current, True, ())
    paths = {}
    for method, (source, python, posterior, options) in methods.items():
        out = output / phase
        out.mkdir(parents=True, exist_ok=True)
        target = out / method
        path = out / f'{method}_22.hdf5'
        if not path.is_file():
            hard_input = out / f'{method}.ibd'
            segments = ibd_outputs / source / 'chr_22.ibd.segments.gz'
            link = Path(str(hard_input) + '.segments.gz')
            if not link.is_symlink():
                link.symlink_to(segments.resolve())
            command = [str(python), '-m', 'snipar.scripts.impute', *genotype, '--pedigree',
                       str(fixture / 'pedigree.txt'), '--ibd', str(hard_input), '--chr_range', '22',
                       '--threads', str(threads), '--out', str(target) + '_@', '-silent_progress', *options]
            if posterior:
                command += ['--ibd-posterior', str(ibd_outputs / source / 'chr_22.ibd.posterior.hdf5')]
            with (out / f'{method}.log').open('w') as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        paths[method] = path
    scores = _family_score(paths, fixture)
    return {mask: _mask_results(mask, *moments, 399, 94821) for mask, moments in scores.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--current-python', type=Path, required=True)
    parser.add_argument('--legacy-python', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--phases', nargs='+', choices=('unphased', 'beagle'), default=('unphased', 'beagle'))
    args = parser.parse_args()
    results = {}
    for mask in ('one_parent', 'zero_parents'):
        fixture = args.fixture.resolve() / mask
        output = args.output.resolve() / mask
        ibd = output / 'ibd'
        parameters = SimpleNamespace(fixture=fixture, output=ibd, current_python=args.current_python,
                                     legacy_python=args.legacy_python, rerun=False, p_error=.001,
                                     seed=20261127, threads=args.threads)
        for method in ('baseline-pair', 'current-pair'):
            run_cli(parameters, method)
        if 'beagle' in args.phases:
            run_benchmark_ibd([str(fixture), str(ibd), '--stage', 'phased', '--phase-arms',
                               'phased-statistical', '--threads', str(args.threads), '--skip-scoring'])
        results[mask] = {'ibd': score_ibd(args.fixture, fixture, ibd), 'imputation': {}}
        for phase in args.phases:
            results[mask]['imputation'][phase] = impute(fixture, output, ibd, current=args.current_python,
                                                        legacy=args.legacy_python, phase=phase, threads=args.threads)
            print(f'Empirical {mask} {phase} dosage estimates verified', flush=True)
    (args.output / 'metrics.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
