#!/usr/bin/env python3
"""Measure sibling IBD calibration against held-out Mendelian transmissions in 1kGP quads."""

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path

import numpy as np
from numba import config, set_num_threads

from prepare_1000g_eur import BASES, read_diploid_gt, vcf_header
from snipar.ibd import _chromosome_hmm, _infer_block, _posterior_tables
from snipar.ld import compute_ld_scores_packed, pack_ld_hardcalls
from snipar.map import _linear_endpoints
from snipar.resources import available_cpus
from snipar.vcf import _info_maf


def load_panel(pedigree, sample_ids, reference_size, seed):
    rows = np.loadtxt(pedigree, dtype=str, skiprows=1)
    available = set(sample_ids)
    offspring = defaultdict(list)
    related_families = set(rows[(rows[:, 2:4] != '0').any(axis=1), 0])
    parents = set(rows[:, 2:4].flat)
    unrelated = defaultdict(dict)
    for family, iid, father, mother, _, _, ancestry in rows:
        if iid in available and father in available and mother in available:
            offspring[father, mother].append((iid, ancestry))
        if iid in available and family not in related_families and iid not in parents:
            unrelated[ancestry].setdefault(family, iid)
    quads = [(parents + tuple(row[0] for row in children), children[0][1])
             for parents, children in offspring.items() if len(children) == 2]
    if not quads:
        raise ValueError('No complete parent–sibling quads found in the genotype source')
    rng = np.random.default_rng(seed)
    references = {}
    for ancestry in sorted({population for _, population in quads}):
        candidates = sorted(unrelated[ancestry].values())
        references[ancestry] = rng.choice(candidates, min(reference_size, len(candidates)), replace=False)
        if len(references[ancestry]) < 32:
            raise ValueError(f'Fewer than 32 independent allele-frequency references for {ancestry}')
    return quads, references


def read_calls(input_vcf, pedigree, reference_size, seed, thinning_bp):
    """Separate anchor and HMM markers before evaluating any transmissions."""
    with gzip.open(input_vcf, 'rb') as source:
        sample_ids = vcf_header(source)
        quads, references = load_panel(pedigree, sample_ids, reference_size, seed)
        quad_ids = list(dict.fromkeys(iid for family, _ in quads for iid in family))
        selected = quad_ids + [iid for cohort in references.values() for iid in cohort]
        rows = {iid: index for index, iid in enumerate(selected)}
        source_rows = {iid: index for index, iid in enumerate(sample_ids)}
        indices = np.asarray([source_rows[iid] for iid in selected])
        anchor_positions, anchor_calls, positions, calls = [], [], [], []
        next_site = 0
        for record in source:
            fields = record.split(b'\t', 9)
            if fields[3] not in BASES or fields[4] not in BASES or fields[6] not in (b'.', b'PASS'):
                continue
            position = int(fields[1])
            anchor = position % 5 == 0
            if not anchor and position < next_site:
                continue
            maf = _info_maf(fields[7])
            if maf is not None and maf < .001:
                continue
            selected_rows = indices[:len(quad_ids)] if anchor else indices
            alleles = read_diploid_gt(fields, len(sample_ids), selected_rows)
            dosages = np.where((alleles <= 1).all(axis=1), alleles.sum(axis=1, dtype=np.int16), -1).astype(np.int8)
            if anchor:
                anchor_positions.append(position)
                anchor_calls.append(dosages)
            else:
                positions.append(position)
                calls.append(dosages)
                next_site = position + thinning_bp
    if not positions or not anchor_positions:
        raise ValueError('No qualifying training or held-out variants in the VCF')
    return quads, references, rows, *(np.asarray(a) for a in (positions, calls, anchor_positions, anchor_calls))


def parental_sharing(anchors, dosages, positions, max_gap):
    """Label the sharing bit when two informative adjacent anchors agree."""
    informative = ((dosages[:, 0] == 1) & np.isin(dosages[:, 1], [0, 2]))
    carried = dosages[:, 2:] - dosages[:, 1:2] // 2
    informative &= ((carried >= 0) & (carried <= 1)).all(axis=1) & (dosages[:, 2:] >= 0).all(axis=1)
    direct = carried[informative, 0] == carried[informative, 1]
    anchors = anchors[informative]
    state = np.full(len(positions), -1, dtype=np.int8)
    if len(anchors) < 2:
        return state, len(anchors)
    right = np.searchsorted(anchors, positions, side='right')
    covered = np.flatnonzero((right > 0) & (right < len(anchors)))
    right = right[covered]
    keep = (direct[right - 1] == direct[right]) & (anchors[right] - anchors[right - 1] <= max_gap)
    state[covered[keep]] = direct[right[keep]]
    return state, len(anchors)


def calibration(probability, truth, viterbi):
    winners = probability.argmax(axis=1)
    confidence = probability.max(axis=1)
    correct = winners == truth
    bins = []
    edges = [0, .5, .7, .8, .9, .95, .99, 1.00001]
    for low, high in zip(edges[:-1], edges[1:]):
        selected = (confidence >= low) & (confidence < high)
        if np.any(selected):
            bins.append({'lower': low, 'upper': min(high, 1), 'n': int(selected.sum()),
                         'confidence': float(confidence[selected].mean(dtype=np.float64)),
                         'accuracy': float(correct[selected].mean())})
    n = len(truth)
    true_prob = probability[np.arange(n), truth].astype(np.float64)
    one_hot = np.eye(3)[truth]
    return {'n': n, 'accuracy': float(correct.mean()), 'viterbi_accuracy': float((viterbi == truth).mean()),
            'brier': float(np.square(probability - one_hot).sum(axis=1).mean()),
            'log_loss': float(-np.log(np.maximum(true_prob, np.finfo(np.float32).tiny)).mean()),
            'ece': sum(row['n'] * abs(row['confidence'] - row['accuracy']) for row in bins) / n,
            'observed_state_fractions': one_hot.mean(axis=0).tolist(),
            'expected_state_fractions': probability.mean(axis=0, dtype=np.float64).tolist(), 'bins': bins}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf', type=Path, help='NYGC 3,202-sample GT-only autosomal VCF for one chromosome.')
    parser.add_argument('pedigree', type=Path, help='1kGP pedigree and population table.')
    parser.add_argument('map', type=Path, help='deCODE GRCh38 gzip map with pos/chr/cM columns.')
    parser.add_argument('output', type=Path, help='Output JSON file.')
    parser.add_argument('--reference-size', type=int, default=128)
    parser.add_argument('--seed', type=int, default=20261128)
    parser.add_argument('--thinning-bp', type=int, default=500)
    parser.add_argument('--max-anchor-gap', type=int, default=100_000)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if min(args.reference_size, args.thinning_bp, args.max_anchor_gap, args.threads) < 1:
        parser.error('Reference count, distances and threads must be positive')
    set_num_threads(min(args.threads, available_cpus(), config.NUMBA_NUM_THREADS))
    quads, references, rows, positions, calls, anchors, anchor_calls = read_calls(
        args.vcf, args.pedigree, args.reference_size, args.seed, args.thinning_bp,
    )
    print(f'{len(quads)} quads, {len(positions):,} candidate HMM markers, {len(anchors):,} holdout sites', flush=True)
    bp, cm = np.loadtxt(args.map, skiprows=1, usecols=(0, 2), unpack=True)
    full_map = np.interp(positions, bp, cm)
    _linear_endpoints(positions, full_map, bp, cm)
    frequencies, markers_by_ancestry, ld_by_ancestry = {}, {}, {}
    for ancestry, iids in references.items():
        genotypes = calls[:, [rows[iid] for iid in iids]]
        called = genotypes >= 0
        freq = np.maximum(genotypes, 0).sum(axis=1) / np.maximum(2 * called.sum(axis=1), 1)
        markers = np.flatnonzero((freq > .01) & (freq < .99) & (called.mean(axis=1) >= .95))
        packed = np.empty((len(markers), 3, (len(iids) + 63) // 64), dtype=np.uint64)
        pack_ld_hardcalls(genotypes[markers].T, packed)
        genetic_map = full_map[markers]
        frequencies[ancestry], markers_by_ancestry[ancestry] = freq[markers], markers
        ld_by_ancestry[ancestry] = compute_ld_scores_packed(packed, genetic_map, 1.)
    results, labelled, predictions, paths = [], [], [], []
    for family, ancestry in quads:
        marker = markers_by_ancestry[ancestry]
        pos = positions[marker]
        family_rows = [rows[iid] for iid in family]
        reference_rows = anchor_calls[:, family_rows]
        paternal, n_paternal = parental_sharing(anchors, reference_rows, pos, args.max_anchor_gap)
        maternal, n_maternal = parental_sharing(anchors, reference_rows[:, [1, 0, 2, 3]], pos, args.max_anchor_gap)
        valid = (paternal >= 0) & (maternal >= 0)
        if not np.any(valid):
            raise ValueError(f'No matching paternal and maternal transmission labels in {family}')
        genotypes = np.ascontiguousarray(calls[marker][:, family_rows[2:]].T)
        genetic_map = full_map[marker]
        transitions, emissions = _chromosome_hmm(frequencies[ancestry], genetic_map,
                                                 1 / ld_by_ancestry[ancestry], np.full(len(pos), .001))
        hard, posterior = _infer_block(np.array([[0, 1]]), genotypes, np.arange(len(pos)),
                                       transitions, emissions, 0, _posterior_tables(transitions, emissions))
        truth, probability, viterbi = (paternal + maternal)[valid], posterior[0, valid], hard[0, valid]
        results.append({'ids': family, 'ancestry': ancestry, 'hmm_markers': len(pos),
                        'paternal_anchors': n_paternal, 'maternal_anchors': n_maternal,
                        **calibration(probability, truth, viterbi)})
        labelled.append(truth)
        predictions.append(probability)
        paths.append(viterbi)
    summary = {'config': {'reference_size': {k: len(v) for k, v in references.items()}, 'seed': args.seed,
                          'thinning_bp': args.thinning_bp, 'max_anchor_gap': args.max_anchor_gap},
               'quads': results, 'pooled': calibration(np.concatenate(predictions), np.concatenate(labelled),
                                                        np.concatenate(paths))}
    args.output.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary['pooled'], indent=2))


if __name__ == '__main__':
    main()
