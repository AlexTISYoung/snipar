"""Derive matched chr22 sibship-size cohorts and known phase-switch controls."""

import argparse
from collections import Counter
from contextlib import ExitStack
import csv
import gzip
from itertools import zip_longest
import json
from pathlib import Path
from shutil import copyfile

import numpy as np

from prepare_1000g_eur import (
    PHASED_GT, read_diploid_gt, vcf_header, write_bed_genotypes, write_vcf_genotypes, write_vcf_header,
)


UNORDERED_GT = np.array([b'0/0\t', b'0/1\t', b'0/1\t', b'1/1\t'])


def nested_metadata(base, output, sizes):
    with (base / 'pedigree.txt').open() as stream:
        header = next(stream)
        pedigree = [line.split() for line in stream]
    children = Counter()
    child_number = {}
    for family, iid, father, mother in pedigree:
        if father != '0' and mother != '0':
            children[family] += 1
            child_number[iid] = children[family]
    if not children or set(children.values()) != {5}:
        raise ValueError('Prepare the source fixture with five children per family')
    fam = [line.split() for line in (base / 'chr_22.fam').read_text().splitlines()]
    sample_ids = np.array([row[1] for row in fam])
    for size in sizes:
        destination = output / f'k{size}'
        include = np.array([child_number.get(iid, 0) <= size for iid in sample_ids])
        if destination != base:
            destination.mkdir(parents=True, exist_ok=True)
            (destination / 'truth').mkdir(exist_ok=True)
            for name in ('chr_22.bim', 'map.txt'):
                copyfile(base / name, destination / name)
            for ending in ('bed', 'bim', 'fam'):
                name = f'truth/chr_22.{ending}'
                copyfile(base / name, destination / name)
            with (destination / 'pedigree.txt').open('w') as stream:
                stream.write(header)
                stream.writelines(' '.join(row) + '\n' for row in pedigree if child_number.get(row[1], 0) <= size)
            with (destination / 'chr_22.fam').open('w') as stream:
                stream.writelines(' '.join(fam[i]) + '\n' for i in np.flatnonzero(include))
            with ((base / 'truth/transmissions.tsv').open() as source,
                  (destination / 'truth/transmissions.tsv').open('w') as stream):
                stream.write(next(source))
                stream.writelines(line for line in source if child_number[line.split('\t', 2)[1]] <= size)
            with ((base / 'relationships.seg').open() as source,
                  (destination / 'relationships.seg').open('w') as stream):
                stream.write(next(source))
                for line in source:
                    first, second, _ = line.split('\t', 2)
                    if max(child_number.get(first, 0), child_number.get(second, 0)) <= size:
                        stream.write(line)
        yield size, destination, np.flatnonzero(include), sample_ids


def prepare(base, switch_rate, seed):
    output = base.parent
    cohorts = list(nested_metadata(base, output, (2, 3, 4, 5)))
    sample_ids = cohorts[0][3]
    rng = np.random.default_rng(seed)
    initial = rng.integers(2, size=len(sample_ids), dtype=np.uint8).astype(bool)
    flipped = initial.copy()
    switch_sites = [[] for _ in sample_ids]
    manifest_base = json.loads((base / 'manifest.json').read_text())
    canonical_unphased = base / 'unphased/normalized.tmp.gz'
    with ExitStack() as stack:
        source = stack.enter_context(gzip.open(base / 'chr_22.vcf.gz', 'rb'))
        if vcf_header(source) != sample_ids.tolist():
            raise ValueError('Source VCF sample order differs from FAM order')
        writers = []
        for size, destination, indices, _ in cohorts:
            (destination / 'unphased').mkdir(exist_ok=True)
            (destination / 'switch').mkdir(exist_ok=True)
            swap_vcf = stack.enter_context(gzip.open(destination / 'switch/chr_22.vcf.gz', 'wb', compresslevel=1))
            path = canonical_unphased if destination == base else destination / 'unphased/chr_22.vcf.gz'
            unphased_vcf = stack.enter_context(gzip.open(path, 'wb', compresslevel=1))
            write_vcf_header(swap_vcf, sample_ids[indices])
            write_vcf_header(unphased_vcf, sample_ids[indices])
            phased = bed = None
            if destination != base:
                phased = stack.enter_context(gzip.open(destination / 'chr_22.vcf.gz', 'wb', compresslevel=1))
                bed = stack.enter_context((destination / 'chr_22.bed').open('wb'))
                write_vcf_header(phased, sample_ids[indices])
                bed.write(b'\x6c\x1b\x01')
            writers.append((indices, unphased_vcf, swap_vcf, phased, bed))
        sample_rows = np.arange(len(sample_ids))
        n_markers = 0
        below_one_percent = np.zeros(len(writers), dtype=int)
        monomorphic = np.zeros(len(writers), dtype=int)
        for line in source:
            fields = line.split(b'\t', 9)
            haplotypes = read_diploid_gt(fields, len(sample_ids), sample_rows)
            if n_markers:
                changed = np.flatnonzero(rng.random(len(sample_ids)) < switch_rate)
                flipped[changed] = ~flipped[changed]
                for row in changed:
                    switch_sites[row].append(n_markers)
            switched = haplotypes.copy()
            switched[flipped] = haplotypes[flipped, ::-1]
            for i, (indices, unphased, swap, phased, bed) in enumerate(writers):
                observed = haplotypes[indices]
                count = observed.sum()
                maf = min(count, observed.size - count) / observed.size
                below_one_percent[i] += maf < .01
                monomorphic[i] += maf == 0
                write_vcf_genotypes(unphased, fields, observed, UNORDERED_GT, count)
                write_vcf_genotypes(swap, fields, switched[indices], PHASED_GT, count)
                if phased is not None:
                    write_vcf_genotypes(phased, fields, observed, PHASED_GT, count)
                    write_bed_genotypes(bed, observed)
            n_markers += 1
    canonical_unphased.replace(base / 'unphased/chr_22.vcf.gz')

    if n_markers != manifest_base['family_snps']:
        raise ValueError('Source VCF marker count differs from its manifest')
    for i, (size, destination, indices, _) in enumerate(cohorts):
        with (destination / 'truth/phase_switches.tsv').open('w') as stream:
            stream.write('IID\tinitial_swapped\tswitch_before_zero_based_marker_index\n')
            for row in indices:
                sites = ','.join(map(str, switch_sites[row])) or '.'
                stream.write(f'{sample_ids[row]}\t{int(initial[row])}\t{sites}\n')
        keys = ('chromosome', 'start', 'end', 'seed', 'eur_samples', 'eur_complete_trios', 'eur_full_sibships',
                'founder_samples', 'simulated_families', 'unrelated_samples', 'families_with_both_parents',
                'families_with_one_parent', 'families_without_parents')
        manifest = {key: manifest_base[key] for key in keys} | {
            'children_per_family': size, 'observed_samples': len(indices), 'family_snps': n_markers,
            'retained_panel': 'k5', 'k5_panel_min_maf': manifest_base.get('min_family_maf', .01),
            'monomorphic_retained_snps': int(monomorphic[i]),
            'retained_snps_below_1_percent_maf': int(below_one_percent[i]),
            'phase_switch_seed': seed, 'phase_switch_probability': switch_rate,
            'phase_switch_events': sum(len(switch_sites[row]) for row in indices),
        }
        pending = destination / '.evaluation-manifest.json.tmp'
        pending.write_text(json.dumps(manifest, indent=2) + '\n')
        pending.replace(destination / 'evaluation-manifest.json')
        if destination != base:
            copyfile(destination / 'evaluation-manifest.json', destination / 'manifest.json')
        print(json.dumps(manifest), flush=True)


def write_oracle_segments(fixture, output):
    """Project simulated transmissions onto the common marker grid in pair order."""
    from snipar.ibd import _SEGMENT_HEADER, find_segments
    from snipar.pedigree import get_sibpairs_from_ped

    bim = np.loadtxt(fixture / 'chr_22.bim', usecols=(1, 2, 3), dtype=str)
    snps = bim[:, 0]
    genetic_map, positions = bim[:, 1].astype(float), bim[:, 2].astype(np.int64)
    pedigree = np.loadtxt(fixture / 'pedigree.txt', dtype=str, skiprows=1, ndmin=2)
    pairs, _ = get_sibpairs_from_ped(pedigree)
    with (fixture / 'truth/transmissions.tsv').open() as stream:
        transmissions = {row['IID']: row for row in csv.DictReader(stream, delimiter='\t')}

    def haplotypes(iid):
        row = transmissions[iid]
        result = np.empty((2, len(genetic_map)), dtype=bool)
        for index, parent in enumerate(('paternal', 'maternal')):
            crossovers = row[f'{parent}_crossover_cM']
            events = np.empty(0) if crossovers == '.' else np.fromstring(crossovers, sep=',')
            result[index] = (int(row[f'{parent}_start']) + np.searchsorted(events, genetic_map, side='right')) % 2
        return result

    output.parent.mkdir(parents=True, exist_ok=True)
    family = None
    paths = {}
    with gzip.open(output, 'wt', compresslevel=3) as stream:
        stream.write(_SEGMENT_HEADER)
        for first, second in pairs:
            if transmissions[first]['FID'] != family:
                family = transmissions[first]['FID']
                paths.clear()
            for iid in (first, second):
                if iid not in paths:
                    paths[iid] = haplotypes(iid)
            states = (paths[first] == paths[second]).sum(axis=0, dtype=np.int8)
            segments = find_segments(states, genetic_map, snps, positions)
            stream.writelines(segment.to_text(first, second, 22) for segment in segments)


def phase_quality(fixture, condition):
    """Measure switches between consecutive heterozygotes against simulated phase."""
    relative_path = 'beagle/beagle.vcf.gz' if condition == 'beagle' else 'switch/chr_22.vcf.gz'
    with (fixture / 'truth/transmissions.tsv').open() as stream:
        children = {row['IID'] for row in csv.DictReader(stream, delimiter='\t')}
    with (gzip.open(fixture / 'chr_22.vcf.gz', 'rb') as truth,
          gzip.open(fixture / relative_path, 'rb') as phased):
        ids, phased_ids = vcf_header(truth), vcf_header(phased)
        if len(phased_ids) != len(ids) or set(phased_ids) != set(ids):
            raise ValueError('The phase test must retain the original study samples')
        is_child = np.array([iid in children for iid in ids])
        phased_id_rows = {iid: row for row, iid in enumerate(phased_ids)}
        truth_rows = np.arange(len(ids))
        phased_rows = np.array([phased_id_rows[iid] for iid in ids])
        switches, comparisons = np.zeros((2, len(ids)), dtype=np.int64)
        seen, previous = np.zeros((2, len(ids)), dtype=bool)
        dose_changes = n_markers = 0
        for true_line, test_line in zip_longest(truth, phased):
            if true_line is None or test_line is None:
                raise ValueError('The phase test must retain all original markers')
            true_fields = true_line.split(b'\t', 9)
            test_fields = test_line.split(b'\t', 9)
            if true_fields[1:5] != test_fields[1:5]:
                raise ValueError('Marker order or alleles differ between phase conditions')
            original = read_diploid_gt(true_fields, len(ids), truth_rows)
            tested = read_diploid_gt(test_fields, len(ids), phased_rows)
            dose_changes += np.count_nonzero(original.sum(axis=1) != tested.sum(axis=1))
            heterozygous = (original[:, 0] != original[:, 1]) & (tested[:, 0] != tested[:, 1])
            orientation = original[:, 0] != tested[:, 0]
            comparable = seen & heterozygous
            switches += comparable & (orientation != previous)
            comparisons += comparable
            seen |= heterozygous
            previous[heterozygous] = orientation[heterozygous]
            n_markers += 1

    groups = {'all': np.ones(len(ids), dtype=bool), 'children': is_child, 'nonchildren': ~is_child}
    result = {
        'phase': condition, 'n_markers': n_markers, 'n_samples': len(ids),
        'genotype_dosage_changes': int(dose_changes),
        'switch_counts': {key: int(switches[mask].sum()) for key, mask in groups.items()},
        'comparison_counts': {key: int(comparisons[mask].sum()) for key, mask in groups.items()},
    }
    result['switch_rates'] = {key: result['switch_counts'][key] / count if count else None
                              for key, count in result['comparison_counts'].items()}
    if dose_changes:
        raise ValueError(f'{condition} changed {dose_changes} genotype dosages')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base', type=Path, help='Existing k5 fixture from prepare_1000g_eur.py')
    parser.add_argument('--switch-rate', type=float, default=.001)
    parser.add_argument('--seed', type=int, default=20261129)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--oracle-out', type=Path, help='Write oracle segments for existing cohorts under this directory')
    modes.add_argument('--phase-quality', choices=('beagle', 'switch'), help='Check phase in existing cohorts')
    args = parser.parse_args()
    if not 0 <= args.switch_rate <= .5 or args.base.name != 'k5':
        parser.error('Use a k5 source directory and a phase-switch rate between zero and 0.5')
    if args.phase_quality is not None:
        for size in (2, 3, 4, 5):
            fixture = args.base.parent / f'k{size}'
            result = phase_quality(fixture, args.phase_quality)
            (fixture / f'{args.phase_quality}_phase_quality.json').write_text(json.dumps(result, indent=2) + '\n')
            print(json.dumps({'cohort': fixture.name, **result}), flush=True)
    elif args.oracle_out is None:
        prepare(args.base.resolve(), args.switch_rate, args.seed)
    else:
        for size in (2, 3, 4, 5):
            source = args.base.parent / f'k{size}'
            output = args.oracle_out / f'k{size}/oracle/chr_22.ibd.segments.gz'
            write_oracle_segments(source, output)


if __name__ == '__main__':
    main()
