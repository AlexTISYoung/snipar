"""Prepare repeated masking experiments in the five complete 1kGP sibling quads."""

import argparse
from contextlib import ExitStack
import gzip
import json
from pathlib import Path

import numpy as np

from calibrate_ibd_posteriors import load_panel, parental_sharing
from prepare_1000g_eur import (
    BASES, UNPHASED_GT, read_diploid_gt, vcf_header, write_bed_genotypes,
    write_vcf_genotypes, write_vcf_header,
)
from snipar.vcf import _info_maf


def prepare(source, pedigree, mapfile, output, start, end, reference_size, thinning_bp, seed):
    bp, cm = np.loadtxt(mapfile, skiprows=1, usecols=(0, 2), unpack=True)
    observed_sites, anchors, truth_positions, truth_ids = [], [], [], []
    codes = UNPHASED_GT.copy()
    codes[2] = codes[1]
    prefixes = ('one_parent', 'zero_parents')
    output.mkdir(parents=True)
    with ExitStack() as stack:
        panel = stack.enter_context(gzip.open(source, 'rb'))
        sample_ids = vcf_header(panel)
        quads, references = load_panel(pedigree, sample_ids, reference_size, seed)
        families = [family for family, _ in quads]
        empirical_ids = [iid for family in families for iid in family]
        reference_ids = [iid for cohort in references.values() for iid in cohort]
        iids = empirical_ids + reference_ids
        source_rows = {iid: row for row, iid in enumerate(sample_ids)}
        sample_rows = {iid: row for row, iid in enumerate(iids)}
        indices = np.array([source_rows[iid] for iid in iids])
        parents = np.array([[sample_rows[father], sample_rows[mother]] for father, mother, *_ in families])
        configurations = []
        for mask in prefixes:
            run = output / mask
            (run / 'truth').mkdir(parents=True)
            bed = stack.enter_context((run / 'chr_22.bed').open('wb'))
            truth_bed = stack.enter_context((run / 'truth' / 'chr_22.bed').open('wb'))
            bim = stack.enter_context((run / 'chr_22.bim').open('w'))
            vcf = stack.enter_context(gzip.open(run / 'unphased.vcf.gz', 'wb', compresslevel=1))
            bed.write(b'\x6c\x1b\x01')
            truth_bed.write(b'\x6c\x1b\x01')
            observed = []
            with ((run / 'pedigree.txt').open('w') as ped,
                  (run / 'chr_22.fam').open('w') as fam,
                  (run / 'truth' / 'chr_22.fam').open('w') as truth_fam):
                ped.write('FID\tIID\tFATHER_ID\tMOTHER_ID\n')
                for index, (father, mother, *children) in enumerate(families):
                    fid = f'quad{index + 1}'
                    truth_fam.write(f'{fid} {father} 0 0 1 -9\n{fid} {mother} 0 0 2 -9\n')
                    if mask == 'one_parent':
                        parent = (father, mother)[index % 2]
                        observed.append(parent)
                        fam.write(f'{fid} {parent} 0 0 {index % 2 + 1} -9\n')
                    for child in children:
                        observed.append(child)
                        ped.write(f'{fid}\t{child}\t{father}\t{mother}\n')
                        fam.write(f'{fid} {child} {father} {mother} 0 -9\n')
                # Independent samples enter the allele-frequency fit as singleton
                # candidates with unique, ungenotyped parent identifiers.
                for refindex, iid in enumerate(reference_ids):
                    fid = f'freq{refindex}'
                    observed.append(iid)
                    ped.write(f'{fid}\t{iid}\t{fid}_1\t{fid}_2\n')
                    fam.write(f'{fid} {iid} 0 0 0 -9\n')
            write_vcf_header(vcf, observed)
            rows = np.array([sample_rows[iid] for iid in observed])
            configurations.append((run, bed, truth_bed, bim, vcf, rows))

        next_site = start
        for record in panel:
            fields = record.split(b'\t', 9)
            position = int(fields[1])
            if position < start:
                continue
            if position > end:
                break
            if fields[3] not in BASES or fields[4] not in BASES or fields[6] not in (b'.', b'PASS'):
                continue
            maf = _info_maf(fields[7])
            if maf is not None and maf < .001:
                continue
            if position % 5 == 0:
                family_calls = read_diploid_gt(fields, len(sample_ids), indices[:len(empirical_ids)])
                if np.all(family_calls < 2):
                    anchors.append((position, family_calls.sum(axis=1, dtype=np.int8)))
                continue
            if position < next_site:
                continue
            haplotypes = read_diploid_gt(fields, len(sample_ids), indices)
            if np.any(haplotypes > 1):
                continue
            reference = haplotypes[len(empirical_ids):]
            maf = reference.mean()
            if min(maf, 1 - maf) < .01:
                continue
            sid = b':'.join((b'22', fields[1], fields[3], fields[4]))
            fields[2] = sid
            genetic = np.interp(position, bp, cm)
            for _, bed, truth_bed, bim, vcf, observed in configurations:
                called = haplotypes[observed]
                write_bed_genotypes(bed, called)
                write_bed_genotypes(truth_bed, haplotypes[parents.ravel()])
                write_vcf_genotypes(vcf, fields, called, codes, int(called.sum()))
                bim.write(f'22\t{sid.decode()}\t{genetic:.9f}\t{position}\t{fields[4].decode()}\t{fields[3].decode()}\n')
            observed_sites.append((position, genetic))
            truth_positions.append(position)
            truth_ids.append(sid.decode())
            next_site = position + thinning_bp
    if not observed_sites or not anchors:
        raise ValueError('No evaluation or held-out markers in the interval')

    for run, *_ in configurations:
        (run / 'truth' / 'chr_22.bim').write_bytes((run / 'chr_22.bim').read_bytes())
        np.savetxt(run / 'map.txt', observed_sites, fmt=('%d', '%.9f'), header='pposition gposition', comments='')
    anchor_positions, anchor_calls = zip(*anchors)
    anchor_positions = np.asarray(anchor_positions)
    anchor_calls = np.stack(anchor_calls)
    positions = np.asarray(truth_positions)
    states, child_pairs = [], []
    for father, mother, first, second in families:
        cols = [sample_rows[iid] for iid in (father, mother, first, second)]
        pat, _ = parental_sharing(anchor_positions, anchor_calls[:, cols], positions, 100_000)
        mat, _ = parental_sharing(anchor_positions, anchor_calls[:, np.asarray(cols)[[1, 0, 2, 3]]],
                                  positions, 100_000)
        states.append(np.where((pat >= 0) & (mat >= 0), pat + mat, -1))
        child_pairs.append((first, second))
    np.savez_compressed(output / 'heldout_truth.npz', pairs=np.asarray(child_pairs),
                        variant_ids=np.asarray(truth_ids), positions=positions, states=np.stack(states))
    summary = {'families': len(families), 'masks': list(prefixes), 'start': start, 'end': end,
               'evaluation_markers': len(positions), 'heldout_anchors': len(anchors),
               'labelled_pairs_markers': int((np.stack(states) >= 0).sum()),
               'references': {k: len(v) for k, v in references.items()}, 'seed': seed}
    (output / 'manifest.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'pedigree', 'map', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--start', type=int, default=17_000_000)
    parser.add_argument('--end', type=int, default=25_000_000)
    parser.add_argument('--reference-size', type=int, default=64)
    parser.add_argument('--thinning-bp', type=int, default=300)
    parser.add_argument('--seed', type=int, default=20261128)
    args = parser.parse_args()
    if not 0 < args.start < args.end or min(args.reference_size, args.thinning_bp) < 1:
        parser.error('A valid interval and positive reference and thinning sizes are required')
    prepare(args.source, args.pedigree, args.map, args.output, args.start, args.end,
            args.reference_size, args.thinning_bp, args.seed)


if __name__ == '__main__':
    main()
