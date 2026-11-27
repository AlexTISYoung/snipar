"""Profile the SNIPar commands on repeatable family-genetics workloads.

Run from a Python environment with SNIPar's dependencies and native extensions.
The optional --source checkout needs its own compiled extensions. PLINK is not
required. R profiling requires Rscript, lme4 and numDeriv. JSON and profiler
outputs go in the selected output directory.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed, SnpData


def run(command, cwd, env, log):
    start = perf_counter()
    with log.open('w') as output:
        result = subprocess.run(command, cwd=cwd, env=env, stdout=output, stderr=subprocess.STDOUT)
    return perf_counter() - start, result.returncode


def trio_data(directory):
    """Create a clean shared GRM fixture with 1,500 trios and 1,000 SNPs."""
    rng = np.random.default_rng(1944)
    n, m = 1500, 1000
    freqs = rng.uniform(.05, .45, m)
    paternal = rng.random((n, 2, m)) < freqs
    maternal = rng.random((n, 2, m)) < freqs
    gts = np.concatenate((paternal[:, 0].astype(np.float32) + maternal[:, 0],
                          paternal.sum(axis=1), maternal.sum(axis=1))).astype(np.float32)
    family = np.array([f'F{i:04}' for i in range(n)])
    ids = np.char.add(np.tile(family, 3), np.repeat(['_kid', '_dad', '_mom'], n))
    pos = np.column_stack((np.ones(m), np.arange(m) / 100, np.arange(m) + 1))
    snps = SnpData(iid=np.column_stack((np.tile(family, 3), ids)),
                   sid=np.array([f'rs{i + 1}' for i in range(m)]), pos=pos, val=gts)
    Bed.write(str(directory / 'trios'), snps, count_A1=True)
    np.savetxt(directory / 'trios.ped', np.column_stack((family, ids.reshape(3, n).T)),
               fmt='%s', header='FID IID FATHER_ID MOTHER_ID', comments='')


def measure(label, command, directory, output, env, repetitions, python=True):
    execution = [sys.executable, *command] if python else command
    _, status = run(execution, directory, env, output / f'{label}.warm.log')
    row = {'command': command, 'exit_status': status, 'seconds': []}
    if status:
        return row
    for i in range(repetitions):
        elapsed, status = run(execution, directory, env, output / f'{label}.{i + 1}.log')
        if status:
            row['exit_status'] = status
            return row
        row['seconds'].append(round(elapsed, 6))
    if python:
        profiled = [sys.executable, '-m', 'cProfile', '-o', str(output / f'{label}.prof'), *command]
    else:
        # pgs.R reads the trailing arguments. Pass only its own options to main.
        script = command[1]
        code = ('args <- commandArgs(TRUE); Rprof(args[1]); '
                'source(args[2]); pgs_binary_main(args[-c(1,2)]); Rprof(NULL)')
        profiled = [command[0], '-e', code, str(output / f'{label}.Rprof'), script, *command[2:]]
    _, row['profile_exit_status'] = run(profiled, directory, env, output / f'{label}.profile.log')
    row['median_seconds'] = round(float(np.median(row['seconds'])), 6)
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repetitions', type=int, default=3)
    args = parser.parse_args(argv)
    if args.repetitions < 1:
        parser.error('--repetitions must be positive')
    source, output = args.source.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PYTHONPATH=str(source), OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1',
               OMP_NUM_THREADS='1', NUMBA_NUM_THREADS='4', NUMBA_CACHE_DIR=str(output / 'numba-cache'))
    results = {}

    with TemporaryDirectory(prefix='snipar-profile-') as temp:
        directory = Path(temp)
        # A distinct destination for every archive copy is required by the CLI.
        # Remove it after each run to keep the same CLI arguments in the timing loop.
        command = ['-m', 'snipar.example.snipar_example_data', '--dest', 'example_data']
        results['snipar_example_data'] = {'command': command, 'seconds': []}
        for i in range(args.repetitions + 2):
            profiled = i == args.repetitions + 1
            options = ['-m', 'cProfile', '-o', str(output / 'snipar_example_data.prof')] if profiled else []
            elapsed, status = run([sys.executable, *options, *command], directory, env,
                                  output / f'snipar_example_data.{i}.log')
            if status:
                raise RuntimeError(f'snipar_example_data failed, see {output}')
            if 0 < i <= args.repetitions:
                results['snipar_example_data']['seconds'].append(round(elapsed, 6))
            if not profiled:
                shutil.rmtree(directory / 'example_data')
        example = directory / 'example_data'
        trio_data(directory)
        (directory / 'simulation').mkdir()
        commands = [
            ('ibd', example, ['--bed', 'chr_@', '--king', 'king.kin0', '--agesex', 'agesex.txt',
                              '--out', 'chr_@', '--threads', '4', '--ld_out']),
            ('impute', example, ['--ibd', 'chr_@.ibd', '--bed', 'chr_@', '--king', 'king.kin0',
                                 '--agesex', 'agesex.txt', '--out', 'chr_@', '--threads', '4']),
            ('gwas', example, ['phenotype.txt', '--bed', 'chr_@', '--imp', 'chr_@', '--out', 'chr_@_young', '--cpus', '1']),
            ('correlate', example, ['chr_@_young', 'effect', '--ldscores', 'chr_@', '--save_delete']),
            ('pgs', example, ['direct', '--bed', 'chr_@', '--imp', 'chr_@', '--weights', 'direct_weights.txt',
                              '--SNP', 'sid', '--beta_col', 'ldpred_beta', '--A1', 'nt1', '--A2', 'nt2',
                              '--phenofile', 'phenotype.txt']),
            ('simulate', directory, ['500', '.5', 'simulation/', '--nfam', '1000', '--n_am', '3', '--r_par', '.5',
                                     '--save_par_gts']),
            ('gts_matrix', example, ['--bed', 'chr_@', '--imp', 'chr_@', '--out', 'genotype_matrix']),
            ('make_rdr_grms', directory, ['trios.bed', 'rdr', '--ped', 'trios.ped']),
        ]
        for name, cwd, argv in commands:
            results[name] = measure(name, ['-m', f'snipar.scripts.{name}', *argv], cwd,
                                    output, env, args.repetitions)
            print(name, results[name].get('median_seconds', f"failed ({results[name]['exit_status']})"), flush=True)
        results['estimate_sim_effects'] = measure(
            'estimate_sim_effects', ['estimate_sim_effects.py', 'chr_1_young.sumstats.hdf5', 'phenotype.effects.txt'],
            example, output, env, args.repetitions,
        )
        results['simulate_phenotype'] = measure(
            'simulate_phenotype', ['-m', 'snipar.example.simulate_phenotype', '--bed', 'trios',
                                   '--ncausal', '500', '--seed', '1944'], directory, output, env, args.repetitions,
        )
        if shutil.which('Rscript') and (example / 'direct.pgs.txt').exists():
            scores = pd.read_csv(example / 'direct.pgs.txt', sep=r'\s+', dtype=str)
            phenotype = pd.read_csv(example / 'phenotype.txt', sep=r'\s+', dtype=str)
            phenotype = scores[['FID', 'IID']].merge(phenotype[['IID', 'case_status']], on='IID')
            phenotype.to_csv(example / 'matched_phenotype.txt', sep='\t', index=False)
            results['pgs.R'] = measure('pgs.R', ['Rscript', str(source / 'snipar/scripts/pgs.R'), 'direct_binary',
                                      '--pgs', 'direct.pgs.txt', '--phenofile', 'matched_phenotype.txt',
                                      '--phen', 'case_status'], example, output, env, args.repetitions, python=False)
    (output / 'timings.json').write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
