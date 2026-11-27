"""Time IBD, imputation, GWAS, score calculation and score analysis on one EUR fixture."""

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter

import h5py
import numpy as np
import pandas as pd
import psutil


def timed_command(command, log, environment):
    start = perf_counter()
    peak_rss = 0
    include_children = True
    with log.open('w') as output:
        with subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT, env=environment) as process:
            parent = psutil.Process(process.pid)
            while process.poll() is None:
                rss = 0
                children = []
                if include_children:
                    try:
                        children = parent.children(recursive=True)
                    except psutil.NoSuchProcess:
                        pass
                    except (psutil.AccessDenied, PermissionError):
                        include_children = False
                for child in [parent, *children]:
                    try:
                        rss += child.memory_info().rss
                    except psutil.NoSuchProcess:
                        pass
                peak_rss = max(peak_rss, rss)
                try:
                    process.wait(timeout=.05)
                except subprocess.TimeoutExpired:
                    pass
        if process.returncode:
            raise RuntimeError(f'{log.stem} exited with status {process.returncode}. See {log}')
    return {'seconds': round(perf_counter() - start, 3), 'peak_rss_mib': round(peak_rss / 2**20, 1),
            'rss_scope': 'process_tree' if include_children else 'main_process'}


def pipeline(data, out, genotype, chrom, workers, threads):
    pedigree = ['--pedigree', str(data / 'pedigree.txt')]
    phenotype = str(data / 'phenotype.txt')
    chromosomes = ['--chr_range', str(chrom)]
    native_threads = ['--threads', str(threads)]
    ibd = str(out / 'ibd_@')
    imputation = str(out / 'imputed_@')
    score = str(out / 'scores')
    analysis = str(out / 'score_analysis')
    return {
        'ibd': ('ibd', [*genotype, *chromosomes, *pedigree, *native_threads,
                       '--map', str(data / 'map.txt'), '--p_error', '.001', '--out', ibd]),
        'imputation': ('impute', [*genotype, *chromosomes, *pedigree, *native_threads,
                                  '--ibd', ibd + '.ibd', '--out', imputation, '-silent_progress']),
        'gwas': ('gwas', [phenotype, *genotype, *chromosomes, *native_threads, '--cpus', str(workers),
                          '--imp', imputation, '--out', str(out / 'gwas_@')]),
        'pgs_calculation': ('pgs', [score, *genotype, *chromosomes, *native_threads,
                                    '--imp', imputation, '--weights', str(data / 'weights.txt'), '--no_am_adj']),
        'pgs_analysis': ('pgs', [analysis, '--pgs', score + '.pgs.txt', '--phenofile', phenotype,
                                 *native_threads, '--no_am_adj']),
    }


def gwas_modes(data, out, genotype, chrom, workers, threads):
    common = [str(data / 'phenotype.txt'), *genotype, '--chr_range', str(chrom),
              '--threads', str(threads), '--cpus', str(workers)]
    observed = ['--pedigree', str(data / 'pedigree.txt')]
    imputed = ['--imp', str(out / 'imputed_@')]
    grm = ['--grm', str(data / 'relationships.seg'), '--sparse_thresh', '.1']
    designs = {
        'trio_sib_meta': observed,
        'sib_diff': [*observed, '--sib_diff'],
        'impute_unrel': [*imputed, '--impute_unrel'],
        'impute_unrel_observed': [*observed, '--impute_unrel'],
        'robust': [*imputed, '--robust'],
        'parental_sum': [*imputed, '--parsum'],
        'fit_sib': [*imputed, '--fit_sib'],
        'grm': [*imputed, *grm, '--covar', str(data / 'covariates.txt')],
        'grm_only': [*imputed, *grm, '--no_sib_var', '--vc_list', '.4', '.6'],
    }
    for name, options in designs.items():
        prefix = f'gwas_{name}_@'
        yield name, [*common, *options, '--out', str(out / prefix)], out / prefix.replace('@', str(chrom))


def step_command(module, options, profile=None):
    command = [sys.executable]
    if module == 'ibd':
        seeded = ('import runpy, numpy as np; np.random.seed(20261127); '
                  'runpy.run_module("snipar.scripts.ibd", run_name="__main__", alter_sys=True)')
        if profile is not None:
            seeded = f'import cProfile; cProfile.run({seeded!r}, {str(profile)!r})'
        return [*command, '-c', seeded, *options]
    if profile is not None:
        command += ['-m', 'cProfile', '-o', str(profile)]
    return [*command, '-m', f'snipar.scripts.{module}', *options]


def output_checks(out, chrom):
    segments = pd.read_csv(out / f'ibd_{chrom}.ibd.segments.gz', sep=r'\s+')
    with h5py.File(out / f'imputed_{chrom}.hdf5') as imputation:
        n_families, n_imputed_variants = imputation['imputed_par_gts'].shape
        # Keep the validation memory bounded independently of panel width.
        imputed = imputation['imputed_par_gts']
        finite_imputations = sum(np.isfinite(imputed[:, start:start + 2048]).sum().item()
                                 for start in range(0, n_imputed_variants, 2048))
    stats = pd.read_csv(out / f'gwas_{chrom}.sumstats.gz', sep=r'\s+', usecols=['direct_Beta'])
    scores = pd.read_csv(out / 'scores.pgs.txt', sep=r'\s+')
    effects = pd.read_csv(out / 'score_analysis.effects.txt', sep=r'\s+').set_index('statistic')
    checks = {'ibd_segments': len(segments), 'imputed_families': n_families,
              'imputation_variants': n_imputed_variants, 'finite_imputations': finite_imputations,
              'gwas_variants': len(stats), 'finite_direct_effects': stats.direct_Beta.notna().sum().item(),
              'pgs_rows': len(scores), 'finite_proband_scores': scores.proband.notna().sum().item(),
              'pgs_direct_effect': float(effects.loc['2.proband', 'est'])}
    if not all(checks[key] for key in ('ibd_segments', 'finite_imputations', 'finite_direct_effects',
                                       'finite_proband_scores')) or not np.isfinite(checks['pgs_direct_effect']):
        raise ValueError(f'Pipeline did not produce all expected results in {out}')
    return checks


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture', type=Path)
    parser.add_argument('output', type=Path)
    genotypes = parser.add_mutually_exclusive_group(required=True)
    genotypes.add_argument('--bed', help='BED prefix, with @ for chromosome substitution.')
    genotypes.add_argument('--bgen', help='BGEN prefix, with @ for chromosome substitution.')
    genotypes.add_argument('--vcf', help='VCF filename or @ chromosome pattern.')
    parser.add_argument('--vcf-min-info-maf', type=float, default=0.)
    parser.add_argument('--chromosome', type=int, default=22)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--workers', type=int, default=1, help='GWAS workers.')
    parser.add_argument('--repeat', type=int, default=1)
    parser.add_argument('--profile', action='store_true', help='Save cProfile data and time the profiled run separately.')
    parser.add_argument('--gwas-modes', action='store_true', help='Also run and check the supported GWAS estimator modes.')
    args = parser.parse_args(argv)
    if min(args.chromosome, args.threads, args.workers, args.repeat) < 1:
        parser.error('Chromosome, threads, workers and repetitions must be positive')
    flag, source = next((f'--{kind}', getattr(args, kind)) for kind in ('bed', 'bgen', 'vcf') if getattr(args, kind))
    genotype = [flag, str(Path(source).absolute())]
    if args.vcf_min_info_maf:
        if args.vcf is None:
            parser.error('--vcf-min-info-maf requires VCF input')
        genotype += ['--vcf-min-info-maf', str(args.vcf_min_info_maf)]
    environment = os.environ | dict.fromkeys(('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                               'VECLIB_MAXIMUM_THREADS'), '1')
    environment['NUMBA_NUM_THREADS'] = str(args.threads)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    result = {'python': platform.python_version(), 'platform': platform.platform(),
              'source': genotype, 'workers': args.workers, 'threads': args.threads, 'runs': []}
    for repeat in range(args.repeat + args.profile):
        profiling = repeat == args.repeat
        out = root / ('profile' if profiling else f'run_{repeat + 1}')
        out.mkdir()
        run = {'profiled': profiling, 'steps': {}}
        for step, (module, options) in pipeline(args.fixture.resolve(), out, genotype, args.chromosome,
                                                args.workers, args.threads).items():
            command = step_command(module, options, out / f'{step}.prof' if profiling else None)
            measurement = timed_command(command, out / f'{step}.log', environment)
            run['steps'][step] = measurement | {'command': command}
            print(f'{out.name} {step.ljust(17)} {measurement["seconds"]:8.2f} s '
                  f'{measurement["peak_rss_mib"]:7.1f} MiB ({measurement["rss_scope"]})', flush=True)
        run['checks'] = output_checks(out, args.chromosome)
        result['runs'].append(run)
        (root / 'timings.json').write_text(json.dumps(result, indent=2) + '\n')
        if args.gwas_modes and not profiling:
            methods = {}
            for name, options, prefix in gwas_modes(args.fixture.resolve(), out, genotype, args.chromosome,
                                                    args.workers, args.threads):
                log = out / f'gwas_{name}.log'
                command = step_command('gwas', options)
                measurement = timed_command(command, log, environment)
                stats = pd.read_csv(f'{prefix}.sumstats.gz', sep=r'\s+', usecols=['direct_Beta'])
                estimated = int(stats.direct_Beta.notna().sum())
                if not estimated:
                    raise ValueError(f'No finite direct estimates for GWAS {name}. See {log}')
                methods[name] = measurement | {'command': command, 'finite_direct_effects': estimated}
                print(f'{out.name} {name.ljust(23)} {measurement["seconds"]:8.2f} s '
                      f'{estimated:,} finite direct effects', flush=True)
            run['gwas_modes'] = methods
            (root / 'timings.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
