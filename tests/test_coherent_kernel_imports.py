"""CLI and package imports must share Numba caches outside the project cwd."""

import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1] / "snipar"
pytestmark = pytest.mark.skipif(importlib.util.find_spec('numba') is None,
                                reason='Cross-entry-point cache regression requires Numba')


@pytest.mark.parametrize('first', ['cli', 'package'])
def test_fresh_process_cross_entry_point_cache(first, tmp_path):
    source_root = tmp_path / 'source'
    library = source_root / 'snipar'
    library.mkdir(parents=True)
    for name in ('__init__.py', 'coherent.py', 'coherent_kernels.py', 'vcf_cli.py'):
        shutil.copy2(ROOT / name, library / name)
    scripts = source_root / 'bin'
    scripts.mkdir()
    shutil.copy2(ROOT / 'scripts/correlate.py', scripts / 'correlate.py')
    outside = tmp_path / 'outside'
    outside.mkdir()
    rng = np.random.default_rng(2401)
    n = 40
    sampling = np.array([[.1, -.03], [-.03, .15]])
    signal = np.array([[.8, .1], [.1, .4]])
    effects = rng.normal(size=(n, 2)) @ np.linalg.cholesky(signal+sampling).T
    # Heterogeneous sampling exercises the compiled numerical fit rather
    # than the exact common-covariance solution.
    scales = np.linspace(.8, 1.2, n)
    effects *= scales[:, None]
    data = tmp_path / 'effects.npz'
    np.savez(data, effects=effects, sampling=sampling[None] * scales[:, None, None]**2)
    sumstats = tmp_path / 'sumstats.tsv'
    with sumstats.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['SNP', 'CHR', 'pos', 'freq', 'direct_Beta', 'direct_SE',
                         'avg_NTC_Beta', 'avg_NTC_SE', 'r_direct_avg_NTC'])
        for index, effect in enumerate(effects):
            writer.writerow(['rs'+str(index), 1, index+1, .3, effect[0], scales[index]*np.sqrt(.1),
                             effect[1], scales[index]*np.sqrt(.15), -.03/np.sqrt(.015)])
    package_script = """import json,sys,numpy as np
from snipar.coherent import fit_coherent,fast_kernels
data=np.load(sys.argv[1])
result=fit_coherent(data['effects'],data['sampling'],n_blocks=2)
print(json.dumps({'metric_names':result['metric_names'],
                  'estimates':result['estimates'].tolist(),'se':result['se'].tolist(),
                  'compiled_module':fast_kernels._criterion_jit.py_func.__module__}))
"""
    results = {}
    order = [first, 'package' if first == 'cli' else 'cli']
    for index, mode in enumerate(order):
        environment = os.environ.copy()
        environment.update(PYTHONPATH=str(source_root),
                           NUMBA_CACHE_DIR=str(tmp_path / 'cache'), NUMBA_DEBUG_CACHE='1',
                           # Cache portability does not depend on OpenMP shared memory.
                           NUMBA_THREADING_LAYER='workqueue',
                           OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1',
                           NUMEXPR_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', NUMBA_NUM_THREADS='1')
        if mode == 'cli':
            prefix = tmp_path / 'cli_result'
            command = [sys.executable, str(scripts / 'correlate.py'), str(sumstats), str(prefix),
                       '--ld-weight', 'none', '--mode', 'average', '--n_blocks', '2', '--save_delete']
        else:
            command = [sys.executable, '-c', package_script, str(data)]
        process = subprocess.run(command, cwd=outside, env=environment, text=True,
                                 capture_output=True, timeout=120)
        assert process.returncode == 0, process.stdout + '\n' + process.stderr
        if index == 1:
            assert 'data loaded from' in process.stdout, process.stdout
        if mode == 'cli':
            with Path(str(prefix)+'_corrs.txt').open() as handle:
                records = list(csv.DictReader(handle, delimiter='\t'))
            results[mode] = {row['correlation']: row for row in records}
        else:
            results[mode] = json.loads(process.stdout.strip().splitlines()[-1])
            assert results[mode]['compiled_module'] == 'snipar.coherent_kernels'
    # The CLI additionally exports legacy derived rows. Compare the core
    # metrics by name so compatibility output ordering cannot mask differences.
    for field, column in [('estimates', 'est'), ('se', 'SE')]:
        cli_values = [float(results['cli'][name][column])
                      for name in results['package']['metric_names']]
        np.testing.assert_allclose(cli_values, results['package'][field],
                                   atol=1e-12, rtol=1e-10, equal_nan=True)


def test_missing_cache_locator_uses_uncached_jit(monkeypatch):
    from snipar import coherent_kernels
    from numba import njit

    def njit_without_writable_cache(**options):
        if options['cache']:
            def reject_cache(function):
                raise RuntimeError("cannot cache function: no locator available for file")
            return reject_cache
        return njit(**options)

    monkeypatch.setattr(coherent_kernels, '_numba_njit', njit_without_writable_cache)

    @coherent_kernels.njit(cache=True)
    def increment(value):
        return value + 1

    assert increment(2) == 3
    assert increment.signatures
