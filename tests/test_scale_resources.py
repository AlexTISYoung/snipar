from types import SimpleNamespace

import pytest

from snipar import resources


def test_batch_memory_bounds_include_fixed_and_concurrent_work(monkeypatch):
    monkeypatch.setattr(resources, 'available_cpus', lambda: 12)
    assert resources.plan_batches(1000, 10, 10, fixed_bytes_per_worker=50,
                                  fixed_bytes=100, memory_limit=601) == (8, 1)
    assert resources.plan_batches(1000, 10, 2, fixed_bytes_per_worker=50,
                                  fixed_bytes=100, memory_limit=1300) == (2, 55)
    assert resources.plan_batches(1000, 10, 2, 11, memory_limit=10000) == (2, 11)
    with pytest.raises(MemoryError, match='minimum workspace'):
        resources.plan_batches(10, 100, 2, memory_limit=99)


def test_container_and_scheduler_limits(tmp_path, monkeypatch):
    group = tmp_path / 'child'
    group.mkdir()
    for path, max_mem, usage, cpus in ((tmp_path, 3000, 2500, '200000 100000'),
                                       (group, 2000, 1000, 'max 100000')):
        (path / 'memory.max').write_text(str(max_mem))
        (path / 'memory.current').write_text(str(usage))
        (path / 'cpu.max').write_text(cpus)
    monkeypatch.setattr(resources, '_cgroup_dirs', lambda: iter((('', group), ('', tmp_path))))
    monkeypatch.setattr(resources.psutil, 'virtual_memory', lambda: SimpleNamespace(available=100_000))
    monkeypatch.setattr(resources.os, 'process_cpu_count', lambda: 16, raising=False)
    assert resources.get_memory_budget() == 300
    assert resources.available_cpus() == 2
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '1')
    assert resources.available_cpus() == 1


def test_spawn_environment_restored(monkeypatch):
    monkeypatch.setenv('OPENBLAS_NUM_THREADS', '8')
    monkeypatch.delenv('MKL_NUM_THREADS', raising=False)
    with resources.worker_environment(2):
        assert resources.os.environ['OPENBLAS_NUM_THREADS'] == '2'
        assert resources.os.environ['MKL_NUM_THREADS'] == '2'
    assert resources.os.environ['OPENBLAS_NUM_THREADS'] == '8'
    assert 'MKL_NUM_THREADS' not in resources.os.environ
