"""Plan numerical work against the CPU and memory available to this process."""

from contextlib import contextmanager
import os
from pathlib import Path, PurePosixPath
import sys

import psutil


MIB = 1024 ** 2


def _read_file(path):
    try:
        return path.read_text().strip()
    except OSError:
        return ''


def _cgroup_dirs():
    """Yield each applicable Linux controller and its mounted ancestors."""
    if sys.platform != 'linux':
        return
    groups = {}
    for line in _read_file(Path('/proc/self/cgroup')).splitlines():
        fields = line.split(':', maxsplit=2)
        if len(fields) == 3:
            groups.update((controller, PurePosixPath(fields[2])) for controller in fields[1].split(','))
    for line in _read_file(Path('/proc/self/mountinfo')).splitlines():
        before, separator, after = line.partition(' - ')
        fields, fs = before.split(), after.split()
        if not separator or len(fields) < 5 or len(fs) < 3:
            continue
        if fs[0] == 'cgroup2':
            controllers = ('',)
        elif fs[0] == 'cgroup':
            controllers = tuple(set(fs[2].split(',')) & {'memory', 'cpu'})
        else:
            continue
        root = PurePosixPath(fields[3])
        mount = Path(fields[4].replace(r'\040', ' '))
        for controller in controllers:
            group = groups.get(controller)
            if group is None:
                continue
            # In a cgroup namespace the process can be at the mounted root.
            relative = group.relative_to(root) if group.is_relative_to(root) else PurePosixPath('.')
            location = mount / relative
            while True:
                yield controller, location
                if location == mount:
                    break
                location = location.parent


def _positive_environment(name):
    try:
        return max(0, int(os.environ.get(name, '0')))
    except ValueError:
        return 0


def available_cpus():
    """Return the logical CPU allowance, including affinity and quotas."""
    if hasattr(os, 'process_cpu_count'):
        count = os.process_cpu_count() or 1
    else:
        try:
            count = len(os.sched_getaffinity(0))
        except (AttributeError, OSError):
            count = os.cpu_count() or 1
    slurm_cpus = _positive_environment('SLURM_CPUS_PER_TASK')
    if slurm_cpus:
        count = min(count, slurm_cpus)
    for controller, path in _cgroup_dirs():
        if controller == '':
            limit = _read_file(path / 'cpu.max').split()
        elif controller == 'cpu':
            limit = [_read_file(path / 'cpu.cfs_quota_us'), _read_file(path / 'cpu.cfs_period_us')]
        else:
            continue
        if len(limit) == 2 and all(value.isdecimal() and int(value) > 0 for value in limit):
            count = min(count, (int(limit[0]) + int(limit[1]) - 1) // int(limit[1]))
    return max(1, count)


def get_memory_budget():
    """Estimate memory available for additional work, leaving room for the OS and allocator."""
    available = psutil.virtual_memory().available
    for controller, path in _cgroup_dirs():
        if controller == '':
            limits = (_read_file(path / 'memory.max'), _read_file(path / 'memory.high'))
            usage = _read_file(path / 'memory.current')
        elif controller == 'memory':
            limits = (_read_file(path / 'memory.limit_in_bytes'),)
            usage = _read_file(path / 'memory.usage_in_bytes')
        else:
            continue
        if usage.isdecimal():
            for limit in limits:
                if limit.isdecimal():
                    available = min(available, max(0, int(limit) - int(usage)))
    # Slurm may expose the job request even when the system has no memory cgroup.
    node = _positive_environment('SLURM_MEM_PER_NODE')
    per_cpu = _positive_environment('SLURM_MEM_PER_CPU')
    if node or per_cpu:
        allocated = (node or (per_cpu * available_cpus())) * MIB
        available = min(available, max(0, allocated - psutil.Process().memory_info().rss))
    return available * 3 // 5


def plan_batches(n_items, bytes_per_item, workers=1, batch_size=None, fixed_bytes_per_worker=0,
                 fixed_bytes=0, *, memory_limit=None):
    """Return worker count and maximum chunk length for the estimated peak memory.

    Fixed costs cover new allocations, independent of chunk length. The default
    memory limit reflects usage measured at the time the plan is made. Callers
    that launch several independent tasks can pass an allocated share explicitly.
    """
    if (n_items < 1 or bytes_per_item < 1 or workers < 1 or fixed_bytes_per_worker < 0
            or fixed_bytes < 0 or (batch_size is not None and batch_size < 1)):
        raise ValueError('Batch lengths, item sizes and worker counts must be positive')
    available = (get_memory_budget() if memory_limit is None else memory_limit) - fixed_bytes
    per_worker = fixed_bytes_per_worker + bytes_per_item
    if available < per_worker:
        raise MemoryError(f'Estimated minimum workspace is {per_worker + fixed_bytes:,} bytes, '
                          f'with {max(0, available + fixed_bytes):,} bytes available')
    workers = min(workers, n_items, available_cpus(), available // per_worker)
    items = (available // workers - fixed_bytes_per_worker) // bytes_per_item
    balanced = (n_items + workers - 1) // workers
    return int(workers), int(min(items, balanced, n_items if batch_size is None else batch_size))


@contextmanager
def worker_environment(threads):
    """Cap thread pools when newly spawned workers import their numeric libraries."""
    keys = ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'OMP_NUM_THREADS')
    previous = {key: os.environ.get(key) for key in keys}
    try:
        for key in keys:
            configured = _positive_environment(key)
            os.environ[key] = str(min(configured, threads) if configured else threads)
        yield
    finally:
        for key, old in previous.items():
            if old is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = old


def limit_native_threads(threads):
    """Limit BLAS in the current process while leaving explicit OpenMP teams alone."""
    from threadpoolctl import threadpool_limits

    return threadpool_limits(limits=threads, user_api='blas')
