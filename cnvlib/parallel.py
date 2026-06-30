"""Utilities for multi-core parallel processing."""

from __future__ import annotations

import atexit
import gzip
import os
import tempfile
from concurrent import futures
from contextlib import contextmanager, suppress
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator
    from concurrent.futures.process import ProcessPoolExecutor


class SerialPool:
    """Mimic the concurrent.futures.Executor interface, but run in serial."""

    def __init__(self) -> None:
        pass

    def submit(self, func: Callable, *args) -> SerialFuture:
        """Just call the function on the arguments."""
        try:
            result = func(*args)
            return SerialFuture(result=result)
        except Exception as exc:
            return SerialFuture(exception=exc)

    def map(self, func: Callable, iterable: Iterator[Any]) -> map:
        """Just apply the function to `iterable`."""
        return map(func, iterable)

    def shutdown(self, wait=True) -> None:
        """Do nothing."""


class SerialFuture:
    """Mimic the concurrent.futures.Future interface."""

    def __init__(self, result: Any = None, exception: Exception | None = None) -> None:
        self._result = result
        self._exception = exception

    def result(self) -> Any:
        if self._exception is not None:
            raise self._exception
        return self._result


def available_cpus() -> int:
    """Number of CPUs this process may actually use.

    Honors CPU affinity and cgroup limits (e.g. ``taskset``, container quotas,
    or the kernel ``nr_cpus=1`` boot parameter) rather than the machine's total
    core count, so a constrained single-CPU environment is reported as 1.
    """
    # Python 3.13+: respects affinity, cgroups, and -X cpu_count / PYTHON_CPU_COUNT.
    process_cpu_count = getattr(os, "process_cpu_count", None)
    if process_cpu_count is not None:
        n = process_cpu_count()
        if n:
            return int(n)
    # Linux/Unix: respects affinity set by taskset / sched_setaffinity.
    sched_getaffinity = getattr(os, "sched_getaffinity", None)
    if sched_getaffinity is not None:
        try:
            n = len(sched_getaffinity(0))
        except OSError:
            n = 0
        if n:
            return n
    return os.cpu_count() or 1


@contextmanager
def pick_pool(nprocs: int) -> Iterator[SerialPool | ProcessPoolExecutor]:
    """Yield a process pool, or a serial stand-in when parallelism is moot.

    `nprocs` <= 0 means "use all available CPUs". The requested count is clamped
    to the number of usable CPUs, so asking for more workers than cores never
    oversubscribes. When the effective count is 1 -- which includes every
    single-CPU host -- work runs in-process via `SerialPool` instead of forking
    workers: real multiprocessing offers no speedup there, and on constrained
    single-CPU build environments forking the pool can deadlock (#1103).
    """
    avail = available_cpus()
    nprocs = avail if nprocs < 1 else min(nprocs, avail)
    if nprocs <= 1:
        yield SerialPool()
    else:
        with futures.ProcessPoolExecutor(max_workers=nprocs) as pool:
            yield pool


def rm(path: str) -> None:
    """Safely remove a file."""
    with suppress(OSError):
        os.unlink(path)


def to_chunks(bed_fname: str, chunk_size: int = 5000) -> Iterator[str]:
    """Split a BED file into `chunk_size`-line parts for parallelization."""
    k, chunk = 0, 0
    fd, name = tempfile.mkstemp(suffix=".bed", prefix=f"tmp.{chunk}.")
    outfile = os.fdopen(fd, "w")
    atexit.register(rm, name)
    opener = gzip.open if bed_fname.endswith(".gz") else open
    with opener(bed_fname) as infile:
        for line in infile:
            if isinstance(line, bytes):
                line = line.decode()
            if line[0] == "#":
                continue
            k += 1
            outfile.write(line)
            if k % chunk_size == 0:
                outfile.close()
                yield name
                chunk += 1
                fd, name = tempfile.mkstemp(suffix=".bed", prefix=f"tmp.{chunk}.")
                outfile = os.fdopen(fd, "w")
    outfile.close()
    if k % chunk_size:
        outfile.close()
        yield name
