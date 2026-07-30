"""Utilities for multi-core parallel processing."""

from __future__ import annotations

import atexit
import gzip
import math
import os
import tempfile
from concurrent import futures
from contextlib import contextmanager, suppress
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator
    from concurrent.futures.process import ProcessPoolExecutor

#: Upper bound on BED lines per parallel chunk, to cap each task's peak memory.
MAX_CHUNK_SIZE = 5000
#: Chunks created per worker. More than one chunk per worker lets the pool
#: balance uneven per-region cost dynamically.
CHUNKS_PER_PROCESS = 2


class SerialPool:
    """Mimic the concurrent.futures.Executor interface, but run in serial."""

    def __init__(self) -> None:
        pass

    def submit(self, func: Callable, *args, **kwargs) -> SerialFuture:
        """Just call the function on the arguments."""
        try:
            result = func(*args, **kwargs)
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


def effective_procs(nprocs: int) -> int:
    """Resolve a requested process count to the number of workers to use.

    `nprocs` <= 0 means "use all available CPUs". The requested count is clamped
    to the number of usable CPUs, so asking for more workers than cores never
    oversubscribes.
    """
    avail = available_cpus()
    return avail if nprocs < 1 else min(nprocs, avail)


@contextmanager
def pick_pool(nprocs: int) -> Iterator[SerialPool | ProcessPoolExecutor]:
    """Yield a process pool, or a serial stand-in when parallelism is moot.

    `nprocs` is resolved by `effective_procs`: at most 1 means "use all
    available CPUs", and larger requests are clamped to the usable CPU count.
    When the effective worker count is 1 -- which includes every single-CPU
    host -- work runs in-process via `SerialPool` instead of forking workers:
    real multiprocessing offers no speedup there, and on constrained single-CPU
    build environments forking the pool can deadlock (#1103).
    """
    nprocs = effective_procs(nprocs)
    if nprocs <= 1:
        yield SerialPool()
    else:
        with futures.ProcessPoolExecutor(max_workers=nprocs) as pool:
            yield pool


def rm(path: str) -> None:
    """Safely remove a file."""
    with suppress(OSError):
        os.unlink(path)


def _bed_lines(bed_fname: str) -> Iterator[str]:
    """Yield the content lines of a possibly gzipped BED, skipping comments."""
    opener = gzip.open if bed_fname.endswith(".gz") else open
    with opener(bed_fname) as infile:
        for line in infile:
            if isinstance(line, bytes):
                line = line.decode()
            if line[0] == "#":
                continue
            yield line


def to_chunks(bed_fname: str, *, nprocs: int = 1) -> Iterator[str]:
    """Split a BED file into parts for parallel processing.

    Chunk size is derived from the number of workers that will consume the
    chunks, resolved from `nprocs` by `effective_procs`. Each worker is given
    `CHUNKS_PER_PROCESS` chunks so the pool can balance uneven per-region cost
    dynamically, and chunks hold at most `MAX_CHUNK_SIZE` lines. A single worker
    gets one chunk per `MAX_CHUNK_SIZE` lines, reproducing the fixed-size
    behavior this replaces: that fixed size serialized any BED smaller than one
    chunk, so `coverage -p N` ran single-threaded on small (e.g. WGS access)
    BEDs regardless of N (#526).

    Deriving the size from the worker count makes chunk boundaries a function of
    `nprocs`, so only callers whose per-chunk processing is independent of those
    boundaries may pass it (`cnvlib.cli.guess_baits` deliberately does not).

    Yields
    ------
    str
        Path to a temporary BED file holding the next chunk of regions, in the
        input file's order. Each is registered for removal at interpreter exit;
        callers should ``rm`` it as soon as the chunk has been processed.
    """
    nworkers = effective_procs(nprocs)
    # A lone worker gains nothing from extra chunks, just extra temp files.
    nchunks = nworkers * CHUNKS_PER_PROCESS if nworkers > 1 else 1
    if os.path.isfile(bed_fname):
        nlines = sum(1 for _ in _bed_lines(bed_fname))
        chunk_size = min(MAX_CHUNK_SIZE, max(1, math.ceil(nlines / nchunks)))
    else:
        # A stream (FIFO, process substitution) can only be read once, so skip
        # the counting pass rather than hanging on a second open().
        chunk_size = MAX_CHUNK_SIZE
    k, chunk = 0, 0
    name, outfile = "", None
    for line in _bed_lines(bed_fname):
        if outfile is None:
            fd, name = tempfile.mkstemp(suffix=".bed", prefix=f"tmp.{chunk}.")
            outfile = os.fdopen(fd, "w")
            atexit.register(rm, name)
        k += 1
        outfile.write(line)
        if k % chunk_size == 0:
            outfile.close()
            outfile = None
            chunk += 1
            yield name
    if outfile is not None:
        outfile.close()
        yield name
