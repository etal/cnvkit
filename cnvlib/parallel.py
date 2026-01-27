"""Utilities for multi-core parallel processing."""

from __future__ import annotations
import atexit
import tempfile
import gzip
import os
from contextlib import contextmanager, suppress
from concurrent import futures
from typing import TYPE_CHECKING, Any, Union
from collections.abc import Callable

if TYPE_CHECKING:
    from collections.abc import Iterator
    from concurrent.futures.process import ProcessPoolExecutor


class SerialPool:
    """Mimic the concurrent.futures.Executor interface, but run in serial."""

    def __init__(self) -> None:
        pass

    def submit(self, func: Callable, *args) -> SerialFuture:
        """Just call the function on the arguments."""
        return SerialFuture(func(*args))

    def map(self, func: Callable, iterable: Iterator[Any]) -> map:
        """Just apply the function to `iterable`."""
        return map(func, iterable)

    def shutdown(self, wait=True) -> None:
        """Do nothing."""
        pass


class SerialFuture:
    """Mimic the concurrent.futures.Future interface."""

    def __init__(self, result: str) -> None:
        self._result = result

    def result(self) -> str:
        return self._result


@contextmanager
def pick_pool(nprocs: int) -> Iterator[Union[SerialPool, ProcessPoolExecutor]]:
    if nprocs == 1:
        yield SerialPool()
    else:
        if nprocs < 1:
            nprocs = None
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
