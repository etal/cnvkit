"""Utilities for multi-core parallel processing."""
from __future__ import absolute_import, division, print_function
from builtins import object

import atexit
import tempfile
import gzip
import os
from contextlib import contextmanager
from concurrent import futures
#  from concurrent.futures import wait


class SerialPool(object):
    """Mimic the concurrent.futures.Executor interface, but run in serial."""

    def __init__(self):
        pass

    def submit(self, func, *args):
        """Just call the function on the arguments."""
        return SerialFuture(func(*args))

    def map(self, func, iterable):
        """Just apply the function to `iterable`."""
        return map(func, iterable)

    def shutdown(self, wait=True):
        """Do nothing."""
        pass



class SerialFuture(object):
    """Mimic the concurrent.futures.Future interface."""

    def __init__(self, result):
        self._result = result

    def result(self):
        return self._result



@contextmanager
def pick_pool(nprocs):
    if nprocs == 1:
        yield SerialPool()
    else:
        if nprocs < 1:
            nprocs = None
        with futures.ProcessPoolExecutor(max_workers=nprocs) as pool:
            yield pool


def rm(path):
    """Safely remove a file."""
    try:
        os.unlink(path)
    except OSError:
        pass


def to_chunks(bed_fname, chunk_size=5000):
    """Split a BED file into `chunk_size`-line parts for parallelization."""
    k, chunk = 0, 0
    fd, name = tempfile.mkstemp(suffix=".bed", prefix="tmp.%s." % chunk)
    outfile = os.fdopen(fd, "w")
    atexit.register(rm, name)
    opener = (gzip.open if bed_fname.endswith(".gz") else open)
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
                fd, name = tempfile.mkstemp(suffix=".bed",
                                            prefix="tmp.%s." % chunk)
                outfile = os.fdopen(fd, "w")
    outfile.close()
    if k % chunk_size:
        outfile.close()
        yield name
