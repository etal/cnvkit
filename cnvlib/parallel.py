"""Utilities for multi-core parallel processing."""
from __future__ import absolute_import, division, print_function
import multiprocessing


class SerialPool(object):
    """Mimic the multiprocessing.Pool interface, but run in serial."""

    def __init__(self):
        pass

    def apply_async(self, func, args):
        """Just call the function."""
        func(*args)

    # No-ops to mimic multiprocessing.Pool
    def close(self): pass
    def join(self): pass


def pick_pool(nprocs):
    if nprocs == 1:
        return SerialPool()
    if nprocs < 1:
        nprocs = multiprocessing.cpu_count()
    return multiprocessing.Pool(nprocs)
