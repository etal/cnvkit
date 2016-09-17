"""Utilities for multi-core parallel processing."""
from __future__ import absolute_import, division, print_function
from builtins import object
from contextlib import contextmanager

from concurrent import futures
from concurrent.futures import wait


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
        raise StopIteration

    if nprocs < 1:
        nprocs = None
    with futures.ProcessPoolExecutor(max_workers=nprocs) as pool:
        yield pool
