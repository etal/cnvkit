"""Unit tests for the parallel module."""

import unittest
from concurrent import futures
from unittest import mock

import numpy as np
import pandas as pd

from cnvlib import parallel
from cnvlib.cnary import CopyNumArray


# Module-level functions for multiprocessing tests
# (must be at module level to be picklable)
def _failing_worker():
    """Worker function that raises an exception."""
    raise ValueError("Worker process error")


def _simple_task(x):
    """Simple worker function for testing."""
    return x * 2


class ParallelTests(unittest.TestCase):
    """Tests for the parallel module."""

    def test_serial_pool_exception_propagation(self):
        """Test that SerialPool properly stores and re-raises exceptions."""

        def failing_func():
            raise ValueError("Test error message")

        pool = parallel.SerialPool()
        future = pool.submit(failing_func)

        # The exception should be stored in the future
        self.assertIsNotNone(future._exception)
        self.assertIsInstance(future._exception, ValueError)

        # The exception should be raised when calling result()
        with self.assertRaises(ValueError) as ctx:
            future.result()
        self.assertEqual(str(ctx.exception), "Test error message")

    def test_serial_pool_success(self):
        """Test that SerialPool correctly handles successful execution."""

        def success_func(x, y):
            return x + y

        pool = parallel.SerialPool()
        future = pool.submit(success_func, 2, 3)

        # No exception should be stored
        self.assertIsNone(
            future._exception,
            "SerialFuture should not have exception for successful execution",
        )

        # Result should be available
        result = future.result()
        self.assertEqual(result, 5, "SerialPool should execute success_func(2, 3) = 5")

    def test_serial_future_creation_with_result(self):
        """Test SerialFuture creation with a result value."""
        future = parallel.SerialFuture(result=42)
        self.assertEqual(future.result(), 42)
        self.assertIsNone(future._exception)

    def test_serial_future_creation_with_exception(self):
        """Test SerialFuture creation with an exception."""
        exc = RuntimeError("test exception")
        future = parallel.SerialFuture(exception=exc)

        with self.assertRaises(RuntimeError) as ctx:
            future.result()
        self.assertEqual(str(ctx.exception), "test exception")

    @unittest.skipIf(
        parallel.available_cpus() < 2,
        "real multiprocessing requires more than one usable CPU (#1103)",
    )
    def test_process_pool_exception_propagation(self):
        """Test that ProcessPoolExecutor properly propagates exceptions from workers."""
        # Test that exceptions from worker processes are propagated
        with futures.ProcessPoolExecutor(max_workers=2) as pool:
            future = pool.submit(_failing_worker)

            with self.assertRaises(ValueError) as ctx:
                future.result()
            self.assertEqual(
                str(ctx.exception),
                "Worker process error",
                "ProcessPoolExecutor should propagate ValueError from worker process",
            )

    def test_pick_pool_with_multiple_processes(self):
        """pick_pool(2) yields a real ProcessPoolExecutor when more than one CPU
        is usable, and a SerialPool on a single-CPU host (#1103)."""
        expected = (
            futures.ProcessPoolExecutor
            if parallel.available_cpus() > 1
            else parallel.SerialPool
        )
        with parallel.pick_pool(2) as pool:
            self.assertIsInstance(pool, expected)
            self.assertEqual(pool.submit(_simple_task, 21).result(), 42)

    def test_pick_pool_with_single_process(self):
        """Test pick_pool correctly creates SerialPool for single process."""

        def simple_task(x):
            return x * 2

        # Test with 1 process (should use SerialPool)
        with parallel.pick_pool(1) as pool:
            self.assertIsInstance(
                pool,
                parallel.SerialPool,
                "pick_pool(1) should return SerialPool for single-process execution",
            )
            future = pool.submit(simple_task, 21)
            result = future.result()
            self.assertEqual(
                result, 42, "SerialPool should execute simple_task(21) = 42"
            )

    def test_pick_pool_single_cpu_never_forks(self):
        """On a single-CPU host, pick_pool must run in-process rather than fork
        workers, however many are requested (regression test for #1103)."""
        with mock.patch.object(parallel, "available_cpus", return_value=1):
            for nprocs in (2, 8, 0, -1):
                with parallel.pick_pool(nprocs) as pool:
                    self.assertIsInstance(
                        pool,
                        parallel.SerialPool,
                        f"pick_pool({nprocs}) must be serial on a single CPU",
                    )
                    self.assertEqual(pool.submit(_simple_task, 21).result(), 42)

    def test_pick_pool_clamps_to_available_cpus(self):
        """pick_pool never starts more workers than there are usable CPUs."""
        with mock.patch.object(parallel, "available_cpus", return_value=4):
            with parallel.pick_pool(64) as pool:
                self.assertIsInstance(pool, futures.ProcessPoolExecutor)
                self.assertEqual(pool._max_workers, 4)
