"""Unit tests for the parallel module."""

import glob
import gzip
import os
import tempfile
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


class ToChunksTests(unittest.TestCase):
    """Tests for splitting a BED file into per-worker chunks."""

    def _write_bed(self, nlines, *, comments=0, gzipped=False):
        """Write a BED of `nlines` regions, optionally gzipped, and return it."""
        lines = [f"#comment {i}\n" for i in range(comments)]
        lines += [f"chr1\t{i * 100}\t{i * 100 + 50}\tbin{i}\n" for i in range(nlines)]
        suffix = ".bed.gz" if gzipped else ".bed"
        fd, fname = tempfile.mkstemp(suffix=suffix)
        os.close(fd)
        self.addCleanup(parallel.rm, fname)
        opener = gzip.open if gzipped else open
        with opener(fname, "wt") as fh:
            fh.writelines(lines)
        return fname, lines[comments:]

    @staticmethod
    def _read_chunks(chunk_fnames):
        """Read back the lines of each chunk file, in chunk order."""
        chunks = []
        for fname in chunk_fnames:
            with open(fname) as fh:
                chunks.append(fh.readlines())
            parallel.rm(fname)
        return chunks

    def test_to_chunks_fills_every_worker(self):
        """A BED smaller than the maximum chunk size must still be split across
        the available workers; a fixed chunk size made `coverage -p N` run
        single-threaded on small (e.g. WGS access) BEDs regardless of N (#526).
        """
        bed, lines = self._write_bed(440)
        with mock.patch.object(parallel, "available_cpus", return_value=16):
            chunks = self._read_chunks(parallel.to_chunks(bed, nprocs=16))
        self.assertEqual(
            len(chunks),
            16 * parallel.CHUNKS_PER_PROCESS,
            "a 440-region BED under 16 workers must be split CHUNKS_PER_PROCESS "
            "ways per worker, not into a single serialized chunk",
        )
        self.assertEqual(
            [line for chunk in chunks for line in chunk],
            lines,
            "chunking must preserve every region, in the BED's original order",
        )
        self.assertTrue(all(chunks), "no chunk may be empty")

    def test_to_chunks_serial_yields_one_chunk(self):
        """With one worker, a BED below the size cap stays a single chunk."""
        bed, lines = self._write_bed(12)
        chunks = self._read_chunks(parallel.to_chunks(bed))
        self.assertEqual(len(chunks), 1)
        self.assertEqual(chunks[0], lines)

    def test_to_chunks_respects_max_chunk_size(self):
        """The per-worker size is capped, so a large BED under few workers is
        still split into bounded chunks rather than one huge task."""
        bed, lines = self._write_bed(10)
        with mock.patch.object(parallel, "MAX_CHUNK_SIZE", 3):
            chunks = self._read_chunks(parallel.to_chunks(bed))
        self.assertEqual([len(chunk) for chunk in chunks], [3, 3, 3, 1])
        self.assertEqual([line for chunk in chunks for line in chunk], lines)

    def test_to_chunks_exact_multiple_yields_no_empty_chunk(self):
        """When the line count divides evenly by the chunk size, the last chunk
        is not an empty file -- and no unyielded temp file is left behind."""
        bed, lines = self._write_bed(9)
        before = set(glob.glob(os.path.join(tempfile.gettempdir(), "tmp.*.bed")))
        with mock.patch.object(parallel, "MAX_CHUNK_SIZE", 3):
            chunks = self._read_chunks(parallel.to_chunks(bed))
        after = set(glob.glob(os.path.join(tempfile.gettempdir(), "tmp.*.bed")))
        self.assertEqual([len(chunk) for chunk in chunks], [3, 3, 3])
        self.assertEqual([line for chunk in chunks for line in chunk], lines)
        self.assertEqual(after - before, set(), "left an unyielded chunk file")

    def test_to_chunks_empty_bed_yields_nothing(self):
        """A BED with no regions yields no chunks and creates no temp file, so
        callers see an empty iteration rather than an empty chunk."""
        bed, _lines = self._write_bed(0, comments=2)
        before = set(glob.glob(os.path.join(tempfile.gettempdir(), "tmp.*.bed")))
        with mock.patch.object(parallel, "available_cpus", return_value=4):
            chunks = self._read_chunks(parallel.to_chunks(bed, nprocs=4))
        after = set(glob.glob(os.path.join(tempfile.gettempdir(), "tmp.*.bed")))
        self.assertEqual(chunks, [])
        self.assertEqual(after - before, set())

    def test_to_chunks_gzipped_input_and_comments(self):
        """Comment lines are dropped and don't count toward chunk sizing, for
        gzipped input as well as plain text."""
        for gzipped in (False, True):
            with self.subTest(gzipped=gzipped):
                bed, lines = self._write_bed(8, comments=3, gzipped=gzipped)
                with mock.patch.object(parallel, "available_cpus", return_value=4):
                    chunks = self._read_chunks(parallel.to_chunks(bed, nprocs=4))
                flat = [line for chunk in chunks for line in chunk]
                self.assertEqual(flat, lines)
                self.assertFalse(any(line.startswith("#") for line in flat))
                self.assertEqual([len(chunk) for chunk in chunks], [1] * 8)

    def test_to_chunks_rejects_positional_worker_count(self):
        """`nprocs` is keyword-only: the parameter used to be a line count, so a
        stale positional call must fail loudly rather than request 5000 workers.
        """
        bed, _lines = self._write_bed(4)
        with self.assertRaises(TypeError):
            list(parallel.to_chunks(bed, 5000))


if __name__ == "__main__":
    unittest.main(verbosity=2)
