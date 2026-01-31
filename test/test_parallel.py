"""Unit tests for the parallel module."""

import unittest

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
        self.assertIsNone(future._exception, "SerialFuture should not have exception for successful execution")

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

    def test_process_pool_exception_propagation(self):
        """Test that ProcessPoolExecutor properly propagates exceptions from workers."""
        from concurrent import futures

        # Test that exceptions from worker processes are propagated
        with futures.ProcessPoolExecutor(max_workers=2) as pool:
            future = pool.submit(_failing_worker)

            with self.assertRaises(ValueError) as ctx:
                future.result()
            self.assertEqual(str(ctx.exception), "Worker process error",
                           "ProcessPoolExecutor should propagate ValueError from worker process")

    def test_pick_pool_with_multiple_processes(self):
        """Test pick_pool correctly creates ProcessPoolExecutor for multiple processes."""
        # Test with 2 processes
        with parallel.pick_pool(2) as pool:
            future = pool.submit(_simple_task, 21)
            result = future.result()
            self.assertEqual(result, 42, "ProcessPoolExecutor should execute _simple_task(21) = 42")

    def test_pick_pool_with_single_process(self):
        """Test pick_pool correctly creates SerialPool for single process."""

        def simple_task(x):
            return x * 2

        # Test with 1 process (should use SerialPool)
        with parallel.pick_pool(1) as pool:
            self.assertIsInstance(pool, parallel.SerialPool,
                                "pick_pool(1) should return SerialPool for single-process execution")
            future = pool.submit(simple_task, 21)
            result = future.result()
            self.assertEqual(result, 42, "SerialPool should execute simple_task(21) = 42")


class CNAryByGeneTests(unittest.TestCase):
    """Tests for the by_gene method with various index types."""

    def test_by_gene_with_numeric_index(self):
        """Test by_gene with standard numeric index (default)."""
        # Create test data with numeric index
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 6,
                "start": [100, 200, 300, 400, 500, 600],
                "end": [150, 250, 350, 450, 550, 650],
                "gene": ["GeneA", "GeneA", "Antitarget", "GeneB", "GeneB", "GeneB"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2, -0.3],
            }
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        # Group by gene and check for duplicates
        gene_groups = list(cnarr.by_gene())
        all_coords = []

        for gene_name, gene_cnarr in gene_groups:
            for row in gene_cnarr.data.itertuples():
                coord = (row.chromosome, row.start, row.end)
                self.assertNotIn(
                    coord,
                    all_coords,
                    f"Duplicate coordinate found in {gene_name}: {coord}",
                )
                all_coords.append(coord)

        # Should have processed all 6 rows
        self.assertEqual(len(all_coords), 6)

    def test_by_gene_with_non_sequential_index(self):
        """Test by_gene with non-sequential numeric index."""
        # Create test data with gaps in the index
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "start": [100, 200, 300, 400, 500],
                "end": [150, 250, 350, 450, 550],
                "gene": ["GeneA", "GeneA", "Antitarget", "GeneB", "GeneB"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2],
            },
            index=[0, 5, 10, 15, 20],  # Non-sequential index
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        # Group by gene and check for duplicates
        gene_groups = list(cnarr.by_gene())
        all_coords = []

        for gene_name, gene_cnarr in gene_groups:
            for row in gene_cnarr.data.itertuples():
                coord = (row.chromosome, row.start, row.end)
                self.assertNotIn(
                    coord,
                    all_coords,
                    f"Duplicate coordinate found in {gene_name}: {coord}",
                )
                all_coords.append(coord)

        # Should have processed all 5 rows
        self.assertEqual(len(all_coords), 5)

    def test_by_gene_preserves_boundaries(self):
        """Test that by_gene correctly handles gene boundaries without overlap."""
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "start": [100, 200, 300, 400, 500],
                "end": [150, 250, 350, 450, 550],
                "gene": ["GeneA", "GeneA", "Antitarget", "GeneB", "GeneB"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2],
            }
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        gene_groups = list(cnarr.by_gene())

        # Check that we have the expected groups
        gene_names = [name for name, _ in gene_groups]
        self.assertIn("GeneA", gene_names)
        self.assertIn("GeneB", gene_names)
        self.assertIn("Antitarget", gene_names)

        # Check GeneA has exactly 2 rows
        gene_a = next(ga for name, ga in gene_groups if name == "GeneA")
        self.assertEqual(len(gene_a), 2)
        self.assertEqual(gene_a.data.iloc[0]["start"], 100)
        self.assertEqual(gene_a.data.iloc[1]["end"], 250)

        # Check GeneB has exactly 2 rows
        gene_b = next(ga for name, ga in gene_groups if name == "GeneB")
        self.assertEqual(len(gene_b), 2)
        self.assertEqual(gene_b.data.iloc[0]["start"], 400)
        self.assertEqual(gene_b.data.iloc[1]["end"], 550)

    def test_by_gene_with_duplicate_index_labels(self):
        """Test by_gene with duplicate index labels (edge case for get_loc)."""
        # Create test data with duplicate index labels
        # This simulates the case where get_loc() returns a slice or boolean array
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 6,
                "start": [100, 200, 300, 400, 500, 600],
                "end": [150, 250, 350, 450, 550, 650],
                "gene": ["GeneA", "GeneA", "Antitarget", "GeneB", "GeneB", "GeneB"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2, -0.3],
            },
            index=[0, 1, 1, 2, 3, 3],  # Duplicate labels at positions 1-2 and 5-6
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        # Group by gene and verify no duplicates
        gene_groups = list(cnarr.by_gene())
        all_coords = []

        for gene_name, gene_cnarr in gene_groups:
            for row in gene_cnarr.data.itertuples():
                coord = (row.chromosome, row.start, row.end)
                self.assertNotIn(
                    coord,
                    all_coords,
                    f"Duplicate coordinate found in {gene_name}: {coord}",
                )
                all_coords.append(coord)

        # Should have processed all 6 rows
        self.assertEqual(len(all_coords), 6)

    def test_by_gene_multiple_genes_per_bin(self):
        """Test by_gene with bins spanning multiple genes (comma-separated)."""
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "start": [100, 200, 300, 400, 500],
                "end": [150, 250, 350, 450, 550],
                "gene": ["GeneA", "GeneA,GeneB", "GeneB", "Antitarget", "GeneC"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2],
            }
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        # Group by gene
        gene_groups = list(cnarr.by_gene())
        gene_dict = {name: ga for name, ga in gene_groups}

        # The bin at position 1 with "GeneA,GeneB" should appear in both genes
        self.assertIn("GeneA", gene_dict)
        self.assertIn("GeneB", gene_dict)
        self.assertIn("GeneC", gene_dict)

        # GeneA should include the shared bin
        gene_a_coords = [
            (row.chromosome, row.start, row.end)
            for row in gene_dict["GeneA"].data.itertuples()
        ]
        self.assertIn(("chr1", 200, 250), gene_a_coords)

        # GeneB should also include the shared bin
        gene_b_coords = [
            (row.chromosome, row.start, row.end)
            for row in gene_dict["GeneB"].data.itertuples()
        ]
        self.assertIn(("chr1", 200, 250), gene_b_coords)


if __name__ == "__main__":
    unittest.main()
