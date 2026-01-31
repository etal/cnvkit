#!/usr/bin/env python
"""Unit tests for RNA import functionality (cnvlib.rna)."""
import logging
import unittest

import numpy as np
import pandas as pd

from cnvlib import rna

logging.basicConfig(level=logging.ERROR, format="%(message)s")


class RNAImportTests(unittest.TestCase):
    """Tests for RNA import and processing functions."""

    def setUp(self):
        """Create sample test data."""
        # Sample gene counts DataFrame (genes x samples)
        self.sample_counts = pd.DataFrame(
            {
                "sample1": [100, 200, 50, 0, 150],
                "sample2": [120, 180, 45, 0, 140],
                "sample3": [90, 210, 55, 0, 160],
            },
            index=["ENSG00000001", "ENSG00000002", "ENSG00000003", "ENSG00000004", "ENSG00000005"],
        )

        # Gene info DataFrame
        self.gene_info = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1", "chr2", "chr2", "chr3"],
                "start": [1000, 5000, 10000, 15000, 20000],
                "end": [2000, 6000, 11000, 16000, 21000],
                "gene": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
                "gc": [0.5, 0.6, 0.4, 0.55, 0.45],
                "tx_length": [1000, 1200, 800, 1100, 900],
            },
            index=["ENSG00000001", "ENSG00000002", "ENSG00000003", "ENSG00000004", "ENSG00000005"],
        )

        # Transcript lengths Series
        self.tx_lengths = pd.Series(
            [1000, 1200, 800, 1100, 900],
            index=["ENSG00000001", "ENSG00000002", "ENSG00000003", "ENSG00000004", "ENSG00000005"],
        )

    def test_align_gene_info_empty_intersection(self):
        """Test error when no genes match between sample data and gene resource."""
        # Create gene_info with completely different gene IDs
        mismatched_gene_info = self.gene_info.copy()
        mismatched_gene_info.index = ["ENSG99990001", "ENSG99990002", "ENSG99990003", "ENSG99990004", "ENSG99990005"]

        with self.assertRaises(ValueError) as cm:
            rna.align_gene_info_to_samples(
                mismatched_gene_info, self.sample_counts, self.tx_lengths, []
            )

        error_msg = str(cm.exception)
        self.assertIn("No genes in common", error_msg)
        self.assertIn("Sample data has 5 genes", error_msg)
        self.assertIn("gene resource has 5 genes", error_msg)
        self.assertIn("ENSG00000001", error_msg)  # Sample gene ID
        self.assertIn("ENSG99990001", error_msg)  # Gene resource ID

    def test_align_gene_info_invalid_tx_lengths(self):
        """Test filtering of genes with invalid transcript lengths."""
        # Create gene_info with some invalid tx_lengths
        # Pass tx_lengths=None so gene_info values are used
        invalid_gene_info = self.gene_info.copy()
        invalid_gene_info.loc["ENSG00000003", "tx_length"] = 0
        invalid_gene_info.loc["ENSG00000004", "tx_length"] = -100

        result_gi, result_sc, result_log2 = rna.align_gene_info_to_samples(
            invalid_gene_info, self.sample_counts, tx_lengths=None, normal_ids=[]
        )

        # Should have filtered out 2 genes with invalid tx_length
        self.assertEqual(len(result_gi), 3)
        self.assertNotIn("ENSG00000003", result_gi.index)
        self.assertNotIn("ENSG00000004", result_gi.index)
        # Valid genes should remain
        self.assertIn("ENSG00000001", result_gi.index)
        self.assertIn("ENSG00000002", result_gi.index)
        self.assertIn("ENSG00000005", result_gi.index)

    def test_align_gene_info_all_invalid_tx_lengths(self):
        """Test error when all genes have invalid transcript lengths."""
        # Create gene_info with all invalid tx_lengths
        # Pass tx_lengths=None so gene_info values are used
        invalid_gene_info = self.gene_info.copy()
        invalid_gene_info["tx_length"] = 0

        with self.assertRaises(ValueError) as cm:
            rna.align_gene_info_to_samples(
                invalid_gene_info, self.sample_counts, tx_lengths=None, normal_ids=[]
            )

        error_msg = str(cm.exception)
        self.assertIn("All genes have invalid transcript lengths", error_msg)

    def test_normalize_read_depths_empty(self):
        """Test error when sample depths are empty."""
        empty_depths = pd.DataFrame()

        with self.assertRaises(ValueError) as cm:
            rna.normalize_read_depths(empty_depths, [])

        error_msg = str(cm.exception)
        self.assertIn("Sample read depths sum to zero or less", error_msg)
        self.assertIn("Input shape:", error_msg)

    def test_normalize_read_depths_all_zeros(self):
        """Test error when all sample depths are zero."""
        zero_depths = pd.DataFrame(
            {
                "sample1": [0, 0, 0],
                "sample2": [0, 0, 0],
            },
            index=["gene1", "gene2", "gene3"],
        )

        with self.assertRaises(ValueError) as cm:
            rna.normalize_read_depths(zero_depths, [])

        error_msg = str(cm.exception)
        self.assertIn("Sample read depths sum to zero or less", error_msg)
        self.assertIn("all values zero: True", error_msg)

    def test_normalize_read_depths_valid(self):
        """Test successful normalization with valid data."""
        valid_depths = pd.DataFrame(
            {
                "sample1": [10.0, 20.0, 30.0],
                "sample2": [15.0, 25.0, 35.0],
            },
            index=["gene1", "gene2", "gene3"],
        )

        result = rna.normalize_read_depths(valid_depths, [])

        # Should return normalized log2 values
        self.assertEqual(result.shape, valid_depths.shape)
        self.assertFalse(result.isna().any().any())
        # Values should be in reasonable log2 range
        self.assertTrue((result > -10).all().all())
        self.assertTrue((result < 10).all().all())

    def test_attach_gene_info_nan_handling(self):
        """Test that NaN values in gene_minima are properly handled."""
        # Create sample data with all NaN for one gene
        sample_data_log2 = pd.DataFrame(
            {
                "sample1": [1.0, np.nan, 2.0],
                "sample2": [1.5, np.nan, 2.5],
                "sample3": [0.5, np.nan, 1.5],
            },
            index=["gene1", "gene2", "gene3"],
        )

        sample_counts = pd.DataFrame(
            {
                "sample1": [100, 0, 200],
                "sample2": [120, 0, 180],
                "sample3": [80, 0, 220],
            },
            index=["gene1", "gene2", "gene3"],
        )

        gene_info = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1", "chr2"],
                "start": [1000, 5000, 10000],
                "end": [2000, 6000, 11000],
                "gene": ["GENE1", "GENE2", "GENE3"],
                "gc": [0.5, 0.6, 0.4],
                "tx_length": [1000, 1200, 800],
                "weight": [1.0, 1.0, 1.0],
            },
            index=["gene1", "gene2", "gene3"],
        )

        # Should not raise an assertion error
        cnrs = list(rna.attach_gene_info_to_cnr(sample_counts, sample_data_log2, gene_info))

        # Should have one CNR per sample
        self.assertEqual(len(cnrs), 3)

        # Check that gene2 (all NaN) has been filled with NULL_LOG2_COVERAGE
        for cnr in cnrs:
            gene2_log2 = cnr.data[cnr.data["gene"] == "GENE2"]["log2"].values[0]
            self.assertEqual(gene2_log2, rna.NULL_LOG2_COVERAGE)

    def test_safe_log2_zero_handling(self):
        """Test that safe_log2 handles zeros correctly."""
        values = np.array([0, 1, 10, 100])
        min_log2 = -5

        result = rna.safe_log2(values, min_log2)

        # Zero should be converted to approximately min_log2
        self.assertAlmostEqual(result[0], min_log2, places=2)
        # Other values should be reasonable
        self.assertTrue(result[1] > min_log2)
        self.assertTrue(result[2] > result[1])
        self.assertTrue(result[3] > result[2])


if __name__ == "__main__":
    unittest.main()
