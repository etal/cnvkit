#!/usr/bin/env python
"""Unit tests for RNA import functionality (cnvlib.rna)."""

import logging
import os
import tempfile
import unittest
import warnings

import numpy as np
import pandas as pd

from cnvlib import import_rna, rna

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
            index=[
                "ENSG00000001",
                "ENSG00000002",
                "ENSG00000003",
                "ENSG00000004",
                "ENSG00000005",
            ],
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
            index=[
                "ENSG00000001",
                "ENSG00000002",
                "ENSG00000003",
                "ENSG00000004",
                "ENSG00000005",
            ],
        )

        # Transcript lengths Series
        self.tx_lengths = pd.Series(
            [1000, 1200, 800, 1100, 900],
            index=[
                "ENSG00000001",
                "ENSG00000002",
                "ENSG00000003",
                "ENSG00000004",
                "ENSG00000005",
            ],
        )

    def test_align_gene_info_empty_intersection(self):
        """Test error when no genes match between sample data and gene resource."""
        # Create gene_info with completely different gene IDs
        mismatched_gene_info = self.gene_info.copy()
        mismatched_gene_info.index = [
            "ENSG99990001",
            "ENSG99990002",
            "ENSG99990003",
            "ENSG99990004",
            "ENSG99990005",
        ]

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

        result_gi, _result_sc, _result_log2 = rna.align_gene_info_to_samples(
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

    def test_align_gene_info_zero_spread_cohort(self):
        """Degenerate cohort (uniform counts across samples) yields uniform weights without a RuntimeWarning.

        Regression for the divide-by-zero at ``weight / weight.max()``: when every
        sample has identical counts, ``gene_spreads = std(axis=1)`` collapses to
        zero, so ``gmean(weights)`` is all-zero and ``weight.max() == 0`` would
        produce ``0 / 0`` → NaN with ``RuntimeWarning: invalid value encountered
        in divide``. The guard must short-circuit to uniform 1.0 weights.
        """
        uniform_counts = pd.DataFrame(
            {
                "sample1": [100, 200, 50, 150, 75],
                "sample2": [100, 200, 50, 150, 75],
                "sample3": [100, 200, 50, 150, 75],
            },
            index=[
                "ENSG00000001",
                "ENSG00000002",
                "ENSG00000003",
                "ENSG00000004",
                "ENSG00000005",
            ],
        )

        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            result_gi, _result_sc, _result_log2 = rna.align_gene_info_to_samples(
                self.gene_info, uniform_counts, self.tx_lengths, normal_ids=[]
            )

        # No NaN weights, and the degenerate case falls back to uniform 1.0.
        self.assertFalse(result_gi["weight"].isna().any())
        self.assertTrue(np.isfinite(result_gi["weight"]).all())
        self.assertTrue((result_gi["weight"] == 1.0).all())

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
        cnrs = list(
            rna.attach_gene_info_to_cnr(sample_counts, sample_data_log2, gene_info)
        )

        # Should have one CNR per sample
        self.assertEqual(len(cnrs), 3)

        # Check that gene2 (all NaN) has been filled with NULL_LOG2_COVERAGE
        for cnr in cnrs:
            gene2_log2 = cnr.data[cnr.data["gene"] == "GENE2"]["log2"].to_numpy()[0]
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


class ImportRnaIntegrationTests(unittest.TestCase):
    """End-to-end `import-rna` from count files + gene resource to .cnr."""

    COUNT_FILES = (
        "formats/rna-sample-A.counts.txt",
        "formats/rna-sample-B.counts.txt",
        "formats/rna-sample-C.counts.txt",
    )
    GENE_RESOURCE = "formats/rna-gene-resource.tsv"

    def test_do_import_rna_counts_to_cnr(self):
        """do_import_rna turns per-gene counts into one finite-log2 .cnr/sample."""
        all_data, cnrs = import_rna.do_import_rna(
            self.COUNT_FILES, "counts", self.GENE_RESOURCE
        )
        cnrs = list(cnrs)
        self.assertEqual(len(cnrs), len(self.COUNT_FILES))
        for cnr in cnrs:
            # Valid CopyNumArray: expected columns, finite log2, valid coords
            for col in ("chromosome", "start", "end", "gene", "log2", "depth"):
                self.assertIn(col, cnr.data.columns)
            self.assertGreater(len(cnr), 0)
            self.assertTrue(np.isfinite(cnr["log2"]).all())
            self.assertTrue((cnr.start < cnr.end).all())
            self.assertEqual(list(cnr.chromosome.unique()), ["1"])
        # Summary table has one row per retained gene
        self.assertEqual(len(all_data), len(cnrs[0]))
        self.assertEqual(cnrs[0].sample_id, "rna-sample-A")

    def test_do_import_rna_unknown_format_raises(self):
        with self.assertRaises(RuntimeError):
            import_rna.do_import_rna(self.COUNT_FILES[:1], "bogus", self.GENE_RESOURCE)


class NormalizeReadDepthsNormalAnchorTests(unittest.TestCase):
    """Pin the current --normal anchoring invariants surfaced in cnvkit-ccy.

    These tests document what ``normalize_read_depths`` does *today* with
    ``normal_ids`` set, as a regression baseline for any future redesign
    (tracked in cnvkit-9wj). They exercise the suspected bug from GH #352
    and confirm it does not reproduce at the unit level: the issue lies in
    the design (small-normal-cohort information loss) rather than a numerical
    error in the implementation.
    """

    @staticmethod
    def _make_cohort(
        n_tumor, n_normal, n_genes, altered_mask, tumor_alt_factor=2.0, seed=0
    ):
        rng = np.random.default_rng(seed)
        baseline = rng.lognormal(mean=4.0, sigma=1.0, size=n_genes)
        cols = {}
        normal_ids = []
        for i in range(n_normal):
            sid = f"normal{i}"
            normal_ids.append(sid)
            cols[sid] = baseline * rng.lognormal(0, 0.05, n_genes)
        for i in range(n_tumor):
            v = baseline.copy()
            v[altered_mask] *= tumor_alt_factor
            cols[f"tumor{i}"] = v * rng.lognormal(0, 0.05, n_genes)
        return (pd.DataFrame(cols, index=[f"g{i}" for i in range(n_genes)]), normal_ids)

    def test_single_normal_is_self_divide(self):
        """One normal => that normal's log2 is ~0 at every gene by construction.

        With one normal sample, the per-gene "median of normals" reduces to
        that single column, so the anchor divide is a self-divide. The normal
        sample's resulting log2 row collapses to ``safe_log2(1.0)`` = a small
        positive offset from the NULL_LOG2_COVERAGE shift, with vanishing
        per-gene variance. This is the design behavior to pin until the
        cnvkit-9wj redesign lands.
        """
        n_genes = 200
        altered = np.zeros(n_genes, dtype=bool)
        altered[:50] = True
        depths, normal_ids = self._make_cohort(10, 1, n_genes, altered, seed=1)
        log2 = rna.normalize_read_depths(depths.copy(), normal_ids)

        normal_log2 = log2[normal_ids[0]]
        # All values are within numerical noise of the safe_log2 floor offset.
        expected_offset = np.log2(1.0 + 2**rna.NULL_LOG2_COVERAGE)
        self.assertLess(abs(normal_log2.median() - expected_offset), 1e-3)
        self.assertLess(normal_log2.std(), 1e-3)
        # And the chrX-altered region in tumors is recovered at ~+1.0 (2x gain)
        tumor_cols = [c for c in log2.columns if c.startswith("tumor")]
        tumor_alt_med = log2[tumor_cols].loc[altered].median().median()
        self.assertGreater(tumor_alt_med, 0.8)
        self.assertLess(tumor_alt_med, 1.2)

    def test_multi_normal_carries_residual_qc_signal(self):
        """With >= 3 normals, individual normals' .cnr files retain non-zero spread.

        This documents the design rationale for the >=3-normals recommendation:
        the median-of-normals anchor leaves each normal with residual deviation
        from peer normals (interpretable as QC signal), rather than collapsing
        to a tautological zero.
        """
        n_genes = 200
        altered = np.zeros(n_genes, dtype=bool)
        altered[:50] = True
        depths, normal_ids = self._make_cohort(10, 3, n_genes, altered, seed=2)
        log2 = rna.normalize_read_depths(depths.copy(), normal_ids)

        # Each normal retains a non-degenerate spread across genes; not flat-zero.
        for nid in normal_ids:
            self.assertGreater(log2[nid].std(), 1e-3)

    def test_import_rna_warns_for_few_normals(self):
        """do_import_rna emits a warning for 0 < n_normals < 3."""
        with tempfile.TemporaryDirectory() as tmp:
            gene_res = os.path.join(tmp, "genes.tsv")
            rng = np.random.default_rng(7)
            n_genes = 80
            gene_ids = [f"ENSG{i:05d}.1" for i in range(n_genes)]
            with open(gene_res, "w") as fh:
                fh.write("# fixture\n")
                fh.write(
                    "Gene stable ID\tGC\tChr\tStart\tEnd\tName\tNCBI\tTxLen\tTSL\n"
                )
                for i, g in enumerate(gene_ids):
                    fh.write(
                        f"{g}\t{int(rng.integers(30, 65))}\t1\t"
                        f"{100000 + i * 5000}\t{102000 + i * 5000}\t"
                        f"G{i}\t{1000 + i}\t1500\ttsl1\n"
                    )
            baseline = rng.lognormal(mean=5.0, sigma=1.0, size=n_genes)
            fnames = []
            normal_fnames = []
            # 1 normal + 5 tumors -> should warn
            for tag, is_normal in [("normal-0", True)] + [
                (f"tumor-{i}", False) for i in range(5)
            ]:
                vals = baseline * rng.lognormal(0, 0.1, n_genes)
                f = os.path.join(tmp, f"{tag}.counts.txt")
                with open(f, "w") as fh:
                    for g, v in zip(gene_ids, vals, strict=True):
                        fh.write(f"{g}\t{int(v)}\n")
                fnames.append(f)
                if is_normal:
                    normal_fnames.append(f)

            with self.assertLogs(level="WARNING") as cm:
                _, cnrs = import_rna.do_import_rna(
                    fnames,
                    "counts",
                    gene_res,
                    normal_fnames=normal_fnames,
                )
                list(cnrs)  # exhaust generator so all logging fires
            joined = "\n".join(cm.output)
            self.assertIn("import-rna --normal", joined)
            self.assertIn("3 normals", joined)


if __name__ == "__main__":
    unittest.main()
