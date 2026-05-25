#!/usr/bin/env python
"""Tests for segment post-processing filters (cnvlib.segfilters)."""

import logging
import unittest

import pytest

logging.basicConfig(level=logging.ERROR, format="%(message)s")

import numpy as np
import pandas as pd
from skgenome import tabio, GenomicArray as GA

import cnvlib
from cnvlib import (
    cnary,
    fix,
    params,
    plots,
    reports,
    segfilters,
    segmentation,
    segmetrics,
    vary,
)


class SegfiltersTests(unittest.TestCase):
    """Tests for squash/bic segment filters."""

    def test_squash_region_nan_weights(self):
        """squash_region handles NaN weights without NaN in output."""
        df = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1", "chr1"],
                "start": [0, 1000, 2000],
                "end": [1000, 2000, 3000],
                "gene": ["G1", "G2", "G3"],
                "log2": [0.1, 0.2, 0.3],
                "probes": [5, 5, 5],
                "weight": [1.0, np.nan, 1.0],
            }
        )
        result = segfilters.squash_region(df)
        self.assertFalse(np.isnan(result["weight"].iat[0]))
        self.assertAlmostEqual(result["weight"].iat[0], 2.0)
        # log2 should be weighted average of bins with valid weights
        self.assertAlmostEqual(
            result["log2"].iat[0], np.average([0.1, 0.3], weights=[1.0, 1.0])
        )

        # All-NaN weights: should fall through to unweighted mean
        df_allnan = df.copy()
        df_allnan["weight"] = np.nan
        result2 = segfilters.squash_region(df_allnan)
        self.assertAlmostEqual(result2["log2"].iat[0], np.mean([0.1, 0.2, 0.3]))

    def test_bic_nan_weights(self):
        """BIC filter handles NaN segment weights without NaN RSS."""
        df = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "start": [0, 1000, 2000, 3000, 4000],
                "end": [1000, 2000, 3000, 4000, 5000],
                "gene": ["A", "B", "C", "D", "E"],
                "log2": [0.0, 0.1, 0.0, -0.1, 0.0],
                "probes": [10, 10, 10, 10, 10],
                "weight": [1.0, np.nan, 2.0, np.nan, 1.0],
            }
        )
        segarr = cnary.CopyNumArray(df)
        # Should not crash; NaN-weight segments get fallback RSS
        result = segfilters.bic(segarr)
        self.assertGreater(len(result), 0)
        self.assertFalse(np.isnan(result["weight"]).any())

    def test_squash_region_nan_gene(self):
        """squash_region handles NaN gene names without crashing (issue #900)."""
        df = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1"],
                "start": [0, 1000],
                "end": [1000, 2000],
                "gene": ["GeneA", float("nan")],
                "log2": [0.1, 0.2],
                "probes": [5, 5],
                "weight": [1.0, 1.0],
            }
        )
        result = segfilters.squash_region(df)
        gene_val = result["gene"].iat[0]
        self.assertIsInstance(gene_val, str)
        self.assertNotIn("nan", gene_val.lower())
        self.assertEqual(gene_val, "GeneA")

        # All-NaN genes should produce the placeholder "-"
        df_allnan = df.copy()
        df_allnan["gene"] = [float("nan"), float("nan")]
        result2 = segfilters.squash_region(df_allnan)
        self.assertEqual(result2["gene"].iat[0], "-")

    def test_squash_by_groups_no_cross_region_merge(self):
        """squash_by_groups merges only contiguous runs (#677).

        The old code combined a run-id cumsum with an integer chromosome
        index by *addition*, which collided on position-unsorted input and
        merged non-consecutive bins -- yielding segments with start > end and
        a bogus end coordinate (the constant "10488" in the ticket).  This is
        the shared segment-construction path used by HMM, CBS, and all segment
        filters, so it is not specific to HMM or to the segmetrics command that
        first surfaced it.  Boundaries are now detected by adjacency.
        """
        # chr1 appears in two runs (300-400, then 100-200) with a chr2 bin
        # between them; all three share the same state/level.
        df = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr2", "chr1"],
                "start": [300, 50, 100],
                "end": [400, 90, 200],
                "gene": ["-", "-", "-"],
                "log2": [0.0, 0.0, 0.0],
                "probes": [1, 1, 1],
                "weight": [1.0, 1.0, 1.0],
            }
        )
        cn = cnary.CopyNumArray(df)
        levels = pd.Series([2, 2, 2], index=cn.data.index)
        out = segfilters.squash_by_groups(cn, levels)
        # No segment may have start >= end ...
        self.assertTrue((out.start < out.end).all(), out.data)
        # ... and the two chr1 runs stay separate (3 segments, not 2) with all
        # probes conserved and no cross-chromosome absorption.
        self.assertEqual(len(out), 3)
        self.assertEqual(out["probes"].sum(), 3)

    def test_squash_region_unsorted_span(self):
        """squash_region spans [min start, max end] (#677).

        A run whose bins are not in ascending order must not yield start > end.
        """
        df = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1"],
                "start": [300, 100],
                "end": [400, 200],
                "gene": ["A", "A"],
                "log2": [0.0, 0.0],
                "probes": [1, 1],
                "weight": [1.0, 1.0],
            }
        )
        result = segfilters.squash_region(df)
        self.assertEqual(result["start"].iat[0], 100)
        self.assertEqual(result["end"].iat[0], 400)

    def test_squash_by_groups_cn1_cn2_nan(self):
        """Allele-specific squash splits on known cn1/cn2 changes but treats
        NaN (missing BAF) as 'no change', matching legacy grouping (#677).

        Without this, the adjacency check would split every NaN row off on its
        own (NaN != NaN), changing .cns output for `call`/`ampdel` on segments
        that lack BAF data (call sets cn1/cn2 = NaN there).
        """
        base = {
            "chromosome": ["chr1"] * 4,
            "start": [0, 100, 200, 300],
            "end": [100, 200, 300, 400],
            "gene": ["-"] * 4,
            "log2": [0.0] * 4,
            "probes": [1, 1, 1, 1],
            "weight": [1.0] * 4,
            "cn": [2, 2, 2, 2],
        }
        # A known cn1/cn2 change in the middle splits into 2 segments.
        known = cnary.CopyNumArray(
            pd.DataFrame({**base, "cn1": [1, 1, 2, 2], "cn2": [1, 1, 0, 0]})
        )
        out = segfilters.squash_by_groups(known, known["cn"])
        self.assertEqual(len(out), 2)
        # All-NaN cn1/cn2 (no BAF) stays merged into 1 segment, not 4.
        allnan = cnary.CopyNumArray(
            pd.DataFrame({**base, "cn1": [np.nan] * 4, "cn2": [np.nan] * 4})
        )
        out_nan = segfilters.squash_by_groups(allnan, allnan["cn"])
        self.assertEqual(len(out_nan), 1)
        self.assertTrue((out_nan.start < out_nan.end).all())

    def test_squash_by_groups_empty_and_single(self):
        """squash_by_groups handles empty and single-row input (#677)."""
        cols = ["chromosome", "start", "end", "gene", "log2", "probes", "weight"]
        empty = cnary.CopyNumArray(pd.DataFrame({c: [] for c in cols}))
        self.assertEqual(len(segfilters.squash_by_groups(empty, empty["log2"])), 0)
        one = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"],
                    "start": [10],
                    "end": [20],
                    "gene": ["G"],
                    "log2": [0.1],
                    "probes": [1],
                    "weight": [1.0],
                }
            )
        )
        out = segfilters.squash_by_groups(one, one["log2"])
        self.assertEqual(len(out), 1)
        self.assertEqual((out["start"].iat[0], out["end"].iat[0]), (10, 20))
