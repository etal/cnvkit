#!/usr/bin/env python
"""Tests for the reference command (and clustering used to build references)."""

import logging
import os
import shutil
import tempfile
import unittest
import warnings

import pytest

logging.basicConfig(level=logging.ERROR, format="%(message)s")
warnings.filterwarnings("ignore", category=DeprecationWarning)

import numpy as np
import pandas as pd
import pysam
from skgenome import tabio, GenomicArray as GA

import cnvlib
from conftest import linecount
from cnvlib import (
    access,
    antitarget,
    autobin,
    batch,
    bintest,
    call,
    cluster,
    cmdutil,
    cnary,
    commands,
    core,
    coverage,
    diagram,
    export,
    fix,
    heatmap,
    import_rna,
    importers,
    metrics,
    parallel,
    params,
    plots,
    purity,
    reference,
    reports,
    samutil,
    scatter,
    segfilters,
    segmentation,
    segmetrics,
    smoothing,
    vary,
)


class ReferenceTests(unittest.TestCase):
    """Tests for the reference and fix commands."""

    def test_reference(self):
        """The 'reference' command."""
        # Empty/unspecified antitargets
        nlines = linecount("formats/amplicon.cnr") - 1
        ref = commands.do_reference(["formats/amplicon.cnr"], ["formats/empty"])
        self.assertEqual(len(ref), nlines)
        ref = commands.do_reference(["formats/amplicon.cnr"])
        self.assertEqual(len(ref), nlines)
        # Empty/unspecified antitargets, flat reference
        nlines = linecount("formats/amplicon.bed")
        ref = commands.do_reference_flat("formats/amplicon.bed", "formats/empty")
        self.assertEqual(len(ref), nlines)
        ref = commands.do_reference_flat("formats/amplicon.bed")
        self.assertEqual(len(ref), nlines)
        # Misc
        ref = cnvlib.read("formats/reference-tr.cnn")
        targets, antitargets = reference.reference2regions(ref)
        self.assertLess(0, len(antitargets))
        self.assertEqual(len(antitargets), (ref["gene"] == "Background").sum())
        self.assertEqual(len(targets), len(ref) - len(antitargets))

    def test_reference_gender_input(self):
        """Test whether correct log2-ratios are calculated for sex chromosomes in reference command"""
        expected_haploid_log2 = [-1.0]
        expected_diploid_log2 = [0.0]

        # Test case when sample gender is provided by user
        ref_male = commands.do_reference(
            ["formats/ref_test_male.cnn"], female_samples=False
        )  # is_haploid_x_reference=False
        ref_female = commands.do_reference(
            ["formats/ref_test_female.cnn"], female_samples=True
        )  # is_haploid_x_reference=False
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_x_filter()].to_list(), expected_diploid_log2
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_y_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_female["log2"][ref_female.chr_x_filter()].to_list(),
            expected_diploid_log2,
        )
        ref_male = commands.do_reference(
            ["formats/ref_test_male.cnn"],
            female_samples=False,
            is_haploid_x_reference=True,
        )
        ref_female = commands.do_reference(
            ["formats/ref_test_female.cnn"],
            female_samples=True,
            is_haploid_x_reference=True,
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_x_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_y_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_female["log2"][ref_female.chr_x_filter()].to_list(),
            expected_haploid_log2,
        )

        # Test case when sample gender is guessed from input
        ref_male = commands.do_reference(
            ["formats/ref_test_male.cnn"]
        )  # is_haploid_x_reference=False
        ref_female = commands.do_reference(
            ["formats/ref_test_female.cnn"]
        )  # is_haploid_x_reference=False
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_x_filter()].to_list(), expected_diploid_log2
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_y_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_female["log2"][ref_female.chr_x_filter()].to_list(),
            expected_diploid_log2,
        )
        ref_male = commands.do_reference(
            ["formats/ref_test_male.cnn"], is_haploid_x_reference=True
        )
        ref_female = commands.do_reference(
            ["formats/ref_test_female.cnn"], is_haploid_x_reference=True
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_x_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_y_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_female["log2"][ref_female.chr_x_filter()].to_list(),
            expected_haploid_log2,
        )

    def test_fix(self):
        """The 'fix' command."""
        # Extract fake target/antitarget bins from a combined file
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg]
        rng = np.random.default_rng(42)  # Use fixed seed for reproducible tests
        tgt_bins.log2 += rng.standard_normal(len(tgt_bins)) / 5
        anti_bins = ref[is_bg]
        anti_bins.log2 += rng.standard_normal(len(anti_bins)) / 5
        blank_bins = cnary.CopyNumArray([])
        # Typical usage (hybrid capture)
        cnr = commands.do_fix(tgt_bins, anti_bins, ref)
        self.assertTrue(0 < len(cnr) <= len(ref))
        # Blank antitargets (WGS or amplicon)
        cnr = commands.do_fix(tgt_bins, blank_bins, ref[~is_bg])
        self.assertTrue(0 < len(cnr) <= len(tgt_bins))

    def test_fix_degenerate_no_nan(self):
        """do_fix never emits NaN log2/weight from degenerate input (#521/#524).

        Zero-coverage sentinel and NaN bins (from malformed inputs or dead
        regions in WGS) must not survive into the .cnr, where they would crash
        smoothing/CBS downstream.
        """
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg].copy()
        anti_bins = ref[is_bg].copy()
        # Inject degenerate values into a few sample bins
        tgt_log2 = tgt_bins["log2"].copy()
        tgt_log2.iloc[0] = np.nan
        tgt_log2.iloc[1] = params.NULL_LOG2_COVERAGE
        tgt_bins["log2"] = tgt_log2
        # Inject a degenerate value into the reference, too
        ref_log2 = ref["log2"].copy()
        ref_log2.iloc[2] = np.nan
        ref["log2"] = ref_log2
        cnr = commands.do_fix(tgt_bins, anti_bins, ref)
        self.assertFalse(np.isnan(cnr["log2"]).any())
        self.assertFalse(np.isnan(cnr["weight"]).any())


class ClusterTests(unittest.TestCase):
    """Tests for clustering used in reference building (cnvlib.cluster)."""

    def test_cluster_hierarchical(self):
        """Test hierarchical clustering on synthetic data with 2 clear groups."""
        from cnvlib.cluster import hierarchical

        rng = np.random.default_rng(42)
        n_bins = 100
        # Group A: 5 samples with similar pattern
        base_a = rng.standard_normal(n_bins)
        group_a = np.array([base_a + rng.normal(0, 0.1, n_bins) for _ in range(5)])
        # Group B: 5 samples with a different pattern
        base_b = rng.standard_normal(n_bins)
        group_b = np.array([base_b + rng.normal(0, 0.1, n_bins) for _ in range(5)])
        samples = np.vstack([group_a, group_b])

        clusters = hierarchical(samples, min_cluster_size=3)
        # Should produce exactly 2 clusters; all samples must be assigned
        self.assertEqual(len(clusters), 2)
        all_indices = sorted(idx for c in clusters for idx in c)
        self.assertEqual(all_indices, list(range(10)))
        # Check they separate the groups correctly
        c0 = set(clusters[0])
        c1 = set(clusters[1])
        group_a_set = set(range(5))
        group_b_set = set(range(5, 10))
        self.assertTrue(
            (c0 == group_a_set and c1 == group_b_set)
            or (c0 == group_b_set and c1 == group_a_set)
        )

    def test_cluster_kmedoids(self):
        """Test k-medoids clustering produces valid clusters."""
        from cnvlib.cluster import kmedoids

        rng = np.random.default_rng(42)
        n_bins = 100
        samples = np.vstack(
            [
                rng.standard_normal((5, n_bins)),
                rng.standard_normal((5, n_bins)) + 3,
            ]
        )
        clusters = kmedoids(samples, k=2)
        self.assertEqual(len(clusters), 2)
        all_indices = sorted(idx for c in clusters for idx in c)
        self.assertEqual(all_indices, list(range(10)))

    def test_create_clusters_column_names(self):
        """Test that create_clusters produces expected column names."""
        rng = np.random.default_rng(42)
        n_bins = 50
        n_samples = 6
        # Create logr_matrix with pseudocount row + sample rows
        logr = np.vstack(
            [
                np.zeros(n_bins),  # pseudocount
                rng.standard_normal((n_samples, n_bins)),
            ]
        )
        depths = np.vstack(
            [
                np.zeros(n_bins),  # pseudocount
                rng.uniform(10, 100, (n_samples, n_bins)),
            ]
        )
        sample_ids = [f"sample_{i}" for i in range(n_samples)]
        cols = reference.create_clusters(
            logr, depths, min_cluster_size=2, sample_ids=sample_ids
        )
        # Should have at least one cluster with log2_, depth_, spread_ columns
        self.assertTrue(len(cols) > 0)
        log2_keys = [k for k in cols if k.startswith("log2_")]
        depth_keys = [k for k in cols if k.startswith("depth_")]
        spread_keys = [k for k in cols if k.startswith("spread_")]
        self.assertTrue(len(log2_keys) > 0)
        self.assertEqual(len(log2_keys), len(depth_keys))
        self.assertEqual(len(log2_keys), len(spread_keys))
        # All column arrays should have n_bins elements
        for val in cols.values():
            self.assertEqual(len(val), n_bins)
