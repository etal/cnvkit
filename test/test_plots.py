#!/usr/bin/env python
"""Tests for plotting commands: scatter, heatmap, diagram."""

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
from conftest import linecount

import cnvlib
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
from skgenome import GenomicArray as GA
from skgenome import tabio


class PlotTests(unittest.TestCase):
    """Smoke tests for plotting commands."""

    def test_scatter(self):
        """The 'scatter' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        fig = scatter.do_scatter(cnarr, segarr)
        self.assertIsNotNone(fig)
        # With a gene zoom
        fig = scatter.do_scatter(cnarr, segarr, show_gene="BRAF")
        self.assertIsNotNone(fig)

    def test_scatter_genome_y_floor(self):
        """Genome-wide autoscale floors y_min so a single deep homozygous
        deletion can't compress the whole plot (gh#385)."""
        # lazy: defer matplotlib import to keep headless test collection fast
        from matplotlib import pyplot  # noqa: PLC0415

        probes = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 100, 200, "A", 0.0],
                ["chr1", 200, 300, "A", -0.1],
                ["chr1", 300, 400, "A", -12.0],
                ["chr2", 100, 200, "B", 0.05],
                ["chr2", 200, 300, "B", -0.2],
            ],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        segments = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 100, 400, "A", -12.0],
                ["chr2", 100, 300, "B", 0.0],
            ],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        segments.data["probes"] = [3, 2]
        _fig, ax = pyplot.subplots()
        scatter.cnv_on_genome(ax, probes, segments)
        y_lo, _y_hi = ax.get_ylim()
        pyplot.close(_fig)
        self.assertGreaterEqual(
            y_lo,
            scatter.AUTO_Y_MIN_FLOOR,
            "Deep deletion must not pull the genome-wide y-axis below the floor",
        )
        # An explicit --y-min still overrides the floor
        _fig, ax = pyplot.subplots()
        scatter.cnv_on_genome(ax, probes, segments, y_min=-15.0, y_max=2.0)
        y_lo, _y_hi = ax.get_ylim()
        pyplot.close(_fig)
        self.assertAlmostEqual(y_lo, -15.0)

    def test_heatmap(self):
        """The 'heatmap' command."""
        cnarrs = [cnvlib.read("formats/amplicon.cnr")]
        ax = heatmap.do_heatmap(cnarrs)
        self.assertIsNotNone(ax)
        # With desaturation
        ax = heatmap.do_heatmap(cnarrs, do_desaturate=True)
        self.assertIsNotNone(ax)


class GeneCoordsTests(unittest.TestCase):
    """Tests for plots.gene_coords_by_name gene-label selection."""

    def test_gene_coords_by_name(self):
        """`-g` labels only requested genes, not co-binned neighbors (gh#458).

        A bin whose `gene` column packs several names (e.g. "ERBB2,MIR4728")
        must not surface the unrequested neighbor (MIR4728) when only ERBB2
        is selected.
        """
        cnarr = cnary.CopyNumArray.from_rows(
            [
                ["chr17", 37800000, 37850000, "STARD3", 0.0],
                ["chr17", 37850000, 37860000, "ERBB2,MIR4728", 0.0],
                ["chr17", 37860000, 37870000, "ERBB2", 0.0],
                ["chr17", 37880000, 37890000, "GRB7", 0.0],
            ]
        )
        # Single requested gene: co-binned MIR4728 must be hidden
        coords = plots.gene_coords_by_name(cnarr, ["ERBB2"])
        self.assertEqual(list(coords.keys()), ["chr17"])
        self.assertEqual(len(coords["chr17"]), 1)
        start, end, name = coords["chr17"][0]
        self.assertEqual(name, "ERBB2")
        self.assertEqual((start, end), (37850000, 37870000))
        # When MIR4728 *is* requested, it should appear
        names_seen = set()
        for _s, _e, label in plots.gene_coords_by_name(cnarr, ["ERBB2", "MIR4728"])[
            "chr17"
        ]:
            names_seen.update(label.split(","))
        self.assertEqual(names_seen, {"ERBB2", "MIR4728"})
