#!/usr/bin/env python
"""Tests for preprocessing commands: access, antitarget, autobin, target, coverage."""

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


class PreprocessingTests(unittest.TestCase):
    """Tests for preprocessing commands: access, antitarget, autobin, target, coverage."""

    def test_access(self):
        fasta = "formats/chrM-Y-trunc.hg19.fa"
        for min_gap_size, expect_nrows in ((None, 7), (500, 3), (1000, 2)):
            with self.subTest(min_gap_size=min_gap_size):
                acc = commands.do_access(
                    fasta, [], min_gap_size, skip_noncanonical=False
                )
                self.assertEqual(len(acc), expect_nrows)
        excludes = ["formats/dac-my.bed", "formats/my-targets.bed"]
        for min_gap_size, expect_nrows in (
            (None, 12),
            (2, 10),
            (20, 5),
            (200, 3),
            (2000, 2),
        ):
            with self.subTest(min_gap_size=min_gap_size, excludes=True):
                acc = commands.do_access(
                    fasta, excludes, min_gap_size, skip_noncanonical=False
                )
                self.assertEqual(len(acc), expect_nrows)
        # Dropping chrM, keeping only chrY
        acc = commands.do_access(fasta, excludes, 10, skip_noncanonical=True)
        self.assertEqual(len(acc), 5)

    @pytest.mark.slow
    def test_antitarget(self):
        """The 'antitarget' command."""
        baits = tabio.read_auto("formats/nv2_baits.interval_list")
        acs = tabio.read_auto("../data/access-5k-mappable.hg19.bed")
        self.assertLess(0, len(commands.do_antitarget(baits)))
        self.assertLess(0, len(commands.do_antitarget(baits, acs)))
        self.assertLess(0, len(commands.do_antitarget(baits, acs, 200000)))
        self.assertLess(0, len(commands.do_antitarget(baits, acs, 10000, 5000)))

    def test_autobin(self):
        """The 'autobin' command."""
        bam_fname = "formats/na12878-chrM-Y-trunc.bam"
        target_bed = "formats/my-targets.bed"
        targets = tabio.read(target_bed, "bed")
        access_bed = "../data/access-5k-mappable.hg19.bed"
        accessible = tabio.read(access_bed, "bed").filter(chromosome="chrY")
        for method in ("amplicon", "wgs", "hybrid"):
            (cov, bs), _ = autobin.do_autobin(
                bam_fname, method, targets=targets, access=accessible
            )
            self.assertGreater(cov, 0)
            self.assertGreater(bs, 0)

    def test_autobin_chrom_name_mismatch(self):
        """WGS autobin gives a clear, actionable error -- not the cryptic
        'cannot convert float NaN to integer' -- when the BAM and access
        regions share no chromosome names (#421)."""
        bam_fname = "formats/na12878-chrM-Y-trunc.bam"
        # Access regions on a contig absent from the BAM => no shared chroms,
        # so the estimated read depth is NaN.
        access_arr = GA.from_columns(
            {"chromosome": ["nonexistent_contig"], "start": [0], "end": [100000]}
        )
        with self.assertRaises(ValueError) as cm:
            autobin.do_autobin(bam_fname, "wgs", access=access_arr)
        msg = str(cm.exception)
        self.assertNotIn("NaN", msg)
        self.assertIn("chromosome", msg.lower())

    def test_bedcov_filters_absent_contigs(self):
        """bedcov keeps regions on BAM contigs and drops those on absent
        contigs instead of erroring out entirely (#620)."""
        bam = "formats/na12878-chrM-Y-trunc.bam"
        samutil.ensure_bam_index(bam)  # self-contained: don't rely on test order
        bam_chroms = samutil.get_bam_chroms(bam)
        with tempfile.NamedTemporaryFile("w+t", suffix=".bed", delete=False) as f:
            f.write("chrM\t251\t277\tfoo\n")
            f.write("absent_contig\t100\t200\tbar\n")
            bed = f.name
        try:
            table = coverage.bedcov(bed, bam, 0, bam_chroms=bam_chroms)
            self.assertEqual(list(table.chromosome.unique()), ["chrM"])
            self.assertEqual(len(table), 1)
        finally:
            os.unlink(bed)

    def test_bedcov_all_absent_contigs(self):
        """bedcov on an all-absent BED: raise a clear error by default, but
        return empty when allow_empty=True (the per-chunk parallel path)."""
        bam = "formats/na12878-chrM-Y-trunc.bam"
        samutil.ensure_bam_index(bam)  # self-contained: don't rely on test order
        bam_chroms = samutil.get_bam_chroms(bam)
        with tempfile.NamedTemporaryFile("w+t", suffix=".bed", delete=False) as f:
            f.write("absent_contig\t100\t200\tbar\n")
            bed = f.name
        try:
            with self.assertRaises(ValueError) as cm:
                coverage.bedcov(bed, bam, 0, bam_chroms=bam_chroms)
            self.assertIn("don't match", str(cm.exception))
            empty = coverage.bedcov(
                bed, bam, 0, bam_chroms=bam_chroms, allow_empty=True
            )
            self.assertEqual(len(empty), 0)
        finally:
            os.unlink(bed)

    def test_coverage_skips_marked_duplicates(self):
        """Coverage ignores reads flagged as duplicates in both the depth
        (bedcov) and by-count paths, so marking duplicates (e.g. Picard
        MarkDuplicates) is equivalent to physically removing them (#689)."""
        seqlen = 50
        n_total, n_dup = 10, 6  # 4 reads survive the duplicate filter
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
        }
        tmpdir = tempfile.mkdtemp()
        bam = os.path.join(tmpdir, "dup.bam")
        with pysam.AlignmentFile(bam, "wb", header=header) as out:
            for i in range(n_total):
                a = pysam.AlignedSegment()
                a.query_name = f"r{i}"
                a.query_sequence = "A" * seqlen
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 100
                a.mapping_quality = 60
                a.cigarstring = f"{seqlen}M"
                a.query_qualities = pysam.qualitystring_to_array("I" * seqlen)
                a.is_duplicate = i >= (n_total - n_dup)
                out.write(a)
        pysam.index(bam)
        bed = os.path.join(tmpdir, "region.bed")
        with open(bed, "w") as fh:
            fh.write("chr1\t90\t160\tregionA\n")
        n_kept = n_total - n_dup
        try:
            # Depth path (default): basecount counts only non-duplicate reads
            table = coverage.bedcov(bed, bam, 0)
            self.assertEqual(float(table["basecount"].iloc[0]), n_kept * seqlen)
            # By-count path: same surviving read count
            bamfile = pysam.AlignmentFile(bam, "rb")
            count, _row = coverage.region_depth_count(
                bamfile, "chr1", 90, 160, "regionA", 0
            )
            self.assertEqual(count, n_kept)
        finally:
            shutil.rmtree(tmpdir)

    def test_coverage_partial_chrom_mismatch(self):
        """coverage tolerates a BED mixing present and absent contigs,
        returning only the present-contig regions, for any process count and
        either depth/count method (#620: failure was process-count- and
        method-dependent)."""
        bam = "formats/na12878-chrM-Y-trunc.bam"
        with tempfile.NamedTemporaryFile("w+t", suffix=".bed", delete=False) as f:
            f.write("chrM\t251\t277\tfoo\n")
            f.write("chrY\t11170\t11213\tbar\n")
            f.write("absent_contig\t100\t200\tbaz\n")
            bed = f.name
        try:
            for by_count in (False, True):
                for nprocs in (1, 2):
                    with self.subTest(by_count=by_count, processes=nprocs):
                        cna = commands.do_coverage(
                            bed, bam, by_count=by_count, processes=nprocs
                        )
                        self.assertEqual(
                            sorted(cna.chromosome.unique()), ["chrM", "chrY"]
                        )
                        self.assertEqual(len(cna), 2)
        finally:
            os.unlink(bed)

    def test_coverage_all_chrom_mismatch(self):
        """coverage raises a clear error when no BED chromosome matches the
        BAM, for any process count or method (#620)."""
        bam = "formats/na12878-chrM-Y-trunc.bam"
        with tempfile.NamedTemporaryFile("w+t", suffix=".bed", delete=False) as f:
            f.write("absent_contig\t100\t200\tbaz\n")
            bed = f.name
        try:
            for by_count in (False, True):
                for nprocs in (1, 2):
                    with self.subTest(by_count=by_count, processes=nprocs):
                        with self.assertRaises(ValueError) as cm:
                            commands.do_coverage(
                                bed, bam, by_count=by_count, processes=nprocs
                            )
                        self.assertIn("don't match", str(cm.exception))
        finally:
            os.unlink(bed)

    @pytest.mark.slow
    def test_coverage(self):
        """The 'coverage' command."""
        bed = "formats/my-targets.bed"
        bam = "formats/na12878-chrM-Y-trunc.bam"
        # Cover distinct code paths: depth-based default, count-based with
        # mapq filtering, and multiprocessing.
        for by_count, min_mapq, nprocs in [
            (False, 0, 1),
            (True, 30, 1),
            (False, 0, 2),
        ]:
            cna = commands.do_coverage(
                bed, bam, by_count=by_count, min_mapq=min_mapq, processes=nprocs
            )
            self.assertEqual(len(cna), 4)
            self.assertTrue((cna.log2 != 0).any())
            self.assertGreater(cna.log2.nunique(), 1)

    def test_coverage_cram(self):
        """The 'coverage' command with CRAM input."""
        bed = "formats/my-targets.bed"
        cram = "formats/na12878-chrM-Y-trunc.cram"
        fasta = "formats/chrM-Y-trunc.hg19.fa"
        bam = "formats/na12878-chrM-Y-trunc.bam"
        # Verify CRAM coverage matches BAM coverage
        cna_cram = commands.do_coverage(bed, cram, fasta=fasta)
        cna_bam = commands.do_coverage(bed, bam)
        self.assertEqual(len(cna_cram), len(cna_bam))
        np.testing.assert_array_almost_equal(
            cna_cram["depth"].values, cna_bam["depth"].values, decimal=1
        )

    def test_coverage_bedgraph(self):
        """The 'coverage' command with bedGraph input."""
        bed = "formats/my-targets.bed"
        bedgraph = "formats/na12878-chrM-Y-trunc.bed.gz"

        # Test basic bedGraph coverage
        cna = commands.do_coverage(bed, bedgraph)
        self.assertEqual(len(cna), 4)
        self.assertTrue((cna.log2 != 0).any())
        self.assertGreater(cna.log2.nunique(), 1)

        # Verify sample_id is set correctly
        self.assertEqual(cna.sample_id, "na12878-chrM-Y-trunc")

        # Verify depth values are calculated
        self.assertTrue((cna["depth"] >= 0).all())

    def test_coverage_bedgraph_vs_bam(self):
        """Compare bedGraph and BAM coverage outputs."""
        bed = "formats/my-targets.bed"
        bam = "formats/na12878-chrM-Y-trunc.bam"
        bedgraph = "formats/na12878-chrM-Y-trunc.bed.gz"

        # Get coverage from both inputs
        cna_bam = commands.do_coverage(bed, bam)
        cna_bedgraph = commands.do_coverage(bed, bedgraph)

        # Should have same number of regions
        self.assertEqual(len(cna_bam), len(cna_bedgraph))

        # Convert to DataFrames for easier comparison
        df_bam = cna_bam.data.copy()
        df_bedgraph = cna_bedgraph.data.copy()

        # Sort both by chromosome, start, end
        df_bam = df_bam.sort_values(["chromosome", "start", "end"]).reset_index(
            drop=True
        )
        df_bedgraph = df_bedgraph.sort_values(
            ["chromosome", "start", "end"]
        ).reset_index(drop=True)

        # Chromosomes and positions should match exactly
        pd.testing.assert_series_equal(
            df_bam["chromosome"], df_bedgraph["chromosome"], check_names=False
        )
        pd.testing.assert_series_equal(
            df_bam["start"], df_bedgraph["start"], check_names=False
        )
        pd.testing.assert_series_equal(
            df_bam["end"], df_bedgraph["end"], check_names=False
        )

        # Depth values should be reasonably close
        # Note: Small differences expected due to different depth calculation methods
        # (pysam.bedcov vs bedtools genomecov)
        np.testing.assert_allclose(
            df_bam["depth"].values,
            df_bedgraph["depth"].values,
            rtol=0.05,  # Allow 5% relative difference
            atol=1.0,  # Allow 1.0 absolute difference for low-depth regions
            err_msg="bedGraph and BAM coverage depths should be reasonably close",
        )

    def test_coverage_bedgraph_missing_index(self):
        """Test error handling when tabix index is missing."""
        bed = "formats/my-targets.bed"

        # Create a temporary bedGraph without index
        with tempfile.NamedTemporaryFile(suffix=".bed.gz", delete=False) as tmp:
            tmp_path = tmp.name
            # Copy bedGraph content but not the index
            with open("formats/na12878-chrM-Y-trunc.bed.gz", "rb") as src:
                tmp.write(src.read())

        try:
            # Should raise FileNotFoundError for missing index
            with self.assertRaises(FileNotFoundError) as cm:
                commands.do_coverage(bed, tmp_path)

            self.assertIn("tabix index", str(cm.exception))
        finally:
            os.unlink(tmp_path)

    def test_coverage_bedgraph_chr_naming(self):
        """Test chromosome name mismatch handling (chr1 vs 1)."""
        bed = "formats/my-targets.bed"
        # This bedGraph has chromosomes named without 'chr' prefix (M, Y)
        # while the BED file uses 'chr' prefix (chrM, chrY)
        bedgraph_nochr = "formats/na12878-M-Y-trunc-nochr.bed.gz"

        # Should still work by auto-matching chromosome names
        cna = commands.do_coverage(bed, bedgraph_nochr)

        # Should get same number of regions
        self.assertEqual(len(cna), 4)

        # Should have valid coverage values
        self.assertTrue((cna.log2 != 0).any())
        self.assertGreater(cna.log2.nunique(), 1)

    @pytest.mark.slow
    def test_target_annotate_refflat_assigns_genes(self):
        """'target --annotate refFlat' actually populates gene names (#688).

        The refFlat reader parses the gene column, tabio sorts the annotation,
        and into_ranges overlaps it onto the baits. Guards against silent
        annotation failures: a real bait set should come back almost entirely
        gene-named, with names drawn from the refFlat's first column.
        """
        annot_fname = "formats/refflat-mini.txt"
        baits = tabio.read_auto("formats/baits-funky.bed")
        out = commands.do_target(
            baits, do_split=True, avg_size=200, annotate=annot_fname
        )
        named = [g for g in out["gene"] if g != "-"]
        # Nearly all bins on these gene-dense baits should be annotated
        self.assertGreater(len(named), 0.9 * len(out))
        # Names come from the refFlat first column (gene symbols)
        all_names = {n for g in named for n in g.split(",")}
        self.assertIn("DDX11L1", all_names)
        self.assertIn("WASH7P", all_names)

    def test_target_annotate_wrong_build_warns(self):
        """Warn when annotation shares chrom names but assigns no genes (#688).

        Mimics a wrong-genome-build refFlat: chromosome names match (so
        compare_chrom_names passes), but no coordinates overlap, leaving every
        region unnamed. Without the warning this failure is silent.
        """
        annot_fname = "formats/refflat-mini.txt"  # chr1 genes end well before 210 Mb
        baits = GA.from_columns(
            {
                "chromosome": ["chr1", "chr1"],
                "start": [210_000_000, 210_010_000],
                "end": [210_001_000, 210_011_000],
                "gene": ["-", "-"],
            }
        )
        with self.assertLogs(level="WARNING") as cm:
            out = commands.do_target(baits, annotate=annot_fname)
        self.assertTrue(all(g == "-" for g in out["gene"]))
        self.assertTrue(any("different genome build" in m for m in cm.output))

    def test_target(self):
        """The 'target' command."""
        annot_fname = "formats/refflat-mini.txt"
        for bait_fname in (
            "formats/nv2_baits.interval_list",
            "formats/amplicon.bed",
            "formats/baits-funky.bed",
        ):
            baits = tabio.read_auto(bait_fname)
            bait_len = len(baits)
            # No splitting: w/o and w/ re-annotation
            r1 = commands.do_target(baits)
            self.assertEqual(len(r1), bait_len)
            r1a = commands.do_target(baits, do_short_names=True, annotate=annot_fname)
            self.assertEqual(len(r1a), len(r1))
            # Splitting, w/o and w/ re-annotation
            r2 = commands.do_target(
                baits, do_short_names=True, do_split=True, avg_size=100
            )
            self.assertGreater(len(r2), len(r1))
            for _c, subarr in r2.by_chromosome():
                self.assertTrue(subarr.start.is_monotonic_increasing, bait_fname)
                self.assertTrue(subarr.end.is_monotonic_increasing, bait_fname)
                # Bins are non-overlapping; next start >= prev. end
                self.assertTrue(
                    (
                        (subarr.start.to_numpy()[1:] - subarr.end.to_numpy()[:-1]) >= 0
                    ).all()
                )
            r2a = commands.do_target(
                baits,
                do_short_names=True,
                do_split=True,
                avg_size=100,
                annotate=annot_fname,
            )
            self.assertEqual(len(r2a), len(r2))
            # Original regions object should be unmodified
            self.assertEqual(len(baits), bait_len)
