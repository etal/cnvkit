#!/usr/bin/env python
"""Tests for export commands."""

import logging
import math
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


class ExportTests(unittest.TestCase):
    """Tests for export commands."""

    def test_export_bed_vcf(self):
        """The 'export' command for formats with absolute copy number."""
        for fname, ploidy, is_f in [
            ("tr95t.cns", 2, True),
            ("cl_seq.cns", 6, True),
            ("amplicon.cns", 2, False),
        ]:
            with self.subTest(fname=fname, ploidy=ploidy):
                cns = cnvlib.read("formats/" + fname)
                # BED
                for show in ("ploidy", "variant", "all"):
                    tbl_bed = export.export_bed(
                        cns, ploidy, True, None, is_f, cns.sample_id, show
                    )
                    if show == "all":
                        self.assertEqual(len(tbl_bed), len(cns), f"{fname} {ploidy}")
                    else:
                        self.assertLess(len(tbl_bed), len(cns))
                # VCF
                _vheader, vcf_body = export.export_vcf(cns, ploidy, True, None, is_f)
                self.assertTrue(0 < len(vcf_body.splitlines()) < len(cns))

    def test_export_seg_values_match_input(self):
        """SEG output reproduces input segment values, not just row counts.

        Guards the clinical deliverable against value-scrambling regressions.
        """
        cns = cnvlib.read("formats/tr95t.cns")
        seg = export.export_seg(["formats/tr95t.cns"])
        self.assertEqual(len(seg), len(cns))
        self.assertTrue((seg["ID"] == cns.sample_id).all())
        self.assertTrue((seg["chrom"].to_numpy() == cns.chromosome.to_numpy()).all())
        np.testing.assert_allclose(
            seg["seg.mean"].to_numpy(), cns["log2"].to_numpy(), atol=1e-6
        )
        np.testing.assert_array_equal(
            seg["num.mark"].to_numpy(), cns["probes"].to_numpy()
        )
        # SEG uses 1-based start, inclusive end
        np.testing.assert_array_equal(
            seg["loc.start"].to_numpy(), cns.start.to_numpy() + 1
        )
        np.testing.assert_array_equal(seg["loc.end"].to_numpy(), cns.end.to_numpy())
        self.assertTrue((seg["loc.start"] < seg["loc.end"]).all())

    def test_export_bed_values(self):
        """BED export emits finite non-negative integer copy numbers."""
        cns = cnvlib.read("formats/tr95t.cns")
        bed = export.export_bed(cns, 2, True, None, True, cns.sample_id, "all")
        self.assertEqual(len(bed), len(cns))
        self.assertTrue(
            (bed["chromosome"].to_numpy() == cns.chromosome.to_numpy()).all()
        )
        self.assertTrue((bed["start"] < bed["end"]).all())
        self.assertEqual(bed["ncopies"].dtype.kind, "i")
        self.assertTrue((bed["ncopies"] >= 0).all())
        # A copy-number gain segment (max log2) must export as > ploidy copies
        gain_idx = int(np.argmax(cns["log2"].to_numpy()))
        self.assertGreater(bed["ncopies"].iloc[gain_idx], 2)

    def test_export_vcf_values(self):
        """VCF export emits well-formed, non-degenerate SV records.

        Note: <DUP>/<DEL> is decided per chromosome's expected copies (sex-aware),
        so it does not track the ploidy-relative FOLD_CHANGE_LOG sign on the sex
        chromosomes; we assert structure and consistency, not that cross-field tie.
        """
        cns = cnvlib.read("formats/tr95t.cns")
        _vheader, vcf_body = export.export_vcf(cns, 2, True, None, True)
        var_lines = [ln for ln in vcf_body.splitlines() if not ln.startswith("#")]
        self.assertTrue(0 < len(var_lines) < len(cns))  # only non-neutral segments
        svtypes, fold_changes = set(), []
        for ln in var_lines:
            fields = ln.split("\t")
            self.assertEqual(len(fields), 10)  # CHROM..FORMAT + 1 sample
            info = dict(kv.split("=", 1) for kv in fields[7].split(";") if "=" in kv)
            self.assertLess(int(fields[1]), int(info["END"]))  # POS < END
            self.assertTrue(math.isfinite(float(info["FOLD_CHANGE_LOG"])))
            self.assertEqual(fields[4], f"<{info['SVTYPE']}>")  # ALT matches SVTYPE
            svtypes.add(info["SVTYPE"])
            fold_changes.append(float(info["FOLD_CHANGE_LOG"]))
            # CN, when present in FORMAT, parses as a non-negative integer
            fmt = fields[8].split(":")
            if "CN" in fmt:
                cn = int(fields[9].split(":")[fmt.index("CN")])
                self.assertGreaterEqual(cn, 0)
        self.assertEqual(svtypes, {"DUP", "DEL"})  # tr95t has both gains and losses
        self.assertGreater(len(set(fold_changes)), 1)  # not a stuck/degenerate value

    def test_export_cdt_jtv(self):
        """The 'export' command for CDT and Java TreeView formats."""
        fnames = ["formats/p2-20_1.cnr", "formats/p2-20_2.cnr"]
        sample_ids = list(map(core.fbase, fnames))
        nrows = linecount(fnames[0]) - 1
        for fmt_key, header2 in (("cdt", 2), ("jtv", 0)):
            table = export.merge_samples(fnames)
            formatter = export.EXPORT_FORMATS[fmt_key]
            _oh, outrows = formatter(sample_ids, table)
            self.assertEqual(len(list(outrows)), nrows + header2)

    def test_export_nexus(self):
        """The 'export nexus-basic' and 'nexus-ogt' commands."""
        cnr = cnvlib.read("formats/amplicon.cnr")
        table_nb = export.export_nexus_basic(cnr)
        self.assertEqual(len(table_nb), len(cnr))
        varr = commands.load_het_snps(
            "formats/na12878_na12882_mix.vcf", None, None, 15, None
        )
        table_ogt = export.export_nexus_ogt(cnr, varr, 0.05)
        self.assertEqual(len(table_ogt), len(cnr))

    def test_export_seg(self):
        """The 'export seg' command."""
        seg_rows = export.export_seg(["formats/tr95t.cns"])
        self.assertGreater(len(seg_rows), 0)
        seg2_rows = export.export_seg(["formats/tr95t.cns", "formats/cl_seq.cns"])
        self.assertGreater(len(seg2_rows), len(seg_rows))

    def test_export_theta(self):
        """The 'export theta' command."""
        segarr = cnvlib.read("formats/tr95t.cns")
        len_seg_auto = len(segarr.autosomes())
        table_theta = export.export_theta(segarr, None)
        self.assertEqual(len(table_theta), len_seg_auto)
        ref = cnvlib.read("formats/reference-tr.cnn")
        table_theta = export.export_theta(segarr, ref)
        self.assertEqual(len(table_theta), len_seg_auto)
        varr = commands.load_het_snps(
            "formats/na12878_na12882_mix.vcf", "NA12882", "NA12878", 15, None
        )
        tumor_snps, normal_snps = export.export_theta_snps(varr)
        self.assertLess(len(tumor_snps), len(varr))
        self.assertGreater(len(tumor_snps), 0)
        self.assertLess(len(normal_snps), len(varr))
        self.assertGreater(len(normal_snps), 0)
