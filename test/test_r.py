#!/usr/bin/env python
"""Unit tests for CNVkit that require an R installation."""

import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np
import pytest

logging.basicConfig(level=logging.ERROR, format="%(message)s")

import cnvlib
from cnvlib import core, segmentation
from cnvlib.segmentation.cbs import CBS_RSCRIPT


class RTests(unittest.TestCase):
    """Tests that depend on the R statistical environment."""

    def setUp(self):
        self.tas_cnr = cnvlib.read("formats/amplicon.cnr")

    @pytest.mark.slow
    def test_cbs(self):
        _test_method(self, "cbs")

    @pytest.mark.slow
    def test_cbs_weights_align_with_dropped_bins(self):
        """Regression test for issue #868.

        R's ``read.delim`` reads an empty ``start`` or a literal ``NA``
        chromosome as ``NA``, and DNAcopy's ``CNA()`` silently drops those
        markers. The weights vector passed to ``segment()`` must be filtered
        the same way, or R raises "length of weights should be the same as the
        number of probes".
        """
        rows = [
            "\t".join(("chromosome", "start", "end", "gene", "log2", "depth", "weight"))
        ]
        # A clean two-level signal so CBS has something to segment ...
        for i in range(20):
            start = 100 * (i + 1)
            log2 = 0.5 if i < 10 else -0.5
            rows.append(f"chr1\t{start}\t{start + 50}\tG\t{log2}\t30\t0.6")
        # ... plus two bins CNA() drops: an empty start and an "NA" chromosome.
        rows.insert(6, "chr1\t\t450\tG\t0.5\t30\t0.6")
        rows.insert(16, "NA\t9999\t10050\tG\t-0.5\t30\t0.6")
        with tempfile.NamedTemporaryFile(suffix=".cnr", mode="w+t") as tmp:
            tmp.write("\n".join(rows) + "\n")
            tmp.flush()
            with core.temp_write_text(CBS_RSCRIPT, mode="w+t") as script_fname:
                seg_out = core.call_quiet(
                    "Rscript",
                    "--no-restore",
                    "--no-environ",
                    script_fname,
                    tmp.name,
                    "test868",
                    "0.0001",
                    "FALSE",
                )
        # No exception means the weights stayed aligned with the kept probes.
        self.assertIn("test868", seg_out.decode())

    @pytest.mark.slow
    def test_cbs_sample_id_with_quote(self):
        """Regression: a sample ID containing a single quote must not break R.

        The ID used to be interpolated into a single-quoted R string literal,
        so an apostrophe (e.g. ``O'Brien_001``) produced invalid R source and a
        subprocess parse error. It is now passed as a command-line argument and
        needs no escaping. Exercised end-to-end through ``do_segmentation`` so
        the Python caller's argument plumbing is covered too.
        """
        cnr = self.tas_cnr.copy()
        cnr.meta["sample_id"] = "O'Brien_001"
        cns, raw_str = segmentation.do_segmentation(
            cnr, "cbs", processes=1, save_dataframe=True
        )
        self.assertGreater(len(cns), 0)
        self.assertIn("O'Brien_001", raw_str)

    @pytest.mark.slow
    def test_cbs_probes_path_with_awkward_chars(self):
        """Regression: a probes-file path with a quote, space, or backslash must
        not break R.

        The path used to be interpolated into a double-quoted R string literal;
        it is now passed as a command-line argument to ``read.delim``.
        """
        rows = [
            "\t".join(("chromosome", "start", "end", "gene", "log2", "depth", "weight"))
        ]
        for i in range(20):
            start = 100 * (i + 1)
            log2 = 0.5 if i < 10 else -0.5
            rows.append(f"chr1\t{start}\t{start + 50}\tG\t{log2}\t30\t0.6")
        tmpdir = tempfile.mkdtemp()
        try:
            probes_fname = os.path.join(tmpdir, "O'Brien \\sample.cnr")
            with open(probes_fname, "w") as fh:
                fh.write("\n".join(rows) + "\n")
            with core.temp_write_text(CBS_RSCRIPT, mode="w+t") as script_fname:
                seg_out = core.call_quiet(
                    "Rscript",
                    "--no-restore",
                    "--no-environ",
                    script_fname,
                    probes_fname,
                    "sample",
                    "0.0001",
                    "FALSE",
                )
        finally:
            shutil.rmtree(tmpdir)
        self.assertIn("sample", seg_out.decode())

    @pytest.mark.slow
    def test_cbs_nan_log2_in_memory(self):
        """CBS segments an in-memory .cnr with NaN log2 values (#881).

        The reporter ran ``segment`` (no --drop-low-coverage) on a .cnr with NaN
        values and hit ``subprocess.CalledProcessError``. The file-read path is
        guarded by read_tab, but an in-memory CopyNumArray (e.g. from batch)
        reaches segmentation with the NaN bins intact. They must be dropped
        before the Savitzky-Golay outlier filter and DNAcopy, not crash.
        """
        cnr = self.tas_cnr.copy()
        log2 = cnr["log2"].to_numpy().copy()
        # Scatter a few NaN values across the bins, as a degenerate .cnr would
        log2[[3, 50, 120, len(log2) - 2]] = np.nan
        cnr["log2"] = log2
        cns = segmentation.do_segmentation(cnr, "cbs", processes=1)
        self.assertGreater(len(cns), 0)
        self.assertFalse(np.isnan(cns["log2"].to_numpy()).any())

    def test_cbs_invocation_omits_vanilla(self):
        """Regression test for issue #491.

        ``--vanilla`` bundles ``--no-init-file``/``--no-site-file``, so R skips
        the user's ``.Rprofile``/``Rprofile.site`` -- breaking custom installs
        that set ``.libPaths()`` there. The segmentation R call must pass only
        ``--no-restore --no-environ``. We mock ``call_quiet`` to capture the
        argv before R runs, so this test needs no R installation (hence it is
        not marked ``slow``).
        """
        captured = []

        class _StopBeforeR(Exception):
            pass

        def fake_call_quiet(*args):
            captured.append(args)
            raise _StopBeforeR

        with (
            mock.patch.object(core, "call_quiet", fake_call_quiet),
            self.assertRaises(_StopBeforeR),
        ):
            segmentation._do_segmentation(
                self.tas_cnr, "cbs", None, 0.0001, skip_outliers=0
            )

        argv = captured[0]
        self.assertNotIn("--vanilla", argv)
        self.assertNotIn("--no-init-file", argv)
        self.assertNotIn("--no-site-file", argv)
        self.assertIn("--no-restore", argv)
        self.assertIn("--no-environ", argv)

    @pytest.mark.slow
    def test_cbs_all_bins_unsegmentable(self):
        """Issue #868: an arm whose bins are all unsegmentable yields no segments.

        When every bin has a chromosome R reads as ``NA``, ``CNA()`` drops them
        all and DNAcopy emits a header-only SEG table. ``_do_segmentation`` must
        return an empty array, not crash the SEG parser. Tested directly (not via
        ``do_segmentation``), whose process pool would otherwise swallow the error.
        """
        cnr = self.tas_cnr[self.tas_cnr["chromosome"] == "chr2"]
        df = cnr.data.head(6).reset_index(drop=True).copy()
        df["chromosome"] = "NA"
        arm = cnr.as_dataframe(df)
        segarr = segmentation._do_segmentation(
            arm, "cbs", None, 0.0001, skip_outliers=0
        )
        self.assertEqual(len(segarr), 0)

    @pytest.mark.slow
    def test_cbs_empty_arm_merges(self):
        """Issue #868: an empty (unsegmentable) arm merges cleanly via the pool.

        ``do_segmentation`` concatenates per-arm results and stitches the
        ``save_dataframe`` SEG strings in the parent process (not inside the
        exception-swallowing pool), so an arm that yields no segments must not
        crash either step.
        """
        cnr = self.tas_cnr[self.tas_cnr["chromosome"] == "chr2"]
        df = cnr.data.reset_index(drop=True).copy()
        df.loc[5, "chromosome"] = "NA"  # one normal chr2 arm + one empty NA arm
        mixed = cnr.as_dataframe(df)
        cns, rstr = segmentation.do_segmentation(
            mixed, "cbs", processes=1, skip_outliers=0, save_dataframe=True
        )
        self.assertGreater(len(cns), 0)
        self.assertEqual(cns.chromosome.unique().tolist(), ["chr2"])
        self.assertIn("chrom", rstr)  # valid stitched SEG string, no crash

    @pytest.mark.slow
    def test_cbs_smooth_single_bin_arm(self):
        """Issue #594: --smooth-cbs must not crash on a single-bin arm.

        DNAcopy's ``smooth.CNA`` aborts with "NA/NaN/Inf in foreign function
        call (arg 7)" when an arm has only one bin (e.g. a lone surviving chrY
        bin). Segmentation runs the R script once per arm, so the guard skips
        smoothing for arms too short to smooth; the lone bin must still come
        back as exactly one segment carrying its own log2. Tested directly (not
        via ``do_segmentation``), whose process pool would swallow the crash.
        """
        cnr = self.tas_cnr[self.tas_cnr["chromosome"] == "chr1"]
        arm = cnr.as_dataframe(cnr.data.head(1).reset_index(drop=True).copy())
        segarr = segmentation._do_segmentation(
            arm, "cbs", None, 0.0001, skip_outliers=0, smooth_cbs=True
        )
        self.assertEqual(len(segarr), 1)
        self.assertAlmostEqual(segarr["log2"].iloc[0], arm["log2"].iloc[0], places=6)


def _test_method(self, method):
    for cnr in (
        self.tas_cnr,
        # self.wgs_cnr
    ):
        cns, raw_str = segmentation.do_segmentation(
            cnr, method, processes=1, save_dataframe=True
        )
        self.assertGreater(len(cns), 0)
        self.assertGreater(len(raw_str), 0)
        # Parallel should produce the same results
        p_cns, p_raw_str = segmentation.do_segmentation(
            cnr, method, processes=2, save_dataframe=True
        )
        self.assertEqual(cns.data.shape, p_cns.data.shape)
        self.assertEqual(len(cns.meta), len(p_cns.meta))
        self.assertEqual(raw_str, p_raw_str)


if __name__ == "__main__":
    unittest.main(verbosity=2)
