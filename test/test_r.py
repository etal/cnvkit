#!/usr/bin/env python
"""Unit tests for CNVkit that require an R installation."""

import tempfile
import unittest
import logging

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
        self.wgs_cnr = cnvlib.read("formats/wgs-chr17.cnr")

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
            script = CBS_RSCRIPT % {
                "probes_fname": tmp.name,
                "sample_id": "test868",
                "threshold": 0.0001,
                "smooth_cbs": False,
            }
            with core.temp_write_text(script, mode="w+t") as script_fname:
                seg_out = core.call_quiet(
                    "Rscript", "--no-restore", "--no-environ", script_fname
                )
        # No exception means the weights stayed aligned with the kept probes.
        self.assertIn("test868", seg_out.decode())

    @pytest.mark.slow
    def test_cbs_nan_log2_in_memory(self):
        """CBS segments an in-memory .cnr with NaN log2 values (gh#881).

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
