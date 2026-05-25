#!/usr/bin/env python
"""Tests for the VariantArray class (cnvlib.vary)."""

import logging
import math
import unittest
import warnings

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


class VariantArrayTests(unittest.TestCase):
    """Tests for the VariantArray class."""

    def test_read(self):
        """Instantiate from a VCF file."""
        variants = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        self.assertGreater(len(variants), 0)

    def test_mirrored_baf_all_nan(self):
        """All-NaN frequencies mirror to all-NaN without warning (gh#407).

        np.nanmedian warns 'Mean of empty slice' on an all-NaN slice under
        older numpy; with pytest's filterwarnings=error that breaks the build
        (the 3.11-min job). _mirrored_baf must short-circuit before calling
        .median() on such a slice.
        """
        vals = pd.Series([np.nan, np.nan, np.nan])
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = vary._mirrored_baf(vals)
        self.assertEqual(len(out), len(vals))
        self.assertTrue(out.isna().all())
        # A partially-NaN slice still mirrors normally (median over real values)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            mixed = vary._mirrored_baf(pd.Series([np.nan, 0.2, 0.3]))
        self.assertTrue(math.isnan(mixed.iloc[0]))
        # median 0.25 < 0.5 -> below-half branch (0.5 - shift); both already minor
        self.assertAlmostEqual(mixed.iloc[1], 0.2)
        self.assertAlmostEqual(mixed.iloc[2], 0.3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
