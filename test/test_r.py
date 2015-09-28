#!/usr/bin/env python

"""Unit tests for CNVkit that require an R installation."""
from __future__ import absolute_import, division, print_function

import unittest

from cnvlib import segmentation


class RTests(unittest.TestCase):
    """Tests that depend on the R statistical environment."""

    def test_segment(self):
        # Each method
        for method in ("cbs", "flasso"):
                cns = segmentation.do_segmentation("formats/amplicon.cnr",
                                                   False, method)
                self.assertTrue(len(cns) > 0)
        # With the R dataframe
        cns, dframe = segmentation.do_segmentation("formats/amplicon.cnr", True,
                                                   "flasso", 0.01)
        self.assertTrue(len(cns) > 0)
        self.assertTrue(len(dframe) > 0)



if __name__ == '__main__':
    unittest.main(verbosity=2)
