#!/usr/bin/env python

"""Unit tests for CNVkit that require an R installation."""
from __future__ import absolute_import, division, print_function

import unittest

import cnvlib
from cnvlib import segmentation


class RTests(unittest.TestCase):
    """Tests that depend on the R statistical environment."""

    def test_segment(self):
        cnarr = cnvlib.read("formats/amplicon.cnr")
        # Each method
        for method in ("cbs", "flasso"):
            cns = segmentation.do_segmentation(cnarr, method)
            self.assertGreater(len(cns), 0)
            # With the R dataframe
            cns, dframe = segmentation.do_segmentation(cnarr, "flasso", 0.01,
                                                       save_dataframe=True)
            self.assertGreater(len(cns), 0)
            self.assertGreater(len(dframe), 0)



if __name__ == '__main__':
    unittest.main(verbosity=2)
