#!/usr/bin/env python
"""Unit tests for CNVkit that require an R installation."""
from __future__ import absolute_import, division, print_function

import unittest

import cnvlib
from cnvlib import segmentation


class RTests(unittest.TestCase):
    """Tests that depend on the R statistical environment."""

    def setUp(self):
        self.tas_cnr = cnvlib.read('formats/amplicon.cnr')
        self.wgs_cnr = cnvlib.read('formats/wgs-chr17.cnr')

    def test_cbs(self):
        _test_method(self, "cbs")

    def test_flasso(self):
        _test_method(self, "flasso")


def _test_method(self, method):
    for cnr in (self.tas_cnr,
                # self.wgs_cnr
               ):
        cns, raw_str = segmentation.do_segmentation(cnr, method, processes=1,
                                                    save_dataframe=True)
        self.assertGreater(len(cns), 0)
        self.assertGreater(len(raw_str), 0)
        # Parallel should produce the same results
        p_cns, p_raw_str = segmentation.do_segmentation(cnr, method,
                                                        processes=2,
                                                        save_dataframe=True)
        self.assertEqual(cns.data.shape, p_cns.data.shape)
        self.assertEqual(len(cns.meta), len(p_cns.meta))
        self.assertEqual(raw_str, p_raw_str)


if __name__ == '__main__':
    unittest.main(verbosity=2)
