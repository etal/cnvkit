#!/usr/bin/env python

"""Unit tests for the supporting library ngfrills."""

import unittest

import numpy
# from Bio._py3k import StringIO

import cnvlib
from cnvlib import ngfrills


NV2_TGT_BED = 'interval-nv2.target-267.bed'
NV2_ANTI_BED = 'interval-nv2.antitarget-15-150kb.bed'
EXO_TGT_BED = 'interval-exome.target-267.bed'
EXO_ANTI_BED = 'interval-exome.antitarget-9-90kb.bed'
ACCESS_BED = 'access-10kb.hg19.bed'


class NGTests(unittest.TestCase):
    """Tests for ngfrills functions."""


if __name__ == '__main__':
    unittest.main(verbosity=2)
