#!/usr/bin/env python

"""Command-line interface for CNVkit, the Copy Number Variation toolkit."""
from future import standard_library
standard_library.install_aliases()

import logging
from cnvlib import commands

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = commands.parse_args()
    args.func(args)
