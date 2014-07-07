#!/usr/bin/env python

"""Command-line interface for CNVkit, the Copy Number Variation toolkit."""

from cnvlib import commands

if __name__ == '__main__':
    args = commands.parse_args()
    args.func(args)
