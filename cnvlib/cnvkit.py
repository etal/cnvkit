#!/usr/bin/env python3
"""Command-line interface for CNVkit, the Copy Number Variation toolkit."""
import logging
from . import commands


def main():
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = commands.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
