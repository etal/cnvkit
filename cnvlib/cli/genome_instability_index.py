#!/usr/bin/env python3
"""Report the aneuploid fraction and segment count in each sample's genome.

For the reported CNA fraction, the denominator is the sequencing-accessible
footprint of the autosomes, and the numerator is the footprint of the regions
with non-neutral copy number (per the 'call' command).

The reported count is the absolute number of segments (or bins) with
non-neutral copy number.

These two values have been described as the Genome Instability Index (G2I),
where segmentation is performed only at the level of chromosome arms:
https://doi.org/10.1186/1755-8794-5-54

In this implementation, .cnr, .cns, or .call.cns files are accepted; .call.cns
is preferred, but if absolute copy number has not already been called, it will
be automatically inferred using the same thresholds as in the 'call' command.
"""
import argparse
import logging
import sys

from ..cmdutil import read_cna
from ..call import do_call


def argument_parsing():
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('cnv_files', nargs='+',
                    help="CNVkit .cnr or .cns filenames.")
    AP.add_argument("--diploid-parx-genome",
                    type=str,
                    help="Considers the given human genome's PAR of chromosome X as autosomal. Example: 'grch38'")
    AP.add_argument('-o', '--output',
                    type=argparse.FileType("w"), default=sys.stdout,
                    help="Output filename.")
    return AP.parse_args()


def cna_stats(cnarr, diploid_parx_genome):
    """Fraction of the mapped genome with altered copy number."""
    cnarr = cnarr.autosomes(diploid_parx_genome=diploid_parx_genome)
    denom = cnarr.total_range_size()
    if "cn" not in cnarr:
        cnarr = do_call(cnarr)
    aneuploid = cnarr[cnarr["cn"] != 2]
    numer = aneuploid.total_range_size()
    frac = numer / denom
    count = len(aneuploid)
    return frac, count


def genome_instability_index(args):
    print("Sample", "CNA_Fraction", "CNA_Count", sep='\t', file=args.output)
    for fname in args.cnv_files:
        cnarr = read_cna(fname)
        frac, count = cna_stats(cnarr, args.diploid_parx_genome)
        print(fname, frac, count, sep='\t', file=args.output)


def main():
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    arguments = argument_parsing()
    genome_instability_index(arguments)


if __name__ == '__main__':
    main()
