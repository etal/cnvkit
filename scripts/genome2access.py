#!/usr/bin/env python

"""List the locations of accessible sequence regions in a FASTA file.

Inaccessible regions, e.g. telomeres and centromeres, are masked out with N in
the reference genome sequence; this script scans those to identify the
coordinates of the accessible regions (those between the long spans of N's).
"""

import argparse
import sys

from cnvlib.ngfrills import echo
# from cnvlib import sorter_chrom_at


def get_regions(fasta_fname):
    """Find accessible sequence regions (those not masked out with 'N')."""
    with open(fasta_fname) as infile:
        chrom = i = cursor = run_start = None
        for line in infile:
            if line.startswith('>'):
                # Emit the last chromosome's last run, if any
                if run_start is not None:
                    yield logme(chrom, run_start, cursor)
                else:
                    echo("\tChromosome ended with a telomere (N's)")
                # Start new chromosome
                chrom = line.split(None, 1)[0][1:]
                run_start = None
                cursor = 0
                echo("Scanning", chrom)
            else:
                for i, char in enumerate(line.rstrip()):
                    if char == 'N':
                        if run_start is not None:
                            # End of a run; emit the current region
                            yield logme(chrom, run_start, cursor + i)
                            run_start = None
                    else:
                        if run_start is None:
                            # Start of a new run of non-N characters
                            run_start = cursor + i
                            echo("\tStarted new run at", run_start)
                cursor += i + 1
        # Emit the last run if it's accessible (i.e. not a telomere)
        if run_start is not None:
            yield logme(chrom, run_start, cursor)


def logme(chrom, run_start, run_end):
    echo("\tEnded run at", run_end, "==>", run_end - run_start, 'b')
    return (chrom, run_start, run_end)


def join_regions(regions, min_gap_size):
    """Filter regions, joining those separated by small gaps."""
    regions = iter(regions)
    prev_chrom, prev_start, prev_end = next(regions)
    for chrom, start, end in regions:
        if chrom != prev_chrom:
            # New chromosome -- emit the remainder & reset things
            yield (prev_chrom, prev_start, prev_end)
            prev_chrom, prev_start, prev_end = chrom, start, end
        else:
            gap = start - prev_end
            assert gap > 0, ("Impossible gap between %s %d-%d and %d-%d (=%d)"
                             % (chrom, prev_start, prev_end, start, end, gap))
            if gap < min_gap_size:
                # Join with the previous region
                echo("\tJoining %s %d-%d and %d-%d (gap %d)"
                     % (chrom, prev_start, prev_end,
                        start, end, gap))
                prev_end = end
            else:
                # Keep the gap; emit the previous region as-is
                echo("\tKeeping gap between %s %d and %d (size %d)"
                     % (chrom, prev_end, start, gap))
                yield (prev_chrom, prev_start, prev_end)
                prev_chrom, prev_start, prev_end = chrom, start, end
    # If the last chromosome had no gaps, emit it too
    if prev_start == 0:
        yield (prev_chrom, prev_start, prev_end)


if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("fa_fname",
                    help="Genome FASTA file name")
    AP.add_argument("-s", "--min-gap-size", type=int, default=10000,
                    help="""Minimum gap size between accessible sequence
                    regions.  Regions separated by less than this distance will
                    be joined together.""")
    AP.add_argument("-o", "--output",
                    type=argparse.FileType('w'), default=sys.stdout,
                    help="Output file name")
    args = AP.parse_args()

    # Closes over args.output
    def write_row(chrom, run_start, run_end):
        args.output.write("%s\t%s\t%s\n"
                        % (chrom, run_start, run_end))
        args.output.flush()

    # for row in sorted(join_regions(get_regions(args.fa_fname),
    #                                args.min_gap_size),
    #                   key=sorter_chrom_at(0)):
    for row in join_regions(get_regions(args.fa_fname), args.min_gap_size):
        write_row(*row)

