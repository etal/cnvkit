#!/usr/bin/env python

"""Convert tumor segments and normal .cnr or reference .cnn to THetA input.

Follows the THetA segmentation import script but avoid repeating the pileups,
since we already have the mean depth of coverage in each target bin.

The options for average depth of coverage and read length do not matter
crucially for proper operation of THetA; increased read counts per bin simply
increase the confidence of THetA's results.
"""
from __future__ import division, print_function

import csv
import sys

import cnvlib


def calculate_theta_fields(seg, ref_rows, chrom_id, expect_depth, read_length):
    """Convert a segment's info to a row of THetA input.

    For the normal/reference bin count, take the mean of the bin values within
    each segment so that segments match between tumor and normal.
    """
    segment_size = seg["end"] - seg["start"]
    def logratio2count(log2_ratio):
        """Calculate a segment's read count from log2-ratio.

        Math:
            nbases = read_length * read_count
        and
            nbases = segment_size * read_depth
        where
            read_depth = read_depth_ratio * expect_depth

        So:
            read_length * read_count = segment_size * read_depth
            read_count = segment_size * read_depth / read_length
        """
        read_depth = (2 ** log2_ratio) * expect_depth
        read_count = segment_size * read_depth / read_length
        return int(round(read_count))

    tumor_count = logratio2count(seg["coverage"])
    ref_count = logratio2count(ref_rows["coverage"].mean())

    # e.g. "start_1_93709:end_1_19208166"
    row_id = "start_%d_%d:end_%d_%d" % (chrom_id, seg["start"],
                                        chrom_id, seg["end"])
    return {"#ID": row_id,
            "chrm": chrom_id,
            "start": seg["start"],
            "end": seg["end"],
            "tumorCount": tumor_count,
            "normalCount": ref_count}


def main(args):
    """Serialize tumor and reference segment coverages to THetA's input format.

    Columns: ID, chrm, start, end, tumorCount, normalCount.

    Chromosome IDs are integers 1 through 24.
    """
    tumor_segs = cnvlib.read(args.tumor_cns)
    ref_vals = cnvlib.read(args.reference)
    if not args.output:
        args.output = tumor_segs.sample_id + ".input"
    with open(args.output, 'w') as outfile:
        dwriter = csv.DictWriter(outfile,
                                 ["#ID", "chrm", "start", "end", "tumorCount",
                                  "normalCount"],
                                 dialect="excel-tab")
        dwriter.writeheader()
        prev_chrom = None
        chrom_id = 0
        for seg, ref_rows in ref_vals.by_segment(tumor_segs):
            if seg["chromosome"] != prev_chrom:
                chrom_id += 1
                prev_chrom = seg["chromosome"]
            fields = calculate_theta_fields(seg, ref_rows, chrom_id,
                                            args.coverage, args.read_length)
            dwriter.writerow(fields)
    print("Wrote", args.output, file=sys.stderr)


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("tumor_cns",
                    help="Tumor-sample segments from CNVkit (.cns)")
    AP.add_argument("reference",
                    help="""Normal-sample bin-level read depths (.cnr), or a
                    corresponding CNVkit reference (.cnn).""")
    AP.add_argument("-c", "--coverage", type=int, default=500,
                    help="Expected genome-wide depth of coverage.")
    AP.add_argument("-r", "--read-length", type=int, default=100,
                    help="Read length.")
    AP.add_argument("-o", "--output",
                    help="Output filename (THetA input file, *.input).")
    main(AP.parse_args())
