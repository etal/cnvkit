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

import cnvlib


def calculate_theta_fields(seg, ref_rows, chrom_id, expect_depth, read_length):
    # Convert log2-ratio to read-count-per-bin
    # calculated as:
    #   normalized_coverage = 2 ^ log2_ratio_coverage
    #   count = (segment_size * read_length)
    #           / (normalized_coverage * expected_depth)
    seg_size = seg["end"] - seg["start"]
    tumor_ratio = 2 ** seg["coverage"]
    tumor_count = int(seg_size * read_length
                      / (tumor_ratio * expect_depth))
    # Take the mean of reference bins within each segment
    ref_ratio = 2 ** ref_rows["coverage"].mean()
    ref_count = int(seg_size * read_length
                    / (ref_ratio * expect_depth))

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
    """Serialize tumor and reference segment coverages in THetA's input format.

    Columns: ID, chrm, start, end, tumorCount, normalCount.

    Chromosome IDs are integers 1 through 24.
    """
    tumor_segs = cnvlib.read(args.tumor_cns)
    ref_vals = cnvlib.read(args.reference)

    # Generate the THetA output file
    with open(args.output or tumor_segs.sample_id + ".input", 'w') as outfile:
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



if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("tumor_cns")
    AP.add_argument("reference",
                    help="""Normal-sample bin-level read depths, or a
                    corresponding CNVkit reference.""")
    AP.add_argument("-c", "--coverage", type=int, default=500,
                    help="Expected genome-wide depth of coverage.")
    AP.add_argument("-r", "--read-length", type=int, default=100,
                    help="Read length.")
    AP.add_argument("-o", "--output", help="Output filename.")
    main(AP.parse_args())
