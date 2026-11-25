#!/usr/bin/env python
"""Convert between tabular formats using scikit-genome I/O."""

import argparse
import logging
import sys

from .. import tabio


def argument_parsing():
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument(
        "infile",
        metavar="FILE",
        default=sys.stdin,
        nargs="?",
        help="Input filename. [Default: stdin]",
    )
    AP.add_argument(
        "-o", "--output", metavar="FILE", help="Output filename. [Default: stdout]"
    )
    AP.add_argument(
        "-f",
        "--from",
        metavar="FORMAT",
        default="auto",
        dest="in_fmt",
        help="Input format. [Default: auto-detect]",
    )
    AP.add_argument(
        "-t",
        "--to",
        metavar="FORMAT",
        default="bed4",
        dest="out_fmt",
        help="Output format. [Required]",
    )

    AP.add_argument(
        "--flatten",
        action="store_true",
        help="""Split overlapping regions at every boundary, giving each output
                    region the fields of the input regions covering it. Partly
                    overlapping regions yield fragments narrower than either
                    input, and any 'weight' or 'probes' field is divided among
                    them, so --merge is usually the better choice for bait
                    files. Not recommended with --exons.""",
    )
    AP.add_argument(
        "--merge",
        metavar="BASEPAIRS",
        nargs="?",
        type=int,
        const=1,
        help="""Merge overlapping regions with different names.
                    Recommended with --exons. The optional argument value is
                    the gap to bridge: 0 also merges bookended regions, and
                    negative values merge regions up to that many bases apart.
                    Any positive value merges only regions that genuinely
                    overlap. [Default: %(const)s]""",
    )
    # ENH combine --gff-type, --refflat-type

    AP_fmt = AP.add_argument_group("Format-specific options")
    AP_fmt.add_argument(
        "--gff-tag",
        help=f"""GFF attribute key to take gene names from, e.g. 'gene_name'.
                This is a literal key, not a regular expression. [Default: the
                first present of {", ".join(tabio.gff.DEFAULT_TAGS)}]""",
    )
    AP_fmt.add_argument(
        "--gff-type",
        help="""GFF record type (column 3) to use exclusively.
                    E.g. for GFF3 files, 'exon', 'gene', and other SOFA terms
                    can be used here.""",
    )
    AP_fmt.add_argument(
        "--refflat-type",
        choices=("exon", "cds", "transcript"),
        default="transcript",
        help="""Emit each exon, or the coding region, instead of the whole
                    transcribed region. (Applies to the UCSC refFlat, refGene,
                    genePred and genePredExt input formats.)""",
    )

    return AP.parse_args()


def skg_convert(args) -> None:
    logging.info(
        "Converting %s%s to %s",
        "input" if args.infile is sys.stdin else args.infile,
        " from " + args.in_fmt if args.in_fmt != "auto" else "",
        args.out_fmt,
    )

    if args.in_fmt == "auto":
        args.in_fmt = tabio.sniff_region_format(args.infile)
    # Format-specific input options
    kwargs = {}
    if args.in_fmt == "gff":
        if args.gff_tag:
            kwargs["tag"] = args.gff_tag
        if args.gff_type:
            kwargs["keep_type"] = args.gff_type
    elif args.in_fmt in ("genepred", "genepredext", "refflat", "refgene"):
        if args.refflat_type == "exon":
            kwargs["exons"] = True
        elif args.refflat_type == "cds":
            kwargs["cds"] = True
    regions = tabio.read(args.infile, args.in_fmt, **kwargs)

    # Post-processing
    if args.flatten:
        regions = regions.flatten()
    elif args.merge is not None:
        regions = regions.merge(bp=args.merge)

    tabio.write(regions, args.output, args.out_fmt)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    arguments = argument_parsing()
    skg_convert(arguments)


if __name__ == "__main__":
    main()
