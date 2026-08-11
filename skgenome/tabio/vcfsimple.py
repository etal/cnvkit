"""Simple VCF I/O.

Read only coordinate info & store the remaining columns as unparsed strings.
Just enough functionality to extract a subset of samples and/or perform
bedtools-like operations on VCF records.
"""

import logging

import numpy as np
import pandas as pd
from Bio.File import as_handle

from . import vcfspan


# TODO save VCF header (as string, the whole text block) in meta{header=}
def read_vcf_simple(infile):
    """Read VCF file without samples."""
    # ENH: Make all readers return a tuple (header_string, body_table)
    # ENH: usecols -- need to trim dtypes dict to match?
    header_lines = []
    with as_handle(infile, "r") as handle:
        for line in handle:
            if line.startswith("##"):
                header_lines.append(line)
            else:
                assert line.startswith("#CHR")
                header_line = line
                header_lines.append(line)
                break

        # Extract sample names from VCF header, keep as column names
        header_fields = header_line.split("\t")
        sample_ids = header_fields[9:]
        colnames = [
            "chromosome",
            "start",
            "id",
            "ref",
            "alt",
            "qual",
            "filter",
            "info",
            "format",
            *sample_ids,
        ]
        dtypes: dict[str, type] = dict.fromkeys(colnames, str)
        dtypes["start"] = int
        del dtypes["qual"]
        table = pd.read_csv(
            handle,
            sep="\t",
            header=None,
            na_filter=False,
            names=colnames,
            converters={"qual": parse_qual},
            dtype=dtypes,
        )
    # ENH: do things with filter, info
    table["start"] -= 1
    table["end"] = table["info"].apply(parse_end_from_info)
    set_ends(table)
    logging.info("Loaded %d plain records", len(table))
    return table


def read_vcf_sites(infile):
    """Read VCF contents into a DataFrame."""
    # The INFO column is read whole rather than converted straight into 'end':
    # the span needs SVLEN out of the same field, and a converter yields one
    # value per column. It is dropped again below, so the columns returned are
    # unchanged.
    colnames = ["chromosome", "start", "id", "ref", "alt", "qual", "filter", "info"]
    dtypes = {
        "chromosome": str,
        "start": int,
        "id": str,
        "ref": str,
        "alt": str,
        "filter": str,
        "info": str,
    }
    table = pd.read_csv(
        infile,
        sep="\t",
        comment="#",
        header=None,
        na_filter=False,
        names=colnames,
        usecols=colnames,
        converters={"qual": parse_qual},
        dtype=dtypes,
    )
    table["start"] -= 1
    table["end"] = table["info"].apply(parse_end_from_info)
    set_ends(table)
    del table["info"]
    logging.info("Loaded %d plain records", len(table))
    return table


def parse_end_from_info(info):
    """Parse the END position, if present, from an INFO field.

    Only a field keyed exactly "END" counts.  Keys merely ending in those
    characters are other fields entirely: CIEND, the confidence interval
    around END, is standard on imprecise structural variants.  Return -1
    where the record declares no usable END: it loses the maximum `set_ends`
    takes, leaving the record's span to its other terms.

    The INFO header is not consulted.  END is a reserved key, fixed by the
    spec at Number=1, Type=Integer, so its type is knowable from the key
    alone; htslib nonetheless declines to type an undeclared one ("assuming
    Type=String") and reports the footprint instead, so the two reader
    families part company on a file that omits the declaration.  That
    divergence is deliberate: reading what a careless writer meant, on
    files htslib will not, is the reason these readers exist.
    """
    for field in info.split(";"):
        key, _, value = field.partition("=")
        if key == "END" and value not in ("", "."):
            return int(value)
    return -1


def svlen_span_from_info(alt, info):
    """Reference bases INFO/SVLEN reaches for a record's ALT alleles, or 0.

    The sibling of `parse_end_from_info`, matching the key exactly for the
    same reason, and deciding which alleles SVLEN may speak for in the one
    place both reader families share.
    """
    for field in info.split(";"):
        key, _, value = field.partition("=")
        if key == "SVLEN":
            return vcfspan.svlen_span(alt.split(","), value.split(","))
    return 0


def parse_qual(qual):
    """Parse a QUAL value as a number or NaN."""
    # ENH: only appy na_filter to this column
    if qual == ".":
        return np.nan
    return float(qual)


def set_ends(table) -> None:
    """Set each record's 'end' to the furthest into the reference it reaches.

    Three things can say how far a record reaches, and the specification
    takes whichever says furthest: the reference allele it quotes, a declared
    INFO/END, and an INFO/SVLEN against an allele that replaces reference.
    The alternate alleles' own lengths do not enter into it -- they describe
    what replaces the span, not how far it runs -- and a record carries one
    reference allele however many alternates it lists.

    Taking the maximum settles the malformed cases without a branch for each.
    A record declaring no END carries the -1 sentinel, which loses to the
    reference allele; so does an END below the record's own POS, which states
    no span at all and which `read` could not repair even in principle, since
    for VCF an end below the start impugns the declaration rather than the
    orientation.  Only the warning below distinguishes that case, because a
    reader is entitled to know a declared value was unusable.

    An END shorter than the reference allele is not warned about: htslib is
    silent there too, and the maximum exists precisely so that such a record
    reads sensibly rather than exceptionally.

    Reads the 'start', 'ref', 'alt' and 'info' columns and rewrites 'end'.
    """
    declared = table.end != -1
    backwards = declared & (table.end <= table.start)
    if backwards.any():
        first = int(backwards.to_numpy().argmax())
        # POS and END as the file spells them, 1-based: `start` has already
        # been shifted by the time this runs, and `_describe_row`'s 0-based
        # rendering would name a coordinate absent from the user's own VCF.
        logging.warning(
            "Discarding INFO/END in %d of %d records for declaring an end "
            "below the record's own POS, first %s:%d with END=%d; using the "
            "record's own span instead",
            int(backwards.sum()),
            len(table),
            table.chromosome.iat[first],
            int(table.start.iat[first]) + 1,
            int(table.end.iat[first]),
        )
    starts = table.start.to_numpy()
    reach = starts + table.ref.str.len().to_numpy(dtype=int)
    # One vectorized pass decides whether any row needs the per-row parse.
    # The test is for a symbolic allele rather than for SVLEN itself: only a
    # symbolic allele can let SVLEN reach the reference, while SVLEN alone is
    # common on ordinary indels -- FreeBayes writes it on every one, so a
    # test for SVLEN fires on files that cannot contain a single span and
    # walks every row to learn it. Subscript rather than attribute: `.info`
    # is the DataFrame's own method.
    if table["alt"].str.contains("<", regex=False).any():
        spans = np.fromiter(
            (
                svlen_span_from_info(alt, info)
                for alt, info in zip(table["alt"], table["info"], strict=True)
            ),
            dtype=int,
            count=len(table),
        )
        reach = np.maximum(reach, starts + spans)
    table["end"] = np.maximum(table.end.to_numpy(), reach)
