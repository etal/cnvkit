"""Simple VCF I/O.

Read only coordinate info & store the remaining columns as unparsed strings.
Just enough functionality to extract a subset of samples and/or perform
bedtools-like operations on VCF records.
"""

import logging

import numpy as np
import pandas as pd
from Bio.File import as_handle


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
    colnames = ["chromosome", "start", "id", "ref", "alt", "qual", "filter", "end"]
    dtypes = {
        "chromosome": str,
        "start": int,
        "id": str,
        "ref": str,
        "alt": str,
        "filter": str,
    }
    table = pd.read_csv(
        infile,
        sep="\t",
        comment="#",
        header=None,
        na_filter=False,
        names=colnames,
        usecols=colnames,
        converters={"end": parse_end_from_info, "qual": parse_qual},
        dtype=dtypes,
    )
    table["start"] -= 1
    set_ends(table)
    logging.info("Loaded %d plain records", len(table))
    return table


def parse_end_from_info(info):
    """Parse the END position, if present, from an INFO field.

    Only a field keyed exactly "END" counts.  Keys merely ending in those
    characters are other fields entirely: CIEND, the confidence interval
    around END, is standard on imprecise structural variants.  Return -1
    where the record declares no usable END, leaving `set_ends` to fall
    back to the reference footprint.

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


def parse_qual(qual):
    """Parse a QUAL value as a number or NaN."""
    # ENH: only appy na_filter to this column
    if qual == ".":
        return np.nan
    return float(qual)


def set_ends(table) -> None:
    """Fill in each unusable 'end' from the reference allele's length.

    A record without INFO/END spans exactly the reference bases it quotes,
    starting at 'start'.  The alternate alleles do not enter into it: they
    describe what replaces that span, not how far it reaches, and a record
    carries one reference allele however many alternates it lists.

    An END below the record's own POS is unusable in the same way.  END is
    defined as the last position of a span running from POS, so a smaller
    value describes nothing, while the reference allele is mandatory and
    still says how far the record reaches.  Discarding it is what htslib
    does -- it warns "INFO/END=... is smaller than POS" and falls back to
    the same footprint -- so the text readers and `vcfio` agree here rather
    than each being self-consistent.  Honouring it would emit an interval
    ending before it starts, which every half-open overlap predicate in
    skgenome mis-answers silently, and which `read` deliberately does not
    repair for VCF: the coordinates are not written backwards, the
    declaration is wrong.  An END equal to POS is a one-base span and stands.
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
            "reference allele's length instead",
            int(backwards.sum()),
            len(table),
            table.chromosome.iat[first],
            int(table.start.iat[first]) + 1,
            int(table.end.iat[first]),
        )
    need_end_idx = ~declared | backwards
    ref_sz = table.loc[need_end_idx, "ref"].str.len()
    table.loc[need_end_idx, "end"] = table.loc[need_end_idx, "start"] + ref_sz
