"""I/O for Generic Feature Format (GFF).

Columns:

1. seqid/reference/seqname/chromosome
2. source
3. type/method/feature
4. start: in 1-based integer coordinates
5. end: in 1-based integer coordinates
6. score: float or '.' (for NA)
7. strand: [+-.?]
8. phase/frame: [012.]
9. attribute/group: string

Specs:

- http://gmod.org/wiki/GFF3
- http://gmod.org/wiki/GFF2
- http://mblab.wustl.edu/GTF2.html
- http://www.ensembl.org/info/website/upload/gff.html
- https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/subsets/SOFA.obo
"""

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pandas.core.frame import DataFrame
    from pandas.core.series import Series

# Symbol-bearing tags first: NCBI GFF3 puts the transcript or protein accession
# in Name and the symbol in gene, while GTF and GENCODE GFF3 put the stable ID
# in gene_id and the symbol in gene_name.
DEFAULT_TAGS = ("gene_name", "gene", "Name", "gene_id")


def read_gff(
    infile: str,
    tag: str | Sequence[str] = DEFAULT_TAGS,
    keep_type: str | None = None,
) -> DataFrame:
    """Read a GFF3/GTF/GFF2 file into a DataFrame.

    Works for all three formats because we only try to extract the gene name,
    at most, from column 9.

    Parameters
    ----------
    infile : filename or open handle
        Source file.
    tag : str or sequence of str
        GFF attribute tag, or tags in order of priority, from which to take
        gene names. Each tag is a literal attribute key, not a regular
        expression. Each row is labeled with the value of the first tag
        present in its column 9, wherever that tag appears in the column.
        Neither the GFF3 convention ("Name") nor the GTF one ("gene_id")
        reliably holds a gene symbol, so the default prefers "gene_name" and
        "gene", then falls back to "Name" and "gene_id".
    keep_type : str
        If specified, only keep rows with this value in the 'type' field
        (column 3). In GFF3, these terms are standardized in the Sequence
        Ontology Feature Annotation (SOFA).
    """
    colnames = [
        "chromosome",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attribute",
    ]
    coltypes = ["str", "str", "str", "int", "int", "str", "str", "str", "str"]
    dframe = pd.read_csv(
        infile,
        sep="\t",
        comment="#",
        header=None,
        na_filter=False,
        names=colnames,
        dtype=dict(zip(colnames, coltypes, strict=True)),
    )
    dframe = (
        dframe.assign(
            start=dframe.start - 1,
            score=dframe.score.replace(".", "nan").astype("float"),
        )
        .sort_values(["chromosome", "start", "end"])
        .reset_index(drop=True)
    )
    if keep_type:
        ok_type = dframe["type"] == keep_type
        logging.info(
            "Keeping %d '%s' / %d total records", ok_type.sum(), keep_type, len(dframe)
        )
        dframe = dframe[ok_type]
    tags = (tag,) if isinstance(tag, str) else tag
    genes = _extract_tags(dframe["attribute"], tags)
    if len(genes) and genes.isna().all():
        logging.warning(
            "No GFF records carry the attribute tag(s) %s; every record is labeled '-'",
            ", ".join(tags),
        )
    dframe["gene"] = genes.fillna("-").astype("str")
    return dframe


def _extract_tags(attributes: Series, tags: Sequence[str]) -> Series:
    """Extract the value of the first tag present in each attribute string."""
    values = pd.Series(None, index=attributes.index, dtype="object")
    for tag in tags:
        missing = values.isna()
        if not missing.any():
            break
        key = re.escape(tag)
        # A tag must begin column 9 (after at most one whitespace character) or
        # follow ";", "; " or two or more whitespace characters, so neither
        # "havana_gene=" nor "pseudogene 1" is the tag "gene", while GTF
        # attributes padded with extra spaces after ";" still parse. The
        # cost is that a double space inside a value also counts as a boundary:
        # Python's re has no variable-width lookbehind to tell the two apart.
        # The lookbehinds follow the tag rather than precede it so that the
        # regex engine can scan for the tag as a literal prefix, which is
        # several times faster.
        # A value cannot begin with a quote or ";", so an empty value such as
        # 'gene_name "";' or "gene=;" falls through to the next tag.
        rx = re.compile(
            rf"{key}(?<![^;\s]{key})(?<![^;\s]\s{key})"
            r'[= ]"?([^\s";]\S*?)"?(?:;|$)'
        )
        values.loc[missing] = attributes[missing].str.extract(rx, expand=False)
    return values
