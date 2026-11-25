"""I/O for UCSC 'genePred' formats.

The formats are tabular, with columns defined by the SQL table definitions
shown with each function. In alternative-splicing situations, each transcript
has a row in these tables. A leading '#'-prefixed line of field names, as
written by the UCSC Table Browser, is skipped. Coordinates are 0-based and
half-open, as in BED.

Every reader returns the columns gene, accession, chromosome, strand, start and
end, with one row per transcript, per exon, or per coding region (omitting
noncoding transcripts), and does not merge or deduplicate rows.

Generally the more-standard GFF or GTF would be used for this information; these
formats are essentially UCSC Genome Browser database dumps. Of these, refGene
is the one UCSC publishes per genome assembly (e.g. refGene.txt.gz).

See:

- https://genome.ucsc.edu/FAQ/FAQformat.html#format9
- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from pandas.core.frame import DataFrame


def read_genepred(infile: str, cds: bool = False, exons: bool = False) -> DataFrame:
    """Gene Predictions.

    ::

        table genePred
        "A gene prediction."
            (
            string  name;               "Name of gene"
            string  chrom;              "Chromosome name"
            char[1] strand;             "+ or - for strand"
            uint    txStart;            "Transcription start position"
            uint    txEnd;              "Transcription end position"
            uint    cdsStart;           "Coding region start"
            uint    cdsEnd;             "Coding region end"
            uint    exonCount;          "Number of exons"
            uint[exonCount] exonStarts; "Exon start positions"
            uint[exonCount] exonEnds;   "Exon end positions"
            )

    This table has no gene symbol, so the transcript name serves as both the
    gene and the accession. The `cds` and `exons` options are as in
    `read_refflat`.
    """
    return _read_genepred(infile, cds=cds, exons=exons)


def read_genepred_ext(infile: str, cds: bool = False, exons: bool = False) -> DataFrame:
    """Gene Predictions (Extended).

    The refGene table is an example of the genePredExt format.

    ::

        table genePredExt
        "A gene prediction with some additional info."
            (
            string name;        	"Name of gene (usually transcript_id from GTF)"
            string chrom;       	"Chromosome name"
            char[1] strand;     	"+ or - for strand"
            uint txStart;       	"Transcription start position"
            uint txEnd;         	"Transcription end position"
            uint cdsStart;      	"Coding region start"
            uint cdsEnd;        	"Coding region end"
            uint exonCount;     	"Number of exons"
            uint[exonCount] exonStarts; "Exon start positions"
            uint[exonCount] exonEnds;   "Exon end positions"
            int score;            	"Score"
            string name2;       	"Alternate name (e.g. gene_id from GTF)"
            string cdsStartStat; 	"enum('none','unk','incmpl','cmpl')"
            string cdsEndStat;   	"enum('none','unk','incmpl','cmpl')"
            lstring exonFrames; 	"Exon frame offsets {0,1,2}"
            )

    The gene is taken from name2 and the accession from name. The `cds` and
    `exons` options are as in `read_refflat`.
    """
    return _read_genepred(infile, after=_GENEPRED_EXT_TAIL, cds=cds, exons=exons)


def read_refgene(infile: str, cds: bool = False, exons: bool = False) -> DataFrame:
    """Gene predictions (extended) plus a "bin" column (e.g. refGene.txt)

    Same as genePredExt, but an additional first column of integers with the
    label "bin", which UCSC Genome Browser uses for optimization.
    """
    return _read_genepred(
        infile, before=("_bin",), after=_GENEPRED_EXT_TAIL, cds=cds, exons=exons
    )


def read_refflat(infile: str, cds: bool = False, exons: bool = False) -> DataFrame:
    """Gene predictions and RefSeq genes with gene names (e.g. refFlat.txt).

    This version of genePred associates the gene name with the gene prediction
    information. For example, the UCSC "refFlat" database lists HGNC gene names
    and RefSeq accessions for each gene, alongside the gene model coordinates
    for transcription region, coding region, and exons.

    ::

        table refFlat
        "A gene prediction with additional geneName field."
            (
            string  geneName;           "Name of gene as it appears in Genome Browser."
            string  name;               "Name of gene"
            string  chrom;              "Chromosome name"
            char[1] strand;             "+ or - for strand"
            uint    txStart;            "Transcription start position"
            uint    txEnd;              "Transcription end position"
            uint    cdsStart;           "Coding region start"
            uint    cdsEnd;             "Coding region end"
            uint    exonCount;          "Number of exons"
            uint[exonCount] exonStarts; "Exon start positions"
            uint[exonCount] exonEnds;   "Exon end positions"
            )

    Parameters
    ----------
    cds : bool
        Emit each gene's CDS region (coding and introns, but not UTRs) instead
        of the full transcript region (default). Noncoding transcripts, whose
        CDS is empty, are omitted.
    exons : bool
        Emit individual exonic regions for each gene instead of the full
        transcribed genomic region (default). Mutually exclusive with `cds`.

    """
    return _read_genepred(infile, before=("gene",), cds=cds, exons=exons)


#: The five columns genePredExt (and so refGene) adds after exonEnds
_GENEPRED_EXT_TAIL = (
    "_score",
    "gene",
    "_cds_start_stat",
    "_cds_end_stat",
    "_exon_frames",
)


def _read_genepred(
    infile: str,
    *,
    before: tuple[str, ...] = (),
    after: tuple[str, ...] = (),
    cds: bool = False,
    exons: bool = False,
) -> DataFrame:
    """Read one of UCSC's genePred-family tables.

    Every format in the family is genePred's ten columns with extra columns
    before, after, or both: `before` names those preceding the transcript
    name, `after` those following exonEnds. A name starting with "_" is parsed
    but not returned. Where neither supplies a "gene" column, the transcript
    name serves as the gene.
    """
    # ENH: choice of regions=('transcript', 'cds', 'exons') instead of flags?
    if cds and exons:
        raise ValueError("Arguments 'cds' and 'exons' are mutually exclusive")

    # Each mode takes the emitted region from a different pair of columns, or
    # from the exon lists when `exons` is set
    if exons:
        region = ("_start_tx", "_end_tx", "_start_cds", "_end_cds")
    elif cds:
        region = ("_start_tx", "_end_tx", "start", "end")
    else:
        region = ("start", "end", "_start_cds", "_end_cds")
    colnames = [
        *before,
        "accession",
        "chromosome",
        "strand",
        *region,
        "_exon_count",
        "_exon_starts",
        "_exon_ends",
        *after,
    ]
    cols_shared = ["gene", "accession", "chromosome", "strand"]
    # Parse every column as text, and missing trailing fields as empty strings
    # rather than NaN, so that a header line no wider than the table reads as
    # a row and the width check below sees the gap. Coordinates are converted
    # once the header is dropped.
    dframe = pd.read_csv(
        infile, sep="\t", header=None, na_filter=False, names=colnames, dtype=str
    )
    # Given surplus fields, pandas makes the leading ones an index; given too
    # few, it leaves the trailing columns empty. Either way fields would land
    # under the wrong names, so refuse a table of the wrong width.
    if not isinstance(dframe.index, pd.RangeIndex):
        raise ValueError(f"Expected {len(colnames)} tab-separated columns, found more")
    # The UCSC Table Browser writes a '#'-prefixed line of field names first.
    # Only leading lines are skipped; a '#' anywhere else is data.
    is_header = dframe[colnames[0]].str.startswith("#").cummin()
    dframe = dframe[~is_header]
    if (dframe[colnames[-1]] == "").any():
        raise ValueError(f"Expected {len(colnames)} tab-separated columns, found fewer")
    if "gene" not in dframe.columns:
        dframe = dframe.assign(gene=dframe["accession"])
    if exons:
        dframe = pd.DataFrame.from_records(
            _split_exons(dframe[[*cols_shared, "_exon_starts", "_exon_ends"]]),
            columns=[*cols_shared, "start", "end"],
        )
    dframe = dframe.astype({"start": "int64", "end": "int64"})
    if cds:
        # UCSC marks a noncoding transcript with cdsStart == cdsEnd
        is_coding = dframe["start"] != dframe["end"]
        if not is_coding.all():
            logging.info(
                "Skipping %d noncoding transcripts of %d",
                int((~is_coding).sum()),
                len(dframe),
            )
            dframe = dframe[is_coding]

    # UCSC stores these coordinates 0-based and half-open, as in BED
    return (
        dframe[[*cols_shared, "start", "end"]]
        .sort_values(["chromosome", "start", "end"])
        .reset_index(drop=True)
    )


def _split_commas(field):
    return field.rstrip(",").split(",")


def _split_exons(dframe):
    """Split exons into individual rows.

    The columns are the shared output columns, then the comma-separated exon
    starts and ends.
    """
    for row in dframe.itertuples(index=False):
        *shared, starts, ends = row
        for start, end in zip(_split_commas(starts), _split_commas(ends), strict=True):
            yield (*shared, start, end)
