"""I/O for UCSC 'genePred' formats.

The formats are tabular, with no header row, and columns defined by the SQL
table definitions shown with each function. In alternative-splicing situations,
each transcript has a row in these tables.

Note: The parsers here load the gene information in each row and deduplicate
identical rows, but do not merge non-identical rows.

Generally the more-standard GFF or GTF would be used for this information; these
formats are essentially UCSC Genome Browser database dumps.

See:

- https://genome.ucsc.edu/FAQ/FAQformat.html#format9
- ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/
"""
from __future__ import absolute_import, division, print_function

import pandas as pd


def read_genepred(infile, exons=False):
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

    """
    raise NotImplementedError


def read_genepred_ext(infile, exons=False):
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

    """
    raise NotImplementedError


def read_refgene(infile, exons=False):
    """Gene predictions (extended) plus a "bin" column (e.g. refGene.txt)

    Same as genePredExt, but an additional first column of integers with the
    label "bin", which UCSC Genome Browser uses for optimization.
    """
    raise NotImplementedError


def read_refflat(infile, cds=False, exons=False):
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
        of the full transcript region (default).
    exons : bool
        Emit individual exonic regions for each gene instead of the full
        transcribed genomic region (default). Mutually exclusive with `cds`.

    """
    # ENH: choice of regions=('transcript', 'cds', 'exons') instead of flags?
    if cds and exons:
        raise ValueError("Arguments 'cds' and 'exons' are mutually exclusive")

    cols_shared = ['gene', 'accession', 'chromosome', 'strand']
    converters = None
    if exons:
        cols_rest = ['_start_tx', '_end_tx',  # Transcription
                     '_start_cds', '_end_cds',  # Coding region
                     '_exon_count', 'exon_starts', 'exon_ends']
        converters = {'exon_starts': _split_commas, 'exon_ends': _split_commas}
    elif cds:
        # Use CDS instead of transcription region
        cols_rest = ['_start_tx', '_end_tx',
                     'start', 'end',
                     '_exon_count', '_exon_starts', '_exon_ends']
    else:
        cols_rest = ['start', 'end',
                     '_start_cds', '_end_cds',
                     '_exon_count', '_exon_starts', '_exon_ends']
    colnames = cols_shared + cols_rest
    usecols = [c for c in colnames if not c.startswith('_')]
    # Parse the file contents
    dframe = pd.read_table(infile,  header=None, na_filter=False,
                           names=colnames, usecols=usecols,
                           dtype={c: str for c in cols_shared},
                           converters=converters)

    # Calculate values for output columns
    if exons:
        dframe = pd.DataFrame.from_records(_split_exons(dframe),
                                           columns=cols_shared + ['start', 'end'])
        dframe['start'] = dframe['start'].astype('int')
        dframe['end'] = dframe['end'].astype('int')

    return (dframe.assign(start=dframe.start - 1)
            .sort_values(['chromosome', 'start', 'end'])
            .reset_index(drop=True))


def _split_commas(field):
    return field.rstrip(',').split(',')


def _split_exons(dframe):
    """Split exons into individual rows."""
    for row in dframe.itertuples(index=False):
        shared = row[:4]
        for start, end in zip(row.exon_starts, row.exon_ends):
            yield shared + (start, end)
