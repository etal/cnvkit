"""I/O for UCSC 'genePred' formats.

The formats are tabular, with no header row, and columns defined by the SQL
table definitions shown with each function. In alternative-splicing situations,
each transcript has a row in these tables.

Generally the more-standard GFF or GTF would be used for this information; these
formats are essentially UCSC Genome Browser database dumps.

See: https://genome.ucsc.edu/FAQ/FAQformat.html#format9
"""
from __future__ import absolute_import, division, print_function
from builtins import next
from past.builtins import basestring

import csv
import logging
import math

import numpy as np
import pandas as pd
from Bio.File import as_handle


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


def read_refflat(infile, cds=False, exons=False):
    """UCSC "flat" gene annotation format, e.g. refFlat.txt.

    Parse genes; merge those with same name and overlapping regions.

    "Gene Predictions and RefSeq Genes with Gene Names"

    UCSC refFlat: RefSeq genes, with additional name (as it appears in Genome
    Browser) field.  A version of genePred that associates the gene name with
    the gene prediction information. 


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

    """
    # ENH: choice of regions=('transcript', 'cds', 'exon') instead of bool args?
    if cds and exons:
        raise ValueError("Arguments 'cds' and 'exons' are mutually exclusive")

    cols_shared = ['gene', 'accession', 'chromosome', 'strand']
    if exons:
        cols_rest = ['_start_tx', '_end_tx',  # Transcription
                     '_start_cds', '_end_cds',  # Coding region
                     'exon_count', 'exon_starts', 'exon_ends']
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
    try:
        dframe = pd.read_table(infile,  header=None,
                               names=colnames, usecols=usecols,
                               dtype={c: str for c in cols_shared})
    except (pd.parser.CParserError, csv.Error) as err:
        raise ValueError("Unexpected dataframe contents:\n%s\n%s" %
                            (line, next(infile)))

    # Calculate values for output columns
    if exons:
        # TODO -- pandas trickery; then _merge_overlapping
        # es = dframe['exon_starts'].str.rstrip(',').str.split(',')
        # ee = dframe['exon_ends'].str.rstrip(',').str.split(',')
        # exons = list(zip(es, ee))
        # dframe.apply,applymap,...
        raise NotImplementedError
    else:
        dframe = dframe.sort_values(['chromosome', 'start', 'end'])
        dframe['start'] -= 1

    # NB: same gene name can appear on alt. contigs
    dframe = (dframe.groupby(by=['chromosome', 'strand', 'gene'],
                             as_index=False, group_keys=False, sort=False)
              .apply(_merge_overlapping))
    return dframe


def _merge_overlapping(dframe):
    """Merge overlapping regions within a group."""
    # Short-circuit the simple, common cases
    if len(dframe) == 1:
        return dframe
    if dframe['start'].nunique() == dframe['end'].nunique() == 1:
        return _squash_rows(dframe, accession=_join_strings)
    # Identify & enumerate (non)overlapping groups of rows
    overlap_sizes = dframe.end.cummax().values[:-1] - dframe.start.values[1:]
    non_overlapping = np.r_[False, (overlap_sizes <= 0)]
    # Squash rows within each non-overlapping group
    return (dframe.groupby(non_overlapping * np.arange(len(non_overlapping)),
                           as_index=False, group_keys=False, sort=False)
            .apply(_squash_rows,
                   accession=_join_strings,
                   start=pd.Series.min,
                   end=pd.Series.max))


def _squash_rows(dframe, **kwargs):
    """Reduce multiple rows into one, combining 'accession' field."""
    i = dframe.first_valid_index()
    row = dframe.loc[i:i, :]
    for key, combiner in kwargs.viewitems():
        row[key] = combiner(dframe[key])
    return row


def _join_strings(ser):
    return ','.join(ser.unique())

