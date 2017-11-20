"""I/O for SEG format.

This is the output of DNAcopy segmentation, widely used to serialize segment
data.

The format is BED-like, but with a header row included and the
columns:

    - ID, "sampleName"
    - chrom, "chromosome"
    - loc.start, "start"
    - loc.end, "end"
    - num.mark, "nbrOfLoci" (optional)
    - seg.mean, "mean"

See: https://software.broadinstitute.org/software/igv/SEG
"""
from __future__ import absolute_import, division, print_function
from builtins import next
from past.builtins import basestring
#  from itertools import zip_longest
from future.moves.itertools import zip_longest

import collections
import csv
import logging
import math

import pandas as pd
from Bio.File import as_handle

LOG2_10 = math.log(10, 2)   # To convert log10 values to log2

# To catch exceptions from pandas versions 0.18 -- 0.20
CSV_ERRORS = (
    # Raised by the pandas 'python' CSV parser, at some point, I think
    csv.Error,
)
if hasattr(pd, 'errors'):
    # New in pandas 0.20
    CSV_ERRORS += (pd.errors.ParserError,)
if hasattr(pd.io.common, 'CParserError'):
    # Deprecated in pandas 0.20
    # Same as pandas.parser.CParserError in <0.20
    CSV_ERRORS += (pd.io.common.CParserError,)


def read_seg(infile, sample_id=None,
             chrom_names=None, chrom_prefix=None, from_log10=False):
    """Read one sample from a SEG file.

    Parameters
    ----------
    sample_id : string, int or None
        If a string identifier, return the sample matching that ID.  If a
        positive integer, return the sample at that index position, counting
        from 0. If None (default), return the first sample in the file.
    chrom_names : dict
        Map (string) chromosome IDs to names. (Applied before chrom_prefix.)
        e.g. {'23': 'X', '24': 'Y', '25': 'M'}
    chrom_prefix : str
        Prepend this string to chromosome names. (Usually 'chr' or None)
    from_log10 : bool
        Convert values from log10 to log2.

    Returns
    -------
    DataFrame of the selected sample's segments.
    """
    results = parse_seg(infile, chrom_names, chrom_prefix, from_log10)
    if isinstance(sample_id, int):
        # Select sample by index number
        for i, (_sid, dframe) in enumerate(results):
            if i == sample_id:
                return dframe
        else:
            raise IndexError("No sample index %d found in SEG file" % sample_id)

    elif isinstance(sample_id, basestring):
        # Select sample by name
        for sid, dframe in results:
            if sid == sample_id:
                return dframe
        else:
            raise IndexError("No sample ID '%s' found in SEG file" % sample_id)
    else:
        # Select the first sample
        sid, dframe = next(results)
        try:
            next(results)
        except StopIteration:
            pass
        else:
            logging.warning("WARNING: SEG file contains multiple samples; "
                            "returning the first sample '%s'", sid)
        return dframe


def parse_seg(infile, chrom_names=None, chrom_prefix=None, from_log10=False):
    """Parse a SEG file as an iterable of samples.

    Coordinates are automatically converted from 1-indexed to half-open
    0-indexed (Python-style indexing).

    Parameters
    ----------
    chrom_names : dict
        Map (string) chromosome IDs to names. (Applied before chrom_prefix.)
        e.g. {'23': 'X', '24': 'Y', '25': 'M'}
    chrom_prefix : str
        Prepend this string to chromosome names. (Usually 'chr' or None)
    from_log10 : bool
        Convert values from log10 to log2.

    Yields
    ------
    Tuple of (string sample ID, DataFrame of segments)
    """
    # Scan through any leading garbage to find the header
    with as_handle(infile) as handle:
        n_tabs = None
        for line in handle:
            n_tabs = line.count('\t')
            if n_tabs == 0:
                # Skip misc. R output (e.g. "WARNING...") before the header
                continue
            if n_tabs == 5:
                col_names = ["sample_id", "chromosome", "start", "end", "probes",
                             "log2"]
            elif n_tabs == 4:
                col_names = ["sample_id", "chromosome", "start", "end", "log2"]
            else:
                raise ValueError("SEG format expects 5 or 6 columns; found {}: {}"
                                .format(n_tabs + 1, line))
            break
        else:
            raise ValueError("SEG file contains no data")
        # Parse the SEG file contents
        try:
            dframe = pd.read_table(handle, names=col_names, header=None,
                                # * pandas.io.common.CParserError: Error
                                #   tokenizing data. C error: Calling
                                #   read(nbytes) on source failed. Try
                                #   engine='python'.
                                engine='python',
                                # * engine='c' only:
                                # na_filter=False,
                                # dtype={
                                #     'sample_id': 'str',
                                #     'chromosome': 'str',
                                #     'start': 'int',
                                #     'end': 'int',
                                #     'log2': 'float'
                                # },
                                )
            dframe['sample_id'] = dframe['sample_id'].astype("str")
            dframe['chromosome'] = dframe['chromosome'].astype("str")
        except CSV_ERRORS as err:
            raise ValueError("Unexpected dataframe contents:\n%s\n%s" %
                             (err, next(handle)))

    # Calculate values for output columns
    if chrom_names:
        dframe['chromosome'] = dframe['chromosome'].replace(chrom_names)
    if chrom_prefix:
        dframe['chromosome'] = dframe['chromosome'].apply(lambda c:
                                                          chrom_prefix + c)
    if from_log10:
        dframe['log2'] *= LOG2_10
    dframe['gene'] = "-"
    dframe['start'] -= 1
    keep_columns = dframe.columns.drop(['sample_id'])
    for sid, sample in dframe.groupby(by='sample_id', sort=False):
        yield sid, sample.loc[:, keep_columns]


def write_seg(dframe, sample_id=None, chrom_ids=None):
    """Format a dataframe or list of dataframes as SEG.

    To put multiple samples into one SEG table, pass `dframe` and `sample_id` as
    equal-length lists of data tables and sample IDs in matching order.
    """
    assert sample_id is not None
    if isinstance(dframe, pd.DataFrame):
        first = dframe
        first_sid = sample_id
        sids = dframes = None
    else:
        assert not isinstance(sample_id, basestring)
        dframes = iter(dframe)
        sids = iter(sample_id)
        first = next(dframes)
        first_sid = next(sids)

    if chrom_ids in (None, True):
        chrom_ids = create_chrom_ids(first)
    results = [format_seg(first, first_sid, chrom_ids)]
    if dframes is not None:
        # Unpack matching lists of data and sample IDs
        results.extend(
            format_seg(subframe, sid, chrom_ids)
            for subframe, sid in zip_longest(dframes, sids))
    return pd.concat(results)


def format_seg(dframe, sample_id, chrom_ids):
    assert dframe is not None
    assert sample_id is not None
    chroms = (dframe.chromosome.replace(chrom_ids) if chrom_ids
              else dframe.chromosome)
    rename_cols = {"log2": "seg.mean",
                   "start": "loc.start",
                   "end": "loc.end"}
    # NB: in some programs the "sampleName" column is labeled "ID"
    reindex_cols = ["ID", "chrom", "loc.start", "loc.end", "seg.mean"]
    if "probes" in dframe:
        rename_cols["probes"] = "num.mark" # or num_probes
        reindex_cols.insert(-1, "num.mark")
    return (dframe.assign(ID=sample_id,
                          chrom=chroms,
                          start=dframe.start + 1)
            .rename(columns=rename_cols)
            .reindex(columns=reindex_cols))


def create_chrom_ids(segments):
    """Map chromosome names to integers in the order encountered."""
    mapping = collections.OrderedDict(
        (chrom, i+1)
        for i, chrom in enumerate(segments.chromosome.drop_duplicates())
        if str(i + 1) != chrom)
    return mapping
