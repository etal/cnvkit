"""Base class for an array of annotated genomic regions."""
from __future__ import print_function, absolute_import, division
from builtins import next, object, zip
from past.builtins import basestring

import logging
import sys
import warnings
from collections import OrderedDict

import numpy as np
import pandas as pd

from . import core


class GenomicArray(object):
    """An array of genomic intervals. Base class for genomic data structures.

    Can represent most BED-like tabular formats with arbitrary additional
    columns.
    """
    _required_columns = ("chromosome", "start", "end")
    _required_dtypes = (str, int, int)

    def __init__(self, data_table, meta_dict=None):
        # Validation
        if (data_table is None or
            (isinstance(data_table, (list, tuple)) and not len(data_table)) or
            (isinstance(data_table, pd.DataFrame) and not len(data_table.columns))
           ):
            data_table = self._make_blank()
        else:
            if not isinstance(data_table, pd.DataFrame):
                # Rarely if ever needed -- prefer from_rows, from_columns, etc.
                data_table = pd.DataFrame(data_table)
            if not all(c in data_table.columns for c in self._required_columns):
                raise ValueError("data table must have at least columns %r;"
                                 "got %r" % (self._required_columns,
                                             data_table.columns))
            # Ensure columns are the right type
            # (in case they've been automatically converted to the wrong type,
            # e.g. chromosome names as integers; genome coordinates as floats)
            if len(data_table):
                def ok_dtype(col, dt):
                    return isinstance(data_table[col].iat[0], dt)
            else:
                def ok_dtype(col, dt):
                    return data_table[col].dtype == np.dtype(dt)
            for col, dtype in zip(self._required_columns, self._required_dtypes):
                if not ok_dtype(col, dtype):
                    data_table[col] = data_table[col].astype(dtype)

        self.data = data_table
        self.meta = (dict(meta_dict)
                     if meta_dict is not None and len(meta_dict)
                     else {})

    @staticmethod
    def row2label(row):
        return "{}:{}-{}".format(row.chromosome, row.start, row.end)

    @classmethod
    def _make_blank(cls):
        """Create an empty dataframe with the columns required by this class."""
        spec = list(zip(cls._required_columns, cls._required_dtypes))
        try:
            arr = np.zeros(0, dtype=spec)
            return pd.DataFrame(arr)
        except TypeError as exc:
            raise TypeError("{}: {}".format(exc, spec))

    @classmethod
    def from_columns(cls, columns, meta_dict=None):
        """Create a new instance from column arrays, given as a dict."""
        table = pd.DataFrame.from_dict(columns)
        ary = cls(table, meta_dict)
        ary.sort_columns()
        return ary

    @classmethod
    def from_rows(cls, rows, columns=None, meta_dict=None):
        """Create a new instance from a list of rows, as tuples or arrays."""
        if columns is None:
            columns = cls._required_columns
        table = pd.DataFrame.from_records(rows, columns=columns)
        return cls(table, meta_dict)

    def as_columns(self, **columns):
        """Wrap the named columns in this instance's metadata."""
        return self.__class__.from_columns(columns, self.meta)
        # return self.__class__(self.data.loc[:, columns], self.meta.copy())

    def as_dataframe(self, dframe):
        """Wrap the given pandas dataframe in this instance's metadata."""
        return self.__class__(dframe.reset_index(drop=True), self.meta.copy())

    # def as_index(self, index):
    #     """Subset with fancy/boolean indexing; reuse this instance's metadata."""
    #     """Extract rows by indices, reusing this instance's metadata."""
    #     if isinstance(index, (int, slice)):
    #         return self.__class__(self.data.iloc[index], self.meta.copy())
    #     else:
    #         return self.__class__(self.data[index], self.meta.copy())

    def as_rows(self, rows):
        """Wrap the given rows in this instance's metadata."""
        return self.from_rows(rows,
                              columns=self.data.columns,
                              meta_dict=self.meta)

    # Container behaviour

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                self.data.equals(other.data))

    def __len__(self):
        return len(self.data)

    def __contains__(self, key):
        return key in self.data.columns

    def __getitem__(self, index):
        """Access a portion of the data.

        Cases:

        - single integer: a row, as pd.Series
        - string row name: a column, as pd.Series
        - a boolean array: masked rows, as_dataframe
        - tuple of integers: selected rows, as_dataframe
        """
        if isinstance(index, int):
            # A single row
            return self.data.iloc[index]
            # return self.as_dataframe(self.data.iloc[index:index+1])
        elif isinstance(index, basestring):
            # A column, by name
            return self.data[index]
        elif (isinstance(index, tuple) and
              len(index) == 2 and
              index[1] in self.data.columns):
            # Row index, column index -> cell value
            return self.data.loc[index]
        elif isinstance(index, slice):
            # return self.as_dataframe(self.data.take(index))
            return self.as_dataframe(self.data[index])
        else:
            # Iterable -- selected row indices or boolean array, probably
            try:
                if isinstance(index, type(None)) or len(index) == 0:
                    empty = pd.DataFrame(columns=self.data.columns)
                    return self.as_dataframe(empty)
            except TypeError:
                raise TypeError("object of type %r " % type(index) +
                                "cannot be used as an index into a " +
                                self.__class__.__name__)
            return self.as_dataframe(self.data[index])
            # return self.as_dataframe(self.data.take(index))

    def __setitem__(self, index, value):
        """Assign to a portion of the data."""
        if isinstance(index, int):
            self.data.iloc[index] = value
        elif isinstance(index, basestring):
            self.data[index] = value
        elif (isinstance(index, tuple) and
              len(index) == 2 and
              index[1] in self.data.columns):
            self.data.loc[index] = value
        else:
            assert isinstance(index, slice) or len(index) > 0
            self.data[index] = value

    def __delitem__(self, index):
        return NotImplemented

    def __iter__(self):
        return self.data.itertuples(index=False)

    __next__ = next

    @property
    def chromosome(self):
        return self.data['chromosome']

    @property
    def start(self):
        return self.data['start']

    @property
    def end(self):
        return self.data['end']

    @property
    def sample_id(self):
        return self.meta.get('sample_id')

    # Traversal

    def autosomes(self, also=()):
        """Select chromosomes w/ integer names, ignoring any 'chr' prefixes."""
        with warnings.catch_warnings():
            # NB: We're not using the deprecated part of this pandas method
            warnings.simplefilter("ignore", UserWarning)
            is_auto = self.chromosome.str.match(r"(chr)?\d+$",
                                                as_indexer=True, na=False)
        if not is_auto.any():
            # The autosomes, if any, are not named with plain integers
            return self
        if also:
            if isinstance(also, basestring):
                also = [also]
            for a_chrom in also:
                is_auto |= (self.chromosome == a_chrom)
        return self[is_auto]

    def by_chromosome(self):
        """Iterate over bins grouped by chromosome name."""
        # Workaround for pandas 0.18.0 bug:
        # https://github.com/pydata/pandas/issues/13179
        self.data.chromosome = self.data.chromosome.astype(str)
        for chrom, subtable in self.data.groupby("chromosome", sort=False):
            yield chrom, self.as_dataframe(subtable)

    def by_ranges(self, other, mode='inner', keep_empty=True):
        """Group rows by another GenomicArray's bin coordinate ranges.

        For example, this can be used to group SNVs by CNV segments.

        Bins in this array that fall outside the other array's bins are skipped.

        Parameters
        ----------
        other : GenomicArray
            Another GA instance.
        mode : string
            Determines what to do with bins that overlap a boundary of the
            selection. Possible values are:

            - ``inner``: Drop the bins on the selection boundary, don't emit them.
            - ``outer``: Keep/emit those bins as they are.
            - ``trim``: Emit those bins but alter their boundaries to match the
              selection; the bin start or end position is replaced with the
              selection boundary position.
        keep_empty : bool
            Whether to also yield `other` bins with no overlapping bins in
            `self`, or to skip them when iterating.

        Yields
        ------
        tuple
            (other bin, GenomicArray of overlapping rows in self)
        """
        chrom_lookup = dict(self.by_chromosome())
        for chrom, bin_rows in other.by_chromosome():
            if chrom in chrom_lookup:
                subranges = chrom_lookup[chrom]._iter_ranges(
                    None, bin_rows['start'], bin_rows['end'], mode)
                for bin_row, subrange in zip(bin_rows, subranges):
                    yield bin_row, subrange
            else:
                if keep_empty:
                    for bin_row in bin_rows:
                        yield bin_row, self.as_rows([])

    def coords(self, also=()):
        """Iterate over plain coordinates of each bin: chromosome, start, end.

        Parameters
        ----------
        also : str, or iterable of strings
            Also include these columns from `self`, in addition to chromosome,
            start, and end.

        Example, yielding rows in BED format:

        >>> probes.coords(also=["gene", "strand"])
        """
        cols = list(GenomicArray._required_columns)
        if also:
            if isinstance(also, basestring):
                cols.append(also)
            else:
                cols.extend(also)
        coordframe = self.data.loc[:, cols]
        return coordframe.itertuples(index=False)

    def labels(self):
        return self.data.apply(self.row2label, axis=1)

    def in_range(self, chrom=None, start=None, end=None, mode='inner'):
        """Get the GenomicArray portion within the given genomic range.

        Parameters
        ----------
        chrom : str or None
            Chromosome name to select. Use None if `self` has only one
            chromosome.
        start : int or None
            Start coordinate of range to select, in 0-based coordinates.
            If None, start from 0.
        end : int or None
            End coordinate of range to select. If None, select to the end of the
            chromosome.
        mode : str
            As in `by_ranges`: ``outer`` includes bins straddling the range
            boundaries, ``trim`` additionally alters the straddling bins'
            endpoints to match the range boundaries, and ``inner`` excludes
            those bins.

        Returns
        -------
        GenomicArray
            The subset of `self` enclosed by the specified range.
        """
        if isinstance(start, (int, np.int64, float, np.float64)):
            start = [int(start)]
        if isinstance(end, (int, np.int64, float, np.float64)):
            end = [int(end)]
        results = self._iter_ranges(chrom, start, end, mode)
        return next(results)

    def in_ranges(self, chrom=None, starts=None, ends=None, mode='inner'):
        """Get the GenomicArray portion within the specified ranges.

        Similar to `in_ranges`, but concatenating the selections of all the
        regions specified by the `starts` and `ends` arrays.

        Parameters
        ----------
        chrom : str or None
            Chromosome name to select. Use None if `self` has only one
            chromosome.
        starts : int array, or None
            Start coordinates of ranges to select, in 0-based coordinates.
            If None, start from 0.
        ends : int array, or None
            End coordinates of ranges to select. If None, select to the end of the
            chromosome. If `starts` and `ends` are both specified, they must be
            arrays of equal length.
        mode : str
            As in `by_ranges`: ``outer`` includes bins straddling the range
            boundaries, ``trim`` additionally alters the straddling bins'
            endpoints to match the range boundaries, and ``inner`` excludes
            those bins.

        Returns
        -------
        GenomicArray
            Concatenation of all the subsets of `self` enclosed by the specified
            ranges.
        """
        return self.concat(self._iter_ranges(chrom, starts, ends, mode))

    def _iter_ranges(self, chrom, starts, ends, mode):
        """Iterate through sub-ranges."""
        assert mode in ('inner', 'outer', 'trim')
        if chrom:
            assert isinstance(chrom, basestring)  # ENH: accept array?
            try:
                table = self.data[self.data['chromosome'] == chrom]
            except KeyError:
                raise KeyError("Chromosome %s is not in this probe set" % chrom)
        else:
            # Unsafe, but faster if we've already subsetted by chromosome
            table = self.data

        # Edge cases
        if not len(table):
            yield self.as_rows([])
            raise StopIteration

        if starts is None and ends is None:
            yield self.as_dataframe(table)
            raise StopIteration

        if starts is not None and len(starts):
            if mode == 'inner':
                # Only rows entirely after the start point
                start_idxs = table.start.searchsorted(starts)
            else:
                # Include all rows overlapping the start point
                start_idxs = table.end.searchsorted(starts, 'right')
        else:
            starts = np.zeros(len(ends) if ends is not None else 1,
                              dtype=np.int_)
            start_idxs = starts.copy()

        if ends is not None and len(ends):
            if mode == 'inner':
                end_idxs = table.end.searchsorted(ends, 'right')
            else:
                end_idxs = table.start.searchsorted(ends)
        else:
            end_idxs = np.repeat(len(table), len(starts))
            ends = [None] * len(starts)

        for start_idx, start_val, end_idx, end_val in zip(start_idxs, starts,
                                                          end_idxs, ends):
            subtable = table[start_idx:end_idx]
            if mode == 'trim':
                subtable = subtable.copy()
                # Update 5' endpoints to the boundary
                if start_val:
                    subtable.start = subtable.start.clip_lower(start_val)
                # Update 3' endpoints to the boundary
                if end_val:
                    subtable.end = subtable.end.clip_upper(end_val)
            yield self.as_dataframe(subtable)

    def match_to_bins(self, other, key, default=0.0, fill=False,
                      summary_func=np.median):
        """Take values of the other array at each of this array's bins.

        For example, group SNVs (other) by CNV segments (self) and calculate the
        median of each SNV group's allele frequencies.

        Parameters
        ----------
        other : GenomicArray
        key : str
            Column name containing the values to be summarized within each bin.
        default : object
            Value assigned to indices that fall outside the `other` array's bins,
            or chromosomes that appear in `self` but not `other`. Should be the
            same type as the column specified by `key`, or compatible.
        fill : bool
            Not used.
        summary_func : callable
            Function to summarize the `key` values in `other` within each of
            the ranges in `self`. A descriptive statistic; default median.

        Returns
        -------
        array
            The `key` column values in `other` corresponding to this array's bin
            locations, the same length as `self`.
        """
        def rows2value(rows):
            if len(rows) == 0:
                return default
            elif len(rows) == 1:
                return rows[0, key]
            else:
                return summary_func(rows[key])

        all_out_vals = [rows2value(other_rows) for _bin, other_rows in
                        other.by_ranges(self, mode='outer', keep_empty=True)]
        return np.asarray(all_out_vals)

    # Modification

    def add(self, other):
        """Combine this array's data with another GenomicArray (in-place).

        Any optional columns must match between both arrays.
        """
        if not isinstance(other, self.__class__):
            raise ValueError("Argument (type %s) is not a %s instance"
                             % (type(other), self.__class__))
        if not len(other.data):
            return self.copy()
        self.data = pd.concat([self.data, other.data])
        self.sort()

    def concat(self, others):
        """Concatenate several GenomicArrays, keeping this array's metadata.

        This array's data table is not implicitly included in the result.
        """
        result = self.as_dataframe(pd.concat([otr.data for otr in others]))
        result.sort()
        return result

    def copy(self):
        """Create an independent copy of this object."""
        return self.as_dataframe(self.data.copy())

    def add_columns(self, **columns):
        """Add the given columns to a copy of this GenomicArray.

        Parameters
        ----------
        **columns : array
            Keyword arguments where the key is the new column's name and the
            value is an array of the same length as `self` which will be the new
            column's values.

        Returns
        -------
        GenomicArray or subclass
            A new instance of `self` with the given columns included in the
            underlying dataframe.
        """
        # return self.as_dataframe(self.data.assign(**columns))
        result = self.copy()
        for key, values in columns.items():
            result[key] = values
        return result

    def keep_columns(self, columns):
        """Extract a subset of columns, reusing this instance's metadata."""
        return self.__class__(self.data.loc[:, columns], self.meta.copy())

    def drop_extra_columns(self):
        """Remove any optional columns from this GenomicArray.

        Returns
        -------
        GenomicArray or subclass
            A new copy with only the minimal set of columns required by the
            class (e.g. chromosome, start, end for GenomicArray; may be more for
            subclasses).
        """
        table = self.data.loc[:, self._required_columns]
        return self.as_dataframe(table)

    def select(self, selector=None, **kwargs):
        """Take a subset of rows where the given condition is true.

        Parameters
        ----------
        selector : callable
            A boolean function which will be applied to each row to select rows
            where the result is True.
        **kwargs : string
            Keyword arguments like ``chromosome="chr7"`` or
            ``gene="Background"``, which will select rows where the keyed field
            equals the specified value.

        Return
        ------
        GenomicArray
            Subset of `self` where the specified condition is True.
        """
        table = self.data
        if selector is not None:
            table = table[table.apply(selector, axis=1)]
        for key, val in list(kwargs.items()):
            assert key in self
            table = table[table[key] == val]
        return self.as_dataframe(table)

    def shuffle(self):
        """Randomize the order of bins in this array (in-place)."""
        np.random.seed(0xA5EED)  # For reproducible results
        order = np.arange(len(self.data))
        np.random.shuffle(order)
        self.data = self.data.iloc[order]
        return order

    def sort(self):
        """Sort this array's bins in-place, with smart chromosome ordering."""
        table = self.data.copy()
        table['SORT_KEY'] = self.chromosome.apply(core.sorter_chrom)
        table.sort_values(by=['SORT_KEY', 'start'], inplace=True)
        del table['SORT_KEY']
        self.data = table.reset_index(drop=True)

    def sort_columns(self):
        """Sort this array's columns in-place, per class definition."""
        extra_cols = []
        for col in self.data.columns:
            if col not in self._required_columns:
                extra_cols.append(col)
        sorted_colnames = list(self._required_columns) + sorted(extra_cols)
        assert len(sorted_colnames) == len(self.data.columns)
        self.data = self.data.reindex(columns=sorted_colnames)

    def _get_gene_map(self):
        """Map unique gene names to their indices in `self`.

        Returns
        -------
        OrderedDict
            An (ordered) dictionary of unique gene names and the data indices of
            their segments in the order of occurrence (genomic order).
        """
        if 'gene' not in self.data:
            return OrderedDict()

        genes = OrderedDict()
        for ix, row in self.data.iterrows():
            if pd.isnull(row.gene):
                continue
            for gene in row.gene.split(','):
                if gene not in genes:
                    genes[gene] = []
                genes[gene].append(ix)
        return genes
