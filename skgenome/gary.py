"""Base class for an array of annotated genomic regions."""
from __future__ import print_function, absolute_import, division
from builtins import next, object, zip
from past.builtins import basestring

import logging
import warnings
from collections import OrderedDict

import numpy as np
import pandas as pd

from .chromsort import sorter_chrom
from .intersect import by_ranges, into_ranges, iter_ranges
from .merge import flatten, merge
from .rangelabel import to_label
from .subtract import subtract
from .subdivide import subdivide


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
                raise ValueError("data table must have at least columns %r; "
                                 "got %r" % (self._required_columns,
                                             tuple(data_table.columns)))
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
        try:
            out = self.from_rows(rows,
                                 columns=self.data.columns,
                                 meta_dict=self.meta)
        except AssertionError:
            columns = self.data.columns.tolist()
            firstrow = next(iter(rows))
            raise RuntimeError("Passed %d columns %r, but "
                               "%d elements in first row: %s",
                               len(columns), columns, len(firstrow), firstrow)
        return out

    # Container behaviour

    def __bool__(self):
        return bool(len(self.data))

    __nonzero__ = __bool__  # Py2.7

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
            # (as_indexer introduced before 0.18.1, deprecated 0.20.1)
            warnings.simplefilter("ignore", UserWarning)
            kwargs = dict(na=False)
            if pd.__version__ < "0.20.1":
                kwargs["as_indexer"] = True
            is_auto = self.chromosome.str.match(r"(chr)?\d+$", **kwargs)
        if not is_auto.any():
            # The autosomes, if any, are not named with plain integers
            return self
        if also:
            if isinstance(also, basestring):
                also = [also]
            for a_chrom in also:
                is_auto |= (self.chromosome == a_chrom)
        return self[is_auto]

    def by_arm(self, min_gap_size=1e5, min_arm_bins=50):
        """Iterate over bins grouped by chromosome arm (inferred)."""
        # ENH:
        # - Accept GArray of actual centromere regions as input
        #   -> find largest gap (any size) within cmere region, split there
        # - Cache centromere locations once found
        self.data.chromosome = self.data.chromosome.astype(str)
        for chrom, subtable in self.data.groupby("chromosome", sort=False):
            margin = max(min_arm_bins, int(round(.1 * len(subtable))))
            if len(subtable) > 2 * margin + 1:
                # Found a candidate centromere
                gaps = (subtable.start.values[margin+1:-margin] -
                        subtable.end.values[margin:-margin-1])
                cmere_idx = gaps.argmax() + margin + 1
                cmere_size = gaps[cmere_idx - margin - 1]
            else:
                cmere_idx = 0
                cmere_size = 0
            if cmere_idx and cmere_size >= min_gap_size:
                logging.debug("%s centromere at %d of %d bins (size %s)",
                             chrom, cmere_idx, len(subtable), cmere_size)
                p_arm = subtable.index[:cmere_idx]
                yield chrom, self.as_dataframe(subtable.loc[p_arm,:])
                q_arm = subtable.index[cmere_idx:]
                yield chrom, self.as_dataframe(subtable.loc[q_arm,:])
            else:
                # No centromere found -- emit the whole chromosome
                if cmere_idx:
                    logging.debug("%s: Ignoring centromere at %d of %d bins (size %s)",
                                  chrom, cmere_idx, len(subtable), cmere_size)
                else:
                    logging.debug("%s: Skipping centromere search, too small",
                                  chrom)
                yield chrom, self.as_dataframe(subtable)

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
        for bin_row, subrange in by_ranges(self.data, other.data,
                                           mode, keep_empty):
            if len(subrange):
                yield bin_row, self.as_dataframe(subrange)
            elif keep_empty:
                yield bin_row, self.as_rows(subrange)

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
        return self.data.apply(to_label, axis=1)

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
        results = iter_ranges(self.data, chrom, start, end, mode)
        return self.as_dataframe(next(results))

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
        table = pd.concat(iter_ranges(self.data, chrom, starts, ends, mode))
        return self.as_dataframe(table)

    def into_ranges(self, other, column, default, summary_func=None):
        """Re-bin values from `column` into the corresponding ranges in `other`.

        Match overlapping/intersecting rows from `other` to each row in `self`.
        Then, within each range in `other`, extract the value(s) from `column`
        in `self`, using the function `summary_func` to produce a single value
        if multiple bins in `self` map to a single range in `other`.

        For example, group SNVs (self) by CNV segments (other) and calculate the
        median (summary_func) of each SNV group's allele frequencies.

        Parameters
        ----------
        other : GenomicArray
            Bins to
        column : string
            Column name in `self` to extract values from.
        default
            Value to assign to indices in `other` that do not overlap any bins in
            `self`. Type should be the same as or compatible with the output
            field specified by `column`, or the output of `summary_func`.
        summary_func : callable, dict of string-to-callable, or None
            Specify how to reduce 1 or more `other` rows into a single value for
            the corresponding row in `self`.

                - If callable, apply to the `column` field each group of rows in
                  `other` column.
                - If a single-element dict of column name to callable, apply to that
                  field in `other` instead of `column`.
                - If None, use an appropriate summarizing function for the datatype
                  of the `column` column in `other` (e.g. median of numbers,
                  concatenation of strings).
                - If some other value, assign that value to `self` wherever there is
                  an overlap.

        Returns
        -------
        pd.Series
            The extracted and summarized values from `self` corresponding to
            other's genomic ranges, the same length as `other`.
        """
        if column not in self:
            logging.warning("No '%s' column available for summary calculation",
                            column)
            return pd.Series(np.repeat(default, len(other)))
        return into_ranges(self.data, other.data, column, default, summary_func)

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
        self.data = self.data.append(other.data, ignore_index=True)
        self.sort()

    def concat(self, others):
        """Concatenate several GenomicArrays, keeping this array's metadata.

        This array's data table is not implicitly included in the result.
        """
        table = pd.concat([otr.data for otr in others], ignore_index=True)
        result = self.as_dataframe(table)
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

    def keep_columns(self, colnames):
        """Extract a subset of columns, reusing this instance's metadata."""
        colnames = self.data.columns.intersection(colnames)
        return self.__class__(self.data.loc[:, colnames], self.meta.copy())

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

    def filter(self, func=None, **kwargs):
        """Take a subset of rows where the given condition is true.

        Parameters
        ----------
        func : callable
            A boolean function which will be applied to each row to keep rows
            where the result is True.
        **kwargs : string
            Keyword arguments like ``chromosome="chr7"`` or
            ``gene="Antitarget"``, which will keep rows where the keyed field
            equals the specified value.

        Return
        ------
        GenomicArray
            Subset of `self` where the specified condition is True.
        """
        table = self.data
        if func is not None:
            table = table[table.apply(func, axis=1)]
        for key, val in list(kwargs.items()):
            assert key in self
            table = table[table[key] == val]
        return self.as_dataframe(table)

    def shuffle(self):
        """Randomize the order of bins in this array (in-place)."""
        order = np.arange(len(self.data))
        np.random.seed(0xA5EED)
        np.random.shuffle(order)
        self.data = self.data.iloc[order]
        return order

    def sort(self):
        """Sort this array's bins in-place, with smart chromosome ordering."""
        sort_key = self.data.chromosome.apply(sorter_chrom)
        self.data = (self.data.assign(_sort_key_=sort_key)
                     .sort_values(by=['_sort_key_', 'start', 'end'],
                                  kind='mergesort')
                     .drop('_sort_key_', axis=1)
                     .reset_index(drop=True))

    def sort_columns(self):
        """Sort this array's columns in-place, per class definition."""
        extra_cols = []
        for col in self.data.columns:
            if col not in self._required_columns:
                extra_cols.append(col)
        sorted_colnames = list(self._required_columns) + sorted(extra_cols)
        assert len(sorted_colnames) == len(self.data.columns)
        self.data = self.data.reindex(columns=sorted_colnames)

    # Genome arithmetic

    def cut(self, other, combine=None):
        """Split this array's regions at the boundaries in `other`."""
        # TODO
        return NotImplemented

    def flatten(self, combine=None, split_columns=None):
        """Split this array's regions where they overlap."""
        return self.as_dataframe(flatten(self.data, combine=combine,
                                         split_columns=split_columns))

    def merge(self, bp=0, stranded=False, combine=None):
        """Merge adjacent or overlapping regions into single rows.

        Similar to 'bedtools merge'.
        """
        return self.as_dataframe(merge(self.data, bp, stranded, combine))

    def resize_ranges(self, bp, chrom_sizes=None):
        """Resize each genomic bin by a fixed number of bases at each end.

        Bin 'start' values have a minimum of 0, and `chrom_sizes` can
        specify each chromosome's maximum 'end' value.

        Similar to 'bedtools slop'.

        Parameters
        ----------
        bp : int
            Number of bases in each direction to expand or shrink each bin.
            Applies to 'start' and 'end' values symmetrically, and may be
            positive (expand) or negative (shrink).
        chrom_sizes : dict of string-to-int
            Chromosome name to length in base pairs. If given, all chromosomes
            in `self` must be included.
        """
        table = self.data
        limits = dict(lower=0)
        if chrom_sizes:
            limits['upper'] = self.chromosome.replace(chrom_sizes)
        table = table.assign(start=(table['start'] - bp).clip(**limits),
                             end=(table['end'] + bp).clip(**limits))
        if bp < 0:
            # Drop any bins that now have zero or negative size
            ok_size = table['end'] - table['start'] > 0
            logging.debug("Dropping %d bins with size <= 0", (~ok_size).sum())
            table = table[ok_size]
        # Don't modify the original
        return self.as_dataframe(table.copy())

    def squash(self, combine=None):
        """Combine some groups of rows, by some criteria, into single rows."""
        # TODO
        return NotImplemented

    def subdivide(self, avg_size, min_size=0, verbose=False):
        """Split this array's regions into roughly equal-sized sub-regions."""
        return self.as_dataframe(subdivide(self.data, avg_size, min_size,
                                           verbose))

    def subtract(self, other):
        """Remove the overlapping regions in `other` from this array."""
        return self.as_dataframe(subtract(self.data, other.data))

    def total_range_size(self):
        """Total number of bases covered by all (merged) regions."""
        if not len(self):
            return 0
        regions = merge(self.data, bp=1)
        return regions.end.sum() - regions.start.sum()

    def _get_gene_map(self):
        """Map unique gene names to their indices in this array.

        Returns
        -------
        OrderedDict
            An (ordered) dictionary of unique gene names and the data indices of
            their segments in the order of occurrence (genomic order).
        """
        if 'gene' not in self.data:
            return OrderedDict()

        genes = OrderedDict()
        for idx, genestr in self.data['gene'].iteritems():
            if pd.isnull(genestr):
                continue
            for gene in genestr.split(','):
                if gene not in genes:
                    genes[gene] = []
                genes[gene].append(idx)
        return genes
