"""Definitions for a generic array of genomic positions."""
from __future__ import print_function, absolute_import

import sys

import numpy as np
import pandas as pd

from . import core, ngfrills


def uniq(arr):
    """Because I don't know how to do this with Pandas yet."""
    # XXX see: pd.Categorical
    # return pd.Categorical(arr, ordered=True)
    prev = None
    for elem in arr:
        if elem != prev:
            yield elem
            prev = elem


# NB: Start by implementing all CNVkit features here, then split later
class GenomicArray(object):
    """An array of genomic intervals.

    Can represent most BED-like tabular formats with arbitrary additional
    columns: SEG, interval list, ...

    Required columns: chromosome, start
    """
    _required_columns = ("chromosome", "start", "end")
    # Extra columns, in order:
    # "gene", "log2", "gc", "rmask", "spread", "weight", "probes"

    def __init__(self, data_table, meta_dict=None):
        if not all(c in data_table.columns for c in
                   self._required_columns):
            raise ValueError("data table must have at least columns "
                             + repr(self._required_columns))
        self.data = data_table
        self.meta = (dict(meta_dict)
                     if meta_dict is not None and len(meta_dict)
                     else {})

    @staticmethod
    def row2label(row):
        return "{}:{}-{}".format(row['chromosome'], row['start'], row['end'])

    @classmethod
    def from_columns(cls, columns, meta_dict=None):
        """Create a new instance from column arrays, given by name."""
        # TODO - ensure columns are sorted properly
        table = pd.DataFrame.from_dict(columns)
        return cls(table, meta_dict)

    @classmethod
    def from_rows(cls, rows, columns=None, meta_dict=None):
        """Create a new instance from a list of rows, as tuples or arrays."""
        if columns is None:
            columns = cls._required_columns
        table = pd.DataFrame.from_records(rows, columns=columns)
        return cls(table, meta_dict)

    def as_columns(self, **columns):
        """Extract a subset of columns, reusing this instance's metadata."""
        return self.__class__.from_columns(columns, self.meta)
        # return self.__class__(self.data.loc[:, columns], self.meta.copy())

    def as_dataframe(self, dframe):
        return self.__class__(dframe.reset_index(drop=True), self.meta.copy())

    # def as_index(self, index):
    #     """Subset with fancy/boolean indexing; reuse this instance's metadata."""
    #     if isinstance(index, (int, slice)):
    #         return self.__class__(self.data.iloc[index], self.meta.copy())
    #     else:
    #         return self.__class__(self.data[index], self.meta.copy())

    def as_rows(self, rows):
        """Extract rows by indices, reusing this instance's metadata."""
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
              isinstance(index[0], (int, slice, tuple)) and
              isinstance(index[1], basestring)):
            # Row index, column index -> cell value
            return self.data.loc[index]
        elif isinstance(index, slice):
            # return self.as_dataframe(self.data.take(index))
            return self.as_dataframe(self.data[index])
        else:
            # Selected row indices or boolean array, probably
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
        """Assign to a portion of the data.
        """
        # self.data[index] = value
        if isinstance(index, int):
            self.data.iloc[index] = value
        elif isinstance(index, basestring):
            self.data[index] = value
        elif (isinstance(index, tuple) and
              len(index) == 2 and
              isinstance(index[0], (int, slice, tuple)) and
              isinstance(index[1], basestring)):
            self.data.loc[index] = value
        else:
            assert isinstance(index, slice) or len(index) > 0
            self.data[index] = value

    def __delitem__(self, index):
        return NotImplemented

    def __iter__(self):
        # return iter(self.data)
        return (row for idx, row in self.data.iterrows())

    __next__ = next
    # def __next__(self):
    #     return next(iter(self))

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

    # XXX troubled
    # def by_coords(self, other, mode='trim'):
    def by_bin(self, bins, mode='trim'):
        """Group rows by another CopyNumArray; trim row start/end to bin edges.

        Returns an iterable of (bin, GenomicArray of overlapping rows))

        modes are:  exclusive, trim, inclusive
            (when coordinates are on a boundary, what happens to the overlapped
            bins? drop, trim to size, include whole)

        default = 'trim': If a probe overlaps with a bin boundary, the probe
        start or end position is replaced with the bin boundary position. Probes
        outside any segments are skipped. This is appropriate for most other
        comparisons between CopyNumArray objects.
        """
        for chrom, bin_rows in bins.by_chromosome():
            try:
                cn_rows = self[self.chromosome == chrom]
            except KeyError:
                continue
            # Traverse rows and bins together, matching up start/end points
            for bin_row in bin_rows:
                binned_rows = cn_rows.in_range(chrom, bin_row['start'],
                                               bin_row['end'], trim=True)
                yield bin_row, self.as_rows(binned_rows)

    def by_chromosome(self):
        """Iterate over bins grouped by chromosome name."""
        for chrom in uniq(self.chromosome):
            yield chrom, self[self.chromosome == chrom]

    def coords(self, also=()):
        """Get plain coordinates of each bin: chromosome, start, end.

        With `also`, also include those columns.

        Example:

        >>> probes.coords(also=["name", "strand"])
        """
        cols = list(self._required_columns)
        if also:
            cols.extend(also)
        return self.data.loc[:, cols]

    def labels(self):
        return self.data.apply(self.row2label, axis=1)

    def in_range(self, chrom, start=0, end=None, trim=False):
        """Get the CopyNumArray portion within the given genomic range.

        If trim=True, include bins straddling the range boundaries, and trim
        the bins endpoints to the boundaries.
        """
        try:
            table = self.data[self.data['chromosome'] == chrom]
        except KeyError:
            raise KeyError("Chromosome %s is not in this probe set" % chrom)
        if start:
            if trim:
                # Include all rows overlapping the start point
                table = table[table.end.searchsorted(start, 'right'):]
                # Update 5' endpoints to the boundary
                table.start[table.start < start] = start
            else:
                # Only rows entirely after the start point
                table = table[table.start.searchsorted(start):]
        if end:
            if trim:
                table = table[:table.start.searchsorted(end)]
                # Update 3' endpoints to the boundary
                table.end[table.end > end] = end
            else:
                table = table[:table.end.searchsorted(end, 'right')]
        return self.as_dataframe(table)

    # Modification

    def add_array(self, other):
        """Combine this array's data with another CopyNumArray (in-place).

        Any optional columns must match between both arrays.
        """
        if not isinstance(other, self.__class__):
            raise ValueError("Argument (type %s) is not a %s instance"
                             % (type(other), self.__class__))
        # if not (self.data.index.names == other.data.index.names and
        #         list(self.data) == list(other.data)):
        #     raise ValueError("DataFrame indices or columns do not match")
        table = pd.concat([self.data.set_index(['chromosome', 'start']),
                           other.data.set_index(['chromosome', 'start'])])
        self.data = table.sort().reset_index()

    def copy(self):
        """Create an independent copy of this object."""
        return self.as_dataframe(self.data.copy())

    def add_columns(self, **columns):
        """Create a new CNA, adding the specified extra columns to this CNA."""
        # return self.as_dataframe(self.data.assign(**columns))
        result = self.copy()
        for key, values in columns.iteritems():
            result[key] = values
        return result

    def keep_columns(self, columns):
        """Extract a subset of columns, reusing this instance's metadata."""
        # XXX
        # required_cols = self.data.columns[:5]
        return self.__class__(self.data.loc[:, columns], self.meta.copy())

    def drop_extra_columns(self):
        """Remove any optional columns from this CopyNumArray.

        Returns a new copy with only the core columns retained:
            log2 value, chromosome, start, end, bin name.
        """
        table = self.data.loc[:, self._required_columns]
        return self.as_dataframe(table)

    # def reorder(self, key):
    #     """Apply a different ordering of chromosomes to this array.

    #     Key is either:

    #     - A function returning comparable (Python-sortable) values; or
    #     - A list of all chromosome values in the array, in the desired order.

    #     """
    #     if callable(key):
    #         chrom_keys = np.apply_along_axis(key, self.data["chromosome"], 0)
    #         new_order = np.argsort(chrom_keys)
    #         # TODO - sort self.data by new_order; reindex
    #     else:
    #         assert len(key) <= len(uniq(self.data["chromosome"]))
    #         # TODO - get indices of each listed chrom; sort into that order
    #         # NB: put any missing chroms at the end of the array, in their
    #         # current order

    def select(self, selector=None, **kwargs):
        """Take a subset of rows where the given condition is true.

        Arguments can be a function (lambda expression) returning a bool, which
        will be used to select True rows, and/or keyword arguments like
        gene="Background" or chromosome="chr7", which will select rows where the
        keyed field equals the specified value.
        """
        table = self.data
        if selector is not None:
            table = table[table.apply(selector, axis=1)]
        for key, val in kwargs.items():
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

    def sort(self, key=None):
        """Sort the bins in this array (in-place). Leaves chromosomes alone.
        """
        # if key is None:
        #     # Sort by chrom, then by start position
        #     chrom_keys = list(map(core.sorter_chrom, self.data['chromosome']))
        #     order = np.lexsort((self.data['start'], chrom_keys))
        # else:
        #     # Sort by the given key, using a stable sort algorithm
        #     if isinstance(key, basestring):
        #         keys = self.data[key]
        #     elif callable(key):
        #         keys = list(map(key, self.data))
        #     else:
        #         if not len(key) == len(self):
        #             raise ValueError("Sort key, as an array, must have the "
        #                              "same length as the CopyNumArray to sort "
        #                              "(%d vs. %d)." % (len(key), len(self)))
        #         keys = key
        #     order = np.argsort(keys, kind='mergesort')
        # self.data = self.data.take(order)
        table = self.data.set_index(['chromosome', 'start'])
        table.sort(inplace=True)
        self.data = table.reset_index()

    # I/O

    @classmethod
    def read(cls, infile, sample_id=None):
        if sample_id is None:
            if isinstance(infile, basestring):
                sample_id = core.fbase(infile)
            else:
                sample_id = '<unknown>'
        # Create a multi-index of genomic coordinates (like GRanges)
        table = pd.read_table(infile, na_filter=False,
                              # index_col=['chromosome', 'start']
        )
        # OR: Replace chromosome names with integers
        # table = pd.read_table(infile, na_filter=False)
        # chrom_names = uniq(table['chromosome'])
        # chrom_ids = np.arange(len(chrom_names))
        # chrom_id_col = np.zeros(len(table), dtype=np.int_)
        # for cn, ci in zip(chrom_names, chrom_ids):
        #     chrom_id_col[table['chromosome'] == cn] = ci
        # table['chromosome'] = chrom_id_col
        # table.set_index(['chromosome', 'start'], inplace=True)
        return cls(table, {"sample_id": sample_id})

    def write(self, outfile=sys.stdout):
        """Write coverage data to a file or handle in tabular format.

        This is similar to BED or BedGraph format, but with extra columns.

        To combine multiple samples in one file and/or convert to another
        format, see the 'export' subcommand.
        """
        with ngfrills.safe_write(outfile) as handle:
            self.data.to_csv(handle, index=False, sep='\t', float_format='%.6g')

