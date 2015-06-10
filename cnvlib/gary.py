"""Definitions for a generic array of genomic positions."""
from __future__ import print_function, absolute_import

import sys

import numpy as np
import pandas as pd

from cnvlib import core, metrics, ngfrills, params


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
    # ENH: categorical index? currently MultiIndex

    def __init__(self, data_table, meta_dict=None):
        # if not (hasattr(data_table, "index") and
        #         hasattr(data_table.index, "names") and
        #         data_table.index.names == ["chromosome", "start"]):
        #     raise ValueError("data table must be indexed on "
        #                      "['chromosome', 'start']")
        if not all(c in data_table.columns for c in
                   ("chromosome", "start", "end")):
            raise ValueError("data table must have at least columns "
                             "['chromosome', 'start', 'end']")
        self.data = data_table
        self.meta = (dict(meta_dict)
                     if meta_dict is not None and len(meta_dict)
                     else {})

    @staticmethod
    def row2label(row):
        # return "{1}:{0}-{2}".format(row['end'], *row.name)
        return "{}:{}-{}".format(row['chromosome'], row['start'], row['end'])

    @classmethod
    def from_columns(cls, columns, meta_dict=None):
        """Create a new instance from column arrays, given by name."""
        # TODO - ensure columns are sorted properly
        table = pd.DataFrame.from_dict(columns)
        return cls(table, meta_dict)

    @classmethod
    def from_rows(cls, rows, extra_cols=(), meta_dict=None):
        """Create a new instance from a list of rows, as tuples or arrays."""
        cols = ['chromosome', 'start', 'end', 'gene', 'log2']
        if extra_cols:
            cols.extend(extra_cols)
        table = pd.DataFrame.from_records(rows, columns=cols)
        return cls(table, meta_dict)

    def as_columns(self, columns):
        """Extract a subset of columns, reusing this instance's metadata."""
        return self.__class__.from_columns(columns, self.meta)
        # return self.__class__(self.data.loc[:, columns], self.meta.copy())

    def as_dataframe(self, dframe):
        return self.__class__(dframe, self.meta.copy())

    # def as_index(self, index):
    #     """Subset with fancy/boolean indexing; reuse this instance's metadata."""
    #     if isinstance(index, (int, slice)):
    #         return self.__class__(self.data.iloc[index], self.meta.copy())
    #     else:
    #         return self.__class__(self.data[index], self.meta.copy())

    def as_rows(self, rows):
        """Extract rows by indices, reusing this instance's metadata."""
        return self.from_rows(rows,
                              extra_cols=self.data.columns[5:],
                              meta_dict=self.meta)

    # Container behaviour

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                (self.data.shape == other.data.shape) and
                (self.data == other.data).all().all())

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
            # Row index, column index
            return self.as_dataframe(self.data.loc[index])
        elif isinstance(index, slice):
            return self.as_dataframe(self.data.take(index))
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
        # return self.data.reset_index()['chromosome']
        return self.data['chromosome']

    @property
    def start(self):
        return self.data['start']

    @property
    def end(self):
        return self.data['end']

    @property
    def sample_id(self):
        # return self.data.reset_index()['start']
        return self.meta.get('sample_id')

    # Traversal

    # XXX troubled
    # def by_coords(self, other, mode='trim'):
    def by_bin(self, bins, mode='trim'):
        """Group rows by another CopyNumArray; trim row start/end to bin edges.

        Returns an iterable of (bin, CopyNumArray of overlapping cnarray rows))

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
        # table = self.data.reset_index()
        for chrom in uniq(self.chromosome):
            yield chrom, self[self.chromosome == chrom]

    # XXX hair: some genes overlap; some bins cover multiple genes
    #   -> option: whether to split gene names on commas
    def by_gene(self, ignore=('-', 'CGH', '.')):
        """Iterate over probes grouped by gene name.

        Emits pairs of (gene name, CNA of rows with same name)

        Groups each series of intergenic bins as a 'Background' gene; any
        'Background' bins within a gene are grouped with that gene.
        Bins with names in `ignore` are treated as 'Background' bins, but retain
        their name.
        """
        for gene in uniq(self.data['gene']):
            yield gene, self[self.data['gene'] == gene]

    # XXX superseded by by_genome_array? by_neighbors?
    def by_segment(self, segments):
        """Group cnarray rows by the segments that row midpoints land in.

        Returns an iterable of segments and rows grouped by overlap with each
        segment.

        Note that segments don't necessarily cover all probes (some near
        telo/centromeres may have been dropped as outliers during segmentation).
        These probes are grouped with the nearest segment, so the endpoint of
        the first/last probe may not match the corresponding segment endpoint.
        This is appropriate if the segments were obtained from this probe array.
        """
        curr_probes_idx = []
        segments = iter(segments)
        curr_segment = next(segments)
        next_segment = None
        for i, row in enumerate(self):
            probe_midpoint = 0.5 * (row['start'] + row['end'])
            if (row['chromosome'] == curr_segment['chromosome'] and
                curr_segment['start'] <= probe_midpoint <= curr_segment['end']):
                # Probe is within the current segment
                curr_probes_idx.append(i)

            elif row['chromosome'] != curr_segment['chromosome']:
                # Probe should be at the start of the next chromosome.
                # Find the matching segment.
                if next_segment is None:
                    next_segment = next(segments)

                # Skip any extra segments on the chromosome after the current
                # probes (e.g. probes are targets, trailing segments are based
                # on antitargets alone)
                while row['chromosome'] != next_segment['chromosome']:
                    try:
                        next_segment = next(segments)
                    except StopIteration:
                        raise ValueError("Segments are missing chromosome %r"
                                            % row['chromosome'])

                # Emit the current (completed) group
                yield (curr_segment,
                       self.as_dataframe(self.data.take(curr_probes_idx)))
                # Begin a new group of probes
                curr_segment, next_segment = next_segment, None
                curr_probes_idx = [i]

            elif row['start'] < curr_segment['start']:
                # Probe is near the start of the current chromosome, but we've
                # already seen another outlier here (rare/nonexistent case).
                # Group with the current (upcoming) segment.
                curr_probes_idx.append(i)

            elif row['end'] > curr_segment['end']:
                # Probe is at the end of an accessible region (e.g. p or q arm)
                # on the current chromosome.
                # Group this probe with whichever of the current or next
                # segments is closer.
                if next_segment is None:
                    next_segment = next(segments)
                if (next_segment['chromosome'] != row['chromosome']
                    or (next_segment['start'] - probe_midpoint) >
                    (probe_midpoint - curr_segment['end'])):
                    # The current segment is closer than the next. Group here.
                    curr_probes_idx.append(i)
                else:
                    # The next segment is closer. Emit the current group
                    # Begin a new group of probes
                    yield (curr_segment,
                           self.as_dataframe(self.data.take(curr_probes_idx)))
                    # Reset/update trackers for the next group of probes
                    curr_segment, next_segment = next_segment, None
                    curr_probes_idx = [i]
            else:
                raise ValueError("Mismatch between probes and segments\n" +
                                    "Probe: %s\nSegment: %s"
                                    % (self.row2label(row),
                                       self.row2label(curr_segment)))
        # Emit the remaining probes
        yield curr_segment, self.as_dataframe(self.data.take(curr_probes_idx))

    # s.data.loc['chr1'] -> all where chrom. index == "chr1"
    # s.data.loc['chr1', 93709] -> first row, as a Series
    # s.data.iloc[0], s.data.ix[0] -> the same
    # s.data.iloc[1:5] -> rows 2-6, keeping the index
    # s.data.sort(inplace=True); s.data.loc["chr15":"chr19"]
    # t.data[t.data['gene'].isin(["BRAF", "MET"])] -> select multiple genes

    # def coords(self, also=()):
    # def coords(self, strand=False, name=False):
    def coords(self, also=()):
        """Get plain coordinates of each bin: chromosome, start, end.

        With `also`, also include those columns.

        Example:

        >>> probes.coords(also=["name", "strand"])
        """
        cols = ["chromosome", "start", "end"]
        if also:
            cols.extend(also)
        # return self.data.reset_index().loc[:, cols]
        return self.data.loc[:, cols]

    def labels(self):
        return self.data.apply(self.row2label, axis=1)

    def in_range(self, chrom, start=0, end=None, trim=False):
        """Get the CopyNumArray portion within the given genomic range.

        If trim=True, include bins straddling the range boundaries, and trim
        the bins endpoints to the boundaries.
        """
        try:
            # table = self.data.loc[chrom]
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
                table = table[:table['end'].searchsorted(end, 'right')]
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
        self.data = pd.concat([self.data, other.data]).sort()

    # XXX mostly ok up to here
    def center_all(self, peak=False):
        """Recenter coverage values to the autosomes' average (in-place)."""
        # ideal: normalize to the total number of reads in this sample
        chr_x = core.guess_chr_x(self)
        chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
        mask_autosome = ((self.chromosome != chr_x) &
                            (self.chromosome != chr_y))
        mid = self.data['log2'][mask_autosome].median()
        mask_cvg = (mask_autosome &
                    (self.data['log2'] >= mid - 1.1) &
                    (self.data['log2'] <= mid + 1.1))
        if peak and sum(mask_cvg) > 210:
            # Estimate the location of peak density
            # hack: from a smoothed histogram -- enh: kernel density estimate
            x = self.data['log2'][mask_cvg]
            w = self['weight'][mask_cvg] if 'weight' in self else None
            resn = int(round(numpy.sqrt(len(x))))
            x_vals, x_edges = numpy.histogram(x, bins=8*resn, weights=w)
            xs = smoothing.smoothed(x_vals, resn)
            mid = x_edges[numpy.argmax(xs)]
        self.data['log2'] -= mid

    def copy(self):
        """Create an independent copy of this object."""
        return self.as_dataframe(self.data.copy())

    # XXX
    def add_columns(self, **columns):
        """Create a new CNA, adding the specified extra columns to this CNA."""
        cols = {key: self.data[key]
                for key in self.data.dtype.names}
        cols.update(columns)
        return self.__class__.from_columns(self.sample_id, **cols)

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
        table = self.data.loc[:, ('chromosome', 'start', 'end', 'gene',
                                  'log2')]
        return self.as_dataframe(table)

    def drop_low_coverage(self):
        """Drop bins with extremely low log2 coverage values.

        These are generally bins that had no reads mapped, and so were
        substituted with a small dummy log2 value to avoid divide-by-zero
        errors.
        """
        return self.as_dataframe(self.data[
                self.data['log2'] > params.NULL_LOG2_COVERAGE])

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

    def sort(self):
        """Sort the bins in this array (in-place). Leaves chromosomes alone.
        """
        self.data.sort(inplace=True)

    def squash_genes(self, ignore=('-', 'CGH', '.'), squash_background=False,
                     summary_stat=metrics.biweight_location):
        """Combine consecutive bins with the same targeted gene name.

        The `ignore` parameter lists bin names that not be counted as genes to
        be output.

        Parameter `summary_stat` is a function that summarizes an array of
        coverage values to produce the "squashed" gene's coverage value. By
        default this is the biweight location, but you might want median, mean,
        max, min or something else in some cases.

        Optional columns, if present, are dropped.
        """
        def squash_rows(name, rows):
            """Combine multiple rows (for the same gene) into one row."""
            chrom = core.check_unique(rows['chromosome'], 'chromosome')
            start = rows[0]['start']
            end = rows[-1]['end']
            cvg = summary_stat(rows['coverage'])
            outrow = [chrom, start, end, name, cvg]
            # Handle extra fields
            # ENH - no coverage stat; do weighted average as appropriate
            for xfield in ('gc', 'rmask', 'spread', 'weight'):
                if xfield in self:
                    outrow.append(summary_stat(rows[xfield]))
            if 'probes' in self:
                outrow.append(sum(rows['probes']))
            return tuple(outrow)

        outrows = []
        for name, subarr in self.by_gene(ignore):
            if name == 'Background' and not squash_background:
                outrows.extend(map(tuple, subarr))
            else:
                outrows.append(squash_rows(name, subarr))
        return self.as_rows(outrows)

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
        # table = table.pivot_table(index=['chromosome', 'start'])
        # OR?
        # table.set_index(['chromosome', 'start'], inplace=True)
        return cls(table, {"sample_id": sample_id})

    def write(self, outfile=sys.stdout):
        """Write coverage data to a file or handle in tabular format.

        This is similar to BED or BedGraph format, but with extra columns.

        To combine multiple samples in one file and/or convert to another
        format, see the 'export' subcommand.
        """
        with ngfrills.safe_write(outfile) as handle:
            self.data.to_csv(handle, sep='\t', float_format='%.6g')

