"""Definitions for the core data structure, a copy number array."""
from __future__ import print_function

import sys

import numpy
from numpy.lib import recfunctions as rfn
from Bio.File import as_handle
from Bio._py3k import basestring, map, zip

from . import core, metrics, ngfrills, smoothing


class CopyNumArray(object):
    """An array of genomic intervals, treated like aCGH probes."""
    # http://docs.scipy.org/doc/numpy/user/basics.rec.html

    _dtype = (('chromosome', 'O'),
              ('start', 'i4'),
              ('end', 'i4'),
              ('gene', 'O'),
              ('coverage', 'f4'))
    _dtype_gc = ('gc', 'f4')
    _dtype_rmask = ('rmask', 'f4')
    _dtype_spread = ('spread', 'f4')
    _dtype_weight = ('weight', 'f4')
    _dtype_probes = ('probes', 'i4')

    def __init__(self, sample_id, extra_column_names=()):
        """Initialize an empty CopyNumArray with the requested columns."""
        xtra = []
        dtype = list(self._dtype)
        for col_name, col_dtype in (
            ('gc', self._dtype_gc),
            ('rmask', self._dtype_rmask),
            ('spread', self._dtype_spread),
            ('weight', self._dtype_weight),
            ('probes', self._dtype_probes)):
            if col_name in extra_column_names:
                xtra.append(col_name)
                dtype.append(col_dtype)
        self._xtra = tuple(xtra)
        self.data = numpy.empty(0, dtype=dtype)
        self.sample_id = sample_id

    @classmethod
    def from_array(cls, sample_id, array):
        cnarr = cls(sample_id, array.dtype.names[5:])
        cnarr.data = numpy.asarray(array, dtype=cnarr.data.dtype)
        return cnarr

    @classmethod
    def from_columns(cls, sample_id, **columns):
        these_cols = set(columns)
        required_cols = {'chromosome', 'start', 'end', 'gene', 'coverage'}
        optional_cols = {'gc', 'rmask', 'spread', 'weight', 'probes'}
        missing_cols = required_cols - these_cols
        if missing_cols:
            raise ValueError("Missing required column(s): " +
                             " ".join(sorted(missing_cols)))
        extra_cols = these_cols - required_cols
        bogus_cols = extra_cols - optional_cols
        if bogus_cols:
            raise ValueError("Unrecognized column name(s): " +
                             " ".join(sorted(bogus_cols)))
        cnarr = cls(sample_id, list(extra_cols))
        # XXX better way?
        table = list(zip(*[columns[key] for key in cnarr.data.dtype.names]))
        cnarr.data = numpy.asarray(table, dtype=cnarr.data.dtype)
        return cnarr

    @classmethod
    def from_rows(cls, sample_id, row_data, extra_keys=()):
        new_cna = cls(sample_id, extra_keys)
        if len(row_data) == 1:
            row_data = [tuple(row_data[0])]
        try:
            # Rows might be plain tuples
            new_array = numpy.asarray(row_data, dtype=new_cna.data.dtype)
        except ValueError:
            # "Setting void-array with object members using buffer"
            # All rows are numpy.ndarray
            new_array = rfn.stack_arrays(row_data, usemask=False,
                                         asrecarray=True, autoconvert=False)
            # print(new_array.dtype)

        new_cna.data = new_array
        return new_cna

    # Container-like behavior

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                (self.data == other.data).all())

    def __len__(self):
        return len(self.data)

    def __contains__(self, key):
        return (key in self.data.dtype.fields)

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        if isinstance(value, numpy.void) and isinstance(index, int):
            # Avoid "Setting void-array with object members using buffer"
            value = tuple(value)
        self.data[index] = value

    def __delitem__(self, index):
        return NotImplemented

    def __iter__(self):
        return iter(self.data)

    def __next__(self):
        return next(self.data)

    @property
    def coverage(self):
        return self.data['coverage']

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
    def gene(self):
        return self.data['gene']

    def labels(self):
        return (row2label(row) for row in self.data)

    def by_bin(self, bins):
        """Group rows by another CopyNumArray; trim row start/end to bin edges.

        Returns an iterable of (bin, CopyNumArray of overlapping cnarray rows))

        If a probe overlaps with a bin boundary, the probe start or end position
        is replaced with the bin boundary position. Probes outside any segments
        are skipped. This is appropriate for most other comparisons between
        CopyNumArray objects.
        """
        cn_by_chrom = dict(self.by_chromosome())
        for chrom, bin_rows in bins.by_chromosome():
            cn_rows = cn_by_chrom.get(chrom)
            if not cn_rows:
                continue

            # Traverse rows and bins together, matching up start/end points
            # cn_rows = iter(cn_rows)
            # curr_row = next(cn_rows)
            # row_start, row_end = curr_row['start'], curr_row['end']
            for bin_row in bin_rows:
                # Iterate over rows until we're clear of the bin
                # while True:
                #     if row_start < bin_row['start']:
                #         # Skip a gap in bin1
                #         row_start = bin_row['start']
                #     elif row_start > bin_row['start']:
                #         # Skip a gap in bin2
                #         bin_row['start'] = row_start
                #     # Row start points match now: emit up to the next endpoint
                #     bin_row['end'] = min(bin_row['end'], row_end)
                #     selected_rows.append(bin_row)
                #     if row_end > bin_row['end']:
                #         # Trim & hold bin2; take next bin1
                #         row_start = bin_row['end']
                #         break
                #     elif row_end < bin_row['end']:
                #         try:
                #             # Trim & hold bin1; take next bin2
                #             bin_row['start'] = row_end
                #             curr_row = next(bin_rows)
                #             row_start, row_end = curr_row['start'], curr_row['end']
                #             continue
                #         except StopIteration:
                #             # At the 3' telomere
                #             break
                #     elif row_end == bin_row['end']:
                #         try:
                #             # Take next of each of bin1 and bin2
                #             curr_row = next(bin_rows)
                #             row_start, row_end = curr_row['start'], curr_row['end']
                #         except StopIteration:
                #             # At the 3' telomere
                #             pass
                #         break

                # ENH - use a faster 1-pass implementation
                binned_rows = cn_rows.in_range(chrom, bin_row['start'],
                                               bin_row['end'], trim=True)

                yield bin_row, self.to_rows(binned_rows)

    def by_chromosome(self):
        """Iterate over probes grouped by chromosome name."""
        uniq_chrom, idx = numpy.unique(self.chromosome, True)
        sort_idx = idx.argsort()
        idx = idx.take(sort_idx)
        uniq_chrom = uniq_chrom.take(sort_idx)
        idx2 = numpy.concatenate((idx[1:], [len(self.data)]))
        for chrom, start, end in zip(uniq_chrom, idx, idx2):
            subarr = self.to_array(self.data[start:end])
            yield chrom, subarr

    def by_gene(self, ignore=('-', 'CGH', '.')):
        """Iterate over probes grouped by gene name.

        Emits pairs of (gene name, CNA of rows with same name)

        Groups each series of intergenic bins as a 'Background' gene; any
        'Background' bins within a gene are grouped with that gene.
        Bins with names in `ignore` are treated as 'Background' bins, but retain
        their name.
        """
        curr_chrom = None
        curr_bg_rows_idx = []
        curr_gene_name = None
        curr_gene_rows_idx = []

        for i, row in enumerate(self.data):
            gene = str(row['gene'])
            if gene == 'Background' or gene in ignore:
                # This background *may* be in an intergenic region
                if curr_chrom != row['chromosome']:
                    # New chromosome (not intergenic): emit the BG & reset
                    if curr_bg_rows_idx:
                        yield ('Background',
                               self.to_array(self.data.take(curr_bg_rows_idx)))
                        curr_bg_rows_idx = []
                # Include this row in the current background
                curr_chrom = row['chromosome']
                curr_bg_rows_idx.append(i)
            elif gene == curr_gene_name:
                # Continue the current gene
                # Any "background" was intronic; include in the current gene
                if curr_bg_rows_idx:
                    curr_gene_rows_idx.extend(curr_bg_rows_idx)
                    curr_bg_rows_idx = []
                # Add this row to the current gene
                curr_gene_rows_idx.append(i)
            else:
                # New gene
                # Emit the last gene, if any
                if curr_gene_rows_idx:
                    yield (curr_gene_name,
                           self.to_array(self.data.take(curr_gene_rows_idx)))
                # Emit the subsequent background, if any
                if curr_bg_rows_idx:
                    yield ('Background',
                           self.to_array(self.data.take(curr_bg_rows_idx)))
                    # Reset BG trackers
                    curr_bg_rows_idx = []
                # Start tracking the new gene
                curr_gene_rows_idx = [i]
                curr_chrom = row['chromosome']
                curr_gene_name = gene

        # Remainder
        if curr_gene_rows_idx:
            yield (curr_gene_name,
                   self.to_array(self.data.take(curr_gene_rows_idx)))
        if curr_bg_rows_idx:
            yield 'Background', self.to_array(self.data.take(curr_bg_rows_idx))

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
        for i, row in enumerate(self.data):
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
                       self.to_array(self.data.take(curr_probes_idx)))
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
                           self.to_array(self.data.take(curr_probes_idx)))
                    # Reset/update trackers for the next group of probes
                    curr_segment, next_segment = next_segment, None
                    curr_probes_idx = [i]
            else:
                raise ValueError("Mismatch between probes and segments\n" +
                                    "Probe: %s\nSegment: %s"
                                    % (row2label(row), row2label(curr_segment)))
        # Emit the remaining probes
        yield curr_segment, self.to_array(self.data.take(curr_probes_idx))

    def center_all(self, peak=False):
        """Recenter coverage values to the autosomes' average (in-place)."""
        chr_x = core.guess_chr_x(self)
        chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
        mask_autosome = ((self.chromosome != chr_x) &
                            (self.chromosome != chr_y))
        mid = numpy.median(self.coverage[mask_autosome])
        mask_cvg = (mask_autosome &
                    (self.coverage >= mid - 1.1) &
                    (self.coverage <= mid + 1.1))
        if peak and sum(mask_cvg) > 210:
            # Estimate the location of peak density
            # hack: from a smoothed histogram -- enh: kernel density estimate
            x = self.coverage[mask_cvg]
            w = self['weight'][mask_cvg] if 'weight' in self else None
            resn = int(round(numpy.sqrt(len(x))))
            x_vals, x_edges = numpy.histogram(x, bins=8*resn, weights=w)
            xs = smoothing.smoothed(x_vals, resn)
            # DBG: Check the fit
            # ngfrills.echo("Centering: median", mid,
            #               ", mode", x_edges[numpy.argmax(xs)],
            #               ", resolution", x_edges[2] - x_edges[1])
            # from matplotlib import pyplot
            # _fig, ax = pyplot.subplots()
            # ax.plot(x_vals, c='k', alpha=.5)
            # ax.plot(xs, c='r', lw=2)
            # median_idx = x_edges.searchsorted(mid)
            # ax.axvspan(median_idx - 1, median_idx, color='teal', alpha=.3)
            # ax.vlines(numpy.argmax(xs), *ax.get_ylim(), colors='salmon')
            # pyplot.show()
            # --
            mid = x_edges[numpy.argmax(xs)]
        self.data['coverage'] -= mid

    def copy(self):
        """Create an independent copy of this object."""
        return self.to_rows(numpy.copy(self.data))

    def add_columns(self, **columns):
        """Create a new CNA, adding the specified extra columns to this CNA."""
        cols = {key: self.data[key]
                for key in self.data.dtype.names}
        cols.update(columns)
        return self.__class__.from_columns(self.sample_id, **cols)

    def drop_extra_columns(self):
        """Remove any optional columns from this CopyNumArray.

        Returns a new copy with only the core columns retained:
            log2 value, chromosome, start, end, bin name.
        """
        result = self.__class__(self.sample_id)
        result.data = rfn.drop_fields(self.data, self._xtra)
        return result

    def merge(self, other):
        """Combine this array's data with another CopyNumArray (in-place).

        Any optional columns must match between both arrays.
        """
        assert isinstance(other, CopyNumArray), (
            "Argument (type %s) is not a CopyNumArray instance" % type(other))
        self.data = numpy.concatenate((self.data, other.data))
        self.sort()

    def in_range(self, chrom, start=0, end=None, trim=False):
        """Get the CopyNumArray portion within the given genomic range.

        If trim=True, include bins straddling the range boundaries, and trim
        the bins endpoints to the boundaries.
        """
        rows = self.data[self.chromosome == chrom]
        if not len(rows):
            raise ValueError("Chromosome %s is not in this probe set" % chrom)
        if start:
            if trim:
                # Include all rows overlapping the start point
                rows = rows[rows['end'].searchsorted(start, 'right'):]
                # Update 5' endpoints to the boundary
                rows['start'][rows['start'] < start] = start
            else:
                # Only rows entirely after the start point
                rows = rows[rows['start'].searchsorted(start):]
        if end:
            if trim:
                rows = rows[:rows['start'].searchsorted(end)]
                # Update 3' endpoints to the boundary
                rows['end'][rows['end'] > end] = end
            else:
                rows = rows[:rows['end'].searchsorted(end, 'right')]
        return self.to_rows(rows)

    def select(self, selector=None, **kwargs):
        """Take a subset of rows where the given condition is true.

        Arguments can be a function (lambda expression) returning a bool, which
        will be used to select True rows, and/or keyword arguments like
        gene="Background" or chromosome="chr7", which will select rows where the
        keyed field equals the specified value.
        """
        outrows = self.data
        if selector is not None:
            assert callable(selector)
            condition = numpy.apply_along_axis(selector, 0, outrows)
            outrows = numpy.extract(condition, outrows)
        for key, val in kwargs.items():
            assert key in self
            outrows = outrows[outrows[key] == val]
        return self.to_rows(outrows)

    def shuffle(self):
        """Randomize the order of bins in this array (in-place)."""
        numpy.random.seed(0xA5EED)  # For reproducible results
        order = numpy.arange(len(self.data))
        numpy.random.shuffle(order)
        self.data = self.data[order]
        return order

    def sort(self, key=None):
        """Sort the bins in this array (in-place).

        Optional argument 'key' is one of:

            - a function that computes a sorting key from a CopyNumArray row
            - a string identifier for an existing data column
            - a list/array/iterable of precomputed keys equal in length to the
              number of rows in this CopyNumArray.

        By default, bins are sorted by chromosomal coordinates.
        """
        if key is None:
            # Sort by chrom, then by start position
            chrom_keys = list(map(core.sorter_chrom, self.data['chromosome']))
            order = numpy.lexsort((self.data['start'], chrom_keys))
        else:
            # Sort by the given key, using a stable sort algorithm
            if isinstance(key, basestring):
                keys = self.data[key]
            elif callable(key):
                keys = list(map(key, self.data))
            else:
                if not len(key) == len(self):
                    raise ValueError("Sort key, as an array, must have the "
                                     "same length as the CopyNumArray to sort "
                                     "(%d vs. %d)." % (len(key), len(self)))
                keys = key
            order = numpy.argsort(keys, kind='mergesort')
        self.data = self.data.take(order)

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
        return self.to_rows(outrows)

    def to_array(self, array):
        return self.__class__.from_array(self.sample_id, array)

    def to_rows(self, rows):
        """Like from_rows, reusing this instance's metadata."""
        return self.__class__.from_rows(self.sample_id, rows, self._xtra)

    # I/O

    @classmethod
    def read(cls, infile, sample_id=None):
        """Parse a tabular table of coverage data from a handle or filename."""
        if sample_id is None:
            if isinstance(infile, basestring):
                sample_id = core.fbase(infile)
            else:
                sample_id = '<unknown>'
        with as_handle(infile) as handle:
            try:
                header = next(handle)
            except StopIteration:
                # Don't crash on empty files
                return cls(sample_id)
            # Build CNA...
            xtra = _sniff_xtra(header)
            cnarr = cls(sample_id, xtra)
            arr = numpy.loadtxt(handle, delimiter="\t", dtype=cnarr.data.dtype, ndmin=1)
        cnarr.data = arr
        return cnarr

    def write(self, outfile=sys.stdout):
        """Write coverage data to a file or handle in tabular format.

        This is similar to BED or BedGraph format, but with extra columns.

        To combine multiple samples in one file and/or convert to another
        format, see the 'export' subcommand.
        """
        colnames = ['chromosome', 'start', 'end', 'gene', 'log2']
        colnames.extend(self._xtra)
        rows = (list(map(str, row)) for row in self.data)
        with ngfrills.safe_write(outfile) as handle:
            header = '\t'.join(colnames) + '\n'
            handle.write(header)
            handle.writelines('\t'.join(row) + '\n' for row in rows)


def _sniff_xtra(header_line):
    colnames = tuple(map(str.strip, header_line.split('\t')))
    assert colnames[:5] == ('chromosome', 'start', 'end', 'gene', 'log2'), \
            colnames[:5]
    xtra = colnames[5:]
    assert all(x in ('gc', 'rmask', 'spread', 'weight', 'probes')
               for x in xtra)
    return xtra


def row2label(row):
    return "{}:{}-{}".format(row['chromosome'], row['start'], row['end'])
