"""Definitions for CNVkit's core data structure, a copy number array."""
from __future__ import print_function, absolute_import

import sys

import numpy as np
import pandas as pd

from . import core, gary, metrics, ngfrills
from .ngfrills import echo


class CopyNumArray(gary.GenomicArray):
    """An array of genomic intervals, treated like aCGH probes.

    Required columns: chromosome, start, end, gene, log2

    Optional columns: gc, rmask, spread, weight, probes
    """
    _required_columns = ("chromosome", "start", "end", "gene", "log2")
    # ENH: make gene optional
    # Extra columns, in order:
    # "gc", "rmask", "spread", "weight", "probes"

    def __init__(self, data_table, meta_dict=None):
        gary.GenomicArray.__init__(self, data_table, meta_dict)

    # Traversal

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
            # XXX TODO - include Background/ignore probes within a gene
            # see cnarray.py
            if not (gene == 'Background' or gene in ignore):
                yield gene, self[self.data['gene'] == gene]

    # XXX superseded by by_neighbors?
    def by_segment(self, segments):
        """Group rows by the segments that row midpoints land in.

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

    # Manipulation

    # XXX mostly ok up to here
    def center_all(self, peak=False):
        """Recenter coverage values to the autosomes' average (in-place)."""
        # ideal: normalize to the total number of reads in this sample
        chr_x = self.guess_chr_x()
        chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
        mask_autosome = ((self.chromosome != chr_x) &
                         (self.chromosome != chr_y))
        mid = self.data['log2'][mask_autosome].median()
        # mask_cvg = (mask_autosome &
        #             (self.data['log2'] >= mid - 1.1) &
        #             (self.data['log2'] <= mid + 1.1))
        # if peak and sum(mask_cvg) > 210:
        #     # Estimate the location of peak density
        #     # hack: from a smoothed histogram -- enh: kernel density estimate
        #     x = self.data['log2'][mask_cvg]
        #     w = self['weight'][mask_cvg] if 'weight' in self else None
        #     resn = int(round(np.sqrt(len(x))))
        #     x_vals, x_edges = np.histogram(x, bins=8*resn, weights=w)
        #     xs = smoothing.smoothed(x_vals, resn)
        #     mid = x_edges[np.argmax(xs)]
        self.data['log2'] -= mid

    def drop_low_coverage(self):
        """Drop bins with extremely low log2 coverage values.

        These are generally bins that had no reads mapped, and so were
        substituted with a small dummy log2 value to avoid divide-by-zero
        errors.
        """
        return self.as_dataframe(self.data[
                self.data['log2'] > params.NULL_LOG2_COVERAGE])

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
            cvg = summary_stat(rows['log2'])
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

    # Chromosomal gender
    # XXX refactor all this

    def shift_xx(self, male_reference=False, chr_x=None):
        """Adjust chrX coverages (divide in half) for apparent female samples."""
        if chr_x is None:
            chr_x = self.guess_chr_x()
        outprobes = self.copy()
        is_xx = self.guess_xx(chr_x=chr_x, male_reference=male_reference)
        if is_xx and male_reference:
            # Female: divide X coverages by 2 (in log2: subtract 1)
            outprobes['log2'][outprobes.chromosome == chr_x] -= 1.0
            # Male: no change
        elif not is_xx and not male_reference:
            # Male: multiply X coverages by 2 (in log2: add 1)
            outprobes['log2'][outprobes.chromosome == chr_x] += 1.0
            # Female: no change
        return outprobes

    def guess_chr_x(self):
        return ('chrX' if self[0, 'chromosome'].startswith('chr')
                else 'X')

    def guess_xx(self, male_reference=False, chr_x=None, verbose=True):
        """Guess whether a sample is female from chrX relative coverages.

        Recommended cutoff values:
            -0.5 -- raw target data, not yet corrected
            +0.5 -- probe data already corrected on a male profile
        """
        cutoff = 0.5 if male_reference else -0.5
        # ENH - better coverage approach: take Z-scores or rank of +1,0 or 0,-1
        # based on the available probes, then choose which is more probable
        rel_chrx_cvg = self.get_relative_chrx_cvg(chr_x=chr_x)
        is_xx = (rel_chrx_cvg >= cutoff)
        if verbose:
            echo("Relative log2 coverage of X chromosome:", rel_chrx_cvg,
                "(assuming %s)" % ('male', 'female')[is_xx])
        return is_xx

    def get_relative_chrx_cvg(self, chr_x=None):
        """Get the relative log-coverage of chrX in a sample."""
        if chr_x is None:
            chr_x = self.guess_chr_x()
        chromosome_x = self[self.chromosome == chr_x]
        if not len(chromosome_x):
            echo("*WARNING* No", chr_x, "found in probes; check the input")
            return
        chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
        autosomes = self[(self.chromosome != chr_x) &
                        (self.chromosome != chr_y)]
        auto_cvgs = autosomes['log2']
        x_cvgs = chromosome_x['log2']
        if 'probes' in self:
            # Weight segments by number of probes to ensure good behavior
            auto_sizes = autosomes['probes']
            x_sizes = chromosome_x['probes']
            # ENH: weighted median
            rel_chrx_cvg = (np.average(x_cvgs, weights=x_sizes) -
                            np.average(auto_cvgs, weights=auto_sizes))
        else:
            rel_chrx_cvg = np.median(x_cvgs) - np.median(auto_cvgs)
        return rel_chrx_cvg

    def expect_flat_cvg(self, is_male_reference=None, chr_x=None):
        """Get the uninformed expected copy ratios of each bin.

        Create an array of log2 coverages like a "flat" reference.

        This is a neutral copy ratio at each autosome (log2 = 0.0) and sex
        chromosomes based on whether the reference is male (XX or XY).
        """
        if chr_x is None:
            chr_x = self.guess_chr_x()
        if is_male_reference is None:
            is_male_reference = not self.guess_xx(chr_x=chr_x, verbose=False)
        chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
        cvg = np.zeros(len(self), dtype=np.float_)
        if is_male_reference:
            # Single-copy X, Y
            idx = np.asarray((self.chromosome == chr_x) |
                             (self.chromosome == chr_y))
        else:
            # Y will be all noise, so replace with 1 "flat" copy
            idx = np.asarray(self.chromosome == chr_y)
        cvg[idx] = -1.0
        return cvg

