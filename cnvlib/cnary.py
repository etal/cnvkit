"""CNVkit's core data structure, a copy number array."""
from __future__ import print_function, absolute_import, division
from builtins import map
from past.builtins import basestring

import logging

import numpy as np
import pandas as pd

from . import core, gary, metrics, params, smoothing


class CopyNumArray(gary.GenomicArray):
    """An array of genomic intervals, treated like aCGH probes.

    Required columns: chromosome, start, end, gene, log2

    Optional columns: gc, rmask, spread, weight, probes
    """
    _required_columns = ("chromosome", "start", "end", "gene", "log2")
    _required_dtypes = (str, int, int, str, float)
    # ENH: make gene optional
    # Extra columns, in order:
    # "gc", "rmask", "spread", "weight", "probes"

    def __init__(self, data_table, meta_dict=None):
        gary.GenomicArray.__init__(self, data_table, meta_dict)

    @property
    def log2(self):
        return self.data["log2"]

    @log2.setter
    def log2(self, value):
        self.data["log2"] = value

    @property
    def _chr_x_label(self):
        if 'chr_x' in self.meta:
            return self.meta['chr_x']
        chr_x = ('chrX' if self[0, 'chromosome'].startswith('chr')
                 else 'X')
        self.meta['chr_x'] = chr_x
        return chr_x

    @property
    def _chr_y_label(self):
        if 'chr_y' in self.meta:
            return self.meta['chr_y']
        chr_y = ('chrY' if self._chr_x_label.startswith('chr') else 'Y')
        self.meta['chr_y'] = chr_y
        return chr_y

    # More meta to store:
    #   is_sample_male = bool
    #   is_reference_male = bool
    #   file_path
    #   * invalidate 'chr_x' if .chromosome/['chromosome'] is set

    # Traversal

    # XXX hair: some genes overlap; some bins cover multiple genes
    #   -> option: whether to split gene names on commas
    def by_gene(self, ignore=params.IGNORE_GENE_NAMES):
        """Iterate over probes grouped by gene name.

        Emits pairs of (gene name, CNA of rows with same name)

        Groups each series of intergenic bins as a 'Background' gene; any
        'Background' bins within a gene are grouped with that gene.
        Bins with names in `ignore` are treated as 'Background' bins, but retain
        their name.
        """
        start_idx = end_idx = None
        for _chrom, subgary in self.by_chromosome():
            prev_idx = 0
            for gene in pd.unique(subgary.data['gene']):
                if not (gene == 'Background' or gene in ignore):
                    gene_idx = (subgary.data['gene'] == gene).nonzero()[0]
                    if not len(gene_idx):
                        logging.warn("Specified gene name somehow missing: %s",
                                     gene)
                        continue
                    start_idx = gene_idx[0]
                    end_idx = gene_idx[-1] + 1
                    if prev_idx < start_idx:
                        # Include intergenic regions
                        yield "Background", subgary.as_dataframe(
                                subgary.data.iloc[prev_idx:start_idx])
                    yield gene, subgary.as_dataframe(
                            subgary.data.iloc[start_idx:end_idx])
                    prev_idx = end_idx
            if prev_idx < len(subgary) - 1:
                # Include the telomere
                yield "Background", subgary.as_dataframe(
                        subgary.data.iloc[prev_idx:])

    # Manipulation

    def center_all(self, estimator=np.median, skip_low=False):
        """Recenter coverage values to the autosomes' average (in-place)."""
        est_funcs = {
            "mean": np.mean,
            "median": np.median,
            "mode": metrics.modal_location,
            "biweight": metrics.biweight_location,
        }
        if isinstance(estimator, basestring):
            if estimator in est_funcs:
                estimator = est_funcs[estimator]
            else:
                raise ValueError("Estimator must be a function or one of: %s"
                                 % ", ".join(map(repr, est_funcs)))
        table = (self.drop_low_coverage() if skip_low else self).autosomes()
        if table:
            self.data['log2'] -= estimator(table['log2'])

    def drop_low_coverage(self):
        """Drop bins with extremely low log2 coverage values.

        These are generally bins that had no reads mapped, and so were
        substituted with a small dummy log2 value to avoid divide-by-zero
        errors.
        """
        return self.as_dataframe(self.data[self.data['log2'] >
                                           params.NULL_LOG2_COVERAGE -
                                           params.MIN_REF_COVERAGE])

    def squash_genes(self, summary_func=metrics.biweight_location,
                     squash_background=False, ignore=params.IGNORE_GENE_NAMES):
        """Combine consecutive bins with the same targeted gene name.

        The `ignore` parameter lists bin names that not be counted as genes to
        be output.

        Parameter `summary_func` is a function that summarizes an array of
        coverage values to produce the "squashed" gene's coverage value. By
        default this is the biweight location, but you might want median, mean,
        max, min or something else in some cases.
        """
        def squash_rows(name, rows):
            """Combine multiple rows (for the same gene) into one row."""
            if len(rows) == 1:
                return tuple(rows.iloc[0])
            chrom = core.check_unique(rows.chromosome, 'chromosome')
            start = rows.iloc[0]['start']
            end = rows.iloc[-1]['end']
            cvg = summary_func(rows.log2)
            outrow = [chrom, start, end, name, cvg]
            # Handle extra fields
            # ENH - no coverage stat; do weighted average as appropriate
            for xfield in ('gc', 'rmask', 'spread', 'weight'):
                if xfield in self:
                    outrow.append(summary_func(rows[xfield]))
            if 'probes' in self:
                outrow.append(sum(rows['probes']))
            return tuple(outrow)

        outrows = []
        for name, subarr in self.by_gene(ignore):
            if name == 'Background' and not squash_background:
                outrows.extend(subarr.data.itertuples(index=False))
            else:
                outrows.append(squash_rows(name, subarr.data))
        return self.as_rows(outrows)

    # Chromosomal gender

    def shift_xx(self, male_reference=False):
        """Adjust chrX coverages (divide in half) for apparent female samples."""
        outprobes = self.copy()
        is_xx = self.guess_xx(male_reference=male_reference)
        if is_xx and male_reference:
            # Female: divide X coverages by 2 (in log2: subtract 1)
            outprobes[outprobes.chromosome == self._chr_x_label, 'log2'] -= 1.0
            # Male: no change
        elif not is_xx and not male_reference:
            # Male: multiply X coverages by 2 (in log2: add 1)
            outprobes[outprobes.chromosome == self._chr_x_label, 'log2'] += 1.0
            # Female: no change
        return outprobes

    def guess_xx(self, male_reference=False, verbose=True):
        """Guess whether a sample is female from chrX relative coverages.

        Recommended cutoff values:
            -0.5 -- raw target data, not yet corrected
            +0.5 -- probe data already corrected on a male profile
        """
        cutoff = 0.5 if male_reference else -0.5
        # ENH - better coverage approach: take Z-scores or rank of +1,0 or 0,-1
        # based on the available probes, then choose which is more probable
        rel_chrx_cvg = self.get_relative_chrx_cvg()
        if rel_chrx_cvg is None:
            return
        is_xx = (rel_chrx_cvg >= cutoff)
        if verbose:
            logging.info("Relative log2 coverage of X chromosome: %g "
                         "(assuming %s)",
                         rel_chrx_cvg, ('male', 'female')[is_xx])
        return is_xx

    def get_relative_chrx_cvg(self):
        """Get the relative log-coverage of chrX in a sample."""
        chromosome_x = self[self.chromosome == self._chr_x_label]
        if not len(chromosome_x):
            logging.warn("*WARNING* No %s found in probes; check the input",
                         self._chr_x_label)
            return
        autosomes = self.autosomes()
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

    def expect_flat_cvg(self, is_male_reference=None):
        """Get the uninformed expected copy ratios of each bin.

        Create an array of log2 coverages like a "flat" reference.

        This is a neutral copy ratio at each autosome (log2 = 0.0) and sex
        chromosomes based on whether the reference is male (XX or XY).
        """
        if is_male_reference is None:
            is_male_reference = not self.guess_xx(verbose=False)
        cvg = np.zeros(len(self), dtype=np.float_)
        if is_male_reference:
            # Single-copy X, Y
            idx = np.asarray((self.chromosome == self._chr_x_label) |
                             (self.chromosome == self._chr_y_label))
        else:
            # Y will be all noise, so replace with 1 "flat" copy
            idx = np.asarray(self.chromosome == self._chr_y_label)
        cvg[idx] = -1.0
        return cvg

    # Reporting

    def residuals(self, segments=None):
        """Difference in log2 value of each bin from its segment mean.

        If segments are just regions (e.g. RegionArray) with no log2 values
        precalculated, subtract the median of this array's log2 values within
        each region. If no segments are given, subtract each chromosome's
        median.

        """
        if not segments:
            resids = [subcna.log2 - subcna.log2.median()
                      for _chrom, subcna in self.by_chromosome()]
        elif "log2" in segments:
            resids = [subcna.log2 - seg.log2
                      for seg, subcna in self.by_ranges(segments)]
        else:
            resids = [subcna.log2 - subcna.log2.median()
                      for _seg, subcna in self.by_ranges(segments)]
        return np.concatenate(resids) if resids else np.array([])

    def guess_average_depth(self, segments=None, window=100):
        """Estimate the effective average read depth from variance.

        Assume read depths are Poisson distributed, converting log2 values to
        absolute counts. Then the mean depth equals the variance , and the average
        read depth is the estimated mean divided by the estimated variance.
        Use robust estimators (Tukey's biweight location and midvariance) to
        compensate for outliers and overdispersion.

        With `segments`, take the residuals of this array's log2 values from
        those of the segments to remove the confounding effect of real CNVs.

        If `window` is an integer, calculate and subtract a smoothed trendline
        to remove the effect of CNVs without segmentation (skipped if `segments`
        are given).

        See: http://www.evanmiller.org/how-to-read-an-unlabeled-sales-chart.html
        """
        # Try to drop allosomes
        cnarr = self.autosomes()
        if not len(cnarr):
            cnarr = self
        # Remove variations due to real/likely CNVs
        y_log2 = cnarr.residuals(segments)
        if window:
            y_log2 -= smoothing.smoothed(y_log2, window)
        # Guess Poisson parameter from absolute-scale values
        y = np.exp2(y_log2)
        # ENH: use weight argument to these stats
        loc = metrics.biweight_location(y)
        spread = metrics.biweight_midvariance(y, loc)
        scale = loc / spread**2
        return scale
