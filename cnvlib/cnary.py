"""CNVkit's core data structure, a copy number array."""
from __future__ import print_function, absolute_import, division
from builtins import map
from past.builtins import basestring

import logging

import numpy as np
import pandas as pd

from . import core, gary, descriptives, params, smoothing


class CopyNumArray(gary.GenomicArray):
    """An array of genomic intervals, treated like aCGH probes.

    Required columns: chromosome, start, end, gene, log2

    Optional columns: gc, rmask, spread, weight, probes
    """
    _required_columns = ('chromosome', 'start', 'end', 'gene',
                         'log2', 'ratio')
    _required_dtypes = (str, int, int, str, float, float)
    # Extra columns:
    # "depth", "ratio",
    # "gc", "rmask", "spread", "weight", "probes"

    def __init__(self, data_table, meta_dict=None):
        if not isinstance(data_table, pd.DataFrame):
            # Empty but conformant table
            data_table = self._make_blank()
        # NB: 'log2' values are deprecated in favor of 'ratio' (absolute scale)
        if 'log2' in data_table.columns:
            # Every bin needs a log2 value; the others can be NaN
            d2 = data_table.dropna(subset=['log2'])
            if len(d2) < len(data_table):
                logging.warn("Dropped %d rows with missing log2 values",
                             len(data_table) - len(d2))
                data_table = d2
        else:
            # Shim to provide log2 column for functions that need it
            if 'ratio' in data_table.columns:
                key = 'ratio'
            elif 'depth' in data_table.columns:
                key = 'depth'
            else:
                raise ValueError("Missing 'log2'/'ratio'/'depth' column")
            data_table['log2'] = (np.log2(data_table[key])
                                  .replace(-np.inf, params.NULL_LOG2_COVERAGE))
        if 'ratio' not in data_table.columns:
            # Shim to load files created before 'ratio' was introduced
            if 'log2' in data_table.columns:
                data_table['ratio'] = np.exp2(data_table['log2'])
            elif 'depth' in data_table.columns:
                data_table['ratio'] = data_table['depth'] #/ data_table['depth'].median()
                # Likely .cnn: handle targets & antitargets like 'reference'
                bg_mask = data_table['gene'] == 'Background'
                bg_cna = CopyNumArray(data_table[bg_mask])
                fg_cna = CopyNumArray(data_table[~bg_mask])
                bg_cna.center_all(skip_low=False)
                fg_cna.center_all(skip_low=True)
                full_cna = fg_cna.concat([fg_cna, bg_cna])
                data_table = full_cna.data

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
        chr_x = ('chrX' if self.chromosome.iat[0].startswith('chr') else 'X')
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
    #   * invalidate 'chr_x' if .chromosome/['chromosome'] is set

    # Traversal

    def by_gene(self, ignore=params.IGNORE_GENE_NAMES):
        """Iterate over probes grouped by gene name.

        Groups each series of intergenic bins as a 'Background' gene; any
        'Background' bins within a gene are grouped with that gene.
        Bins with names in `ignore` are treated as 'Background' bins, but retain
        their name.

        Bins' gene names are split on commas to accommodate overlapping genes
        and bins that cover multiple genes.

        Return an iterable of pairs of: (gene name, CNA of rows with same name)
        """
        start_idx = end_idx = None
        for _chrom, subgary in self.by_chromosome():
            prev_idx = 0
            for gene, gene_idx in subgary._get_gene_map().items():
                if not (gene == 'Background' or gene in ignore):
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

    def center_all(self, estimator=np.median, skip_low=False, by_chrom=True):
        """Recenter coverage values to the autosomes' average (in-place)."""
        est_funcs = {
            "mean": np.mean,
            "median": np.median,
            "mode": descriptives.modal_location,
            "biweight": descriptives.biweight_location,
        }
        if isinstance(estimator, basestring):
            if estimator in est_funcs:
                estimator = est_funcs[estimator]
            else:
                raise ValueError("Estimator must be a function or one of: %s"
                                 % ", ".join(map(repr, est_funcs)))
        cnarr = (self.drop_low_coverage() if skip_low else self).autosomes()
        if cnarr:
            if by_chrom:
                values = np.array([estimator(subarr['log2'])
                                   for _c, subarr in cnarr.by_chromosome()
                                   if len(subarr)])
            else:
                values = cnarr['log2']
            self.data['log2'] -= estimator(values)

    def drop_low_coverage(self):
        """Drop bins with extremely low log2 coverage values.

        These are generally bins that had no reads mapped, and so were
        substituted with a small dummy log2 value to avoid divide-by-zero
        errors.
        """
        return self.as_dataframe(self.data[self.data['log2'] >
                                           params.NULL_LOG2_COVERAGE -
                                           params.MIN_REF_COVERAGE])

    def squash_genes(self, summary_func=descriptives.biweight_location,
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
            start = rows.start.iat[0]
            end = rows.end.iat[-1]
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

    def shift_xx(self, male_reference=False, is_xx=None):
        """Adjust chrX coverages (divide in half) for apparent female samples."""
        outprobes = self.copy()
        if is_xx is None:
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
        if is_male_reference:
            # Single-copy X, Y
            idx = np.asarray((self.chromosome == self._chr_x_label) |
                             (self.chromosome == self._chr_y_label))
        else:
            # Y will be all noise, so replace with 1 "flat" copy
            idx = np.asarray(self.chromosome == self._chr_y_label)
        # ENH: do math with ratio, then calculate log2 from it
        ratio = np.ones(len(self), dtype=np.float_)
        log2 = np.zeros(len(self), dtype=np.float_)
        ratio[idx] = 0.5
        log2[idx] = -1.0
        depth = ratio * self.guess_average_depth()
        return pd.DataFrame({'depth': depth,
                             'log2': log2,
                             'ratio': ratio})

    # Reporting

    def residuals(self, segments=None):
        """Difference in log2 value of each bin from its segment mean.

        If segments are just regions (e.g. GenomicArray) with no log2 values
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
        if segments is None and window:
            y_log2 -= smoothing.smoothed(y_log2, window)
        # Guess Poisson parameter from absolute-scale values
        y = np.exp2(y_log2)
        # ENH: use weight argument to these stats
        loc = descriptives.biweight_location(y)
        spread = descriptives.biweight_midvariance(y, loc)
        scale = loc / spread**2
        return scale
