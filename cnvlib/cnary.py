"""CNVkit's core data structure, a copy number array."""
from __future__ import print_function, absolute_import, division
from builtins import map
from past.builtins import basestring

import logging

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata

from . import core, gary, descriptives, params, smoothing
from .metrics import segment_mean


class CopyNumArray(gary.GenomicArray):
    """An array of genomic intervals, treated like aCGH probes.

    Required columns: chromosome, start, end, gene, log2

    Optional columns: gc, rmask, spread, weight, probes
    """
    _required_columns = ("chromosome", "start", "end", "gene", "log2")
    _required_dtypes = (str, int, int, str, float)
    # ENH: make gene optional
    # Extra columns, in order:
    # "depth", "gc", "rmask", "spread", "weight", "probes"

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

        Group each series of intergenic bins as a "Background" gene; any
        "Background" bins within a gene are grouped with that gene.

        Bins' gene names are split on commas to accommodate overlapping genes
        and bins that cover multiple genes.

        Parameters
        ----------
        ignore : list or tuple of str
            Gene names to treat as "Background" bins instead of real genes,
            grouping these bins with the surrounding gene or background region.
            These bins will still retain their name in the output.

        Yields
        ------
        tuple
            Pairs of: (gene name, CNA of rows with same name)
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

    def center_all(self, estimator=pd.Series.median, skip_low=False, by_chrom=True):
        """Re-center log2 values to the autosomes' average (in-place).

        Parameters
        ----------
        estimator : str or callable
            Function to estimate central tendency. If a string, must be one of
            'mean', 'median', 'mode', 'biweight' (for biweight location). Median
            by default.
        skip_low : bool
            Whether to drop very-low-coverage bins (via `drop_low_coverage`)
            before estimating the center value.
        by_chrom : bool
            If True, first apply `estimator` to each chromosome separately, then
            apply `estimator` to the per-chromosome values, to reduce the impact
            of uneven targeting or extreme aneuploidy. Otherwise, apply
            `estimator` to all log2 values directly.
        """
        est_funcs = {
            "mean": pd.Series.mean,
            "median": pd.Series.median,
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
                values = pd.Series([estimator(subarr['log2'])
                                    for _c, subarr in cnarr.by_chromosome()
                                    if len(subarr)])
            else:
                values = cnarr['log2']
            self.data['log2'] -= estimator(values)

    def drop_low_coverage(self):
        """Drop bins with extremely low log2 coverage or copy ratio values.

        These are generally bins that had no reads mapped due to sample-specific
        issues. A very small log2 ratio or coverage value may have been
        substituted to avoid domain or divide-by-zero errors.
        """
        return self.as_dataframe(self.data[self.data['log2'] >
                                           params.NULL_LOG2_COVERAGE -
                                           params.MIN_REF_COVERAGE])

    def squash_genes(self, summary_func=descriptives.biweight_location,
                     squash_background=False, ignore=params.IGNORE_GENE_NAMES):
        """Combine consecutive bins with the same targeted gene name.

        Parameters
        ----------
        summary_func : callable
            Function to summarize an array of log2 values to produce a
            new log2 value for a "squashed" (i.e. reduced) region. By default
            this is the biweight location, but you might want median, mean, max,
            min or something else in some cases.
        squash_background : bool
            If True, also reduce consecutive "Background" bins into a single
            bin. Otherwise, keep "Background" and ignored bins as they are in
            the output.
        ignore : list or tuple of str
            Bin names to be treated as "Background" instead of as unique genes.

        Return
        ------
        CopyNumArray
            Another, usually smaller, copy of `self` with each gene's bins
            reduced to a single bin with appropriate values.
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

    # Chromosomal sex (gender)

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
        """Detect chromosomal sex; return True if a sample is probably female.

        Uses `compare_sex_chromosomes` to calculate coverage ratios of the X and
        Y chromosomes versus autosomes.

        Parameters
        ----------
        male_reference : bool
            Was this sample normalized to a male reference copy number profile?
        verbose : bool
            If True, print (i.e. log to console) the ratios of the log2
            coverages of the X and Y chromosomes versus autosomes, the
            "maleness" ratio of male vs. female expectations for each sex
            chromosome, and the inferred chromosomal sex.

        Returns
        -------
        bool
            True if the coverage ratios indicate the sample is female.
        """
        is_xy, stats = self.compare_sex_chromosomes(male_reference)
        if is_xy is None:
            return None
        elif verbose:
            logging.info("Relative log2 coverage of %s=%.3g, %s=%.3g "
                         "(maleness=%.3g x %.3g = %.3g) --> assuming %s",
                         self._chr_x_label, stats['chrx_ratio'],
                         self._chr_y_label, stats['chry_ratio'],
                         stats['chrx_male_lr'], stats['chry_male_lr'],
                         stats['chrx_male_lr'] * stats['chry_male_lr'],
                         ('female', 'male')[is_xy])
        return ~is_xy

    def compare_sex_chromosomes(self, male_reference=False, skip_low=False):
        """Compare coverage ratios of sex chromosomes versus autosomes.

        Perform 4 Mann-Whitney tests of the log2 coverages on chromosomes X and
        Y, separately shifting for assumed male and female chromosomal sex.
        Compare the U values obtained to infer whether the male or female
        assumption fits the data better.

        Parameters
        ----------
        male_reference : bool
            Whether a male reference copy number profile was used to normalize
            the data. If so, a male sample should have log2 values of 0 on X and
            Y, and female +1 on X, deep negative (below -3) on Y. Otherwise, a
            male sample should have log2 values of -1 on X and 0 on
            Y, and female 0 on X, deep negative (below -3) on Y.
        skip_low : bool
            If True, drop very-low-coverage bins (via `drop_low_coverage`)
            before comparing log2 coverage ratios. Included for completeness,
            but shouldn't affect the result much since the M-W test is
            nonparametric and p-values are not used here.

        Returns
        -------
        bool
            True if the sample appears male.
        dict
            Calculated values used for the inference: relative log2 ratios of
            chromosomes X and Y versus the autosomes; the Mann-Whitney U values
            from each test; and ratios of U values for male vs. female
            assumption on chromosomes X and Y.
        """
        cnarr = self.drop_low_coverage() if skip_low else self
        chrx = cnarr[cnarr.chromosome == self._chr_x_label]
        if not len(chrx):
            logging.warn("*WARNING* No %s found in probes; check the input",
                         self._chr_x_label)
            return None, {}

        chrx_l = chrx['log2'].values
        auto = cnarr.autosomes()
        auto_l = auto['log2'].values

        # ENH: use scipy.stats.ttest_ind(a, b, equal_var=False)
        def compare_to_auto(vals):
            # u = mannwhitneyu(auto_l, vals, alternative='two-sided')[0]
            # From scipy.stats.mannwhitneyu -- skip the p-value calculation
            # (the tiebreaker there can raise an irrelevant error)
            n1 = len(auto_l)
            n2 = len(vals)
            ranked = rankdata(np.concatenate((auto_l, vals)))
            rankx = ranked[0:n1]  # get the x-ranks
            u1 = n1*n2 + (n1*(n1+1))/2.0 - np.sum(rankx, axis=0)  # calc U for x
            return n1*n2 - u1 # remainder is U for y

        if male_reference:
            female_chrx_u = compare_to_auto(chrx_l - 1)
            male_chrx_u = compare_to_auto(chrx_l)
        else:
            female_chrx_u = compare_to_auto(chrx_l)
            male_chrx_u = compare_to_auto(chrx_l + 1)
        # Mann-Whitney U score is greater for similar-mean sets
        chrx_male_lr = male_chrx_u / female_chrx_u
        # Similar for chrY if it's present
        chry = cnarr[cnarr.chromosome == self._chr_y_label]
        if len(chry):
            chry_l = chry['log2']
            male_chry_u = compare_to_auto(chry_l)
            female_chry_u = compare_to_auto(chry_l + 3)
            chry_male_lr = male_chry_u / female_chry_u
        else:
            # If chrY is missing, don't sabotage the inference
            male_chry_u = female_chry_u = np.nan
            chry_male_lr = 1.0
        # Relative log2 values, for convenient reporting
        auto_mean = segment_mean(auto, skip_low=skip_low)
        chrx_mean = segment_mean(chrx, skip_low=skip_low)
        chry_mean = segment_mean(chry, skip_low=skip_low)
        return (chrx_male_lr * chry_male_lr > 1.0,
                dict(chrx_ratio=chrx_mean - auto_mean,
                     chry_ratio=chry_mean - auto_mean,
                     # For debugging, mainly
                     chrx_male_lr=chrx_male_lr,
                     chry_male_lr=chry_male_lr,
                     female_chrx_u=female_chrx_u,
                     male_chrx_u=male_chrx_u,
                     female_chry_u=female_chry_u,
                     male_chry_u=male_chry_u,
                    ))

    def expect_flat_log2(self, is_male_reference=None):
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

        Parameters
        ----------
        segments : GenomicArray, CopyNumArray, or None
            Determines the "mean" value to which `self` log2 values are relative:

            - If CopyNumArray, use the log2 values as the segment means to
              subtract.
            - If GenomicArray with no log2 values, group `self` by these ranges
              and subtract each group's median log2 value.
            - If None, subtract each chromosome's median.

        Returns
        -------
        array
            Residual log2 values from `self` relative to `segments`; same length
            as `self`.
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

    def _guess_average_depth(self, segments=None, window=100):
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
        if spread > 0:
            return loc / spread**2
        return loc
