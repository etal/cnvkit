"""CNVkit's core data structure, a copy number array."""

from __future__ import annotations
import logging
from typing import TYPE_CHECKING, Any, Optional, Union
from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy.stats import median_test

from skgenome import GenomicArray
from . import core, descriptives, params, smoothing
from .segmetrics import segment_mean

if TYPE_CHECKING:
    from collections.abc import Iterator
    from numpy import bool_, float64, ndarray
    from pandas.core.frame import DataFrame
    from pandas.core.series import Series


class CopyNumArray(GenomicArray):
    """An array of genomic intervals, treated like aCGH probes.

    Required columns: chromosome, start, end, gene, log2

    Optional columns: gc, rmask, spread, weight, probes
    """

    _required_columns = ("chromosome", "start", "end", "gene", "log2")  # type: ignore
    _required_dtypes = (str, int, int, str, float)  # type: ignore
    # ENH: make gene optional
    # Extra columns, in order:
    # "depth", "gc", "rmask", "spread", "weight", "probes"

    def __init__(
        self, data_table: DataFrame, meta_dict: Optional[dict[str, str]] = None
    ) -> None:
        GenomicArray.__init__(self, data_table, meta_dict)

    @property
    def log2(self):
        return self.data["log2"]

    @log2.setter
    def log2(self, value) -> None:
        self.data["log2"] = value

    @property
    def chr_x_label(self) -> str:
        """The name of the X chromosome.

        This is either "X" or "chrX".
        """
        key = "chr_x"
        if key in self.meta:
            return self.meta[key]
        if len(self):
            chr_x_label = "chrX" if self.chromosome.iat[0].startswith("chr") else "X"
            self.meta[key] = chr_x_label
            return chr_x_label
        return ""

    def chr_x_filter(self, diploid_parx_genome: Optional[str] = None) -> Series:
        """All regions on X, potentially without PAR1/2."""
        x = self.chromosome == self.chr_x_label
        if diploid_parx_genome is not None:
            # Exclude PAR since they are expected to be diploid (i.e. autosomal).
            x &= ~self.parx_filter(genome_build=diploid_parx_genome)
        return x

    def parx_filter(self, genome_build: str) -> Series:
        """All PAR1/2 regions on X."""
        genome_build = genome_build.lower()
        assert genome_build in params.SUPPORTED_GENOMES_FOR_PAR_HANDLING
        f = self.chromosome == self.chr_x_label
        par1_start, par1_end = params.PSEUDO_AUTSOMAL_REGIONS[genome_build]["PAR1X"]
        par2_start, par2_end = params.PSEUDO_AUTSOMAL_REGIONS[genome_build]["PAR2X"]
        f &= ((self.start >= par1_start) & (self.end <= par1_end)) | (
            (self.start >= par2_start) & (self.end <= par2_end)
        )
        return f

    @property
    def chr_y_label(self) -> str:
        """The name of the Y chromosome."""
        if "chr_y" in self.meta:
            return self.meta["chr_y"]
        if len(self):
            chr_y = "chrY" if self.chr_x_label.startswith("chr") else "Y"
            self.meta["chr_y"] = chr_y
            return chr_y
        return ""

    def pary_filter(self, genome_build: str) -> Series:
        """All PAR1/2 regions on Y."""
        genome_build = genome_build.lower()
        assert genome_build in params.SUPPORTED_GENOMES_FOR_PAR_HANDLING
        f = self.chromosome == self.chr_y_label
        par1_start, par1_end = params.PSEUDO_AUTSOMAL_REGIONS[genome_build]["PAR1Y"]
        par2_start, par2_end = params.PSEUDO_AUTSOMAL_REGIONS[genome_build]["PAR2Y"]
        f &= ((self.start >= par1_start) & (self.end <= par1_end)) | (
            (self.start >= par2_start) & (self.end <= par2_end)
        )
        return f

    def chr_y_filter(self, diploid_parx_genome: Optional[str] = None) -> Series:
        """All regions on Y, potentially without PAR1/2."""
        y = self.chromosome == self.chr_y_label
        if diploid_parx_genome is not None:
            # Exclude PAR on Y since they cannot be covered (everything is mapped to X).
            y &= ~self.pary_filter(genome_build=diploid_parx_genome)
        return y

    def autosomes(
        self, diploid_parx_genome: Optional[str] = None, also: Optional[Series] = None
    ) -> CopyNumArray:
        """Overrides GenomeArray.autosomes()."""
        if diploid_parx_genome is not None:
            if also is None:
                also = self.parx_filter(diploid_parx_genome)
            elif isinstance(also, pd.Series):
                also |= self.parx_filter(diploid_parx_genome)
            else:
                raise NotImplementedError("Cannot combine pd.Series with non-Series.")
        return super().autosomes(also=also)

    # More meta to store:
    #   is_sample_male = bool
    #   is_haploid_x_reference = bool
    #   * invalidate 'chr_x' if .chromosome/['chromosome'] is set

    # Traversal

    def by_gene(
        self, ignore: tuple[str, str, str] = params.IGNORE_GENE_NAMES
    ) -> Iterator[tuple[str, CopyNumArray]]:
        """Iterate over probes grouped by gene name.

        Group each series of intergenic bins as an "Antitarget" gene; any
        "Antitarget" bins within a gene are grouped with that gene.

        Bins' gene names are split on commas to accommodate overlapping genes
        and bins that cover multiple genes.

        Parameters
        ----------
        ignore : list or tuple of str
            Gene names to treat as "Antitarget" bins instead of real genes,
            grouping these bins with the surrounding gene or intergenic region.
            These bins will still retain their name in the output.

        Yields
        ------
        tuple
            Pairs of: (gene name, CNA of rows with same name)
        """
        ignore += params.ANTITARGET_ALIASES
        for _chrom, subgary in self.by_chromosome():
            prev_idx = 0
            for gene, gene_idx in subgary._get_gene_map().items():
                if gene not in ignore:
                    if not len(gene_idx):
                        logging.warning("Specified gene name somehow missing: %s", gene)
                        continue
                    start_idx = gene_idx[0]
                    end_idx = gene_idx[-1] + 1
                    if prev_idx < start_idx:
                        # Include intergenic regions
                        yield (
                            params.ANTITARGET_NAME,
                            subgary.as_dataframe(subgary.data.loc[prev_idx:start_idx]),
                        )
                    yield (
                        gene,
                        subgary.as_dataframe(subgary.data.loc[start_idx:end_idx]),
                    )
                    prev_idx = end_idx
            if prev_idx < len(subgary) - 1:
                # Include the telomere
                yield (
                    params.ANTITARGET_NAME,
                    subgary.as_dataframe(subgary.data.loc[prev_idx:]),
                )

    # Manipulation

    def center_all(
        self,
        estimator: Union[Callable, str] = pd.Series.median,
        by_chrom: bool = True,
        skip_low: bool = False,
        verbose: bool = False,
        diploid_parx_genome: Optional[str] = None,
    ) -> None:
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
        diploid_parx_genome : String
             Whether to include the PAR1/2 on chr X from the given genome (build)
             as part of the autosomes
        """
        est_funcs = {
            "mean": pd.Series.mean,
            "median": pd.Series.median,
            "mode": descriptives.modal_location,
            "biweight": descriptives.biweight_location,
        }
        if isinstance(estimator, str):
            if estimator in est_funcs:
                estimator = est_funcs[estimator]
            else:
                raise ValueError(
                    "Estimator must be a function or one of: "
                    + ", ".join(map(repr, est_funcs))
                )
        cnarr = (
            self.drop_low_coverage(verbose=verbose) if skip_low else self
        ).autosomes(diploid_parx_genome=diploid_parx_genome)
        if cnarr:
            if by_chrom:
                values = pd.Series(
                    [
                        estimator(subarr["log2"])
                        for _c, subarr in cnarr.by_chromosome()
                        if len(subarr)
                    ]
                )
            else:
                values = cnarr["log2"]
            shift = -estimator(values)
            if verbose:
                logging.info("Shifting log2 values by %f", shift)
            self.data["log2"] += shift

    def drop_low_coverage(self, verbose: bool = False) -> CopyNumArray:
        """Drop bins with extremely low log2 coverage or copy ratio values.

        These are generally bins that had no reads mapped due to sample-specific
        issues. A very small log2 ratio or coverage value may have been
        substituted to avoid domain or divide-by-zero errors.
        """
        min_cvg = params.NULL_LOG2_COVERAGE - params.MIN_REF_COVERAGE
        drop_idx = self.data["log2"] < min_cvg
        if "depth" in self:
            drop_idx |= self.data["depth"] == 0
        if verbose and drop_idx.any():
            logging.info("Dropped %d low-coverage bins", drop_idx.sum())
        return self[~drop_idx]

    def squash_genes(
        self,
        summary_func: Callable = descriptives.biweight_location,
        squash_antitarget: bool = False,
        ignore: tuple[str, str, str] = params.IGNORE_GENE_NAMES,
    ) -> CopyNumArray:
        """Combine consecutive bins with the same targeted gene name.

        Parameters
        ----------
        summary_func : callable
            Function to summarize an array of log2 values to produce a
            new log2 value for a "squashed" (i.e. reduced) region. By default
            this is the biweight location, but you might want median, mean, max,
            min or something else in some cases.
        squash_antitarget : bool
            If True, also reduce consecutive "Antitarget" bins into a single
            bin. Otherwise, keep "Antitarget" and ignored bins as they are in
            the output.
        ignore : list or tuple of str
            Bin names to be treated as "Antitarget" instead of as unique genes.

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
            # Build row matching self.data.columns order exactly
            outrow = []
            for col in self.data.columns:
                if col == "chromosome":
                    outrow.append(core.check_unique(rows.chromosome, "chromosome"))
                elif col == "start":
                    outrow.append(rows.start.iat[0])
                elif col == "end":
                    outrow.append(rows.end.iat[-1])
                elif col == "gene":
                    outrow.append(name)
                elif col == "log2":
                    outrow.append(summary_func(rows.log2))
                elif col == "probes":
                    # Special case: sum probes rather than average
                    outrow.append(sum(rows[col]))
                else:
                    # All other fields: use summary function
                    outrow.append(summary_func(rows[col]))
            return tuple(outrow)

        outrows = []
        for name, subarr in self.by_gene(ignore):
            if not len(subarr):
                continue
            if name in params.ANTITARGET_ALIASES and not squash_antitarget:
                outrows.extend(subarr.data.itertuples(index=False))
            else:
                outrows.append(squash_rows(name, subarr.data))
        return self.as_rows(outrows)

    # Chromosomal sex

    def shift_xx(
        self,
        is_haploid_x_reference: bool = False,
        is_xx: Optional[bool_] = None,
        diploid_parx_genome: Optional[str] = None,
    ) -> CopyNumArray:
        """Adjust chrX log2 ratios to match the ploidy of the reference sex.

        I.e. add 1 to chrX log2 ratios for a male sample vs. female reference,
        or subtract 1 for a female sample vs. male reference, so that chrX log2
        values are comparable across samples with different chromosomal sex.
        """
        outprobes = self.copy()
        if is_xx is None:
            is_xx = self.guess_xx(
                is_haploid_x_reference=is_haploid_x_reference,
                diploid_parx_genome=diploid_parx_genome,
            )
        if is_xx and is_haploid_x_reference:
            # Female: divide X coverages by 2 (in log2: subtract 1)
            outprobes[outprobes.chromosome == self.chr_x_label, "log2"] -= 1.0
            # Male: no change
        elif not is_xx and not is_haploid_x_reference:
            # Male: multiply X coverages by 2 (in log2: add 1)
            outprobes[outprobes.chromosome == self.chr_x_label, "log2"] += 1.0
            # Female: no change
        return outprobes

    def guess_xx(
        self,
        is_haploid_x_reference: bool = False,
        diploid_parx_genome: Optional[str] = None,
        verbose: bool = True,
    ) -> Optional[bool_]:
        """Detect chromosomal sex; return True if a sample is probably female.

        Uses `compare_sex_chromosomes` to calculate coverage ratios of the X and
        Y chromosomes versus autosomes.

        Parameters
        ----------
        is_haploid_x_reference : bool
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
        is_xy, stats = self.compare_sex_chromosomes(
            is_haploid_x_reference, diploid_parx_genome
        )
        if is_xy is None:
            return None
        if verbose:
            logging.info(
                "Relative log2 coverage of %s=%.3g, %s=%.3g "
                "(maleness=%.3g x %.3g = %.3g) --> assuming %s",
                self.chr_x_label,
                stats["chrx_ratio"],
                self.chr_y_label,
                stats["chry_ratio"],
                stats["chrx_male_lr"],
                stats["chry_male_lr"],
                stats["chrx_male_lr"] * stats["chry_male_lr"],
                "male" if is_xy else "female",
            )
        return ~is_xy

    def compare_sex_chromosomes(
        self,
        is_haploid_x_reference: bool = False,
        diploid_parx_genome: Optional[str] = None,
        skip_low: bool = False,
    ) -> Union[
        tuple[bool_, dict[str, Union[float64, float]]],
        tuple[bool_, dict[str, float64]],
        tuple[None, dict[Any, Any]],
    ]:
        """Compare coverage ratios of sex chromosomes versus autosomes.

        Perform 4 Mood's median tests of the log2 coverages on chromosomes X and
        Y, separately shifting for assumed male and female chromosomal sex.
        Compare the chi-squared values obtained to infer whether the male or
        female assumption fits the data better.

        Parameters
        ----------
        is_haploid_x_reference : bool
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
        if not len(self):
            return None, {}

        chrx = self[self.chr_x_filter(diploid_parx_genome)]
        if not len(chrx):
            logging.warning(
                "No %s found in sample; is the input truncated?", self.chr_x_label
            )
            return None, {}

        auto = self.autosomes(diploid_parx_genome=diploid_parx_genome)
        if skip_low:
            chrx = chrx.drop_low_coverage()
            auto = auto.drop_low_coverage()
        auto_l = auto["log2"].to_numpy()
        use_weight = "weight" in self
        auto_w = auto["weight"].to_numpy() if use_weight else None

        def compare_to_auto(vals, weights):
            # Mood's median test stat is chisq -- near 0 for similar median
            try:
                stat, _p, _med, cont = median_test(
                    auto_l, vals, ties="ignore", lambda_="log-likelihood"
                )
            except ValueError:
                # "All values are below the grand median (0.0)"
                stat = None
            else:
                if stat == 0 and 0 in cont:
                    stat = None
            # In case Mood's test failed for either sex
            if use_weight:
                med_diff = abs(
                    descriptives.weighted_median(auto_l, auto_w)
                    - descriptives.weighted_median(vals, weights)
                )
            else:
                med_diff = abs(np.median(auto_l) - np.median(vals))
            return (stat, med_diff)

        def compare_chrom(vals, weights, female_shift, male_shift):
            """Calculate "maleness" ratio of test statistics.

            The ratio is of the female vs. male chi-square test statistics from
            the median test. If the median test fails for either sex, (due to
            flat/trivial input), use the ratio of the absolute difference in
            medians.
            """
            female_stat, f_diff = compare_to_auto(vals + female_shift, weights)
            male_stat, m_diff = compare_to_auto(vals + male_shift, weights)
            # Statistic is smaller for similar-median sets
            if female_stat is not None and male_stat is not None:
                return female_stat / max(male_stat, 0.01)
            # Difference in medians is also smaller for similar-median sets
            return f_diff / max(m_diff, 0.01)

        female_x_shift, male_x_shift = (-1, 0) if is_haploid_x_reference else (0, +1)
        chrx_male_lr = compare_chrom(
            chrx["log2"].to_numpy(),
            (chrx["weight"].to_numpy() if use_weight else None),
            female_x_shift,
            male_x_shift,
        )
        combined_score = chrx_male_lr
        # Similar for chrY if it's present
        chry = self[self.chr_y_filter(diploid_parx_genome)]
        if len(chry):
            if skip_low:
                chry = chry.drop_low_coverage()
            chry_male_lr = compare_chrom(
                chry["log2"].to_numpy(),
                (chry["weight"].to_numpy() if use_weight else None),
                +3,
                0,
            )
            if np.isfinite(chry_male_lr):
                combined_score *= chry_male_lr
        else:
            # If chrY is missing, don't sabotage the inference
            chry_male_lr = np.nan
        # Relative log2 values, for convenient reporting
        auto_mean = segment_mean(auto, skip_low=skip_low)
        chrx_mean = segment_mean(chrx, skip_low=skip_low)
        chry_mean = segment_mean(chry, skip_low=skip_low)
        return (
            combined_score > 1.0,
            dict(
                chrx_ratio=chrx_mean - auto_mean,
                chry_ratio=chry_mean - auto_mean,
                combined_score=combined_score,
                # For debugging, mainly
                chrx_male_lr=chrx_male_lr,
                chry_male_lr=chry_male_lr,
            ),
        )

    def expect_flat_log2(
        self,
        is_haploid_x_reference: Optional[bool] = None,
        diploid_parx_genome: Optional[str] = None,
    ) -> ndarray:
        """Get the uninformed expected copy ratios of each bin.

        Create an array of log2 coverages like a "flat" reference.

        This is a neutral copy ratio at each autosome (log2 = 0.0) and sex
        chromosomes based on whether the reference is male (XX or XY).
        """
        if is_haploid_x_reference is None:
            is_haploid_x_reference = not self.guess_xx(
                diploid_parx_genome=diploid_parx_genome, verbose=False
            )
        cvg = np.zeros(len(self), dtype=np.float64)
        if is_haploid_x_reference:
            # Single-copy X, Y
            idx = (
                self.chr_x_filter(diploid_parx_genome).to_numpy()
                | (self.chr_y_filter(diploid_parx_genome)).to_numpy()
            )
        else:
            # Y will be all noise, so replace with 1 "flat" copy, including PAR1/2.
            idx = (self.chr_y_filter()).to_numpy()
        cvg[idx] = -1.0
        return cvg

    # Reporting

    def residuals(
        self, segments: Optional[Union[CopyNumArray, GenomicArray]] = None
    ) -> Series:
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
            resids = [
                subcna.log2 - subcna.log2.median()
                for _chrom, subcna in self.by_chromosome()
            ]
        elif "log2" in segments:
            resids = [
                bins_lr - seg_lr
                for seg_lr, bins_lr in zip(
                    segments["log2"],
                    self.iter_ranges_of(
                        segments, "log2", mode="inner", keep_empty=True
                    ), strict=False,
                )
                if len(bins_lr)
            ]
        else:
            resids = [
                lr - lr.median()
                for lr in self.iter_ranges_of(segments, "log2", keep_empty=False)
            ]
        return pd.concat(resids) if resids else pd.Series([])

    def smooth_log2(self, bandwidth: None = None, by_arm: bool = True) -> ndarray:
        """Smooth log2 values with a sliding window.

        Account for chromosome and (optionally) centromere boundaries. Use bin
        weights if present.

        Returns
        -------
        array
            Smoothed log2 values from `self`, the same length as `self`.
        """
        if bandwidth is None:
            bandwidth = smoothing.guess_window_size(
                self.log2, weights=(self["weight"] if "weight" in self else None)
            )

        parts = self.by_arm() if by_arm else self.by_chromosome()
        if "weight" in self:
            out = [
                smoothing.savgol(
                    subcna["log2"].to_numpy(),
                    bandwidth,
                    weights=subcna["weight"].to_numpy(),
                )
                for _chrom, subcna in parts
            ]
        else:
            out = [
                smoothing.savgol(subcna["log2"].to_numpy(), bandwidth)
                for _chrom, subcna in parts
            ]
        return np.concatenate(out)

    def _guess_average_depth(self, segments=None, window=100, diploid_parx_genome=None):
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
        cnarr = self.autosomes(diploid_parx_genome=diploid_parx_genome)
        if not len(cnarr):
            cnarr = self
        # Remove variations due to real/likely CNVs
        y_log2 = cnarr.residuals(segments)
        if segments is None and window:
            y_log2 -= smoothing.savgol(y_log2, window)
        # Guess Poisson parameter from absolute-scale values
        y = np.exp2(y_log2)
        # ENH: use weight argument to these stats
        loc = descriptives.biweight_location(y)
        spread = descriptives.biweight_midvariance(y, loc)
        if spread > 0:
            return loc / spread**2
        return loc
