"""CNVkit's core data structure, a copy number array."""

from __future__ import annotations
import logging
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from skgenome import GenomicArray
from skgenome.chromnames import infer_sex_chrom_labels
from skgenome.genomebuild import get_genome_build
from . import core, descriptives, params, smoothing

if TYPE_CHECKING:
    from collections.abc import Callable
    from collections.abc import Iterator
    from numpy import bool_, float64, ndarray
    from pandas.core.frame import DataFrame
    from pandas.core.series import Series


def is_female_default(is_xx: bool | bool_ | None) -> bool:
    """Resolve an indeterminate chromosomal-sex guess to the safe default.

    `CopyNumArray.guess_xx` returns None when sex can't be determined (no X
    chromosome present, or a degenerate coverage comparison). For copy-number
    purposes the only effect of a "male" call is that a haploid X is treated as
    the expected baseline; assuming male for a sample whose X is actually
    diploid would spuriously inflate chrX by +1 (the gh#360 failure family). So
    when there isn't positive evidence either way, default to female (diploid-X
    baseline). A real (non-None) guess is returned unchanged as a plain ``bool``.
    """
    return True if is_xx is None else bool(is_xx)


class CopyNumArray(GenomicArray):
    """An array of genomic intervals, treated like aCGH probes.

    Required columns: chromosome, start, end, gene, log2

    Optional columns: gc, rmask, spread, weight, probes
    """

    _required_columns = ("chromosome", "start", "end", "gene", "log2")  # type: ignore[assignment]
    _required_dtypes = (str, int, int, str, float)  # type: ignore[assignment]
    # ENH: make gene optional
    # Extra columns, in order:
    # "depth", "gc", "rmask", "spread", "weight", "probes"

    def __init__(
        self, data_table: DataFrame, meta_dict: dict[str, str] | None = None
    ) -> None:
        GenomicArray.__init__(self, data_table, meta_dict)

    @property
    def log2(self):
        return self.data["log2"]

    @log2.setter
    def log2(self, value) -> None:
        self.data["log2"] = value

    def _ensure_sex_chrom_labels(self) -> None:
        """Populate ``meta['chr_x']`` and ``meta['chr_y']`` if not already set."""
        if "chr_x" in self.meta and "chr_y" in self.meta:
            return
        if len(self):
            x, y = infer_sex_chrom_labels(self.chromosome.unique())
        else:
            x, y = None, None
        self.meta.setdefault("chr_x", x)
        self.meta.setdefault("chr_y", y)

    @property
    def chr_x_label(self) -> str | None:
        """The name of the X chromosome in this dataset, if present.

        Returns ``"chrX"`` or ``"X"`` if one of those is present in the data,
        or ``None`` if no X chromosome is detected. For yeast and other
        Roman-numeral genomes where ``chrX`` is autosome 10, returns ``None``.
        """
        self._ensure_sex_chrom_labels()
        label = self.meta["chr_x"]
        return label if label is None else str(label)

    def chr_x_filter(self, diploid_parx_genome: str | None = None) -> Series:
        """All regions on X, potentially without PAR1/2.

        Returns an all-False Series when no X chromosome is detected.
        """
        label = self.chr_x_label
        if label is None:
            return pd.Series(False, index=self.data.index)
        x = self.chromosome == label
        if diploid_parx_genome is not None:
            # Exclude PAR since they are expected to be diploid (i.e. autosomal).
            x &= ~self.parx_filter(genome_build=diploid_parx_genome)
        return x

    def parx_filter(self, genome_build: str) -> Series:
        """All PAR1/2 regions on X."""
        build = get_genome_build(genome_build)
        label = self.chr_x_label
        if label is None:
            return pd.Series(False, index=self.data.index)
        f = self.chromosome == label
        par1_start, par1_end = build.par_regions["PAR1X"]
        par2_start, par2_end = build.par_regions["PAR2X"]
        f &= ((self.start >= par1_start) & (self.end <= par1_end)) | (
            (self.start >= par2_start) & (self.end <= par2_end)
        )
        return f

    @property
    def chr_y_label(self) -> str | None:
        """The name of the Y chromosome in this dataset, if present.

        Returns ``None`` when no Y chromosome is detected.
        """
        self._ensure_sex_chrom_labels()
        label = self.meta["chr_y"]
        return label if label is None else str(label)

    def pary_filter(self, genome_build: str) -> Series:
        """All PAR1/2 regions on Y."""
        build = get_genome_build(genome_build)
        label = self.chr_y_label
        if label is None:
            return pd.Series(False, index=self.data.index)
        f = self.chromosome == label
        par1_start, par1_end = build.par_regions["PAR1Y"]
        par2_start, par2_end = build.par_regions["PAR2Y"]
        f &= ((self.start >= par1_start) & (self.end <= par1_end)) | (
            (self.start >= par2_start) & (self.end <= par2_end)
        )
        return f

    def chr_y_filter(self, diploid_parx_genome: str | None = None) -> Series:
        """All regions on Y, potentially without PAR1/2.

        Returns an all-False Series when no Y chromosome is detected.
        """
        label = self.chr_y_label
        if label is None:
            return pd.Series(False, index=self.data.index)
        y = self.chromosome == label
        if diploid_parx_genome is not None:
            # Exclude PAR on Y since they cannot be covered (everything is mapped to X).
            y &= ~self.pary_filter(genome_build=diploid_parx_genome)
        return y

    def autosomes(  # type: ignore[override]
        self, diploid_parx_genome: str | None = None, also: Series | None = None
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
        self, ignore: tuple[str, ...] = params.IGNORE_GENE_NAMES
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
            prev_pos = 0
            for gene, gene_idx_labels in subgary._get_gene_map().items():
                if gene not in ignore:
                    if not len(gene_idx_labels):
                        logging.warning("Specified gene name somehow missing: %s", gene)
                        continue
                    # Convert index labels to positional indices
                    start_label = gene_idx_labels[0]
                    end_label = gene_idx_labels[-1]
                    start_loc = subgary.data.index.get_loc(start_label)
                    end_loc = subgary.data.index.get_loc(end_label)

                    # get_loc() can return int, slice, or boolean array depending on index type:
                    # - int: unique label -> position
                    # - slice: duplicate consecutive labels -> range of positions
                    # - ndarray (bool): duplicate non-consecutive labels -> mask of positions
                    # For duplicates, use the first occurrence for start, last for end
                    if isinstance(start_loc, slice):
                        start_pos = start_loc.start
                    elif isinstance(start_loc, np.ndarray):
                        # Boolean mask: convert to positions and take first
                        start_pos = np.where(start_loc)[0][0]
                    else:
                        start_pos = start_loc

                    if isinstance(end_loc, slice):
                        end_pos = end_loc.stop
                    elif isinstance(end_loc, np.ndarray):
                        # Boolean mask: convert to positions, take last, +1 for exclusive end
                        end_pos = np.where(end_loc)[0][-1] + 1
                    else:
                        end_pos = end_loc + 1

                    if prev_pos < start_pos:
                        # Include intergenic regions
                        yield (
                            params.ANTITARGET_NAME,
                            subgary.as_dataframe(subgary.data.iloc[prev_pos:start_pos]),
                        )
                    yield (
                        gene,
                        subgary.as_dataframe(subgary.data.iloc[start_pos:end_pos]),
                    )
                    prev_pos = end_pos
            if prev_pos < len(subgary):
                # Include the telomere
                yield (
                    params.ANTITARGET_NAME,
                    subgary.as_dataframe(subgary.data.iloc[prev_pos:]),
                )

    # Manipulation

    def center_all(
        self,
        estimator: Callable | str = pd.Series.median,
        by_chrom: bool = True,
        skip_low: bool = False,
        verbose: bool = False,
        diploid_parx_genome: str | None = None,
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
        assert callable(estimator)
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

        NaN log2/depth values (from degenerate or malformed inputs) are also
        dropped here: a bare ``log2 < min_cvg`` comparison is False for NaN, so
        without this NaN bins would slip past into segmentation and trigger an
        'invalid value encountered in greater' warning or a CBS crash (#521).
        """
        min_cvg = params.NULL_LOG2_COVERAGE - params.MIN_REF_COVERAGE
        log2 = self.data["log2"]
        drop_idx = (log2 < min_cvg) | log2.isna()
        if "depth" in self:
            depth = self.data["depth"]
            drop_idx |= (depth == 0) | depth.isna()
        if verbose and drop_idx.any():
            logging.info("Dropped %d low-coverage bins", drop_idx.sum())
        return self[~drop_idx]

    def squash_genes(
        self,
        summary_func: Callable = descriptives.biweight_location,
        squash_antitarget: bool = False,
        ignore: tuple[str, ...] = params.IGNORE_GENE_NAMES,
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
                match col:
                    case "chromosome":
                        outrow.append(core.check_unique(rows.chromosome, "chromosome"))
                    case "start":
                        outrow.append(rows.start.iat[0])
                    case "end":
                        outrow.append(rows.end.iat[-1])
                    case "gene":
                        outrow.append(name)
                    case "log2":
                        outrow.append(summary_func(rows.log2))
                    case "probes":
                        # Special case: sum probes rather than average
                        outrow.append(sum(rows[col]))
                    case _:
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
        is_xx: bool | bool_ | None = None,
        diploid_parx_genome: str | None = None,
    ) -> CopyNumArray:
        """Adjust chrX log2 ratios to match the ploidy of the reference sex.

        I.e. add 1 to chrX log2 ratios for a male sample vs. female reference,
        or subtract 1 for a female sample vs. male reference, so that chrX log2
        values are comparable across samples with different chromosomal sex.
        """
        outprobes = self.copy()
        if is_xx is None:
            # guess_xx may still return None (no chrX); default to female so an
            # indeterminate sample is never silently given the male +1 shift.
            is_xx = is_female_default(
                self.guess_xx(
                    is_haploid_x_reference=is_haploid_x_reference,
                    diploid_parx_genome=diploid_parx_genome,
                )
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
        diploid_parx_genome: str | None = None,
        verbose: bool = True,
    ) -> bool_ | None:
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
            # Decision is AND of the two maleness ratios -- show each ratio
            # separately so the log accurately reflects which axis carried
            # the call (a multiplicative product would be misleading here).
            logging.info(
                "Relative log2 coverage of %s=%.3g, %s=%.3g "
                "(maleness chrX=%.3g, chrY=%.3g) --> assuming %s",
                self.chr_x_label or "chrX",
                stats["chrx_ratio"],
                self.chr_y_label or "chrY",
                stats["chry_ratio"],
                stats["chrx_male_lr"],
                stats["chry_male_lr"],
                "male" if is_xy else "female",
            )
        return ~is_xy

    def compare_sex_chromosomes(
        self,
        is_haploid_x_reference: bool = False,
        diploid_parx_genome: str | None = None,
        skip_low: bool = False,
    ) -> (
        tuple[bool_, dict[str, float64 | float]]
        | tuple[bool_, dict[str, float64]]
        | tuple[None, dict[Any, Any]]
    ):
        """Compare coverage ratios of sex chromosomes versus autosomes.

        For each sex chromosome, compute a "maleness" ratio that asks how
        close the observed median log2 ratio (chromosome minus autosomes)
        sits to the male expectation versus the female expectation. The
        observed quantity is a *median difference* -- robust to noisy
        low-coverage bins on its own without further test-statistic
        machinery -- and the expected positions follow the standard
        haploid/diploid conventions for chrX (depending on
        ``is_haploid_x_reference``) and an absent-versus-present convention
        for chrY (female chrY at ``params.NULL_LOG2_COVERAGE``, male chrY
        at autosome-median). A ratio > 1 means the male hypothesis fits
        better; the 1.0 boundary is the midpoint between the two expected
        positions, so the decision is threshold-free up to that geometry.

        The two ratios combine with an AND-gate: the sample is male only
        when *both* chrX and chrY independently look more male than female.
        This is monotonic toward the safe female default (the design
        principle's "female default when no positive male evidence") and
        treats chrX and chrY symmetrically rather than letting one
        multiplicatively dominate the other -- chrY routinely beats chrX
        by 30-300x in real data, so a multiplicative combination let chrY
        decide the call (#954) and let representation noise on either axis
        flip the inference between ``.cnr`` and ``.cns`` (#785).

        Parameters
        ----------
        is_haploid_x_reference : bool
            Whether a male reference copy number profile was used to normalize
            the data. If so, a male sample should have log2 values of 0 on X
            and Y, and female +1 on X, deep negative on Y. Otherwise, a male
            sample should have log2 values of -1 on X and 0 on Y, and female
            0 on X, deep negative on Y.
        skip_low : bool
            If True, drop very-low-coverage bins (via `drop_low_coverage`)
            before comparing log2 coverage ratios. Included for completeness,
            but shouldn't affect the result much since the maleness ratios
            use medians.

        Returns
        -------
        bool
            True if the sample appears male.
        dict
            Calculated values used for the inference: relative log2 ratios of
            chromosomes X and Y versus the autosomes, and the per-chromosome
            maleness ratios.
        """
        if not len(self):
            return None, {}

        chrx = self[self.chr_x_filter(diploid_parx_genome)]
        if not len(chrx):
            # No usable chrX rows for inference. Two shapes both land here:
            # (a) the context-aware classifier
            # (``skgenome.chromnames.infer_sex_chrom_labels``) returned
            # ``chr_x_label is None`` -- the assembly has no chrX (yeast
            # and other Roman-numeral genomes, where ``chrX`` is autosome
            # 10; ZW-sex species; custom non-human assemblies, #669); and
            # (b) ``chr_x_label`` is set but the filter still yielded no
            # rows -- e.g. all chrX bins fall inside PAR under
            # ``--diploid-parx-genome``, or the data has been pre-filtered
            # to autosomes only. Both correctly map to "sex inference is
            # not applicable here"; the honest ``None`` propagates to
            # ``is_female_default`` at decision consumers. The pre-#669
            # WARNING ("is the input truncated?") falsely framed the
            # structurally-correct cases as data corruption.
            return None, {}

        auto = self.autosomes(diploid_parx_genome=diploid_parx_genome)
        if skip_low:
            chrx = chrx.drop_low_coverage()
            auto = auto.drop_low_coverage()
        use_weight = "weight" in self

        def _median(values, weights):
            if use_weight and len(values):
                return float(descriptives.weighted_median(values, weights))
            return float(np.median(values))

        auto_l = auto["log2"].to_numpy()
        auto_w = auto["weight"].to_numpy() if use_weight else None
        auto_med = _median(auto_l, auto_w)

        def _maleness(observed_med, female_shift, male_shift):
            """Ratio of residuals to the female vs. male expected positions.

            ``female_shift`` / ``male_shift`` are what would be added to a
            true-female / true-male observation to bring it to the autosome
            median, so the implied expected positions on the observed log2
            ratio scale are ``-female_shift`` and ``-male_shift``. Returns
            ``f_resid / m_resid``; > 1 ⇔ closer to the male expected position.
            """
            ratio = observed_med - auto_med
            f_resid = abs(ratio + female_shift)
            m_resid = abs(ratio + male_shift)
            # Tiny floor so the ratio doesn't explode to inf when an observed
            # median sits exactly at the male expected position.
            return f_resid / max(m_resid, 1e-6)

        chrx_l = chrx["log2"].to_numpy()
        chrx_w = chrx["weight"].to_numpy() if use_weight else None
        chrx_med = _median(chrx_l, chrx_w)
        female_x_shift, male_x_shift = (-1, 0) if is_haploid_x_reference else (0, +1)
        chrx_male_lr = _maleness(chrx_med, female_x_shift, male_x_shift)

        # Female chrY sits at NULL_LOG2_COVERAGE (the no-reads sentinel; that
        # negation flips it to a positive ``female_shift`` because _maleness
        # adds the shift to the observed ratio to align it with autosomes).
        # Male chrY sits at autosome-median (mapping-suppressed) or above. The
        # wide absence/presence separation makes the chrY maleness check act
        # as a presence/absence detector while keeping the 1.0 boundary
        # principled (midpoint between the two expectations).
        chry = self[self.chr_y_filter(diploid_parx_genome)]
        if skip_low and len(chry):
            chry = chry.drop_low_coverage()
        if len(chry):
            chry_l = chry["log2"].to_numpy()
            chry_w = chry["weight"].to_numpy() if use_weight else None
            chry_med = _median(chry_l, chry_w)
            chry_male_lr = _maleness(chry_med, -params.NULL_LOG2_COVERAGE, 0.0)
            is_male = chrx_male_lr > 1.0 and chry_male_lr > 1.0
            chry_ratio = chry_med - auto_med
        else:
            # No chrY data at all (or all dropped as low-coverage) -- can't
            # gate on Y, so fall back to chrX-only, same as the pre-redesign
            # behavior for this edge case.
            chry_male_lr = np.nan
            chry_ratio = np.nan
            is_male = chrx_male_lr > 1.0

        return (
            np.bool_(is_male),
            dict(
                chrx_ratio=chrx_med - auto_med,
                chry_ratio=chry_ratio,
                chrx_male_lr=chrx_male_lr,
                chry_male_lr=chry_male_lr,
            ),
        )

    def expect_flat_log2(
        self,
        is_haploid_x_reference: bool | None = None,
        diploid_parx_genome: str | None = None,
    ) -> ndarray:
        """Get the uninformed expected copy ratios of each bin.

        Create an array of log2 coverages like a "flat" reference.

        This is a neutral copy ratio at each autosome (log2 = 0.0) and sex
        chromosomes based on whether the reference is male (XX or XY).
        """
        if is_haploid_x_reference is None:
            # Default to a female (diploid-X) reference when sex is unknown,
            # rather than letting a None guess become a male (haploid-X) one.
            is_haploid_x_reference = not is_female_default(
                self.guess_xx(diploid_parx_genome=diploid_parx_genome, verbose=False)
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

    def residuals(self, segments: CopyNumArray | GenomicArray | None = None) -> Series:
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
                    ),
                    strict=True,
                )
                if len(bins_lr)
            ]
        else:
            resids = [
                lr - lr.median()
                for lr in self.iter_ranges_of(segments, "log2", keep_empty=False)
            ]
        return pd.concat(resids) if resids else pd.Series([])

    def smooth_log2(self, bandwidth: int | None = None, by_arm: bool = True) -> ndarray:
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
