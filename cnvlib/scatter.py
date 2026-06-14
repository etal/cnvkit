"""The 'scatter' command for rendering copy number as scatter plots."""

from __future__ import annotations

import collections
import logging
from typing import TYPE_CHECKING

import numpy as np
from matplotlib import pyplot

from skgenome.chromnames import detect_chr_prefix, diagnose_missing_chromosome
from skgenome.rangelabel import unpack_range

from . import core, params, plots
from .cnary import CopyNumArray as CNA
from .plots import MB

if TYPE_CHECKING:
    from matplotlib.axes._axes import Axes
    from matplotlib.figure import Figure

    from cnvlib.cnary import CopyNumArray

HIGHLIGHT_COLOR = "gold"
POINT_COLOR = "#606060"
SEG_COLOR = "darkorange"
TREND_COLOR = "#A0A0A0"
# Distinct colors and markers for ``--show-snvs`` overlays drawn on top of
# the het SNP scatter. Colors picked from Wong's colorblind-friendly palette
# so they remain discriminable from each other, from the gray het dots, and
# from the darkorange segment overlay (Wong 2011, Nature Methods 8:441).
# Markers are also distinct so the overlays read in grayscale prints. LOH
# gets the higher-contrast reddish-purple + 'x' to ensure the typically rare
# and low-VAF LOH dots are not visually buried by the bulk het scatter. (#290)
LOH_SNV_COLOR = "#CC79A7"  # reddish purple
LOH_SNV_MARKER = "x"
SOMATIC_SNV_COLOR = "#0072B2"  # blue
SOMATIC_SNV_MARKER = "+"

# Floor for the auto-scaled y-axis lower limit, so a single deep homozygous
# deletion (log2 far below 0) can't compress the rest of the plot into a sliver.
# Segments below this are flagged; override with the --y-min option. (#385)
AUTO_Y_MIN_FLOOR = -5.0


def _sex_labels(arr) -> tuple[str | None, str | None]:
    """Best-effort ``(x_label, y_label)`` from any genomic array.

    Returns ``(None, None)`` for arrays without the labels (e.g. plain
    ``GenomicArray`` or ``VariantArray``) so callers can pass the result
    to ``choose_segment_color`` unconditionally.
    """
    return (
        getattr(arr, "chr_x_label", None),
        getattr(arr, "chr_y_label", None),
    )


def do_scatter(
    cnarr: CopyNumArray,
    segments: CopyNumArray | None = None,
    variants: CopyNumArray | None = None,
    show_range: str | None = None,
    show_gene: str | None = None,
    do_trend: bool = False,
    by_bin: bool = False,
    window_width: float = 1e6,
    y_min: float | None = None,
    y_max: float | None = None,
    fig_size: tuple[float, float] | None = None,
    antitarget_marker: str | None = None,
    segment_color: str = SEG_COLOR,
    title: str | None = None,
    loh_variants: CopyNumArray | None = None,
    somatic_variants: CopyNumArray | None = None,
) -> Figure:
    """Plot probe log2 coverages and segmentation calls together.

    Parameters
    ----------
    cnarr : CopyNumArray
        Bin-level copy number data to plot.
    segments : CopyNumArray, optional
        Segmented copy number data to overlay.
    variants : VariantArray, optional
        Variant allele frequency data to plot.
    show_range : str, optional
        Genomic range to display (e.g., "chr1:1000000-2000000").
    show_gene : str, optional
        Gene name to zoom into.
    do_trend : bool, optional
        Plot a smoothed trendline. Default is False.
    by_bin : bool, optional
        Plot by bin index rather than genomic coordinates. Default is False.
    window_width : float, optional
        Smoothing window width in bases (or bins if by_bin=True).
        Default is 1e6.
    y_min : float, optional
        Minimum y-axis value.
    y_max : float, optional
        Maximum y-axis value.
    fig_size : tuple of float, optional
        Figure size as (width, height) in inches.
    antitarget_marker : str, optional
        Marker style for antitarget bins.
    segment_color : str, optional
        Color for segment lines. Default is SEG_COLOR.
    title : str, optional
        Plot title.
    loh_variants : VariantArray, optional
        Tumor-homozygous loci to overlay as LOH evidence in the VAF panel
        (#290). Plotted with a distinct color; excluded from the BAF
        segment-overlay trend so the trend reflects only the het subset.
    somatic_variants : VariantArray, optional
        Somatic-SNV loci (VCF SOMATIC flag or T/N-inferred) to overlay in
        the VAF panel with a distinct color (#290). Also excluded from the
        BAF segment-overlay trend.

    Returns
    -------
    matplotlib.figure.Figure
        The scatter plot figure.
    """
    if by_bin:
        bp_per_bin = sum(c.end.iat[-1] for _, c in cnarr.by_chromosome()) / len(cnarr)
        window_width /= bp_per_bin
        show_range_bins = plots.translate_region_to_bins(show_range, cnarr)
        (
            cnarr,
            segments,
            variants,
            (
                loh_variants,
                somatic_variants,
            ),
        ) = plots.update_binwise_positions(
            cnarr, segments, variants, extra_variants=[loh_variants, somatic_variants]
        )
        global MB
        orig_mb = MB
        MB = 1

    if not show_gene and not show_range:
        fig = genome_scatter(
            cnarr,
            segments,
            variants,
            do_trend,
            y_min,
            y_max,
            title,
            segment_color,
            loh_variants=loh_variants,
            somatic_variants=somatic_variants,
        )
    else:
        if by_bin:
            show_range = show_range_bins
        fig = chromosome_scatter(
            cnarr,
            segments,
            variants,
            show_range,
            show_gene,
            antitarget_marker,
            do_trend,
            by_bin,
            window_width,
            y_min,
            y_max,
            title,
            segment_color,
            loh_variants=loh_variants,
            somatic_variants=somatic_variants,
        )

    if by_bin:
        # Reset to avoid permanently altering the value of cnvlib.scatter.MB
        MB = orig_mb
    if fig_size:
        width, height = fig_size
        fig.set_size_inches(w=width, h=height)
    return fig


# === Genome-level scatter plots ===


def genome_scatter(
    cnarr: CopyNumArray,
    segments: CopyNumArray | None = None,
    variants: CopyNumArray | None = None,
    do_trend: bool = False,
    y_min: float | None = None,
    y_max: float | None = None,
    title: str | None = None,
    segment_color: str = SEG_COLOR,
    loh_variants: CopyNumArray | None = None,
    somatic_variants: CopyNumArray | None = None,
) -> Figure:
    """Plot all chromosomes, concatenated on one plot."""
    if (cnarr or segments) and variants:
        # Lay out top 3/5 for the CN scatter, bottom 2/5 for SNP plot
        axgrid = pyplot.GridSpec(5, 1, hspace=0.85)
        axis = pyplot.subplot(axgrid[:3])
        axis2 = pyplot.subplot(axgrid[3:], sharex=axis)
        # Place chromosome labels between the CNR and SNP plots
        axis2.tick_params(labelbottom=False)
        chrom_sizes = plots.chromosome_sizes(cnarr or segments)  # type: ignore[arg-type]
        axis2 = snv_on_genome(
            axis2,
            variants,
            chrom_sizes,
            segments,
            do_trend,
            segment_color,
            loh_variants=loh_variants,
            somatic_variants=somatic_variants,
        )
    else:
        _fig, axis = pyplot.subplots()
    if title is None:
        title = (cnarr or segments or variants).sample_id  # type: ignore[union-attr]
    if cnarr or segments:
        axis.set_title(title)
        axis = cnv_on_genome(
            axis,
            cnarr,  # type: ignore[arg-type]
            segments,  # type: ignore[arg-type]
            do_trend,
            y_min,
            y_max,
            segment_color,  # type: ignore[arg-type]
        )
    else:
        axis.set_title(f"Variant allele frequencies: {title}")
        chrom_sizes = collections.OrderedDict(
            (chrom, subarr["end"].max())
            for chrom, subarr in variants.by_chromosome()  # type: ignore[union-attr]
        )
        axis = snv_on_genome(
            axis,
            variants,
            chrom_sizes,
            segments,
            do_trend,
            segment_color,
            loh_variants=loh_variants,
            somatic_variants=somatic_variants,
        )
    return axis.get_figure()  # type: ignore[return-value, no-any-return]


def cnv_on_genome(
    axis: Axes,
    probes: CopyNumArray,
    segments: CopyNumArray,
    do_trend: bool = False,
    y_min: float | None = None,
    y_max: float | None = None,
    segment_color: str = SEG_COLOR,
) -> Axes:
    """Plot bin ratios and/or segments for all chromosomes on one plot."""
    # Configure axes etc.
    axis.axhline(color="k")
    axis.set_ylabel("Copy ratio (log2)")
    if not (y_min and y_max):
        if segments:
            # Auto-scale y-axis according to segment mean-coverage values,
            # excluding chrY (low ploidy) and chr6 (HLA region's noisy log2)
            # to avoid pulling the scale to spuriously low values. For non-human
            # data without sex chromosomes these exclusions become a no-op.
            prefix = detect_chr_prefix(segments.chromosome.unique())
            exclude_for_scaling = {f"{prefix}6"}
            chr_y_for_scaling = segments.chr_y_label
            if chr_y_for_scaling is not None:
                exclude_for_scaling.add(chr_y_for_scaling)
            low_chroms = segments.chromosome.isin(exclude_for_scaling)
            seg_auto_vals = segments[~low_chroms]["log2"].dropna()
            if not y_min:
                y_min = (
                    max(AUTO_Y_MIN_FLOOR, np.nanmin([seg_auto_vals.min() - 0.2, -1.5]))
                    if len(seg_auto_vals)
                    else -2.5
                )
            if not y_max:
                y_max = (
                    np.nanmax([seg_auto_vals.max() + 0.2, 1.5])
                    if len(seg_auto_vals)
                    else 2.5
                )
        else:
            if not y_min:
                y_min = -2.5
            if not y_max:
                y_max = 2.5
    axis.set_ylim(y_min, y_max)
    if segments:
        # Flag segments clipped below the (possibly floored) y-axis, so a deep
        # homozygous deletion isn't silently dropped off the bottom (#385).
        hidden_seg = segments["log2"] < y_min
        if hidden_seg.sum():
            logging.warning(
                "With y_min=%s, %d genome-wide segment(s) are hidden below the "
                "plot --> add '--y-min %s' to see them",
                y_min,
                int(hidden_seg.sum()),
                int(np.floor(segments["log2"].min())),
            )

    # Group probes by chromosome (to calculate plotting coordinates)
    if probes:
        chrom_sizes = plots.chromosome_sizes(probes)
        chrom_probes = dict(probes.by_chromosome())
    else:
        chrom_sizes = plots.chromosome_sizes(segments)
    # Same for segment calls
    chrom_segs = dict(segments.by_chromosome()) if segments else {}

    # Plot points & segments
    x_starts = plots.plot_chromosome_dividers(axis, chrom_sizes)
    sex_labels = _sex_labels(segments)
    for chrom, x_offset in x_starts.items():
        if probes and chrom in chrom_probes:
            subprobes = chrom_probes[chrom]
            x = 0.5 * (subprobes["start"] + subprobes["end"]) + x_offset
            axis.scatter(
                x,
                subprobes["log2"],
                marker=".",
                color=POINT_COLOR,
                edgecolor="none",
                alpha=0.2,
            )
            if do_trend:
                # ENH break trendline by chromosome arm boundaries?
                # Here and in subsequent occurrences, it's important to use snap=False
                # to avoid short lines/segment disappearing when saving as PNG.
                # See also: https://github.com/etal/cnvkit/issues/604
                axis.plot(
                    x,
                    subprobes.smooth_log2(),
                    color=POINT_COLOR,
                    linewidth=2,
                    zorder=-1,
                    snap=False,
                )

        if chrom in chrom_segs:
            for seg in chrom_segs[chrom]:
                color = choose_segment_color(
                    seg, segment_color, sex_chrom_labels=sex_labels
                )
                axis.plot(
                    (seg.start + x_offset, seg.end + x_offset),
                    (seg.log2, seg.log2),
                    color=color,
                    linewidth=3,
                    solid_capstyle="round",
                    snap=False,
                )
    return axis


def snv_on_genome(
    axis,
    variants,
    chrom_sizes,
    segments,
    do_trend,
    segment_color,
    loh_variants=None,
    somatic_variants=None,
):
    """Plot a scatter-plot of SNP chromosomal positions and shifts.

    Optional ``loh_variants`` and ``somatic_variants`` are overlaid in
    distinct colors and excluded from the segment-overlay trend so the trend
    continues to reflect only the het subset (#290).
    """
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("VAF")
    x_starts = plots.plot_chromosome_dividers(axis, chrom_sizes)

    # Calculate the coordinates of plot components
    chrom_snvs = dict(variants.by_chromosome())
    chrom_loh = dict(loh_variants.by_chromosome()) if loh_variants else {}
    chrom_som = dict(somatic_variants.by_chromosome()) if somatic_variants else {}
    if segments:
        chrom_segs = dict(segments.by_chromosome())
    elif do_trend:
        # Pretend a single segment covers each chromosome
        chrom_segs = dict.fromkeys(chrom_snvs)
    else:
        chrom_segs = {}

    sex_labels = _sex_labels(segments)
    for chrom, x_offset in x_starts.items():
        if chrom not in chrom_snvs:
            continue

        snvs = chrom_snvs[chrom]
        # Plot the points
        axis.scatter(
            snvs["start"].values + x_offset,
            snvs["alt_freq"].values,
            color=POINT_COLOR,
            edgecolor="none",
            alpha=0.2,
            marker=".",
        )
        # Overlay LOH-evidence and somatic markers on top of the het dots.
        # Distinct color AND marker shape (set in module-level constants) so
        # the overlays read in grayscale prints as well as in color. Higher
        # alpha than het dots so they remain visible through dense het
        # regions. Legend handles (``label=``) are attached only on the
        # first per-chromosome draw to avoid N duplicate legend entries; the
        # axis-level legend call at the end of the loop emits the legend.
        if chrom in chrom_loh:
            loh_snvs = chrom_loh[chrom]
            axis.scatter(
                loh_snvs["start"].values + x_offset,
                loh_snvs["alt_freq"].values,
                color=LOH_SNV_COLOR,
                alpha=0.7,
                marker=LOH_SNV_MARKER,
                label=_overlay_label_once(axis, "LOH"),
            )
        if chrom in chrom_som:
            som_snvs = chrom_som[chrom]
            axis.scatter(
                som_snvs["start"].values + x_offset,
                som_snvs["alt_freq"].values,
                color=SOMATIC_SNV_COLOR,
                alpha=0.7,
                marker=SOMATIC_SNV_MARKER,
                label=_overlay_label_once(axis, "somatic"),
            )
        # Trend bars: always calculated, only shown on request.
        # NB: trend uses only the het subset (snvs), never the overlays --
        # the BAF math that drives the segment overlay must be unchanged.
        if chrom in chrom_segs:
            # Draw average VAF within each segment
            segs = chrom_segs[chrom]
            for seg, v_freq in get_segment_vafs(snvs, segs):
                if seg:
                    posn = [seg.start + x_offset, seg.end + x_offset]
                    color = choose_segment_color(
                        seg,
                        segment_color,
                        default_bright=False,
                        sex_chrom_labels=sex_labels,
                    )
                else:
                    posn = [snvs.start.iat[0] + x_offset, snvs.start.iat[-1] + x_offset]
                    color = TREND_COLOR
                axis.plot(
                    posn,
                    [v_freq, v_freq],
                    color=color,
                    linewidth=2,
                    zorder=-1,
                    solid_capstyle="round",
                    snap=False,
                )
    if loh_variants is not None or somatic_variants is not None:
        # Legend lives outside the plot (anchored to the right edge of the
        # VAF panel) so it doesn't occlude the dot field. Marker keys are
        # de-duplicated by _overlay_label_once. ``bbox_inches="tight"`` at
        # savefig captures the legend even though it sits past the axes.
        axis.legend(
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0,
            fontsize="small",
            frameon=False,
        )
    return axis


def _overlay_label_once(axis: Axes, label: str) -> str | None:
    """Return ``label`` the first time it's seen for this axis, else None.

    Prevents matplotlib's legend from accumulating N duplicate entries when
    the same overlay kind is plotted once per chromosome in genome view.
    Uses the axis itself as the dedupe scope so different subplots get
    independent legends.
    """
    seen = getattr(axis, "_cnvkit_overlay_labels_seen", None)
    if seen is None:
        seen = set()
        axis._cnvkit_overlay_labels_seen = seen  # type: ignore[attr-defined]
    if label in seen:
        return None
    seen.add(label)
    return label


# === Chromosome-level scatter plots ===


def chromosome_scatter(
    cnarr,
    segments,
    variants,
    show_range,
    show_gene,
    antitarget_marker,
    do_trend,
    by_bin,
    window_width,
    y_min,
    y_max,
    title,
    segment_color,
    loh_variants=None,
    somatic_variants=None,
):
    """Plot a specified region on one chromosome.

    Possibilities::

             Options | Shown
        ------------ | --------
        -c      | -g | Genes | Region
        ------- | -- | ----- | ------
        -       | +  | given | auto: gene(s) + margin
        chr     | -  | none  | whole chrom
        chr     | +  | given | whole chrom
        chr:s-e | -  | all   | given
        chr:s-e | +  | given | given

    """
    (
        sel_probes,
        sel_segs,
        sel_snvs,
        sel_loh,
        sel_som,
        window_coords,
        genes,
        chrom,
    ) = select_range_genes(
        cnarr,
        segments,
        variants,
        show_range,
        show_gene,
        window_width,
        loh_variants=loh_variants,
        somatic_variants=somatic_variants,
    )
    # Create plots
    if cnarr or segments:
        # Plot CNVs at chromosome level
        if variants:
            # Lay out top 3/5 for the CN scatter, bottom 2/5 for SNP plot
            axgrid = pyplot.GridSpec(5, 1, hspace=0.5)
            axis = pyplot.subplot(axgrid[:3])
            axis2 = pyplot.subplot(axgrid[3:], sharex=axis)
            # Plot allele freqs for only the selected region
            snv_on_chromosome(
                axis2,
                sel_snvs,
                sel_segs,
                genes,
                do_trend,
                by_bin,
                segment_color,
                loh_variants=sel_loh,
                somatic_variants=sel_som,
            )
        else:
            _fig, axis = pyplot.subplots()
            if by_bin:
                axis.set_xlabel("Position (bin)")
            else:
                axis.set_xlabel("Position (Mb)")
        axis = cnv_on_chromosome(
            axis,
            sel_probes,
            sel_segs,
            genes,
            antitarget_marker=antitarget_marker,
            do_trend=do_trend,
            x_limits=window_coords,
            y_min=y_min,
            y_max=y_max,
            segment_color=segment_color,
        )
    elif variants:
        # Only plot SNVs in a single-panel layout
        _fig, axis = pyplot.subplots()
        axis = snv_on_chromosome(
            axis,
            sel_snvs,
            sel_segs,
            genes,
            do_trend,
            by_bin,
            segment_color,
            loh_variants=sel_loh,
            somatic_variants=sel_som,
        )

    if title is None:
        title = "%s %s" % ((cnarr or segments or variants).sample_id, chrom)
    axis.set_title(title)
    return axis.get_figure()


def select_range_genes(
    cnarr,
    segments,
    variants,
    show_range,
    show_gene,
    window_width,
    *,
    loh_variants=None,
    somatic_variants=None,
):
    """Determine which datapoints to show based on the given options.

    Behaviors::

        start/end   show_gene
           +           +       given region + genes; err if any gene outside it
           -           +       window +/- around genes
           +           -       given region, highlighting any genes within it
           -           -       whole chromosome, no genes

    If `show_range` is a chromosome name only, no start/end positions, then the
    whole chromosome will be shown.

    If region start/end coordinates are given and `show_gene` is '' or ',' (or
    all commas, etc.), then instead of highlighting all genes in the selection,
    no genes will be highlighted.
    """
    chrom, start, end = unpack_range(show_range)
    window_coords: tuple = ()
    if start is None and end is None:
        # Either the specified range is only chrom, no start-end, or gene names
        # were given
        pass
    else:
        # Viewing region coordinates were specified -- take them as given
        # Fill in open-ended ranges' endpoints
        if start is None or start < 0:
            start = 0
        if not end:
            # Default selection endpoint to the maximum chromosome position
            end = (cnarr or segments or variants).filter(chromosome=chrom).end.iat[-1]
        if end <= start:
            raise ValueError(
                f"Coordinate range {chrom}:{start}-{end} (from {show_range}) "
                + "has size <= 0"
            )
        window_coords = (start, end)

    gene_ranges = []
    if show_gene is None:
        if window_coords:
            if cnarr:
                # Highlight all genes within the given range
                gene_ranges = plots.gene_coords_by_range(cnarr, chrom, start, end)[
                    chrom
                ]
            if not gene_ranges and (end - start) < 10 * window_width:
                # No genes in the selected region, so if the selection is small
                # (i.e. <80% of the displayed window, <10x window padding),
                # highlight the selected region itself.
                # (To prevent this, use show_gene='' or window_width=0)
                logging.info(
                    "No genes found in selection; will highlight the "
                    "selected region itself instead"
                )
                gene_ranges = [(start, end, "Selection")]
                window_coords = (max(0, start - window_width), end + window_width)

    else:
        gene_names = filter(None, show_gene.split(","))
        if gene_names:
            # Scan for probes matching the specified gene(s)
            gene_coords = plots.gene_coords_by_name(cnarr or segments, gene_names)
            if len(gene_coords) > 1:
                raise ValueError(
                    f"Genes {show_gene} are split across chromosomes "
                    f"{list(gene_coords.keys())}"
                )
            g_chrom, gene_ranges = gene_coords.popitem()
            if chrom:
                # Confirm that the selected chromosomes match
                core.assert_equal(
                    "Chromosome also selected by region (-c) does not match",
                    **{"chromosome": chrom, "gene(s)": g_chrom},
                )
            else:
                chrom = g_chrom

            gene_ranges.sort()
            if window_coords:
                # Verify all genes fit in the given window
                for gene_start, gene_end, gene_name in gene_ranges:
                    if not (start <= gene_start and gene_end <= end):
                        raise ValueError(
                            f"Selected gene {gene_name} "
                            + f"({chrom}:{gene_start}-{gene_end}) "
                            + f"is outside specified region {show_range}"
                        )
            elif not show_range:
                # Set the display window to the selected genes +/- a margin
                window_coords = (
                    max(0, gene_ranges[0][0] - window_width),
                    gene_ranges[-1][1] + window_width,
                )

    # Prune plotted elements to the selected region
    sel_probes = cnarr.in_range(chrom, *window_coords) if cnarr else CNA([])
    sel_segs = (
        segments.in_range(chrom, *window_coords, mode="trim") if segments else CNA([])
    )
    sel_snvs = variants.in_range(chrom, *window_coords) if variants else None
    sel_loh = loh_variants.in_range(chrom, *window_coords) if loh_variants else None
    sel_som = (
        somatic_variants.in_range(chrom, *window_coords) if somatic_variants else None
    )
    logging.info(
        "Showing %d probes and %d selected genes in region %s",
        len(sel_probes),
        len(gene_ranges),
        (chrom + ":{}-{}".format(*window_coords) if window_coords else chrom),
    )
    if cnarr and len(sel_probes) == 0:
        # Help the user diagnose why the selection has no probes
        # (e.g. yeast Roman-numeral names not present in the .cnr file,
        # a "chr" prefix mismatch, or a window outside any probe).
        logging.warning(
            "%s",
            diagnose_missing_chromosome(chrom, cnarr.chromosome.unique()),
        )

    return (
        sel_probes,
        sel_segs,
        sel_snvs,
        sel_loh,
        sel_som,
        window_coords,
        gene_ranges,
        chrom,
    )


def cnv_on_chromosome(
    axis,
    probes,
    segments,
    genes,
    antitarget_marker=None,
    do_trend=False,
    x_limits=None,
    y_min=None,
    y_max=None,
    segment_color=SEG_COLOR,
):
    """Draw a scatter plot of probe values with optional segments overlaid.

    Parameters
    ----------
    genes : list
        Of tuples: (start, end, gene name)
    """
    # TODO - allow plotting just segments without probes
    # Get scatter plot coordinates
    x = 0.5 * (probes["start"] + probes["end"]) * MB  # bin midpoints
    y = probes["log2"]
    w = 46 * probes["weight"] ** 2 + 2 if "weight" in probes else np.repeat(30, len(x))

    # Configure axes
    if not y_min:
        y_min = max(AUTO_Y_MIN_FLOOR, min(y.min() - 0.1, -0.3)) if len(y) else -1.1
    if not y_max:
        y_max = max(0.3, y.max() + (0.25 if genes else 0.1)) if len(y) else 1.1
    if x_limits:
        x_min, x_max = x_limits
        axis.set_xlim(x_min * MB, x_max * MB)
    else:
        set_xlim_from(axis, probes, segments)
    setup_chromosome(axis, y_min, y_max, "Copy ratio (log2)")
    if genes:
        highlight_genes(axis, genes, min(2.4, y.max() + 0.1) if len(y) else 0.1)

    if antitarget_marker in (None, "o"):
        # Plot targets and antitargets with the same marker
        axis.scatter(x, y, w, color=POINT_COLOR, alpha=0.4, marker="o")
    else:
        # Use the given marker to plot antitargets separately
        x_fg = []
        y_fg = []
        w_fg = []
        x_bg = []
        y_bg = []
        # w_bg = []
        is_bg = probes["gene"].isin(params.ANTITARGET_ALIASES)
        for x_pt, y_pt, w_pt, is_bg_pt in zip(x, y, w, is_bg, strict=True):
            if is_bg_pt:
                x_bg.append(x_pt)
                y_bg.append(y_pt)
                # w_bg.append(w_pt)
            else:
                x_fg.append(x_pt)
                y_fg.append(y_pt)
                w_fg.append(w_pt)
        axis.scatter(x_fg, y_fg, w_fg, color=POINT_COLOR, alpha=0.4, marker="o")
        axis.scatter(x_bg, y_bg, color=POINT_COLOR, alpha=0.5, marker=antitarget_marker)

    # Add a local trend line
    if do_trend:
        axis.plot(
            x,
            probes.smooth_log2(),  # .1),
            color=POINT_COLOR,
            linewidth=2,
            zorder=-1,
            snap=False,
        )

    # Draw segments as horizontal lines
    if segments:
        sex_labels = _sex_labels(segments)
        for row in segments:
            color = choose_segment_color(
                row, segment_color, sex_chrom_labels=sex_labels
            )
            axis.plot(
                (row.start * MB, row.end * MB),
                (row.log2, row.log2),
                color=color,
                linewidth=4,
                solid_capstyle="round",
                snap=False,
            )
        # Warn about segments masked by default pruning at 'y_min=-5':
        hidden_seg = segments.log2 < y_min
        if hidden_seg.sum():
            logging.warning(
                "WARNING: With 'y_min=%s' %s segments are hidden"
                " --> Add parameter '--y-min %s' to see them",
                y_min,
                hidden_seg.sum(),
                np.floor(segments.log2.min()).astype(int),
            )
            # Signal them as triangles crossing y-axis:
            x_hidden = segments.start[hidden_seg] * MB
            y_hidden = np.array([y_min] * len(x_hidden))
            axis.scatter(
                x_hidden,
                y_hidden,
                marker="^",
                linewidth=3,
                snap=False,
                color=segment_color,
                edgecolor="none",
                clip_on=False,
                zorder=10,
            )
    return axis


def snv_on_chromosome(
    axis,
    variants,
    segments,
    genes,
    do_trend,
    by_bin,
    segment_color,
    loh_variants=None,
    somatic_variants=None,
):
    """Draw VAFs for a single-chromosome region.

    Optional ``loh_variants`` and ``somatic_variants`` are overlaid in
    distinct colors and excluded from the segment-overlay trend so the trend
    continues to reflect only the het subset (#290).
    """
    # TODO set x-limits if not already done for probes/segments
    # set_xlim_from(axis, None, segments, variants)
    # setup_chromosome(axis, 0.0, 1.0, "VAF")
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("VAF")
    if by_bin:
        axis.set_xlabel("Position (bin)")
    else:
        axis.set_xlabel("Position (Mb)")
    axis.get_yaxis().tick_left()
    axis.get_xaxis().tick_top()
    axis.tick_params(which="both", direction="out", labelbottom=False, labeltop=False)

    x_mb = variants["start"].values * MB
    y = variants["alt_freq"].values
    axis.scatter(x_mb, y, color=POINT_COLOR, alpha=0.3)
    overlay_present = False
    if loh_variants is not None and len(loh_variants):
        axis.scatter(
            loh_variants["start"].values * MB,
            loh_variants["alt_freq"].values,
            color=LOH_SNV_COLOR,
            alpha=0.8,
            marker=LOH_SNV_MARKER,
            label="LOH",
        )
        overlay_present = True
    if somatic_variants is not None and len(somatic_variants):
        axis.scatter(
            somatic_variants["start"].values * MB,
            somatic_variants["alt_freq"].values,
            color=SOMATIC_SNV_COLOR,
            alpha=0.8,
            marker=SOMATIC_SNV_MARKER,
            label="somatic",
        )
        overlay_present = True
    if overlay_present:
        axis.legend(
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0,
            fontsize="small",
            frameon=False,
        )
    if segments or do_trend:
        # Draw average VAF within each segment.
        # NB: trend uses only the het subset (variants), never the overlays
        # -- the BAF math driving the segment overlay must stay unchanged.
        sex_labels = _sex_labels(segments) if segments else (None, None)
        for seg, v_freq in get_segment_vafs(variants, segments):
            if seg:
                posn = [seg.start * MB, seg.end * MB]
                color = choose_segment_color(
                    seg,
                    segment_color,
                    default_bright=False,
                    sex_chrom_labels=sex_labels,
                )
            else:
                posn = [variants.start.iat[0] * MB, variants.start.iat[-1] * MB]
                color = TREND_COLOR
            axis.plot(
                posn,
                [v_freq, v_freq],
                color=color,
                linewidth=2,
                zorder=1,
                solid_capstyle="round",
                snap=False,
            )

    if genes:
        highlight_genes(axis, genes, 0.9)
    return axis


def set_xlim_from(axis, probes=None, segments=None, variants=None) -> None:
    """Configure axes for plotting a single chromosome's data.

    Parameters
    ----------
    probes : CopyNumArray
    segments : CopyNumArray
    variants : VariantArray
        All should already be subsetted to the region that will be plotted.
    """
    min_x = np.inf
    max_x = 0
    for arr in (probes, segments, variants):
        if arr and len(arr):
            max_x = max(max_x, arr.end.iat[-1])
            min_x = min(min_x, arr.start.iat[0])
    if max_x <= min_x:
        if min_x != np.inf:
            logging.warning(
                "WARNING: selection start %s > end %s; did you "
                "correctly sort the input file by genomic "
                "location?",
                min_x,
                max_x,
            )
        raise ValueError(
            "No usable data points to plot out of "
            f"{len(probes) if probes else 0} probes, "
            f"{len(segments) if segments else 0} segments, "
            f"{len(variants) if variants else 0} variants"
        )
    axis.set_xlim(min_x * MB, max_x * MB)


def setup_chromosome(axis, y_min=None, y_max=None, y_label=None) -> None:
    """Configure axes for plotting a single chromosome's data."""
    if y_min and y_max:
        axis.set_ylim(y_min, y_max)
        if y_min < 0 < y_max:
            axis.axhline(color="k")
    if y_label:
        axis.set_ylabel(y_label)
    axis.tick_params(which="both", direction="out")
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()


# === Shared ===


def choose_segment_color(
    segment,
    highlight_color,
    default_bright=True,
    sex_chrom_labels: tuple[str | None, str | None] = (None, None),
):
    """Choose a display color based on a segment's CNA status.

    Uses the fields added by the 'call' command. If these aren't present, use
    `highlight_color` for everything.

    For sex chromosomes, some single-copy deletions or gains might not be
    highlighted, since sample sex isn't used to infer the neutral ploidies.

    *sex_chrom_labels* is ``(x_label, y_label)`` for the genome being plotted;
    either may be None when the genome has no detected sex chromosome (e.g.
    yeast). When both are None, every chromosome is treated as autosomal
    (ploidy 2 expected). The labels must use the same naming convention as
    ``segment.chromosome`` (both ``"chrX"`` or both ``"X"``) so the dict
    lookup matches; in CNVkit's own pipeline this is guaranteed because the
    labels are derived from the array via ``infer_sex_chrom_labels``.
    """
    neutral_color = TREND_COLOR
    if "cn" not in segment._fields:
        # No 'call' info
        return highlight_color if default_bright else neutral_color

    x_label, y_label = sex_chrom_labels
    expected_ploidies: dict[str, tuple[int, int]] = {}
    if x_label is not None:
        expected_ploidies[x_label] = (1, 2)
    if y_label is not None:
        expected_ploidies[y_label] = (0, 1)

    # Detect copy number alteration
    if segment.cn not in expected_ploidies.get(segment.chromosome, (2,)):
        return highlight_color

    # Detect regions of allelic imbalance / LOH on autosomes
    if (
        segment.chromosome not in expected_ploidies
        and "cn1" in segment._fields
        and "cn2" in segment._fields
        and (segment.cn1 != segment.cn2)
    ):
        return highlight_color

    return neutral_color


def get_segment_vafs(variants, segments):
    """Group SNP allele frequencies by segment.

    Assume variants and segments were already subset to one chromosome.

    Yields
    ------
    tuple
        (segment, value)
    """
    # No segments -> one fake segment covering the whole region.
    chunks = variants.by_ranges(segments) if segments else [(None, variants)]
    for seg, seg_snvs in chunks:
        # ENH: seg_snvs.tumor_boost()
        freqs = seg_snvs["alt_freq"].values
        # Separately emit VAFs above and below .5 for plotting
        idx_above_mid = freqs > 0.5
        for idx_vaf in (idx_above_mid, ~idx_above_mid):
            if sum(idx_vaf) > 1:
                yield (seg, np.median(freqs[idx_vaf]))


def highlight_genes(axis, genes, y_posn) -> None:
    """Show gene regions with background color and a text label."""
    # Rotate text in proportion to gene density
    ngenes = len(genes)
    text_size = "small" if ngenes <= 6 else "x-small"
    text_rot: str | int
    if ngenes <= 3:
        text_rot = "horizontal"
    elif ngenes <= 6:
        text_rot = 30
    elif ngenes <= 10:
        text_rot = 45
    elif ngenes <= 20:
        text_rot = 60
    else:
        text_rot = "vertical"
    for gene in genes:
        gene_start, gene_end, gene_name = gene
        # Highlight and label gene region
        # (rescale positions from bases to megabases)
        axis.axvspan(
            gene_start * MB, gene_end * MB, alpha=0.5, color=HIGHLIGHT_COLOR, zorder=-1
        )
        axis.text(
            0.5 * (gene_start + gene_end) * MB,
            y_posn,
            gene_name,
            horizontalalignment="center",
            rotation=text_rot,
            size=text_size,
        )
