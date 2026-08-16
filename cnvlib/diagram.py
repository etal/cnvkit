"""Chromosome diagram drawing functions.

This uses and abuses Biopython's BasicChromosome module. It depends on
ReportLab, too, so we isolate this functionality here so that the rest of
CNVkit will run without it. (And also to keep the codebase tidy.)
"""

from __future__ import annotations

import collections
import math
import warnings
from typing import TYPE_CHECKING, Any, NamedTuple

# Reportlab triggers a DeprecationWarning via load_module on import, which can
# become an error under `-W error`. Silence it for just these imports rather
# than muzzling DeprecationWarning process-wide (which would hide every other
# deprecation, including CNVkit's own).
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    from Bio.Graphics import BasicChromosome as BC
    from reportlab.graphics import renderPDF
    from reportlab.lib import colors
    from reportlab.lib.units import inch
    from reportlab.pdfgen import canvas

import numpy as np

from skgenome.rangelabel import unpack_range

from . import descriptives, params, plots, reports
from .segmetrics import segment_mean

if TYPE_CHECKING:
    from collections.abc import Iterator

    from reportlab.graphics.shapes import Drawing

    from cnvlib.cnary import CopyNumArray

TELOMERE_LENGTH = 6e6  # For illustration only
CHROM_FATNESS = 0.3
PAGE_SIZE = (11.0 * inch, 8.5 * inch)
# Page layout, shared by `bc_organism_draw` and by the resolution the bin half
# is aggregated to. Every feature is stroked a point wide whatever its span,
# so the points a chromosome gets are also the features it can distinguish.
MARGIN_TOP = 1.25 * inch
MARGIN_BOTTOM = 0.1 * inch
MARGIN_SIDE = 0.5 * inch
CHROM_WRAP = 12  # Chromosomes per row
# Gene-column values that name no gene: placeholders and off-target bins
_UNNAMED = frozenset(params.IGNORE_GENE_NAMES + params.ANTITARGET_ALIASES)


def _bp_per_point(chrom_sizes: collections.OrderedDict) -> float:
    """Base pairs the page gives one drawn point, on the tightest row.

    Biopython divides a chromosome's vertical space by ``scale_num``, which
    ``build_chrom_diagram`` sets to the longest chromosome plus its two
    telomeres, and applies a 0.95 factor to the row's top coordinate. The rows
    come from ``bc_organism_draw`` above, so the arithmetic belongs here rather
    than in a measured constant: a whole genome gets ~220 points per scale and
    a single contig ~490, and a constant fitted to one of those aggregates the
    other twice as coarsely as the page can draw.
    """
    nrows = math.ceil(len(chrom_sizes) / CHROM_WRAP)
    row_height = (PAGE_SIZE[1] - MARGIN_TOP - MARGIN_BOTTOM) / nrows
    points = min(
        (PAGE_SIZE[1] - MARGIN_TOP - row_height * row) * 0.95
        - (MARGIN_BOTTOM + row_height * (nrows - row - 1))
        for row in range(nrows)
    )
    scale = max(chrom_sizes.values()) + 2 * TELOMERE_LENGTH
    # Biopython applies its 0.95 to an absolute page coordinate, so beyond ~20
    # rows the space it computes goes negative; draw everything as one point
    # rather than dividing by it.
    return scale / points if points > 1 else float(scale)


def create_diagram(
    cnarr: CopyNumArray,
    segarr: CopyNumArray,
    threshold: float,
    min_probes: int,
    outfname: str,
    show_range: str | None = None,
    title: str | None = None,
    show_labels: bool = True,
    *,
    threshold_low: float | None = None,
    threshold_high: float | None = None,
    gene_names: list[str] | None = None,
    squash_genes: bool = True,
) -> str:
    """Create the diagram.

    Segments (``segarr``) are drawn on the left half of each chromosome, one
    feature apiece. The bin-level array (``cnarr``) is drawn on the right, but
    not one feature per bin: Biopython strokes every feature a point wide
    whatever its span, and a chromosome body is only a couple of hundred
    points tall, so on a comprehensive panel or an exome the bins overdraw
    each other and the page shows whichever bin was painted last. The bin half
    is therefore aggregated by ``_aggregate_bins``, to genes where the .cnr
    names them and to the page's own resolution everywhere else. Passing
    ``squash_genes=False`` draws one feature per bin instead, which suits a
    small genome or an assembly whose contigs are short enough for individual
    bins to be legible.

    Gene labels are drawn for segments whose log2 ratio passes a threshold.
    By default the symmetric ``threshold`` is used (label when
    ``abs(log2) >= threshold``). Passing ``threshold_high`` and/or
    ``threshold_low`` switches to directional thresholds: label gains at or
    above ``threshold_high`` and losses at or below ``threshold_low``; a side
    left as ``None`` disables labeling in that direction. ``gene_names``, if
    given, restricts labels to those genes (among the ones passing the
    threshold), dropping co-binned neighbors.
    """
    if cnarr and segarr:
        do_both = True  # Draw segments on left, probes on right.
        cnarr_is_seg = False  # Are probes actually the segmented values?
    else:
        if cnarr:
            cnarr_is_seg = False
        elif segarr:
            cnarr = segarr
            cnarr_is_seg = True
        else:
            raise ValueError("Must provide argument cnarr or segarr, or both. ")
        do_both = False

    if show_range:  # type: ignore[unreachable]
        chrom, start, end = unpack_range(show_range)
        if not (start is None and end is None):
            raise ValueError(
                "Must provide chromosome only "
                "(genomic-range not allowed for 'diagram')."
            )
        if cnarr:
            cnarr = cnarr.in_range(chrom=chrom, start=None, end=None)
        if segarr:
            segarr = segarr.in_range(chrom=chrom, start=None, end=None)
    # The symmetric -t/--threshold maps to the pair (-threshold, +threshold),
    # for which "log2 >= high or log2 <= low" is exactly "abs(log2) >= threshold"
    # -- so the default output is unchanged. Directional flags override one or
    # both sides; a None side disables labeling in that direction.
    low: float | None
    high: float | None
    if threshold_low is None and threshold_high is None:
        low, high = -threshold, threshold
    else:
        low, high = threshold_low, threshold_high
    # Selecting labels means grouping every bin by gene, seconds of work on a
    # WGS array -- skip it when nothing will be labeled. An empty set then
    # leaves every row unlabeled below.
    gene_labels = (
        _get_gene_labels(cnarr, segarr, cnarr_is_seg, low, high, min_probes)
        if show_labels
        else set()
    )
    requested = set(gene_names) if gene_names else None

    # A gene that spans several segments, or is co-binned with its neighbors,
    # reaches this loop once per row it appears in; label it only at the first.
    # TODO - label the gene's most significant row instead of its first
    seen_genes: set[str] = set()

    # The segment half (below) is drawn separately, one feature per segment
    features = collections.defaultdict(list)
    strand = 1 if do_both else None  # Draw on the chr. right half or full width
    chrom_sizes = plots.chromosome_sizes(cnarr)
    drawn = (
        cnarr
        if cnarr_is_seg or not squash_genes
        else _aggregate_bins(cnarr, _bp_per_point(chrom_sizes))
    )
    for row in drawn:
        span = _feature_span(row.start, row.end, chrom_sizes[row.chromosome])
        if span is not None:
            lo, hi = span
            feat_name = _gene_feature_label(
                row.gene, gene_labels, requested, seen_genes
            )
            features[row.chromosome].append(
                (
                    lo,
                    hi,
                    strand,
                    feat_name,
                    colors.Color(*plots.cvg2rgb(row.log2, not cnarr_is_seg)),
                )
            )
    if do_both:
        # Draw segments in the left half of each chromosome (strand -1)
        for chrom, segrows in segarr.by_chromosome():
            for srow in segrows:
                span = _feature_span(srow.start, srow.end, chrom_sizes[chrom])
                if span is not None:
                    lo, hi = span
                    features[chrom].append(
                        (
                            lo,
                            hi,
                            -1,
                            None,
                            colors.Color(*plots.cvg2rgb(srow.log2, False)),
                        )
                    )

    # Generate the diagram PDF
    if not outfname:
        outfname = cnarr.sample_id + "-diagram.pdf"
    drawing = build_chrom_diagram(features, chrom_sizes, cnarr.sample_id, title)
    cvs = canvas.Canvas(outfname, pagesize=PAGE_SIZE)
    renderPDF.draw(drawing, cvs, 0, 0)
    cvs.showPage()
    cvs.save()
    return outfname


class _Feature(NamedTuple):
    """One drawn row of the bin half, in the shape the feature loop reads."""

    chromosome: str
    start: int
    end: int
    log2: float
    gene: str


def _aggregate_bins(cnarr: CopyNumArray, resolution: float) -> Iterator[_Feature]:
    """Reduce bins to the features a page can actually distinguish.

    Each targeted gene becomes one feature spanning its bins, carrying the
    log2 ratio ``genemetrics`` reports for it, so the ideogram and the table
    agree. Every stretch the .cnr does not name -- off-target bins, or all of
    them in an unannotated file -- is cut on a grid of ``resolution`` base
    pairs, one drawn point, since nothing finer survives the page. The two
    rules apply per stretch, not per file: an annotated whole-genome .cnr gets
    gene features at its genes and gridded features between them.

    Whatever then still shares a drawn point is combined, because the page
    cannot show it separately and drawing it anyway means the last feature
    painted decides the color. On a panel this changes nothing -- genes are
    far apart in page terms -- but an annotated whole-genome .cnr puts
    thousands of gene features on a couple of hundred points, and there the
    gene-by-gene reading is already lost.
    """
    for chrom, subarr in cnarr.by_chromosome():
        cells: dict[int, list[tuple[str, CopyNumArray]]] = {}
        for name, rows in subarr.by_gene():
            if not len(rows):
                continue
            if name in _UNNAMED:
                for cell, cell_rows in _split_on_grid(rows, resolution):
                    cells.setdefault(cell, []).append(("", cell_rows))
            else:
                cell = int(rows.start.min() // resolution)
                cells.setdefault(cell, []).append((name, rows))
        for _cell, members in sorted(cells.items()):
            yield _combine(chrom, members)


def _split_on_grid(
    rows: CopyNumArray, resolution: float
) -> Iterator[tuple[int, CopyNumArray]]:
    """Cut a run of unnamed bins at the boundaries between drawn points."""
    # Bins arrive in genomic order, so equal cells are consecutive and the
    # boundaries between them are where the cell index steps up.
    cell_of_bin = (rows.start.to_numpy() // resolution).astype(int)
    edges = np.flatnonzero(np.diff(cell_of_bin)) + 1
    for lo, hi in zip(
        np.concatenate(([0], edges)),
        np.concatenate((edges, [len(cell_of_bin)])),
        strict=True,
    ):
        yield int(cell_of_bin[lo]), rows.as_dataframe(rows.data.iloc[lo:hi])


def _combine(chrom: str, members: list[tuple[str, CopyNumArray]]) -> _Feature:
    """One drawn feature from everything that landed on one point.

    A gene with a point to itself keeps the weighted mean ``genemetrics``
    publishes for it, so the plot and the table agree wherever the page has
    room to say so.

    Where several things share a point, targeted genes decide the color and
    the off-target bins around them are left out of it: on a capture panel a
    gene sits inside a much larger off-target bin, and averaging the two
    reports the background instead of the gene -- measured on
    test/formats/p2-20_1.cnr, METTL8 reads -1.17 in the table and -0.11 if its
    neighbors are included.

    The summary is then the biweight location of the bins in question, not
    their mean. A mean cannot be used across many bins here because ``fix``
    floors zero-coverage bins near -20, which is a censoring marker rather
    than a measurement: on test/formats/wgs-chr17.cnr, where 21,361 of 165,626
    bins sit at that floor, a weighted mean drives 288 of the 371 combined
    cells past the color scale's -1.33 cutoff, against 12 for the biweight
    location. Which estimator each summarizing site in CNVkit should use is a
    wider question than this plot.
    """
    names = list(dict.fromkeys(name for name, _rows in members if name))
    start = min(rows.start.min() for _name, rows in members)
    end = max(rows.end.max() for _name, rows in members)
    if len(members) == 1 and names:
        return _Feature(chrom, start, end, segment_mean(members[0][1]), names[0])
    speaking = [rows for name, rows in members if name] or [
        rows for _name, rows in members
    ]
    log2s = np.concatenate([rows.log2.to_numpy() for rows in speaking])
    return _Feature(
        chrom, start, end, descriptives.biweight_location(log2s), ",".join(names)
    )


def _split_gene_names(gene: Any) -> list[str]:
    """Individual gene names in a row's ``gene`` field, in order, deduplicated.

    A bin or segment names every gene it covers, comma-joined; placeholders
    ("-", "Antitarget", ...) name none. The field may also be missing: the BED
    and interval-list readers default it to "-", but ``read_tab`` passes a
    blank column through as NaN, and an array built in memory need not set it
    at all.
    """
    if not isinstance(gene, str):
        return []
    names = (name.strip() for name in gene.split(","))
    return list(dict.fromkeys(name for name in names if name and name not in _UNNAMED))


def _get_gene_labels(
    cnarr: CopyNumArray,
    segarr: CopyNumArray,
    cnarr_is_seg: bool,
    threshold_low: float | None,
    threshold_high: float | None,
    min_probes: int,
) -> set[str]:
    """Gene names whose copy ratio passes a directional threshold.

    A gene qualifies when its log2 value is at or above ``threshold_high``
    (a gain) or at or below ``threshold_low`` (a loss). Either bound may be
    ``None`` to disable labeling in that direction.

    Segments name every gene they cover, so the qualifying rows are split into
    individual names: the drawn rows are bins and segments, not genes, and a
    whole comma-joined string would match only the rows spelling it the same
    way.
    """
    if cnarr_is_seg:
        # Only segments (.cns): build the directional mask directly.
        mask = cnarr.data["log2"].map(
            lambda v: _passes(v, threshold_low, threshold_high)
        )
        rows = cnarr.data[mask].itertuples(index=False)
        probes_attr = "probes"
    elif segarr:
        # Both segments and bin-level ratios. gene_metrics_by_segment filters on
        # abs(log2) >= t, so feed it the least-restrictive magnitude covering
        # both bounds and refine the direction afterward.
        rows = (
            row
            for row in reports.gene_metrics_by_segment(
                cnarr, segarr, _min_magnitude(threshold_low, threshold_high)
            )
            if _passes(row.log2, threshold_low, threshold_high)
        )
        probes_attr = "segment_probes"
    else:
        # Only bin-level ratios (.cnr)
        rows = (
            row
            for row in reports.gene_metrics_by_gene(
                cnarr, _min_magnitude(threshold_low, threshold_high)
            )
            if _passes(row.log2, threshold_low, threshold_high)
        )
        probes_attr = "probes"
    return {
        name
        for row in rows
        if getattr(row, probes_attr) >= min_probes
        for name in _split_gene_names(row.gene)
    }


def _passes(
    log2: float, threshold_low: float | None, threshold_high: float | None
) -> bool:
    """True if ``log2`` is a gain (>= high) or loss (<= low); None disables a side."""
    return (threshold_high is not None and log2 >= threshold_high) or (
        threshold_low is not None and log2 <= threshold_low
    )


def _min_magnitude(threshold_low: float | None, threshold_high: float | None) -> float:
    """Least-restrictive symmetric cutoff covering both directional bounds."""
    mags = [abs(t) for t in (threshold_low, threshold_high) if t is not None]
    return min(mags) if mags else 0.0


def _feature_span(start: int, end: int, chrom_size: int) -> tuple[int, int] | None:
    """0-based half-open span ``[lo, hi)`` for a feature, or None if out of range.

    Normalizes reverse-oriented intervals to span ``[min, max]``.
    ``skgenome.tabio.read`` keeps such a row out of every region and
    copy-number array it returns -- reversing it out of a BED, refusing it in
    CNVkit's own formats -- but the VCF readers are exempt and an array built
    in memory is not covered at all. Biopython's chromosome renderer asserts
    ``0 <= start <= end <= length``, so an un-normalized reversed interval
    would crash the diagram.

    A feature reaching the start of the chromosome is clamped there rather
    than dropped. The renderer's lower bound is what the subtraction can
    violate, and it can only be violated by one base; dropping instead cost
    the whole first feature, which since the bin half is aggregated is not one
    bin but every bin up to the first drawn point -- 474 of them, three named
    genes among them, on test/formats/wgs-chr17.cnr.
    """
    lo = max(0, min(start, end) - 1)
    hi = max(start, end)
    # `lo` is clamped, `hi` is not, so an interval lying entirely below zero
    # would come back inverted; Biopython asserts against exactly that.
    if lo <= hi <= chrom_size:
        return lo, hi
    return None


def _gene_feature_label(
    gene: Any,
    gene_labels: set[str],
    requested: set[str] | None,
    seen: set[str],
) -> str | None:
    """Build the on-plot label for a drawn row, or None if it earns none.

    Keep the row's gene names that passed the threshold, that ``--gene``
    requested (when given), and that no earlier row has already been labeled
    with; the kept names join into one label and are recorded in ``seen``.
    """
    names = [
        name
        for name in _split_gene_names(gene)
        if name in gene_labels
        and name not in seen
        and (requested is None or name in requested)
    ]
    if not names:
        return None
    seen.update(names)
    # TODO - line-wrap a multi-gene label (reportlab won't do \n)
    return ", ".join(names)


def build_chrom_diagram(
    features: collections.defaultdict[
        str, list[tuple[int, int, int | None, str | None, colors.Color]]
    ],
    chr_sizes: collections.OrderedDict,
    sample_id: str,
    title: str | None = None,
) -> Drawing:
    """Create a PDF of color-coded features on chromosomes."""
    max_chr_len = max(chr_sizes.values())

    chr_diagram = BC.Organism()
    chr_diagram.page_size = PAGE_SIZE
    chr_diagram.title_size = 18

    for chrom, length in list(chr_sizes.items()):
        chrom_features = features.get(chrom)
        if not chrom_features:
            continue
        body = BC.AnnotatedChromosomeSegment(length, chrom_features)
        body.label_size = 4
        body.scale = length
        body.chr_percent = CHROM_FATNESS

        # Create opening and closing telomeres
        tel_start = BC.TelomereSegment()
        tel_start.scale = TELOMERE_LENGTH  # type: ignore[assignment]
        tel_start.chr_percent = CHROM_FATNESS
        tel_end = BC.TelomereSegment(inverted=True)
        tel_end.scale = TELOMERE_LENGTH  # type: ignore[assignment]
        tel_end.chr_percent = CHROM_FATNESS

        # Assemble the chromosome diagram in order
        cur_chromosome = BC.Chromosome(chrom)
        cur_chromosome.title_size = 14
        # Set the scale to the MAXIMUM length plus the two telomeres in bp,
        # want the same scale used on all chromosomes so they can be
        # compared to each other
        cur_chromosome.scale_num = max_chr_len + 2 * TELOMERE_LENGTH
        cur_chromosome.add(tel_start)
        cur_chromosome.add(body)
        cur_chromosome.add(tel_end)
        chr_diagram.add(cur_chromosome)

    if not title:
        title = "Sample " + sample_id
    return bc_organism_draw(chr_diagram, title)


def bc_organism_draw(org: BC.Organism, title: str, wrap: int = CHROM_WRAP) -> Drawing:
    """Modified copy of Bio.Graphics.BasicChromosome.Organism.draw.

    Instead of stacking chromosomes horizontally (along the x-axis), stack rows
    vertically, then proceed with the chromosomes within each row.

    Parameters
    ----------
    org :
        The chromosome diagram object being modified.
    title : str
        The output title of the produced document.
    wrap : int
        Maximum number of chromosomes per row; the remainder will be wrapped to
        the next row(s).
    """
    margin_top = MARGIN_TOP
    margin_bottom = MARGIN_BOTTOM
    margin_side = MARGIN_SIDE

    width, height = org.page_size
    cur_drawing = BC.Drawing(width, height)

    # Draw the title text
    title_string = BC.String(width / 2, height - margin_top + 0.5 * inch, title)
    title_string.fontName = "Helvetica-Bold"
    title_string.fontSize = org.title_size
    title_string.textAnchor = "middle"
    cur_drawing.add(title_string)

    # Layout subcomponents (individual chromosomes), wrapping into rows
    if len(org._sub_components) > 0:
        nrows = math.ceil(len(org._sub_components) / wrap)
        x_pos_change = (width - 2 * margin_side) / wrap
        y_pos_change = (height - margin_top - margin_bottom) / nrows

    cur_x_pos = margin_side
    cur_row = 0
    for i, sub_component in enumerate(org._sub_components):
        if i % wrap == 0 and i != 0:
            cur_row += 1
            cur_x_pos = margin_side
        # Set the page coordinates of the chromosome drawing
        sub_component.start_x_position = cur_x_pos + 0.05 * x_pos_change
        sub_component.end_x_position = cur_x_pos + 0.95 * x_pos_change
        sub_component.start_y_position = height - margin_top - y_pos_change * cur_row
        sub_component.end_y_position = margin_bottom + y_pos_change * (
            nrows - cur_row - 1
        )
        # Render the chromosome drawing
        sub_component.draw(cur_drawing)
        # Update the locations for the next chromosome
        cur_x_pos += x_pos_change

    # Draw a legend
    # (Rect coordinates are: left, bottom, width, height)
    # Bounding box -- near-bottom, center
    cur_drawing.add(
        BC.Rect(
            width / 2 - 0.8 * inch,
            0.5 * inch,
            1.6 * inch,
            0.4 * inch,
            fillColor=colors.white,
        )
    )
    # Red box & label -- in left half of bounding box
    cur_drawing.add(
        BC.Rect(
            width / 2 - 0.7 * inch,
            0.6 * inch,
            0.2 * inch,
            0.2 * inch,
            fillColor=colors.Color(0.8, 0.2, 0.2),
        )
    )
    cur_drawing.add(
        BC.String(
            width / 2 - 0.42 * inch,
            0.65 * inch,
            "Gain",
            fontName="Helvetica",
            fontSize=12,
        )
    )
    # Blue box & label -- in right half of bounding box
    cur_drawing.add(
        BC.Rect(
            width / 2 + 0.07 * inch,
            0.6 * inch,
            0.2 * inch,
            0.2 * inch,
            fillColor=colors.Color(0.2, 0.2, 0.8),
        )
    )
    cur_drawing.add(
        BC.String(
            width / 2 + 0.35 * inch,
            0.65 * inch,
            "Loss",
            fontName="Helvetica",
            fontSize=12,
        )
    )

    # Let the caller take care of writing to the file...
    return cur_drawing


def bc_chromosome_draw_label(
    self: BC.Chromosome, cur_drawing: Drawing, label_name: str
) -> None:
    """Monkeypatch to Bio.Graphics.BasicChromosome.Chromosome._draw_label.

    Draw a label for the chromosome. Mod: above the chromosome, not below.
    """
    # Center on chromosome image
    x_position = 0.5 * (self.start_x_position + self.end_x_position)
    # Place at the bottom of the diagram?
    y_position = self.start_y_position + 0.1 * inch  # was: self.end_y_position
    label_string = BC.String(x_position, y_position, label_name)
    label_string.fontName = "Times-BoldItalic"
    label_string.fontSize = self.title_size
    label_string.textAnchor = "middle"
    cur_drawing.add(label_string)


BC.Chromosome._draw_label = bc_chromosome_draw_label  # type: ignore[method-assign]
