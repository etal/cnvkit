"""Chromosome diagram drawing functions.

This uses and abuses Biopython's BasicChromosome module. It depends on
ReportLab, too, so we isolate this functionality here so that the rest of
CNVkit will run without it. (And also to keep the codebase tidy.)
"""

from __future__ import annotations

import collections
import math
import warnings
from typing import TYPE_CHECKING, Any

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

from skgenome.rangelabel import unpack_range

from . import params, plots, reports

if TYPE_CHECKING:
    from reportlab.graphics.shapes import Drawing

    from cnvlib.cnary import CopyNumArray

TELOMERE_LENGTH = 6e6  # For illustration only
CHROM_FATNESS = 0.3
PAGE_SIZE = (11.0 * inch, 8.5 * inch)


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
) -> str:
    """Create the diagram.

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
    gene_labels = _get_gene_labels(cnarr, segarr, cnarr_is_seg, low, high, min_probes)
    requested = set(gene_names) if gene_names else None

    # NB: If multiple segments cover the same gene (gene contains breakpoints),
    # all those segments are marked as "hits".  We'll uniquify them.
    # TODO - use different logic to only label the gene's signficant segment(s)
    seen_genes = set()

    # Consolidate genes & coverage values as chromosome features
    features = collections.defaultdict(list)
    strand = 1 if do_both else None  # Draw on the chr. right half or full width
    chrom_sizes = plots.chromosome_sizes(cnarr)
    if not cnarr_is_seg:
        cnarr = cnarr.squash_genes()
    for row in cnarr:
        span = _feature_span(row.start, row.end, chrom_sizes[row.chromosome])
        if span is not None:
            lo, hi = span
            if show_labels and row.gene in gene_labels and row.gene not in seen_genes:
                seen_genes.add(row.gene)
                # TODO - line-wrap multi-gene labels (reportlab won't do \n)
                feat_name = _gene_feature_label(row.gene, requested)
            else:
                feat_name = None
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
    drawing = build_chrom_diagram(features, chrom_sizes, cnarr.sample_id, title)  # type: ignore[arg-type]
    cvs = canvas.Canvas(outfname, pagesize=PAGE_SIZE)
    renderPDF.draw(drawing, cvs, 0, 0)
    cvs.showPage()
    cvs.save()
    return outfname


def _get_gene_labels(
    cnarr: CopyNumArray,
    segarr: CopyNumArray,
    cnarr_is_seg: bool,
    threshold_low: float | None,
    threshold_high: float | None,
    min_probes: int,
) -> list[Any]:
    """Label genes whose copy ratio passes a directional threshold.

    A gene qualifies when its log2 value is at or above ``threshold_high``
    (a gain) or at or below ``threshold_low`` (a loss). Either bound may be
    ``None`` to disable labeling in that direction.
    """
    if cnarr_is_seg:
        # Only segments (.cns): build the directional mask directly.
        mask = cnarr.data["log2"].map(
            lambda v: _passes(v, threshold_low, threshold_high)
        )
        sel = cnarr.data[mask & ~cnarr.data.gene.isin(params.IGNORE_GENE_NAMES)]
        rows = sel.itertuples(index=False)
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
    return [row.gene for row in rows if getattr(row, probes_attr) >= min_probes]


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

    Normalizes reverse-oriented intervals: reverse-direction PCR primers are
    sometimes stored with ``start > end`` in the input BED, so the interval is
    taken to span ``[min, max]`` regardless of column order. Biopython's
    chromosome renderer asserts ``0 <= start <= end <= length``, so an
    un-normalized reversed interval would otherwise crash the diagram.
    """
    lo = min(start, end) - 1
    hi = max(start, end)
    if lo >= 0 and hi <= chrom_size:
        return lo, hi
    return None


def _gene_feature_label(gene: str, requested: set[str] | None) -> str | None:
    """Build the on-plot label for a qualifying segment or gene.

    Without ``requested`` (no --gene), expand the bin's comma-joined gene
    string for display, preserving legacy output. With ``requested``, show only
    the requested genes -- dropping co-binned neighbors the user did not ask for
    and collapsing duplicate names.
    """
    if requested is None:
        return gene.replace(",", ", ")
    names = list(
        dict.fromkeys(n.strip() for n in gene.split(",") if n.strip() in requested)
    )
    return ", ".join(names) if names else None


def build_chrom_diagram(
    features: collections.defaultdict[
        str, list[tuple[int, int, int, None, colors.Color]]
    ],
    chr_sizes: collections.OrderedDict,
    sample_id: str,
    title: None = None,
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
        title = "Sample " + sample_id  # type: ignore[assignment]
    return bc_organism_draw(chr_diagram, title)  # type: ignore[arg-type]


def bc_organism_draw(org: BC.Organism, title: str, wrap: int = 12) -> Drawing:
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
    margin_top = 1.25 * inch
    margin_bottom = 0.1 * inch
    margin_side = 0.5 * inch

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
