"""Chromosome diagram drawing functions.

This uses and abuses Biopython's BasicChromosome module. It depends on
ReportLab, too, so we isolate this functionality here so that the rest of CNVkit
will run without it. (And also to keep the codebase tidy.)
"""
from __future__ import absolute_import, division
import collections
import math
import warnings

from Bio.Graphics import BasicChromosome as BC
from reportlab.graphics import renderPDF
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas

from . import params, plots, reports

# Silence Biopython's whinging
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

TELOMERE_LENGTH = 6e6   # For illustration only
CHROM_FATNESS = 0.3
PAGE_SIZE = (11.0*inch, 8.5*inch)


def create_diagram(cnarr, segarr, threshold, min_probes, outfname, title=None):
    """Create the diagram."""
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
    gene_labels = _get_gene_labels(cnarr, segarr, cnarr_is_seg, threshold,
                                   min_probes)

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
        if row.start - 1 >= 0 and row.end <= chrom_sizes[row.chromosome]:  # Sanity check
            if row.gene in gene_labels and row.gene not in seen_genes:
                seen_genes.add(row.gene)
                feat_name = row.gene
                if "," in feat_name:
                    # TODO - line-wrap multi-gene labels (reportlab won't do \n)
                    feat_name = feat_name.replace(",", ", ")
            else:
                feat_name = None
            features[row.chromosome].append(
                (row.start - 1, row.end, strand, feat_name,
                 colors.Color(*plots.cvg2rgb(row.log2, not cnarr_is_seg))))
    if do_both:
        # Draw segments in the left half of each chromosome (strand -1)
        for chrom, segrows in segarr.by_chromosome():
            for srow in segrows:
                if srow.start - 1 >= 0 and srow.end <= chrom_sizes[chrom]:  # Sanity check
                    features[chrom].append(
                        (srow.start - 1, srow.end, -1, None,
                         colors.Color(*plots.cvg2rgb(srow.log2, False))))

    # Generate the diagram PDF
    if not outfname:
        outfname = cnarr.sample_id + '-diagram.pdf'
    drawing = build_chrom_diagram(features, chrom_sizes, cnarr.sample_id, title)
    cvs = canvas.Canvas(outfname, pagesize=PAGE_SIZE)
    renderPDF.draw(drawing, cvs, 0, 0)
    cvs.showPage()
    cvs.save()
    return outfname


def _get_gene_labels(cnarr, segarr, cnarr_is_seg, threshold, min_probes):
    """Label genes where copy ratio value exceeds threshold."""
    if cnarr_is_seg:
        # Only segments (.cns)
        sel = cnarr.data[(cnarr.data.log2.abs() >= threshold) &
                          ~cnarr.data.gene.isin(params.IGNORE_GENE_NAMES)]
        rows = sel.itertuples(index=False)
        probes_attr = 'probes'
    elif segarr:
        # Both segments and bin-level ratios
        rows = reports.gene_metrics_by_segment(cnarr, segarr, threshold)
        probes_attr = 'segment_probes'
    else:
        # Only bin-level ratios (.cnr)
        rows = reports.gene_metrics_by_gene(cnarr, threshold)
        probes_attr = 'n_bins'
    return [row.gene for row in rows if getattr(row, probes_attr) >= min_probes]


def build_chrom_diagram(features, chr_sizes, sample_id, title=None):
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
        tel_start.scale = TELOMERE_LENGTH
        tel_start.chr_percent = CHROM_FATNESS
        tel_end = BC.TelomereSegment(inverted=True)
        tel_end.scale = TELOMERE_LENGTH
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


def bc_organism_draw(org, title, wrap=12):
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
    margin_top = 1.25*inch
    margin_bottom = 0.1*inch
    margin_side = 0.5*inch

    width, height = org.page_size
    cur_drawing = BC.Drawing(width, height)

    # Draw the title text
    title_string = BC.String(width / 2, height - margin_top + .5*inch, title)
    title_string.fontName = 'Helvetica-Bold'
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
        sub_component.start_y_position = (height - margin_top
                                          - y_pos_change * cur_row)
        sub_component.end_y_position = (margin_bottom
                                        + y_pos_change * (nrows - cur_row - 1))
        # Render the chromosome drawing
        sub_component.draw(cur_drawing)
        # Update the locations for the next chromosome
        cur_x_pos += x_pos_change

    # Draw a legend
    # (Rect coordinates are: left, bottom, width, height)
    # Bounding box -- near-bottom, center
    cur_drawing.add(BC.Rect(width/2 - .8*inch, .5*inch, 1.6*inch, .4*inch,
                            fillColor=colors.white))
    # Red box & label -- in left half of bounding box
    cur_drawing.add(BC.Rect(width/2 - .7*inch, .6*inch, .2*inch, .2*inch,
                            fillColor=colors.Color(.8, .2, .2)))
    cur_drawing.add(BC.String(width/2 - .42*inch, .65*inch, "Gain",
                              fontName='Helvetica', fontSize=12))
    # Blue box & label -- in right half of bounding box
    cur_drawing.add(BC.Rect(width/2 + .07*inch, .6*inch, .2*inch, .2*inch,
                            fillColor=colors.Color(.2, .2, .8)))
    cur_drawing.add(BC.String(width/2 + .35*inch, .65*inch, "Loss",
                              fontName='Helvetica', fontSize=12))

    # Let the caller take care of writing to the file...
    return cur_drawing


def bc_chromosome_draw_label(self, cur_drawing, label_name):
    """Monkeypatch to Bio.Graphics.BasicChromosome.Chromosome._draw_label.

    Draw a label for the chromosome. Mod: above the chromosome, not below.
    """
    # Center on chromosome image
    x_position = 0.5 * (self.start_x_position + self.end_x_position)
    # Place at the bottom of the diagram?
    y_position = self.start_y_position + 0.1*inch  # was: self.end_y_position
    label_string = BC.String(x_position, y_position, label_name)
    label_string.fontName = 'Times-BoldItalic'
    label_string.fontSize = self.title_size
    label_string.textAnchor = 'middle'
    cur_drawing.add(label_string)

BC.Chromosome._draw_label = bc_chromosome_draw_label
