"""Transform bait intervals into targets more suitable for CNVkit."""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

import numpy as np

from skgenome import tabio

from . import antitarget

if TYPE_CHECKING:
    from collections.abc import Iterator

    from pandas.core.series import Series

    from skgenome.gary import GenomicArray


#: Bases added to each side of a genic region before its fine bins are drawn,
#: so that a breakpoint just outside a gene is resolved at the genic bin size.
GENIC_FLANK = 1000

#: GFF feature types that describe a whole sequence rather than a gene. NCBI
#: RefSeq GFF3 files, for example, open each chromosome with a 'region' row.
SEQUENCE_REGION_TYPES = ("region", "chromosome", "contig", "supercontig", "scaffold")


def do_target(
    bait_arr: GenomicArray,
    annotate: str | None = None,
    do_short_names: bool = False,
    do_split: bool = False,
    avg_size: int | float = 200 / 0.75,
    intergenic_bin_scale: float | None = None,
    genic: GenomicArray | None = None,
) -> GenomicArray:
    """Transform bait intervals into targets more suitable for CNVkit.

    Parameters
    ----------
    bait_arr : GenomicArray
        Bait intervals from a BED or interval file.
    annotate : str, optional
        Path to annotation file (BED, GFF, refFlat, etc.) to assign gene names
        to target regions.
    do_short_names : bool, optional
        Reduce multi-accession bait labels to be short and consistent.
        Default is False.
    do_split : bool, optional
        Normalize target sizes: join contiguous intervals and divide them into
        pieces of about `avg_size`. Bin boundaries already present in the bait
        file are not preserved. Default is False.
    avg_size : float, optional
        Average target size when splitting large intervals.
        Default is 200/0.75 (~267 bp).
    intergenic_bin_scale : float, optional
        With `do_split`, split the intervals with `split_genic`: about
        `avg_size` in and near genic regions, and this many times larger
        between them. By default the whole of each interval is split evenly.
    genic : GenomicArray, optional
        The genic regions for `intergenic_bin_scale`. Default: the genes in
        `annotate`.

    Returns
    -------
    GenomicArray
        Processed target intervals ready for CNVkit analysis: zero-width bait
        intervals have been dropped and overlapping ones merged, so the targets
        are disjoint and each genomic position falls in at most one bin.
        Without `do_split`, bait intervals that merely abut are kept apart.

    Raises
    ------
    ValueError
        If `intergenic_bin_scale` is given without `do_split`, or with neither
        `genic` nor `annotate`, or as `split_genic` raises it.
    """
    if intergenic_bin_scale is not None and not do_split:
        raise ValueError("Intergenic bin scaling applies only with do_split.")
    tgt_arr = bait_arr.copy()
    # Drop zero-width regions
    tgt_arr = tgt_arr[tgt_arr.start != tgt_arr.end]
    unmerged_baits = tgt_arr  # kept for re-labeling the split bins below
    # Merge duplicate and overlapping baits, joining their gene labels, so that
    # each position is covered by one target. Left in place, duplicated
    # coordinates reach `fix` as "Duplicated genomic coordinates" (#567), and
    # overlapping baits count the shared reads twice. Bookended baits (the
    # consecutive tiles of a capture kit) are kept separate, so a target BED
    # without overlaps is passed through unchanged. `subdivide` would merge
    # these anyway, but doing it here holds the postcondition on both paths.
    n_baits = len(tgt_arr)
    tgt_arr = tgt_arr.merge(bp=1)
    if len(tgt_arr) < n_baits:
        logging.info(
            "Merged overlapping or duplicated baits: %d intervals -> %d targets",
            n_baits,
            len(tgt_arr),
        )
    annotation = tabio.read_auto(annotate) if annotate else None
    if do_split:
        if intergenic_bin_scale is not None:
            if genic is None:
                if annotation is None:
                    raise ValueError(
                        "Intergenic bin scaling needs genic regions or an annotation."
                    )
                genic = annotation
            logging.info("Splitting targets, more coarsely between genes")
            tgt_arr = split_genic(tgt_arr, genic, avg_size, intergenic_bin_scale)
        else:
            logging.info("Splitting large targets")
            tgt_arr = tgt_arr.subdivide(avg_size, 0)
        if annotation is None:
            # Unless an annotation is about to replace every label: splitting
            # re-merges bookended regions before dividing them, so the joined
            # label of a whole contiguous run would otherwise be stamped on
            # every bin the run is cut into. Take each bin's label from the
            # baits it actually covers instead.
            tgt_arr["gene"] = unmerged_baits.into_ranges(tgt_arr, "gene", "-")
    if annotation is not None:
        logging.info("Applying annotations as target names")
        antitarget.compare_chrom_names(tgt_arr, annotation)
        tgt_arr["gene"] = annotation.into_ranges(tgt_arr, "gene", "-")
        # compare_chrom_names only catches *fully* disjoint chromosome names
        # (e.g. 'chr1' vs '1'). If the names overlap but no coordinates do --
        # most often an annotation file for a different genome build -- every
        # region is left unnamed. Flag that rather than failing silently (#688).
        if len(tgt_arr) and not (tgt_arr["gene"] != "-").any():
            logging.warning(
                "No target regions were assigned a gene name from %s, even "
                "though its chromosome names match the targets. This usually "
                "means the annotation file is for a different genome build -- "
                "check that its assembly and chromosome naming match your "
                "targets.",
                annotate,
            )
    if do_short_names:
        logging.info("Shortening target interval labels")
        tgt_arr["gene"] = list(shorten_labels(tgt_arr["gene"]))
    return tgt_arr  # type: ignore[no-any-return]


def split_genic(
    regions: GenomicArray,
    genic: GenomicArray,
    avg_size: float,
    intergenic_bin_scale: float,
    flank: int = GENIC_FLANK,
) -> GenomicArray:
    """Divide regions into fine bins near genes and coarse bins between them.

    The parts of `regions` within `flank` bases of a genic region are divided
    into bins of about `avg_size`, and the remaining, intergenic parts into
    bins of about ``avg_size * intergenic_bin_scale``. Together the bins tile
    `regions` exactly, as uniform binning with `GenomicArray.subdivide` does,
    and contiguous regions are likewise joined before they are divided.

    Two rules keep bins from becoming much narrower than their neighbors,
    since such bins are noisy enough to produce spurious segments:

    - A genic region narrower than `avg_size`, flanks included, is widened
      to `avg_size`, so that it fills at least one fine bin.
    - An intergenic piece narrower than ``2 * avg_size`` is binned at the
      genic size, joining the genic bins it abuts. Conversely, a genic piece
      that a region boundary clips to less than half of `avg_size` is binned
      with the intergenic bins it abuts.

    Parameters
    ----------
    regions : GenomicArray
        The regions to divide into bins, e.g. the sequencing-accessible
        regions of the genome.
    genic : GenomicArray
        The regions to bin finely, e.g. the transcribed spans of a gene
        annotation, or a chosen subset of genes. GFF rows describing a whole
        sequence (`SEQUENCE_REGION_TYPES`) are ignored.
    avg_size : float
        Average size of the bins in and near genic regions.
    intergenic_bin_scale : float
        Ratio of the average intergenic bin size to `avg_size`.
    flank : int, optional
        Bases added to each side of each genic region.

    Returns
    -------
    GenomicArray
        Sorted bins, with coordinate columns only.

    Raises
    ------
    ValueError
        If `intergenic_bin_scale` is not a finite number greater than 1, if
        `genic` holds no regions apart from whole-sequence GFF rows, or if no
        genic region overlaps `regions`, which usually means the two files use
        different chromosome names or genome builds.
    """
    if not 1 < intergenic_bin_scale < math.inf:
        raise ValueError(
            "Intergenic bin scale must be a finite number greater than 1, "
            f"not {intergenic_bin_scale}."
        )
    if "type" in genic:
        genic = genic[~genic["type"].isin(SEQUENCE_REGION_TYPES)]
    if not len(genic):
        raise ValueError(
            "No genic regions were given, apart from any GFF rows spanning "
            "a whole sequence."
        )
    loci = genic.drop_extra_columns().resize_ranges(flank)
    pad = np.ceil((avg_size - (loci.end - loci.start)).clip(lower=0) / 2)
    loci = loci.add_columns(
        start=(loci.start - pad).clip(lower=0).astype(int),
        end=(loci.end + pad).astype(int),
    ).merge()
    space = regions.drop_extra_columns().merge()
    near_genes = space.intersection(loci, mode="trim")
    if not len(near_genes):
        raise ValueError(
            "No genic region overlaps the regions to be binned; check that "
            "both use the same genome build and chromosome names."
        )
    intergenic = space.subtract(loci)
    is_sliver = (intergenic.end - intergenic.start) < 2 * avg_size
    near_genes = near_genes.concat([near_genes, intergenic[is_sliver]]).merge()
    is_stub = (near_genes.end - near_genes.start) < avg_size / 2
    intergenic = intergenic.concat([intergenic[~is_sliver], near_genes[is_stub]])
    fine = near_genes[~is_stub].subdivide(avg_size)
    coarse = intergenic.subdivide(avg_size * intergenic_bin_scale)
    logging.info(
        "Binned %d bp near genes into %d bins and %d bp between genes into %d bins",
        fine.total_range_size(),
        len(fine),
        coarse.total_range_size(),
        len(coarse),
    )
    if not len(fine):
        logging.warning(
            "No bins have the genic size, because the genic regions, extended by "
            "their flanks and to at least one genic bin, overlap the regions to be "
            "binned only in pieces narrower than half a genic bin. Check that the "
            "genic regions are the intended genes."
        )
    elif not len(coarse):
        logging.warning(
            "No bins have the intergenic size, because the genic regions, extended "
            "by their flanks and to at least one genic bin, leave no intergenic "
            "stretch as wide as two genic bins. Check that the genic regions are "
            "the intended genes."
        )
    return fine.concat([fine, coarse])  # type: ignore[no-any-return]


def shorten_labels(gene_labels: Series) -> Iterator[str]:
    """Reduce multi-accession interval labels to the minimum consistent.

    So: BED or interval_list files have a label for every region. We want this
    to be a short, unique string, like the gene name. But if an interval list is
    instead a series of accessions, including additional accessions for
    sub-regions of the gene, we can extract a single accession that covers the
    maximum number of consecutive regions that share this accession.

    e.g.::

        ...
        mRNA|JX093079,ens|ENST00000342066,mRNA|JX093077,ref|SAMD11,mRNA|AF161376,mRNA|JX093104
        ens|ENST00000483767,mRNA|AF161376,ccds|CCDS3.1,ref|NOC2L
        ...

    becomes::

        ...
        mRNA|AF161376
        mRNA|AF161376
        ...
    """
    longest_name_len = 0
    curr_names: set[str] = set()
    curr_gene_count = 0

    for label in gene_labels:
        next_names = set(label.rstrip().split(","))
        assert len(next_names)
        overlap = curr_names.intersection(next_names)
        if overlap:
            # Continuing the same gene; update shared accessions
            curr_names = filter_names(overlap)
            curr_gene_count += 1
        else:
            # End of the old gene -- emit shared name(s)
            for _i in range(curr_gene_count):
                out_name = shortest_name(curr_names)
                yield out_name
                longest_name_len = max(longest_name_len, len(out_name))

            # Start of a new gene
            curr_gene_count = 1
            curr_names = next_names
    # Final emission
    for _i in range(curr_gene_count):
        out_name = shortest_name(curr_names)
        yield out_name
        longest_name_len = max(longest_name_len, len(out_name))

    logging.info("Longest name length: %d", longest_name_len)


def filter_names(names: set[str], exclude: tuple[str] = ("mRNA",)) -> set[str]:
    """Remove less-meaningful accessions from the given set."""
    if len(names) > 1:
        ok_names = {n for n in names if not any(n.startswith(ex) for ex in exclude)}
        if ok_names:
            return ok_names
    # Names are not filter-worthy; leave them as they are for now
    return names


def shortest_name(names: set[str]) -> str:
    """Return the shortest trimmed name from the given set."""
    name = min(filter_names(names), key=len)
    if len(name) > 2 and "|" in name[1:-1]:
        # Split 'DB|accession' and extract the accession sans-DB
        name = name.split("|")[-1]
    return name
