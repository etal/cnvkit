"""Supporting functions for the 'antitarget' command."""
from __future__ import absolute_import, division, print_function
from builtins import map

import logging
import re

from skgenome import GenomicArray as GA

from .params import INSERT_SIZE, MIN_REF_COVERAGE, ANTITARGET_NAME


def do_antitarget(targets, access=None, avg_bin_size=150000,
                  min_bin_size=None):
    """Derive off-targt ("antitarget") bins from target regions."""
    if not min_bin_size:
        min_bin_size = 2 * int(avg_bin_size * (2 ** MIN_REF_COVERAGE))
    return get_antitargets(targets, access, avg_bin_size, min_bin_size)


def get_antitargets(targets, accessible, avg_bin_size, min_bin_size):
    """Generate antitarget intervals between/around target intervals.

    Procedure:

    - Invert target intervals
    - Subtract the inverted targets from accessible regions
    - For each of the resulting regions:

        - Shrink by a fixed margin on each end
        - If it's smaller than min_bin_size, skip
        - Divide into equal-size (region_size/avg_bin_size) portions
        - Emit the (chrom, start, end) coords of each portion
    """
    if accessible:
        # Chromosomes' accessible sequence regions are given -- use them
        accessible = drop_noncanonical_contigs(accessible, targets)
    else:
        # Chromosome accessible sequence regions not known -- use heuristics
        # (chromosome length is endpoint of last probe; skip initial
        # <magic number> of bases that are probably telomeric)
        TELOMERE_SIZE = 150000
        accessible = guess_chromosome_regions(targets, TELOMERE_SIZE)
    pad_size = 2 * INSERT_SIZE
    bg_arr = (accessible.resize_ranges(-pad_size)
              .subtract(targets.resize_ranges(pad_size))
              .subdivide(avg_bin_size, min_bin_size))
    bg_arr['gene'] = ANTITARGET_NAME
    return bg_arr


def drop_noncanonical_contigs(accessible, targets, verbose=True):
    """Drop contigs that are not targeted or canonical chromosomes.

    Antitargets will be binned over chromosomes that:

    - Appear in the sequencing-accessible regions of the reference genome
      sequence, and
    - Contain at least one targeted region, or
    - Are named like a canonical chromosome (1-22,X,Y for human)

    This allows antitarget binning to pick up canonical chromosomes that do not
    contain any targets, as well as non-canonical or oddly named chromosomes
    that were targeted.
    """
    # TODO - generalize: (garr, by="name", verbose=True):

    access_chroms, target_chroms = compare_chrom_names(accessible, targets)
    # Filter out untargeted alternative contigs and mitochondria
    untgt_chroms = access_chroms - target_chroms
    # Autosomes typically have numeric names, allosomes are X and Y
    if any(is_canonical_contig_name(c) for c in target_chroms):
        chroms_to_skip = [c for c in untgt_chroms
                            if not is_canonical_contig_name(c)]
    else:
        # Alternative contigs have longer names -- skip them
        max_tgt_chr_name_len = max(map(len, target_chroms))
        chroms_to_skip = [c for c in untgt_chroms
                            if len(c) > max_tgt_chr_name_len]
    if chroms_to_skip:
        logging.info("Skipping untargeted chromosomes %s",
                     ' '.join(sorted(chroms_to_skip)))
        skip_idx = accessible.chromosome.isin(chroms_to_skip)
        accessible = accessible[~skip_idx]
    return accessible


def compare_chrom_names(a_regions, b_regions):
    a_chroms = set(a_regions.chromosome.unique())
    b_chroms = set(b_regions.chromosome.unique())
    if a_chroms and a_chroms.isdisjoint(b_chroms):
        msg = "Chromosome names do not match between files"
        a_fname = a_regions.meta.get('filename')
        b_fname = b_regions.meta.get('filename')
        if a_fname and b_fname:
            msg += " {} and {}".format(a_fname, b_fname)
        msg += ": {} vs. {}".format(', '.join(map(repr, sorted(a_chroms)[:3])),
                                    ', '.join(map(repr, sorted(b_chroms)[:3])))
        raise ValueError(msg)
    return a_chroms, b_chroms


def guess_chromosome_regions(targets, telomere_size):
    """Determine (minimum) chromosome lengths from target coordinates."""
    endpoints = [subarr.end.iat[-1] for _c, subarr in targets.by_chromosome()]
    whole_chroms = GA.from_columns({
        'chromosome': targets.chromosome.drop_duplicates(),
        'start': telomere_size,
        'end': endpoints})
    return whole_chroms


# TODO - move to skgenome.chromsort

# CNVkit's original inclusion regex
re_canonical = re.compile(r"(chr)?(\d+|[XYxy])$")
# goleft indexcov's exclusion regex
# TODO drop chrM, MT
re_noncanonical = re.compile(r"^chrEBV$|^NC|_random$|Un_|^HLA\-|_alt$|hap\d$")

def is_canonical_contig_name(name):
    return bool(re_canonical.match(name))
    #  return not re_noncanonical.search(name))


def _drop_short_contigs(garr):
    """Drop contigs that are much shorter than the others.

    Cutoff is where a contig is less than half the size of the next-shortest
    contig.
    """
    from .plots import chromosome_sizes
    from skgenome.chromsort import detect_big_chroms

    chrom_sizes = chromosome_sizes(garr)
    n_big, thresh = detect_big_chroms(chromosome_sizes.values())
    chrom_names_to_keep = {c for c, s in chrom_sizes.items()
                           if s >= thresh}
    assert len(chrom_names_to_keep) == n_big
    return garr[garr.chromosome.isin(chrom_names_to_keep)]
