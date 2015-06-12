"""CNV utilities."""
from __future__ import absolute_import, division, print_function
import sys
import os.path

import numpy
from Bio.File import as_handle

from .ngfrills import echo, safe_write

# __________________________________________________________________________
# I/O helpers

def parse_tsv(infile, keep_header=False):
    """Parse a tabular data table into an iterable of lists.

    Rows are split on tabs.  Header row is optionally included in the output.
    """
    with as_handle(infile) as handle:
        lines = iter(handle)
        header = next(lines)
        if keep_header:
            yield header.rstrip().split('\t')
        for line in lines:
            yield line.rstrip().split('\t')


def write_tsv(outfname, table, colnames=None):
    """Write the CGH file."""
    if not outfname:
        outfname = sys.stdout
    with safe_write(outfname) as handle:
        if colnames:
            header = '\t'.join(colnames) + '\n'
            handle.write(header)
        handle.writelines('\t'.join(map(str, row)) + '\n'
                           for row in table)


# __________________________________________________________________________
# Sorting key functions


def sorter_chrom(label):
    """Create a sorting key from chromosome label.

    Sort by integers first, then letters or strings. The prefix "chr"
    (case-insensitive), if present, is stripped automatically for sorting.

    E.g. chr1 < chr2 < chr10 < chrX < chrY
    """
    chrom = (label[3:] if label.lower().startswith('chr')
             else label)
    return (chr(int(chrom)) if chrom.isdigit()
            else chrom)


def sorter_chrom_at(index):
    """Create a sort key function that gets chromosome label at a list index."""
    return lambda row: sorter_chrom(row[index])


# __________________________________________________________________________
# Guess gender: XX or XY?

# XXX refactor all this

def shift_xx(probes, male_reference=False, chr_x=None):
    """Adjust chrX coverages (divide in half) for apparent female samples."""
    if chr_x is None:
        chr_x = guess_chr_x(probes)
    outprobes = probes.copy()
    is_xx = guess_xx(probes, chr_x=chr_x, male_reference=male_reference)
    if is_xx and male_reference:
        # Female: divide X coverages by 2 (in log2: subtract 1)
        outprobes['log2'][outprobes.chromosome == chr_x] -= 1.0
        # Male: no change
    elif not is_xx and not male_reference:
        # Male: multiply X coverages by 2 (in log2: add 1)
        outprobes['log2'][outprobes.chromosome == chr_x] += 1.0
        # Female: no change
    return outprobes


def guess_chr_x(probes):
    return ('chrX' if probes[0, 'chromosome'].startswith('chr')
            else 'X')


def guess_xx(probes, male_reference=False, chr_x=None, verbose=True):
    """Guess whether a sample is female from chrX relative coverages.

    Recommended cutoff values:
        -0.5 -- raw target data, not yet corrected
        +0.5 -- probe data already corrected on a male profile
    """
    cutoff = 0.5 if male_reference else -0.5
    # ENH - better coverage approach: take Z-scores or rank of +1,0 or 0,-1
    # based on the available probes, then choose which is more probable
    rel_chrx_cvg = get_relative_chrx_cvg(probes, chr_x=chr_x)
    is_xx = (rel_chrx_cvg >= cutoff)
    if verbose:
        echo("Relative log2 coverage of X chromosome:", rel_chrx_cvg,
             "(assuming %s)" % ('male', 'female')[is_xx])
    return is_xx


def get_relative_chrx_cvg(probes, chr_x=None):
    """Get the relative log-coverage of chrX in a sample."""
    if chr_x is None:
        chr_x = guess_chr_x(probes)
    chromosome_x = probes[probes.chromosome == chr_x]
    if not len(chromosome_x):
        echo("*WARNING* No", chr_x, "found in probes; check the input")
        return
    chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
    autosomes = probes[(probes.chromosome != chr_x) &
                       (probes.chromosome != chr_y)]
    auto_cvgs = autosomes['log2']
    x_cvgs = chromosome_x['log2']
    if 'probes' in probes:
        # Weight segments by number of probes to ensure good behavior
        auto_sizes = autosomes['probes']
        x_sizes = chromosome_x['probes']
        # ENH: weighted median
        rel_chrx_cvg = (numpy.average(x_cvgs, weights=x_sizes) -
                        numpy.average(auto_cvgs, weights=auto_sizes))
    else:
        rel_chrx_cvg = numpy.median(x_cvgs) - numpy.median(auto_cvgs)
    return rel_chrx_cvg


def expect_flat_cvg(cnarr, is_male_reference=None, chr_x=None):
    """Get the uninformed expected copy ratios of each bin.

    Create an array of log2 coverages like a "flat" reference.

    This is a neutral copy ratio at each autosome (log2 = 0.0) and sex
    chromosomes based on whether the reference is male (XX or XY).
    """
    if chr_x is None:
        chr_x = guess_chr_x(cnarr)
    if is_male_reference is None:
        is_male_reference = not guess_xx(cnarr, chr_x=chr_x, verbose=False)
    chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
    cvg = numpy.zeros(len(cnarr), dtype=numpy.float_)
    if is_male_reference:
        # Single-copy X, Y
        cvg[(cnarr.chromosome == chr_x) |
            (cnarr.chromosome == chr_y)] = -1.0
    else:
        # Y will be all noise, so replace with 1 "flat" copy
        cvg[cnarr.chromosome == chr_y] = -1.0
    return cvg


# __________________________________________________________________________
# More helpers

def assert_equal(msg, **values):
    """Evaluate and compare two or more values for equality.

    Sugar for a common assertion pattern. Saves re-evaluating (and retyping)
    the same values for comparison and error reporting.

    Example:

    >>> assert_equal("Mismatch", expected=1, saw=len(['xx', 'yy']))
    ...
    ValueError: Mismatch: expected = 1, saw = 2

    """
    ok = True
    key1, val1 = values.popitem()
    msg += ": %s = %r" % (key1, val1)
    for okey, oval in values.iteritems():
        msg += ", %s = %r" % (okey, oval)
        if oval != val1:
            ok = False
    if not ok:
        raise ValueError(msg)


def check_unique(items, title):
    """Ensure all items in an iterable are identical; return that one item."""
    its = set(items)
    assert len(its) == 1, ("Inconsistent %s keys: %s"
                           % (title, ' '.join(map(str, sorted(its)))))
    return its.pop()


def fbase(fname):
    """Strip directory and all extensions from a filename."""
    return os.path.basename(fname).split('.', 1)[0]


def rbase(fname):
    """Strip directory and final extension from a filename."""
    return os.path.basename(fname).rsplit('.', 1)[0]
