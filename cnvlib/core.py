"""CNV utilities."""
from __future__ import absolute_import, division, print_function
from builtins import map
from past.builtins import basestring

import contextlib
import logging
import os
import subprocess
import sys
import tempfile
from itertools import takewhile


# __________________________________________________________________________
# I/O helpers

def call_quiet(*args):
    """Safely run a command and get stdout; print stderr if there's an error.

    Like subprocess.check_output, but silent in the normal case where the
    command logs unimportant stuff to stderr. If there is an error, then the
    full error message(s) is shown in the exception message.
    """
    # args = map(str, args)
    if not len(args):
        raise ValueError("Must supply at least one argument (the command name)")
    try:
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except OSError as exc:
        raise RuntimeError("Could not find the executable %r" % args[0]
                           + " -- is it installed correctly?"
                           + "\n(Original error: %s)" % exc)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError("Subprocess command failed:\n$ %s\n\n%s"
                           % (' '.join(args), err))
    return out


def ensure_path(fname):
    """Create dirs and move an existing file to avoid overwriting, if necessary.

    If a file already exists at the given path, it is renamed with an integer
    suffix to clear the way.
    """
    if '/' in os.path.normpath(fname):
        # Ensure the output directory exists
        dname = os.path.dirname(os.path.abspath(fname))
        if dname and not os.path.isdir(dname):
            try:
                os.makedirs(dname)
            except OSError as exc:
                raise OSError("Output path " + fname +
                              " contains a directory " + dname +
                              " that cannot be created: %s" % exc)
    if os.path.isfile(fname):
        # Add an integer suffix to the existing file name
        cnt = 1
        bak_fname = "%s.%d" % (fname, cnt)
        while os.path.isfile(bak_fname):
            cnt += 1
            bak_fname = "%s.%d" % (fname, cnt)
        os.rename(fname, bak_fname)
        logging.info("Moved existing file %s -> %s", fname, bak_fname)
    return True


@contextlib.contextmanager
def safe_write(outfile, verbose=True):
    """Write to a filename or file-like object with error handling.

    If given a file name, open it. If the path includes directories that don't
    exist yet, create them.  If given a file-like object, just pass it through.
    """
    if isinstance(outfile, basestring):
        dirname = os.path.dirname(outfile)
        if dirname and not os.path.isdir(dirname):
            os.mkdir(dirname)
            logging.info("Created directory %s", dirname)
        with open(outfile, 'w') as handle:
            yield handle
    else:
        yield outfile

    # Log the output path, if possible
    if verbose:
        if isinstance(outfile, basestring):
            outfname = outfile
        elif hasattr(outfile, 'name') and outfile not in (sys.stdout,
                                                          sys.stderr):
            outfname = outfile.name
        else:
            # Probably stdout or stderr -- don't ruin the pipeline
            return
        logging.info("Wrote %s", outfname)


@contextlib.contextmanager
def temp_write_text(text, mode="w+b"):
    """Save text to a temporary file.

    NB: This won't work on Windows b/c the file stays open.
    """
    with tempfile.NamedTemporaryFile(mode=mode) as tmp:
        tmp.write(text)
        tmp.flush()
        yield tmp.name


def write_tsv(outfname, rows, colnames=None):
    """Write rows, with optional column header, to tabular file."""
    with safe_write(outfname or sys.stdout) as handle:
        if colnames:
            header = '\t'.join(colnames) + '\n'
            handle.write(header)
        handle.writelines('\t'.join(map(str, row)) + '\n'
                           for row in rows)


def write_text(outfname, text, *more_texts):
    """Write one or more strings (blocks of text) to a file."""
    with safe_write(outfname or sys.stdout) as handle:
        handle.write(text)
        if more_texts:
            for mtext in more_texts:
                handle.write(mtext)


def write_dataframe(outfname, dframe, header=True):
    """Write a pandas.DataFrame to a tabular file."""
    with safe_write(outfname or sys.stdout) as handle:
        dframe.to_csv(handle, header=header,
                      index=False, sep='\t', float_format='%.6g')


# __________________________________________________________________________
# Sorting key functions

def sorter_chrom(label):
    """Create a sorting key from chromosome label.

    Sort by integers first, then letters or strings. The prefix "chr"
    (case-insensitive), if present, is stripped automatically for sorting.

    E.g. chr1 < chr2 < chr10 < chrX < chrY < chrM
    """
    # Strip "chr" prefix
    chrom = (label[3:] if label.lower().startswith('chr')
             else label)
    if chrom in ('X', 'Y'):
        key = (1000, chrom)
    else:
        # Separate numeric and special chromosomes
        nums = ''.join(takewhile(str.isdigit, chrom))
        chars = chrom[len(nums):]
        nums = int(nums) if nums else 0
        if not chars:
            key = (nums, '')
        elif len(chars) == 1:
            key = (2000 + nums, chars)
        else:
            key = (3000 + nums, chars)
    return key


def sorter_chrom_at(index):
    """Create a sort key function that gets chromosome label at a list index."""
    return lambda row: sorter_chrom(row[index])


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
    for okey, oval in values.items():
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
