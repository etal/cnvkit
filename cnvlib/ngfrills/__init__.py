"""NGS utilities."""
from __future__ import absolute_import, division, print_function

import contextlib
import os
import subprocess
import sys
import tempfile

from Bio._py3k import basestring, map

from .faidx import *
from .regions import *
from .samutil import *
from .shared import *


# __________________________________________________________
# Shell

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
                              " that cannot be created: " + str(exc))
    if os.path.isfile(fname):
        # Add an integer suffix to the existing file name
        cnt = 1
        bak_fname = "%s.%d" % (fname, cnt)
        while os.path.isfile(bak_fname):
            cnt += 1
            bak_fname = "%s.%d" % (fname, cnt)
        os.rename(fname, bak_fname)
        echo("Moved existing file", fname, "->", bak_fname)
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
            echo("Created directory", dirname)
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
        echo("Wrote", outfname)


@contextlib.contextmanager
def temp_write_text(text):
    """Save text to a temporary file.

    NB: This won't work on Windows b/c the file stays open.
    """
    with tempfile.NamedTemporaryFile() as tmp:
        tmp.write(text)
        tmp.flush()
        yield tmp.name
