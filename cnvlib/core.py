"""CNV utilities."""

import contextlib
import logging
import os
import subprocess
import tempfile


# __________________________________________________________________________
# I/O helpers


def call_quiet(*args):
    """Safely run a command and get stdout; print stderr if there's an error.

    Like subprocess.check_output, but silent in the normal case where the
    command logs unimportant stuff to stderr. If there is an error, then the
    full error message(s) is shown in the exception message.
    """
    if not len(args):
        raise ValueError("Must supply at least one argument (the command name)")
    try:
        proc = subprocess.run(args, check=True, capture_output=True)
    except OSError as exc:
        raise RuntimeError(
            f"Could not find the executable {args[0]!r} -- is it installed correctly?"
            f"\n(Original error: {exc})"
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Subprocess command failed:\n$ {' '.join(args)}\n\n{exc}"
        ) from exc
    return proc.stdout


def ensure_path(fname) -> bool:
    """Create dirs and move an existing file to avoid overwriting, if necessary.

    If a file already exists at the given path, it is renamed with an integer
    suffix to clear the way.
    """
    if "/" in os.path.normpath(fname):
        # Ensure the output directory exists
        dname = os.path.dirname(os.path.abspath(fname))
        if dname and not os.path.isdir(dname):
            try:
                os.makedirs(dname)
            except OSError as exc:
                raise OSError(
                    f"Output path {fname} contains a directory {dname} "
                    f"that cannot be created: {exc}"
                ) from exc
    if os.path.isfile(fname):
        # Add an integer suffix to the existing file name
        cnt = 1
        bak_fname = f"{fname}.{cnt}"
        while os.path.isfile(bak_fname):
            cnt += 1
            bak_fname = f"{fname}.{cnt}"
        os.rename(fname, bak_fname)
        logging.info("Moved existing file %s -> %s", fname, bak_fname)
    return True


@contextlib.contextmanager
def temp_write_text(text, mode="w+b"):
    """Save text to a temporary file.

    NB: This won't work on Windows b/c the file stays open.
    """
    with tempfile.NamedTemporaryFile(mode=mode) as tmp:
        tmp.write(text)
        tmp.flush()
        yield tmp.name


# __________________________________________________________________________
# More helpers


def assert_equal(msg, **values) -> None:
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
    msg += f": {key1} = {val1!r}"
    for okey, oval in values.items():
        msg += f", {okey} = {oval!r}"
        if oval != val1:
            ok = False
    if not ok:
        raise ValueError(msg)


def check_unique(items, title):
    """Ensure all items in an iterable are identical; return that one item."""
    its = set(items)
    assert len(its) == 1, "Inconsistent %s keys: %s" % (
        title,
        " ".join(map(str, sorted(its))),
    )
    return its.pop()


def fbase(fname: str) -> str:
    """Strip directory and all extensions from a filename."""
    base = os.path.basename(fname)
    # Gzip extension usually follows another extension
    if base.endswith(".gz"):
        base = base[:-3]
    # Cases to drop more than just the last dot
    known_multipart_exts = (
        ".antitargetcoverage.cnn",
        ".targetcoverage.cnn",
        ".antitargetcoverage.csv",
        ".targetcoverage.csv",
        # Pipeline suffixes
        ".recal.bam",
        ".deduplicated.realign.bam",
    )
    for ext in known_multipart_exts:
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    else:
        base = base.rsplit(".", 1)[0]
    return base
