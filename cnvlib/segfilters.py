"""Filter copy number segments."""
from __future__ import absolute_import, division, print_function

import functools
import logging


def require_column(*colnames):
    """Ensure these columns are in the CopyNumArray the wrapped function takes.
    """
    if len(colnames) == 1:
        msg = "'{}' filter requires column '{}'"
    else:
        msg = "'{}' filter requires columns " + \
                ", ".join(["'{}'"] * len(colnames))
    def wrap(func):
        @functools.wraps(func)
        def wrapped_f(segarr):
            filtname = func.func_name
            if any(c not in segarr for c in colnames):
                raise ValueError(msg.format(filtname, *colnames))
            result = func(segarr)
            logging.info("Filtered by '%s' from %d to %d rows",
                         filtname, len(segarr), len(result))
            return result
        return wrapped_f
    return wrap


@require_column('cn')
def ampdel(segarr):
    """Merge segments by amplified/deleted/neutral copy number status.

    Follow the clinical reporting convention:

    - 5+ copies (2.5-fold gain) is amplification
    - 0 copies is homozygous/deep deletion
    - CNAs of lesser degree are not reported

    This is recommended only for highlighting segments corresponding to
    actionable, usually focal, CNAs. Real and potentially informative but
    lower-level CNAs will be merged together.
    """
    return segarr


@require_column('depth')
def bic(segarr):
    """Merge segments by Bayesian Information Criterion.

    See: BIC-seq
    """
    return segarr


@require_column('ci_lo', 'ci_hi')
def ci(segarr):
    """Merge segments by confidence interval (overlapping 0)."""
    return segarr


@require_column('cn')
def cn(segarr):
    """Merge segments by integer copy number."""
    return segarr


@require_column('sem')
def sem(segarr):
    """Merge segments by Standard Error of the Mean (SEM)."""
    return segarr
