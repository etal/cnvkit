"""Filter copy number segments."""
from __future__ import absolute_import, division, print_function

import functools
import logging

import numpy as np
import pandas as pd

from .metrics import segment_mean

def require_column(*colnames):
    """Ensure these columns are in the CopyNumArray the wrapped function takes.

    Also log the number of rows in the array before and after filtration.
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


def squash_region(cnarr):
    """Reduce a CopyNumArray to 1 row, keeping fields sensible."""
    assert 'weight' in cnarr and 'probes' in cnarr
    out = {'chromosome': [cnarr['chromosome'].iat[0]],
           'start': cnarr['start'].iat[0],
           'end': cnarr['end'].iat[-1],
          }
    out['log2'] = np.average(cnarr['log2'], weights=cnarr['weight'])
    out['gene'] = ','.join(cnarr['gene'].drop_duplicates())
    out['probes'] = cnarr['probes'].sum()
    out['weight'] = cnarr['weight'].sum()
    if 'depth' in cnarr:
        out['depth'] = np.average(cnarr['depth'], weights=cnarr['weight'])
    return pd.DataFrame(out)


def enumerate_chroms(chroms):
    names = chroms.unique()
    nums = np.arange(len(names))
    return chroms.replace(names, nums)


@require_column('cn')
def ampdel(segarr):
    """Merge segments by amplified/deleted/neutral copy number status.

    Follow the clinical reporting convention:

    - 5+ copies (2.5-fold gain) is amplification
    - 0 copies is homozygous/deep deletion
    - CNAs of lesser degree are not reported

    This is recommended only for selecting segments corresponding to
    actionable, usually focal, CNAs. Real and potentially informative but
    lower-level CNAs will be merged together.
    """
    levels = pd.Series(np.zeros(len(segarr)))
    levels[segarr['cn'] == 0] = -1
    levels[segarr['cn'] >= 5] = 1
    # or: segarr['log2'] >= np.log2(2.5)
    # TODO - handle a/b allele amplifications separately
    #   i.e. don't merge amplified segments if cn1, cn2 are not the same

    groups = levels.diff().fillna(0).abs().cumsum().astype(int)
    groups += enumerate_chroms(segarr['chromosome'])

    squashed = (segarr.data.assign(group=groups)
                .groupby('group', as_index=False, group_keys=False, sort=False)
                .apply(squash_region))
    return segarr.as_dataframe(squashed)


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
