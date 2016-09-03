"""Filter copy number segments."""
from __future__ import absolute_import, division, print_function

import functools
import logging

import numpy as np
import pandas as pd

from .descriptives import weighted_median


def require_column(*colnames):
    """Wrapper to coordinate the segment-filtering functions.

    Verify that the given columns are in the CopyNumArray the wrapped function
    takes. Also log the number of rows in the array before and after filtration.
    """
    if len(colnames) == 1:
        msg = "'{}' filter requires column '{}'"
    else:
        msg = "'{}' filter requires columns " + \
                ", ".join(["'{}'"] * len(colnames))
    def wrap(func):
        @functools.wraps(func)
        def wrapped_f(segarr):
            filtname = func.__name__
            if any(c not in segarr for c in colnames):
                raise ValueError(msg.format(filtname, *colnames))
            result = func(segarr)
            logging.info("Filtered by '%s' from %d to %d rows",
                         filtname, len(segarr), len(result))
            return result
        return wrapped_f
    return wrap


def squash_by_groups(cnarr, levels):
    """Reduce CopyNumArray rows to a single row within each given level."""
    # Enumerate runs of identical values
    change_levels = enumerate_changes(levels)
    # Enumerate chromosomes
    chrom_names = cnarr['chromosome'].unique()
    change_levels += cnarr['chromosome'].replace(chrom_names,
                                                 np.arange(len(chrom_names)))
    data = cnarr.data.assign(_group=change_levels)
    groupkey = ['_group']
    if 'cn1' in cnarr:
        # Keep allele-specific CNAs separate
        data['_g1'] = enumerate_changes(cnarr['cn1'])
        data['_g2'] = enumerate_changes(cnarr['cn2'])
        groupkey.extend(['_g1', '_g2'])
    data = data.groupby(groupkey, sort=False).apply(squash_region)
    return cnarr.as_dataframe(data)


def enumerate_changes(levels):
    """Assign a unique integer to each run of identical values.

    Repeated but non-consecutive values will be assigned different integers.
    """
    return levels.diff().fillna(0).abs().cumsum().astype(int)


def squash_region(cnarr):
    """Reduce a CopyNumArray to 1 row, keeping fields sensible.

    Most fields added by the `segmetrics` command will be dropped.
    """
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
    if 'baf' in cnarr:
        out['baf'] = np.average(cnarr['baf'], weights=cnarr['weight'])
    if 'cn' in cnarr:
        out['cn'] = weighted_median(cnarr['cn'], cnarr['weight'])
        if 'cn1' in cnarr:
            out['cn1'] = weighted_median(cnarr['cn1'], cnarr['weight'])
            out['cn2'] = out['cn'] - out['cn1']
    return pd.DataFrame(out)


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
    return squash_by_groups(segarr, levels)


@require_column('depth')
def bic(segarr):
    """Merge segments by Bayesian Information Criterion.

    See: BIC-seq (Xi 2011), doi:10.1073/pnas.1110574108
    """
    return NotImplemented


@require_column('ci_lo', 'ci_hi')
def ci(segarr):
    """Merge segments by confidence interval (overlapping 0).

    Segments with lower CI above 0 are kept as gains, upper CI below 0 as
    losses, and the rest with CI overlapping zero are collapsed as neutral.
    """
    levels = pd.Series(np.zeros(len(segarr)))
    levels[segarr['ci_lo'] > 0] = 1
    levels[segarr['ci_hi'] < 0] = -1
    return squash_by_groups(segarr, levels)


@require_column('cn')
def cn(segarr):
    """Merge segments by integer copy number."""
    return squash_by_groups(segarr, segarr['cn'])


@require_column('sem')
def sem(segarr, zscore=1.96):
    """Merge segments by Standard Error of the Mean (SEM).

    Use each segment's SEM value to estimate a 95% confidence interval (via
    `zscore`). Segments with lower CI above 0 are kept as gains, upper CI below
    0 as losses, and the rest with CI overlapping zero are collapsed as neutral.
    """
    margin = segarr['sem'] * zscore
    levels = pd.Series(np.zeros(len(segarr)))
    levels[segarr['log2'] - margin > 0] = 1
    levels[segarr['log2'] + margin < 0] = -1
    return squash_by_groups(segarr, levels)
