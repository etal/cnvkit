"""Filter copy number segments."""
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


def squash_by_groups(cnarr, levels, by_arm=False):
    """Reduce CopyNumArray rows to a single row within each given level."""
    # Enumerate runs of identical values
    change_levels = enumerate_changes(levels)
    assert (change_levels.index == levels.index).all()
    assert cnarr.data.index.is_unique
    assert levels.index.is_unique
    assert change_levels.index.is_unique
    if by_arm:
        # Enumerate chromosome arms
        arm_levels = []
        for i, (_chrom, cnarm) in enumerate(cnarr.by_arm()):
            arm_levels.append(np.repeat(i, len(cnarm)))
        change_levels += np.concatenate(arm_levels)
    else:
        # Enumerate chromosomes
        chrom_names = cnarr['chromosome'].unique()
        chrom_col = (cnarr['chromosome']
                     .replace(chrom_names, np.arange(len(chrom_names))))
        change_levels += chrom_col
    data = cnarr.data.assign(_group=change_levels)
    groupkey = ['_group']
    if 'cn1' in cnarr:
        # Keep allele-specific CNAs separate
        data['_g1'] = enumerate_changes(cnarr['cn1'])
        data['_g2'] = enumerate_changes(cnarr['cn2'])
        groupkey.extend(['_g1', '_g2'])
    data = (data.groupby(groupkey, as_index=False, group_keys=False, sort=False)
            .apply(squash_region)
            .reset_index(drop=True))
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
    assert 'weight' in cnarr
    out = {'chromosome': [cnarr['chromosome'].iat[0]],
           'start': cnarr['start'].iat[0],
           'end': cnarr['end'].iat[-1],
          }
    region_weight = cnarr['weight'].sum()
    if region_weight > 0:
        out['log2'] = np.average(cnarr['log2'], weights=cnarr['weight'])
    else:
        out['log2'] = np.mean(cnarr['log2'])
    out['gene'] = ','.join(cnarr['gene'].drop_duplicates())
    out['probes'] = cnarr['probes'].sum() if 'probes' in cnarr else len(cnarr)
    out['weight'] = region_weight
    if 'depth' in cnarr:
        if region_weight > 0:
            out['depth'] = np.average(cnarr['depth'], weights=cnarr['weight'])
        else:
            out['depth'] = np.mean(cnarr['depth'])
    if 'baf' in cnarr:
        if region_weight > 0:
            out['baf'] = np.average(cnarr['baf'], weights=cnarr['weight'])
        else:
            out['baf'] = np.mean(cnarr['baf'])
    if 'cn' in cnarr:
        if region_weight > 0:
            out['cn'] = weighted_median(cnarr['cn'], cnarr['weight'])
        else:
            out['cn'] = np.median(cnarr['cn'])
        if 'cn1' in cnarr:
            if region_weight > 0:
                out['cn1'] = weighted_median(cnarr['cn1'], cnarr['weight'])
            else:
                out['cn1'] = np.median(cnarr['cn1'])
            out['cn2'] = out['cn'] - out['cn1']
    if 'p_bintest' in cnarr:
        # Only relevant for single-bin segments, but this seems safe/conservative
        out['p_bintest'] = cnarr['p_bintest'].max()
    return pd.DataFrame(out)


@require_column('cn')
def ampdel(segarr):
    """Merge segments by amplified/deleted/neutral copy number status.

    Follow the clinical reporting convention:

    - 5+ copies (2.5-fold gain) is amplification
    - 0 copies is homozygous/deep deletion
    - CNAs of lesser degree are not reported

    This is recommended only for selecting segments corresponding to
    actionable, usually focal, CNAs. Any real and potentially informative but
    lower-level CNAs will be dropped.
    """
    levels = np.zeros(len(segarr))
    levels[segarr['cn'] == 0] = -1
    levels[segarr['cn'] >= 5] = 1
    # or: segarr['log2'] >= np.log2(2.5)
    cnarr = squash_by_groups(segarr, pd.Series(levels))
    return cnarr[(cnarr['cn'] == 0) | (cnarr['cn'] >= 5)]


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
    levels = np.zeros(len(segarr))
    # if len(segarr) < 10:
    #     logging.warning("* segarr :=\n%s", segarr)
    #     logging.warning("* segarr['ci_lo'] :=\n%s", segarr['ci_lo'])
    #     logging.warning("* segarr['ci_lo']>0 :=\n%s", segarr['ci_lo'] > 0)
    levels[segarr['ci_lo'].values > 0] = 1
    levels[segarr['ci_hi'].values < 0] = -1
    return squash_by_groups(segarr, pd.Series(levels))


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
    levels = np.zeros(len(segarr))
    levels[segarr['log2'] - margin > 0] = 1
    levels[segarr['log2'] + margin < 0] = -1
    return squash_by_groups(segarr, pd.Series(levels))
