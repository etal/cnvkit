"""Z-test for single-bin copy number alterations."""
import logging

import numpy as np
import pandas as pd
from scipy.stats import norm

from . import params, segfilters


def do_bintest(cnarr, segments=None, alpha=0.005, target_only=False):
    """Get a probability for each bin based on its Z-score.

    Adds a column w/ p-values to the input .cnr. With `segments`, the Z-score is
    relative to the enclosing segment's mean, otherwise it is relative to 0.

    Bin p-values are corrected for multiple hypothesis testing by the
    Benjamini-Hochberg method.

    Returns: bins where the probability < `alpha`.
    """
    cnarr = cnarr.copy()
    # Subtract segment means, if given, to report only the CNA bins that
    # weren't already detected (including exon-size CNAs within a 
    # larger-scale, smaller-amplitude CNA)
    resid = cnarr.residuals(segments)
    if not resid.index.is_unique:
        # Overlapping segments, maybe?
        dup_idx = resid.index.duplicated(keep=False)
        logging.warning("Segments may overlap at %d bins; dropping duplicate values",
                        dup_idx.sum())
        logging.debug("Duplicated indices: %s", " ".join(map(str, resid[dup_idx].head(50))))
        resid = resid[~resid.index.duplicated()]
        cnarr = cnarr.as_dataframe(cnarr.data.loc[resid.index])
    if len(cnarr) != len(resid):
        logging.info("Segments do not cover all bins (%d), only %d of them",
                     len(cnarr), len(resid))
        cnarr = cnarr.as_dataframe(cnarr.data.loc[resid.index])

    cnarr['log2'] = resid
    cnarr['probes'] = 1

    if target_only:
        antitarget_idx = cnarr['gene'].isin(params.ANTITARGET_ALIASES)
        if antitarget_idx.any():
            logging.info("Ignoring %d off-target bins", antitarget_idx.sum())
            # NB: bins no longer match the original input
            cnarr = cnarr[~antitarget_idx]

    cnarr['p_bintest'] = z_prob(cnarr)
    is_sig = cnarr['p_bintest'] < alpha
    logging.info("Significant hits in {}/{} bins ({:.3g}%)"
                 .format(is_sig.sum(), len(is_sig),
                         100 * is_sig.sum() / len(is_sig)))
    # if segments:
    #     return spike_into_segments(cnarr, segments, is_sig)

    # May be empty
    hits = cnarr[is_sig]
    return hits


def z_prob(cnarr):
    """Calculate z-test p-value at each bin."""
    # Bin weights ~ 1-variance; bin log2 values already centered at 0.0
    sd = np.sqrt(1 - cnarr['weight'])
    # Convert to Z-scores
    z = cnarr['log2'] / sd
    # Two-sided survival function (1-CDF) probability
    p = 2. * norm.cdf(-np.abs(z))
    # Similar to the above -- which is better?
    # p2 = 2 * norm.pdf(cnarr['log2'], loc=0, scale=sd)
    # if not np.allclose(p, p2):
    #     print("Max diff:", np.abs(p - p2).max())
    #     print("Median diff:", np.median(np.abs(p - p2)))
    #     print("Ratio:", (p / p2).mean())
    # Correct for multiple hypothesis tests
    return p_adjust_bh(p)


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def spike_into_segments(cnarr, segments, is_sig):
    if is_sig.any():
        # Splice significant hits into the given segments
        # NB: residuals() above ensures hits all occur within segments
        cnarr['is_sig'] = is_sig
        chunks = []
        for segment, seghits in cnarr.by_ranges(segments, keep_empty=True):
            if seghits['is_sig'].any():
                # Merge each run of adjacent non-significant bins within this
                # segment, leaving the significant hits as single-bin segments
                levels = seghits['is_sig'].cumsum() * seghits['is_sig']
                chunks.append(seghits.data
                                .assign(_levels=levels)
                                .groupby('_levels', sort=False)
                                .apply(segfilters.squash_region)
                                .reset_index(drop=True))
            else:
                # Keep this segment as-is
                chunks.append(pd.DataFrame.from_records([segment],
                                                        columns=segments.data.columns))
        return cnarr.as_dataframe(pd.concat(chunks,
                                            # pandas 0.23+
                                            #sort=False
                                            ))
    else:
        # Nothing to do
        return segments
