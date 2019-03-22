"""Segmentation by Hidden Markov Model."""
import collections
import logging

import numpy as np
import pandas as pd
import scipy.special
import pomegranate as pom

from ..descriptives import biweight_midvariance
from ..segfilters import squash_by_groups


def segment_hmm(cnarr, method, window=None, processes=1):
    """Segment bins by Hidden Markov Model.

    Use Viterbi method to infer copy number segments from sequential data.

    With b-allele frequencies ('baf' column in `cnarr`), jointly segment
    log-ratios and b-allele frequencies across a chromosome.

    Parameters
    ----------
    cnarr : CopyNumArray
        The bin-level data to segment.
    method : string
        One of 'hmm' (3 states, flexible means), 'hmm-tumor' (5 states, flexible
        means), 'hmm-germline' (3 states, fixed means).

    Results
    -------
    segarr : CopyNumArray
        The segmented data.
    """
    # NB: Incorporate weights into smoothed log2 estimates
    # (Useful kludge until weighted HMM is in place)
    orig_log2 = cnarr['log2'].values.copy()
    cnarr['log2'] = cnarr.smoothed()#window)

    logging.info("Building model from observations")
    model = hmm_get_model(cnarr, method, processes)

    logging.info("Predicting states from model")
    observations = as_observation_matrix(cnarr)
    states = np.concatenate([np.array(model.predict(obs, algorithm='map'))
                             for obs in observations])

    logging.info("Done, now finalizing")
    logging.debug("Model states: %s", model.states)
    logging.debug("Predicted states: %s", states[:100])
    logging.debug(str(collections.Counter(states)))
    logging.debug("Observations: %s", observations[0][:100])
    logging.debug("Edges: %s", model.edges)

    # Merge adjacent bins with the same state to create segments
    cnarr['log2'] = orig_log2
    cnarr['probes'] = 1
    segarr = squash_by_groups(cnarr, pd.Series(states), by_arm=True)
    return segarr


def hmm_get_model(cnarr, method, processes):
    """

    Parameters
    ----------
    cnarr : CopyNumArray
        The normalized bin-level values to be segmented.
    method : string
        One of 'hmm', 'hmm-tumor', 'hmm-germline'.
    processes : int
        Number of parallel jobs to run.

    Returns
    -------
    model :
        A pomegranate HiddenMarkovModel trained on the given dataset.
    """
    assert method in ('hmm-tumor', 'hmm-germline', 'hmm')
    observations = as_observation_matrix(cnarr.autosomes())

    # Estimate standard deviation from the full distribution, robustly
    stdev = biweight_midvariance(np.concatenate(observations), initial=0)
    if method == 'hmm-germline':
        state_names = ["loss", "neutral", "gain"]
        distributions = [
            pom.NormalDistribution(-1.0, stdev, frozen=True),
            pom.NormalDistribution(0.0, stdev, frozen=True),
            pom.NormalDistribution(0.585, stdev, frozen=True),
        ]
    elif method == 'hmm-tumor':
        state_names = ["del", "loss", "neutral", "gain", "amp"]
        distributions = [
            pom.NormalDistribution(-2.0, stdev, frozen=False),
            pom.NormalDistribution(-0.5, stdev, frozen=False),
            pom.NormalDistribution(0.0, stdev, frozen=True),
            pom.NormalDistribution(0.3, stdev, frozen=False),
            pom.NormalDistribution(1.0, stdev, frozen=False),
        ]
    else:
        state_names = ["loss", "neutral", "gain"]
        distributions = [
            pom.NormalDistribution(-1.0, stdev, frozen=False),
            pom.NormalDistribution(0.0, stdev, frozen=False),
            pom.NormalDistribution(0.585, stdev, frozen=False),
        ]

    n_states = len(distributions)
    # Starts -- prefer neutral
    binom_coefs = scipy.special.binom(n_states - 1, range(n_states))
    start_probabilities = binom_coefs / binom_coefs.sum()
    # Ends -- equally likely
    #end_probabilities = np.ones(n_states) / n_states

    # Prefer to keep the current state in each transition
    # All other transitions are equally likely, to start
    transition_matrix = (np.identity(n_states) * 100
                        + np.ones((n_states, n_states)) / n_states)

    model = pom.HiddenMarkovModel.from_matrix(transition_matrix, distributions,
        start_probabilities, state_names=state_names, name=method)

    model.fit(sequences=observations,
              weights=[len(obs) for obs in observations],
              distribution_inertia = .8,  # Allow updating dists, but slowly
              edge_inertia=0.1,
              # lr_decay=.75,
              pseudocount=5,
              use_pseudocount=True,
              max_iterations=100000,
              n_jobs=processes,
              verbose=False)
    return model


def as_observation_matrix(cnarr, variants=None):
    """Extract HMM fitting values from `cnarr`.

    For each chromosome arm, extract log2 ratios as a numpy array.

    Future: If VCF of variants is given, or 'baf' column has already been
    added to `cnarr` from the same, then the BAF values are a second row/column
    in each numpy array.

    Returns: List of numpy.ndarray, one per chromosome arm.
    """
    # TODO incorporate weights -- currently handled by smoothing
    # TODO incorporate inter-bin distances
    observations = [arm.log2.values
                    for _c, arm in cnarr.by_arm()]
    return observations
    # --- TODO incorporate variant BAF/zygosity/LOH/ROH ---
    if variants is not None:
        print("Variants!")
        cnarr['baf'] = variants.baf_by_ranges(cnarr)
    if 'baf' in cnarr:
        print("BAF!")
        observations = [np.hstack([arm.log2.values, arm['baf'].values])
                        for _c, arm in cnarr.by_arm()]
    # /---
