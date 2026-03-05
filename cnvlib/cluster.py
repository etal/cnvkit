#!/usr/bin/env python
"""Cluster control samples by coverage correlation.

Supports hierarchical clustering (default) and k-medoids.

References:

    Guo et al. 2019, "Comparative study of exome copy number variant detection
    tools using the K-means clustering algorithm" (PMC6537193):
    - k-means with Pearson correlation of read depth is ~200x faster than kNN
      with equivalent accuracy; optimal k=4-5; clusters correspond to
      sequencing batches/centers
"""

import logging

import numpy as np


def hierarchical(samples, min_cluster_size=4):
    """Cluster samples using hierarchical (UPGMA) clustering on Pearson correlation.

    Uses the inconsistency statistic to automatically determine the number of
    clusters, which adapts to the data without requiring a pre-specified k.

    Parameters
    ----------
    samples : 2D numpy array
        Matrix of samples' log2 values. Rows are samples, columns are bins.
    min_cluster_size : int
        Minimum number of samples per cluster. Clusters smaller than this are
        absorbed into the nearest larger cluster.

    Returns
    -------
    list of numpy arrays
        Each array contains the indices of samples in that cluster.
    """
    from scipy.cluster.hierarchy import cophenet, fcluster, inconsistent, linkage
    from scipy.spatial.distance import pdist, squareform

    n_samples = len(samples)
    if n_samples < 2:
        return [np.array(range(n_samples))]

    dists = pdist(samples, metric="correlation")  # 1 - Pearson r
    # Replace NaN distances (from constant-coverage bins) with max distance
    dists = np.nan_to_num(dists, nan=1.0)
    Z = linkage(dists, method="average")  # UPGMA

    # Log cophenetic correlation coefficient as a clustering quality metric
    coph_corr, _ = cophenet(Z, dists)
    logging.info("Cophenetic correlation coefficient: %.4f", coph_corr)

    depth = 2
    incon = inconsistent(Z, d=depth)
    threshold = _find_inconsistency_threshold(
        Z, incon, min_cluster_size, n_samples, depth
    )
    if threshold is None:
        return [np.arange(n_samples)]
    labels = fcluster(Z, t=threshold, criterion="inconsistent", depth=depth)
    logging.info(
        "Hierarchical clustering: %d samples -> %d clusters (threshold=%.3f)",
        n_samples,
        len(set(labels)),
        threshold,
    )
    dist_matrix = squareform(dists)
    return _labels_to_clusters(labels, min_cluster_size, dist_matrix)


def _find_inconsistency_threshold(Z, incon, min_cluster_size, n_samples, depth):
    """Find an inconsistency threshold that yields valid clusters.

    Search for the lowest threshold where every cluster has at least
    min_cluster_size samples. Returns None if no valid split is found,
    indicating all samples should be in a single cluster.
    """
    from scipy.cluster.hierarchy import fcluster

    # Use actual inconsistency values from the linkage as candidate thresholds
    candidates = np.unique(incon[:, 3])
    candidates = candidates[candidates > 0]
    if len(candidates) == 0:
        logging.info(
            "All samples are highly correlated; using single reference cluster"
        )
        return None

    for threshold in candidates:
        labels = fcluster(Z, t=threshold, criterion="inconsistent", depth=depth)
        n_clusters = len(set(labels))
        if n_clusters <= 1:
            continue
        # Check that all clusters are large enough
        label_counts = np.bincount(labels)  # index 0 unused (labels are 1-based)
        min_size = label_counts[1:].min()
        if min_size >= min_cluster_size:
            return threshold

    logging.info("All samples are highly correlated; using single reference cluster")
    return None


def _labels_to_clusters(labels, min_cluster_size, dist_matrix):
    """Convert scipy's 1-indexed label array to a list of index arrays.

    Clusters smaller than min_cluster_size are merged into the nearest
    larger cluster, measured by average distance between cluster members.
    """
    unique_labels = np.unique(labels)
    large_clusters: list[np.ndarray] = []
    small_clusters: list[np.ndarray] = []
    for label in unique_labels:
        indices = np.where(labels == label)[0]
        if len(indices) >= min_cluster_size:
            large_clusters.append(indices)
            logging.info("Cluster %d: %d samples", label, len(indices))
        else:
            logging.info(
                "Cluster %d: %d samples (below min size %d, merging)",
                label,
                len(indices),
                min_cluster_size,
            )
            small_clusters.append(indices)

    if not large_clusters:
        # All clusters too small -> return everything as one cluster
        return [np.arange(len(labels))]

    # Merge each small cluster into its nearest large cluster
    for small_idx in small_clusters:
        # Average distance from small cluster members to each large cluster
        best_target = 0
        best_dist = np.inf
        for i, large_idx in enumerate(large_clusters):
            avg_dist = dist_matrix[np.ix_(small_idx, large_idx)].mean()
            if avg_dist < best_dist:
                best_dist = avg_dist
                best_target = i
        large_clusters[best_target] = np.concatenate(
            [large_clusters[best_target], small_idx]
        )
        logging.info(
            "Merged %d small-cluster samples into cluster %d (avg dist=%.4f)",
            len(small_idx),
            best_target + 1,
            best_dist,
        )

    return large_clusters


def kmedoids(samples, k=None, min_cluster_size=4):
    """Cluster samples using k-medoids (PAM) on Pearson correlation distance.

    When k is not specified, it is selected by maximizing the silhouette score over a
    range of candidate values.

    Parameters
    ----------
    samples : 2D numpy array
        Matrix of samples' log2 values. Rows are samples, columns are bins.
    k : int, optional
        Number of clusters. If None, selected via silhouette score.
    min_cluster_size : int
        Minimum number of samples per cluster. Clusters smaller than this are
        absorbed into the nearest larger cluster.

    Returns
    -------
    list of numpy arrays
        Each array contains the indices of samples in that cluster.
    """
    from scipy.spatial.distance import pdist, squareform

    n_samples = len(samples)
    if n_samples < 2:
        return [np.array(range(n_samples))]

    dists = pdist(samples, metric="correlation")  # 1 - Pearson r
    dists = np.nan_to_num(dists, nan=1.0)
    dist_matrix = squareform(dists)

    if k is None:
        k = _select_k(dist_matrix, n_samples)

    logging.info(
        "Clustering %d samples by k-medoids (k=%d) on correlation distance",
        n_samples,
        k,
    )
    labels = _pam(dist_matrix, k)

    return _labels_to_clusters(labels, min_cluster_size, dist_matrix)


def _pam(dist_matrix, k):
    """Partitioning Around Medoids (PAM) on a precomputed distance matrix.

    A simple implementation of the PAM algorithm:
    1. Initialize medoids using greedy selection (BUILD phase).
    2. Iteratively swap medoids with non-medoids to minimize total cost (SWAP).

    Parameters
    ----------
    dist_matrix : 2D numpy array
        Symmetric distance matrix (n x n).
    k : int
        Number of clusters.

    Returns
    -------
    labels : numpy array
        1-indexed cluster labels for each sample.
    """
    n = len(dist_matrix)
    if k >= n:
        return np.arange(1, n + 1)

    # BUILD: greedy medoid initialization
    medoids = [int(np.argmin(dist_matrix.sum(axis=1)))]
    for _ in range(1, k):
        # Current distance from each point to nearest medoid
        current_dists = dist_matrix[:, medoids].min(axis=1)
        # Gain for each candidate: reduction in total distance
        gains = np.maximum(0, current_dists[:, np.newaxis] - dist_matrix).sum(axis=0)
        # Exclude existing medoids
        gains[medoids] = -1.0
        medoids.append(int(np.argmax(gains)))

    medoid_arr = np.array(medoids)

    # SWAP: iteratively improve medoids
    for _iteration in range(100):
        # Assign each point to nearest medoid (vectorized)
        dists_to_medoids = dist_matrix[:, medoid_arr]
        labels = np.argmin(dists_to_medoids, axis=1)
        total_cost = dists_to_medoids[np.arange(n), labels].sum()

        improved = False
        non_medoid_mask = np.ones(n, dtype=bool)
        non_medoid_mask[medoid_arr] = False
        for m_idx in range(k):
            for candidate in np.where(non_medoid_mask)[0]:
                # Try swapping (vectorized cost)
                new_medoid_arr = medoid_arr.copy()
                new_medoid_arr[m_idx] = candidate
                new_dists = dist_matrix[:, new_medoid_arr]
                new_cost = new_dists.min(axis=1).sum()
                if new_cost < total_cost:
                    medoid_arr = new_medoid_arr
                    total_cost = new_cost
                    improved = True
                    break
            if improved:
                break
        if not improved:
            break

    # Final assignment (1-indexed labels)
    labels = np.argmin(dist_matrix[:, medoid_arr], axis=1) + 1
    return labels


def _select_k(dist_matrix, n_samples):
    """Select k for k-medoids by maximizing silhouette score.

    Tries k from 2 to min(10, n_samples-1) and picks the k with the highest
    silhouette score. Short-circuits if the score declines for 2 consecutive
    values of k. Falls back to k=2 if all scores are negative.
    """
    from sklearn.metrics import silhouette_score

    max_k = min(10, n_samples - 1)
    if max_k < 2:
        return 1

    best_k = 2
    best_score = -1.0
    consecutive_declines = 0

    for k in range(2, max_k + 1):
        labels = _pam(dist_matrix, k)
        # Need at least 2 distinct labels for silhouette
        if len(set(labels)) < 2:
            continue
        score = silhouette_score(dist_matrix, labels, metric="precomputed")
        logging.info("  k=%d: silhouette=%.4f", k, score)
        if score > best_score:
            best_score = score
            best_k = k
            consecutive_declines = 0
        else:
            consecutive_declines += 1
            if consecutive_declines >= 2:
                break

    logging.info("Selected k=%d (silhouette=%.4f)", best_k, best_score)
    return best_k
