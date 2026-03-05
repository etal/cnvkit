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


def hierarchical(samples, min_cluster_size=5):
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
    from scipy.spatial.distance import pdist

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
    labels = fcluster(Z, t=threshold, criterion="inconsistent", depth=depth)
    logging.info(
        "Hierarchical clustering: %d samples -> %d clusters (threshold=%.3f)",
        n_samples,
        len(set(labels)),
        threshold,
    )
    return _labels_to_clusters(labels, min_cluster_size)


def _find_inconsistency_threshold(Z, incon, min_cluster_size, n_samples, depth):
    """Find an inconsistency threshold that yields valid clusters.

    Search for the lowest threshold where every cluster has at least
    min_cluster_size samples.
    """
    from scipy.cluster.hierarchy import fcluster

    best_threshold = 10.0  # Default: single cluster (very permissive)

    for threshold in np.arange(0.5, 3.1, 0.1):
        labels = fcluster(Z, t=threshold, criterion="inconsistent", depth=depth)
        n_clusters = len(set(labels))
        if n_clusters <= 1:
            continue
        # Check that all clusters are large enough
        label_counts = np.bincount(labels)  # index 0 unused (labels are 1-based)
        min_size = label_counts[1:].min()
        if min_size >= min_cluster_size:
            best_threshold = threshold
            break

    if best_threshold >= 10.0:
        logging.info(
            "All samples are highly correlated; using single reference cluster"
        )

    return best_threshold


def _labels_to_clusters(labels, min_cluster_size):
    """Convert scipy's 1-indexed label array to a list of index arrays.

    Clusters smaller than min_cluster_size are merged into the nearest
    larger cluster.
    """
    unique_labels = np.unique(labels)
    clusters = []
    small_indices: list[int] = []
    for label in unique_labels:
        indices = np.where(labels == label)[0]
        if len(indices) >= min_cluster_size:
            clusters.append(indices)
            logging.info("Cluster %d: %d samples", label, len(indices))
        else:
            logging.info(
                "Cluster %d: %d samples (below min size %d, merging)",
                label,
                len(indices),
                min_cluster_size,
            )
            small_indices.extend(indices)

    if not clusters:
        # All clusters too small -> return everything as one cluster
        return [np.arange(len(labels))]

    if small_indices:
        # Merge small-cluster samples into the first (largest) cluster
        largest_idx = max(range(len(clusters)), key=lambda i: len(clusters[i]))
        clusters[largest_idx] = np.concatenate(
            [clusters[largest_idx], np.array(small_indices)]
        )
        logging.info(
            "Merged %d small-cluster samples into cluster %d",
            len(small_indices),
            largest_idx + 1,
        )

    return clusters


def kmedoids(samples, k=None, min_cluster_size=1):
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
        Minimum number of samples per cluster.

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

    return _labels_to_clusters(labels, min_cluster_size)


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
    medoids = []
    # First medoid: point with smallest total distance to all others
    medoids.append(int(np.argmin(dist_matrix.sum(axis=1))))
    for _ in range(1, k):
        # Next medoid: point that reduces total distance the most
        non_medoids = [i for i in range(n) if i not in medoids]
        # Current distance from each point to nearest medoid
        current_dists = dist_matrix[:, medoids].min(axis=1)
        best_gain = -1.0
        best_candidate = non_medoids[0]
        for candidate in non_medoids:
            # Gain: reduction in distance if candidate becomes a medoid
            gain = np.sum(np.maximum(0, current_dists - dist_matrix[:, candidate]))
            if gain > best_gain:
                best_gain = gain
                best_candidate = candidate
        medoids.append(best_candidate)

    # SWAP: iteratively improve medoids
    max_iter = 100
    for _iteration in range(max_iter):
        # Assign each point to nearest medoid
        labels = np.argmin(dist_matrix[:, medoids], axis=1)
        total_cost = sum(dist_matrix[i, medoids[labels[i]]] for i in range(n))

        improved = False
        for m_idx in range(k):
            cluster_members = np.where(labels == m_idx)[0]
            for candidate in cluster_members:
                if candidate == medoids[m_idx]:
                    continue
                # Try swapping
                new_medoids = medoids.copy()
                new_medoids[m_idx] = candidate
                new_labels = np.argmin(dist_matrix[:, new_medoids], axis=1)
                new_cost = sum(
                    dist_matrix[i, new_medoids[new_labels[i]]] for i in range(n)
                )
                if new_cost < total_cost:
                    medoids = new_medoids
                    total_cost = new_cost
                    improved = True
                    break
            if improved:
                break
        if not improved:
            break

    # Final assignment (1-indexed labels)
    labels = np.argmin(dist_matrix[:, medoids], axis=1) + 1
    return labels


def _select_k(dist_matrix, n_samples):
    """Select k for k-medoids by maximizing silhouette score.

    Tries k from 2 to min(10, n_samples-1) and picks the k with the highest
    silhouette score. Falls back to k=2 if all scores are negative.
    """
    from sklearn.metrics import silhouette_score

    max_k = min(10, n_samples - 1)
    if max_k < 2:
        return 1

    best_k = 2
    best_score = -1.0

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

    logging.info("Selected k=%d (silhouette=%.4f)", best_k, best_score)
    return best_k
