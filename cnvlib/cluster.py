#!/usr/bin/env python
"""Cluster control samples by coverage correlation.

Supports hierarchical clustering (default), k-means, and Markov clustering (MCL).

References:

    Guo et al. 2019, "Comparative study of exome copy number variant detection
    tools using the K-means clustering algorithm" (PMC6537193):
    - k-means with Pearson correlation of read depth is ~200x faster than kNN
      with equivalent accuracy; optimal k=4-5; clusters correspond to
      sequencing batches/centers

    Stijn van Dongen, Graph Clustering by Flow Simulation,
    PhD thesis, University of Utrecht, May 2000.
    https://micans.org/mcl/
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
    from scipy.cluster.hierarchy import fcluster, inconsistent, linkage
    from scipy.spatial.distance import pdist

    n_samples = len(samples)
    if n_samples < 2:
        return [np.array(range(n_samples))]

    dists = pdist(samples, metric="correlation")  # 1 - Pearson r
    # Replace NaN distances (from constant-coverage bins) with max distance
    dists = np.nan_to_num(dists, nan=1.0)
    Z = linkage(dists, method="average")  # UPGMA

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


def kmeans(samples, k=None):
    """Cluster samples using k-means on PCA-reduced data.

    Parameters
    ----------
    samples : 2D numpy array
        Matrix of samples' log2 values. Rows are samples, columns are bins.
    k : int, optional
        Number of clusters. If None, estimated from sample count.

    Returns
    -------
    list of lists
        Each inner list contains the indices of samples in that cluster.
    """
    from collections import defaultdict

    from scipy.cluster import vq

    if k is None:
        from math import log

        k = max(1, round(log(len(samples), 3)))

    logging.info("Clustering %d samples by k-means, where k = %d", len(samples), k)
    obs = pca_sk(samples, 3)
    obs = vq.whiten(obs)
    _centroids, labels = vq.kmeans2(obs, k, minit="++")
    clusters = defaultdict(list)
    for idx, label in enumerate(labels):
        clusters[label].append(idx)
    return list(clusters.values())


def markov(samples, inflation=5, max_iterations=100, by_pca=True):
    """Markov-cluster control samples by their read depths' correlation.

    Each of the matrices in the resulting iterable (list) can be processed the
    same as the input to calculate average log2 and spread values for that
    cluster.

    Parameters
    ----------
    samples : array
        Matrix of samples' read depths or normalized log2 values, as columns.
    inflation : float
        Inflation parameter for MCL. Must be >1; higher more granular clusters.
    by_pca : bool
        If true, similarity is by PCA; otherwise, by Pearson correlation.

    Return
    ------
    results : list
        A list of matrices representing non-overlapping column-subsets of the
        input, where each set of samples represents a cluster.
    """
    if inflation <= 1:
        raise ValueError("inflation must be > 1")

    if by_pca:
        pca_matrix = pca_sk(samples, 2)
        from scipy.spatial import distance

        dists = distance.squareform(distance.pdist(pca_matrix))
        M = 1 - (dists / dists.max())
    else:
        M = np.corrcoef(samples)

    M, clusters = mcl(M, max_iterations, inflation)
    return clusters


# https://github.com/koteth/python_mcl/blob/master/mcl/mcl_clustering.py
# https://github.com/GuyAllard/markov_clustering/blob/master/markov_clustering/mcl.py
# https://stackoverflow.com/questions/44243525/mcl-clustering-implementation-in-python-deal-with-overlap
def mcl(M, max_iterations, inflation, expansion=2):
    """Markov cluster algorithm."""
    logging.debug("MCL input matrix:\n%s", M)
    M = normalize(M)
    for i in range(max_iterations):
        M_prev = M
        M = inflate(expand(M, expansion), inflation)
        if converged(M, M_prev):
            logging.debug("Converged at iteration %d", i)
            break
        M = prune(M)

    clusters = get_clusters(M)
    return M, clusters


def normalize(A):
    """Normalize matrix columns."""
    return A / A.sum(axis=0)


def inflate(A, inflation):
    """Apply cluster inflation with the given element-wise exponent.

    From the mcl manual:

    This value is the main handle for affecting cluster granularity.
    This parameter is the usually only one that may require tuning.

    By default it is set to 2.0 and this is a good way to start. If you want to
    explore cluster structure in graphs with MCL, vary this parameter to obtain
    clusterings at different levels of granularity.  It is usually chosen
    somewhere in the range [1.2-5.0]. -I 5.0 will tend to result in fine-grained
    clusterings, and -I 1.2 will tend to result in very coarse grained
    clusterings. A good set of starting values is 1.4, 2, 4, and 6.
    Your mileage will vary depending on the characteristics of your data.

    Low values for -I, like -I 1.2, will use more CPU/RAM resources.

    Use mcl's cluster validation tools 'clm dist' and 'clm info' to test the
    quality and coherency of your clusterings.
    """
    return normalize(np.power(A, inflation))


def expand(A, expansion):
    """Apply cluster expansion with the given matrix power."""
    return np.linalg.matrix_power(A, expansion)


def converged(M, M_prev):
    """Test convergence.

    Criterion: homogeneity(??) or no change from previous round.
    """
    return np.allclose(M, M_prev)


# https://stackoverflow.com/questions/17772506/markov-clustering
def get_clusters(M):
    """Extract clusters from the matrix.

    Interpretation: "Attractors" are the non-zero elements of the matrix
    diagonal. The nodes in the same row as each attractor form a cluster.

    Overlapping clusterings produced by MCL are extremely rare, and always a
    result of symmetry in the input graph.

    Returns
    -------
    result : list
        A list of arrays of sample indices. The indices in each list item
        indicate the elements of that cluster; the length of the list is the
        number of clusters.
    """
    attractors_idx = M.diagonal().nonzero()[0]
    clusters_idx = [M[idx].nonzero()[0] for idx in attractors_idx]
    return clusters_idx


# https://github.com/GuyAllard/markov_clustering/blob/master/markov_clustering/mcl.py
def prune(M, threshold=0.001):
    """Remove many small entries while retaining most of M's stochastic mass.

    After pruning, vectors are rescaled to be stochastic again.
    (stochastic: values are all non-negative and sum to 1.)

    This step is purely to keep computation tractable in mcl by making the
    matrix more sparse (i.e. full of zeros), enabling sparse-matrix tricks to
    work.

    ----

    mcl:
        The default setting is something like -P 4000 -S 500 -R 600, where:

      -P <int> (1/cutoff)
      -S <int> (selection number)
      -R <int> (recover number)
      ---
      -pct <pct> (recover percentage)
      -p <num> (cutoff)

    After computing a new (column stochastic) matrix vector during expansion
    (which  is  matrix  multiplication c.q.  squaring), the vector is
    successively exposed to different pruning strategies. Pruning effectively
    perturbs the MCL process a little in order to obtain matrices that are
    genuinely sparse, thus keeping the computation tractable.

    mcl proceeds as follows:

    First, entries that are smaller than cutoff are
    removed, resulting in a vector with  at most 1/cutoff entries.

        * The cutoff can be supplied either by -p, or as the inverse value by
        -P.  The latter is more intuitive, if your intuition is like mine (P
        stands for precision or pruning).

    Second, if the remaining stochastic mass (i.e. the sum of all remaining
    entries) is less than <pct>/100 and the number of remaining entries is
    less than <r> (as specified by the -R flag), mcl will try to regain ground
    by recovering the largest discarded entries. If recovery was not necessary,
    mcl tries to prune the vector further down to at most s entries (if
    applicable), as specified by the -S flag. If this results in a vector that
    satisfies the recovery condition then recovery is attempted, exactly as
    described above. The latter will not occur of course if <r> <= <s>.

    """
    pruned = M.copy()
    pruned[pruned < threshold] = 0
    return pruned


def pca_sk(data, n_components=None):
    """Principal component analysis using scikit-learn.

    Parameters
    ----------
    data : 2D NumPy array
    n_components : int

    Returns: PCA-transformed data with `n_components` columns.
    """
    from sklearn.decomposition import PCA

    return PCA(n_components=n_components).fit_transform(data)


def pca_plain(data, n_components=None):
    """Principal component analysis using numpy eigenvalues.

    Source:
    https://stackoverflow.com/questions/13224362/principal-component-analysis-pca-in-python

    Parameters
    ----------
    data : 2D NumPy array
    n_components : int

    Returns: PCA-transformed data with `n_components` columns.
    """
    # Standarize the data
    data = data.copy()
    data -= data.mean(axis=0)
    data /= data.std(axis=0)
    # Covariance matrix
    C = np.cov(data)
    # Eigenvectors & eigenvalues of covariance matrix
    # ('eigh' rather than 'eig' since C is symmetric, for performance)
    E, V = np.linalg.eigh(C)
    # Sort eigenvalues in decreasing order, & eigenvectors to match,
    # keeping only the first `n_components` components
    key = np.argsort(E)[::-1][:n_components]
    E, V = E[key], V[:, key]
    # Transformation the data using eigenvectors
    U = np.dot(data, V)  # or: np.dot(V.T, data.T).T
    # U = np.dot(V.T, data.T).T
    # Return the re-scaled data, eigenvalues, and eigenvectors
    return U  # , E, V


def plot_clusters(M, cluster_indices) -> None:
    """Scatter plot first 2 components, colorized by cluster.

    Parameters
    ----------
    M : np.array
        PCA'd matrix. Rows are samples.
    cluster_indices : iterable of np.array
        Indices of samples in each cluster, as present in M.
    """
    from matplotlib import pyplot as plt

    _fig, ax = plt.subplots(1, 1)
    for cl_idx in cluster_indices:
        ax.scatter(M[cl_idx, 0], M[cl_idx, 1])
    plt.show()
