#!/usr/bin/env python
"""Markov-cluster control samples' correlation matrix.

See:

    Stijn van Dongen, Graph Clustering by Flow Simulation,
    PhD thesis, University of Utrecht, May 2000.
    https://micans.org/mcl/

"""

import logging

import numpy as np


def kmeans(samples, k=None):
    from scipy.cluster import vq

    if k is None:
        from math import log

        k = max(1, round(log(len(samples), 3)))
        # E.g. n=66 -> k=2, 16 -> 3, 47 -> 4, 141 -> 5, 421 -> 6, 1263 -> 7

    print("Clustering", len(samples), "samples by k-means, where k =", k)
    obs = pca_sk(samples, 3)
    obs = vq.whiten(obs)  # Needed?
    _centroids, labels = vq.kmeans2(obs, k, minit="++")
    # XXX shorter way to do this grouping?
    from collections import defaultdict

    clusters = defaultdict(list)
    for idx, label in enumerate(labels):
        clusters[label].append(idx)
    clusters = clusters.values()
    # plot_clusters(obs, clusters)
    return clusters


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
        pca_matrix = pca_sk(samples, 2)  # pca_plain
        # Convert to similarity matrix
        from scipy.spatial import distance

        dists = distance.squareform(distance.pdist(pca_matrix))
        M = 1 - (dists / dists.max())
    else:
        M = np.corrcoef(samples)

    M, clusters = mcl(M, max_iterations, inflation)
    # plot_clusters(M, clusters)
    return clusters


# https://github.com/koteth/python_mcl/blob/master/mcl/mcl_clustering.py
# https://github.com/GuyAllard/markov_clustering/blob/master/markov_clustering/mcl.py
# https://stackoverflow.com/questions/44243525/mcl-clustering-implementation-in-python-deal-with-overlap
def mcl(M, max_iterations, inflation, expansion=2):
    """Markov cluster algorithm."""
    print("M_init:\n", M)
    M = normalize(M)
    # print("M_norm:\n", M)
    for i in range(max_iterations):
        M_prev = M
        M = inflate(expand(M, expansion), inflation)
        # print("M_inflate_%d:\n" % i, M)
        if converged(M, M_prev):
            logging.debug("Converged at iteration %d", i)
            break
        M = prune(M)
        # print("M_prune_%d:\n" % i, M)

    # print("M_final:\n", M)
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


def plot_clusters(M, cluster_indices):
    """Scatter plot first 2 components, colorized by cluster.

    Parameters
    ----------
    M : np.array
        PCA'd matrix. Rows are samples.
    cluster_indices : iterable of np.array
        Indices of samples in each cluster, as present in M.
    """
    from matplotlib import pyplot as plt

    # colors = list("krgbo")[:len(cluster_indices)]
    _fig, ax = plt.subplots(1, 1)
    for cl_idx in cluster_indices:
        ax.scatter(M[cl_idx, 0], M[cl_idx, 1])  # c=color
    plt.show()
