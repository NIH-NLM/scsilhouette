"""
Exact silhouette scores for any range of cells.

Same definition as sklearn.metrics.silhouette_samples (exact, no sampling or approximation). The cells
are sorted by cluster once, so the sum of distances from a cell to every cluster is a single numpy
reduceat over the distances of a block of rows. Which cells to score is the caller's choice: all of
them, or a range, so that a workflow can score the ranges as separate tasks.
"""
import numpy as np
from sklearn.metrics import pairwise_distances


def _sorted_reference(X, y):
    """The cells sorted by cluster, with the size and first position of each cluster."""
    n = X.shape[0]
    k = len(np.unique(y))
    if k < 2 or k > n - 1:
        raise ValueError("Number of labels must be between 2 and n_samples - 1")
    order = np.argsort(y, kind="stable")
    freq = np.bincount(y[order]).astype(np.float64)
    starts = np.concatenate(([0], np.cumsum(freq[:-1]).astype(np.int64)))
    return X[order], freq, starts


def _score_block(X_rows, y_rows, X_sorted, freq, starts, metric):
    """Silhouette of a block of rows against all cells."""
    D = pairwise_distances(X_rows, X_sorted, metric=metric)                   # rows x n
    sums = np.add.reduceat(D, starts, axis=1, dtype=np.float64)               # rows x K: sum of distances to each cluster (float64, as sklearn)
    rows = np.arange(len(y_rows))
    own_freq = freq[y_rows]
    with np.errstate(divide="ignore", invalid="ignore"):
        intra = sums[rows, y_rows] / (own_freq - 1)                           # exclude the cell itself (distance 0)
        intra[own_freq == 1] = 0.0
        mean_to = sums / freq                                                 # rows x K
        mean_to[rows, y_rows] = np.inf                                        # leave out the own cluster
        inter = mean_to.min(axis=1)
        s = (inter - intra) / np.maximum(intra, inter)
    s[own_freq == 1] = 0.0                                                    # sklearn: singleton clusters score 0
    return np.nan_to_num(s, copy=False)


def silhouette_rows(X, y, start, end, metric="euclidean", block_memory_mb=256):
    """
    Exact silhouette of cells start:end (original order) against all cells.

    Parameters
    ----------
    X : array (n_cells, n_features), all cells
    y : integer cluster codes (n_cells,)
    start, end : the range of cells to score
    block_memory_mb : size of the distance block held at one time (rows x n_cells x 8 bytes)
    """
    X = np.asarray(X)
    y = np.asarray(y)
    X_sorted, freq, starts = _sorted_reference(X, y)
    step = max(1, int(block_memory_mb * 1024 * 1024 // (X.shape[0] * 8)))
    parts = [
        _score_block(X[lo:hi], y[lo:hi], X_sorted, freq, starts, metric)
        for lo, hi in ((lo, min(lo + step, end)) for lo in range(start, end, step))
    ]
    return np.concatenate(parts) if parts else np.empty(0)
