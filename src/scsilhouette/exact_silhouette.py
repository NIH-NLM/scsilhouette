"""
Exact silhouette scores for any range of cells, identical to sklearn.metrics.silhouette_samples.

sklearn scores all the cells in one call. To let a workflow score ranges of cells as separate tasks,
this module runs the same steps on the rows start:end only. The steps, the data types and the order of
the additions are those of sklearn, so each cell gets exactly the score sklearn gives it. This holds
bit for bit for a float32 embedding (what CellxGene stores). For a float64 embedding the distances come
from a matrix product whose result depends on how the rows are grouped, so a range can differ from
a whole-file run at about the 10th decimal.

_silhouette_reduce is the logic of sklearn.metrics.cluster._unsupervised._silhouette_reduce
(BSD 3-Clause, scikit-learn developers), copied here because that function is private.
"""
import numpy as np
from sklearn.metrics import pairwise_distances


def _silhouette_reduce(D_chunk, start, labels, label_freqs):
    """For a block of rows, the sum of the distances to the own cluster and the smallest mean distance to another."""
    n_chunk_samples = D_chunk.shape[0]
    cluster_distances = np.zeros((n_chunk_samples, len(label_freqs)), dtype=D_chunk.dtype)
    for i in range(n_chunk_samples):
        cluster_distances[i] += np.bincount(labels, weights=D_chunk[i], minlength=len(label_freqs))

    end = start + n_chunk_samples
    intra_index = (np.arange(n_chunk_samples), labels[start:end])
    intra_cluster_distances = cluster_distances[intra_index]
    cluster_distances[intra_index] = np.inf
    cluster_distances /= label_freqs
    inter_cluster_distances = cluster_distances.min(axis=1)
    return intra_cluster_distances, inter_cluster_distances


def silhouette_rows(X, y, start, end, metric="euclidean", block_memory_mb=256):
    """
    Silhouette of the cells start:end (original order) against all cells.

    Parameters
    ----------
    X : array (n_cells, n_features), all cells
    y : integer cluster codes (n_cells,), as sklearn's LabelEncoder gives
    start, end : the range of cells to score
    block_memory_mb : size of the distance block held at one time (rows x n_cells x 8 bytes)
    """
    X = np.asarray(X)
    y = np.asarray(y)
    n_cells = len(y)
    label_freqs = np.bincount(y)
    if len(label_freqs) < 2 or len(label_freqs) > n_cells - 1:
        raise ValueError("Number of labels must be between 2 and n_samples - 1")

    step = max(1, int(block_memory_mb * 1024 * 1024 // (n_cells * 8)))
    parts = []
    for lo in range(start, end, step):
        hi = min(lo + step, end)
        D = pairwise_distances(X[lo:hi], X, metric=metric)
        intra, inter = _silhouette_reduce(D, lo, y, label_freqs)

        denom = (label_freqs - 1).take(y[lo:hi], mode="clip")
        with np.errstate(divide="ignore", invalid="ignore"):
            intra /= denom
            sil = inter - intra
            sil /= np.maximum(intra, inter)
        parts.append(np.nan_to_num(sil))                       # nan is a cluster of size 1: the score is 0
    return np.concatenate(parts) if parts else np.empty(0)
