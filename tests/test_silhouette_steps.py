"""
The silhouette steps (silhouette-inputs, silhouette-cells, silhouette-summary) must give exactly the
scores and files of compute-silhouette, so a workflow can run them as separate tasks.
"""
import glob

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from sklearn.metrics import silhouette_samples
from sklearn.preprocessing import LabelEncoder

from scsilhouette import compute
from scsilhouette.exact_silhouette import silhouette_rows

COMMON = dict(organ="neocortex", first_author="Jorstad", journal="Science", year="2023",
              dataset_version_id="abcdef123456")


def _embedding_and_labels(n=3000, k=7, seed=0):
    rng = np.random.default_rng(seed)
    labels = rng.integers(0, k, n)
    centers = rng.normal(scale=3, size=(k, 2))
    embedding = (centers[labels] + rng.normal(size=(n, 2))).astype("float32")
    return embedding, np.array([f"type {i}" for i in labels])


def test_rows_equal_sklearn_for_any_ranges():
    embedding, names = _embedding_and_labels()
    names[5] = "a cluster of one cell"            # sklearn scores it 0
    codes = LabelEncoder().fit_transform(names)
    expected = silhouette_samples(embedding, names)

    n = len(names)
    pieces = [silhouette_rows(embedding, codes, s, min(s + 700, n), block_memory_mb=1) for s in range(0, n, 700)]
    assert np.array_equal(np.concatenate(pieces), expected)


def test_steps_equal_compute_silhouette(tmp_path, monkeypatch):
    embedding, names = _embedding_and_labels(n=2000)
    embedding[17] = np.nan                         # a cell without a position is dropped by both
    obs = pd.DataFrame({"ct": names}, index=[f"cell{i}" for i in range(len(names))])
    obs.index.name = "sample_id"
    adata = ad.AnnData(sp.csr_matrix(np.ones((len(names), 5), dtype="float32")), obs=obs)
    adata.obsm["X_UMAP"] = embedding
    h5ad = tmp_path / "small.h5ad"
    adata.write_h5ad(h5ad)

    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir()
    new.mkdir()

    monkeypatch.chdir(old)
    compute.run_silhouette(h5ad_path=str(h5ad), cluster_header="ct", embedding_key="X_UMAP",
                           disease="normal", filter_normal=False, **COMMON)

    monkeypatch.chdir(new)
    compute.extract_silhouette_inputs(str(h5ad), "ct", "X_UMAP", **COMMON)
    emb, lab = glob.glob("silhouette_embedding_*.npy")[0], glob.glob("silhouette_labels_*.npy")[0]
    n_cells = len(np.load(lab))
    assert n_cells == len(names) - 1
    cut = n_cells // 3                             # the workflow decides the ranges; here three
    compute.score_cells(emb, lab, 0, cut, "scores_0.csv")
    compute.score_cells(emb, lab, cut, 2 * cut, "scores_1.csv")
    compute.score_cells(emb, lab, 2 * cut, n_cells, "scores_2.csv")
    parts = [pd.read_csv(f"scores_{i}.csv") for i in range(3)]
    pd.concat(parts, ignore_index=True).to_csv("all_scores.csv", index=False)
    compute.summarize_scores(glob.glob("silhouette_cells_*.csv")[0], glob.glob("silhouette_info_*.json")[0],
                             "all_scores.csv", disease="normal", **COMMON)

    for pattern in ("silhouette_scores_*.csv", "scsilhouette_cluster_summary_*.csv", "annotation_*.json"):
        a, b = glob.glob(f"{old}/{pattern}")[0], glob.glob(f"{new}/{pattern}")[0]
        assert open(a, "rb").read() == open(b, "rb").read(), pattern
