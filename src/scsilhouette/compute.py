# src/scsilhouette/compute.py

from pathlib import Path
from typing import List, Optional
import pandas as pd
import anndata as ad
from sklearn.metrics import silhouette_samples
from . import viz


def run_silhouette(
    h5ad_path: str,
    label_keys: List[str],
    embedding_key: str,
    output_dir: str,
    show_obs: bool = False,
    save_scores: bool = False,
    save_cluster_summary: bool = False,
    save_csv: bool = False,
    save_plots: bool = False,
    qc_correlations: bool = False,
    log_pca_dims: bool = False,
) -> None:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(h5ad_path)

    if show_obs:
        obs_preview = adata.obs.head(10)
        obs_preview.to_csv(output_path / "obs_preview.csv")
        print("[✓] Saved obs preview:", output_path / "obs_preview.csv")

    if log_pca_dims:
        dims = adata.obsm[embedding_key].shape[1]
        dim_path = output_path / "embedding_dims.txt"
        with open(dim_path, "w") as f:
            f.write(f"Embedding: {embedding_key}, PCA dimensions: {dims}\n")
        print("[✓] Saved PCA dimension info:", dim_path)

    for label_key in label_keys:
        embedding = adata.obsm[embedding_key]
        labels = adata.obs[label_key].values

        silhouette_vals = silhouette_samples(embedding, labels)
        cell_scores = adata.obs[[label_key]].copy()
        cell_scores["silhouette_score"] = silhouette_vals

        if save_scores:
            score_path = output_path / f"{label_key}_scores.csv"
            cell_scores.to_csv(score_path, index=False)
            print(f"[✓] Saved: {score_path.name}")

        cluster_summary = (
            cell_scores.groupby(label_key)["silhouette_score"]
            .agg(["mean", "median", "count"])
            .reset_index()
            .rename(
                columns={
                    "mean": "mean_silhouette_score",
                    "median": "median_silhouette_score",
                    "count": "n_cells",
                }
            )
        )

        if save_cluster_summary:
            cluster_path = output_path / f"{label_key}_cluster_summary.csv"
            cluster_summary.to_csv(cluster_path, index=False)
            print(f"[✓] Saved: {cluster_path.name}")

        if save_plots:
            viz.plot_score_distribution(cell_scores, output_path, label_key)
            viz.plot_cluster_summary(cluster_summary, output_path, label_key)
            viz.plot_cluster_size_vs_score(cluster_summary, output_path, label_key)

        if qc_correlations:
            viz.plot_qc_boxplots(cell_scores, adata.obs, output_path, label_key)

        if save_csv:
            cluster_summary.to_csv(output_path / f"{label_key}_summary.csv", index=False)

        viz.plot_all(cell_scores, cluster_summary, output_path, label_key)


