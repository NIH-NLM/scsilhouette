# src/scsilhouette/compute.py
import os
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics import silhouette_samples
from sklearn.metrics.pairwise import cosine_distances
from . import viz, nsforest


def compute_silhouette_scores(adata, label_key, embedding_key, metric="euclidean"):
    X = adata.obsm[embedding_key]
    labels = adata.obs[label_key]
    if metric == "cosine":
        dist_matrix = cosine_distances(X)
        scores = silhouette_samples(dist_matrix, labels, metric="precomputed")
    else:
        scores = silhouette_samples(X, labels, metric=metric)
    return pd.DataFrame({"cell_id": adata.obs_names, label_key: labels, f"silhouette_score_{metric}": scores})


def run_silhouette(
    h5ad_path: str,
    label_keys: list[str],
    embedding_key: str,
    output_dir: str,
    show_obs: bool = False,
    save_scores: bool = False,
    save_cluster_summary: bool = False,
    save_csv: bool = False,
    save_plots: bool = False,
    qc_correlations: bool = False,
    nsforest_path: str = None,
):
    adata = sc.read_h5ad(h5ad_path)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if show_obs:
        obs_preview = adata.obs.head()
        obs_preview.to_csv(Path(output_dir) / "obs_preview.csv")
        print(f"[✓] Saved obs preview: {output_dir}/obs_preview.csv")

    for label_key in label_keys:
        for metric in ["euclidean", "cosine"]:
            cell_scores = compute_silhouette_scores(adata, label_key, embedding_key, metric=metric)

            if save_scores:
                score_path = Path(output_dir) / f"{label_key}_scores_{metric}.csv"
                cell_scores.to_csv(score_path, index=False)
                print(f"[✓] Saved: {score_path.name}")

            if save_cluster_summary:
                cluster_summary = (
                    cell_scores.groupby(label_key)[f"silhouette_score_{metric}"]
                    .agg(["mean", "median", "std", "count"])
                    .reset_index()
                )
                summary_path = Path(output_dir) / f"{label_key}_cluster_summary_{metric}.csv"
                cluster_summary.to_csv(summary_path, index=False)
                print(f"[✓] Saved: {summary_path.name}")

                if save_plots:
                    viz.plot_all(cluster_summary, output_dir, label_key, suffix=metric)

            if qc_correlations:
                viz.plot_qc_boxplots(adata, cell_scores, output_dir, label_key, metric=metric)

            if nsforest_path:
                fscore_df = nsforest.read_nsforest_csv(nsforest_path)
                silhouette_df = cell_scores[[label_key, f"silhouette_score_{metric}"]].rename(
                    columns={f"silhouette_score_{metric}": "silhouette_score"}
                )
                viz.plot_fscore_vs_silhouette(
                    fscore_df, silhouette_df, output_dir, label_key, score_col="silhouette_score"
                )

