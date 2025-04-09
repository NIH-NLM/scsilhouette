
# src/scsilhouette/compute.py

import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics import silhouette_samples
from pathlib import Path
from . import viz


def compute_silhouette_scores(adata, label_key, embedding_key, metric="euclidean"):
    X = adata.obsm[embedding_key]
    labels = adata.obs[label_key]
    scores = silhouette_samples(X, labels, metric=metric)
    return pd.DataFrame({
        label_key: labels.values,
        f"silhouette_score_{metric}": scores
    })


def run_silhouette(
    h5ad_path,
    label_keys,
    embedding_key,
    output_dir,
    save_scores=False,
    save_cluster_summary=False,
    save_csv=False,
    save_plots=False,
    show_obs=False,
    qc_correlations=False,
    nsforest_path=None,
    use_binary_genes=False,
):
    import os
    import json

    adata = sc.read_h5ad(h5ad_path)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if show_obs:
        adata.obs.to_csv(Path(output_dir) / "obs_preview.csv")

    nsforest_df = None
    if nsforest_path:
        nsforest_df = pd.read_csv(nsforest_path)
        print(f"[✓] Loaded NSForest file: {nsforest_path}")

    for label_key in label_keys:
        for metric in ["euclidean", "cosine"]:
            suffix = f"{'binary_' if use_binary_genes else 'all'}_{embedding_key}_{metric}"

            # Subset to binary genes if applicable
            if use_binary_genes and nsforest_df is not None:
                genes = nsforest_df["binary_genes"].dropna().str.split(",").explode().str.strip().unique()
                common_genes = [g for g in genes if g in adata.raw.var_names]
                if not common_genes:
                    print(f"[!] No common binary genes found in dataset for {label_key}")
                    continue
                adata_subset = adata[:, common_genes].copy()
                X = adata.obsm[embedding_key]
            else:
                adata_subset = adata

            cell_scores = compute_silhouette_scores(adata_subset, label_key, embedding_key, metric=metric)

            if save_scores:
                cell_scores.to_csv(Path(output_dir) / f"{label_key}_{suffix}_scores.csv", index=False)

            cluster_summary = (
                cell_scores.groupby(label_key)
                .agg(
                    n_cells=(f"silhouette_score_{metric}", "count"),
                    mean_score=(f"silhouette_score_{metric}", "mean"),
                )
                .reset_index()
            )

            if save_cluster_summary:
                cluster_summary.to_csv(Path(output_dir) / f"{label_key}_{suffix}_cluster_summary.csv", index=False)

            if save_plots:
                viz.plot_all(
                    cluster_summary=cluster_summary,
                    cell_scores=cell_scores,
                    output_dir=output_dir,
                    label=label_key,
                    score_col=f"silhouette_score_{metric}",
                    suffix=suffix,
                    file_format="svg"
                )

            if qc_correlations:
                if nsforest_df is not None:
                    viz.plot_fscore_vs_silhouette(
                        nsforest_df,
                        cluster_summary.rename(columns={"mean_score": f"silhouette_score_{metric}"}),
                        output_dir,
                        label_key,
                        score_col=f"silhouette_score_{metric}",
                        suffix=suffix
                    )
