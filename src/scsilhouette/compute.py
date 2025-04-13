# src/scsilhouette/compute.py

import os
import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics import silhouette_samples
from typing import Optional, List
from .nsforest import load_nsforest_binary_genes


def compute_silhouette_scores(adata, label_key, embedding_key, metric="euclidean", gene_set: Optional[List[str]] = None):
    if gene_set:
        adata = adata[:, gene_set].copy()
        print(f"[i] Using {len(gene_set)} binary genes for silhouette score calculation")

    X = adata.obsm[embedding_key]
    labels = adata.obs[label_key].astype(str)

    silhouette_vals = silhouette_samples(X, labels, metric=metric)
    return pd.DataFrame({
        "cell_id": adata.obs_names,
        label_key: labels,
        f"silhouette_score_{metric}": silhouette_vals
    })


def run_silhouette(
    h5ad_path: str,
    label_keys: List[str],
    embedding_key: str,
    output_dir: str,
    nsforest_path: Optional[str] = None,
    metric: str = "euclidean",
    use_binary_genes: bool = False,
    save_scores: bool = False,
    save_cluster_summary: bool = False,
    save_csv: bool = False,
    qc_correlations: bool = False,
    show_obs: bool = False,
):
    os.makedirs(output_dir, exist_ok=True)
    adata = sc.read(h5ad_path)

    if show_obs:
        obs_preview = adata.obs.head()
        print("[✓] Saved obs preview:", os.path.join(output_dir, "obs_preview.csv"))
        obs_preview.to_csv(os.path.join(output_dir, "obs_preview.csv"))

    binary_gene_dict = {}
    fscore_df = None
    if nsforest_path:
        binary_gene_dict, fscore_df = load_nsforest_binary_genes(nsforest_path)
        print("[✓] Loaded NSForest file:", nsforest_path)

    for label_key in label_keys:
        gene_set = None
        if use_binary_genes and label_key in binary_gene_dict:
            gene_set = binary_gene_dict[label_key]
            gene_set = [g for g in gene_set if g in adata.var_names]
            print(f"[i] Using {len(gene_set)} binary genes for cluster '{label_key}'")

        cell_scores = compute_silhouette_scores(
            adata, label_key, embedding_key, metric=metric, gene_set=gene_set
        )

        cluster_summary = (
            cell_scores.groupby(label_key)[f"silhouette_score_{metric}"]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(columns={
                "mean": f"mean_silhouette_{metric}",
                "std": f"std_silhouette_{metric}",
                "count": "n_cells"
            })
        )

        suffix = "binary_genes_" + embedding_key if use_binary_genes else embedding_key

        if save_scores:
            cell_scores.to_csv(os.path.join(output_dir, f"{label_key}_silhouette_scores_{suffix}.csv"), index=False)
            print(f"[✓] Saved: {label_key}_silhouette_scores_{suffix}.csv")

        if save_cluster_summary:
            cluster_summary.to_csv(os.path.join(output_dir, f"{label_key}_cluster_summary_{metric}.csv"), index=False)
            print(f"[✓] Saved: {label_key}_cluster_summary_{metric}.csv")

        if qc_correlations and fscore_df is not None:
            fscore_out = os.path.join(output_dir, f"{label_key}_fscore_vs_silhouette_{metric}_corr.csv")
            print(f"[✓] QC Correlation saved to {fscore_out}")
