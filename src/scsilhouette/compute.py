# src/scsilhouette/compute.py

import os
import pickle
from typing import Optional, List
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_samples
import numpy as np


def load_binary_gene_list(path: str) -> List[str]:
    genes = []
    with open(path) as f:
        genes = [line.strip() for line in f if line.strip()]
    print(f"[PASS] Loaded {len(genes)} binary genes from file.")
    return genes


def validate_adata_for_save(adata: sc.AnnData):
    for attr in [adata.obs, adata.var, getattr(adata.raw, "var", pd.DataFrame())]:
        if "_index" in attr.columns:
            attr.drop(columns=["_index"], inplace=True)


def compute_silhouette(
    adata: sc.AnnData,
    label_key: str,
    embedding_key: str,
    gene_subset: Optional[List[str]] = None,
    metric: str = "euclidean"
) -> pd.DataFrame:

    if embedding_key not in adata.obsm:
        raise ValueError(f"{embedding_key} not found in adata.obsm")
    if label_key not in adata.obs:
        raise ValueError(f"{label_key} not found in adata.obs")

    if gene_subset:
        adata_subset = adata[:, adata.var_names.isin(gene_subset)].copy()
        X = adata_subset.X
        print(f"[PASS] silhouette_samples calculated on {len(gene_subset)} binary genes")
    else:
        X = adata.obsm[embedding_key]
        print(f"[PASS] silhouette_samples calculated on {embedding_key}")

    y = adata.obs[label_key].values
    scores = silhouette_samples(X, y, metric=metric)

    df = adata.obs[[label_key]].copy()
    df["silhouette_score_" + metric] = scores
    return df


def summarize_silhouette(df: pd.DataFrame, label_key: str, score_key: str) -> pd.DataFrame:
    summary = (
        df.groupby(label_key)[score_key]
        .agg(["mean", "std", "median", "count"])
        .reset_index()
    )
    print(f"[PASS] Generated summary for {summary.shape[0]} clusters.")
    return summary


def run_silhouette(
    h5ad_path: str,
    label_key: str,
    embedding_key: str,
    output_dir: str,
    use_binary_genes: bool = False,
    gene_list_path: Optional[str] = None,
    metric: str = "euclidean",
    save_scores: bool = False,
    save_cluster_summary: bool = False,
    save_csv: bool = False,
    show_obs: bool = False
):
    adata = sc.read(h5ad_path)

    if show_obs:
        print(adata.obs.head())

    gene_subset = load_binary_gene_list(gene_list_path) if use_binary_genes and gene_list_path else None

    df_scores = compute_silhouette(
        adata,
        label_key,
        embedding_key,
        gene_subset=gene_subset,
        metric=metric
    )

    os.makedirs(output_dir, exist_ok=True)

    if save_scores:
        pkl_path = os.path.join(output_dir, f"{label_key}_silhouette_scores.pkl")
        with open(pkl_path, "wb") as f:
            pickle.dump(df_scores, f)
        print(f"[PASS] Saved silhouette scores to {pkl_path}")

    if save_csv:
        csv_path = os.path.join(output_dir, f"{label_key}_silhouette_scores_binary_genes_{embedding_key}.csv")
        df_scores.to_csv(csv_path, index=False)
        print(f"[PASS] Saved silhouette scores CSV to {csv_path}")

    if save_cluster_summary:
        score_col = df_scores.columns[-1]
        summary = summarize_silhouette(df_scores, label_key, score_col)
        summary_path = os.path.join(output_dir, f"{label_key}_cluster_summary_{metric}.csv")
        summary.to_csv(summary_path, index=False)
        print(f"[PASS] Saved cluster summary to {summary_path}")


