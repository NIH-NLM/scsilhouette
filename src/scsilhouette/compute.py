# src/scsilhouette/compute.py

import os
from pathlib import Path
import pandas as pd
import scanpy as sc
import json
from .utils import map_gene_symbols_to_ensembl

def run_silhouette(
    h5ad_path: Path,
    label_key: str,
    embedding_key: str = None,
    organism: str = None,
    disease: str = None,
    tissue: str = None,
    cell_count: int = None,
    use_binary_genes: bool = False,
    gene_list_path: Path = None,
    mapping_output_path: Path = None,
    metric: str = "euclidean",
    pca_components: int = None,
    filter_normal: bool = True,
    save_scores: bool = False,
    save_cluster_summary: bool = False,
    save_annotation: bool = False,
    show_annotation: bool = False,
):
    from sklearn.metrics import silhouette_samples
    from sklearn.decomposition import PCA
    import numpy as np

    adata = sc.read(h5ad_path)
    prefix = Path(h5ad_path).stem

    print(f" will save output with {prefix}")
    print(f"[PASS] Loaded AnnData with {adata.shape}")

    # Handle gene subsetting if binary genes used
    if use_binary_genes and gene_list_path is not None:
        genes_df = pd.read_csv(gene_list_path)
        genes = genes_df.iloc[:, 0].tolist()
        print(f"[PASS] Loaded {len(genes)} binary genes from {gene_list_path}")

        if adata.raw is not None and all(g in adata.raw.var_names for g in genes):
            adata_use = adata.raw[:, genes]
            print("[PASS] Subsetting using raw and binary genes")
        else:
            adata_use = adata[:, genes]
            print("[WARN] Subsetting using var and binary genes")
    else:
        if embedding_key is None:
            adata_use = adata
        else:
            adata_use = adata.obsm[embedding_key]

    # PCA optionally
    if pca_components is not None:
        print(f"[INFO] Running PCA with {pca_components} components")
        pca = PCA(n_components=pca_components)
        adata_use = pca.fit_transform(adata_use)


    # Interpret flag string from CLI
    filter_normal = str(filter_normal).lower() == "true"

    # Get the labels first
    labels = adata.obs[label_key].copy()
    
    # Apply filter if needed
    if filter_normal and "disease" in adata.obs.columns:
        normal_mask = adata.obs["disease"] == "normal"
        adata_use = adata.obsm[embedding_key][normal_mask]
        labels = labels[normal_mask]
    else:
        adata_use = adata.obsm[embedding_key]

    scores = silhouette_samples(adata_use, labels, metric=metric)

    # Step 1: Init full column with NaNs
    adata.obs["silhouette_score"] = np.nan

    # Step 2: Assign only where used
    if filter_normal:
        adata.obs.loc[normal_mask, "silhouette_score"] = scores
    else:
        adata.obs["silhouette_score"] = scores


    if save_scores:
        silhouette_scores_csv = f"silhouette_scores_{organism}_{disease}_{tissue}_{prefix}_{label_key}_{embedding_key}_{metric}.csv"
        silhouette_scores_json = f"silhouette_scores_{organism}_{disease}_{tissue}_{prefix}_{label_key}_{embedding_key}_{metric}.json"
        adata.obs[[label_key, "silhouette_score"]].to_csv(silhouette_scores_csv)
        adata.obs[[label_key, "silhouette_score"]].to_json(silhouette_scores_json)
        print(f"[PASS] Saved silhouette scores to {silhouette_scores_csv}")
        print(f"[PASS] Saved silhouette scores to {silhouette_scores_json}")

    if save_cluster_summary:
        cluster_summary_df = (
            adata.obs.groupby(label_key)
            .agg(
                mean_silhouette=("silhouette_score", "mean"),
                std_silhouette=("silhouette_score", "std"),
                median_silhouette=("silhouette_score", "median"),
                count=("silhouette_score", "count"),
            )
            .reset_index()
        )
        cluster_summary_csv = f"cluster_summary_{organism}_{disease}_{tissue}_{prefix}_{label_key}_{embedding_key}_{metric}.csv"
        cluster_summary_json = f"cluster_summary_{organism}_{disease}_{tissue}_{prefix}_{label_key}_{embedding_key}_{metric}.json"
        cluster_summary_df.to_csv(cluster_summary_csv)
        cluster_summary_df.to_json(cluster_summary_json)
        print(f"[PASS] Saved cluster summary to {cluster_summary_csv}")
        print(f"[PASS] Saved cluster summary to {cluster_summary_json}")
    
    # Assume previous context where adata and other variables are defined
    if save_annotation:
        annotation_json = f"annotation_{organism}_{disease}_{tissue}_{prefix}_{label_key}_{embedding_key}_{metric}.json"
        annotation_output = {
            "available_obs_keys": list(adata.obs.columns),
            "available_obsm_keys": list(adata.obsm.keys()),
            "label_key": label_key,
            "embedding_key": embedding_key,
            "dataset_summary": {
                "median_of_medians": cluster_summary_df['median_silhouette'].median(),
                "median_of_means": cluster_summary_df['mean_silhouette'].median(),
                "median_of_stds": cluster_summary_df['std_silhouette'].median(),
                "n_cells": int(adata.n_obs),
                "n_clusters": int(cluster_summary_df.shape[0]),
                "organism": organism,
                "tissue": tissue,
                "disease": disease,
                "filter_normal": filter_normal
            }
        }
        with open(annotation_json, "w") as f:
            json.dump(annotation_output, f, indent=2)
            print(f"[PASS] Saved annotation metadata to {annotation_json}")

    return adata

