# src/scsilhouette/compute.py

import os
from pathlib import Path
import pandas as pd
import scanpy as sc
import json
from typing import Optional
from .utils import map_gene_symbols_to_ensembl
from .logging_config import setup_logger

logger = setup_logger()


def run_silhouette(
    h5ad_path: str,
    cluster_header: str,
    embedding_key: str,
    organ: str,
    first_author: str,
    year: str,
    organism: str = "human",
    disease: str = "normal",
    tissue: str = None,
    cell_count: str = "0",
    output_dir: str = None,
    use_binary_genes: bool = False,
    gene_list_path: str = None,
    metric: str = "euclidean",
    pca_components: int = None,
    filter_normal: bool = True,
    save_scores: bool = True,
    save_cluster_summary: bool = True,
    save_annotation: bool = True,
):
    """
    Compute silhouette scores for single-cell clusters
    
    Parameters
    ----------
    h5ad_path : str
        Path to input h5ad file
    cluster_header : str
        Column name for cell type clusters
    embedding_key : str
        Embedding key (e.g., X_umap)
    organ : str
        Organ/tissue (e.g., kidney)
    first_author : str
        First author (e.g., Lake)
    year : str
        Publication year (e.g., 2023)
    organism : str
        Organism (default: human)
    disease : str
        Disease state (default: normal)
    tissue : str
        Tissue type
    cell_count : str
        Cell count
    output_dir : str
        Output directory (auto-generated if None)
    use_binary_genes : bool
        Use binary genes from NSForest
    gene_list_path : str
        Path to gene list file
    metric : str
        Distance metric for silhouette (default: euclidean)
    pca_components : int
        Number of PCA components (optional)
    filter_normal : bool
        Filter to normal cells only
    save_scores : bool
        Save per-cell silhouette scores
    save_cluster_summary : bool
        Save cluster summary statistics
    save_annotation : bool
        Save annotation metadata
        
    Returns
    -------
    adata : AnnData
        Annotated data with silhouette scores
    """
    from sklearn.metrics import silhouette_samples
    from sklearn.decomposition import PCA
    import numpy as np

    logger.info("Loading data...")
    adata = sc.read(h5ad_path)
    
    stem = Path(h5ad_path).stem
    logger.info(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Use provided output_dir or auto-generate matching NSForest pattern
    if output_dir is None:
        output_dir = f"outputs_{organ}_{first_author}_{year}"
    
    os.makedirs(output_dir, exist_ok=True)
    
    # File prefix pattern matches NSForest: outputs_{organ}_{author}_{year}/{cluster_header}_
    prefix = f"{output_dir}/{cluster_header}"
    
    logger.info(f"Output prefix: {prefix}")
    
    # Handle gene subsetting if binary genes used
    if use_binary_genes and gene_list_path is not None:
        genes_df = pd.read_csv(gene_list_path)
        genes = genes_df.iloc[:, 0].tolist()
        logger.info(f"Loaded {len(genes)} binary genes from {gene_list_path}")

        if adata.raw is not None and all(g in adata.raw.var_names for g in genes):
            adata_use = adata.raw[:, genes]
            logger.info("Subsetting using raw and binary genes")
        else:
            adata_use = adata[:, genes]
            logger.info("Subsetting using var and binary genes")
    else:
        adata_use = adata.obsm[embedding_key] if embedding_key else adata

    # PCA optionally
    if pca_components is not None:
        logger.info(f"Running PCA with {pca_components} components")
        pca = PCA(n_components=pca_components)
        adata_use = pca.fit_transform(adata_use)

    # Interpret filter_normal flag
    filter_normal = str(filter_normal).lower() == "true"

    # Get cluster labels
    labels = adata.obs[cluster_header].copy()

    # Apply filter normal
    if filter_normal and "disease" in adata.obs.columns:
        normal_mask = adata.obs["disease"] == "normal"
        n_normal = sum(normal_mask)
        adata = adata[normal_mask].copy()
        labels = labels[normal_mask]
        logger.info(f"Filtered to {n_normal} normal cells")

    # Get embedding AFTER filtering
    adata_use = adata.obsm[embedding_key]

    # Drop rows with NaNs in embedding
    if np.isnan(adata_use).any():
        nan_mask = ~np.isnan(adata_use).any(axis=1)
        n_nan = np.sum(~nan_mask)
        adata = adata[nan_mask].copy()
        adata_use = adata_use[nan_mask]
        labels = labels[nan_mask]
        logger.info(f"Dropped {n_nan} cells with NaNs in '{embedding_key}'")

    logger.info(f"Computing silhouette scores (metric: {metric})...")
    scores = silhouette_samples(adata_use, labels, metric=metric)

    # Initialize silhouette_score column
    adata.obs["silhouette_score"] = np.nan

    # Assign scores
    if filter_normal:
        adata.obs.loc[normal_mask, "silhouette_score"] = scores
    else:
        adata.obs["silhouette_score"] = scores

    if save_scores:
        silhouette_scores_csv = f"{prefix}_silhouette_scores.csv"
        silhouette_scores_json = f"{prefix}_silhouette_scores.json"
        adata.obs[[cluster_header, "silhouette_score"]].to_csv(silhouette_scores_csv)
        adata.obs[[cluster_header, "silhouette_score"]].to_json(silhouette_scores_json)
        logger.info(f"Saved silhouette scores to {silhouette_scores_csv}")

    if save_cluster_summary:
        logger.info("Computing cluster summary statistics...")
        cluster_summary_df = (
            adata.obs.groupby(cluster_header)
            .agg(
                mean_silhouette=("silhouette_score", "mean"),
                std_silhouette=("silhouette_score", "std"),
                median_silhouette=("silhouette_score", "median"),
                count=("silhouette_score", "count"),
            )
            .reset_index()
        )
        cluster_summary_csv = f"{prefix}_cluster_summary.csv"
        cluster_summary_json = f"{prefix}_cluster_summary.json"
        cluster_summary_df.to_csv(cluster_summary_csv, index=False)
        cluster_summary_df.to_json(cluster_summary_json)
        logger.info(f"Saved cluster summary to {cluster_summary_csv}")
    
    if save_annotation:
        annotation_json = f"{prefix}_annotation.json"
        annotation_output = {
            "available_obs_keys": list(adata.obs.columns),
            "available_obsm_keys": list(adata.obsm.keys()),
            "cluster_header": cluster_header,
            "embedding_key": embedding_key,
            "dataset_summary": {
                "median_of_medians": float(cluster_summary_df['median_silhouette'].median()),
                "median_of_means": float(cluster_summary_df['mean_silhouette'].median()),
                "median_of_stds": float(cluster_summary_df['std_silhouette'].median()),
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
        logger.info(f"Saved annotation metadata to {annotation_json}")

    logger.info("Silhouette analysis complete!")
    return adata
