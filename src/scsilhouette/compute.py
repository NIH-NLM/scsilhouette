"""
Core computation functions for silhouette analysis.

This module contains the main function for computing silhouette scores
from single-cell data.
"""
# src/scsilhouette/compute.py
import json
import os
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc

from .utils import map_gene_symbols_to_ensembl
from .logging_config import setup_logger

logger = setup_logger()


# =============================================================================
# Ontology JSON loaders (mirror cellxgene-harvester helpers)
# =============================================================================

def load_obo_ids(json_path: str, label: str) -> set:
    """Load obo_ids from a resolve_uberon / resolve_disease JSON file."""
    with open(json_path) as f:
        data = json.load(f)
    obo_ids = set(data["obo_ids"])
    roots   = [t["label"] for t in data["root_terms"]]
    logger.info(f"  Loaded {label} JSON : {json_path}")
    logger.info(f"  Root terms          : {', '.join(roots)}")
    logger.info(f"  Total obo_ids       : {len(obo_ids):,}")
    return obo_ids


def load_obo_labels(json_path: str) -> dict:
    """Load obo_id -> label mapping from any resolve JSON file."""
    with open(json_path) as f:
        data = json.load(f)
    return {t["obo_id"]: t["label"] for t in data["terms"] if t.get("obo_id")}


def filter_by_age(adata, hsapdv_json: str):
    """Filter cells using development_stage_ontology_term_id matched against HsapDv obo_ids.

    Age threshold is encoded in the JSON at resolve time
    (cellxgene-harvester resolve-hsapdv --min-age N).
    No numeric comparison — identical pattern to tissue and disease filters.
    """
    id_col = "development_stage_ontology_term_id"
    if id_col not in adata.obs.columns:
        logger.warning(f"  '{id_col}' not found in obs - skipping age filter")
        return adata

    obo_ids    = load_obo_ids(hsapdv_json, "HsapDv")
    obo_labels = load_obo_labels(hsapdv_json)
    mask       = adata.obs[id_col].isin(obo_ids)
    n_before   = adata.n_obs

    all_counts     = adata.obs[id_col].value_counts()
    kept_terms     = [(tid, cnt) for tid, cnt in all_counts.items() if tid in obo_ids]
    excluded_terms = [(tid, cnt) for tid, cnt in all_counts.items() if tid not in obo_ids]

    logger.info(f"  Kept     ({len(kept_terms)} terms, "
                f"{sum(c for _, c in kept_terms):,} cells):")
    for term_id, count in kept_terms:
        label = obo_labels.get(term_id, "unknown label")
        logger.info(f"    KEPT     {term_id}  {label}: {count:,} cells")

    if excluded_terms:
        logger.info(f"  Excluded ({len(excluded_terms)} terms, "
                    f"{sum(c for _, c in excluded_terms):,} cells):")
        for term_id, count in excluded_terms:
            label = obo_labels.get(term_id, "unknown label — not in HsapDv JSON")
            logger.info(f"    EXCLUDED {term_id}  {label}: {count:,} cells")
    else:
        logger.info("  No cells excluded by age filter")

    adata = adata[mask].copy()
    logger.info(f"  Age filter (HsapDv IDs): {n_before:,} -> {adata.n_obs:,} cells "
                f"({n_before - adata.n_obs:,} removed)")
    return adata


# =============================================================================
# Main silhouette function
# =============================================================================

def run_silhouette(
    h5ad_path: str,
    cluster_header: str,
    embedding_key: str,
    organ: str,
    first_author: str,
    year: str,
    organism: str = "human",
    disease: str = "normal",
    use_binary_genes: bool = False,
    gene_list_path: str = None,
    metric: str = "euclidean",
    pca_components: int = None,
    filter_normal: bool = True,
    uberon_json: str = None,
    disease_json: str = None,
    hsapdv_json: str = None,
    save_scores: bool = True,
    save_cluster_summary: bool = True,
    save_annotation: bool = True,
):
    """
    Compute silhouette scores for single-cell clusters.

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
        Disease state label for annotation output (default: normal)
    use_binary_genes : bool
        Use binary genes from NSForest
    gene_list_path : str
        Path to gene list file
    metric : str
        Distance metric for silhouette (default: euclidean)
    pca_components : int
        Number of PCA components (optional)
    filter_normal : bool
        If True, apply tissue + disease + age filters
    uberon_json : str
        Path to UBERON JSON from cellxgene-harvester resolve-uberon.
        Required when filter_normal=True.
    disease_json : str
        Path to disease JSON from cellxgene-harvester resolve-disease.
        Required when filter_normal=True.
    hsapdv_json : str
        Path to HsapDv JSON from cellxgene-harvester resolve-hsapdv --min-age N.
        Age threshold is encoded in the JSON at resolve time.
        Required when filter_normal=True.
    save_scores : bool
        Save per-cell silhouette scores
    save_cluster_summary : bool
        Save cluster summary statistics
    save_annotation : bool
        Save annotation metadata

    Returns
    -------
    adata : AnnData
        Filtered, annotated data with silhouette scores
    """
    from sklearn.metrics import silhouette_samples
    from sklearn.decomposition import PCA

    # ------------------------------------------------------------------
    # Validate
    # ------------------------------------------------------------------
    filter_normal = bool(filter_normal)

    if filter_normal and not uberon_json:
        raise ValueError(
            "--uberon JSON file is required when --filter-normal is True. "
            "Generate it with: cellxgene-harvester resolve-uberon <tissue>"
        )
    if filter_normal and not disease_json:
        raise ValueError(
            "--disease JSON file is required when --filter-normal is True. "
            "Generate it with: cellxgene-harvester resolve-disease normal"
        )
    if filter_normal and not hsapdv_json:
        raise ValueError(
            "--hsapdv JSON file is required when --filter-normal is True. "
            "Generate it with: cellxgene-harvester resolve-hsapdv --min-age N"
        )

    # ------------------------------------------------------------------
    # Load
    # ------------------------------------------------------------------
    logger.info("Loading data...")
    adata = sc.read(h5ad_path)
    logger.info(f"Loaded AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    output_dir = f"outputs_{organ}_{first_author}_{year}"
    os.makedirs(output_dir, exist_ok=True)
    prefix = f"{output_dir}/{cluster_header}"
    logger.info(f"Output prefix: {prefix}")

    # ------------------------------------------------------------------
    # Optional gene subsetting (binary genes path — unchanged)
    # ------------------------------------------------------------------
    if use_binary_genes and gene_list_path is not None:
        genes_df = pd.read_csv(gene_list_path)
        genes    = genes_df.iloc[:, 0].tolist()
        logger.info(f"Loaded {len(genes)} binary genes from {gene_list_path}")
        if adata.raw is not None and all(g in adata.raw.var_names for g in genes):
            adata_use = adata.raw[:, genes]
            logger.info("Subsetting using raw and binary genes")
        else:
            adata_use = adata[:, genes]
            logger.info("Subsetting using var and binary genes")

    # ------------------------------------------------------------------
    # Filtering chain
    # ------------------------------------------------------------------
    if filter_normal:
        logger.info("\n=== Applying Filters ===")

        # 1. Tissue filter - UBERON obo_ids
        logger.info("\n[1/3] Tissue filter")
        if 'tissue_ontology_term_id' not in adata.obs.columns:
            raise ValueError(
                "'tissue_ontology_term_id' not found in adata.obs. "
                "This column is required for UBERON-based tissue filtering. "
                "Ensure the h5ad file originates from CellxGene."
            )
        obo_ids      = load_obo_ids(uberon_json, "UBERON")
        tissue_mask  = adata.obs['tissue_ontology_term_id'].isin(obo_ids)
        n_before     = adata.n_obs
        adata        = adata[tissue_mask].copy()
        logger.info(f"  Tissue filter (UBERON IDs): {n_before:,} -> {adata.n_obs:,} cells "
                    f"({n_before - adata.n_obs:,} removed)")

        if adata.n_obs == 0:
            logger.warning("  WARNING: 0 cells remain after tissue filter.")

        # 2. Disease filter - disease obo_ids
        logger.info("\n[2/3] Disease filter")
        if 'disease_ontology_term_id' not in adata.obs.columns:
            raise ValueError(
                "'disease_ontology_term_id' not found in adata.obs. "
                "This column is required for disease filtering. "
                "Ensure the h5ad file originates from CellxGene."
            )
        disease_ids  = load_obo_ids(disease_json, "disease")
        disease_mask = adata.obs['disease_ontology_term_id'].isin(disease_ids)
        n_before     = adata.n_obs
        adata        = adata[disease_mask].copy()
        logger.info(f"  Disease filter (ontology IDs): {n_before:,} -> {adata.n_obs:,} cells "
                    f"({n_before - adata.n_obs:,} removed)")

        # 3. Age filter - HsapDv ontology IDs
        logger.info("\n[3/3] Age filter")
        adata = filter_by_age(adata, hsapdv_json)

    logger.info(f"\nCells after filtering: {adata.n_obs:,}")

    # ------------------------------------------------------------------
    # Get cluster labels and embedding AFTER filtering
    # ------------------------------------------------------------------
    labels   = adata.obs[cluster_header].copy()
    adata_use = adata.obsm[embedding_key]

    # Optional PCA
    if pca_components is not None:
        logger.info(f"Running PCA with {pca_components} components")
        pca       = PCA(n_components=pca_components)
        adata_use = pca.fit_transform(adata_use)

    # Drop NaN rows in embedding
    if np.isnan(adata_use).any():
        nan_mask  = ~np.isnan(adata_use).any(axis=1)
        n_nan     = np.sum(~nan_mask)
        adata     = adata[nan_mask].copy()
        adata_use = adata_use[nan_mask]
        labels    = labels[nan_mask]
        logger.info(f"Dropped {n_nan:,} cells with NaNs in '{embedding_key}'")

    # ------------------------------------------------------------------
    # Silhouette scores
    # ------------------------------------------------------------------
    logger.info(f"Computing silhouette scores (metric: {metric})...")
    scores = silhouette_samples(adata_use, labels, metric=metric)

    # Assign scores directly — adata is already the filtered subset
    adata.obs["silhouette_score"] = scores

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    if save_scores:
        silhouette_scores_csv  = f"{prefix}_silhouette_scores.csv"
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
        cluster_summary_csv  = f"{prefix}_cluster_summary.csv"
        cluster_summary_json = f"{prefix}_cluster_summary.json"
        cluster_summary_df.to_csv(cluster_summary_csv, index=False)
        cluster_summary_df.to_json(cluster_summary_json)
        logger.info(f"Saved cluster summary to {cluster_summary_csv}")

    if save_annotation:
        annotation_json   = f"{prefix}_annotation.json"
        annotation_output = {
            "available_obs_keys":  list(adata.obs.columns),
            "available_obsm_keys": list(adata.obsm.keys()),
            "cluster_header":      cluster_header,
            "embedding_key":       embedding_key,
            "dataset_summary": {
                "median_of_medians": float(cluster_summary_df['median_silhouette'].median()),
                "median_of_means":   float(cluster_summary_df['mean_silhouette'].median()),
                "median_of_stds":    float(cluster_summary_df['std_silhouette'].median()),
                "n_cells":           int(adata.n_obs),
                "n_clusters":        int(cluster_summary_df.shape[0]),
                "organism":          organism,
                "organ":             organ,
                "disease":           disease,
                "filter_normal":     filter_normal,
                "uberon_json":       uberon_json,
                "disease_json":      disease_json,
                "hsapdv_json":       hsapdv_json,
            }
        }
        with open(annotation_json, "w") as f:
            json.dump(annotation_output, f, indent=2)
        logger.info(f"Saved annotation metadata to {annotation_json}")

        dataset_summary_csv = f"{prefix}_dataset_summary.csv"
        dataset_summary_df  = pd.DataFrame([{
            "organ":             organ,
            "first_author":      first_author,
            "year":              year,
            "cluster_header":    cluster_header,
            "n_cells":           int(adata.n_obs),
            "n_clusters":        int(cluster_summary_df.shape[0]),
            "median_of_medians": float(cluster_summary_df['median_silhouette'].median()),
            "median_of_means":   float(cluster_summary_df['mean_silhouette'].median()),
            "median_of_stds":    float(cluster_summary_df['std_silhouette'].median()),
            "organism":          organism,
            "disease":           disease,
            "filter_normal":     filter_normal,
            "uberon_json":       uberon_json,
            "disease_json":      disease_json,
            "hsapdv_json":       hsapdv_json,
        }])
        dataset_summary_df.to_csv(dataset_summary_csv, index=False)
        logger.info(f"Saved dataset summary to {dataset_summary_csv}")

    logger.info("Silhouette analysis complete!")
    return adata
