"""
Core computation functions for silhouette analysis.

This module contains the main function for computing silhouette scores
from single-cell data.
"""
# src/scsilhouette/compute.py
import json
import os
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc

from .utils import map_gene_symbols_to_ensembl
from .logging_config import setup_logger

logger = setup_logger()

# Terms that definitively exclude a cell from adult filtering.
# Copied verbatim from cellxgene-harvester count_normal_cells.py.
_EXCLUDE_STAGE_TERMS = [
    'fetal', 'embryo', 'newborn', 'prenatal', 'lmp',
    'post-fertilization', 'week post', 'Carnegie stage',
    'trimester', 'gestational'
]


# =============================================================================
# Shared helpers (identical logic to cellxgene-harvester)
# =============================================================================

def load_uberon_obo_ids(uberon_json: str) -> set:
    """
    Load UBERON obo_ids from a resolved JSON file produced by
    cellxgene-harvester resolve-uberon.

    Args:
        uberon_json: Path to uberon_<tissue>.json

    Returns:
        Set of UBERON obo_id strings e.g. {'UBERON:0002113', ...}
    """
    with open(uberon_json) as f:
        data = json.load(f)
    obo_ids = set(data["obo_ids"])
    roots   = [t["label"] for t in data["root_terms"]]
    logger.info(f"  Loaded UBERON JSON : {uberon_json}")
    logger.info(f"  Root terms         : {', '.join(roots)}")
    logger.info(f"  Total obo_ids      : {len(obo_ids):,}")
    return obo_ids


def extract_age_from_stage(stage_label: str) -> Optional[int]:
    """
    Parse numeric age from a development_stage label string.
    e.g. '25-year-old stage' → 25, '3rd decade' → None

    Copied verbatim from cellxgene-harvester count_normal_cells.py.
    """
    if not stage_label or not isinstance(stage_label, str):
        return None
    match = re.search(r'(\d+)[- ]?(?:year|yr)', stage_label.lower())
    return int(match.group(1)) if match else None


def filter_adult_cells(adata, min_age: int):
    """
    Filter cells to adults only using development_stage text parsing.

    Logic (same as cellxgene-harvester count_normal_cells.py):
      - Exclude cells whose stage matches any _EXCLUDE_STAGE_TERMS
      - Include cells whose stage contains 'adult'
      - Include cells where parsed numeric age >= min_age
      - Exclude everything else (unparseable = conservative exclusion)

    Args:
        adata:   AnnData object
        min_age: Minimum age in years (default 15)

    Returns:
        Filtered AnnData copy, labels Series filtered to match
    """
    if 'development_stage' not in adata.obs.columns:
        logger.warning("  'development_stage' not found in obs — skipping age filter")
        return adata

    if min_age == 0:
        logger.info("  min_age=0 — skipping age filter")
        return adata

    adult_mask           = []
    cells_adult          = 0
    cells_child          = 0
    cells_excluded_fetal = 0
    cells_unparseable    = 0

    for stage_val in adata.obs['development_stage'].astype(str):
        stage_lower = stage_val.lower()

        if not stage_val or stage_val in ('nan', 'None') or stage_val.strip() == '':
            adult_mask.append(False)
            cells_unparseable += 1
            continue

        if any(term in stage_lower for term in _EXCLUDE_STAGE_TERMS):
            adult_mask.append(False)
            cells_excluded_fetal += 1
            continue

        if 'adult' in stage_lower:
            adult_mask.append(True)
            cells_adult += 1
            continue

        age = extract_age_from_stage(stage_val)
        if age is not None:
            if age >= min_age:
                adult_mask.append(True)
                cells_adult += 1
            else:
                adult_mask.append(False)
                cells_child += 1
        else:
            adult_mask.append(False)
            cells_unparseable += 1

    n_before       = adata.n_obs
    adata_filtered = adata[adult_mask].copy()
    n_after        = adata_filtered.n_obs

    logger.info(f"  Age filter (>= {min_age} years): {n_before:,} → {n_after:,} cells "
                f"({n_before - n_after:,} removed)")
    logger.info(f"    Adult      : {cells_adult:,}")
    logger.info(f"    Child      : {cells_child:,}")
    logger.info(f"    Fetal/emb  : {cells_excluded_fetal:,}")
    logger.info(f"    Unparseable: {cells_unparseable:,}")

    return adata_filtered


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
    min_age: int = 15,
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
        If True, apply tissue + disease + age filters (requires uberon_json)
    uberon_json : str
        Path to UBERON JSON from cellxgene-harvester resolve-uberon.
        Required when filter_normal=True.
    min_age : int
        Minimum donor age in years for adult cell filtering (default 15)
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

        # 1. Tissue filter — UBERON obo_id precision
        logger.info("\n[1/3] Tissue filter")
        if 'tissue_ontology_term_id' not in adata.obs.columns:
            raise ValueError(
                "'tissue_ontology_term_id' not found in adata.obs. "
                "This column is required for UBERON-based tissue filtering. "
                "Ensure the h5ad file originates from CellxGene."
            )
        obo_ids      = load_uberon_obo_ids(uberon_json)
        tissue_mask  = adata.obs['tissue_ontology_term_id'].isin(obo_ids)
        n_before     = adata.n_obs
        adata        = adata[tissue_mask].copy()
        logger.info(f"  Tissue filter (UBERON IDs): {n_before:,} → {adata.n_obs:,} cells "
                    f"({n_before - adata.n_obs:,} removed)")

        if adata.n_obs == 0:
            present_ids = adata.obs['tissue_ontology_term_id'].value_counts().head(10) \
                if n_before > 0 else pd.Series()
            logger.warning("  WARNING: 0 cells remain after tissue filter.")

        # 2. Disease filter — PATO:0000461 only
        logger.info("\n[2/3] Disease filter")
        if 'disease_ontology_term_id' not in adata.obs.columns:
            raise ValueError(
                "'disease_ontology_term_id' not found in adata.obs. "
                "This column is required for disease filtering. "
                "Ensure the h5ad file originates from CellxGene."
            )
        disease_mask = adata.obs['disease_ontology_term_id'] == 'PATO:0000461'
        n_before     = adata.n_obs
        adata        = adata[disease_mask].copy()
        logger.info(f"  Disease filter (PATO:0000461): {n_before:,} → {adata.n_obs:,} cells "
                    f"({n_before - adata.n_obs:,} removed)")

        # 3. Age filter — development_stage text parsing
        logger.info("\n[3/3] Age filter")
        adata = filter_adult_cells(adata, min_age)

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
                "min_age":           min_age,
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
            "min_age":           min_age,
        }])
        dataset_summary_df.to_csv(dataset_summary_csv, index=False)
        logger.info(f"Saved dataset summary to {dataset_summary_csv}")

    logger.info("Silhouette analysis complete!")
    return adata
