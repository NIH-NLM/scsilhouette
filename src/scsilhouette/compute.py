"""
Core computation functions for silhouette analysis.

This module contains the main function for computing silhouette scores
from single-cell data.
"""
# src/scsilhouette/compute.py
import json
import sys
from pathlib import Path
from typing import Optional

import anndata
import numpy as np
import pandas as pd
import scanpy as sc

from .exact_silhouette import silhouette_rows
from .utils import map_gene_symbols_to_ensembl, get_output_prefix
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
    journal: str,
    year: str,
    dataset_version_id: str,
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
    journal: str
        Journal (e.g., Cell)
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
    if not filter_normal and not use_binary_genes:
        # Only obs and the embeddings are used on this path, so the (large) expression
        # matrix X is never read: open backed, keep obs and obsm in memory, close the file.
        adata_backed = anndata.read_h5ad(h5ad_path, backed="r")
        adata = anndata.AnnData(
            obs=adata_backed.obs.copy(),
            obsm={k: np.asarray(v) for k, v in adata_backed.obsm.items()},
        )
        n_vars_file = adata_backed.n_vars
        adata_backed.file.close()
        logger.info(f"Loaded obs and obsm only (X not read): {adata.n_obs:,} cells x {n_vars_file:,} genes in file")
    else:
        adata = sc.read(h5ad_path)
        logger.info(f"Loaded AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    prefix = get_output_prefix(organ, first_author, journal, year, cluster_header, embedding_key, dataset_version_id)
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
    logger.info(f"Available obsm keys: {list(adata.obsm.keys())}")
    if embedding_key in adata.obsm:
        resolved_key = embedding_key
    elif f'X_{embedding_key}' in adata.obsm:
        resolved_key = f'X_{embedding_key}'
    else:
        matches = [k for k in adata.obsm.keys() if embedding_key in k]
        if len(matches) == 1:
            resolved_key = matches[0]
        elif len(matches) > 1:
            raise KeyError(f"Ambiguous: {matches}")
        else:
            raise KeyError(f"Not found. Available: {list(adata.obsm.keys())}")

    logger.info(f"Embedding resolved: '{embedding_key}' resolves to '{resolved_key}'")
    adata_use = adata.obsm[resolved_key]

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

    write_silhouette_outputs(
        obs=adata.obs,
        obs_keys=list(adata.obs.columns),
        obsm_keys=list(adata.obsm.keys()),
        cluster_header=cluster_header,
        embedding_key=embedding_key,
        prefix=prefix,
        organism=organism,
        organ=organ,
        disease=disease,
        filter_normal=filter_normal,
        uberon_json=uberon_json,
        disease_json=disease_json,
        hsapdv_json=hsapdv_json,
        save_scores=save_scores,
        save_cluster_summary=save_cluster_summary,
        save_annotation=save_annotation,
    )

    logger.info("Silhouette analysis complete!")
    return adata


def write_silhouette_outputs(
    obs,
    obs_keys,
    obsm_keys,
    cluster_header,
    embedding_key,
    prefix,
    organism="human",
    organ="",
    disease="normal",
    filter_normal=False,
    uberon_json=None,
    disease_json=None,
    hsapdv_json=None,
    save_scores=True,
    save_cluster_summary=True,
    save_annotation=True,
):
    """
    Write the per-cell scores, the cluster summary and the annotation JSON.

    ``obs`` is a DataFrame indexed by cell with the cluster column and ``silhouette_score``.
    Shared by run_silhouette and summarize_scores so both write identical files.
    """
    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    if save_scores:
        silhouette_scores_csv  = f"silhouette_scores_{prefix}.csv"
        silhouette_scores_json = f"silhouette_scores_{prefix}.json"
        obs[[cluster_header, "silhouette_score"]].to_csv(silhouette_scores_csv)
        obs[[cluster_header, "silhouette_score"]].to_json(silhouette_scores_json)
        logger.info(f"Saved silhouette scores to {silhouette_scores_csv}")

    if save_cluster_summary or save_annotation:
        logger.info("Computing cluster summary statistics...")
        cluster_summary_df = (
            obs.groupby(cluster_header)
            .agg(
                mean_silhouette=("silhouette_score", "mean"),
                std_silhouette=("silhouette_score", "std"),
                median_silhouette=("silhouette_score", "median"),
                count=("silhouette_score", "count"),
            )
            .reset_index()
        )

    if save_cluster_summary:
        cluster_summary_csv  = f"scsilhouette_cluster_summary_{prefix}.csv"
        cluster_summary_json = f"scsilhouette_cluster_summary_{prefix}.json"
        cluster_summary_df.to_csv(cluster_summary_csv, index=False)
        cluster_summary_df.to_json(cluster_summary_json)
        logger.info(f"Saved cluster summary to {cluster_summary_csv}")

    if save_annotation:
        annotation_json   = f"annotation_{prefix}.json"
        annotation_output = {
            "available_obs_keys":  list(obs_keys),
            "available_obsm_keys": list(obsm_keys),
            "cluster_header":      cluster_header,
            "embedding_key":       embedding_key,
            "dataset_summary": {
                "median_of_medians": float(cluster_summary_df['median_silhouette'].median()),
                "median_of_means":   float(cluster_summary_df['mean_silhouette'].median()),
                "median_of_stds":    float(cluster_summary_df['std_silhouette'].median()),
                "n_cells":           int(len(obs)),
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


def extract_silhouette_inputs(h5ad_path, cluster_header, embedding_key, organ, first_author, journal, year,
                              dataset_version_id):
    """
    Write the small files that silhouette scoring needs, reading only obs and the embedding (not X):
    the embedding (.npy), the cluster codes (.npy), the cell table (.csv) and an info file (.json)
    that holds the number of cells. Cells with NaN in the embedding are dropped, as in run_silhouette.
    """
    prefix = get_output_prefix(organ, first_author, journal, year, cluster_header, embedding_key, dataset_version_id)

    logger.info(f"Reading obs and obsm only (X is not read): {h5ad_path}")
    backed = anndata.read_h5ad(h5ad_path, backed="r")
    obs = backed.obs.copy()
    obsm_keys = list(backed.obsm.keys())
    if embedding_key in obsm_keys:
        resolved = embedding_key
    elif f"X_{embedding_key}" in obsm_keys:
        resolved = f"X_{embedding_key}"
    else:
        matches = [k for k in obsm_keys if embedding_key in k]
        if len(matches) != 1:
            raise KeyError(f"Embedding '{embedding_key}' not found or ambiguous. Available: {obsm_keys}")
        resolved = matches[0]
    embedding = np.asarray(backed.obsm[resolved])
    backed.file.close()
    logger.info(f"Embedding '{embedding_key}' resolves to '{resolved}': {embedding.shape}, {embedding.dtype}")
    if embedding.dtype != np.float32:
        logger.warning(
            f"The embedding is {embedding.dtype}, not float32: scores from silhouette-cells can differ from "
            f"compute-silhouette at about the 10th decimal (see exact_silhouette.py)"
        )

    if cluster_header not in obs.columns:
        raise ValueError(f"Cluster column '{cluster_header}' not found in obs: {list(obs.columns)}")

    nan_mask = np.isnan(embedding).any(axis=1)
    if nan_mask.any():
        logger.info(f"Dropped {int(nan_mask.sum()):,} cells with NaNs in '{embedding_key}'")
        obs = obs[~nan_mask]
        embedding = embedding[~nan_mask]

    from sklearn.preprocessing import LabelEncoder
    codes = LabelEncoder().fit_transform(obs[cluster_header].to_numpy()).astype(np.int32)

    np.save(f"silhouette_embedding_{prefix}.npy", embedding)
    np.save(f"silhouette_labels_{prefix}.npy", codes)
    obs[[cluster_header]].to_csv(f"silhouette_cells_{prefix}.csv")
    info = {
        "obs_keys": list(obs.columns),
        "obsm_keys": obsm_keys,
        "cluster_header": cluster_header,
        "embedding_key": embedding_key,
        "n_cells": int(embedding.shape[0]),
        "n_clusters": int(codes.max()) + 1,
        "embedding_dtype": str(embedding.dtype),     # the scores have this type, as in sklearn
    }
    with open(f"silhouette_info_{prefix}.json", "w") as f:
        json.dump(info, f, indent=2)
    logger.info(f"{info['n_cells']:,} cells, {info['n_clusters']} clusters")


def score_cells(embedding_npy, labels_npy, start, end, output, metric="euclidean"):
    """Exact silhouette of the cells start:end against all cells; the scores are written to ``output``."""
    X = np.load(embedding_npy)
    y = np.load(labels_npy)
    logger.info(f"Scoring cells {start:,}-{end:,} of {X.shape[0]:,} (metric: {metric})")
    scores = silhouette_rows(X, y, start, end, metric=metric)
    pd.DataFrame({"silhouette_score": scores}).to_csv(output, index=False)
    logger.info(f"Saved {output}")


def summarize_scores(cells_csv, info_json, scores_csv, organ, first_author, journal, year, dataset_version_id,
                     organism="human", disease="normal", save_scores=True, save_cluster_summary=True,
                     save_annotation=True):
    """Write the per-cell scores, cluster summary and annotation from a finished table of scores (one per cell)."""
    with open(info_json) as f:
        info = json.load(f)
    cluster_header = info["cluster_header"]
    embedding_key = info["embedding_key"]
    prefix = get_output_prefix(organ, first_author, journal, year, cluster_header, embedding_key, dataset_version_id)

    obs = pd.read_csv(cells_csv, index_col=0, dtype=str)
    # the scores are read back in the type they were computed in, so the averages match compute-silhouette
    scores = pd.read_csv(scores_csv)["silhouette_score"].to_numpy(dtype=info.get("embedding_dtype", "float64"))
    if len(scores) != len(obs) or len(obs) != info["n_cells"]:
        raise ValueError(f"{len(scores)} scores for {len(obs)} cells (expected {info['n_cells']})")
    obs["silhouette_score"] = scores

    write_silhouette_outputs(
        obs=obs,
        obs_keys=info["obs_keys"] + ["silhouette_score"],
        obsm_keys=info["obsm_keys"],
        cluster_header=cluster_header,
        embedding_key=embedding_key,
        prefix=prefix,
        organism=organism,
        organ=organ,
        disease=disease,
        filter_normal=False,
        save_scores=save_scores,
        save_cluster_summary=save_cluster_summary,
        save_annotation=save_annotation,
    )


def compute_summary_stats(
    cluster_summary_path: str,
    nsforest_results_path: str = None,
    metadata_path: str = None,
    cluster_header: str = "",
    organ: str = "",
    first_author: str = "",
    journal: str = "",
    year: str = "",
    embedding: str = "",
):
    """
    Compute dataset-level summary statistics from cluster summaries.

    Computes median-of-medians and other aggregate metrics across all clusters.
    Optionally includes NSForest F-score statistics if results are available.

    When metadata_path is provided, reads a JSON file containing all dataset
    fields from the cellxgene-harvester CSV. All metadata fields are included
    as columns in the output summary CSV. JSON fields override individual params.
    """
    # Load metadata from JSON if provided
    metadata = {}
    if metadata_path:
        with open(metadata_path) as f:
            metadata = json.load(f)

    # JSON fields override individual params
    organ = metadata.get('organ', organ)
    first_author = metadata.get('first_author', first_author)
    journal = metadata.get('journal', journal)
    year = metadata.get('year', year)
    cluster_header = metadata.get('author_cell_type', cluster_header)
    embedding = metadata.get('embedding', embedding)
    dataset_version_id = metadata.get('dataset_version_id', '')

    cluster_summary = pd.read_csv(cluster_summary_path)

    # Compute median of medians (primary QC metric)
    median_of_medians = cluster_summary['median_silhouette'].median()
    mean_of_medians = cluster_summary['median_silhouette'].mean()
    std_of_medians = cluster_summary['median_silhouette'].std()

    # Count clusters by quality
    high_quality = (cluster_summary['median_silhouette'] >= 0.5).sum()
    medium_quality = ((cluster_summary['median_silhouette'] >= 0.25) &
                      (cluster_summary['median_silhouette'] < 0.5)).sum()
    low_quality = (cluster_summary['median_silhouette'] < 0.25).sum()

    # Read NSForest results if available
    median_fscore = None
    mean_fscore = None
    if nsforest_results_path and nsforest_results_path != "NO_FILE":
        try:
            nsforest_df = pd.read_csv(nsforest_results_path)
            median_fscore = nsforest_df['f_score'].median()
            mean_fscore = nsforest_df['f_score'].mean()
        except Exception:
            logger.warning(f"Could not read NSForest results: {nsforest_results_path}")

    # Build prefix with embedding and dataset_version_id suffix for uniqueness
    prefix = get_output_prefix(organ, first_author, journal, year, cluster_header, embedding, dataset_version_id)
    
    # Start with ALL metadata fields (harvester columns), then add computed fields
    summary_data = dict(metadata)
    
    # Normalize ontology summary fields: replace commas with pipe separator,
    # strip numeric commas (122,000 -> 122000). Pipe is bioinformatics standard for
    # multi-value fields and avoids CSV quoting issues.
    
    for _key in ['tissue_ontology_summary',
                 'assay_ontology_summary',
                 'cell_type_ontology_summary',
                 'disease_ontology_summary',
                 'sex_ontology_summary',
                 'development_stage_summary']:
        _val = summary_data.get(_key)
        
        if not _val or (isinstance(_val, float) and pd.isna(_val)):
            summary_data[_key] = 'NA'
        elif isinstance(_val, str) and _val.strip():
            _val = _val.replace('; ', ' | ').replace(', ', ' | ')
            _val = _val.replace(',', '')   # remove numeric commas: 122,000 -> 122000
            summary_data[_key] = _val
        else:
            summary_data[_key] = 'NA'

    summary_data.update({
        'dataset': f"{prefix}",
        'organ': organ,
        'first_author': first_author,
        'journal': journal,
        'year': year,
        'cluster_header': cluster_header,
        'embedding': embedding,
        'n_clusters': len(cluster_summary),
        'median_silhouette': median_of_medians,
        'mean_silhouette': mean_of_medians,
        'std_silhouette': std_of_medians,
        'filtered_cell_count': int(cluster_summary['count'].sum()),
        'high_quality_clusters': int(high_quality),
        'medium_quality_clusters': int(medium_quality),
        'low_quality_clusters': int(low_quality),
        'median_fscore': median_fscore,
        'mean_fscore': mean_fscore,
    })

    # Convert to single-row DataFrame
    summary = pd.DataFrame({k: [v] for k, v in summary_data.items()})

    # Save
    output_path = f"master_dataset_summary_{prefix}.csv"
    summary.to_csv(output_path, index=False)

    # Log
    logger.info(f"Dataset Summary: {organ}_{first_author}_{year}")
    logger.info(f"  Clusters: {len(cluster_summary)} total")
    logger.info(f"  Cells: {int(cluster_summary['count'].sum())}")
    logger.info(f"  Median of median silhouette scores: {median_of_medians:.3f}")
    logger.info(f"    High quality (>= 0.5): {int(high_quality)}")
    logger.info(f"    Medium quality (0.25-0.5): {int(medium_quality)}")
    logger.info(f"    Low quality (< 0.25): {int(low_quality)}")
    if median_fscore is not None:
        logger.info(f"  Median F-score: {median_fscore:.3f}")
        logger.info(f"  Mean F-score: {mean_fscore:.3f}")
    logger.info(f"  Columns in output: {len(summary.columns)}")

    return output_path

