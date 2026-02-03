# src/scsilhouette/cli.py
"""
Command-line interface for scsilhouette.

This module provides CLI commands for computing silhouette scores
and generating visualizations.
"""

import typer
import os
from pathlib import Path
from typing import Optional

from . import compute, viz
from .logging_config import setup_logger

logger = setup_logger()
app = typer.Typer(add_completion=False)


@app.command("compute-silhouette")
def compute_silhouette_command(
    h5ad_path: Path = typer.Option(..., help="Path to input h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for cell type clusters"),
    embedding_key: str = typer.Option(..., help="Embedding key (e.g., X_umap)"),
    organ: str = typer.Option(..., help="Organ/tissue (e.g., kidney)"),
    first_author: str = typer.Option(..., help="First author (e.g., Lake)"),
    year: str = typer.Option(..., help="Publication year (e.g., 2023)"),
    organism: str = typer.Option("human", help="Organism"),
    disease: str = typer.Option("normal", help="Disease state"),
    use_binary_genes: bool = typer.Option(False, help="Use binary genes from NSForest"),
    gene_list_path: Optional[Path] = typer.Option(None, help="Path to gene list"),
    metric: str = typer.Option("euclidean", help="Distance metric for silhouette"),
    filter_normal: str = typer.Option("True", help="Filter to normal cells only"),
    save_scores: bool = typer.Option(True, help="Save per-cell silhouette scores"),
    save_cluster_summary: bool = typer.Option(True, help="Save cluster summary statistics"),
    save_annotation: bool = typer.Option(True, help="Save annotation metadata"),
):
    """
    Compute silhouette scores for single-cell clusters.
    
    This command calculates silhouette scores to assess clustering quality
    in single-cell RNA-seq data.
    
    Args:
        h5ad_path: Path to input h5ad file
        cluster_header: Column name for cell type clusters
        embedding_key: Embedding key (e.g., X_umap)
        organ: Organ/tissue (e.g., kidney)
        first_author: First author (e.g., Lake)
        year: Publication year (e.g., 2023)
    
    Returns:
        None. Writes output files to outputs_{organ}_{first_author}_{year}/
    """
    
    
    logger.info("="*80)
    logger.info("Silhouette Score Analysis Pipeline")
    logger.info("="*80)
    logger.info(f"Input file: {h5ad_path}")
    logger.info(f"Dataset: {organ} / {first_author} / {year}")
    logger.info(f"Cluster header: {cluster_header}")
    logger.info(f"Embedding: {embedding_key}")
    logger.info("="*80)
    
    compute.run_silhouette(
        h5ad_path=str(h5ad_path),
        cluster_header=cluster_header,
        embedding_key=embedding_key,
        organ=organ,
        first_author=first_author,
        year=year,
        organism=organism,
        disease=disease,
        use_binary_genes=use_binary_genes,
        gene_list_path=str(gene_list_path) if gene_list_path else None,
        metric=metric,
        filter_normal=filter_normal,
        save_scores=save_scores,
        save_cluster_summary=save_cluster_summary,
        save_annotation=save_annotation,
    )


@app.command("viz-summary")
def viz_summary_command(
    silhouette_score_path: Path = typer.Option(..., help="Path to silhouette scores CSV"),
    silhouette_score_col: str = typer.Option("silhouette_score", help="Column name for scores"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    fscore_path: Optional[Path] = typer.Option(None, help="Path to NSForest results with F-scores"),
    sort_by: str = typer.Option("median", help="Sort by mean|median|std"),
):
    """Generate silhouette summary visualization with F-scores"""
    
    output_dir = f"outputs_{organ}_{first_author}_{year}"
    os.makedirs(output_dir, exist_ok=True)
    
    viz.plot_silhouette_summary(
        silhouette_score_path=str(silhouette_score_path),
        silhouette_score_col=silhouette_score_col,
        cluster_header=cluster_header,
        fscore_path=str(fscore_path) if fscore_path else None,
        sort_by=sort_by,
        output_dir=output_dir,
    )


@app.command("viz-dotplot")
def viz_dotplot_command(
    h5ad_path: Path = typer.Option(..., help="Path to input .h5ad file"),
    embedding_key: str = typer.Option(..., help="Embedding key (X_umap, X_tsne)"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Generate UMAP/embedding dotplot colored by cluster"""
    
    output_dir = f"outputs_{organ}_{first_author}_{year}"
    
    viz.plot_dotplot(
        h5ad_path=str(h5ad_path),
        cluster_header=cluster_header,
        embedding_key=embedding_key,
        organ=organ,
        first_author=first_author,
        year=year,
        output_dir=output_dir,
    )


@app.command("viz-distribution")
def viz_distribution_command(
    cluster_summary_path: Path = typer.Option(..., help="CSV with cluster summary"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Generate distribution plots of cluster sizes vs silhouette"""
    
    output_dir = f"outputs_{organ}_{first_author}_{year}"
    
    viz.plot_distribution(
        cluster_summary_path=str(cluster_summary_path),
        cluster_header=cluster_header,
        organ=organ,
        first_author=first_author,
        year=year,
        output_dir=output_dir,
    )


def main():
    app()
