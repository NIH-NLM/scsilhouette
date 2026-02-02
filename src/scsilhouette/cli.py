# src/scsilhouette/cli.py

import typer
import os
from pathlib import Path
from typing import Optional, List

from . import compute, download, viz, nsforest
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
    tissue: str = typer.Option(..., help="Tissue type"),
    cell_count: str = typer.Option("0", help="Cell count"),
    output_dir: Optional[Path] = typer.Option(None, help="Output directory (auto-generated if not provided)"),
    use_binary_genes: bool = typer.Option(False, help="Use binary genes from NSForest"),
    gene_list_path: Optional[Path] = typer.Option(None, help="Path to gene list"),
    metric: str = typer.Option("euclidean", help="Distance metric for silhouette"),
    filter_normal: str = typer.Option("True", help="Filter to normal cells only"),
    save_scores: bool = typer.Option(True, help="Save per-cell silhouette scores"),
    save_cluster_summary: bool = typer.Option(True, help="Save cluster summary statistics"),
    save_annotation: bool = typer.Option(True, help="Save annotation metadata"),
):
    """Compute silhouette scores for single-cell clusters"""
    
    # Auto-generate output directory matching NSForest pattern
    if output_dir is None:
        output_dir = Path(f"outputs_{organ}_{first_author}_{year}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("="*80)
    logger.info("Silhouette Score Analysis Pipeline")
    logger.info("="*80)
    logger.info(f"Input file: {h5ad_path}")
    logger.info(f"Dataset: {organ} / {first_author} / {year}")
    logger.info(f"Cluster header: {cluster_header}")
    logger.info(f"Embedding: {embedding_key}")
    logger.info(f"Output directory: {output_dir}")
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
        tissue=tissue,
        cell_count=cell_count,
        output_dir=str(output_dir),
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
    output_dir: Optional[Path] = typer.Option(None, help="Output directory"),
    fscore_path: Optional[Path] = typer.Option(None, help="Path to NSForest results with F-scores"),
    sort_by: str = typer.Option("median", help="Sort by mean|median|std"),
):
    """Generate silhouette summary visualization with F-scores"""
    
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        os.chdir(str(output_dir))
    
    viz.plot_silhouette_summary(
        silhouette_score_path=str(silhouette_score_path),
        silhouette_score_col=silhouette_score_col,
        cluster_header=cluster_header,
        fscore_path=str(fscore_path) if fscore_path else None,
        sort_by=sort_by,
    )

    
@app.command("viz-correlation")
def viz_correlation_command(
    cluster_summary_path: Path = typer.Option(..., help="CSV with cluster summary metrics"),
    x_metric: str = typer.Option(..., help="X-axis metric (e.g., f_score)"),
    y_metrics: str = typer.Option(..., help="Comma-separated Y-axis metrics (e.g., mean_silhouette,median_silhouette)"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    fscore_path: Optional[Path] = typer.Option(None, help="Path to NSForest F-scores"),
    mapping_path: Optional[Path] = typer.Option(None, help="Path to cluster mapping file"),
):
    """Generate correlation plots between metrics"""
    
    viz.plot_correlation_summary(
        cluster_summary_path=str(cluster_summary_path),
        x_metric=x_metric,
        y_metrics=[x.strip() for x in y_metrics.split(",")],
        cluster_header=cluster_header,
        fscore_path=str(fscore_path) if fscore_path else None,
        mapping_path=str(mapping_path) if mapping_path else None,
    )


@app.command("nsforest-genes")
def nsforest_genes_command(
    nsforest_path: Path = typer.Option(..., help="Path to NSForest results CSV"),
    output_path: Path = typer.Option("binary_genes.txt", help="Output file for gene list"),
):
    """Extract binary genes from NSForest results"""
    
    nsforest.extract_binary_genes(
        nsforest_path=str(nsforest_path),
        output_path=str(output_path),
    )


@app.command("viz-dotplot")
def viz_dotplot_command(
    h5ad_path: str = typer.Option(..., help="Path to input .h5ad file"),
    embedding_key: str = typer.Option(..., help="Embedding key (X_umap, X_tsne)"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
):
    """Generate UMAP/embedding dotplot colored by cluster"""
    
    viz.plot_dotplot(
        h5ad_path=h5ad_path,
        embedding_key=embedding_key,
        cluster_header=cluster_header,
    )


@app.command("viz-distribution")
def viz_distribution_command(
    cluster_summary_path: Path = typer.Option(..., help="CSV with cluster summary"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
):
    """Generate distribution plots of cluster sizes vs silhouette"""
    
    viz.plot_distribution(
        cluster_summary_path=str(cluster_summary_path),
        cluster_header=cluster_header,
    )


def main():
    app()
