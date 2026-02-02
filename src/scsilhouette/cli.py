import typer
from pathlib import Path
from typing import Optional, List

from . import compute, download, viz, nsforest

app = typer.Typer(add_completion=False)

@app.command("compute-silhouette")
def compute_silhouette_command(
    h5ad_path: Path = typer.Option(...),
    label_key: str = typer.Option(...),
    embedding_key: str = typer.Option(...),
    organism: str = typer.Option(...),
    disease: str = typer.Option(...),
    tissue: str = typer.Option(...),
    cell_count: str = typer.Option(...),
    use_binary_genes: bool = typer.Option(False),
    gene_list_path: Optional[Path] = typer.Option(None),
    metric: str = typer.Option("euclidean"),
    filter_normal: str = typer.Option(True),
    save_scores: bool = typer.Option(False),
    save_cluster_summary: bool = typer.Option(False),
    save_annotation: bool = typer.Option(False),
):
    compute.run_silhouette(
        h5ad_path=str(h5ad_path),
        label_key=label_key,
        embedding_key=embedding_key,
        organism=organism,
        disease=disease,
        tissue=tissue,
        cell_count=cell_count,
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
    silhouette_score_path: Path = typer.Option(...),
    silhouette_score_col: str = typer.Option(...),
    label_key: str = typer.Option(...),
    fscore_path: Optional[Path] = typer.Option(None),
    sort_by: str = typer.Option("median", help="Sort by mean|median|std"),
):
    viz.plot_silhouette_summary(
        silhouette_score_path=str(silhouette_score_path),
        silhouette_score_col=silhouette_score_col,
        label_key=label_key,
        fscore_path=str(fscore_path) if fscore_path else None,
        sort_by=sort_by,
    )
    
@app.command("viz-correlation")
def viz_correlation_command(
    cluster_summary_path: Path = typer.Option(..., help="CSV with summary metrics (e.g., mean, fscore, count)"),
    x_metric: str = typer.Option(..., help="X-axis metric for correlation (e.g., fscore)"),
    y_metrics: str = typer.Option(..., help="Comma-separated list of Y-axis metrics (e.g., mean,median,std)"),
    label: str = typer.Option(..., help="Label column used for filenames and groupings"),
    fscore_path: Optional[Path] = typer.Option(None, help="Optional path to fscore CSV"),
    mapping_path: Optional[Path] = typer.Option(None, help="Optional mapping file to match cluster labels")
):
    viz.plot_correlation_summary(
        cluster_summary_path=str(cluster_summary_path),
        x_metric=x_metric,
        y_metrics=[x.strip() for x in y_metrics.split(",")],
        label=label,
        fscore_path=str(fscore_path) if fscore_path else None,
        mapping_path=str(mapping_path) if mapping_path else None,
    )

@app.command("nsforest-genes")
def nsforest_genes_command(
    nsforest_path: Path = typer.Option(...),
):
    nsforest.extract_binary_genes(
        nsforest_path=str(nsforest_path),
    )

@app.command("viz-dotplot")
def viz_dotplot_command(
    h5ad_path: str = typer.Option(..., help="Path to input .h5ad file"),
    embedding_key: str = typer.Option(..., help="Embedding key (X_umap, X_scanvi_emb)"),
    label_key: str = typer.Option(..., help="Label key for color"),
):
    viz.plot_dotplot(
        h5ad_path=h5ad_path,
        embedding_key=embedding_key,
        label_key=label_key,
    )

@app.command("viz-distribution")
def viz_distribution_command(
    cluster_summary_path: Path = typer.Option(..., help="CSV file from compute-silhouette with mean, median, std, count"),
    label_key: str = typer.Option(..., help="Column name used as x-axis (e.g., ann_finest_level)"),
):
    viz.plot_distribution(
        cluster_summary_path=str(cluster_summary_path),
        label_key=label_key,
    )


def main():
    app()

