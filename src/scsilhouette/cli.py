# src/scsilhouette/cli.py

import typer
from pathlib import Path
from typing import Optional

from . import compute, download, viz, nsforest

app = typer.Typer(add_completion=False)

@app.command("compute-silhouette")
def compute_silhouette_command(
    h5ad_path: Path = typer.Option(...),
    label_key: str = typer.Option(...),
    embedding_key: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
    use_binary_genes: bool = typer.Option(False),
    gene_list_path: Optional[Path] = typer.Option(None),
    metric: str = typer.Option("euclidean"),
    save_scores: bool = typer.Option(False),
    save_cluster_summary: bool = typer.Option(False),
    save_csv: bool = typer.Option(False),
    show_obs: bool = typer.Option(False)
):
    compute.run_silhouette(
        h5ad_path=str(h5ad_path),
        label_key=label_key,
        embedding_key=embedding_key,
        output_dir=str(output_dir),
        use_binary_genes=use_binary_genes,
        gene_list_path=str(gene_list_path) if gene_list_path else None,
        metric=metric,
        save_scores=save_scores,
        save_cluster_summary=save_cluster_summary,
        save_csv=save_csv,
        show_obs=show_obs
    )

@app.command("viz-summary")
def viz_summary_command(
    silhouette_score_path: Path = typer.Option(...),
    output_dir: Path = typer.Option(...),
    label: str = typer.Option(...),
    score_col: str = typer.Option(...),
    fscore_path: Optional[Path] = typer.Option(None),
    mapping_path: Optional[Path] = typer.Option(None),
    suffix: str = typer.Option(""),
    show: bool = typer.Option(False),
    sort_by: str = typer.Option("median", help="Sort by mean|median|std")
):
    viz.plot_silhouette_summary(
        silhouette_score_path=str(silhouette_score_path),
        output_dir=str(output_dir),
        label=label,
        score_col=score_col,
        fscore_path=str(fscore_path) if fscore_path else None,
        mapping_path=str(mapping_path) if mapping_path else None,
        suffix=suffix,
        show=show,
        sort_by=sort_by,
    )

@app.command("viz-correlation")
def viz_correlation_command(
    cluster_summary_path: Path = typer.Option(...),
    output_dir: Path = typer.Option(...),
    silhouette_metrics: str = typer.Option(...),
    score_col: str = typer.Option(...),
    suffix: str = typer.Option(""),
    show: bool = typer.Option(False)
):
    viz.plot_correlation_summary(
        cluster_summary_path=str(cluster_summary_path),
        output_dir=str(output_dir),
        silhouette_metrics=[x.strip() for x in silhouette_metrics.split(",")],
        score_col=score_col,
        suffix=suffix,
        show=show,
    )

@app.command("nsforest-genes")
def nsforest_genes_command(
    nsforest_path: Path = typer.Option(...),
    output_path: Path = typer.Option(..., help="Path to write unique binary gene list")
):
    nsforest.extract_binary_genes(
        nsforest_path=str(nsforest_path),
        output_path=str(output_path)
    )

def main():
    app()

