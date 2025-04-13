# src/scsilhouette/cli.py

import typer
from pathlib import Path
from typing import List, Optional

from scsilhouette import compute, download, viz

app = typer.Typer(add_completion=False)


@app.command("compute-silhouette")
def compute_silhouette_command(
    h5ad_path: Path = typer.Option(..., help="Path to the .h5ad file"),
    label_keys: List[str] = typer.Option(..., help="List of label keys to use"),
    embedding_key: str = typer.Option("X_umap", help="Embedding key to use"),
    output_dir: Path = typer.Option("results", help="Directory to store results"),
    nsforest_path: Optional[Path] = typer.Option(None, help="Path to NSForest result CSV"),
    metric: str = typer.Option("euclidean", help="Distance metric to use"),
    use_binary_genes: bool = typer.Option(False, help="Use binary genes from NSForest"),
    save_scores: bool = typer.Option(False, help="Save cell-level silhouette scores"),
    save_cluster_summary: bool = typer.Option(False, help="Save cluster-level silhouette summaries"),
    save_csv: bool = typer.Option(False, help="Save CSV files"),
    save_plots: bool = typer.Option(False, help="Save cluster plots"),
    show_obs: bool = typer.Option(False, help="Print and save preview of adata.obs"),
    qc_correlations: bool = typer.Option(False, help="Perform QC correlation analysis"),
):
    compute.run_silhouette(
        h5ad_path=str(h5ad_path),
        label_keys=label_keys,
        embedding_key=embedding_key,
        output_dir=str(output_dir),
        nsforest_path=str(nsforest_path) if nsforest_path else None,
        metric=metric,
        use_binary_genes=use_binary_genes,
        save_scores=save_scores,
        save_cluster_summary=save_cluster_summary,
        save_csv=save_csv,
        save_plots=save_plots,
        show_obs=show_obs,
        qc_correlations=qc_correlations,
    )


@app.command("download")
def download_data():
    download.run()


@app.command("viz")
def visualize():
    typer.echo("✨ Visualization commands go here.")
    # You can expand this with subcommands later for plot_all, etc.


def main():
    app()


if __name__ == "__main__":
    main()

