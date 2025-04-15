# src/scsilhouette/cli.py

import typer
from pathlib import Path
from typing import Optional
import pandas as pd
from . import viz, compute, download, nsforest

app = typer.Typer()

@app.command("download")
def download_command(
    url: str = typer.Option(..., help="URL to download .h5ad dataset"),
    output_dir: Path = typer.Option(..., help="Where to save the file")
):
    path = download.download_h5ad(url, str(output_dir))
    print(f"[\u2713] Downloaded: {path}")

@app.command("compute-silhouette")
def compute_command(
    h5ad_path: Path = typer.Option(...),
    label_keys: str = typer.Option(...),
    embedding_key: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
    nsforest_path: Optional[Path] = typer.Option(None),
    use_binary_genes: bool = typer.Option(False),
    metric: str = typer.Option("euclidean"),
    save_scores: bool = typer.Option(True),
    save_cluster_summary: bool = typer.Option(True),
    save_csv: bool = typer.Option(False),
    qc_correlations: bool = typer.Option(False),
    show_obs: bool = typer.Option(False),
):
    compute.run_silhouette(
        h5ad_path=str(h5ad_path),
        label_keys=label_keys.split(","),
        embedding_key=embedding_key,
        output_dir=str(output_dir),
        nsforest_path=str(nsforest_path) if nsforest_path else None,
        use_binary_genes=use_binary_genes,
        metric=metric,
        save_scores=save_scores,
        save_cluster_summary=save_cluster_summary,
        save_csv=save_csv,
        qc_correlations=qc_correlations,
        show_obs=show_obs,
    )

@app.command("viz-summary")
def viz_summary_command(
    silhouette_score_path: Path = typer.Option(...),
    output_dir: Path = typer.Option(...),
    label: str = typer.Option(...),
    score_col: str = typer.Option(...),
    suffix: str = typer.Option(""),
    show: bool = typer.Option(False),
    fscore_path: Optional[Path] = typer.Option(None),
    mapping_path: Optional[Path] = typer.Option(None)
):
    viz.plot_silhouette_summary(
        silhouette_score_path=str(silhouette_score_path),
        output_dir=str(output_dir),
        label=label,
        score_col=score_col,
        suffix=suffix,
        show=show,
        fscore_path=str(fscore_path) if fscore_path else None,
        mapping_path=str(mapping_path) if mapping_path else None,
    )

@app.command("viz-fscore")
def viz_fscore_command(
    fscore_path: Path = typer.Option(...),
    cluster_summary_path: Path = typer.Option(...),
    output_dir: Path = typer.Option(...),
    label: str = typer.Option(...),
    score_col: str = typer.Option(...),
    suffix: str = typer.Option(""),
    show: bool = typer.Option(False),
    export_csv: bool = typer.Option(False),
    mapping_path: Optional[Path] = typer.Option(None),
    summary_path: Optional[Path] = typer.Option(None),
    silhouette_stat: str = typer.Option("mean")
):
    silhouette_col = f"{silhouette_stat}_silhouette_euclidean"
    viz.plot_fscore_vs_silhouette_files(
        fscore_path=str(fscore_path),
        cluster_summary_path=str(cluster_summary_path),
        output_dir=str(output_dir),
        label=label,
        score_col=silhouette_col,
        suffix=suffix,
        show=show,
        export_csv=export_csv,
        mapping_path=str(mapping_path) if mapping_path else None,
        summary_path=str(summary_path) if summary_path else None,
        silhouette_stat=silhouette_stat
    )

@app.command("nsforest-validate")
def nsforest_validate_command(
    silhouette_csv: Path = typer.Option(...),
    nsforest_csv: Path = typer.Option(...),
    label_key: str = typer.Option(...)
):
    silhouette_df = pd.read_csv(silhouette_csv)
    _, fscore_df = nsforest.load_nsforest_binary_genes(nsforest_csv)
    unmatched = nsforest.validate_cluster_alignment(silhouette_df, fscore_df, label_key)
    if unmatched:
        print(f"[!] Unmatched labels: {unmatched}")
    else:
        print("[\u2713] All clusters match!")

main = app

if __name__ == "__main__":
    app()

