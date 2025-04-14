# src/scsilhouette/cli.py

import typer
from pathlib import Path
from typing import Optional
from . import viz

app = typer.Typer()

@app.command("viz-summary")
def viz_summary_command(
    silhouette_score_path: Path = typer.Option(..., help="Path to silhouette score CSV"),
    output_dir: Path = typer.Option(..., help="Path to save plots and summary CSV"),
    label: str = typer.Option(..., help="Column to group by (e.g. cell_type)"),
    score_col: str = typer.Option(..., help="Silhouette score column name"),
    fscore_path: Optional[Path] = typer.Option(None, help="Optional path to F-score CSV"),
    mapping_path: Optional[Path] = typer.Option(None, help="Optional path to clusterName-to-cell_type mapping CSV"),
    show: bool = typer.Option(False, help="Show plot interactively")
):
    viz.plot_silhouette_summary(
        silhouette_score_path=str(silhouette_score_path),
        output_dir=str(output_dir),
        label_col=label,
        score_col=score_col,
        suffix="",
        show=show,
        fscore_path=str(fscore_path) if fscore_path else None,
        mapping_path=str(mapping_path) if mapping_path else None,
    )


@app.command("viz-fscore")
def viz_fscore_command(
    fscore_path: Path = typer.Option(..., help="Path to f-score CSV"),
    cluster_summary_path: Path = typer.Option(..., help="Path to silhouette summary CSV"),
    output_dir: Path = typer.Option(..., help="Output directory"),
    label: str = typer.Option(..., help="Label column name to join on (e.g. cell_type)"),
    score_col: str = typer.Option(..., help="Silhouette score column to correlate (e.g. mean or median)"),
    suffix: str = typer.Option("", help="Suffix for output file names"),
    show: bool = typer.Option(False, help="Show plot interactively"),
    export_csv: bool = typer.Option(False, help="Export merged and stats CSVs"),
    mapping_path: Optional[Path] = typer.Option(None, help="Optional mapping CSV (clusterName -> cell_type)"),
    summary_path: Optional[Path] = typer.Option(None, help="Optional path to precomputed silhouette summary CSV")
):
    viz.plot_fscore_vs_silhouette_files(
        fscore_path=str(fscore_path),
        cluster_summary_path=str(cluster_summary_path),
        output_dir=str(output_dir),
        label=label,
        score_col=score_col,
        suffix=suffix,
        show=show,
        export_csv=export_csv,
        mapping_path=str(mapping_path) if mapping_path else None,
        summary_path=str(summary_path) if summary_path else None,
    )


# 🔧 RECOMMENDED: Explicit CLI entry point for setuptools or pyproject
main = app

if __name__ == "__main__":
    app()
