
from pathlib import Path
import typer
from scsilhouette import compute, viz, download

app = typer.Typer()

@app.command("compute-silhouette")
def compute_silhouette_command(
    h5ad_path: Path = typer.Option(...),
    label_keys: list[str] = typer.Option(...),
    embedding_key: str = typer.Option("X_umap"),
    output_dir: Path = typer.Option(Path("results")),
    nsforest_path: Path = typer.Option(None),
    metric: str = typer.Option("euclidean"),
    use_binary_genes: bool = typer.Option(False),
    save_scores: bool = typer.Option(False),
    save_cluster_summary: bool = typer.Option(False),
    save_csv: bool = typer.Option(False),
    show_obs: bool = typer.Option(False),
    qc_correlations: bool = typer.Option(False)
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
        show_obs=show_obs,
        qc_correlations=qc_correlations
    )

@app.command("viz")
def viz_command(
    cluster_summary_file: Path = typer.Option(...),
    cell_scores_file: Path = typer.Option(...),
    output_dir: Path = typer.Option("results"),
    label_key: str = typer.Option(...),
    score_col: str = typer.Option(...),
    suffix: str = typer.Option(""),
):
    viz.plot_all(
        cluster_summary_path=str(cluster_summary_file),
        cell_scores_path=str(cell_scores_file),
        output_dir=str(output_dir),
        label=label_key,
        score_col=score_col,
        suffix=suffix,
    )

@app.command("viz-fscore")
def viz_fscore_command(
    fscore_path: Path = typer.Option(...),
    cluster_summary_path: Path = typer.Option(...),
    output_dir: Path = typer.Option("results"),
    label: str = typer.Option(...),
    score_col: str = typer.Option(...),
    suffix: str = typer.Option(""),
):
    viz.plot_fscore_vs_silhouette(
        fscore_path=str(fscore_path),
        cluster_summary_path=str(cluster_summary_path),
        output_dir=str(output_dir),
        label=label,
        score_col=score_col,
        suffix=suffix,
    )

@app.command("download")
def download_command():
    download.fetch_data()

def main():
    app()
