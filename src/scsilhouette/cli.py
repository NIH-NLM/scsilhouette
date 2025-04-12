# src/scsilhouette/cli.py

import typer
from pathlib import Path
from scsilhouette import compute, download, viz

app = typer.Typer(help="scSilhouette CLI for computing and visualizing silhouette scores")

@app.command("compute-silhouette")
def compute_silhouette(
    h5ad_path: Path,
    label_keys: str,
    embedding_key: str = "X_umap",
    output_dir: Path = Path("results"),
    nsforest_path: Path = None,
    metric: str = "euclidean",
    use_binary_genes: bool = False,
    save_scores: bool = False,
    save_cluster_summary: bool = False,
    save_csv: bool = False,
    save_plots: bool = False,
    show_obs: bool = False,
    qc_correlations: bool = False
):
    compute.run_silhouette(
        h5ad_path=h5ad_path,
        label_key=label_keys,
        embedding_key=embedding_key,
        output_dir=output_dir,
        nsforest_path=nsforest_path,
        metric=metric,
        use_binary_genes=use_binary_genes,
        save_scores=save_scores,
        save_cluster_summary=save_cluster_summary,
        save_csv=save_csv,
        save_plots=save_plots,
        show_obs=show_obs,
        qc_correlations=qc_correlations
    )

@app.command("viz")
def run_viz(
    cluster_summary_path: Path,
    cell_scores_path: Path,
    output_dir: Path,
    label: str,
    score_col: str,
    suffix: str = ""
):
    viz.plot_all(
        cluster_summary_path=cluster_summary_path,
        cell_scores_path=cell_scores_path,
        output_dir=output_dir,
        label=label,
        score_col=score_col,
        suffix=suffix
    )

@app.command("download")
def download_test_file(
    dataset_id: str,
    output_dir: Path = Path("data")
):
    download.fetch_dataset(dataset_id, output_dir)

def main():
    app()

if __name__ == "__main__":
    main()

