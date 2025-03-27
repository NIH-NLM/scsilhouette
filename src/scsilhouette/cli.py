# src/scsilhouette/cli.py

import click
from scsilhouette import compute, download, viz
import pandas as pd

@click.group()
def main():
    """scsilhouette CLI: compute silhouette scores, visualize, and manage data."""
    pass

@main.command("download")
@click.option("--url", required=True, help="URL to download .h5ad file from")
@click.option("--output-dir", default="data", show_default=True, help="Directory to save downloaded file")
def run_download_command(url, output_dir):
    """Download .h5ad dataset from URL."""
    download.download_h5ad(url, output_dir)

@main.command("compute")
@click.option("--h5ad-path", required=True, type=click.Path(exists=True), help="Path to input .h5ad file")
@click.option("--label-keys", required=True, multiple=True, help="Label keys in .obs to compute silhouette")
@click.option("--embedding-key", default="X_pca", help="Key in .obsm for embedding (e.g., X_pca, X_umap)")
@click.option("--output-dir", default="results", show_default=True, help="Directory to save outputs")
@click.option("--show-obs", is_flag=True, help="Print and save .obs preview")
@click.option("--qc-correlations", is_flag=True, help="Compute QC correlations with silhouette scores")
@click.option("--save-scores", is_flag=True, help="Save full silhouette scores per cell")
@click.option("--save-cluster-summary", is_flag=True, help="Save per-cluster summary CSV")
@click.option("--save-csv", is_flag=True, help="Save all cluster summaries and scores")
@click.option("--save-plots", is_flag=True, help="Save all relevant plots")
@click.option("--nsforest-path", type=click.Path(exists=True), help="Path to NSForest results (CSV or Excel)")
def run_compute_command(**kwargs):
    """Compute silhouette scores and analyze results."""
    compute.run_silhouette(**kwargs)

@main.command("viz")
@click.option("--summary-csv", required=True, type=click.Path(exists=True), help="Path to cluster summary CSV")
@click.option("--scores-csv", required=True, type=click.Path(exists=True), help="Path to scores CSV")
@click.option("--label", required=True, help="Label name used in plots")
@click.option("--output-dir", default="results", show_default=True, help="Directory to save figures")
def run_viz_command(summary_csv, scores_csv, label, output_dir):
    """Generate all silhouette visualization plots."""
    cluster_summary = pd.read_csv(summary_csv)
    cell_scores = pd.read_csv(scores_csv)
    viz.plot_all(cluster_summary, cell_scores, output_dir, label)

if __name__ == "__main__":
    main()

