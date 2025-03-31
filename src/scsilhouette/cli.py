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

@click.group()
def viz():
    """Visualization commands."""
    pass


@viz.command('fscore-corr')
@click.option('--fscore-csv', required=True, type=click.Path(exists=True), help='Path to NSForest F-score CSV')
@click.option('--silhouette-csv', required=True, type=click.Path(exists=True), help='Path to silhouette scores CSV')
@click.option('--output-dir', required=True, type=click.Path(), help='Directory to save plot and correlation CSV')
@click.option('--label-key', required=True, type=str, help='Label column name used for clustering')
@click.option('--score-type', default='silhouette_score_euclidean',
              type=click.Choice(['silhouette_score_euclidean', 'silhouette_score_cosine']),
              help='Silhouette score type to use')
def fscore_corr(fscore_csv, silhouette_csv, output_dir, label_key, score_type):
    """Compute and plot correlation between NS-Forest F-score and silhouette score."""
    from scsilhouette.viz import plot_fscore_vs_silhouette

    fscore_df = pd.read_csv(fscore_csv)
    silhouette_df = pd.read_csv(silhouette_csv)

    plot_fscore_vs_silhouette(fscore_df, silhouette_df, output_dir, label_key, score_col=score_type)

    click.echo(f"[✓] Correlation plot and CSV saved to: {output_dir}")

if __name__ == "__main__":
    main()

