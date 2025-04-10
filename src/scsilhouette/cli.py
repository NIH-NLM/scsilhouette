# src/scsilhouette/cli.py

import click
import anndata as ad
from scsilhouette import compute, download
from scsilhouette.viz import viz as viz_commands

@click.group()
def main():
    """SCSilhouette CLI entry point."""
    pass

@main.command("download")
@click.option("--output-dir", required=True, type=click.Path(), help="Directory to save downloaded files")
def run_download_command(output_dir):
    """Download example files."""
    download.run_download(output_dir)

@main.command("compute")
@click.option("--h5ad-path", required=True, type=click.Path(exists=True), help="Path to h5ad file")
@click.option("--label-keys", required=True, type=str, help="Cluster label key in adata.obs")
@click.option("--embedding-key", required=True, type=str, help="Key for 2D embedding in adata.obsm")
@click.option("--output-dir", required=True, type=click.Path(), help="Directory to save output files")
@click.option("--nsforest-path", required=False, type=click.Path(exists=True), help="Optional path to NSForest CSV file")
@click.option("--metric", default="euclidean", type=click.Choice(["euclidean", "cosine"]), help="Distance metric for silhouette")
@click.option("--use-binary-genes", is_flag=True, help="Restrict silhouette calc to binary marker genes from NSForest")
@click.option("--save-scores", is_flag=True, help="Save per-cell silhouette scores")
@click.option("--save-cluster-summary", is_flag=True, help="Save silhouette score summary per cluster")
@click.option("--save-csv", is_flag=True, help="Save merged score/F-score correlation CSV")
@click.option("--save-plots", is_flag=True, help="Save all correlation and cluster plots")
@click.option("--qc-correlations", is_flag=True, help="Run QC plots comparing silhouette to F-score")
@click.option("--show-obs", is_flag=True, help="Save obs head CSV for debugging")
def run_compute_command(**kwargs):
    """Compute silhouette scores and correlations."""
    click.echo("[✓] Arguments received:")
    for k, v in kwargs.items():
        click.echo(f"    --{k}={v}")

    # ✅ Safety check for valid label_key
    h5ad_path = kwargs.get("h5ad_path")
    label_key = kwargs.get("label_keys")
    if h5ad_path and label_key:
        adata = ad.read_h5ad(h5ad_path)
        if label_key not in adata.obs.columns:
            raise click.BadParameter(f"[ERROR] Label key '{label_key}' not found in adata.obs")

    compute.run_silhouette(**kwargs)

# expose viz commands
main.add_command(viz_commands, name="viz")

