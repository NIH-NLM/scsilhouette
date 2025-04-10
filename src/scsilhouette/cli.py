# src/scsilhouette/cli.py

import click
from pathlib import Path
from scsilhouette import compute, download
from scsilhouette.viz import viz as viz_commands


@click.group()
def main():
    """scSilhouette CLI"""
    pass


@main.command("compute")
@click.option("--h5ad-path", required=True, type=click.Path(exists=True), help="Path to H5AD file")
@click.option("--label-keys", required=True, type=str, help="Label keys in .obs")
@click.option("--embedding-key", required=True, type=str, help="Key in .obsm for embedding")
@click.option("--output-dir", required=True, type=click.Path(), help="Directory to save results")
@click.option("--nsforest-path", type=click.Path(exists=True), help="Optional: NSForest CSV file")
@click.option("--metric", default="euclidean", type=click.Choice(["euclidean", "cosine"]), show_default=True)
@click.option("--use-binary-genes", is_flag=True, help="Restrict analysis to binary genes from NSForest")
@click.option("--save-scores", is_flag=True, help="Save silhouette scores CSV")
@click.option("--save-cluster-summary", is_flag=True, help="Save cluster summary CSV")
@click.option("--save-csv", is_flag=True, help="Alias: Save all CSVs")
@click.option("--save-plots", is_flag=True, help="Save plots as SVG")
@click.option("--qc-correlations", is_flag=True, help="Run correlation analysis against NSForest F-score")
@click.option("--show-obs", is_flag=True, help="Preview first rows of .obs")
def run_compute_command(**kwargs):
    """Run silhouette scoring and QC analysis on single-cell data."""
    compute.run_silhouette(**kwargs)


@main.command("download")
@click.argument("dataset_id")
@click.argument("output_dir", type=click.Path())
def run_download_command(dataset_id, output_dir):
    """Download example H5AD file by dataset_version_id."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    download.fetch_h5ad_file(dataset_id, output_dir)


# Add viz group (e.g. scsilhouette viz fscore-corr ...)
main.add_command(viz_commands, name="viz")


if __name__ == "__main__":
    main()

