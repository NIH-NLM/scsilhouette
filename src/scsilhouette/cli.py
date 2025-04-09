# src/scsilhouette/cli.py

import os
import click
from scsilhouette import compute
from scsilhouette.viz import viz as viz_commands

@click.group()
def main():
    pass

@main.command("compute")
@click.option("--h5ad-path", required=True, type=click.Path(exists=True), help="Path to .h5ad file")
@click.option("--label-keys", required=True, type=str, help="Label column(s) in .obs for cluster identity")
@click.option("--embedding-key", required=True, type=str, help="Key in .obsm for embeddings (e.g. X_pca or X_umap)")
@click.option("--output-dir", required=True, type=click.Path(), help="Directory to write results")
@click.option("--show-obs", is_flag=True, help="Save preview of AnnData .obs as CSV")
@click.option("--metric", default="euclidean", type=click.Choice(["euclidean", "cosine"]), help="Distance metric for silhouette score")
@click.option("--save-scores", is_flag=True, help="Save per-cell silhouette scores to CSV")
@click.option("--save-cluster-summary", is_flag=True, help="Save per-cluster silhouette stats to CSV")
@click.option("--save-csv", is_flag=True, help="Save all generated data tables")
@click.option("--save-plots", is_flag=True, help="Save visualizations to PNG and SVG")
@click.option("--qc-correlations", is_flag=True, help="Perform quality control correlations and plots")
@click.option("--nsforest-path", type=click.Path(exists=True), help="Path to NS-Forest output CSV")
@click.option("--use-binary-genes", is_flag=True, help="Subset features using NS-Forest binary genes")
def run_compute_command(**kwargs):
    compute.run_silhouette(**kwargs)

# Include viz subcommands
main.add_command(viz_commands)

if __name__ == "__main__":
    main()

