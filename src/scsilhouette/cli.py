# src/scsilhouette/cli.py

import click
from scsilhouette.compute import run_silhouette
from scsilhouette.viz import plot_all
from scsilhouette.download import download_h5ad


@click.group()
def cli():
    """scSilhouette CLI for computing silhouette scores and generating plots."""
    pass


@cli.command()
@click.option("--h5ad-path", required=True, type=click.Path(exists=True), help="Path to .h5ad file")
@click.option("--label-keys", required=True, multiple=True, help="One or more .obs label fields (e.g. cell_type)")
@click.option("--embedding-key", default="X_pca", help="Key for embedding (e.g. X_pca, X_umap)")
@click.option("--output-dir", default="results", type=click.Path(), help="Directory for saving output files")
def compute(h5ad_path, label_keys, embedding_key, output_dir):
    run_silhouette(h5ad_path, label_keys, embedding_key, output_dir)


@cli.command()
@click.option("--input-dir", default="results", type=click.Path(exists=True), help="Input directory with scores")
def viz(input_dir):
    plot_all(input_dir)


@cli.command()
@click.option("--url", required=True, help="URL to download the H5AD file")
@click.option("--output-dir", default="data", type=click.Path(), help="Download location")
def download(url, output_dir):
    download_h5ad(url, output_dir)


if __name__ == "__main__":
    cli()

