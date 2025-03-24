# src/scsilhouette/cli.py

import click
from scsilhouette import compute, download

@click.group()
def main():
    pass

@main.command()
@click.option("--url", required=True, help="URL of the .h5ad file to download.")
@click.option("--output-dir", required=True, help="Directory to save the downloaded file.")
def download_file(url, output_dir):
    download.download_h5ad(url, output_dir)

@main.command(name="compute")
@click.option("--h5ad-path", required=True, type=click.Path(exists=True), help="Path to input .h5ad file.")
@click.option("--label-keys", required=True, multiple=True, help="One or more label keys from adata.obs.")
@click.option("--embedding-key", default="X_pca", help="Key in adata.obsm for embedding. Default: X_pca.")
@click.option("--output-dir", default="results", help="Directory to save results.")
@click.option("--show-obs", is_flag=True, help="Print preview of obs fields.")
@click.option("--log-pca-dims", is_flag=True, help="Log number of PCA components.")
@click.option("--qc-correlations", is_flag=True, help="Run QC correlation plots.")
@click.option("--save-scores", is_flag=True, help="Save silhouette scores as CSV.")
@click.option("--save-cluster-summary", is_flag=True, help="Save cluster summary CSV.")
@click.option("--save-csv", is_flag=True, help="Save all CSVs.")
@click.option("--save-plots", is_flag=True, help="Save all plots.")
def run_compute_command(
    h5ad_path,
    label_keys,
    embedding_key,
    output_dir,
    show_obs,
    log_pca_dims,
    qc_correlations,
    save_scores,
    save_cluster_summary,
    save_csv,
    save_plots,
):
    compute.run_silhouette(
        h5ad_path=h5ad_path,
        label_keys=label_keys,
        embedding_key=embedding_key,
        output_dir=output_dir,
        show_obs=show_obs,
        log_pca_dims=log_pca_dims,
        qc_correlations=qc_correlations,
        save_scores=save_scores,
        save_cluster_summary=save_cluster_summary,
        save_csv=save_csv,
        save_plots=save_plots,
    )

if __name__ == "__main__":
    main()

