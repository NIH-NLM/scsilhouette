# src/scsilhouette/cli.py
import click
from . import compute, viz, download


@click.group()
def main():
    pass


@main.command("compute")
@click.option("--h5ad-path", required=True, type=click.Path(exists=True))
@click.option("--label-keys", required=True, multiple=True)
@click.option("--embedding-key", required=True)
@click.option("--output-dir", required=True, type=click.Path())
@click.option("--show-obs", is_flag=True, default=False)
@click.option("--save-scores", is_flag=True, default=False)
@click.option("--save-cluster-summary", is_flag=True, default=False)
@click.option("--save-csv", is_flag=True, default=False)
@click.option("--save-plots", is_flag=True, default=False)
@click.option("--qc-correlations", is_flag=True, default=False)
@click.option("--nsforest-path", type=click.Path(exists=True), default=None)
def run_compute_command(**kwargs):
    compute.run_silhouette(**kwargs)


@main.command("viz")
@click.option("--fscore-csv", required=True, type=click.Path(exists=True))
@click.option("--silhouette-csv", required=True, type=click.Path(exists=True))
@click.option("--output-dir", required=True, type=click.Path())
@click.option("--label-key", required=True, type=str)
@click.option(
    "--score-type",
    default="silhouette_score_euclidean",
    type=click.Choice(["silhouette_score_euclidean", "silhouette_score_cosine"]),
)
def fscore_corr(fscore_csv, silhouette_csv, output_dir, label_key, score_type):
    fscore_df = pd.read_csv(fscore_csv)
    silhouette_df = pd.read_csv(silhouette_csv)
    viz.plot_fscore_vs_silhouette(fscore_df, silhouette_df, output_dir, label_key, score_col=score_type)


@main.command("download")
def run_download():
    download.download_all()


if __name__ == "__main__":
    main()

