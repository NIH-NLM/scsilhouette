# src/scsilhouette/viz.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr
import click

def plot_score_distribution(cell_scores: pd.DataFrame, output_dir: str, label: str, score_col: str, suffix: str):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 6))
    sns.violinplot(x=label, y=score_col, data=cell_scores, inner='box')
    plt.title(f"Silhouette Score Distribution ({suffix})")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_score_distribution_{suffix}.svg")
    plt.close()

def plot_cluster_summary(cluster_summary: pd.DataFrame, output_dir: str, label: str, score_col: str, suffix: str):
    plt.figure(figsize=(10, 6))
    cluster_summary.plot.bar(x=label, y=score_col, legend=False)
    plt.title(f"Average Silhouette Score per Cluster ({suffix})")
    plt.ylabel("Silhouette Score")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_cluster_summary_{suffix}.svg")
    plt.close()

def plot_cluster_size_vs_score(cluster_summary: pd.DataFrame, output_dir: str, label: str, score_col: str, suffix: str):
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='n_cells', y=score_col, data=cluster_summary)
    plt.title(f"Cluster Size vs Silhouette Score ({suffix})")
    plt.xlabel("Cluster Size")
    plt.ylabel("Silhouette Score")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_cluster_size_vs_score_{suffix}.svg")
    plt.close()

def plot_fscore_vs_silhouette(fscore_df, silhouette_df, output_dir, label_key, score_col):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    merged = pd.merge(fscore_df, silhouette_df, how='inner', left_on='cluster', right_on=label_key)
    if score_col not in merged.columns:
        raise ValueError(f"{score_col} missing")
    r, p = pearsonr(merged[score_col], merged["F-score"])
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=merged, x=score_col, y="F-score")
    plt.title(f"{score_col} vs F-score\nr={r:.2f}, p={p:.3e}")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"fscore_vs_{score_col}.svg")
    plt.close()
    merged[[label_key, score_col, "F-score"]].to_csv(Path(output_dir) / f"fscore_vs_{score_col}_correlation.csv", index=False)

def plot_fscore_vs_silhouette_all(fscore_df, silhouette_df, output_dir, label_key):
    for score_col in ["silhouette_score_euclidean", "silhouette_score_cosine"]:
        if score_col in silhouette_df.columns:
            plot_fscore_vs_silhouette(fscore_df, silhouette_df, output_dir, label_key, score_col)

def plot_all(cluster_summary, output_dir, label, score_col, suffix, cell_scores):
    plot_score_distribution(cell_scores, output_dir, label, score_col, suffix)
    plot_cluster_summary(cluster_summary, output_dir, label, score_col, suffix)
    plot_cluster_size_vs_score(cluster_summary, output_dir, label, score_col, suffix)

@click.command("fscore-corr")
@click.option('--fscore-csv', required=True, type=click.Path(exists=True))
@click.option('--silhouette-csv', required=True, type=click.Path(exists=True))
@click.option('--output-dir', required=True, type=click.Path())
@click.option('--label-key', required=True, type=str)
@click.option('--score-type', default='silhouette_score_euclidean', type=click.Choice(['silhouette_score_euclidean', 'silhouette_score_cosine']))
def fscore_corr(fscore_csv, silhouette_csv, output_dir, label_key, score_type):
    fscore_df = pd.read_csv(fscore_csv)
    silhouette_df = pd.read_csv(silhouette_csv)
    plot_fscore_vs_silhouette(fscore_df, silhouette_df, output_dir, label_key, score_type)
    click.echo(f"[✓] Saved correlation plot and CSV for {score_type} to {output_dir}")

@click.command("fscore-corr-all")
@click.option('--fscore-csv', required=True, type=click.Path(exists=True))
@click.option('--silhouette-csv', required=True, type=click.Path(exists=True))
@click.option('--output-dir', required=True, type=click.Path())
@click.option('--label-key', required=True, type=str)
def fscore_corr_all(fscore_csv, silhouette_csv, output_dir, label_key):
    fscore_df = pd.read_csv(fscore_csv)
    silhouette_df = pd.read_csv(silhouette_csv)
    plot_fscore_vs_silhouette_all(fscore_df, silhouette_df, output_dir, label_key)
    click.echo(f"[✓] Saved all correlation plots and CSVs to {output_dir}")

viz = click.Group()
viz.add_command(fscore_corr)
viz.add_command(fscore_corr_all)

