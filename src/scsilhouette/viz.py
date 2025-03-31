# src/scsilhouette/viz.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr
from pathlib import Path
from typing import Literal


def plot_score_distribution(cell_scores, output_dir, label):
    plt.figure(figsize=(8, 4))
    sns.histplot(cell_scores["silhouette_score"], bins=50, kde=True)
    plt.title(f"Silhouette Score Distribution - {label}")
    plt.xlabel("Silhouette Score")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_score_distribution.png")
    plt.close()


def plot_cluster_summary(cluster_summary, output_dir, label):
    cluster_summary.plot.bar(x=label, y="mean_silhouette_score", legend=False)
    plt.title(f"Mean Silhouette Score by {label}")
    plt.ylabel("Mean Silhouette Score")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_cluster_summary.png")
    plt.close()


def plot_cluster_size_vs_score(cluster_summary, output_dir, label):
    plt.figure(figsize=(6, 6))
    plt.scatter(cluster_summary["n_cells"], cluster_summary["mean_silhouette_score"], alpha=0.7)
    plt.xlabel("Cluster Size (n_cells)")
    plt.ylabel("Mean Silhouette Score")
    plt.title(f"{label}: Cluster Size vs Silhouette Score")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_cluster_size_vs_score.png")
    plt.close()


def plot_qc_boxplots(cell_scores, output_dir, label):
    for metric in ["nCount_RNA", "nFeature_RNA"]:
        if metric not in cell_scores:
            continue
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=cell_scores, x=label, y=metric)
        sns.stripplot(data=cell_scores, x=label, y=metric, color="black", alpha=0.3, jitter=0.2)
        plt.xticks(rotation=90)
        plt.title(f"{metric} by {label}")
        plt.tight_layout()
        plt.savefig(Path(output_dir) / f"{label}_{metric}_boxplot.png")
        plt.close()


def plot_all(cluster_summary, cell_scores, output_dir, label):
    plot_score_distribution(cell_scores, output_dir, label)
    plot_cluster_summary(cluster_summary, output_dir, label)
    plot_cluster_size_vs_score(cluster_summary, output_dir, label)
    plot_qc_boxplots(cell_scores, output_dir, label)


def plot_fscore_vs_silhouette(merged_df, output_dir, score_col="mean_silhouette_score", fscore_col="F-score"):
    corr = merged_df[score_col].corr(merged_df[fscore_col])
    plt.figure(figsize=(7, 6))
    sns.scatterplot(data=merged_df, x=score_col, y=fscore_col, hue="cluster", palette="tab10")
    plt.xlabel("Silhouette Score")
    plt.ylabel("NSForest F-score")
    plt.title(f"{fscore_col} vs. {score_col} (r={corr:.2f})")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"fscore_vs_{score_col}.png")
    plt.close()

def plot_fscore_vs_silhouette(
    fscore_df: pd.DataFrame,
    silhouette_df: pd.DataFrame,
    output_dir: str,
    label_key: str,
    score_col: Literal['silhouette_score_euclidean', 'silhouette_score_cosine'] = 'silhouette_score_euclidean'
) -> None:

    """Plot F-score vs Silhouette Score with correlation and save the plot."""

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    merged = pd.merge(fscore_df, silhouette_df, how='inner', left_on='cluster', right_on=label_key)

    if score_col not in merged.columns:
        raise ValueError(f"'{score_col}' not found in merged dataframe")

    r, p = pearsonr(merged[score_col], merged['F-score'])

    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=score_col, y='F-score', data=merged)
    plt.title(f'F-score vs {score_col}\nPearson r={r:.2f}, p={p:.2e}')
    plt.xlabel(score_col)
    plt.ylabel('F-score')
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f'fscore_vs_{score_col}.png')
    plt.close()

    merged[[label_key, score_col, 'F-score']].to_csv(
        Path(output_dir) / f'fscore_vs_{score_col}_correlation.csv', index=False
    )
