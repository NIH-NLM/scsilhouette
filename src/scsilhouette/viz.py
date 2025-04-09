
# src/scsilhouette/viz.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr


def plot_score_distribution(cell_scores, output_dir, label, score_col, suffix="", file_format="svg"):
    plt.figure(figsize=(8, 4))
    sns.histplot(cell_scores[score_col], bins=50, kde=True)
    plt.title(f"Silhouette Score Distribution - {label} ({suffix})")
    plt.xlabel(score_col)
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_{suffix}_score_distribution.{file_format}", format=file_format)
    plt.close()


def plot_cluster_summary(cluster_summary, output_dir, label, score_col, suffix="", file_format="svg"):
    cluster_summary.plot.bar(x=label, y=score_col, legend=False)
    plt.title(f"Mean Silhouette Score by {label} ({suffix})")
    plt.ylabel("Mean Silhouette Score")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_{suffix}_cluster_summary.{file_format}", format=file_format)
    plt.close()


def plot_cluster_size_vs_score(cluster_summary, output_dir, label, score_col, suffix="", file_format="svg"):
    plt.figure(figsize=(6, 6))
    plt.scatter(cluster_summary["n_cells"], cluster_summary[score_col], alpha=0.7)
    plt.xlabel("Cluster Size (n_cells)")
    plt.ylabel("Mean Silhouette Score")
    plt.title(f"{label} ({suffix}): Cluster Size vs Silhouette Score")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_{suffix}_cluster_size_vs_score.{file_format}", format=file_format)
    plt.close()


def plot_qc_boxplots(cell_scores, output_dir, label, suffix="", file_format="svg"):
    for metric in ["nCount_RNA", "nFeature_RNA"]:
        if metric not in cell_scores:
            continue
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=cell_scores, x=label, y=metric)
        sns.stripplot(data=cell_scores, x=label, y=metric, color="black", alpha=0.3, jitter=0.2)
        plt.xticks(rotation=90)
        plt.title(f"{metric} by {label} ({suffix})")
        plt.tight_layout()
        plt.savefig(Path(output_dir) / f"{label}_{suffix}_{metric}_boxplot.{file_format}", format=file_format)
        plt.close()


def plot_all(cluster_summary, cell_scores, output_dir, label, score_col, suffix="", file_format="svg"):
    plot_score_distribution(cell_scores, output_dir, label, score_col, suffix, file_format)
    plot_cluster_summary(cluster_summary, output_dir, label, score_col, suffix, file_format)
    plot_cluster_size_vs_score(cluster_summary, output_dir, label, score_col, suffix, file_format)
    plot_qc_boxplots(cell_scores, output_dir, label, suffix, file_format)


def plot_fscore_vs_silhouette(fscore_df, silhouette_df, output_dir, label_key, score_col, suffix=""):
    merged = pd.merge(fscore_df, silhouette_df, how='inner', left_on='cluster', right_on=label_key)

    if score_col not in merged.columns:
        raise ValueError(f"'{score_col}' not found in merged dataframe")

    r, p = pearsonr(merged[score_col], merged["F-score"])
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=score_col, y="F-score", data=merged)
    plt.title(f"F-score vs {score_col} ({suffix})\nPearson r={r:.2f}, p={p:.2e}")
    plt.xlabel(score_col)
    plt.ylabel("F-score")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"fscore_vs_{score_col}_{suffix}.svg", format="svg")
    plt.close()

    merged[[label_key, score_col, "F-score"]].to_csv(
        Path(output_dir) / f"fscore_vs_{score_col}_{suffix}_correlation.csv", index=False
    )
