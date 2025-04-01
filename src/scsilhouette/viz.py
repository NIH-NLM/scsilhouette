# src/scsilhouette/viz.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr


def plot_score_distribution(cell_scores, output_dir, label, score_col, suffix=""):
    plt.figure(figsize=(8, 4))
    sns.histplot(cell_scores[score_col], bins=50, kde=True)
    plt.title(f"Silhouette Score Distribution - {label} ({suffix})")
    plt.xlabel(score_col)
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_{suffix}_score_distribution.png")
    plt.close()


def plot_cluster_summary(cluster_summary, output_dir, label, score_col, suffix=""):
    cluster_summary.plot.bar(x=label, y=score_col, legend=False)
    plt.title(f"Mean Silhouette Score by {label} ({suffix})")
    plt.ylabel("Mean Silhouette Score")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_{suffix}_cluster_summary.png")
    plt.close()


def plot_cluster_size_vs_score(cluster_summary, output_dir, label, score_col, suffix=""):
    plt.figure(figsize=(6, 6))
    plt.scatter(cluster_summary["n_cells"], cluster_summary[score_col], alpha=0.7)
    plt.xlabel("Cluster Size (n_cells)")
    plt.ylabel("Mean Silhouette Score")
    plt.title(f"{label} ({suffix}): Cluster Size vs Silhouette Score")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{label}_{suffix}_cluster_size_vs_score.png")
    plt.close()


def plot_qc_boxplots(cell_scores, output_dir, label, suffix=""):
    for metric in ["nCount_RNA", "nFeature_RNA"]:
        if metric not in cell_scores:
            continue
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=cell_scores, x=label, y=metric)
        sns.stripplot(data=cell_scores, x=label, y=metric, color="black", alpha=0.3, jitter=0.2)
        plt.xticks(rotation=90)
        plt.title(f"{metric} by {label} ({suffix})")
        plt.tight_layout()
        plt.savefig(Path(output_dir) / f"{label}_{suffix}_{metric}_boxplot.png")
        plt.close()


def plot_all(cluster_summary, cell_scores, output_dir, label, score_col, suffix=""):
    plot_score_distribution(cell_scores, output_dir, label, score_col, suffix)
    plot_cluster_summary(cluster_summary, output_dir, label, score_col, suffix)
    plot_cluster_size_vs_score(cluster_summary, output_dir, label, score_col, suffix)
    plot_qc_boxplots(cell_scores, output_dir, label, suffix)

