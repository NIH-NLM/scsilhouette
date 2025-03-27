# src/scsilhouette/viz.py

import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import seaborn as sns

def plot_score_distribution(cell_scores: pd.DataFrame, output_dir: str, label: str):
    plt.figure(figsize=(8, 4))
    sns.histplot(cell_scores["silhouette_score"], bins=50, kde=True)
    plt.title(f"Silhouette Score Distribution - {label}")
    plt.xlabel("Silhouette Score")
    plt.tight_layout()
    path = Path(output_dir) / f"{label}_score_distribution.png"
    plt.savefig(path)
    plt.close()

def plot_cluster_summary(cluster_summary: pd.DataFrame, output_dir: str, label: str):
    cluster_summary.plot.bar(x=label, y="mean_silhouette_score", legend=False)
    plt.title(f"Mean Silhouette Score by {label}")
    plt.ylabel("Mean Silhouette Score")
    plt.tight_layout()
    path = Path(output_dir) / f"{label}_cluster_summary.png"
    plt.savefig(path)
    plt.close()

def plot_qc_boxplots(cell_scores: pd.DataFrame, output_dir: str, label: str):
    for qc_metric in ["nCount_RNA", "nFeature_RNA"]:
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=cell_scores, x=label, y=qc_metric)
        sns.stripplot(data=cell_scores, x=label, y=qc_metric, color="black", alpha=0.3, jitter=0.2)
        plt.xticks(rotation=90)
        plt.title(f"{qc_metric} by {label}")
        plt.tight_layout()
        path = Path(output_dir) / f"{label}_{qc_metric}_boxplot.png"
        plt.savefig(path)
        plt.close()

def plot_cluster_size_vs_score(cluster_summary: pd.DataFrame, output_dir: str, label: str):
    plt.figure(figsize=(6, 6))
    plt.scatter(
        cluster_summary["n_cells"],
        cluster_summary["mean_silhouette_score"],
        alpha=0.7
    )
    plt.xlabel("Cluster Size (n_cells)")
    plt.ylabel("Mean Silhouette Score")
    plt.title(f"{label}: Cluster Size vs Silhouette Score")
    plt.tight_layout()
    path = Path(output_dir) / f"{label}_cluster_size_vs_score.png"
    plt.savefig(path)
    plt.close()

def plot_all(cluster_summary: pd.DataFrame, cell_scores: pd.DataFrame, output_dir: str, label: str):
    plot_score_distribution(cell_scores, output_dir, label)
    plot_cluster_summary(cluster_summary, output_dir, label)
    plot_cluster_size_vs_score(cluster_summary, output_dir, label)
    plot_qc_boxplots(cell_scores, output_dir, label)

def plot_fscore_vs_silhouette(merged_df: pd.DataFrame, output_dir: str):
    plt.figure(figsize=(6, 6))
    sns.scatterplot(
        data=merged_df,
        x="mean_silhouette_score",
        y="F-score",
        hue="cluster",
        palette="tab10"
    )
    plt.xlabel("Mean Silhouette Score")
    plt.ylabel("NSForest F-score")
    plt.title("F-score vs. Silhouette Score")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(Path(output_dir) / "fscore_vs_silhouette.png")
    plt.close()

