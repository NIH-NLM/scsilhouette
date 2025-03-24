# src/scsilhouette/viz.py

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_score_distribution(cell_scores: pd.DataFrame, output_dir: Path, label: str):
    plt.figure()
    cell_scores['silhouette_score'].hist(bins=50)
    plt.title(f'Silhouette Score Distribution - {label}')
    plt.xlabel('Silhouette Score')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(output_dir / f"{label}_score_distribution.png")
    plt.close()


def plot_cluster_summary(cluster_summary: pd.DataFrame, output_dir: Path, label: str):
    plt.figure()
    cluster_summary.plot.bar(
        x=label, y="mean_silhouette_score", legend=False
    )
    plt.title(f'Mean Silhouette Score by {label}')
    plt.ylabel('Mean Silhouette Score')
    plt.tight_layout()
    plt.savefig(output_dir / f"{label}_cluster_summary.png")
    plt.close()


def plot_cluster_size_vs_score(cluster_summary: pd.DataFrame, output_dir: Path, label: str):
    plt.figure()
    plt.scatter(
        cluster_summary["n_cells"],
        cluster_summary["mean_silhouette_score"],
        alpha=0.7
    )
    plt.xlabel("Cluster Size (n_cells)")
    plt.ylabel("Mean Silhouette Score")
    plt.title(f"Cluster Size vs Silhouette Score - {label}")
    plt.tight_layout()
    plt.savefig(output_dir / f"{label}_cluster_size_vs_score.png")
    plt.close()


def plot_qc_boxplots(cell_scores: pd.DataFrame, obs: pd.DataFrame, output_dir: Path, label: str):
    merged = pd.concat([cell_scores, obs[["nCount_RNA", "nFeature_RNA"]]], axis=1)

    for feature in ["nCount_RNA", "nFeature_RNA"]:
        plt.figure(figsize=(10, 6))
        sns.boxplot(x=label, y=feature, data=merged, showfliers=False)
        sns.stripplot(x=label, y=feature, data=merged, color='black', alpha=0.3, jitter=0.2)
        plt.title(f"{feature} by {label}")
        plt.ylabel(feature)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(output_dir / f"{label}_{feature}_by_cluster.png")
        plt.close()

        # Correlation scatter
        plt.figure(figsize=(6, 5))
        sns.regplot(x=feature, y="silhouette_score", data=merged)
        plt.title(f"Correlation: {feature} vs Silhouette Score")
        plt.tight_layout()
        plt.savefig(output_dir / f"{label}_corr_{feature}_vs_score.png")
        plt.close()

def plot_all(cell_scores: pd.DataFrame, cluster_summary: pd.DataFrame, output_dir: Path, label: str):
    plot_score_distribution(cell_scores, output_dir, label)
    plot_cluster_summary(cluster_summary, output_dir, label)
    plot_cluster_size_vs_score(cluster_summary, output_dir, label)
