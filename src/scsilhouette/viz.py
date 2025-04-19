# src/scsilhouette/viz.py
import os
from pathlib import Path
from typing import Optional
from typing import List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr


def plot_silhouette_summary(
    silhouette_score_path: str,
    output_dir: str,
    label: str,
    score_col: str,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
    suffix: str = "",
    show: bool = False,
    sort_by: str = "median",
):
    df = pd.read_csv(silhouette_score_path)
    grouped = df.groupby(label)[score_col].agg(["mean", "std", "median", "count"]).reset_index()

    if fscore_path and mapping_path:
        fscore_df = pd.read_csv(fscore_path)
        mapping = pd.read_csv(mapping_path)
        mapping.columns = [label, "clusterName"]
        fscore_df = fscore_df.merge(mapping, on="clusterName", how="inner")
        grouped = grouped.merge(fscore_df[[label, "f_score"]], on=label, how="left")

    if sort_by not in ["mean", "median", "std"]:
        raise ValueError("--sort-by must be one of 'mean', 'median', or 'std'")

    grouped = grouped.sort_values(sort_by, ascending=False)

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(grouped))
    ax.bar(x - 0.2, grouped["median"], width=0.4, label="Median Silhouette")

    if "f_score" in grouped.columns:
        ax.bar(x + 0.2, grouped["f_score"], width=0.4, label="F-score")
        for i, f in enumerate(grouped["f_score"]):
            ax.text(i + 0.2, f + 0.02, f"F={f:.2f}", ha="center", fontsize=6)

    for i, (mean, std, median) in enumerate(zip(grouped["mean"], grouped["std"], grouped["median"])):
        ax.errorbar(i - 0.2, mean, yerr=std, fmt='o', color='black', capsize=5)
        ax.scatter(i - 0.2, median, color='red', zorder=5)
        ax.text(i - 0.25, mean + std + 0.02, f"\u03BC={mean:.2f}", ha="center", fontsize=6)
        ax.text(i - 0.20, mean + std + 0.06, f"M={median:.2f}"   , ha="center", fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(grouped[label], rotation=90)
    ax.set_ylabel("Silhouette / F-score")
    ax.set_title(f"Silhouette Summary with F-score per {label}")
    ax.legend()

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(output_dir, f"{label}_summary_{score_col}{suffix}.png"), bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)

    grouped.to_csv(os.path.join(output_dir, f"{label}_summary_{score_col}{suffix}.csv"), index=False)

def plot_correlation_summary(
    cluster_summary_path: str,
    output_dir: str,
    silhouette_metrics: List[str],
    score_col: str,
    label: str,
    suffix: str = "",
    show: bool = False,
):
    df = pd.read_csv(cluster_summary_path)
    results = []

    for metric in silhouette_metrics:
        if metric not in df.columns or score_col not in df.columns:
            continue

        x = df[score_col]
        y = df[metric]
        r, p = pearsonr(x, y)
        results.append({"metric": metric, "pearson_r": r, "p_value": p})

        fig, ax = plt.subplots()
        sns.regplot(x=x, y=y, ax=ax)
        ax.set_title(f"{score_col} vs {metric} (r={r:.2f}, p={p:.3f})")
        ax.set_xlabel(score_col)
        ax.set_ylabel(metric)
        fig.tight_layout()

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        fig.savefig(os.path.join(output_dir, f"{label}_fscore_vs_{metric}{suffix}.png"))
        if show:
            plt.show()
        plt.close(fig)

    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(output_dir, f"{label}_fscore_correlations{suffix}.csv"), index=False)
    
