# src/scsilhouette/viz.py

import os
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
):
    df = pd.read_csv(silhouette_score_path)
    grouped = df.groupby(label)[score_col].agg(["mean", "std", "median", "count"]).reset_index()

    if fscore_path and mapping_path:
        fscore_df = pd.read_csv(fscore_path)
        mapping_df = pd.read_csv(mapping_path)
        mapping_df.columns = ["clusterName", label]

        fscore_df = fscore_df.merge(mapping_df, on="clusterName", how="inner")
        fscore_df = fscore_df[[label, "f_score"]].drop_duplicates()

        grouped = grouped.merge(fscore_df, on=label, how="left")

    grouped = grouped.sort_values("median", ascending=False)

    fig, ax = plt.subplots(figsize=(20, 6))
    x = np.arange(len(grouped))

    ax.bar(x - 0.2, grouped["median"], width=0.4, label="Median Silhouette", color="skyblue")
    if "f_score" in grouped.columns:
        ax.bar(x + 0.2, grouped["f_score"], width=0.4, label="F-score", color="orange")

    ax.errorbar(x - 0.2, grouped["mean"], yerr=grouped["std"], fmt='o', color='black', capsize=5)
    ax.scatter(x - 0.2, grouped["median"], color='red', zorder=5)

    for i, row in grouped.iterrows():
        # Mean & Median for silhouette
        ax.text(i - 0.2, row["mean"] + row["std"] + 0.02, f"μ={row['mean']:.2f}", ha="center", fontsize=7)
        ax.text(i - 0.2, row["median"] + 0.02, f"M={row['median']:.2f}", ha="center", fontsize=7)
        # F-score if available
        if "f_score" in row and not pd.isna(row["f_score"]):
            ax.text(i + 0.2, row["f_score"] + 0.02, f"F={row['f_score']:.2f}", ha="center", fontsize=7)

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


def plot_fscore_vs_silhouette_files(
    fscore_path: str,
    cluster_summary_path: str,
    output_dir: str,
    label: str,
    score_col: str,
    suffix: str = "",
    mapping_path: Optional[str] = None,
    summary_path: Optional[str] = None,
    silhouette_stat: str = "mean",
    show: bool = False,
    export_csv: bool = False,
):
    fscore_df = pd.read_csv(fscore_path)
    cluster_summary = pd.read_csv(summary_path) if summary_path else pd.read_csv(cluster_summary_path)

    if mapping_path:
        mapping_df = pd.read_csv(mapping_path)
        mapping_df.columns = ["clusterName", label]
        fscore_df = fscore_df.merge(mapping_df, on="clusterName", how="inner")

    fscore_df = fscore_df[[label, "f_score"]].drop_duplicates()
    merged = pd.merge(fscore_df, cluster_summary, on=label, how="left")

    y_col = f"{silhouette_stat}_{score_col.split('_')[-1]}"
    if y_col not in merged.columns:
        raise ValueError(f"{y_col} not found in merged DataFrame.")

    merged = merged.dropna(subset=["f_score", y_col])
    if merged.empty:
        raise ValueError("No data left after merging. Check mapping or column names.")

    r, _ = pearsonr(merged["f_score"], merged[y_col])

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.regplot(data=merged, x="f_score", y=y_col, ax=ax)
    ax.set_title(f"f_score vs {y_col} (r={r:.2f})")
    ax.axhline(merged[y_col].median(), linestyle="--", color="gray")
    ax.axvline(merged["f_score"].median(), linestyle="--", color="gray")
    ax.set_xlabel("f_score")
    ax.set_ylabel(y_col)

    for _, row in merged.iterrows():
        ax.text(row["f_score"], row[y_col], f"{row['f_score']:.2f}", fontsize=7)

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(output_dir, f"{label}_fscore_vs_{y_col}_{suffix}.png"), bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)

    if export_csv:
        merged.to_csv(
            os.path.join(output_dir, f"{label}_merged_fscore_silhouette_{suffix}.csv"),
            index=False,
        )

