# src/scsilhouette/viz.py
import os
from pathlib import Path
from typing import Optional, List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import scanpy as sc

def plot_silhouette_summary(
    silhouette_score_path: str,
    output_dir: str,
    label: str,
    score_col: str,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
    show: bool = False,
    sort_by: str = "median",
):
    df = pd.read_csv(silhouette_score_path)
    prefix = Path(silhouette_score_path).stem
    suffix = prefix
    
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
        ax.text(i - 0.20, mean + std + 0.06, f"M={median:.2f}", ha="center", fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(grouped[label], rotation=90)
    ax.set_ylabel("Silhouette / F-score")
    ax.set_title(f"Silhouette Summary with F-score per {label}")
    ax.legend()

    fig.savefig(os.path.join(output_dir, f"fscore_silhouette_summary_{suffix}.png"), bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)

    # to do - add the cell count from the fscore file to the summary with the fscore,ppv, tn, fp, fn, tp to this file output
    
    grouped.to_csv(os.path.join(output_dir, f"fscore_silhouette_ummary_{prefix}.csv"), index=False)

def plot_correlation_summary(
    cluster_summary_path: str,
    output_dir: str,
    x_metric: str,
    y_metrics: List[str],
    label: str,
    show: bool = False,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
):
    df = pd.read_csv(cluster_summary_path)
    prefix = Path(cluster_summary_path).stem
    suffix = prefix

    # Merge fscore if provided
    if fscore_path:
        fscore_df = pd.read_csv(fscore_path)
        if mapping_path:
            mapping = pd.read_csv(mapping_path)
            mapping.columns = [label, "clusterName"]
            fscore_df = fscore_df.merge(mapping, on="clusterName", how="inner")
        df = df.merge(fscore_df[[label, "f_score"]], on=label, how="left")

    results = []

    for y_metric in y_metrics:
        if y_metric not in df.columns or x_metric not in df.columns:
            print(f"[WARN] Skipping: {x_metric} or {y_metric} not in dataframe columns")
            continue

        x = df[x_metric]
        y = df[y_metric]
        r, p = pearsonr(x, y)
        results.append({"x_metric": x_metric, "y_metric": y_metric, "pearson_r": r, "p_value": p})

        fig, ax = plt.subplots()
        sns.regplot(x=x, y=y, ax=ax)
        ax.set_title(f"{x_metric} vs {y_metric} (r={r:.2f}, p={p:.3f})")
        ax.set_xlabel(x_metric)
        ax.set_ylabel(y_metric)
        fig.tight_layout()

        fig.savefig(os.path.join(output_dir, f"{x_metric}_vs_{y_metric}_{suffix}.png"))
        if show:
            plt.show()
        plt.close(fig)

    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(output_dir, f"correlation_results_{suffix}.csv"), index=False)

def plot_dotplot(
    h5ad_path: str,
    groupby: str,
    embedding_key: str,
    output_dir: str,
    show: bool = False,
):
    adata = sc.read_h5ad(h5ad_path)
    prefix = Path(h5ad_path).stem
    suffix = prefix
    
    sc.pl.embedding(
        adata,
        basis=embedding_key.replace("X_", ""),
        color=groupby,
        show=show,
        save=f"dotplot_{suffix}_{embedding_key}.png",
    )


def plot_heatmap(
    h5ad_path: str,
    groupby: str,
    embedding_key: str,
    output_dir: str,
    suffix: str = "",
    show: bool = False,
):
    adata = sc.read_h5ad(h5ad_path)
    prefix = Path(h5ad_path).stem
    suffix = prefix

    sc.pl.embedding(
        adata,
        basis=embedding_key.replace("X_", ""),
        color=groupby,
        cmap="viridis",
        show=show,
        save=f"heatmap_{suffix}_{embedding_key}.png",
    )


def plot_distribution(
    summary_csv: str,
    label_key: str,
):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    df = pd.read_csv(summary_csv)
    prefix = Path(summary_csv).stem

    df['count_log10'] = np.log10(df['count'].replace(0, np.nan))

    # === Plot 1: Sorted by log10(cell count) ===
    df_sorted_log = df.sort_values("count_log10", ascending=False).reset_index(drop=True)
    fig_log, ax1 = plt.subplots(figsize=(18, 6))

    ax1.bar(df_sorted_log[label_key], df_sorted_log["count_log10"], color="steelblue", label="log10(Cell Count)")
    ax1.set_ylabel("log10(Cell Count)", color="steelblue")
    ax1.tick_params(axis='y', labelcolor="steelblue")
    ax1.tick_params(axis='x', rotation=90)

    ax2 = ax1.twinx()
    ax2.plot(df_sorted_log[label_key], df_sorted_log["mean_silhouette"], color="blue", marker='o', label="Mean")
    ax2.plot(df_sorted_log[label_key], df_sorted_log["median_silhouette"], color="red", marker='o', label="Median")
    ax2.plot(df_sorted_log[label_key], df_sorted_log["std_silhouette"], color="black", marker='o', label="Std Dev")
    ax2.set_ylabel("Silhouette Scores", color="black")
    ax2.tick_params(axis='y', labelcolor="black")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

    plt.title("Cluster Summary: log10(Cell Count) and Silhouette Metrics")
    plt.tight_layout()
    fig_log.savefig(f"{prefix}_log10.png")
    plt.close(fig_log)

    # === Plot 2: Sorted by raw count ===
    df_sorted_raw = df.sort_values("count", ascending=False).reset_index(drop=True)
    fig_raw, ax1_raw = plt.subplots(figsize=(18, 6))

    ax1_raw.bar(df_sorted_raw[label_key], df_sorted_raw["count"], color="steelblue", label="Cell Count")
    ax1_raw.set_ylabel("Cell Count", color="steelblue")
    ax1_raw.tick_params(axis='y', labelcolor="steelblue")
    ax1_raw.tick_params(axis='x', rotation=90)

    ax2_raw = ax1_raw.twinx()
    ax2_raw.plot(df_sorted_raw[label_key], df_sorted_raw["mean_silhouette"], color="blue", marker='o', label="Mean")
    ax2_raw.plot(df_sorted_raw[label_key], df_sorted_raw["median_silhouette"], color="red", marker='o', label="Median")
    ax2_raw.plot(df_sorted_raw[label_key], df_sorted_raw["std_silhouette"], color="black", marker='o', label="Std Dev")
    ax2_raw.set_ylabel("Silhouette Scores", color="black")
    ax2_raw.tick_params(axis='y', labelcolor="black")

    lines1, labels1 = ax1_raw.get_legend_handles_labels()
    lines2, labels2 = ax2_raw.get_legend_handles_labels()
    ax2_raw.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

    plt.title("Cluster Summary: Raw Cell Count and Silhouette Metrics")
    plt.tight_layout()
    fig_raw.savefig(f"{prefix}_raw.png")
    plt.close(fig_raw)

    print(f"[DONE] Saved: {prefix}_log10.png and {prefix}_raw.png")
def plot_dataset_summary(
    cluster_summary_path: str,
    label: str,
    score_col: str = "silhouette_score",
    output_dir: str = "results",
    show: bool = False
):
    df = pd.read_csv(cluster_summary_path)
    prefix = Path(cluster_summary_path).stem

    df["quartile"] = pd.qcut(df["mean_silhouette"], 4, labels=["Q1", "Q2", "Q3", "Q4"])

    for q in ["Q1", "Q2", "Q3", "Q4"]:
        subset = df[df["quartile"] == q]
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.boxplot(data=subset[["mean_silhouette", "std_silhouette", "median_silhouette", "count"]], ax=ax)
        ax.set_title(f"{prefix} Cluster Summary: {q}")
        fig.tight_layout()

        fig_path = Path(output_dir) / f"{prefix}_summary_{q}.png"
        fig.savefig(fig_path)
        if show:
            plt.show()
        plt.close(fig)

