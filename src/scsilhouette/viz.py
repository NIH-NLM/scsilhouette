# src/scsilhouette/viz.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from typing import Optional
import numpy as np


def plot_silhouette_summary(
    silhouette_score_path,
    output_dir,
    label_col,
    score_col,
    suffix="",
    show=False,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
):
    df = pd.read_csv(silhouette_score_path)
    summary_df = df.groupby(label_col)[score_col].agg(['mean', 'median', 'std', 'count']).reset_index()
    summary_df = summary_df.sort_values(by="median", ascending=False)

    if fscore_path:
        fscore_df = pd.read_csv(fscore_path)
        if mapping_path:
            mapping_df = pd.read_csv(mapping_path)
            fscore_df = pd.merge(fscore_df, mapping_df, on="clusterName", how="left")
        fscore_cols = [col for col in ["F-score", "f_score", "fbeta", "PPV", "recall", "TN", "FP", "FN", "TP"] if col in fscore_df.columns]
        if 'cell_type' not in fscore_df.columns:
            raise ValueError("Expected 'cell_type' column in fscore data.")
        summary_df = pd.merge(summary_df, fscore_df[['cell_type'] + fscore_cols].drop_duplicates(), on="cell_type", how="left")
    else:
        fscore_cols = []

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    plot_path = Path(output_dir) / f"{label_col}_summary_{score_col}_{suffix}.png"
    stats_path = Path(output_dir) / f"{label_col}_summary_{score_col}_{suffix}.csv"

    fig, ax = plt.subplots(figsize=(16, 6))
    x = summary_df[label_col]
    mean_vals = summary_df["mean"]
    std_vals = summary_df["std"]
    medians = summary_df["median"]
    x_pos = np.arange(len(summary_df))
    bar_width = 0.4

    ax.bar(x_pos - bar_width/2, mean_vals, yerr=std_vals, alpha=0.7, align='center', ecolor='black', capsize=5, width=bar_width, label='Mean ± Std')
    ax.scatter(x_pos - bar_width/2, medians, color='red', zorder=5, label='Median')

    if "f_score" in summary_df.columns:
        ax.bar(x_pos + bar_width/2, summary_df["f_score"], width=bar_width, alpha=0.7, label='F-score')

    for i, row in summary_df.iterrows():
        text = f"μ={row['mean']:.2f}\nM={row['median']:.2f}"
        if "f_score" in row:
            text += f"\nF={row['f_score']:.2f}"
        ax.text(i, max(row['mean'] + row['std'], row.get("f_score", 0)) + 0.02, text, ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x, rotation=90)
    ax.set_title(f"Silhouette Summary with F-score per {label_col}")
    ax.legend()
    fig.tight_layout()

    plt.savefig(plot_path)
    if show:
        plt.show()
    plt.close()

    summary_df.to_csv(stats_path, index=False)
    print(f"Summary stats exported to: {stats_path}")

def plot_fscore_vs_silhouette_files(
    fscore_path,
    cluster_summary_path,
    output_dir,
    label,
    score_col,
    suffix="",
    show=False,
    export_csv=False,
    mapping_path: Optional[str] = None,
    summary_path: Optional[str] = None,
):
    fscore_df = pd.read_csv(fscore_path)
    summary_df = pd.read_csv(summary_path) if summary_path else pd.read_csv(cluster_summary_path)

    if mapping_path:
        mapping_df = pd.read_csv(mapping_path)
        if 'clusterName' not in mapping_df.columns or 'cell_type' not in mapping_df.columns:
            raise ValueError("Mapping file must contain 'clusterName' and 'cell_type' columns.")
        original_cluster_count = fscore_df['clusterName'].nunique()
        fscore_df = pd.merge(fscore_df, mapping_df, on="clusterName", how="left")
        unmatched = fscore_df[fscore_df['cell_type'].isna()]["clusterName"].unique()
        if len(unmatched) > 0:
            warning_path = Path(output_dir) / f"unmatched_clusters_{suffix}.csv"
            pd.Series(unmatched, name="unmatched_clusterName").to_csv(warning_path, index=False)
            print(f"[WARNING] {len(unmatched)} cluster(s) in fscore file not found in mapping file. Saved to: {warning_path}")
        fscore_df = fscore_df.dropna(subset=['cell_type'])
        print(f"Matched {fscore_df['cell_type'].nunique()} cell types from {original_cluster_count} clusters.")

    if label not in summary_df.columns:
        raise ValueError(f"'{label}' not found in summary. Found: {summary_df.columns.tolist()}")

    fscore_col = next((col for col in ["F-score", "f_score", "fbeta"] if col in fscore_df.columns), None)
    if not fscore_col:
        raise ValueError("No recognized f-score column found in fscore file.")

    score_source_col = score_col if score_col in summary_df.columns else None
    if not score_source_col:
        raise ValueError(f"'{score_col}' not found in summary data.")

    merged = pd.merge(fscore_df, summary_df, how="inner", left_on="cell_type", right_on=label)
    if merged.empty:
        raise ValueError("Merge resulted in 0 rows. Check if mapped cell types match summary labels.")

    corr = merged[[fscore_col, score_source_col]].corr().iloc[0, 1]
    print(f"Correlation ({fscore_col} vs {score_source_col}) = {corr:.3f}")

    merged["fscore_bin"] = pd.qcut(merged[fscore_col], q=4, labels=["Q1", "Q2", "Q3", "Q4"])
    if "median" in merged.columns:
        cluster_order = merged.groupby("cell_type")["median"].median().sort_values(ascending=False).index
        merged["cell_type"] = pd.Categorical(merged["cell_type"], categories=cluster_order, ordered=True)

    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=merged, x=fscore_col, y=score_source_col, hue="fscore_bin", palette="viridis")

    x_vals = merged[[fscore_col]].values
    y_vals = merged[[score_source_col]].values
    model = LinearRegression().fit(x_vals, y_vals)
    pred_y = model.predict(x_vals)
    plt.plot(merged[fscore_col], pred_y, color='black', linestyle='--', linewidth=1, label='Regression Line')

    grouped = merged.groupby("cell_type")[[fscore_col, score_source_col]].median().reset_index()
    for _, row in grouped.iterrows():
        plt.text(row[fscore_col], row[score_source_col] + 0.02, f"{row[fscore_col]:.2f}", ha='center', va='bottom', fontsize=8, color='black')

    plt.title(f"{fscore_col} vs {score_source_col} (r={corr:.2f})")
    plt.tight_layout()
    plt.legend()

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    plot_path = Path(output_dir) / f"{label}_fscore_vs_{score_source_col}_{suffix}.png"
    plt.savefig(plot_path)
    if show:
        plt.show()
    plt.close()

    if export_csv:
        csv_path = Path(output_dir) / f"{label}_fscore_vs_{score_source_col}_{suffix}_merged.csv"
        merged.to_csv(csv_path, index=False)
        print(f"Merged data exported to: {csv_path}")

        stats_df = merged.groupby("cell_type")[[fscore_col, score_source_col]].agg(['mean', 'std', 'median', 'count']).reset_index()
        stats_path = Path(output_dir) / f"{label}_fscore_vs_{score_source_col}_{suffix}_stats.csv"
        stats_df.to_csv(stats_path, index=False)
        print(f"Stats per cluster exported to: {stats_path}")

