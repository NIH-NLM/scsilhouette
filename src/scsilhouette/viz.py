# src/scsilhouette/viz.py

import pandas as pd
import matplotlib.pyplot as plt
import os


def plot_cluster_summary(cluster_summary_path, output_dir, label, score_col, suffix=""):
    df = pd.read_csv(cluster_summary_path)
    if label not in df.columns or score_col not in df.columns:
        raise ValueError(f"'{label}' or '{score_col}' not found in {cluster_summary_path}")
    ax = df.plot.bar(x=label, y=score_col, legend=False)
    plt.ylabel(score_col.replace("_", " ").capitalize())
    plt.title(f"{label} cluster summary")
    filename = f"{label}_cluster_summary_{suffix}.svg"
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()


def plot_score_distribution(cell_scores_path, output_dir, label, score_col, suffix=""):
    df = pd.read_csv(cell_scores_path)
    if label not in df.columns or score_col not in df.columns:
        raise ValueError(f"'{label}' or '{score_col}' not found in {cell_scores_path}")
    df.boxplot(column=score_col, by=label, grid=False, rot=90)
    plt.title(f"{score_col} distribution by {label}")
    plt.suptitle("")
    plt.xlabel(label)
    plt.ylabel(score_col)
    filename = f"{label}_score_distribution_{suffix}.svg"
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()


def plot_fscore_vs_silhouette(fscore_path, cluster_summary_path, output_dir, label, score_col, suffix=""):
    fscore_df = pd.read_csv(fscore_path)
    summary_df = pd.read_csv(cluster_summary_path)

    if "clusterName" not in fscore_df.columns or "f_score" not in fscore_df.columns:
        raise ValueError("fscore file must contain 'clusterName' and 'f_score' columns")

    merged_df = pd.merge(
        fscore_df[["clusterName", "f_score"]],
        summary_df[[label, score_col]],
        left_on="clusterName",
        right_on=label,
        how="inner"
    )

    if merged_df.empty:
        raise ValueError("Merged DataFrame is empty. Check that cluster names and labels match.")

    ax = merged_df.plot.scatter(x="f_score", y=score_col)
    plt.xlabel("f-score")
    plt.ylabel(score_col)
    plt.title(f"f-score vs {score_col}")
    filename = f"{label}_fscore_vs_silhouette_{suffix}.svg"
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()

