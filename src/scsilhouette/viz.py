
# src/scsilhouette/viz.py

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def plot_cluster_summary(cluster_summary_path, output_dir, label, score_col, suffix=""):
    df = pd.read_csv(cluster_summary_path)
    if label not in df.columns or score_col not in df.columns:
        raise ValueError(f"'{label}' or '{score_col}' not found in the cluster summary file.")
    ax = df.plot.bar(x=label, y=score_col, legend=False)
    ax.set_ylabel(score_col)
    ax.set_title(f"{score_col} per {label}")
    output_file = Path(output_dir) / f"{label}_cluster_summary_{suffix}.svg"
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_score_distribution(cell_scores_path, output_dir, label, score_col, suffix=""):
    df = pd.read_csv(cell_scores_path)
    if score_col not in df.columns:
        raise ValueError(f"{score_col} not found in cell scores.")
    ax = df[score_col].hist(bins=50)
    ax.set_xlabel(score_col)
    ax.set_ylabel("Frequency")
    ax.set_title(f"Distribution of {score_col}")
    output_file = Path(output_dir) / f"{label}_score_distribution_{suffix}.svg"
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_all(cluster_summary_path, cell_scores_path, output_dir, label, score_col, suffix=""):
    plot_cluster_summary(cluster_summary_path, output_dir, label, score_col, suffix)
    plot_score_distribution(cell_scores_path, output_dir, label, score_col, suffix)
