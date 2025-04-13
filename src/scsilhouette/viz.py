import pandas as pd
import os
import matplotlib.pyplot as plt

def plot_cluster_summary(cluster_summary_path, output_dir, label, score_col, suffix=""):
    df = pd.read_csv(cluster_summary_path)
    if label not in df.columns or score_col not in df.columns:
        raise ValueError(f"'{label}' or '{score_col}' not found in the cluster summary file")

    ax = df.plot.bar(x=label, y=score_col, legend=False)
    fig_path = os.path.join(output_dir, f"{label}_summary_{suffix}.svg")
    plt.tight_layout()
    plt.savefig(fig_path)
    plt.close()
    print(f"[✓] Saved: {fig_path}")

def plot_score_distribution(cell_scores_path, output_dir, label, score_col, suffix=""):
    df = pd.read_csv(cell_scores_path)
    ax = df.hist(column=score_col, by=label, figsize=(12, 10))
    fig_path = os.path.join(output_dir, f"{label}_score_dist_{suffix}.svg")
    plt.tight_layout()
    plt.savefig(fig_path)
    plt.close()
    print(f"[✓] Saved: {fig_path}")

def plot_all(cluster_summary_path, cell_scores_path, output_dir, label, score_col, suffix=""):
    plot_cluster_summary(cluster_summary_path, output_dir, label, score_col, suffix)
    plot_score_distribution(cell_scores_path, output_dir, label, score_col, suffix)

