import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, linregress
import numpy as np

def plot_silhouette_summary(
    silhouette_score_path,
    output_dir,
    label,
    score_col,
    suffix="",
    show=False,
    fscore_path=None,
    mapping_path=None
):
    df = pd.read_csv(silhouette_score_path)
    grouped = df.groupby(label)[score_col].agg(['mean', 'std', 'median', 'count']).reset_index()
    grouped.rename(columns={'count': 'n_cells'}, inplace=True)

    fscore_df = None
    if fscore_path and mapping_path:
        fscore_df = pd.read_csv(fscore_path)
        mapping = pd.read_csv(mapping_path)
        mapping.columns = ['label', 'clusterName']
        fscore_df = fscore_df.merge(mapping, on='clusterName', how='inner')
        grouped = grouped.merge(fscore_df[['label', 'f_score']], on='label', how='left')

    grouped.sort_values(by='median', ascending=False, inplace=True)

    fig, ax = plt.subplots(figsize=(14, 8))
    bars = ax.bar(grouped[label], grouped['mean'], yerr=grouped['std'], capsize=4)

    for i, (mean, median, fscore) in enumerate(zip(grouped['mean'], grouped['median'], grouped.get('f_score', [None]*len(grouped)))):
        ax.text(i, mean + 0.02, f"Mean: {mean:.2f}\nMed: {median:.2f}", ha='center', fontsize=9)
        if fscore is not None:
            ax.text(i, mean + 0.15, f"F={fscore:.2f}", ha='center', fontsize=9, fontweight='bold')

    ax.set_title(f"Silhouette Score Summary for {label}")
    ax.set_ylabel("Mean Silhouette Score")
    ax.set_xlabel(label)
    ax.set_xticks(range(len(grouped)))
    ax.set_xticklabels(grouped[label], rotation=90)
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{label}_summary_{score_col}_{suffix}.png")
    fig.savefig(output_path)
    print(f"[✓] Saved silhouette summary plot to {output_path}")

    csv_path = os.path.join(output_dir, f"{label}_summary_{score_col}_{suffix}.csv")
    grouped.to_csv(csv_path, index=False)
    print(f"[✓] Summary stats exported to: {csv_path}")

    if show:
        plt.show()
    plt.close()

def plot_fscore_vs_silhouette_files(
    fscore_path,
    cluster_summary_path,
    output_dir,
    label,
    score_col,
    suffix="",
    show=False,
    export_csv=False,
    mapping_path=None,
    summary_path=None
):
    fscore_df = pd.read_csv(fscore_path)
    cluster_summary = pd.read_csv(cluster_summary_path)

    if mapping_path:
        mapping = pd.read_csv(mapping_path)
        mapping.columns = ['label', 'clusterName']
        fscore_df = fscore_df.merge(mapping, on='clusterName', how='inner')

    if summary_path:
        cluster_summary = pd.read_csv(summary_path)

    merged = pd.merge(fscore_df, cluster_summary, on='label', how='inner')

    if score_col not in merged.columns or 'f_score' not in merged.columns:
        raise ValueError("Missing required columns in merged data.")

    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 7))
    sns.scatterplot(data=merged, x='f_score', y=score_col, ax=ax)

    for _, row in merged.iterrows():
        ax.text(row['f_score'], row[score_col], row['label'], fontsize=8, alpha=0.7)

    r, _ = pearsonr(merged['f_score'], merged[score_col])
    slope, intercept, *_ = linregress(merged['f_score'], merged[score_col])
    x_vals = np.linspace(merged['f_score'].min(), merged['f_score'].max(), 100)
    y_vals = slope * x_vals + intercept
    ax.plot(x_vals, y_vals, linestyle='--', color='gray', label=f"r={r:.2f}")

    ax.set_title(f"F-score vs. {score_col} Correlation (r={r:.2f})")
    ax.set_xlabel("F-score")
    ax.set_ylabel(score_col)
    ax.legend()

    plot_path = os.path.join(output_dir, f"{label}_fscore_vs_{score_col}_{suffix}.png")
    fig.savefig(plot_path)
    print(f"[✓] Saved correlation plot to {plot_path}")

    if export_csv:
        merged.to_csv(os.path.join(output_dir, f"merged_fscore_silhouette_{suffix}.csv"), index=False)
        print(f"[✓] Exported merged data CSV.")

    if show:
        plt.show()
    plt.close()

