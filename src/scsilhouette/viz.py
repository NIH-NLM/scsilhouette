# src/scsilhouette/viz.py (Plotly: auto-export to HTML, SVG, PNG for Nextflow)
import os
from pathlib import Path
from typing import Optional, List
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio


def plot_silhouette_summary(
    silhouette_score_path: str,
    label: str,
    silhouette_score_col: str,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
    show: bool = False,
    sort_by: str = "median",
):
    df = pd.read_csv(silhouette_score_path)
    prefix = Path(silhouette_score_path).stem
    suffix = prefix

    grouped = df.groupby(label)[silhouette_score_col].agg([
        "mean", "std", "median", "count", "min", "max",
        lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75)
    ]).reset_index()
    grouped.columns = [label, "mean", "std", "median", "count", "min", "max", "q1", "q3"]

    has_fscore = False
    if fscore_path and mapping_path:
        fscore_df = pd.read_csv(fscore_path)
        mapping = pd.read_csv(mapping_path)
        mapping.columns = [label, "clusterName"]
        fscore_df = fscore_df.merge(mapping, on="clusterName", how="inner")
        grouped = grouped.merge(fscore_df[[label, "f_score"]], on=label, how="inner")
        has_fscore = True

    grouped = grouped.sort_values(by=sort_by, ascending=False)
    x = grouped[label].tolist()

    fig = go.Figure()

    # IQR error bar
    fig.add_trace(go.Bar(
        x=x,
        y=grouped["median"],
        name="Median Silhouette",
        marker_color="skyblue",
        error_y=dict(
            type='data',
            symmetric=False,
            array=grouped["q3"] - grouped["median"],
            arrayminus=grouped["median"] - grouped["q1"],
            thickness=1.5,
            color='black'
        )
    ))

    # Mean dots
    fig.add_trace(go.Scatter(
        x=x,
        y=grouped["mean"],
        mode="markers",
        name="Mean",
        marker=dict(color="red", size=8, symbol="circle")
    ))

    # Optional F-score bars
    if has_fscore:
        fig.add_trace(go.Bar(
            x=x,
            y=grouped["f_score"],
            name="F-score",
            marker_color="orange",
            opacity=0.6
        ))

    # Global + cluster stats
    global_mean = df[silhouette_score_col].mean()
    global_median = df[silhouette_score_col].median()
    global_std = df[silhouette_score_col].std()
    cluster_mean_of_means = grouped["mean"].mean()
    cluster_median_of_medians = grouped["median"].median()
    cluster_std_of_stds = grouped["std"].mean()

    stats_text = (
        f"Global:\nMean={global_mean:.2f}, Median={global_median:.2f}, Std={global_std:.2f}\n"
        f"Cluster:\nMean of Means={cluster_mean_of_means:.2f},\n"
        f"Median of Medians={cluster_median_of_medians:.2f}, Mean of Stds={cluster_std_of_stds:.2f}"
    )

    fig.add_annotation(
        text=stats_text,
        align="left",
        showarrow=False,
        xref="paper", yref="paper",
        x=0.01, y=0.01,
        bordercolor="black",
        borderwidth=1,
        bgcolor="white",
        opacity=0.9,
        font=dict(size=12)
    )

    fig.update_layout(
        title=f"Silhouette Summary per {label}",
        xaxis=dict(title=label, tickangle=45),
        yaxis=dict(title="Silhouette Score"),
        barmode="group",
        legend=dict(x=1.01, y=1),
        margin=dict(l=50, r=50, t=80, b=150),
        height=600,
    )

    # Output files
    outbase = f"{suffix}_silhouette_summary"
    html_path = f"{outbase}.html"
    svg_path = f"{outbase}.svg"
    png_path = f"{outbase}.png"

    fig.write_html(html_path)
    pio.write_image(fig, svg_path, format='svg')
    pio.write_image(fig, png_path, format='png')

    print(f"Saved interactive HTML to {html_path}")
    print(f"Saved static SVG to {svg_path}")
    print(f"Saved PNG image to {png_path}")

    if show:
        fig.show()

def plot_correlation_summary(
    cluster_summary_path: str,
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

        fig.savefig(f"{x_metric}_vs_{y_metric}_{suffix}.png")
        if show:
            plt.show()
        plt.close(fig)

    results_df = pd.DataFrame(results)
    results_df.to_csv(f"correlation_results_{suffix}.csv", index=False)

def plot_dotplot(
    h5ad_path: str,
    label_key: str,
    embedding_key: str,
    show: bool = False,
):
    adata = sc.read_h5ad(h5ad_path)
    prefix = Path(h5ad_path).stem
    suffix = prefix
    
    sc.pl.embedding(
        adata,
        basis=embedding_key.replace("X_", ""),
        color=label_key,
        show=show,
        save=f"dotplot_{suffix}_{embedding_key}.png",
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
    show: bool = False
):
    df = pd.read_csv(cluster_summary_path)

    prefix = Path(cluster_summary_path).stem

    metrics = {
        "mean_silhouette": "Mean Silhouette",
        "median_silhouette": "Median Silhouette",
        "std_silhouette": "Std Dev Silhouette",
        "log10_count": "log10(Cluster Size)"
    }

    df["log10_count"] = np.log10(df["count"] + 1)

    plots = []
    ylims = {}

    # Get global y-limits per row
    for metric in metrics:
        min_val = df[metric].min()
        max_val = df[metric].max()
        ylims[metric] = (min_val - 0.05, max_val + 0.05)

    df["quartile"] = pd.qcut(df["mean_silhouette"], 4, labels=["Q1", "Q2", "Q3", "Q4"])

    fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(20, 16), sharex=False)
    metric_keys = list(metrics.keys())
    quartiles = ["Q1", "Q2", "Q3", "Q4"]

    for row, metric in enumerate(metric_keys):
        for col, q in enumerate(quartiles):
            ax = axes[row][col]
            subset = df[df["quartile"] == q]
            sns.boxplot(y=subset[metric], ax=ax, color="lightgray")
            sns.stripplot(y=subset[metric], ax=ax, color="blue", size=6, jitter=True, alpha=0.7)
            ax.set_title(f"{metrics[metric]} in {q}")
            ax.set_ylabel(metrics[metric] if col == 0 else "")
            ax.set_ylim(ylims[metric])
            for i, val in enumerate(subset[metric]):
                ax.text(0, val, f"{val:.2f}", ha="left", va="center", fontsize=6)

    fig.suptitle(f"Cluster Summary by Quartile (Group by mean silhouette) — {label}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    outfile = f"{prefix}_dataset_summary_quartiles.png"
    fig.savefig(outfile, dpi=300)
    if show:
        plt.show()
    plt.close(fig)
