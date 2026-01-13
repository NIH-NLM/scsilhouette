from typing import Optional, List

from typing import Optional
import os
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

# src/scsilhouette/viz.py

import os
from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go

def plot_silhouette_summary(
    silhouette_score_path: str,
    silhouette_score_col: str,
    label_key: str,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
    sort_by: str = "median",
):
    df = pd.read_csv(silhouette_score_path)
    grouped = df.groupby(label_key)[silhouette_score_col].agg(["mean", "std", "median"]).reset_index()
    
    # Cluster sizes from score table
    cluster_sizes = df[label_key].value_counts().rename("n").reset_index()
    cluster_sizes.columns = [label_key, "n"]
    grouped = grouped.merge(cluster_sizes, on=label_key, how="left")

    if fscore_path and mapping_path:
        fscore_df = pd.read_csv(fscore_path)
        mapping_df = pd.read_csv(mapping_path)
        mapping_df.columns = [label_key, "clusterName"]
        fscore_df = fscore_df.merge(mapping_df, on="clusterName", how="inner")
        grouped = grouped.merge(fscore_df[[label_key, "f_score"]], on=label_key, how="left")

    if sort_by not in ["mean", "median", "std"]:
        raise ValueError("sort_by must be one of 'mean', 'median', or 'std'")
    
    grouped = grouped.sort_values(by=sort_by, ascending=False)
    ordered_labels = grouped[label_key].tolist()

    # Static plot
    fig, ax1 = plt.subplots(figsize=(14, 6))

    sns.boxplot(
        data=df,
        x=label_key,
        y=silhouette_score_col,
        order=ordered_labels,
        ax=ax1,
        color="lightblue",
        showmeans=True,
        meanprops={"marker": "o", "markerfacecolor": "blue", "markeredgecolor": "black"},
    )

    if "f_score" in grouped.columns:
        ax2 = ax1.twinx()
        ax2.bar(
            x=np.arange(len(grouped)) + 0.3,
            height=grouped["f_score"],
            width=0.4,
            color="orange",
            alpha=0.7,
            label="F-score"
        )
        for i, f in enumerate(grouped["f_score"]):
            ax2.text(i + 0.3, f + 0.02, f"{f:.2f}", ha="center", fontsize=6)

    for i, n in enumerate(grouped["n"]):
        ax1.text(i, ax1.get_ylim()[0] - 0.1, f"n={n}", ha="center", va="top", fontsize=6, rotation=90)

    ax1.set_title("Silhouette Score Summary with F-score")
    ax1.set_ylabel("Silhouette Score")
    ax1.set_xlabel(label_key)
    ax1.set_xticks(np.arange(len(ordered_labels)))
    ax1.set_xticklabels(ordered_labels, rotation=90)

    fig.tight_layout()
    out_prefix = Path(silhouette_score_path).stem
    fig.savefig(f"{out_prefix}_silhouette_fscore_summary.svg", bbox_inches="tight")
    plt.close(fig)

    # Interactive HTML
    fig_html = go.Figure()
    for label in ordered_labels:
        y_vals = df[df[label_key] == label][silhouette_score_col]
        fig_html.add_trace(go.Box(y=y_vals, name=label, boxmean=True, marker_color="lightblue"))

    if "f_score" in grouped.columns:
        fig_html.add_trace(go.Bar(
            x=ordered_labels,
            y=grouped["f_score"],
            name="F-score",
            marker_color="orange",
            opacity=0.6,
            yaxis="y2"
        ))

    fig_html.update_layout(
        title="Silhouette Score Summary with F-score",
        xaxis=dict(title=label_key),
        yaxis=dict(title="Silhouette Score"),
        yaxis2=dict(title="F-score", overlaying='y', side='right'),
        boxmode="group",
        margin=dict(l=40, r=40, t=40, b=150)
    )
    fig_html.write_html(f"{out_prefix}_silhouette_fscore_summary.html")

    # Export CSV summary
    grouped.to_csv(f"{out_prefix}_silhouette_fscore_summary.csv", index=False)
    # Save outputs
#    svg_path = f"{prefix}_silhouette_fscore_summary.svg"
    html_path = f"{prefix}_silhouette_fscore_summary.html"
    csv_path = f"{prefix}_silhouette_fscore_summary.csv"

    # Save SVG
#    fig.write_image(svg_path)

    # Save HTML
    pio.write_html(
        fig,
        html_path,
        include_plotlyjs="cdn",
        config={
            "toImageButtonOptions": {
                "format": "svg",
                "filename": prefix + "_silhouette_fscore_summary",
                "height": 800,
                "width": 1600,
                "scale": 1,
            },
            "displaylogo": False,
        }
    )

    grouped.to_csv(csv_path, index=False)

    print(f"Saved SVG: {svg_path}")
    print(f"Saved HTML: {html_path}")
    print(f"Saved CSV: {csv_path}")


def plot_correlation_summary(
    cluster_summary_path: str,
    x_metric: str,
    y_metrics: List[str],
    label: str,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
):
    import pandas as pd
    import plotly.express as px
    from plotly.io import write_image, write_html
    from pathlib import Path
    from scipy.stats import pearsonr

    df = pd.read_csv(cluster_summary_path)
    prefix = Path(cluster_summary_path).stem

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
        title = f"{y_metric} vs {x_metric} (r={r:.2f}, p={p:.3g})"

        fig = px.scatter(
            df,
            x=x_metric,
            y=y_metric,
            hover_name=label,
            title=title,
            width=1200,
            height=600,
        )
        fig.add_shape(
            type="line",
            x0=x.min(), x1=x.max(),
            y0=y.min(), y1=y.max(),
            line=dict(dash="dot", color="gray")
        )

        output_prefix = f"{prefix}_{x_metric}_vs_{y_metric}"
        try:
            fig.write_html(f"{output_prefix}.html")
#            fig.write_image(f"{output_prefix}.svg")
#            fig.write_image(f"{output_prefix}.png")
        except Exception as e:
            print(f"[WARN] Export failed: {e}")

        results.append({
            "x": x_metric,
            "y": y_metric,
            "pearson_r": r,
            "p_value": p,
            "n": len(df)
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv(f"{prefix}_correlation_summary.csv", index=False)

def plot_dotplot(
    h5ad_path: str,
    label_key: str,
    embedding_key: str,
):
    adata = sc.read_h5ad(h5ad_path)
    prefix = Path(h5ad_path).stem
    suffix = prefix
    
    sc.pl.embedding(
        adata,
        basis=embedding_key.replace("X_", ""),
        color=label_key,
        show=False,
        save=f"dotplot_{suffix}_{embedding_key}.png",
    )

def plot_distribution(
    cluster_summary_path: str,
    label_key: str,
):

    df = pd.read_csv(cluster_summary_path)
    prefix = Path(cluster_summary_path).stem

    df['count_log10'] = np.log10(df['count'].replace(0, np.nan))

    # === Plot 1: log10 count + silhouette
    df_log = df.sort_values("count_log10", ascending=False).reset_index(drop=True)
    x = df_log[label_key]

    fig_log = make_subplots(specs=[[{"secondary_y": True}]])

    fig_log.add_trace(go.Bar(
        x=x, y=df_log["count_log10"],
        name="log10(Cell Count)",
        marker_color="steelblue"
    ), secondary_y=False)

    fig_log.add_trace(go.Scatter(
        x=x, y=df_log["mean_silhouette"],
        name="Mean", mode="lines+markers", line=dict(color="blue")
    ), secondary_y=True)

    fig_log.add_trace(go.Scatter(
        x=x, y=df_log["median_silhouette"],
        name="Median", mode="lines+markers", line=dict(color="red")
    ), secondary_y=True)

    fig_log.update_layout(
        title="Cluster Cell Count (log10) vs Silhouette",
        width=1200, height=600,
        margin=dict(l=60, r=60, t=80, b=60),
        xaxis_tickangle=90
    )

    fig_log.update_yaxes(title_text="log10(Cell Count)", secondary_y=False)
    fig_log.update_yaxes(title_text="Silhouette Score", secondary_y=True)

    try:
        fig_log.write_html(f"{prefix}_log10.html")
#        fig_log.write_image(f"{prefix}_log10.svg")
#        fig_log.write_image(f"{prefix}_log10.png")
    except Exception as e:
        print(f"[WARN] Export log10 failed: {e}")

    # === Plot 2: raw count + silhouette
    df_raw = df.sort_values("count", ascending=False).reset_index(drop=True)
    x = df_raw[label_key]

    fig_raw = make_subplots(specs=[[{"secondary_y": True}]])

    fig_raw.add_trace(go.Bar(
        x=x, y=df_raw["count"],
        name="Cell Count",
        marker_color="gray"
    ), secondary_y=False)

    fig_raw.add_trace(go.Scatter(
        x=x, y=df_raw["mean_silhouette"],
        name="Mean", mode="lines+markers", line=dict(color="blue")
    ), secondary_y=True)

    fig_raw.add_trace(go.Scatter(
        x=x, y=df_raw["median_silhouette"],
        name="Median", mode="lines+markers", line=dict(color="red")
    ), secondary_y=True)

    fig_raw.update_layout(
        title="Cluster Cell Count (raw) vs Silhouette",
        width=1200, height=600,
        margin=dict(l=60, r=60, t=80, b=60),
        xaxis_tickangle=90
    )

    fig_raw.update_yaxes(title_text="Cell Count", secondary_y=False)
    fig_raw.update_yaxes(title_text="Silhouette Score", secondary_y=True)

    try:
        fig_raw.write_html(f"{prefix}_raw.html")
#        fig_raw.write_image(f"{prefix}_raw.svg")
#        fig_raw.write_image(f"{prefix}_raw.png")
    except Exception as e:
        print(f"[WARN] Export raw failed: {e}")

