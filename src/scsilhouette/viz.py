# src/scsilhouette/viz.py (Plotly: auto-export to HTML, SVG, PNG for Nextflow)
import os
from pathlib import Path
import scanpy as sc
from typing import Optional, List
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from plotly.io import write_html, write_image
import plotly.io as pio
from scipy.stats import pearsonr

def plot_silhouette_summary(
    silhouette_score_path: str,
    silhouette_score_col: str,
    label_key: str,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
    sort_by: str = "median",
):
    df = pd.read_csv(silhouette_score_path)
    prefix = Path(silhouette_score_path).stem
    suffix = prefix
    
    # Output files
    outbase = f"{suffix}_silhouette_summary"
    html_path = f"{outbase}.html"
#    svg_path = f"{outbase}.svg"
#    png_path = f"{outbase}.png"


    grouped = df.groupby(label_key)[silhouette_score_col].agg([
        "mean", "std", "median", "count", "min", "max",
        lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75)
    ]).reset_index()
    grouped.columns = [label_key, "mean", "std", "median", "count", "min", "max", "q1", "q3"]

    has_fscore = False
    if fscore_path and mapping_path:
        fscore_df = pd.read_csv(fscore_path)
        mapping = pd.read_csv(mapping_path)
        mapping.columns = [label_key, "clusterName"]
        fscore_df = fscore_df.merge(mapping, on="clusterName", how="inner")
        grouped = grouped.merge(fscore_df[[label_key, "f_score"]], on=label_key, how="inner")
        has_fscore = True

    grouped = grouped.sort_values(by=sort_by, ascending=False)
    sorted_labels = grouped[label_key].tolist()

    # Create base figure
    fig = go.Figure()

    # Add box plots (uniform color)
    for cat in sorted_labels:
        cat_data = df[df[label_key] == cat][silhouette_score_col]
        fig.add_trace(go.Box(
            y=cat_data,
            name=cat,
            boxpoints=False,
            whiskerwidth=0.8,
            marker=dict(size=3, color='grey'),
            line=dict(width=1, color='black'),
            fillcolor='lightgrey',
            showlegend=False
        ))

    # Add red dots for mean
#    fig.add_trace(go.Scatter(
#        x=sorted_labels,
#        y=grouped["mean"],
#        mode="markers",
#        name="Mean",
#        marker=dict(color="red", size=8, symbol="circle")
#    ))

    # Update layout
    fig.update_layout(
        yaxis=dict(title="Silhouette Score", range=[-1, 1]),
        xaxis=dict(title=label_key, tickangle=45),
        title=f"{outbase}",
        legend=dict(x=1.01, y=1),
        margin=dict(l=60, r=60, t=80, b=60),
        height=600,
        width=1200
    )

    # Add stats annotation
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

    fig.write_html(html_path)
#    pio.write_image(fig, svg_path, format='svg')
#    pio.write_image(fig, png_path, format='png')

    print(f"Saved interactive HTML to {html_path}")
#    print(f"Saved static SVG to {svg_path}")
#    print(f"Saved PNG image to {png_path}")
    

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

