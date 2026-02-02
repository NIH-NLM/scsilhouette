# src/scsilhouette/viz.py

import os
from pathlib import Path
from typing import Optional, List
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.stats import pearsonr
from .logging_config import setup_logger

logger = setup_logger()


def plot_silhouette_summary(
    silhouette_score_path: str,
    silhouette_score_col: str,
    cluster_header: str,
    fscore_path: str = None,
    output_prefix: str = None,
    sort_by: str = "median",
):
    """Generate silhouette summary boxplot with optional F-scores"""
    
    logger.info("Loading silhouette scores...")
    df = pd.read_csv(silhouette_score_path)

    if silhouette_score_col not in df.columns or cluster_header not in df.columns:
        raise ValueError(f"Missing required columns in silhouette data")

    logger.info("Computing per-cluster statistics...")
    grouped = (
        df.groupby(cluster_header)[silhouette_score_col]
        .agg(["mean", "std", "median", "count"])
        .reset_index()
    )

    # Add IQR
    q1 = df.groupby(cluster_header)[silhouette_score_col].quantile(0.25)
    q3 = df.groupby(cluster_header)[silhouette_score_col].quantile(0.75)
    grouped["q1"] = q1.values
    grouped["q3"] = q3.values

    # F-score integration
    has_fscore = False
    if fscore_path is not None:
        logger.info(f"Loading F-scores from {fscore_path}...")
        fscore_df = pd.read_csv(fscore_path)
        if "clusterName" in fscore_df.columns:
            fscore_df = fscore_df.rename(columns={"clusterName": cluster_header})

        if cluster_header in fscore_df.columns and "f_score" in fscore_df.columns:
            grouped = grouped.merge(
                fscore_df[[cluster_header, "f_score"]],
                on=cluster_header,
                how="left",
            )
            has_fscore = True

    # Sort
    grouped = grouped.sort_values(sort_by, ascending=False)
    ordered_labels = grouped[cluster_header].tolist()

    logger.info("Building plotly figure...")
    fig = go.Figure()

    # Silhouette boxplots
    for label in ordered_labels:
        y = df[df[cluster_header] == label][silhouette_score_col]
        fig.add_trace(
            go.Box(
                y=y,
                x=[label] * len(y),
                name=label,
                marker_color="lightblue",
                line_color="blue",
                boxmean=True,
                showlegend=False,
                yaxis="y1",
            )
        )

    # F-score bars
    if has_fscore:
        fig.add_trace(
            go.Bar(
                x=grouped[cluster_header],
                y=grouped["f_score"],
                marker_color="orange",
                name="F-score",
                yaxis="y2",
                opacity=0.8,
            )
        )

    # Cluster size annotations
    annotations = []
    for _, row in grouped.iterrows():
        annotations.append(
            dict(
                x=row[cluster_header],
                y=row["q3"] + 0.05,
                text=f"n={int(row['count'])}",
                showarrow=False,
                yanchor="bottom",
                font=dict(size=10),
            )
        )

    prefix = (
        output_prefix
        if output_prefix is not None
        else os.path.splitext(os.path.basename(silhouette_score_path))[0]
    )

    # Layout
    fig.update_layout(
        boxmode="group",
        xaxis=dict(title=cluster_header, tickangle=45),
        yaxis=dict(title="Silhouette Score", range=[-1, 1]),
        yaxis2=dict(
            title="F-score",
            overlaying="y",
            side="right",
            range=[0, 1],
            showgrid=False,
        ),
        title=f"Silhouette Summary with F-score",
        width=1600,
        height=800,
        annotations=annotations,
    )

    # Output paths
    svg_path = f"{prefix}_silhouette_fscore_summary.svg"
    html_path = f"{prefix}_silhouette_fscore_summary.html"
    csv_path = f"{prefix}_silhouette_fscore_summary.csv"

    # Write HTML
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
        },
    )

    # Save SVG
    fig.write_image(svg_path)

    # Save summary table
    grouped.to_csv(csv_path, index=False)

    logger.info(f"Saved HTML: {html_path}")
    logger.info(f"Saved SVG: {svg_path}")
    logger.info(f"Saved CSV: {csv_path}")


def plot_correlation_summary(
    cluster_summary_path: str,
    x_metric: str,
    y_metrics: List[str],
    cluster_header: str,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
):
    """Generate correlation scatter plots between metrics"""
    
    logger.info("Loading cluster summary...")
    df = pd.read_csv(cluster_summary_path)
    prefix = Path(cluster_summary_path).stem

    if fscore_path:
        logger.info(f"Loading F-scores from {fscore_path}...")
        fscore_df = pd.read_csv(fscore_path)
        if mapping_path:
            mapping = pd.read_csv(mapping_path)
            mapping.columns = [cluster_header, "clusterName"]
            fscore_df = fscore_df.merge(mapping, on="clusterName", how="inner")
        df = df.merge(fscore_df[[cluster_header, "f_score"]], on=cluster_header, how="left")

    results = []

    for y_metric in y_metrics:
        if y_metric not in df.columns or x_metric not in df.columns:
            logger.warning(f"Skipping: {x_metric} or {y_metric} not in dataframe columns")
            continue

        x = df[x_metric]
        y = df[y_metric]
        r, p = pearsonr(x, y)
        title = f"{y_metric} vs {x_metric} (r={r:.2f}, p={p:.3g})"

        logger.info(f"Generating correlation plot: {y_metric} vs {x_metric}")
        fig = px.scatter(
            df,
            x=x_metric,
            y=y_metric,
            hover_name=cluster_header,
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
            fig.write_image(f"{output_prefix}.svg")
            logger.info(f"Saved: {output_prefix}.html, {output_prefix}.svg")
        except Exception as e:
            logger.warning(f"Export failed: {e}")

        results.append({
            "x": x_metric,
            "y": y_metric,
            "pearson_r": r,
            "p_value": p,
            "n": len(df)
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv(f"{prefix}_correlation_summary.csv", index=False)
    logger.info(f"Saved correlation summary: {prefix}_correlation_summary.csv")


def plot_dotplot(
    h5ad_path: str,
    cluster_header: str,
    embedding_key: str,
):
    """Generate embedding dotplot colored by cluster"""
    
    import scanpy as sc
    
    logger.info("Loading data...")
    adata = sc.read_h5ad(h5ad_path)
    prefix = Path(h5ad_path).stem
    
    # Get embedding coordinates
    embedding = adata.obsm[embedding_key]
    df_plot = pd.DataFrame({
        'x': embedding[:, 0],
        'y': embedding[:, 1],
        cluster_header: adata.obs[cluster_header]
    })
    
    logger.info("Creating plotly scatter plot...")
    fig = px.scatter(
        df_plot,
        x='x',
        y='y',
        color=cluster_header,
        title=f"{embedding_key} colored by {cluster_header}",
        width=1200,
        height=800
    )
    
    fig.update_layout(
        xaxis_title=f"{embedding_key}_1",
        yaxis_title=f"{embedding_key}_2"
    )
    
    output_prefix = f"dotplot_{prefix}_{embedding_key}"
    fig.write_html(f"{output_prefix}.html")
    fig.write_image(f"{output_prefix}.svg")
    
    logger.info(f"Saved: {output_prefix}.html, {output_prefix}.svg")


def plot_distribution(
    cluster_summary_path: str,
    cluster_header: str,
):
    """Generate distribution plots of cluster sizes vs silhouette"""
    
    logger.info("Loading cluster summary...")
    df = pd.read_csv(cluster_summary_path)
    prefix = Path(cluster_summary_path).stem

    df['count_log10'] = np.log10(df['count'].replace(0, np.nan))

    # Plot 1: log10 count + silhouette
    logger.info("Generating log10 count distribution...")
    df_log = df.sort_values("count_log10", ascending=False).reset_index(drop=True)
    x = df_log[cluster_header]

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
        fig_log.write_image(f"{prefix}_log10.svg")
        logger.info(f"Saved: {prefix}_log10.html, {prefix}_log10.svg")
    except Exception as e:
        logger.warning(f"Export log10 failed: {e}")

    # Plot 2: raw count + silhouette
    logger.info("Generating raw count distribution...")
    df_raw = df.sort_values("count", ascending=False).reset_index(drop=True)
    x = df_raw[cluster_header]

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
        fig_raw.write_image(f"{prefix}_raw.svg")
        logger.info(f"Saved: {prefix}_raw.html, {prefix}_raw.svg")
    except Exception as e:
        logger.warning(f"Export raw failed: {e}")
