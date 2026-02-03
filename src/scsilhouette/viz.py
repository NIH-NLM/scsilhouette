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
    sort_by: str = "median",
    output_dir: str = None,
):
    """Generate silhouette summary boxplot with optional F-scores"""
    
    # File prefix pattern: outputs_{organ}_{author}_{year}/{cluster_header}_
    prefix = f"{output_dir}/{cluster_header}"
    logger.info(f"output prefix for files is {prefix}")

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

    # Calculate summary statistics
    median_of_medians = grouped["median"].median()
    median_of_fscores = grouped["f_score"].median() if has_fscore else None

    logger.info("Building plotly figure...")
    fig = go.Figure()

    # Silhouette boxplots
    for i, label in enumerate(ordered_labels):
        y = df[df[cluster_header] == label][silhouette_score_col]
        fig.add_trace(
            go.Box(
                y=y,
                x=[label] * len(y),
                name=label,
                marker_color="lightblue",
                line_color="blue",
                boxmean=True,
                showlegend=(i == 0), # show legend for first item only
                legendgroup="silhouette",
                yaxis="y1",
                offsetgroup="silhouette",
            )
        )

    # Update the first trace to have proper legend name
    if len(fig.data) > 0:
        fig.data[0].name = "Silhouette Score"
        fig.data[0].showlegend = True

    # F-score bars
    if has_fscore:
        fig.add_trace(
            go.Bar(
                x=grouped[cluster_header],
                y=grouped["f_score"],
                marker_color="orange",
                name="F-score",
                yaxis="y2",
                opacity=0.7,
                offsetgroup="fscore",
                legendgroup="fscore",
                showlegend=True,
            )
        )

    # Cluster size annotations - horizontal text just below y=1
    annotations = []
    for _, row in grouped.iterrows():
        annotations.append(
            dict(
                x=row[cluster_header],
                y=0.92,
                text=f"n={int(row['count'])}",
                showarrow=False,
                yanchor="bottom",
                xanchor="center",
                font=dict(size=9),
                textangle=-90,
                yref="paper"
            )
        )

    # Add summary statistics above legend
    summary_text = f"Median of Silhouette medians: {median_of_medians:.3f}"
    if median_of_fscores is not None:
        summary_text += f"<br>Median of F-scores: {median_of_fscores:.3f}"
    
    annotations.append(
        dict(
            x=0.98,
            y=1.03,
            xref="paper",
            yref="paper",
            text=summary_text,
            showarrow=False,
            xanchor="right",
            yanchor="bottom",
            font=dict(size=11, color="darkblue"),
            bgcolor="white",
            bordercolor="black",
            borderwidth=1,
            borderpad=4
        )
    )

    # Layout
    fig.update_layout(
        boxmode="group",
        xaxis=dict(
            title=cluster_header,
            tickangle=-45,
            tickfont=dict(size=10),
            domain=[0.02, 0.98]
        ),
        yaxis=dict(
            title="Silhouette Score",
            range=[-1, 1],
            domain=[0, 0.92]
        ),
        yaxis2=dict(
            title="F-score",
            overlaying="y",
            side="right",
            range=[0, 1],
            showgrid=False,
        ),
        title=dict(
            text="Silhouette Summary with F-score",
            x=0.5,
            xanchor="center",
            font=dict(size=14)
        ),
        width=1800,
        height=900,
        annotations=annotations,
        margin=dict(t=100, b=100, l=80, r=80),
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=0.98,
            xanchor="right",
            x=0.98
        ),
        bargap=0.15,
        bargroupgap=0.05
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
                "height": 900,
                "width": 1800,
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
    logger.info(f"Summary - Median of Medians: {median_of_medians:.3f}")
    if median_of_fscores:
        logger.info(f"Summary - Median of F-scores: {median_of_fscores:.3f}")


def plot_dotplot(
    h5ad_path: str,
    cluster_header: str,
    embedding_key: str,
    organ: str,
    first_author: str,
    year: str,
    output_dir: str,
):
    """Generate embedding dotplot colored by cluster"""
    
    import scanpy as sc
    
    logger.info("Loading data...")
    adata = sc.read_h5ad(h5ad_path)
    
    # Always auto-construct output directory
    prefix = f"{output_dir}/{cluster_header}"
    logger.info(f"output prefix for files is {prefix}")

    os.makedirs(output_dir, exist_ok=True)
    
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
    
    # Output with standard naming pattern
    output_prefix = f"{prefix}_dotplot_{embedding_key}"
    fig.write_html(f"{output_prefix}.html")
    fig.write_image(f"{output_prefix}.svg")
    
    logger.info(f"Saved: {output_prefix}.html, {output_prefix}.svg")

def plot_distribution(
    cluster_summary_path: str,
    cluster_header: str,
    organ: str,
    first_author: str,
    year: str,
    output_dir: str,
):
    """Generate distribution plots of cluster sizes vs silhouette"""
    
    logger.info("Loading cluster summary...")
    df = pd.read_csv(cluster_summary_path)
    
    # Always auto-construct output directory
    prefix = f"{output_dir}/{cluster_header}"
    logger.info(f"output prefix for files is {prefix}")

    os.makedirs(output_dir, exist_ok=True)

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
        width=1800,  # Match viz-summary
        height=900,  # Match viz-summary
        margin=dict(l=80, r=100, t=80, b=150),  # More bottom margin for angled labels
        xaxis=dict(
            title=cluster_header,
            tickangle=-45,  # 45 degree angle like viz-summary
            tickfont=dict(size=10)
        )
    )

    fig_log.update_yaxes(title_text="log10(Cell Count)", secondary_y=False)
    fig_log.update_yaxes(title_text="Silhouette Score", secondary_y=True)

    output_prefix_log = f"{prefix}_distribution_log10"
    try:
        fig_log.write_html(f"{output_prefix_log}.html")
        fig_log.write_image(f"{output_prefix_log}.svg")
        logger.info(f"Saved: {output_prefix_log}.html, {output_prefix_log}.svg")
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
        width=1800,  # Match viz-summary
        height=900,  # Match viz-summary
        margin=dict(l=80, r=100, t=80, b=150),  # More bottom margin for angled labels
        xaxis=dict(
            title=cluster_header,
            tickangle=-45,  # 45 degree angle like viz-summary
            tickfont=dict(size=10)
        )
    )

    fig_raw.update_yaxes(title_text="Cell Count", secondary_y=False)
    fig_raw.update_yaxes(title_text="Silhouette Score", secondary_y=True)

    output_prefix_raw = f"{prefix}_distribution_raw"
    try:
        fig_raw.write_html(f"{output_prefix_raw}.html")
        fig_raw.write_image(f"{output_prefix_raw}.svg")
        logger.info(f"Saved: {output_prefix_raw}.html, {output_prefix_raw}.svg")
    except Exception as e:
        logger.warning(f"Export raw failed: {e}")

