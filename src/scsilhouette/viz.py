from typing import Optional, List

def plot_silhouette_summary(
    silhouette_score_path: str,
    silhouette_score_col: str,
    label_key: str,
    fscore_path: str = None,
    output_prefix: str = None,
    sort_by: str = "median",
):
    import os
    from pathlib import Path
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    from plotly.io import write_html, write_image
    import plotly.io as pio
    from scipy.stats import pearsonr
    import matplotlib.pyplot as plt
    from plotly.offline import plot as plotly_offline

    # Load silhouette scores
    df = pd.read_csv(silhouette_score_path)
    
    # Calculate summary statistics per label
    grouped = df.groupby(label_key)[silhouette_score_col].agg(['mean', 'std', 'median', 'count', 'min', 'max'])
    grouped['q1'] = df.groupby(label_key)[silhouette_score_col].quantile(0.25)
    grouped['q3'] = df.groupby(label_key)[silhouette_score_col].quantile(0.75)
    grouped = grouped.reset_index()

    # Merge in F-score data if provided
    has_fscore = False
    if fscore_path:
        fscore_df = pd.read_csv(fscore_path)
        if 'clusterName' in fscore_df.columns and label_key in df.columns:
            fscore_df = fscore_df.rename(columns={'clusterName': label_key})
        if label_key in fscore_df.columns:
            grouped = grouped.merge(fscore_df[[label_key, 'f_score']], on=label_key, how='left')
            has_fscore = True

    # Order by silhouette median descending
    grouped = grouped.sort_values(by="median", ascending=False)

    # Create the plot
    fig = go.Figure()

    # Silhouette box (whiskers)
    fig.add_trace(go.Box(
        y=grouped['median'],
        x=grouped[label_key],
        name="Silhouette Score (Median)",
        boxpoints=False,
        line=dict(color="blue"),
        marker_color="blue",
        q1=grouped["q1"],
        median=grouped["median"],
        q3=grouped["q3"],
        customdata=grouped[label_key],
        hovertemplate='%{customdata}<br>Median Silhouette: %{y:.3f}<extra></extra>',
    ))

    # F-score box (if available)
    if has_fscore:
        fig.add_trace(go.Box(
            y=grouped["f_score"],
            x=grouped[label_key],
            name="F-score",
            boxpoints=False,
            marker_color="orange",
            line=dict(color="orange"),
            hovertemplate='%{x}<br>F-score: %{y:.3f}<extra></extra>',
        ))

    # Median of medians
    silhouette_median_of_medians = grouped["median"].median()
    fscore_median_of_medians = grouped["f_score"].median() if has_fscore else None

    summary_text = f"<b>Median of Silhouette Medians: {silhouette_median_of_medians:.2f}"
    if fscore_median_of_medians is not None:
        summary_text += f"<br>Median F-score: {fscore_median_of_medians:.2f}"
    summary_text += "</b>"

    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.95, y=0.95,
        text=summary_text,
        showarrow=False,
        align="right",
        bordercolor="black",
        borderwidth=1
    )

    fig.update_layout(
        title="Silhouette Scores and F-scores by Author Cell Type",
        xaxis_title="Cell Type",
        yaxis_title="Score",
        boxmode='group',
        showlegend=True,
        width=1800,
        height=800
    )

    # Output filenames
    prefix = output_prefix or os.path.splitext(silhouette_score_path)[0]
    svg_path = f"{prefix}_silhouette_fscore_plot.svg"
    html_path = f"{prefix}_silhouette_fscore_plot.html"
    stats_path = f"{prefix}_summary_stats.csv"

    try:
        pio.write_image(fig, svg_path)
    except Exception as e:
        print("SVG export failed (kaleido required):", e)

    pio.write_html(fig, html_path)

    # Save summary statistics
    summary_out = pd.DataFrame({
        'median_of_silhouette_medians': [silhouette_median_of_medians],
        'median_of_f_scores': [fscore_median_of_medians if fscore_median_of_medians is not None else '']
    })
    summary_out.to_csv(stats_path, index=False)

    print(f"Saved: {html_path}")
    print(f"Stats: {stats_path}")

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

