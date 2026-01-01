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

    # Load silhouette score data
    df = pd.read_csv(silhouette_score_path)
    if silhouette_score_col not in df.columns or label_key not in df.columns:
        raise ValueError("Missing required columns in silhouette data")

    # Compute summary statistics for silhouette scores
    grouped = df.groupby(label_key)[silhouette_score_col].agg(
        ["mean", "std", "median", "count", "min", "max"]
    ).reset_index()
    grouped["q1"] = df.groupby(label_key)[silhouette_score_col].quantile(0.25).values
    grouped["q3"] = df.groupby(label_key)[silhouette_score_col].quantile(0.75).values

    # Move this up here
    has_fscore = False
    if fscore_path is not None:
        fscore_df = pd.read_csv(fscore_path)
        if "clusterName" in fscore_df.columns:
            fscore_df = fscore_df.rename(columns={"clusterName": label_key})
        if label_key in fscore_df.columns and "f_score" in fscore_df.columns:
            grouped = grouped.merge(
                fscore_df[[label_key, "f_score"]], on=label_key, how="left"
            )
            has_fscore = True

    # Now sort AFTER merging
    grouped = grouped.sort_values("median", ascending=False)
    labels = grouped[label_key].tolist()

    # Create box plot traces for silhouette scores
    box_traces = []
    for i, label in enumerate(labels):
        y = df[df[label_key] == label][silhouette_score_col].values
        box_traces.append(
            go.Box(
                y=y,
                x=[label] * len(y),
                name="Silhouette Score",
                legendgroup="silhouette",  # Keeps consistent legend entry
                boxpoints="outliers",
                marker_color="blue",
                line_color="blue",
                offsetgroup=label, # ensures side-by-side alignment
                showlegend=(i == 0),
                yaxis="y1",
            )
        )

    # Create bar chart traces for F-scores
    bar_traces = []
    if has_fscore:
        for i, label in enumerate(labels):
            fval = grouped.loc[grouped[label_key] == label, "f_score"].values
            if len(fval) == 1:
                bar_traces.append(
                    go.Bar(
                        x=[label],
                        y=[fval[0]],
                        legendgroup="fscore",
                        showlegend=(i == 0),
                        marker_color="orange",
                        offsetgroup=f"{label}_bar",  # differentiate from box
                        name="F-score",
                        yaxis="y2",
                        offset=0.3,  # <-- Pushes it to the right within the group
                        width=0.4,  # <-- Makes the bar wide enough to be seen
                    )
                )

    # Cell count annotations
    annotations = []
    for i, row in enumerate(grouped.itertuples()):
        annotations.append(
            dict(
                x=row.__getattribute__(label_key),
                y=row.max + 0.05,
                text=f"n={row.count}",
                showarrow=False,
                yanchor="bottom",
                font=dict(size=10),
            )
        )

    # Summary stats box (upper right)
    med_of_meds = round(grouped["median"].median(), 3)
    med_fscore = round(grouped["f_score"].median(), 3) if has_fscore else None
    text_lines = [f"Median of medians: {med_of_meds}"]
    if has_fscore:
        text_lines.append(f"Median F-score: {med_fscore}")

    annotations.append(
        dict(
            xref="paper",
            yref="paper",
            x=0.95,
            y=0.95,
            text="<br>".join(text_lines),
            showarrow=False,
            align="right",
            font=dict(size=12, color="black"),
            bordercolor="black",
            borderwidth=1,
            bgcolor="white",
        )
    )

    # Build figure
    fig = go.Figure()
    for trace in box_traces:
        fig.add_trace(trace)
    if has_fscore:
        for bar_trace in bar_traces:
            fig.add_trace(bar_trace)

    # Generate file prefix
    prefix = os.path.splitext(os.path.basename(silhouette_score_path))[0]

    fig.update_layout(
        barmode="group", # <-- This is what tells Plotly to offset them side-by-side
        xaxis=dict(
            title=label_key,
            tickangle=45),
        yaxis=dict(
            title=f"Silhouette Summary with F-scores",
            side="left",
            range=[-1, 1],
        ),
        yaxis2=dict(
            title="F-score",
            overlaying="y",
            side="right",
            range=[0,1],
            showgrid=False,
        ),
        title=f"Silhouette Summary with F-scores – {prefix.replace('_silhouette_scores', '')}",
        width=1600,
        height=800,
        annotations=annotations,
    )


    # Save outputs
    svg_path = f"{prefix}_silhouette_fscore_summary.svg"
    html_path = f"{prefix}_silhouette_fscore_summary.html"
    csv_path = f"{prefix}_silhouette_fscore_summary.csv"

#    fig.write_image(svg_path)

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

