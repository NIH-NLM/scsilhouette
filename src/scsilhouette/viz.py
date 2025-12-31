from typing import Optional, List

def plot_silhouette_summary(
    silhouette_score_path: str,
    silhouette_score_col: str,
    label_key: str,
    fscore_path: Optional[str] = None,
    mapping_path: Optional[str] = None,
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

    df = pd.read_csv(silhouette_score_path)
    prefix = Path(silhouette_score_path).stem
    suffix = prefix
    
    # Output files
    outbase = f"{suffix}_silhouette_summary"
    html_path = f"{outbase}.html"
    svg_path = f"{outbase}.svg"

    grouped = df.groupby(label_key)[silhouette_score_col].agg([
        "mean", "std", "median", "count", "min", "max",
        lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75)
    ]).reset_index()
    grouped.columns = [label_key, "mean", "std", "median", "count", "min", "max", "q1", "q3"]
    
    has_fscore = False
    if fscore_path:
        fscore_df = pd.read_csv(fscore_path)

        if label_key not in fscore_df.columns and 'clusterName' in fscore_df.columns:
            fscore_df = fscore_df.rename(columns={'clusterName': label_key})
            
        grouped = grouped.merge(fscore_df[[label_key, "f_score"]], on=label_key, how="left")
        has_fscore = True

    # --- BEGIN AUTO-MAPPING PATCH ---
    if mapping_path:
        mapping = pd.read_csv(mapping_path)
        mapping.columns = [label_key, "clusterName"]
        fscore_df = fscore_df.merge(mapping, on="clusterName", how="inner")
    else:
        # Auto-map when clusterName exists but label_key does not
        if label_key not in fscore_df.columns:
            if "clusterName" in fscore_df.columns:
                print(
                    f"[INFO] Auto-mapping '{label_key}' ← 'clusterName' "
                    "(no mapping file provided)"
                )
                fscore_df[label_key] = fscore_df["clusterName"]
            else:
                raise ValueError(
                    f"[ERROR] Cannot merge fscore file: "
                    f"'{label_key}' not found and no 'clusterName' column present."
                )
    # --- END AUTO-MAPPING PATCH ---

    grouped = grouped.sort_values(by=sort_by, ascending=False)
    sorted_labels = grouped[label_key].tolist()

    x = np.arange(len(grouped))

    # Static SVG with matplotlib
    fig, ax = plt.subplots(figsize=(max(14, 0.25 * len(grouped)), 6))
    ax.bar(x - 0.2, grouped["median"], width=0.4, label="Median Silhouette", color="steelblue")

    if "f_score" in grouped.columns:
        ax.bar(x + 0.2, grouped["f_score"], width=0.4, label="F-score", color="orange")

    for i, (mean, std, median) in enumerate(zip(grouped["mean"], grouped["std"], grouped["median"])):
        ax.errorbar(i - 0.2, mean, yerr=std, fmt='o', color='black', capsize=4)
        ax.scatter(i - 0.2, median, color='red', zorder=5)
        ax.text(i - 0.25, mean + std + 0.02, f"μ={mean:.2f}", ha="center", fontsize=6)
        ax.text(i - 0.20, mean + std + 0.06, f"M={median:.2f}", ha="center", fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(grouped[label_key], rotation=90, fontsize=6)
    ax.set_ylabel("Silhouette / F-score")
    ax.set_title(f"Silhouette Summary with F-score per {label_key}")
    ax.legend()

    fig.savefig(f"{prefix}_silhouette_fscore_summary.svg", bbox_inches="tight")
    plt.close(fig)

    # Interactive HTML with plotly
    fig_html = go.Figure()
    fig_html.add_trace(go.Bar(
        x=grouped[label_key],
        y=grouped["median"],
        name="Median Silhouette",
        marker_color="steelblue"
    ))

    if "f_score" in grouped.columns:
        fig_html.add_trace(go.Bar(
            x=grouped[label_key],
            y=grouped["f_score"],
            name="F-score",
            marker_color="orange"
        ))

    fig_html.add_trace(go.Scatter(
        x=grouped[label_key],
        y=grouped["mean"],
        mode="markers+lines",
        name="Mean Silhouette",
        error_y=dict(type="data", array=grouped["std"]),
        marker=dict(color="black", size=6)
    ))

    fig_html.update_layout(
        title=f"Silhouette Summary with F-score per {label_key}",
        xaxis_title=label_key,
        yaxis_title="Silhouette / F-score",
        barmode="group"
    )

    plotly_offline(fig_html, filename=f"{prefix}_silhouette_fscore_summary.html", auto_open=False)

    grouped.to_csv(f"{prefix}_silhouette_fscore_summary.csv", index=False)


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

