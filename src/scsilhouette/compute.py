import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics import silhouette_samples
from pathlib import Path


def run_silhouette(
    h5ad_path: str,
    label_keys: list,
    embedding_key: str,
    output_dir: str,
    show_obs: bool = False,
    save_scores: bool = False,
    save_cluster_summary: bool = False,
    compute_qc: bool = False,
):
    import warnings
    warnings.simplefilter("ignore")

    adata = sc.read_h5ad(h5ad_path)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if show_obs:
        adata.obs.head().to_csv(Path(output_dir) / "obs_preview.csv")
        print(f"[✓] Saved obs preview: {Path(output_dir) / 'obs_preview.csv'}")

    X = adata.obsm[embedding_key]
    results = {}

    for label_key in label_keys:
        if label_key not in adata.obs:
            print(f"[!] Label '{label_key}' not found in obs. Skipping.")
            continue

        labels = adata.obs[label_key]
        cell_scores = pd.DataFrame(index=adata.obs_names)
        cell_scores[label_key] = labels
        cell_scores["silhouette_score"] = silhouette_samples(X, labels)

        cluster_summary = (
            cell_scores.groupby(label_key)["silhouette_score"]
            .agg(["mean", "std", "median", "min", "max", "count"])
            .rename(columns={
                "mean": "mean_silhouette_score",
                "std": "std_silhouette_score",
                "median": "median_silhouette_score",
                "min": "min_silhouette_score",
                "max": "max_silhouette_score",
                "count": "n_cells"
            })
            .reset_index()
        )

        if save_scores:
            cell_scores.to_csv(Path(output_dir) / f"{label_key}_scores.csv", index=False)
            print(f"[✓] Saved: {label_key}_scores.csv")

        if save_cluster_summary:
            cluster_summary.to_csv(Path(output_dir) / f"{label_key}_cluster_summary.csv", index=False)
            print(f"[✓] Saved: {label_key}_cluster_summary.csv")

        results[label_key] = {
            "cell_scores": cell_scores,
            "cluster_summary": cluster_summary,
        }

    return results


def compute_fscore_correlation(
    silhouette_df: pd.DataFrame,
    nsforest_df: pd.DataFrame,
    label_col: str,
    output_dir: str,
    score_col: str = "mean_silhouette_score",
    fscore_col: str = "F-score",
    distance_metric: str = "euclidean"
):
    """
    Computes and saves correlation between silhouette score and NS-Forest F-score.
    Saves a CSV summary: cluster, silhouette score, F-score, distance metric, correlation.
    """

    from pathlib import Path
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Align on cluster name
    merged = silhouette_df[[label_col, score_col]].merge(
        nsforest_df[["cluster", fscore_col]],
        left_on=label_col, right_on="cluster"
    )

    corr = merged[score_col].corr(merged[fscore_col])

    # Save annotated results
    merged["distance_metric"] = distance_metric
    merged["correlation"] = corr

    out_csv = Path(output_dir) / f"fscore_correlation_{score_col}_{distance_metric}.csv"
    merged.to_csv(out_csv, index=False)
    print(f"[✓] Correlation: {corr:.3f} saved to {out_csv.name}")

    return merged, corr
