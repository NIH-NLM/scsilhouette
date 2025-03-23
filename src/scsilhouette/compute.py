from typing import Dict, Optional, List
import pandas as pd
import numpy as np
import anndata as ad
from sklearn.metrics import silhouette_samples

def compute_silhouette_score(
    adata: ad.AnnData,
    label_field: str,
    embedding_key: str = "X_pca",
    metadata_fields: Optional[List[str]] = None
) -> Dict:
    if label_field not in adata.obs:
        raise ValueError(f"Label field '{label_field}' not found in `.obs`.")

    if embedding_key not in adata.obsm:
        raise ValueError(f"Embedding key '{embedding_key}' not found in `.obsm`.")

    X = adata.obsm[embedding_key]
    labels = adata.obs[label_field].astype(str)
    sil_scores = silhouette_samples(X, labels)

    cell_scores = pd.DataFrame({
        "cell_id": adata.obs_names,
        "label": labels,
        "silhouette_score": sil_scores
    }).set_index("cell_id")

    metadata_fields = metadata_fields or []
    for field in metadata_fields:
        if field in adata.obs.columns:
            cell_scores[field] = adata.obs[field]

    cell_scores.reset_index(inplace=True)

    cluster_scores = (
        cell_scores.groupby("label")["silhouette_score"]
        .mean()
        .to_dict()
    )

    cluster_summary = (
        cell_scores
        .groupby("label")
        .agg(
            cluster_size=("silhouette_score", "count"),
            mean_silhouette=("silhouette_score", "mean")
        )
        .sort_values("mean_silhouette", ascending=False)
        .reset_index()
    )

    return {
        "cluster_header": label_field,
        "embedding_key": embedding_key,
        "mean_score": float(np.mean(sil_scores)),
        "cluster_scores": cluster_scores,
        "cell_scores": cell_scores,
        "cluster_summary": cluster_summary
    }

