import anndata as ad
import pandas as pd
import numpy as np
from scsilhouette.compute import compute_silhouette_scores

def test_compute_silhouette():
    obs = pd.DataFrame({
        "cell_type": ["A"] * 10 + ["B"] * 10,
        "author_cell_type": ["X"] * 5 + ["Y"] * 5 + ["Z"] * 10
    })
    X = np.random.rand(20, 5)
    adata = ad.AnnData(X=X, obs=obs)
    adata.obsm["X_pca"] = np.random.rand(20, 3)

    result_df, summary_df = compute_silhouette_scores(
        adata=adata,
        label_fields=["cell_type", "author_cell_type"],
        embedding_key="X_pca"
    )

    assert "silhouette_score" in result_df.columns
    assert summary_df.shape[0] >= 2  # At least 2 clusters

