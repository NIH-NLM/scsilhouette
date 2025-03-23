# tests/test_compute.py

import anndata as ad
from scsilhouette.compute import compute_silhouette_score

def test_silhouette_computation_on_real_data():
    adata = ad.read_h5ad("data/test.h5ad")  # Replace with real path
    result = compute_silhouette_score(
        adata,
        label_field="author_cell_type",
        embedding_key="X_pca"
    )

    assert "cluster_header" in result
    assert "mean_score" in result
    assert "cluster_scores" in result
    assert "cell_scores" in result
    assert result["mean_score"] > 0  # sanity check
    print("Mean silhouette score:", result["mean_score"])
