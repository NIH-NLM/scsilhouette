import pytest
from scsilhouette import compute
import anndata as ad
import pandas as pd

@pytest.fixture
def sample_adata(tmp_path):
    # Replace with actual loading code or mock
    return ad.read_h5ad("tests/data/sample.h5ad")

def test_run_silhouette(sample_adata):
    results = compute.run_silhouette(
        h5ad_path="tests/data/sample.h5ad",
        label_keys=["cell_type"],
        embedding_key="X_pca",
        output_dir="tests/results",
        show_obs=False,
        qc_correlations=False,
        save_scores=False,
        save_cluster_summary=False,
        save_csv=False
    )
    assert "cell_type" in results
    assert isinstance(results["cell_type"]["cell_scores"], pd.DataFrame)

def test_compute_fscore_correlation(tmp_path):
    fscore_df = pd.DataFrame({
        "cluster": ["A", "B"],
        "fscore": [0.9, 0.8]
    })
    silhouette_df = pd.DataFrame({
        "cluster": ["A", "B"],
        "mean_silhouette_score": [0.7, 0.6]
    })
    merged_df, corr = compute.compute_fscore_correlation(silhouette_df, fscore_df)
    assert "fscore" in merged_df
    assert "mean_silhouette_score" in merged_df
    assert isinstance(corr, float)
