import pandas as pd
from scsilhouette.viz import (
    plot_score_distribution,
    plot_cluster_summary,
    plot_cluster_size_vs_score
)

def test_viz(tmp_path):
    df = pd.DataFrame({
        "label": ["A", "B", "C"],
        "mean_silhouette": [0.5, 0.6, 0.7],
        "cluster_size": [10, 15, 20]
    })
    scores_df = pd.DataFrame({
        "silhouette_score": [0.2, 0.4, 0.6, 0.8]
    })

    plot_score_distribution(scores_df, output_dir=tmp_path)
    plot_cluster_summary(df, output_dir=tmp_path)
    plot_cluster_size_vs_score(df, output_dir=tmp_path)

    assert (tmp_path / "silhouette_distribution.png").exists()
    assert (tmp_path / "cluster_summary.png").exists()
    assert (tmp_path / "size_vs_score.png").exists()

