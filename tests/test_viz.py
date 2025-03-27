import pandas as pd
import matplotlib.pyplot as plt
from scsilhouette import viz

def test_plot_score_distribution(tmp_path):
    df = pd.DataFrame({"silhouette_score": [0.1, 0.2, 0.3], "cell_type": ["A", "B", "A"]})
    viz.plot_score_distribution(df, tmp_path, "cell_type")
    assert (tmp_path / "cell_type_score_distribution.png").exists()

def test_plot_cluster_summary(tmp_path):
    df = pd.DataFrame({
        "cell_type": ["A", "B"],
        "mean_silhouette_score": [0.5, 0.6]
    })
    viz.plot_cluster_summary(df, tmp_path, "cell_type")
    assert (tmp_path / "cell_type_cluster_summary.png").exists()

def test_plot_qc_boxplots(tmp_path):
    df = pd.DataFrame({
        "silhouette_score": [0.1, 0.2],
        "cell_type": ["X", "Y"],
        "nCount_RNA": [1000, 2000],
        "nFeature_RNA": [300, 500]
    })
    viz.plot_qc_boxplots(df, tmp_path, "cell_type")
    assert (tmp_path / "cell_type_qc_boxplot.png").exists()
