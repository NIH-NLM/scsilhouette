import anndata as ad
import numpy as np
import pandas as pd

from scsilhouette import viz


def adata_with(labels, colours=None):
    obs = pd.DataFrame({"type": pd.Categorical(labels)}, index=[f"c{i}" for i in range(len(labels))])
    adata = ad.AnnData(X=np.zeros((len(labels), 1)), obs=obs)
    if colours is not None:
        adata.uns["type_colors"] = colours
    return adata


def test_the_cellxgene_colours_are_used_in_the_order_of_the_categories():
    adata = adata_with(["b", "a", "c", "a"], ["#111111", "#222222", "#333333"])
    assert viz.cluster_colors(adata, "type") == {"a": "#111111", "b": "#222222", "c": "#333333"}


def test_without_cellxgene_colours_every_cluster_gets_its_own_colour():
    labels = [f"cluster{i}" for i in range(60)]
    colours = viz.cluster_colors(adata_with(labels), "type")
    assert set(colours) == set(labels) and len(set(colours.values())) == 60


def test_colours_that_do_not_match_the_clusters_are_not_used():
    adata = adata_with(["a", "b", "c"], ["#111111", "#222222"])
    assert "#111111" not in viz.cluster_colors(adata, "type").values()
