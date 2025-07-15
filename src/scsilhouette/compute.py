# src/scsilhouette/compute.py

import os
from pathlib import Path
import pandas as pd
import scanpy as sc
from .utils import map_gene_symbols_to_ensembl

def run_silhouette(
    h5ad_path: Path,
    label_key: str,
    embedding_key: str = None,
    organism: str = None,
    disease: str = None,
    tissue: str = None,
    cell_count: int = None,
    use_binary_genes: bool = False,
    gene_list_path: Path = None,
    mapping_output_path: Path = None,
    metric: str = "euclidean",
    pca_components: int = None,
    filter_normal: bool = True,
    save_scores: bool = False,
    save_cluster_summary: bool = False,
    save_annotation: bool = False,
    show_annotation: bool = False,
):
    from sklearn.metrics import silhouette_samples
    from sklearn.decomposition import PCA
    import numpy as np

    adata = sc.read(h5ad_path)
    prefix = Path(h5ad_path).stem

    print(f" will save output with {prefix}")
    print(f"[PASS] Loaded AnnData with {adata.shape}")

    # Handle gene subsetting if binary genes used
    if use_binary_genes and gene_list_path is not None:
        genes_df = pd.read_csv(gene_list_path)
        genes = genes_df.iloc[:, 0].tolist()
        print(f"[PASS] Loaded {len(genes)} binary genes from {gene_list_path}")

        if adata.raw is not None and all(g in adata.raw.var_names for g in genes):
            adata_use = adata.raw[:, genes]
            print("[PASS] Subsetting using raw and binary genes")
        else:
            adata_use = adata[:, genes]
            print("[WARN] Subsetting using var and binary genes")
    else:
        if embedding_key is None:
            adata_use = adata
        else:
            adata_use = adata.obsm[embedding_key]

    # PCA optionally
    if pca_components is not None:
        print(f"[INFO] Running PCA with {pca_components} components")
        pca = PCA(n_components=pca_components)
        adata_use = pca.fit_transform(adata_use)

    # Compute silhouette scores
#    labels = adata.obs[label_key]
#    scores = silhouette_samples(adata_use, labels, metric=metric)
#    adata.obs["silhouette_score"] = scores
#    labels = adata.obs[label_key]
    if filter_normal and "disease" in adata.obs.columns:
        normal_mask = adata.obs["disease"] == "normal"
        adata_use = adata_use[normal_mask]
        labels = labels[normal_mask]
    scores = silhouette_samples(adata_use, labels, metric=metric)
    adata.obs["silhouette_score"] = scores

    if save_scores:
        scores_path = f"silhouette_scores_{organism}_{disease}_{tissue}_{prefix}_{label_key}_{embedding_key}_{metric}.csv"
        adata.obs[[label_key, "silhouette_score"]].to_csv(scores_path)
        print(f"[PASS] Saved silhouette scores to {scores_path}")

    if save_cluster_summary:
        summary_df = (
            adata.obs.groupby(label_key)
            .agg(
                mean_silhouette=("silhouette_score", "mean"),
                std_silhouette=("silhouette_score", "std"),
                median_silhouette=("silhouette_score", "median"),
                count=("silhouette_score", "count"),
            )
            .reset_index()
        )
        summary_path = f"cluster_summary_{organism}_{disease}_{tissue}_{prefix}_{label_key}_{embedding_key}_{metric}.csv"
        summary_df.to_csv(summary_path)
        print(f"[PASS] Saved cluster summary to {summary_path}")

    if save_annotation:
        annotation_path = f"annotation_{organism}_{disease}_{tissue}_{prefix}_{label_key}_{embedding_key}_{metric}.csv"
        with open(annotation_path, "w") as f:
            f.write(f"\nLoaded file: {h5ad_path}\n")
            f.write(f"Shape: {adata.shape} (cells x genes)\n")
            f.write("Cell Annotations (.obs):\n")
            f.write(", ".join(adata.obs.columns.tolist()) + "\n")
            f.write("Gene Annotations (.var):\n")
            f.write(", ".join(adata.var.columns.tolist()) + "\n")
            f.write("Unstructured Metadata (.uns):\n")
            f.write(", ".join(list(adata.uns.keys())) + "\n")
            f.write("\nExample values from .obs:")
            for col in adata.obs.columns[:5]:  # limit output
                f.write(f"{col}: {adata.obs[col].unique()[:5]}\n")
       
        print(f"[PASS] Saved annotation to {annotation_path}")

    return adata

