Silhouette Score Computation Process
====================================

This document outlines the steps scsilhouette takes to compute silhouette scores for annotated single-cell datasets.

Steps
-----

1. **Input Data**

   An `.h5ad` file is provided containing:

   - `adata.X`: matrix of shape (cells, features)
   - `adata.obs`: with cluster labels (e.g. `cell_type`)
   - `adata.obsm`: with low-dimensional embeddings (e.g. `X_pca`)

2. **Parameters Specified via CLI**

   Example command:

   .. code-block:: bash

      scsilhouette compute \
        --h5ad-path data/my_dataset.h5ad \
        --label-keys cell_type \
        --embedding-key X_pca \
        --output-dir results/ \
        --show-obs \
        --qc-correlations

3. **Silhouette Score Computation**

   For each label key (e.g., `cell_type`, `author_cell_type`), silhouette scores are computed using:

   .. code-block:: python

      from sklearn.metrics import silhouette_samples

      scores = silhouette_samples(embedding, labels)

4. **Cluster Summary Generation**

   After computing cell-level silhouette scores, the following metrics are aggregated per cluster:

   - Mean silhouette score
   - Median
   - Standard deviation
   - Count of cells in cluster

5. **Visualization & QC**

   The following plots are generated:

   - Silhouette score distribution histogram
   - Cluster summary bar plot
   - Correlation with `nCount_RNA` and `nFeature_RNA` (if available)
   - Optional F-Score comparison (from NS-Forest output)

6. **Optional NS-Forest Correlation**

   F-scores can be loaded from `.csv` or `.xlsx` files using:

   .. code-block:: bash

      scsilhouette associate --fscore-file markers.xlsx --output-dir results/

   The tool will then compute the correlation between F-score and silhouette score at the cluster level.

7. **Output Files**

   - `*_scores.csv`: silhouette scores for all cells
   - `*_cluster_summary.csv`: per-cluster summary
   - `obs_preview.csv`: snapshot of `.obs`
   - Plots: saved in the results directory


