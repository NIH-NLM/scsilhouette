Process Overview
================

The `scsilhouette` Python package computes silhouette scores on single-cell RNA-seq data
and performs visualizations and association analysis with external metrics like F-scores.

Usage Workflow
--------------

1. **Download** dataset from a `cellxgene` `.h5ad` URL:
   ```
   scsilhouette download --url <H5AD_URL> --output-dir data/
   ```

2. **Compute** silhouette scores:
   ```
   scsilhouette compute --h5ad-path data/sample.h5ad --label-keys cell_type --embedding-key X_pca --output-dir results/
   ```

3. **Analyze** associations with NS-Forest F-scores:
   ```
   scsilhouette compute --h5ad-path ... --fscore-path ... --output-dir results/
   ```

F-Score Format
--------------

The expected format of the NS-Forest F-score CSV:
- Columns: `cluster`, `fscore`, `marker_count`
- One row per cluster-level result

This file is mapped by cluster name and used for correlation plots and analysis.
