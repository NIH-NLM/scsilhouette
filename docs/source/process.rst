Using scsilhouette from the Command Line
=========================================

This page provides example commands for using `scsilhouette`'s CLI.

Compute Silhouette Scores
--------------------------

.. code-block:: bash

   scsilhouette compute-silhouette \
   --h5ad-path sample.h5ad \
   --label-keys cell_type \
   --embedding-key X_umap \
   --use-binary-genes \
   --gene-list-path \
   --metric euclidean \
   --save-scores \
   --save-cluster-summary \
   --save-annotation

Visualize Summary Statistics
-----------------------------

.. code-block:: bash

   scsilhouette viz-summary \
     --silhouette-score-path cell_type_silhouette_scores.csv \
     --label cell_type \
     --score-col silhouette_score \
     --fscore-path nsforest_scores.csv \
     --mapping-path cell_type_cluster_map.csv \
     --sort-by median

Visualize Correlation
---------------------

.. code-block:: bash

   scsilhouette viz-correlation \
     --cluster-summary-path cell_type_cluster_summary.csv
     --x-metric silhouette_score \
     --y-metrics median \
     --label cell_type \
     --fscore-path nsforest_scores.csv \
     --mapping-path cell_type_cluster_map.csv 

Visualize Dotplot
-----------------

.. code-block:: bash

   scsilhouette viz-dotplot \
   --h5ad-path sample.h5ad \
   --label-keys cell_type \
   --embedding-key X_umap 

Visualize Distribution
-----------------

.. code-block:: bash

   scsilhouette viz-distribution \
   --cluster-summary-path cluster_summary.csv \
   --label-keys cell_type 
