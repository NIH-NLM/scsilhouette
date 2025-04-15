.. _scsilhouette.viz.cli:

viz-fscore
==========

This command visualizes the correlation between NS-Forest F-scores and silhouette-based cluster summary scores.

.. code-block:: bash

   scsilhouette viz-fscore \
     --fscore-path data/nsforest_results.csv \
     --cluster-summary-path results/cell_type_cluster_summary.csv \
     --output-dir results/viz \
     --label cell_type \
     --score-col silhouette_score_euclidean \
     --silhouette-stat median \
     --mapping-path data/cell_type_cluster_map.csv \
     --suffix test \
     --show \
     --export-csv

CLI Help Output
---------------

.. literalinclude:: ../_cli_help/viz-fscore.txt
   :language: text

.. figure:: /_static/cell_type_fscore_vs_median_smoke.png
   :alt: F-score vs. Median silhouette
   :figwidth: 80%
   :align: center

   Scatter plot showing correlation between NS-Forest F-scores and median silhouette scores.

viz-summary
===========

This command summarizes silhouette scores per cluster label and optionally overlays NS-Forest F-scores.

.. code-block:: bash

   scsilhouette viz-summary \
     --silhouette-score-path results/scores.csv \
     --output-dir results/summary \
     --label cell_type \
     --score-col silhouette_score_euclidean \
     --fscore-path data/nsforest_results.csv \
     --mapping-path data/cell_type_cluster_map.csv \
     --show

CLI Help Output
---------------

.. literalinclude:: ../_cli_help/viz-summary.txt
   :language: text

.. figure:: /_static/cell_type_summary_silhouette_score_euclidean_.png
   :alt: Summary of silhouette scores by cluster
   :figwidth: 90%
   :align: center

   Bar chart of mean silhouette scores with error bars and overlaid F-scores for each cluster label.

