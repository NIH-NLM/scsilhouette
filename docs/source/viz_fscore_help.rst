.. _viz-fscore:

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

