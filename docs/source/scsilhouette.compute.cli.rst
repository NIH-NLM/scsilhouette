.. _scsilhouette.compute.cli:

compute-silhouette
===================

Computes silhouette scores for annotated `.h5ad` single-cell datasets using a specified embedding and cluster label.

.. code-block:: bash

   scsilhouette compute-silhouette \
     --h5ad-path data/sample.h5ad \
     --label-keys cell_type \
     --embedding-key X_umap \
     --output-dir results/ \
     --save-scores --save-cluster-summary

CLI Help Output
---------------

.. literalinclude:: ../_cli_help/compute-silhouette.txt
   :language: text
