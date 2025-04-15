.. _scsilhouette.nsforest.cli:

nsforest CLI Usage
===================

CLI tool for validating NS-Forest marker outputs and correlating them with silhouette scores.

Example:

.. code-block:: bash

   scsilhouette nsforest-validate \
     --marker-csv data/nsforest_markers.csv \
     --summary-path results/cell_type_cluster_summary.csv \
     --output-dir results/nsforest/

CLI Help Output
---------------

.. literalinclude:: ../_cli_help/nsforest-validate.txt
   :language: text
