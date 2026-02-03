scsilhouette Documentation
==========================

Silhouette score analysis for single-cell clustering quality control.

.. toctree::
   :maxdepth: 2
   :caption: API Reference:

   modules

Installation
------------

.. code-block:: bash

   pip install scsilhouette

Quick Start
-----------

.. code-block:: bash

   scsilhouette compute-silhouette \
       --h5ad-path data.h5ad \
       --cluster-header "cell_type" \
       --embedding-key "X_umap" \
       --organ "kidney" \
       --first-author "Lake" \
       --year "2023"

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

