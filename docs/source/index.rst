.. scsilhouette documentation master file, created by
   sphinx-quickstart on Mon Apr 15 12:00:00 2025.

Welcome to scsilhouette's documentation!
========================================

.. include:: ../README.md
   :parser: myst_parser.sphinx_

.. image:: _static/summary_silhouette.png
   :width: 600
   :align: center

Features
-----------

- Compute silhouette scores across label sets
- Visualize summary stats and correlation with F-scores
- NS-Forest cluster validation support
- CLI and programmatic API access
- Nextflow, Docker, and Conda-ready

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   process

.. toctree::
   :maxdepth: 2
   :caption: Command Line Interface

   scsilhouette.compute.cli
   scsilhouette.viz.cli
   scsilhouette.download.cli
   scsilhouette.nsforest.cli

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   modules

.. toctree::
   :maxdepth: 1
   :caption: Repository

   GitHub <https://github.com/NIH-NLM/scsilhouette>
   Project Website <https://nih-nlm.github.io/scsilhouette>

