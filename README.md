# scsilhouette

[![Docs Status](https://img.shields.io/badge/docs-online-success)](https://nih-nlm.github.io/scsilhouette.github.io/)
[![CI Docs](https://github.com/NIH-NLM/scsilhouette.github.io/actions/workflows/docs.yml/badge.svg)](https://github.com/NIH-NLM/scsilhouette.github.io/actions/workflows/docs.yml)
[![License](https://img.shields.io/github/license/NIH-NLM/scsilhouette.github.io)](https://github.com/NIH-NLM/scsilhouette.github.io/blob/main/LICENSE)

Silhouette scoring and quality assessment for single-cell RNA-seq datasets.

---

## 🔧 Features

- Download `.h5ad` datasets
- Compute silhouette scores for multiple label keys
- Visualize cluster score distributions, summaries, QC boxplots
- Correlate silhouette scores with NS-Forest F-scores
- CLI-based, modular, and Nextflow-ready
- Auto-built documentation using GitHub Actions

---

## 📦 Install

```bash
git clone https://github.com/NIH-NLM/scsilhouette.github.io
cd scsilhouette.github.io
conda env create -f environment.yml
conda activate scsilhouette
pip install -e .
```

## 🚀 Command Line Examples

### Download

```bash
scsilhouette download \
  --url https://datasets.cellxgene.cziscience.com/example.h5ad \
  --output-dir ./data
```

### Compute

```bash
scsilhouette compute \
  --h5ad-path data/example.h5ad \
  --label-keys cell_type author_cell_type \
  --embedding-key X_pca \
  --output-dir results/ \
  --show-obs \
  --qc-correlations \
  --save-scores \
  --save-cluster-summary \
  --save-csv \
  --save-plots
```

### F-Score Correlation

```bash
scsilhouette correlate-fscore \
  --silhouette-csv results/cell_type_scores.csv \
  --fscore-csv nsforest/cell_type_fscore.csv \
  --output-dir results/
```

## 📄 Documentation

📘 Docs: https://nih-nlm.github.io/scsilhouette.github.io

To build locally:

```bash
sphinx-apidoc -f --separate -o docs/source/ src/scsilhouette
make -C docs html
open docs/build/html/index.html
```

## 🧪 Run Tests

```bash
pytest tests/
```

## 🧠 About

This package supports reproducible evaluation of cell clustering quality in scRNA-seq data, with direct comparison to NS-Forest marker gene-based F-scores.

## 📤 GitHub Pages Deployment
CI deploys docs from docs/build/html via GitHub Actions (.github/workflows/docs.yml) on every push to main.

## 📜 License
MIT License. © National Library of Medicine, NIH.


