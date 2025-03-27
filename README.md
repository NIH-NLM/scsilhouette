# scsilhouette

**scsilhouette** is a Python package to compute silhouette scores for single-cell RNA-seq datasets.  
It supports `.h5ad` inputs, correlation with external F-Scores, integrated plotting, and CLI operations.

---

## 🔧 Installation

```bash
conda env create -f environment.yml
conda activate scrnaseq_silhouette
```

---

## 🚀 Usage

### Download an `.h5ad` Dataset

```bash
scsilhouette download --url <dataset_url> --output-dir data/
```

### Compute Silhouette Scores

```bash
scsilhouette compute \
  --h5ad-path data/my_file.h5ad \
  --label-keys cell_type \
  --embedding-key X_pca \
  --output-dir results/ \
  --save-scores --save-plots --save-csv --qc-correlations
```

### Correlate with NSForest F-Scores

```bash
scsilhouette correlate-fscore \
  --silhouette-csv results/cell_type_cluster_summary.csv \
  --fscore-csv path/to/nsforest.csv \
  --label-col cell_type \
  --output-dir results/
```

---

## 📦 Project Structure

```
scsilhouette/
├── src/scsilhouette/
│   ├── compute.py
│   ├── viz.py
│   ├── download.py
│   ├── cli.py
├── tests/
│   ├── test_compute.py
│   ├── test_viz.py
│   ├── test_cli.py
├── docs/
│   ├── source/
│   │   ├── conf.py
│   │   ├── index.rst
│   │   ├── process.rst  ← Dev/test steps
│   │   ├── modules.rst
```

---

## 🧪 Development and Testing

```bash
pytest tests/
```

---

## 📄 Documentation

Build documentation using Sphinx:

```bash
cd docs
make html
```

Auto-deploys via GitHub Actions to:  
[https://nih-nlm.github.io/scsilhouette.github.io/](https://nih-nlm.github.io/scsilhouette.github.io/)

---
