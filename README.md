# scsilhouette

![Docs Status](https://img.shields.io/badge/docs-online-success)
![License](https://img.shields.io/github/license/NIH-NLM/scsilhouette)

Silhouette scoring and F-score correlation for single-cell RNA-seq cluster validation.

---

## 🔧 Features

- Compute silhouette scores from `.h5ad` files
- Summarize silhouette stats per cluster
- Visualize F-score vs silhouette correlation
- Overlay NS-Forest F-scores using manual mappings
- Export merged summary stats for downstream analysis

---

## 📦 Installation

```bash
git clone https://github.com/NIH-NLM/scsilhouette
cd scsilhouette
conda env create -f environment.yml
conda activate scsilhouette
pip install -e .
```

---

## 🚀 Command Line Usage

### Compute Silhouette

```bash
scsilhouette compute-silhouette \
  --h5ad-path data/sample.h5ad \
  --label-keys cell_type \
  --embedding-key X_umap \
  --output-dir results/ \
  --save-scores --save-cluster-summary
```

### Visualize Summary

```bash
scsilhouette viz-summary \
  --silhouette-score-path results/cell_type_silhouette_scores.csv \
  --output-dir results/ \
  --label cell_type \
  --score-col silhouette_score_euclidean \
  --fscore-path data/nsforest_scores.csv \
  --mapping-path data/cell_type_cluster_map.csv \
  --show
```

![Summary Bar Chart](docs/source/_static/cell_type_summary_silhouette_score_euclidean_.png)

### Visualize F-score Correlation

```bash
scsilhouette viz-fscore \
  --fscore-path data/nsforest_scores.csv \
  --cluster-summary-path results/cell_type_cluster_summary.csv \
  --output-dir results/ \
  --label cell_type \
  --score-col silhouette_score_euclidean \
  --silhouette-stat median \
  --mapping-path data/cell_type_cluster_map.csv \
  --suffix smoke \
  --show --export-csv
```

![F-score Correlation](docs/source/_static/cell_type_fscore_vs_median_smoke.png)

---

## 🧠 Quadrant Interpretation

| Quadrant | Meaning                                                                 |
|----------|-------------------------------------------------------------------------|
| Q1       | 💡 High Silhouette, High F-score → Good clustering & markers            |
| Q2       | ⚠️ High Silhouette, Low F-score → Well-separated clusters, weak markers |
| Q3       | 🚧 Low Silhouette, Low F-score → Poor clustering and weak markers       |
| Q4       | 🤔 Low Silhouette, High F-score → Good markers, weak clustering         |

---

## 📄 Documentation

Full API and CLI documentation is auto-generated with [Sphinx](https://www.sphinx-doc.org/) using `autodoc`, `autosummary`, and `literalinclude` for CLI help injection.

Deployed using [GitHub Pages](https://pages.github.com/) at:

🔗 https://nih-nlm.github.io/scsilhouette/

All rendered figures used in this README are stored under `docs/source/_static/` and referenced within reStructuredText files for Sphinx processing.

---

## 🧪 Testing (Coming Soon)

```bash
pytest tests/
```

---

## 📄 License

MIT License © National Library of Medicine, NIH

