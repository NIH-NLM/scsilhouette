# scsilhouette

**scsilhouette** is a Python package to compute silhouette scores for single-cell RNA-seq datasets.  
It supports `.h5ad` files, flexible label fields, modular plotting, QC metrics, and association studies with external cluster metrics such as NS-Forest F-Scores.

---

## 🔧 Installation

```bash
conda env create -f environment.yml
conda activate scrnaseq_silhouette
pip install -e .
```

## 📥 Download Dataset

```bash
scsilhouette download \
  --url https://datasets.cellxgene.cziscience.com/5daeaafe-c79e-4ee4-a9f0-ddf6649adc21.h5ad \
  --output-dir data/
```

## 📊 Compute Silhouette Scores

```bash
scsilhouette compute \
  --h5ad-path data/my_data.h5ad \
  --label-keys cell_type \
  --embedding-key X_pca \
  --output-dir results/ \
  --save-scores \
  --save-cluster-summary \
  --save-csv \
  --save-plots \
  --qc-correlations \
  --show-obs

```

## 📂 F-Score Input Format

To run correlation studies against external metrics like NS-Forest:

Example F-Score CSV:

| label_col | fscore |
| ----------| ------ |
| Alveolar Type 2 | 0.87 |
| Basal	| 0.91 |
| Ciliated | 0.76 |

label_col: Must match label keys used in --label-keys

fscore: Numeric column used for correlation

## 📈 Visualizations

* Silhouette score distribution per cluster

* Cluster summary (mean silhouette per label)

* QC plots: silhouette vs. nCount_RNA, nFeature_RNA

* F-score vs. Silhouette Score (if provided)

##  📚 Build Sphinx Docs

From the root directory of this repositiory *`scsilhouette`*

```bash
cd docs
sphinx-apidoc -o source ../src/scsilhouette
make html
open build/html/index.html
```

and you see your beautiful documentation

##🧪 Run Unit Tests

```bash
pytest tests/
```

** 🌐 Online Documentation

📖 https://nih-nlm.github.io/scsilhouette.github.io/

