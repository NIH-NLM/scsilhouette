# scsilhouette

[![Build and Push Docker Image](https://github.com/NIH-NLM/scsilhouette/actions/workflows/docker-build.yml/badge.svg)](https://github.com/NIH-NLM/scsilhouette/actions/workflows/docker-build.yml)
[![Documentation Status](https://github.com/NIH-NLM/scsilhouette/actions/workflows/docs.yml/badge.svg)](https://nih-nlm.github.io/scsilhouette/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Silhouette score analysis for single-cell clustering quality control.

`scsilhouette` computes silhouette scores to assess clustering quality in single-cell RNA-seq data and provides integrated visualizations with NSForest marker discovery results.

## Documentation

Full documentation is available at: **[https://nih-nlm.github.io/scsilhouette/](https://nih-nlm.github.io/scsilhouette/)**

## Features

- **Compute silhouette scores** for single-cell clusters using embeddings (UMAP, t-SNE, etc.)
- **Integrated visualizations** with NSForest F-scores for comprehensive QC
- **Standardized output** structure compatible with NSForest workflows
- **Distribution analysis** of cluster sizes vs silhouette scores
- **Publication-ready plots** in SVG and interactive HTML formats
- **Ontology-based filtering** using UBERON, disease, and HsapDv JSON files from cellxgene-harvester

## Installation

### From PyPI (when published)
```bash
pip install scsilhouette
```

### From Source
```bash
git clone https://github.com/NIH-NLM/scsilhouette.git
cd scsilhouette
pip install -e .
```

### Using Conda Environment
```bash
git clone https://github.com/NIH-NLM/scsilhouette.git
cd scsilhouette
conda env create -f environment.yml
conda activate scsilhouette
pip install -e .
```

## Quick Start

### Compute Silhouette Scores (unfiltered)
```bash
scsilhouette compute-silhouette \
    --h5ad-path data.h5ad \
    --cluster-header "cell_type" \
    --embedding-key "X_umap" \
    --organ "kidney" \
    --first-author "Lake" \
    --year "2023"
```

### Compute Silhouette Scores (filtered to normal adult tissue)

First generate the three ontology JSON files using
[cellxgene-harvester](https://github.com/NIH-NLM/cellxgene-harvester):

```bash
# One-time per organ
cellxgene-harvester resolve-uberon kidney  > data/uberon_kidney.json

# One-time (organ-independent)
cellxgene-harvester resolve-disease normal > data/disease_normal.json
cellxgene-harvester resolve-hsapdv --min-age 15 > data/hsapdv_adult_15.json
```

Then pass all three to `compute-silhouette`:

```bash
scsilhouette compute-silhouette \
    --h5ad-path data.h5ad \
    --cluster-header "cell_type" \
    --embedding-key "X_umap" \
    --organ "kidney" \
    --first-author "Lake" \
    --year "2023" \
    --filter-normal \
    --uberon  data/uberon_kidney.json \
    --disease data/disease_normal.json \
    --hsapdv  data/hsapdv_adult_15.json
```

Filter stages applied in order:
1. **Tissue** — `tissue_ontology_term_id` matched against UBERON obo_ids
2. **Disease** — `disease_ontology_term_id` matched against disease obo_ids
3. **Age** — `development_stage_ontology_term_id` matched against HsapDv obo_ids (age threshold is encoded in the JSON at resolve time, not here)

**Output directory:** `outputs_kidney_Lake_2023/`

**Output files:**
- `cell_type_silhouette_scores.csv` — per-cell scores
- `cell_type_cluster_summary.csv` — per-cluster mean / median / std
- `cell_type_annotation.json` — dataset metadata and filter provenance

### Generate Summary Visualization with NSForest F-scores
```bash
scsilhouette viz-summary \
    --silhouette-score-path outputs_kidney_Lake_2023/cell_type_silhouette_scores.csv \
    --cluster-header "cell_type" \
    --organ "kidney" \
    --first-author "Lake" \
    --year "2023" \
    --fscore-path outputs_kidney_Lake_2023/cell_type_results.csv
```

### Create Embedding Dotplot
```bash
scsilhouette viz-dotplot \
    --h5ad-path data.h5ad \
    --embedding-key "X_umap" \
    --cluster-header "cell_type" \
    --organ "kidney" \
    --first-author "Lake" \
    --year "2023"
```

### Generate Distribution Plots
```bash
scsilhouette viz-distribution \
    --cluster-summary-path outputs_kidney_Lake_2023/cell_type_cluster_summary.csv \
    --cluster-header "cell_type" \
    --organ "kidney" \
    --first-author "Lake" \
    --year "2023"
```

## Nextflow Workflow Integration

`scsilhouette` is the scsilhouette component of the
[sc-nsforest-qc-nf](https://github.com/NIH-NLM/sc-nsforest-qc-nf)
Nextflow pipeline.  The pipeline runs NSForest marker discovery and silhouette
QC in parallel across every dataset listed in a harvester CSV, then publishes
results to [cell-kn](https://github.com/NIH-NLM/cell-kn) automatically.

### One-time setup

```bash
# Clone the workflow
git clone https://github.com/NIH-NLM/sc-nsforest-qc-nf.git
cd sc-nsforest-qc-nf

# Generate organ-specific ontology JSON files (kidney example)
cellxgene-harvester resolve-uberon  kidney > data/uberon_kidney.json
cellxgene-harvester resolve-disease normal > data/disease_normal.json
cellxgene-harvester resolve-hsapdv  --min-age 15 > data/hsapdv_adult_15.json
```

### Running the workflow

```bash
nextflow run main.nf \
    --datasets_csv   data/homo_sapiens_kidney_harvester_final.csv \
    --organ          kidney \
    --uberon_json    data/uberon_kidney.json \
    --disease_json   data/disease_normal.json \
    --hsapdv_json    data/hsapdv_adult_15.json \
    --github_token   "$(cat ~/.github_token)" \
    -c               configs/macamd64.config
```

> **Security note — `--github_token`:**
> The token is used by the `publish_results` module to push a results branch
> and open a pull request against `NIH-NLM/cell-kn`.  Never hardcode the token
> in a config file or commit it to version control.  The recommended approaches
> are:
>
> **Option 1 — params JSON file (preferred)**
> ```json
> { "github_token": "ghp_xxxxxxxxxxxxxxxxxxxx" }
> ```
> ```bash
> nextflow run main.nf -params-file params.json ...
> ```
> Add `params.json` to `.gitignore`.
>
> **Option 2 — environment variable**
> ```bash
> export GITHUB_TOKEN="ghp_xxxxxxxxxxxxxxxxxxxx"
> nextflow run main.nf --github_token "$GITHUB_TOKEN" ...
> ```
>
> **Option 3 — Nextflow secrets (Nextflow Tower / Seqera Platform)**
> ```bash
> nextflow secrets set GITHUB_TOKEN ghp_xxxxxxxxxxxxxxxxxxxx
> ```
> Then reference it in `nextflow.config` as `params.github_token = secrets.GITHUB_TOKEN`.

### Datasets CSV columns

The `--datasets_csv` file must have these columns:

| Column | Description |
|--------|-------------|
| `h5ad_file` | Absolute or relative path to the `.h5ad` file |
| `first_author` | Surname of first author (used in output naming) |
| `year` | Publication year |
| `author_cell_type` | Column name in `adata.obs` containing cluster labels |
| `embedding` | Embedding key, e.g. `X_umap` |
| `disease` | Disease label string for annotation output |
| `filter_normal` | `True` or `False` — whether to apply the three-stage ontology filter |
| `reference` | Processing flag (see table below) |

**`reference` column values:**

| Value | Behaviour |
|-------|-----------|
| `yes` | Process this dataset |
| `no` | Process this dataset (not a reference atlas) |
| `unk` | Process this dataset (reference status unknown) |
| `exclude` | Skip — explicitly excluded |
| `delete` | Skip — marked for deletion |
| `merge` | Skip — to be merged with another dataset |
| `question` | Skip — under review |

### Published output structure

After the workflow completes, a pull request is opened against `main` on
`NIH-NLM/cell-kn` with results written to:

```
data/prod/{organ}/
├── nsforest/
│   └── {organ}_{first_author}_{year}/    ← one directory per processed dataset
└── scsilhouette/
    └── {organ}_{first_author}_{year}/    ← one directory per processed dataset
```

The cellxgene-harvester outputs and human-in-the-loop annotation are deposited
manually into the same `data/prod/{organ}/` tree before the PR is merged.

## Use Cases

### Quality Control for Single-Cell Clustering

Silhouette scores measure how similar cells are to their assigned cluster
compared to other clusters:

- **Score near +1** — well-matched to cluster
- **Score near 0** — on the border between clusters
- **Score near −1** — possibly assigned to the wrong cluster

### Integration with NSForest

`scsilhouette` is designed to work seamlessly with
[NSForest](https://github.com/JCVenterInstitute/NSForest) marker discovery:

- **Combined QC** — assess both marker quality (F-scores) and cluster separation (silhouette)
- **Standardized outputs** — matching directory structure `outputs_{organ}_{author}_{year}/`
- **Integrated visualizations** — side-by-side comparison of both metrics

## Output Structure

All outputs follow the pattern `outputs_{organ}_{first_author}_{year}/{cluster_header}_*`:

```
outputs_kidney_Lake_2023/
├── cell_type_silhouette_scores.csv
├── cell_type_cluster_summary.csv
├── cell_type_annotation.json
├── cell_type_dataset_summary.csv
├── cell_type_silhouette_fscore_summary.html
├── cell_type_silhouette_fscore_summary.svg
├── cell_type_dotplot_X_umap.html
├── cell_type_dotplot_X_umap.svg
├── cell_type_distribution_log10.html
├── cell_type_distribution_log10.svg
├── cell_type_distribution_raw.html
└── cell_type_distribution_raw.svg
```

## Docker Container

```bash
docker pull ghcr.io/nih-nlm/scsilhouette:1.0
```

### Using on Apple Silicon (M1/M2/M3)

The container is built for `linux/amd64`. On Apple Silicon Macs:

```bash
docker pull --platform linux/amd64 ghcr.io/nih-nlm/scsilhouette:1.0

docker run --platform linux/amd64 -v $(pwd)/data:/data \
  ghcr.io/nih-nlm/scsilhouette:1.0 \
  compute-silhouette \
  --h5ad-path /data/test.h5ad \
  --cluster-header "cell_type" \
  --embedding-key "X_umap" \
  --organ "kidney" \
  --first-author "Lake" \
  --year "2023"
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License — see the LICENSE file for details.

## Acknowledgments

Developed by the National Library of Medicine (NLM), National Institutes of Health (NIH).
Part of the Cell Knowledge Network project.

## Contact

For questions or issues, please [open an issue](https://github.com/NIH-NLM/scsilhouette/issues) on GitHub.

## Related Projects

- [cell-kn](https://github.com/NIH-NLM/cell-kn) — NIH NLM Cell Knowledge Network
- [NSForest](https://github.com/JCVenterInstitute/NSForest) — Marker gene discovery
- [sc-nsforest-qc-nf](https://github.com/NIH-NLM/sc-nsforest-qc-nf) — Nextflow workflow combining NSForest and scsilhouette
- [cellxgene-harvester](https://github.com/NIH-NLM/cellxgene-harvester) — Single-cell data aggregation from CellxGene

## Citation

If you use scsilhouette in your research, please cite:
```
[Citation information will be added upon publication]
```
