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

### Compute Silhouette Scores
```bash
scsilhouette compute-silhouette \
    --h5ad-path data.h5ad \
    --cluster-header "cell_type" \
    --embedding-key "X_umap" \
    --organ "kidney" \
    --first-author "Lake" \
    --year "2023"
```

**Output directory:** `outputs_kidney_Lake_2023/`

**Output files:**
- `cell_type_silhouette_scores.csv` - Per-cell scores
- `cell_type_cluster_summary.csv` - Cluster statistics
- `cell_type_annotation.json` - Metadata

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

**Output:**
- Interactive HTML plot with side-by-side silhouette boxplots and F-score bars
- SVG for publication
- Summary statistics (median of medians)

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

## Use Cases

### Quality Control for Single-Cell Clustering

Silhouette scores measure how similar cells are to their assigned cluster compared to other clusters:
- **Score near +1**: Well-matched to cluster
- **Score near 0**: On the border between clusters
- **Score near -1**: Possibly assigned to wrong cluster

### Integration with NSForest

`scsilhouette` is designed to work seamlessly with [NSForest](https://github.com/JCVenterInstitute/NSForest) marker discovery:
- **Combined QC**: Assess both marker quality (F-scores) and cluster separation (silhouette)
- **Standardized outputs**: Matching directory structure `outputs_{organ}_{author}_{year}/`
- **Integrated visualizations**: Side-by-side comparison of metrics

### Workflow Integration

Used in the [sc-nsforest-qc-nf](https://github.com/NIH-NLM/sc-nsforest-qc-nf) Nextflow pipeline for scalable, parallelized analysis on CloudOS.

## Output Examples

### Silhouette Summary with F-scores
![Silhouette Summary](docs/images/silhouette_summary_example.png)

### Distribution Analysis
Side-by-side comparison of cluster sizes (log10 and raw) with silhouette scores to identify potential issues with very small or very large clusters.

## Output Structure

All outputs follow the pattern: `outputs_{organ}_{first_author}_{year}/{cluster_header}_*`
```
outputs_kidney_Lake_2023/
├── cell_type_silhouette_scores.csv
├── cell_type_silhouette_scores.json
├── cell_type_cluster_summary.csv
├── cell_type_cluster_summary.json
├── cell_type_annotation.json
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

A Docker container is available for reproducible analysis:
```bash
docker pull ghcr.io/nih-nlm/scsilhouette:1.0
```

## Docker Container

A Docker container is available for reproducible analysis:
```bash
docker pull ghcr.io/nih-nlm/scsilhouette:1.0
```

### Using on Apple Silicon (M1/M2/M3)

The container is built for linux/amd64. On Apple Silicon Macs, use the `--platform` flag:
```bash
# Pull the image
docker pull --platform linux/amd64 ghcr.io/nih-nlm/scsilhouette:1.0

# Run commands
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

Note: The container will run under x86_64 emulation on Apple Silicon, which is slightly slower but fully functional.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

Developed by the National Library of Medicine (NLM), National Institutes of Health (NIH).

Part of the Cell Knowledge Network project

## Contact

For questions or issues, please [open an issue](https://github.com/NIH-NLM/scsilhouette/issues) on GitHub.

## Related Projects

- [cell-kn](https://github.com/NIH-NLM/cell-kn) - NIH NLM Cell Knowledge Network
- [NSForest](https://github.com/JCVenterInstitute/NSForest) - Marker gene discovery
- [sc-nsforest-qc-nf](https://github.com/NIH-NLM/sc-nsforest-qc-nf) - Nextflow workflow combining NSForest and scsilhouette
- [cellxgene-harvester](https://github.com/NIH-NLM/cellxgene-harvester) - Single-cell data aggregation from CellxGene

## Citation

If you use scsilhouette in your research, please cite:
```
[Citation information will be added upon publication]
```
