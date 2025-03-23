# scsilhouette

Silhouette score computation for assessing single-cell dataset cluster quality, using `.h5ad` inputs. Supports automated quality control, visualization, CLI usage, and reproducible integration into Nextflow workflows.


## Installation

Clone and install with conda:

```bash
conda env create -f environment.yml
conda activate scsilhouette
pip install -e .
```

## CLI Usage

The Command Line Interface (CLI) is as follows:

```bash
scsilhouette compute \
  --input path/to/data.h5ad \
  --label_field cell_type \
  --embedding_key X_pca \
  --output_dir results/
```

Optional arguments:

* --label_field: the cell annotation (e.g., cell_type, author_cell_type)
* --embedding_key: e.g., X_pca, X_umap, or another key in .obsm
* --output_dir: folder for results


## Download a file

To download a file provide the download function with the specific URL.
To construct the URL, you add the version id to the end of the API.

The existing API is:
```bash
https://api.cellxgene.cziscience.com/curation/v1/dataset_versions/{dataset_version_id}
```

The dataset version id for one of the smallest datasets at cellxgene is:

```bash
17616a35-1b3d-4754-9553-8ddb71267d63
```

then to append the specific version ID the complete URL is

```bash
https://api.cellxgene.cziscience.com/curation/v1/dataset_versions/17616a35-1b3d-4754-9553-8ddb71267d63
```

so to download this **dataset_version_id**

```bash
python -m scsilhouette.download https://api.cellxgene.cziscience.com/curation/v1/dataset_versions/17616a35-1b3d-4754-9553-8ddb71267d63
```

Providing the function with the complete URL.
This saves the file to the current directory.

## Output Artifacts

Upon execution, the following artifacts are written to the output directory

### CSV

* *`silhouette_scores.csv`*: silhouette scores per cell
* *`cluster_summary.csv`* : per cluster stats including
  * *`label`*
  * *`mean-silhouette`*
  * *`cluster-size`*

### Plots

* *`silhouette_distribution.png`*:  histogram silhouette values
* *`cluster_summary.png'*:  bar plot of mean silhouette score per cluster
* *`size_vs_score.png'* : scatter plot of clustersize vs mean silhouette


## Testing

Run tests with

```bash
pytest tests/
```

Test modules include:

* *`test_compute.py`*
* *`test_viz.py`*
* *`test_download.py`*

Each module checks correctness, plotting, and downloading behavior using a fixed test .h5ad dataset.


## Docker

Build container with :

```bash
docker build -t scsilhouette .
```

## Development Notes

### Project Bootstrapping Steps

1. Created *`src/scsilhouette/`* package structure

2. Added CLI via *`argparse`* in *`cli.py`*

3. Added *`compute.py`*, *`viz.py`*, and *`download.py`* modules

4. Installed *`pooch`* for robust *`URL`* downloading

5. Integrated plotting via *`seaborn`* and *`matplotlib`*

6. Included *`setup.cfg`* and *`pyproject.toml`* for packaging

7. Configured *`console_scripts`* for *`scsilhouette`* CLI

8. Verified CLI with: *`scsilhouette --help`*

9. Pushed to GitHub at *`https://github.com/nih-nlm/scsilhouette`*

10. Planned for:

    * Sphinx documentation
    * GitHub Pages via docs/
    * CI tests with GitHub Actions





