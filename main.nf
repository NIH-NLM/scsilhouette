nextflow.enable.dsl=2

// Parameters
params.h5ad          = null
params.label_key     = null
params.embedding_key = null
params.datasets      = null
params.datasets_csv  = null
params.outdir        = "results"

// Determine input source
def dataset_triples_ch =
    params.datasets_csv ?
        Channel
            .fromPath(params.datasets_csv)
            .splitCsv(header: true)
            .map { row ->
                def name = row.h5ad.tokenize('/').last().replace('.h5ad','')
                tuple(file(row.h5ad), row.label_key, row.embedding_key, name)
            } :
    params.datasets ?
        Channel
            .from(params.datasets)
            .map { h5ad, label, embed ->
                def name = h5ad.tokenize('/').last().replace('.h5ad','')
                tuple(file(h5ad), label, embed, name)
            } :
        Channel
            .value([
                file(params.h5ad),
                params.label_key,
                params.embedding_key,
                params.h5ad.tokenize('/').last().replace('.h5ad','')
            ])

workflow {
  dataset_triples_ch
    .each { h5ad_file, label_keys, embedding_key, name ->

      scores = computeSilhouette(h5ad_file, label_keys, embedding_key, name)

      summary = vizSummary(scores[0], label_keys, name)
      dotplot = vizDotplot(scores[0], scores[1], label_keys, name)
      heatmap = vizHeatmap(scores[0], scores[1], label_keys, name)
      dist    = vizDistribution(scores[0], label_keys, name)
      dataset = vizDatasetSummary(scores[0], name)

      mergeReport(summary, dotplot, heatmap, dist, dataset, name)
    }
}

// --- Processes ---
process computeSilhouette {
  input:
  path h5ad_file
  val  label_keys
  val  embedding_key
  val  name

  output:
  tuple path("*_silhouette_scores.csv"), path("*_cluster_summary.csv")

  publishDir "${params.outdir}/${name}/scores", mode: 'copy'
  container "ghcr.io/nih-nlm/scsilhouette:latest"

  script:
  """
  scsilhouette compute-silhouette \\
    --h5ad-path ${h5ad_file} \\
    --label-keys ${label_keys} \\
    --embedding-key ${embedding_key} \\
    --output-dir . \\
    --save-scores --save-cluster-summary
  """
}

process vizSummary {
  input:
  path silhouette_scores
  val  label_keys
  val  name

  output:
  path("*")

  publishDir "${params.outdir}/${name}/summary", mode: 'copy'
  container "ghcr.io/nih-nlm/scsilhouette:latest"

  script:
  """
  scsilhouette viz-summary \\
    --silhouette-score-path ${silhouette_scores} \\
    --output-dir . \\
    --label ${label_keys} \\
    --score-col silhouette_score_euclidean \\
    --fscore-path data/nsforest_scores.csv \\
    --mapping-path data/cell_type_cluster_map.csv \\
    --show --export-csv
  """
}

process vizDotplot {
  input:
  path silhouette_scores
  path cluster_summary
  val  label_keys
  val  name

  output:
  path("*")

  publishDir "${params.outdir}/${name}/dotplot", mode: 'copy'
  container "ghcr.io/nih-nlm/scsilhouette:latest"

  script:
  """
  scsilhouette viz-dotplot \\
    --silhouette-score-path ${silhouette_scores} \\
    --cluster-summary-path ${cluster_summary} \\
    --output-dir . \\
    --label ${label_keys}
  """
}

process vizHeatmap {
  input:
  path silhouette_scores
  path cluster_summary
  val  label_keys
  val  name

  output:
  path("*")

  publishDir "${params.outdir}/${name}/heatmap", mode: 'copy'
  container "ghcr.io/nih-nlm/scsilhouette:latest"

  script:
  """
  scsilhouette viz-heatmap \\
    --silhouette-score-path ${silhouette_scores} \\
    --cluster-summary-path ${cluster_summary} \\
    --output-dir . \\
    --label ${label_keys}
  """
}

process vizDistribution {
  input:
  path silhouette_scores
  val  label_keys
  val  name

  output:
  path("*")

  publishDir "${params.outdir}/${name}/distribution", mode: 'copy'
  container "ghcr.io/nih-nlm/scsilhouette:latest"

  script:
  """
  scsilhouette viz-distribution \\
    --silhouette-score-path ${silhouette_scores} \\
    --output-dir . \\
    --label ${label_keys}
  """
}

process vizDatasetSummary {
  input:
  path silhouette_scores
  val  name

  output:
  path("*")

  publishDir "${params.outdir}/${name}/dataset_summary", mode: 'copy'
  container "ghcr.io/nih-nlm/scsilhouette:latest"

  script:
  """
  scsilhouette viz-dataset-summary \\
    --silhouette-score-path ${silhouette_scores} \\
    --output-dir .
  """
}

process mergeReport {
  input:
  path summary_files
  path dotplot_files
  path heatmap_files
  path distribution_files
  path dataset_summary_files
  val  name

  output:
  path("${name}_report.html")
  path("${name}_report.pdf")

  publishDir "${params.outdir}/${name}/report", mode: 'copy'
  container "ghcr.io/nih-nlm/scsilhouette:latest"

  script:
  """
  cat <<EOF > ${name}_report.html
  <html><head><title>${name} Report</title></head><body>
  <h1>Dataset: ${name}</h1>
  <img src="dataset_summary.png"/>
  <h2>Summary</h2><img src="summary.png"/>
  <h2>Dotplot</h2><img src="dotplot.png"/>
  <h2>Heatmap</h2><img src="heatmap.png"/>
  <h2>Distribution</h2><img src="distribution.png"/>
  </body></html>
  EOF

  weasyprint ${name}_report.html ${name}_report.pdf
  """
}

// ↓↓↓ All processes remain unchanged from previous full version ↓↓↓

