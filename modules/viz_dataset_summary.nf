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

