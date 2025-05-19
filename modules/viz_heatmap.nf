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

