
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

