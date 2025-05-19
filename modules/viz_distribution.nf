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

