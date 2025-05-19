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
