process COMPUTE_SILHOUETTE {
  input:
  path h5ad_file
  val  label_keys
  val  embedding_key
  val  name

  output:
  tuple path("*_silhouette_scores.csv"), path("*_cluster_summary.csv"), val(name)

  publishDir "${params.outdir}/${name}/scores", mode: 'copy'
  container "ghcr.io/nih-nlm/scsilhouette:latest"

  script:
  """
  scsilhouette compute-silhouette \\
    --h5ad-path ${h5ad_file} \\
    --label-keys ${label_keys} \\
    --embedding-key ${embedding_key} \\
    --output-dir . \\
    --save-scores \\
    --save-cluster-summary
  """
}

