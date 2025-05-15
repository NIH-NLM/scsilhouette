#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.h5ads           = "${baseDir}/data/*.h5ad"
params.label_key       = "celltype_level3_fullname"
params.embedding_key   = "X_umap"
params.metric          = "euclidean"
params.fscore_path     = null
params.mapping_path    = null
params.output_dir      = "results"
params.container       = "ghcr.io/myorg/scsilhouette:latest" // or your Dockerhub

workflow {
  Channel
    .fromPath(params.h5ads)
    .map { file -> 
      def tag = file.baseName.replace('.h5ad','')
      tuple(file, tag)
    }
    | compute_silhouette
    | visualize_all
    | make_report
}

process compute_silhouette {

  publishDir "${params.output_dir}", mode: 'copy'

  container "${params.container}"

  input:
  path h5ad
  val tag

  output:
  path "results/silhouette_scores_${tag}_*.csv"
  path "results/cluster_summary_${tag}_*.csv"

  script:
  """
  scsilhouette compute-silhouette \\
    --h5ad-path $h5ad \\
    --label-key ${params.label_key} \\
    --embedding-key ${params.embedding_key} \\
    --metric ${params.metric} \\
    --output-dir results \\
    --save-scores --save-cluster-summary --save-annotation
  """
}

process visualize_all {

  publishDir "${params.output_dir}", mode: 'copy'

  container "${params.container}"

  input:
  path silhouette_scores
  path cluster_summary

  script:
  def tag = silhouette_scores.name.replaceFirst(/silhouette_scores_/, "").replaceFirst(/.csv$/, "")
  def suffix = "_${tag}"

  """
  scsilhouette viz-summary \\
    --silhouette-score-path $silhouette_scores \\
    --output-dir results \\
    --label ${params.label_key} \\
    --score-col silhouette_score \\
    --sort-by median \\
    --suffix "$suffix"

  scsilhouette viz-correlation \\
    --cluster-summary-path $cluster_summary \\
    --output-dir results \\
    --silhouette-metrics mean,median,std \\
    --score-col f_score \\
    --suffix "$suffix"

  scsilhouette viz-dataset-summary \\
    --cluster-summary-path $cluster_summary \\
    --output-dir results \\
    --label ${params.label_key} \\
    --suffix "$suffix"

  scsilhouette viz-distribution \\
    --input-path $cluster_summary \\
    --score-col mean_silhouette \\
    --label ${params.label_key} \\
    --output-dir results \\
    --suffix "$suffix"

  scsilhouette viz-distribution \\
    --input-path $cluster_summary \\
    --score-col f_score \\
    --label ${params.label_key} \\
    --output-dir results \\
    --suffix "$suffix"

  scsilhouette viz-dotplot \\
    --score-path $silhouette_scores \\
    --label ${params.label_key} \\
    --score-col silhouette_score \\
    --output-dir results \\
    --suffix "$suffix"

  scsilhouette viz-heatmap \\
    --score-path $silhouette_scores \\
    --label ${params.label_key} \\
    --score-col silhouette_score \\
    --output-dir results \\
    --suffix "$suffix"
  """
}

process make_report {

  publishDir "${params.output_dir}/report", mode: 'copy'

  container "${params.container}"

  input:
  path("results/*summary*.png")
  path("results/*heatmap*.png")
  path("results/*dotplot*.png")
  path("results/*distribution*.png")
  path("results/*correlation*.png")
  path("results/*summary*.csv")

  output:
  path "report.html"

  script:
  """
  python src/scsilhouette/report_generator.py --input-dir results --output report.html
  """
}
