#!/usr/bin/env nextflow

params.data_dir         = "data"
params.output_dir       = "results"
params.label_keys       = "cell_type"
params.embedding_key    = "X_umap"
params.metric           = "euclidean"
params.show_obs         = true
params.save_scores      = true
params.save_csv         = true
params.save_plots       = true
params.save_cluster_summary = true
params.qc_correlations  = true

process RunScSilhouette {
    tag "$h5ad_file.baseName"

    input:
    path h5ad_file
    path nsforest_file

    output:
    path "${params.output_dir}/${h5ad_file.baseName}", emit: results

    script:
    def output_subdir = "${params.output_dir}/${h5ad_file.baseName}"

    """
    mkdir -p ${output_subdir}
    
    scsilhouette compute \\
      --h5ad-path ${h5ad_file} \\
      --label-keys ${params.label_keys} \\
      --embedding-key ${params.embedding_key} \\
      --output-dir ${output_subdir} \\
      --nsforest-path ${nsforest_file} \\
      ${params.show_obs ? "--show-obs" : ""} \\
      ${params.save_scores ? "--save-scores" : ""} \\
      ${params.save_csv ? "--save-csv" : ""} \\
      ${params.save_plots ? "--save-plots" : ""} \\
      ${params.save_cluster_summary ? "--save-cluster-summary" : ""} \\
      ${params.qc_correlations ? "--qc-correlations" : ""}
    """
}

workflow {
    Channel
      .fromFilePairs("${params.data_dir}/*_ann_*.csv", size: 1)
      .map { fscore -> 
          def h5ad = file(fscore.name.split("_ann_")[0] + ".h5ad")
          tuple(h5ad, fscore)
      }
      | RunScSilhouette
}

