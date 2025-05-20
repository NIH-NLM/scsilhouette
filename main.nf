#!/usr/bin/env nextflow

params.datasets_csv         = "datasets.csv"
params.outdir               = "results"
params.metric               = "euclidean"
params.save_scores          = true
params.save_cluster_summary = true
params.save_annotation      = true

// include { COMPUTE_SILHOUETTE }  from './modules/computeSilhouette.nf'


workflow {

  csv_rows_ch = Channel
    .fromPath(params.datasets_csv)
    .splitCsv(header: true)

  h5ad_ch          = csv_rows_ch.map { row -> file(row.h5ad) }
  label_key_ch     = csv_rows_ch.map { row -> row.label_key }
  embedding_key_ch = csv_rows_ch.map { row -> row.embedding_key }

  compute_silhouette_process(
                 h5ad_ch,
		 label_key_ch,
		 embedding_key_ch,
		 params.output_dir,
		 params.metric,
		 params.save_scores,
		 params.save_cluster_summary,
		 params.save_annotation)
		 
}

process compute_silhouette_process {

    tag { "compute_silhouette_process" }


    input:
        path h5ad_file
        val  label_key
        val  embedding_key
        path output_dir
        val  metric
        val  save_scores
        val  save_cluster_summary
        val  save_annotation

    output:
        tuple path("silhouette_scores*.csv"), path("cluster_summary*.csv")

    script:
    """
    scsilhouette compute-silhouette \\
	--h5ad-path ${h5ad_file} \\
	--label-keys ${label_key} \\
	--embedding-key ${embedding_key} \\
	--output-dir . \\
	--metric ${metric} \\
	--save-scores \\
	--save-cluster-summary \\
	--save-annotation
    """
}