#!/usr/bin/env nextflow

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
