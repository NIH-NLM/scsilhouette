#!/usr/bin/env nextflow

process viz_dataset_summary_process {

    tag ( "viz_dataset_summary_process" }

    input:
        path cluster_summary_path
        val  label
    
    output:
        path("*")

    script:
    """
    scsilhouette viz-dataset-summary \\
    --cluster-summary-path ${cluster_summary_path} \\
    --label ${label}
    """
}

