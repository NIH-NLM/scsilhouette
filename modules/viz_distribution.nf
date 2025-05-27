#!/usr/bin/env nextflow

process viz_distribution_process {

    tag ( "viz_distribution_process" }

    input:
        path cluster_summary
        val  label

    output:
        path("*")

    script:
    """
    scsilhouette viz-distribution \\
        --silhouette-score-path ${silhouette_scores} \\
        --label label
    """
}

