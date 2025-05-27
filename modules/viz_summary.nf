#!/usr/bin/env nextflow

process viz_summary_process {

    tag { "viz_summary_process" }
  
    input:
        path silhouette_scores_path
        val  label
        val  score_col
        val  sort_by

    output:
        path("*")

    script:
    """
    scsilhouette viz-summary \\
        --silhouette-score-path ${silhouette_scores_path} \\
        --label ${label} \\
        --sort_by ${sort_by} \\
        --score-col ${score_col} \\
        --sort-by ${sort_by}
    """
}


