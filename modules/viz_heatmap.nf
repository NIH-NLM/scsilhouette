#!/usr/bin/env nextflow

process viz_heatmap_process {
    input:
        path h5ad_path
        val  groupby
        val  embedding_key

    output:
        path("*")

    script:
    """
    scsilhouette viz-heatmap \\
        --h5ad-path ${h5ad_path} \\
        --groupby ${groupby} \\
        --embedding-key ${embedding_key}
    """
}

