#!/usr/bin/env nextflow

process viz_dotplot_process {

    tag ( "viz_dotplot_process" }
  
    input:
        path h5ad_path
        val  groupby
        val  embedding_key

    output:
        path("*")

    script:
    """
    scsilhouette viz-dotplot \\
        --h5ad-path h5ad_path \\
        --groupby groupby \\
        --embedding-key embedding_key
    """
}

