nextflow.enable.dsl=2

params.h5ad_path     = "data/sample.h5ad"
params.label_keys    = "cell_type"
params.embedding_key = "X_umap"
params.fscore_path   = "data/nsforest_scores.csv"
params.mapping_path  = "data/cell_type_cluster_map.csv"
params.outdir        = "results"
params.final_html    = "results/final_report.html"
params.final_pdf     = "results/final_report.pdf"

workflow {
    scores_ch = computeSilhouette(params.h5ad_path)

    summary_ch      = vizSummary(scores_ch.silhouette_scores)
    dotplot_ch      = vizDotplot(scores_ch.silhouette_scores, scores_ch.cluster_summary)
    heatmap_ch      = vizHeatmap(scores_ch.silhouette_scores, scores_ch.cluster_summary)
    dist_ch         = vizDistribution(scores_ch.silhouette_scores)
    dataset_ch      = vizDatasetSummary(scores_ch.silhouette_scores)

    mergeReport(
        summary_ch,
        dotplot_ch,
        heatmap_ch,
        dist_ch,
        dataset_ch
    )
}

process computeSilhouette {
    input:
    path h5ad_file from file(params.h5ad_path)

    output:
    path "*_silhouette_scores.csv", emit: silhouette_scores
    path "*_cluster_summary.csv",  emit: cluster_summary

    publishDir "${params.outdir}/scores", mode: 'copy'
    container: "ghcr.io/nih-nlm/scsilhouette:latest"

    script:
    """
    scsilhouette compute-silhouette \\
      --h5ad-path ${h5ad_file} \\
      --label-keys ${params.label_keys} \\
      --embedding-key ${params.embedding_key} \\
      --output-dir . \\
      --save-scores --save-cluster-summary
    """
}

process vizSummary {
    input:
    path silhouette_scores

    output:
    path "*", emit: summary_output

    publishDir "${params.outdir}/summary", mode: 'copy'
    container: "ghcr.io/nih-nlm/scsilhouette:latest"

    script:
    """
    scsilhouette viz-summary \\
      --silhouette-score-path ${silhouette_scores} \\
      --output-dir . \\
      --label ${params.label_keys} \\
      --score-col silhouette_score_euclidean \\
      --fscore-path ${params.fscore_path} \\
      --mapping-path ${params.mapping_path} \\
      --show --export-csv
    """
}

process vizDotplot {
    input:
    path silhouette_scores
    path cluster_summary

    output:
    path "*", emit: dotplot_output

    publishDir "${params.outdir}/dotplot", mode: 'copy'
    container: "ghcr.io/nih-nlm/scsilhouette:latest"

    script:
    """
    scsilhouette viz-dotplot \\
      --silhouette-score-path ${silhouette_scores} \\
      --cluster-summary-path ${cluster_summary} \\
      --output-dir . \\
      --label ${params.label_keys}
    """
}

process vizHeatmap {
    input:
    path silhouette_scores
    path cluster_summary

    output:
    path "*", emit: heatmap_output

    publishDir "${params.outdir}/heatmap", mode: 'copy'
    container: "ghcr.io/nih-nlm/scsilhouette:latest"

    script:
    """
    scsilhouette viz-heatmap \\
      --silhouette-score-path ${silhouette_scores} \\
      --cluster-summary-path ${cluster_summary} \\
      --output-dir . \\
      --label ${params.label_keys}
    """
}

process vizDistribution {
    input:
    path silhouette_scores

    output:
    path "*", emit: distribution_output

    publishDir "${params.outdir}/distribution", mode: 'copy'
    container: "ghcr.io/nih-nlm/scsilhouette:latest"

    script:
    """
    scsilhouette viz-distribution \\
      --silhouette-score-path ${silhouette_scores} \\
      --output-dir . \\
      --label ${params.label_keys}
    """
}

process vizDatasetSummary {
    input:
    path silhouette_scores

    output:
    path "*", emit: dataset_summary_output

    publishDir "${params.outdir}/dataset_summary", mode: 'copy'
    container: "ghcr.io/nih-nlm/scsilhouette:latest"

    script:
    """
    scsilhouette viz-dataset-summary \\
      --silhouette-score-path ${silhouette_scores} \\
      --output-dir .
    """
}

process mergeReport {
    input:
    path summary_files
    path dotplot_files
    path heatmap_files
    path dist_files
    path dataset_summary_files

    output:
    path "${params.final_html}"
    path "${params.final_pdf}"

    publishDir "${params.outdir}", mode: 'copy'
    container: "ghcr.io/nih-nlm/scsilhouette:latest"

    script:
    """
    cat <<EOF > report.html
    <html><head><title>scSilhouette Report</title></head><body>
    <h1>Dataset Summary</h1>
    <img src="dataset_summary.png"/>
    <h2>Silhouette Summary</h2>
    <img src="summary.png"/>
    <h2>Dotplot</h2>
    <img src="dotplot.png"/>
    <h2>Heatmap</h2>
    <img src="heatmap.png"/>
    <h2>Distribution</h2>
    <img src="distribution.png"/>
    </body></html>
    EOF

    cp report.html ${params.final_html}
    weasyprint ${params.final_html} ${params.final_pdf}
    """
}

