#!/usr/bin/env nextflow

/**
 * Regenerate Dataset Summaries
 *
 * Lightweight workflow to regenerate _dataset_summary.csv files with ALL
 * cellxgene-harvester CSV columns. Reads existing cluster_summary and
 * nsforest results from a flat S3 results directory — no need to re-run
 * the full pipeline.
 *
 * Usage:
 *   nextflow run regenerate_summaries.nf \
 *       --datasets_csv s3://bucket/path/homo_sapiens_heart_harvester_final.csv \
 *       --organ heart \
 *       --results_dir s3://bucket/path/results \
 *       --outdir ./regenerated_summaries
 */

nextflow.enable.dsl=2

params.datasets_csv = null
params.organ        = null
params.results_dir  = null
params.outdir       = './results'
params.publish_mode = 'copy'

process regenerate_summary_process {
    tag "${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}", mode: params.publish_mode
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(cluster_summary), path(nsforest_results)

    output:
    path("*_dataset_summary.csv")

    script:
    def metadata_json = groovy.json.JsonOutput.toJson(meta)
    """
    echo '${metadata_json}' > metadata.json
    scsilhouette compute-summary-stats \
        --cluster-summary ${cluster_summary} \
        --nsforest-results ${nsforest_results} \
        --metadata metadata.json
    """
}

workflow {
    log.info "============================================"
    log.info "  Regenerate Dataset Summaries"
    log.info "============================================"
    log.info "datasets_csv : ${params.datasets_csv}"
    log.info "organ        : ${params.organ}"
    log.info "results_dir  : ${params.results_dir}"
    log.info "outdir       : ${params.outdir}"
    log.info "============================================"

    if (!params.datasets_csv) { exit 1, "ERROR: --datasets_csv is required" }
    if (!params.organ)        { exit 1, "ERROR: --organ is required" }
    if (!params.results_dir)  { exit 1, "ERROR: --results_dir is required (S3 path to existing results)" }

    csv_ch = Channel.fromPath(params.datasets_csv)
        .splitCsv(header: true)
        .filter { row ->
            def ref = row.reference?.trim()?.toLowerCase()
            if (ref in ['exclude', 'delete', 'merge', 'question']) {
                log.info "Skipping ${row.first_author} ${row.year} — reference='${row.reference}'"
                return false
            }
            if (!(ref in ['yes', 'no', 'unk'])) {
                log.warn "Skipping ${row.first_author} ${row.year} — unrecognised reference value '${row.reference}'"
                return false
            }
            return true
        }
        .map { row ->
            def meta = [:]
            row.each { k, v -> meta[k] = v }
            meta.organ = params.organ

            // Derive file names from naming convention
            // CloudOS results land under results_dir/results/ — auto-append if needed
            def base = params.results_dir.toString().replaceAll('/+$', '')
            def results_path = base.endsWith('/results') ? base : "${base}/results"
            def act = (row.author_cell_type ?: '').replace(' ', '_')
            def prefix = "${params.organ}_${row.first_author}_${row.year}_${act}"
            def cluster_summary = file("${results_path}/${prefix}_cluster_summary.csv")
            def nsforest_results = file("${results_path}/${prefix}_results.csv")

            log.info "Dataset: ${row.first_author} ${row.year} → ${prefix}"
            log.info "  cluster_summary: ${results_path}/${prefix}_cluster_summary.csv"
            log.info "  nsforest_results: ${results_path}/${prefix}_results.csv"

            tuple(meta, cluster_summary, nsforest_results)
        }

    regenerate_summary_process(csv_ch)
}
