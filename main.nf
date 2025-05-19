nextflow.enable.dsl=2

params.h5ad           = null
params.label_key      = null
params.embedding_key  = null
params.datasets       = null
params.datasets_csv   = null
params.outdir         = "results"

// --- Include process modules ---
include { COMPUTE_SILHOUETTE }        from './modules/compute_silhouette.nf'
include { VIZ_SUMMARY }               from './modules/viz_summary.nf'
include { VIZ_DOTPLOT }               from './modules/viz_dotplot.nf'
include { VIZ_HEATMAP }               from './modules/viz_heatmap.nf'
include { VIZ_DISTRIBUTION }          from './modules/viz_distribution.nf'
include { VIZ_DATASET_SUMMARY }       from './modules/viz_dataset_summary.nf'
include { MERGE_REPORT }              from './modules/merge_report.nf'

// --- Input Channel Setup ---
workflow {

    dataset_triples_ch = params.datasets_csv ?
        Channel
            .fromPath(params.datasets_csv)
            .splitCsv(header: true)
            .map { row -> 
                def name = row.h5ad.tokenize('/').last().replace('.h5ad','')
                tuple(file(row.h5ad), row.label_key, row.embedding_key, name)
            } :
        params.datasets ?
        Channel
            .from(params.datasets)
            .map { h5ad, label, embed ->
                def name = h5ad.tokenize('/').last().replace('.h5ad','')
                tuple(file(h5ad), label, embed, name)
            } :
        Channel
            .value([
                file(params.h5ad),
                params.label_key,
                params.embedding_key,
                params.h5ad.tokenize('/').last().replace('.h5ad','')
            ])

    // --- Run pipeline steps ---
    silhouette_out_ch     = COMPUTE_SILHOUETTE(dataset_triples_ch)
    summary_out_ch        = VIZ_SUMMARY(silhouette_out_ch)
    dotplot_out_ch        = VIZ_DOTPLOT(silhouette_out_ch)
    heatmap_out_ch        = VIZ_HEATMAP(silhouette_out_ch)
    distribution_out_ch   = VIZ_DISTRIBUTION(silhouette_out_ch)
    dataset_summary_ch    = VIZ_DATASET_SUMMARY(silhouette_out_ch)

    report_out_ch = MERGE_REPORT(
        summary_out_ch,
        dotplot_out_ch,
        heatmap_out_ch,
        distribution_out_ch,
        dataset_summary_ch
    )
}

