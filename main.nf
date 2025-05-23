#!/usr/bin/env nextflow

include { compute_silhouette_process }  from './modules/compute_silhouette.nf'


workflow {

//  csv_rows_ch = Channel
//    .fromPath(params.datasets_csv)
//    .splitCsv(header: true)

//  h5ad_ch          = csv_rows_ch.map { row -> file(row.h5ad) }
//  label_key_ch     = csv_rows_ch.map { row -> row.label_key }
//  embedding_key_ch = csv_rows_ch.map { row -> row.embedding_key }

  h5ad_ch          = params.h5ad_ch
  label_key_ch     = params.label_key_ch
  embedding_key_ch = params.embedding_key_ch
  
  compute_silhouette_process(
                 h5ad_ch,
		 label_key_ch,
		 embedding_key_ch,
		 params.outdir,
		 params.metric,
		 params.save_scores,
		 params.save_cluster_summary,
		 params.save_annotation)
		 


}

