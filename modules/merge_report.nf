process mergeReport {
  input:
  path summary_files
  path dotplot_files
  path heatmap_files
  path distribution_files
  path dataset_summary_files
  val  name

  output:
  path("${name}_report.html")
  path("${name}_report.pdf")

  publishDir "${params.outdir}/${name}/report", mode: 'copy'
  container "ghcr.io/nih-nlm/scsilhouette:latest"

  script:
  """
  cat <<EOF > ${name}_report.html
  <html><head><title>${name} Report</title></head><body>
  <h1>Dataset: ${name}</h1>
  <img src="dataset_summary.png"/>
  <h2>Summary</h2><img src="summary.png"/>
  <h2>Dotplot</h2><img src="dotplot.png"/>
  <h2>Heatmap</h2><img src="heatmap.png"/>
  <h2>Distribution</h2><img src="distribution.png"/>
  </body></html>
  EOF

  weasyprint ${name}_report.html ${name}_report.pdf
  """
}
