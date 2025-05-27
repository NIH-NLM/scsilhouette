#!/usr/bin/env nextflow

process merge_report_process {

    tag { "merge_report_process" }

    input:
        path summary
        path dotplot
        path heatmap
        path distribution
        path dataset_summary
        val  name

    output:
        path("${name}_report.html")
        path("${name}_report.pdf")

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

