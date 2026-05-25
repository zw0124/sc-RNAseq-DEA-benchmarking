process BIO_ENRICHGO_COMBINED_PLOT {
    tag "combined_reactome"
    publishDir "${params.output}/bio", mode: 'copy', pattern: '*.png'
    cache false
    input:
    path lfc0_filter0_files
    path lfc0_filter05_files
    path lfc05_filter05_files
    output:
    path "combined_reactome_bubble.png", emit: bubble
    path "combined_reactome_terms.tsv", emit: terms

    script:
    """
    Rscript "${projectDir}/process/bio/enrichgo_combined_plot/templates/enrichgo_combined_plot.R"
    """

    stub:
    """
    touch combined_reactome_bubble.png
    echo -e "strategy\tmethod\tDescription\tp.adjust\tGeneRatio\tCount" > combined_reactome_terms.tsv
    """
}
