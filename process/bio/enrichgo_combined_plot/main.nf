process BIO_ENRICHGO_COMBINED_PLOT {
    tag "combined_reactome"
    publishDir "${params.output}/bio", mode: 'copy', pattern: '*'
    cache false
    input:
    path lfc0_filter0_files
    path lfc0_filter05_files
    path lfc05_filter05_files
    output:
    path "*_reactome_bubble.png", emit: bubble
    path "reactome_terms.tsv", emit: terms

    script:
    """
    Rscript "${projectDir}/process/bio/enrichgo_combined_plot/templates/enrichgo_combined_plot.R"
    """

    stub:
    """
    touch kang_reactome_bubble.png
    touch perez_reactome_bubble.png
    touch GSE135779_CD14_11v11_reactome_bubble.png
    echo -e "dataset\tstrategy\tmethod\tDescription\tp.adjust\tGeneRatio\tCount" > reactome_terms.tsv
    """
}
