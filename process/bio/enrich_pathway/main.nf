process BIO_ENRICH_PATHWAY {
    tag "reactome:${meta.scenario}"

    input:
    tuple val(meta), path(filtered_genes), path(universe_file)
    output:
    tuple val(meta), path("bio_enrich_reactome_${meta.scenario}_${meta.method}.tsv"), emit: enrich

    script:
    """
    Rscript "${projectDir}/process/bio/enrich_pathway/templates/enrich_pathway.R" \
        --method '${meta.method}' \
        --filtered-genes '${filtered_genes}' \
        --universe-file '${universe_file}' \
        --output-file "bio_enrich_reactome_${meta.scenario}_${meta.method}.tsv"
    """

    stub:
    method = meta.method
    """
    echo -e "method\tDescription\tp.adjust\tGeneRatio\tCount" > bio_enrich_reactome_${meta.scenario}_${method}.tsv
    """
}
