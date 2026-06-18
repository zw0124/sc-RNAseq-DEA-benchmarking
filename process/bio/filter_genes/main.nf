process BIO_FILTER_GENES {
    tag "${meta.scenario}"

    input:
    tuple val(meta), path(de_result)
    val lfc_threshold
    val adj_p_cutoff

    output:
    tuple val(meta), path("filtered_genes_${meta.dataset}_${meta.scenario}_${meta.method}.tsv"), emit: filtered

    script:
    """
    Rscript "${projectDir}/process/bio/filter_genes/templates/filter_genes.R" \
        --method '${meta.method}' \
        --de-result '${de_result}' \
        --lfc-threshold '${lfc_threshold}' \
        --adj-p-cutoff '${adj_p_cutoff}' \
        --output-file "filtered_genes_${meta.dataset}_${meta.scenario}_${meta.method}.tsv"
    """

    stub:
    """
    touch "filtered_genes_${meta.dataset}_${meta.scenario}_${meta.method}.tsv"
    """
}
