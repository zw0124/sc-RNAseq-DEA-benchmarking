process DESEQ2 {
    tag "${meta.scenario}_${meta.n_cells ?: 'na'}_${meta.n_genes ?: 'na'}_${meta.run}_deseq2"

    input:
    tuple val(meta), path(sc_rds_file)
    val lfc_threshold

    output:
    tuple val(meta), path("${meta.scenario}.deseq2.tsv")

    script:
    meta = meta + [method: 'deseq2']
    """
        DESEQ2.R \
        --input ${sc_rds_file} \
        --output ${meta.scenario}.deseq2.tsv \
        --lfc_threshold ${lfc_threshold}
    """

    stub:
    meta = meta + [method: 'deseq2']
    """
    touch ${meta.scenario}.deseq2.tsv
    """
}
