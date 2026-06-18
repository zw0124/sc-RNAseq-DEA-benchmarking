process LIMMA {
    tag "${meta.scenario}_${meta.n_cells ?: 'na'}_${meta.n_genes ?: 'na'}_${meta.run}_limma"

    input:
    tuple val(meta), path(sc_rds_file)
    val lfc_threshold

    output:
    tuple val(meta), path("${meta.scenario}.limma.tsv")

    script:
    meta = meta + [method: 'limma']
    """
        LIMMA.R \
        --input ${sc_rds_file} \
        --output ${meta.scenario}.limma.tsv \
        --lfc_threshold ${lfc_threshold}
    """

    stub:
    meta = meta + [method: 'limma']
    """
    touch ${meta.scenario}.limma.tsv
    """
}
