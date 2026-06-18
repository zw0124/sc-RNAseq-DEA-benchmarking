process EDGER {
    tag "${meta.scenario}_${meta.n_cells ?: 'na'}_${meta.n_genes ?: 'na'}_${meta.run}_edger"

    input:
    tuple val(meta), path(sc_rds_file)
    val lfc_threshold

    output:
    tuple val(meta), path("${meta.scenario}.edger.tsv")

    script:
    meta = meta + [method: 'edger']
    """
        EDGER.R \
        --input ${sc_rds_file} \
        --output ${meta.scenario}.edger.tsv \
        --lfc_threshold ${lfc_threshold}
    """

    stub:
    meta = meta + [method: 'edger']
    """
    touch ${meta.scenario}.edger.tsv
    """
}
