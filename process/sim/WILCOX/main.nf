process WILCOX {
    tag "${meta.scenario}_${meta.n_cells ?: 'na'}_${meta.n_genes ?: 'na'}_${meta.run}_wilcox"

    input:
    tuple val(meta), path(sc_rds_file)

    output:
    tuple val(meta), path("${meta.scenario}.wilcox.tsv")

    script:
    meta = meta + [method: 'wilcox']
    """
    WILCOX.R \
        --input ${sc_rds_file} \
        --output ${meta.scenario}.wilcox.tsv
    """

    stub:
    meta = meta + [method: 'wilcox']
    """
    touch ${meta.scenario}.wilcox.tsv
    """
}
