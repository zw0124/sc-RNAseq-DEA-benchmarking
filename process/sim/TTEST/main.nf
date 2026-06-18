process TTEST {
    tag "${meta.scenario}_${meta.n_cells ?: 'na'}_${meta.n_genes ?: 'na'}_${meta.run}_ttest"

    input:
    tuple val(meta), path(sc_rds_file)

    output:
    tuple val(meta), path("${meta.scenario}.ttest.tsv")

    script:
    meta = meta + [method: 'ttest']
    """
    TTEST.R \
        --input ${sc_rds_file} \
        --output ${meta.scenario}.ttest.tsv
    """

    stub:
    meta = meta + [method: 'ttest']
    """
    touch ${meta.scenario}.ttest.tsv
    """
}
