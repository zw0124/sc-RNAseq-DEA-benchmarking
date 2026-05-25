process NEBULA {
    cpus 4
    tag "${meta.scenario}_${meta.n_cells ?: 'na'}_${meta.n_genes ?: 'na'}_${meta.run}_nebula"

    input:
    tuple val(meta), path(rds_file)

    output:
    tuple val(meta_nebula), path("nebula_${meta.scenario}.tsv")

    script:
    meta_nebula = meta + [method: 'nebula']
    """
    nebula.R \
        --input ${rds_file} \
        --scenario ${meta.scenario} \
        --output nebula_${meta.scenario}.tsv \
        --n_cores ${task.cpus}
    """

    stub:
    meta_nebula = meta + [method: 'nebula']
    """
    touch nebula_${meta.scenario}.tsv
    """
}
