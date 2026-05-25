process MEMENTO {
    cpus 4
    tag "${meta.scenario}_${meta.n_cells ?: 'na'}_${meta.n_genes ?: 'na'}_${meta.run}_memento"

    input:
    tuple val(meta), path(sc_rds_file)

    output:
    tuple val(meta_memento), path("${meta.scenario}_memento_${meta.run}.tsv")

    script:
    meta_memento = meta + [method: 'memento']
    """
    runmemento.py \
        --input $sc_rds_file \
        --output ${meta.scenario}_memento_${meta.run}.tsv \
        --scenario ${meta.scenario} \
        --n_cores ${task.cpus}
    """

    stub:
    meta_memento = meta + [method: 'memento']
    """
    touch ${meta.scenario}_memento_${meta.run}.tsv
    """
}
