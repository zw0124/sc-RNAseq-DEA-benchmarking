process SIMULATION {
    tag "${meta.scenario}_${meta.run}"

    input:
    val meta
    val n_genes

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad"), emit: h5ad
    tuple val(meta), path("design_${meta.scenario}_${meta.run}.tsv"), emit: design

    script:
    """
    simulation.R \
        --scenario ${meta.scenario} \
        --n_genes ${n_genes} \
        --output ${meta.scenario}_${meta.run}.h5ad \
        --design_output design_${meta.scenario}_${meta.run}.tsv \
        --seed ${meta.run}
    """

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    echo -e "scenario\trequested_seed\tused_seed\tSample\tBatch\tCondition\tsample_cells\tn_total\tn_case\tn_control\tbatch_n\tbatch_case\tbatch_control\ttarget_case_prob\ttarget_control_prob\ttarget_mean_cells\ttarget_cv_cells\tactual_mean_cells\tactual_cv_cells\ttotal_cells\tlfc" > design_${meta.scenario}_${meta.run}.tsv
    """
}
