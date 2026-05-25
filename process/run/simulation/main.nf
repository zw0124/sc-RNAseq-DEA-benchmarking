process SIMULATION {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    

    input:
    val meta

    output:
    tuple val(meta), path("simulation.h5ad")

    script:
    """
    sim-runtime.R \
        simulation.h5ad \
        ${meta.n_genes} \
        ${meta.n_cells} \
        ${meta.run}
    """

    stub:
    """
    touch simulation.h5ad
    """
}
