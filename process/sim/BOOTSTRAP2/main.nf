process  BOOTSTRAP2 {
    cpus 4
    tag "${meta.scenario}_${meta.n_cells ?: 'na'}_${meta.n_genes ?: 'na'}_${meta.run}_bootstrap2"

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')
    val lfc_threshold

    output:
    tuple val(meta), path("${meta.scenario}_BOOTSTRAP2_${meta.run}.tsv")

    script:
    meta = meta + [method: 'BOOTSTRAP2']
    """
    2hboot.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_BOOTSTRAP2_${meta.run}.tsv \
        --lfc_threshold ${lfc_threshold} \
        --n_jobs ${task.cpus}

    """

    stub:
    meta = meta + [method: 'BOOTSTRAP2']
    """
    touch ${meta.scenario}_BOOTSTRAP2_${meta.run}.tsv
    """
}
