process BOOTSTRAP3 {
    tag "${meta.scenario}_${meta.n_cells ?: 'na'}_${meta.n_genes ?: 'na'}_${meta.run}_bootstrap3"

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_BOOTSTRAP3_${meta.run}.tsv")

    script:
    meta = meta + [method: 'BOOTSTRAP3']
    """
    3hboot.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_BOOTSTRAP3_${meta.run}.tsv
    """

    stub:
    meta = meta + [method: 'BOOTSTRAP3']
    """
    echo -e "gene\tp_val\tlfc" > ${meta.scenario}_BOOTSTRAP3_${meta.run}.tsv
    """
}
