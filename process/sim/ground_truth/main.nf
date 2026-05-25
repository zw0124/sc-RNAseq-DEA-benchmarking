process GROUND_TRUTH {
    tag "${meta.scenario}_${meta.run}"

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("de_${meta.scenario}_${meta.run}.tsv")

    script:
    template('ground_truth.py')

    stub:
    """
    touch de_${meta.scenario}_${meta.run}.tsv
    """
}
