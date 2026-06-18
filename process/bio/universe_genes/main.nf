process BIO_UNIVERSE_GENES {
    tag "${meta.scenario}"

    input:
    tuple val(meta), path(input_anndata)

    output:
    tuple val(meta), path("${meta.dataset}_${meta.scenario}.universe.tsv")

    script:
    template('extract_universe.py')

    stub:
    """
    echo -e "gene\nIL7R" > ${meta.dataset}_${meta.scenario}.universe.tsv
    """
}
