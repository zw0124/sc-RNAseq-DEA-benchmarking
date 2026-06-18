process H5_SEURAT {
    tag "${meta.scenario}_${meta.run}"
    

    input:
    tuple val(meta), path(input_anndata)

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.seurat.rds"), emit: rds
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5_seurat.h5ad"), emit: h5ad

    script:
    template 'h5_seurat.R'

    stub:
    """
    touch ${meta.scenario}_${meta.run}.seurat.rds
    touch ${meta.scenario}_${meta.run}.h5_seurat.h5ad
    """
}
