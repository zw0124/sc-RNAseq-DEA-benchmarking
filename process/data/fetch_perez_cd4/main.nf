process FETCH_PEREZ_CD4 {
    tag 'perez_cd4'
    output:
    path 'perez_cd4_full.h5ad'

    script:
    """
    fetch_perez_cd4.py \
        --output perez_cd4_full.h5ad \
        --dataset-id '${params.cellxgene_dataset_id}'
    """

    stub:
    """
    touch perez_cd4_full.h5ad
    """
}
