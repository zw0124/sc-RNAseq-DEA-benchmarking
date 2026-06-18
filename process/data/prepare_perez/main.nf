process PREPARE_PEREZ {
    tag 'perez'
    input:
    path full_h5ad

    output:
    path 'perez_cd4_full_prepared.h5ad', emit: bio
    path 'perez_european_female.h5ad', emit: negative

    script:
    """
    prepare_perez.py \
        --input ${full_h5ad} \
        --bio-output perez_cd4_full_prepared.h5ad \
        --negative-output perez_european_female.h5ad
    """

    stub:
    """
    touch perez_cd4_full_prepared.h5ad
    touch perez_european_female.h5ad
    """
}
