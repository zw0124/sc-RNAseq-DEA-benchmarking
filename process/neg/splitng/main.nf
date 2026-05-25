process SPLIT_FILES_NEGATIVE {
    tag "perez_negative_${sample_size}v${sample_size}_${run}"

    input:
    tuple path(input_h5ad), val(sample_size), val(run)

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    meta = [scenario: "perez_negative_${sample_size}v${sample_size}", run: run, sample_size: sample_size]
    template 'split_files.py'

    stub:
    meta = [scenario: "perez_negative_${sample_size}v${sample_size}", run: run, sample_size: sample_size]
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}
