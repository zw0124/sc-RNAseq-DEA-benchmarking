process PREPARE_RUNTIME_TRACE {
    cache false

    input:
    val runtime_gate

    output:
    path 'fixed_cells.tsv', emit: fixed_cells
    path 'fixed_genes.tsv', emit: fixed_genes

    script:
    """
    python3 "${projectDir}/process/run/prepare_runtime_trace/templates/prepare_runtime_trace.py" \
        --trace-file "${params.output}/trace.txt" \
        --fixed-cells fixed_cells.tsv \
        --fixed-genes fixed_genes.tsv
    """

    stub:
    """
    touch fixed_cells.tsv
    touch fixed_genes.tsv
    """
}
