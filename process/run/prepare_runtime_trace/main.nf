process PREPARE_RUNTIME_TRACE {
    cache false
    publishDir "${params.output}/runtime", mode: 'copy', pattern: 'memory_*.tsv'

    input:
    val runtime_gate

    output:
    path 'fixed_cells.tsv', emit: fixed_cells
    path 'fixed_genes.tsv', emit: fixed_genes
    path 'memory_runs.tsv', emit: memory_runs
    path 'memory_summary.tsv', emit: memory_summary

    script:
    """
    python3 "${projectDir}/process/run/prepare_runtime_trace/templates/prepare_runtime_trace.py" \
        --trace-file "${params.output}/trace.txt" \
        --fixed-cells fixed_cells.tsv \
        --fixed-genes fixed_genes.tsv \
        --memory-runs memory_runs.tsv \
        --memory-summary memory_summary.tsv
    """

    stub:
    """
    touch fixed_cells.tsv
    touch fixed_genes.tsv
    echo -e "scenario\tn_cells\tn_genes\trun\tmethod\tpeak_rss_gib" > memory_runs.tsv
    echo -e "scenario\tn_cells\tn_genes\tmethod\tmean_peak_rss_gib\tsd_peak_rss_gib\tn_runs" > memory_summary.tsv
    """
}
