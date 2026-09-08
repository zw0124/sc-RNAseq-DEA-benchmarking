process PLOT_MEMORY_BENCHMARK {
    publishDir "${params.output}/runtime", mode: 'copy'
    cache false

    input:
    path memory_summary

    output:
    path 'memory_benchmark.png'

    script:
    """
    Rscript "${projectDir}/process/run/plot_memory_benchmark/templates/plot_memory_benchmark.R" \
        --input "${memory_summary}" \
        --output memory_benchmark.png
    """

    stub:
    """
    touch memory_benchmark.png
    """
}
