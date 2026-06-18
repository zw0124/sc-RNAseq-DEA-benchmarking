process PLOT_RUNTIME_BENCHMARK {

    publishDir "${params.output}/runtime", mode: 'copy'
    cache false

    input:
    path fixed_cells
    path fixed_genes

    output:
    path 'runtime_benchmark.png'

    script:
    """
    Rscript "${projectDir}/process/run/plot_runtime_benchmark/templates/plot_runtime_benchmark.R" \
        --fixed-cells "${fixed_cells}" \
        --fixed-genes "${fixed_genes}" \
        --output runtime_benchmark.png
    """

    stub:
    """
    touch runtime_benchmark.png
    """
}
