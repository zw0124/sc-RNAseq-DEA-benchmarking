process PLOT_FPR_BIN {
    publishDir "${params.output}/neg/fpr_bin", mode: 'copy'

    input:
    tuple val(scenario), path(fpr_files)
    val p_val

    output:
    tuple val(scenario), path("fpr_boxplot_${scenario}_p${p_val}.png"), emit: fpr_bin

    script:
    output_file = "fpr_boxplot_${scenario}_p${p_val}.png"
    """
    Rscript "${projectDir}/process/sim/plot_fpr_bin/templates/plot_fpr_bin.R" \
        --input-dir . \
        --p-val "${p_val}" \
        --output "${output_file}"
    """

    stub:
    """
    touch fpr_boxplot_${scenario}_p${p_val}.png
    """
}
