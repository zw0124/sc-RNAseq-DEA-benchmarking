process PLOT_FPR {
    publishDir "${params.output}/neg/FPR", mode: 'copy'

    input:
    tuple val(scenario), path(fpr_files)

    output:
    tuple val(scenario), path("fpr_curve_${scenario}.png"), emit: fpr_curve

    script:
    output_file = "fpr_curve_${scenario}.png"
    plot_title = "Average FPR Curve - ${scenario}"
    template 'plot_fpr.R'
    stub:
    """
    touch fpr_curve_${scenario}.png
    """
}
