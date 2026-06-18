process PLOT_OBS_FDR_POWER_GRID {

    publishDir "${params.output}/Fig_ObsFDR_Power_Grid", mode: 'copy'

    input:
    path "input_dir/*"

    output:
    path "Fig_ObsFDR_Power_Grid_*.png"

    script:
    """
    Rscript "${projectDir}/process/sim/plot_obs_fdr_power_grid/templates/plot_obs_fdr_power_grid.R" --input_dir input_dir
    """

    stub:
    """
    touch Fig_ObsFDR_Power_Grid_stub.png
    """
}
