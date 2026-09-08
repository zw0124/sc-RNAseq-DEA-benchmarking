process PLOT_REVIEW_SIMULATION {
    tag "reviewer_simulation_plots"
    publishDir "${params.output}/reviewer-analysis/simulation", mode: 'copy'

    input:
    path "metrics/*"
    path "design/*"
    path "fpr/*"

    output:
    path "reviewer_figure1.png"
    path "reviewer_figure1_summary.tsv"
    path "review_sample_size_fdr0p01.png"
    path "review_sample_size_fdr0p05.png"
    path "review_sample_size_fdr0p10.png"
    path "review_cell_imbalance.png"
    path "review_cells_per_subject.png"
    path "review_combined_imbalance.png"
    path "review_simulation_runs.tsv"
    path "review_simulation_summary.tsv"
    path "review_simulation_design.tsv"
    path "review_simulation_fpr_runs.tsv"
    path "review_simulation_fpr_summary.tsv"
    path "review_sample_size_runs.tsv"
    path "review_cell_imbalance_runs.tsv"
    path "review_cells_per_subject_runs.tsv"
    path "review_combined_imbalance_runs.tsv"

    script:
    """
    Rscript "${projectDir}/process/sim/plot_review_simulation/templates/plot_review_simulation.R" \
        --metrics-dir metrics \
        --design-dir design \
        --fpr-dir fpr
    """

    stub:
    """
    touch reviewer_figure1.png
    echo -e "method\tlfc\tsetting_value\tthreshold\tobserved_fdr\tpower\tobserved_fdr_sd\tpower_sd\tn_runs\tlfc_label\tn_label" > reviewer_figure1_summary.tsv
    touch review_sample_size_fdr0p01.png
    touch review_sample_size_fdr0p05.png
    touch review_sample_size_fdr0p10.png
    touch review_cell_imbalance.png
    touch review_cells_per_subject.png
    touch review_combined_imbalance.png
    echo -e "scenario\trun\tmethod\tthreshold\tobserved_fdr\tpower" > review_simulation_runs.tsv
    echo -e "family\tsetting\tlfc\tthreshold\tmethod\tmetric\tmean\tsd\tn_runs" > review_simulation_summary.tsv
    echo -e "scenario\trequested_seed\tused_seed\tSample\tBatch\tCondition\tsample_cells\tn_total\tn_case\tn_control\tbatch_n\tbatch_case\tbatch_control\ttarget_case_prob\ttarget_control_prob\ttarget_mean_cells\ttarget_cv_cells\tactual_mean_cells\tactual_cv_cells\ttotal_cells\tlfc" > review_simulation_design.tsv
    echo -e "scenario\trun\tmethod\traw_p_value\tfpr\tfamily\tsetting\tsetting_value\tlfc" > review_simulation_fpr_runs.tsv
    echo -e "family\tsetting\tsetting_value\tlfc\traw_p_value\tmethod\tmean_fpr\tsd_fpr\tn_runs" > review_simulation_fpr_summary.tsv
    echo -e "scenario\trun\tmethod\tthreshold\tobserved_fdr\tpower\tdiscoveries\tfamily\tsetting\tsetting_value\tlfc" > review_sample_size_runs.tsv
    echo -e "scenario\trun\tmethod\tthreshold\tobserved_fdr\tpower\tdiscoveries\tfamily\tsetting\tsetting_value\tlfc" > review_cell_imbalance_runs.tsv
    echo -e "scenario\trun\tmethod\tthreshold\tobserved_fdr\tpower\tdiscoveries\tfamily\tsetting\tsetting_value\tlfc" > review_cells_per_subject_runs.tsv
    echo -e "scenario\trun\tmethod\tthreshold\tobserved_fdr\tpower\tdiscoveries\tfamily\tsetting\tsetting_value\tlfc" > review_combined_imbalance_runs.tsv
    """
}
