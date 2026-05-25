process CALCULATE_METRICS {
    tag "${meta.scenario}_${meta.run}"

    input:
    tuple val(meta), path(pvalues, stageAs: 'pvalues.tsv'), path(ground_truth, stageAs: 'ground_truth.tsv')

    output:
    tuple val(meta), path("prc_${meta.scenario}_${meta.run}_*.tsv"), emit: 'prc'
    tuple val(meta), path("tpr_${meta.scenario}_${meta.run}_*.tsv"), emit: 'tpr'
    tuple val(meta), path("fpr_${meta.scenario}_${meta.run}_*.tsv"), emit: 'fpr'
    tuple val(meta), path("fdr_power_${meta.scenario}_${meta.run}_*.tsv"), emit: 'fdr_power'
    tuple val(meta), path("auc_${meta.scenario}_${meta.run}.tsv"), emit: 'auc'

    script:
    template 'calculate_metrics.py'

    stub:
    """
    touch prc_${meta.scenario}_${meta.run}_deseq2.tsv
    touch tpr_${meta.scenario}_${meta.run}_deseq2.tsv
    touch fpr_${meta.scenario}_${meta.run}_deseq2.tsv
    touch fdr_power_${meta.scenario}_${meta.run}_deseq2.tsv
    touch prc_${meta.scenario}_${meta.run}_ttest.tsv
    touch tpr_${meta.scenario}_${meta.run}_ttest.tsv
    touch fpr_${meta.scenario}_${meta.run}_ttest.tsv
    touch fdr_power_${meta.scenario}_${meta.run}_ttest.tsv

    touch auc_${meta.scenario}_${meta.run}.tsv
    """
}
