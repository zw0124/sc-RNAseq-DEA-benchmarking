include { SPLIT_FILES_NEGATIVE } from '../process/neg/splitng'
include { GROUND_TRUTH } from '../process/sim/ground_truth'
include { PREPROCESSING } from '../process/sim/preprocessing'
include { PVALUES } from '../process/sim/pvalues/main'
include { CALCULATE_METRICS } from '../process/sim/calculate_metrics/main'
include { DESEQ2 } from '../process/sim/DESEQ2'
include { LIMMA } from '../process/sim/LIMMA'
include { EDGER } from '../process/sim/EDGER'
include { WILCOX } from '../process/sim/WILCOX'
include { TTEST } from '../process/sim/TTEST'
include { MEMENTO } from '../process/sim/MEMENTO'
include { NEBULA } from '../process/sim/nebula/main'
include { PLOT_FPR } from '../process/sim/plot_fpr/main'
include { PLOT_FPR_BIN } from '../process/sim/plot_fpr_bin/main'
include { BOOTSTRAP2 } from '../process/sim/BOOTSTRAP2'
include { BOOTSTRAP3 } from '../process/sim/BOOTSTRAP3'
include { H5_SEURAT } from '../process/sim/h5_seurat'

workflow NEGATIVE{
    take:
        negative_input_h5ad
        negative_n_runs
        preprocessing_threshold
        lfc_threshold
        p_val

    main:
        ch_lfc_threshold = channel.value(lfc_threshold)

        ch_runs = channel.from(1..negative_n_runs)
        ch_sample_sizes = channel.from(5, 10, 15)
        ch_negative_design = ch_sample_sizes.combine(ch_runs)
        ch_negative_input_design = negative_input_h5ad
            .combine(ch_negative_design)
            .map { input_h5ad, sample_size, run -> tuple(input_h5ad, sample_size, run) }

        SPLIT_FILES_NEGATIVE(ch_negative_input_design)

        GROUND_TRUTH(SPLIT_FILES_NEGATIVE.out)
        PREPROCESSING(SPLIT_FILES_NEGATIVE.out, preprocessing_threshold)
        ch_methods_input = PREPROCESSING.out

    H5_SEURAT(ch_methods_input)

    ch_sc_rds = H5_SEURAT.out.rds
    ch_sc_h5ad = H5_SEURAT.out.h5ad

    DESEQ2(ch_sc_rds, ch_lfc_threshold)
    LIMMA(ch_sc_rds, ch_lfc_threshold)
    EDGER(ch_sc_rds, ch_lfc_threshold)

    WILCOX(ch_sc_rds)
    TTEST(ch_sc_rds)

    NEBULA(ch_sc_rds)
    

    MEMENTO(ch_sc_h5ad)
    BOOTSTRAP3(ch_sc_h5ad)
    BOOTSTRAP2(ch_sc_h5ad, ch_lfc_threshold)

    sc_all_unprocessed = DESEQ2.out
        .mix(LIMMA.out)
        .mix(EDGER.out)
        .mix(WILCOX.out)
        .mix(TTEST.out)
        .mix(MEMENTO.out)
        .mix(NEBULA.out)
        .mix(BOOTSTRAP3.out)
        .mix(BOOTSTRAP2.out)

    ch_pvalues_input = sc_all_unprocessed
        .map { meta, path -> [meta.scenario + '_' + meta.run, meta, path] }
        .groupTuple()
        .map { _key, meta_list, path_list -> [meta_list, path_list] }

    PVALUES(ch_pvalues_input)


    ch_meta_pval_groundtruth = (PVALUES.out
        .combine(GROUND_TRUTH.out))
        .filter { pval_meta, _pval_path, truth_meta, _truth_path ->
            pval_meta.scenario == truth_meta.scenario && pval_meta.run == truth_meta.run
        }
        .map { pval_meta, pval_path, _truth_meta, truth_path -> tuple(pval_meta, pval_path, truth_path) }

    CALCULATE_METRICS(ch_meta_pval_groundtruth)

    ch_fpr_all_sample_sizes = CALCULATE_METRICS.out.fpr
        .map { meta, path_list -> tuple('perez_negative_sample_sizes', path_list) }
        .groupTuple()
        .map { scenario, path_lists -> tuple(scenario, path_lists.flatten()) }


    PLOT_FPR(ch_fpr_all_sample_sizes)
    PLOT_FPR_BIN(ch_fpr_all_sample_sizes, p_val)
}
