include { SIMULATION } from '../process/sim/simulation'
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
include { PLOT_OBS_FDR_POWER_GRID } from '../process/sim/plot_obs_fdr_power_grid/main'
include { BOOTSTRAP2 } from '../process/sim/BOOTSTRAP2'
include { BOOTSTRAP3 } from '../process/sim/BOOTSTRAP3'
include { H5_SEURAT } from '../process/sim/h5_seurat'

workflow PERFORMANCE {
    take:
    n_runs
    n_genes
    preprocessing_threshold
    lfc_threshold

    main:
    scenarios_regular = [
        'dataset_n10_lfc1',
        'dataset_n10_lfc1p5',
        'dataset_n10_lfc2',
        'dataset_n20_lfc1',
        'dataset_n20_lfc1p5',
        'dataset_n20_lfc2',
        'dataset_n30_lfc1',
        'dataset_n30_lfc1p5',
        'dataset_n30_lfc2',
    ]

    ch_scenarios_regular = channel.from(scenarios_regular)
    ch_runs = channel.from(1..n_runs)

    ch_meta_regular = ch_scenarios_regular.combine(ch_runs).map { scenario, run -> [scenario: scenario, run: run] }
    SIMULATION(ch_meta_regular, n_genes)
    ch_sim_all = SIMULATION.out

    GROUND_TRUTH(ch_sim_all)

    PREPROCESSING(ch_sim_all, preprocessing_threshold)

    ch_methods_input = PREPROCESSING.out

    H5_SEURAT(ch_methods_input)

    ch_sc_rds = H5_SEURAT.out.rds
    ch_sc_h5ad = H5_SEURAT.out.h5ad

    DESEQ2(ch_sc_rds, lfc_threshold)
    LIMMA(ch_sc_rds, lfc_threshold)
    EDGER(ch_sc_rds, lfc_threshold)

    WILCOX(ch_sc_rds)
    TTEST(ch_sc_rds)

    NEBULA(ch_sc_rds)

    MEMENTO(ch_sc_h5ad)
    BOOTSTRAP3(ch_sc_h5ad)
    BOOTSTRAP2(ch_sc_h5ad, lfc_threshold)

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


    ch_obs_fdr_power_grid = CALCULATE_METRICS.out.fdr_power
        .map { _meta, path_list -> path_list }
        .flatten()
        .collect()
    PLOT_OBS_FDR_POWER_GRID(ch_obs_fdr_power_grid)
}
