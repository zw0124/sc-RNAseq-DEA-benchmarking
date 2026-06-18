include { SIMULATION } from '../process/run/simulation'
include { PREPROCESSING } from '../process/sim/preprocessing'
include { H5_SEURAT } from '../process/sim/h5_seurat'
include { DESEQ2 } from '../process/sim/DESEQ2'
include { LIMMA } from '../process/sim/LIMMA'
include { EDGER } from '../process/sim/EDGER'
include { WILCOX } from '../process/sim/WILCOX'
include { TTEST } from '../process/sim/TTEST'
include { MEMENTO } from '../process/sim/MEMENTO'
include { NEBULA } from '../process/sim/nebula/main'
include { BOOTSTRAP2 } from '../process/sim/BOOTSTRAP2'
include { BOOTSTRAP3 } from '../process/sim/BOOTSTRAP3'
include { PREPARE_RUNTIME_TRACE } from '../process/run/prepare_runtime_trace'
include { PLOT_RUNTIME_BENCHMARK } from '../process/run/plot_runtime_benchmark'


workflow RUNTIME {
    take:
        n_runs
        n_fixed_cells
        n_fixed_genes
        preprocessing_threshold
        lfc_threshold

    main:
        ch_fixed_cells = channel.from([scenario: 'fixed_cells', n_cells: n_fixed_cells])
            .combine(1..n_runs)
                .combine([(100..900).step(100), (1000..10000).step(500)].flatten())
                .map{meta, run, n_genes -> meta + [n_genes: n_genes, run: run]}

        ch_fixed_genes = channel.from([scenario: 'fixed_genes', n_genes: n_fixed_genes])
            .combine(1..n_runs)
                .combine([(100..900).step(100), (1000..10000).step(500)].flatten())
                .map{meta, run, n_cells -> meta + [n_cells: n_cells, run: run]}

        ch_sim_input = ch_fixed_cells.mix(ch_fixed_genes)
        

        SIMULATION(ch_sim_input)

        PREPROCESSING(SIMULATION.out, preprocessing_threshold)

        H5_SEURAT(PREPROCESSING.out)

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

        ch_runtime_gate = DESEQ2.out
            .mix(LIMMA.out)
            .mix(EDGER.out)
            .mix(WILCOX.out)
            .mix(TTEST.out)
            .mix(MEMENTO.out)
            .mix(NEBULA.out)
            .mix(BOOTSTRAP3.out)
            .mix(BOOTSTRAP2.out)

        PREPARE_RUNTIME_TRACE(ch_runtime_gate.collect())

        PLOT_RUNTIME_BENCHMARK(PREPARE_RUNTIME_TRACE.out.fixed_cells, PREPARE_RUNTIME_TRACE.out.fixed_genes)
}
