include { NEGATIVE } from './workflow/neg.nf'
include { RUNTIME } from './workflow/run.nf'
include { PERFORMANCE } from './workflow/sim.nf'
include { BIO } from './workflow/bio.nf'
include { FETCH_PEREZ_CD4 } from './process/data/fetch_perez_cd4/main'
include { PREPARE_PEREZ } from './process/data/prepare_perez/main'

workflow {
    preprocessing_threshold = params.preprocessing_threshold
    lfc_threshold = params.lfc_threshold
    lfc_threshold_lfc05 = params.lfc_threshold_lfc05

    negative_enabled = params.negative_enabled
    negative_n_runs = params.negative_n_runs
    negative_p_val = params.negative_p_val
    runtime_enabled = params.runtime_enabled
    runtime_n_runs = params.runtime_n_runs
    runtime_n_fixed_cells = params.runtime_n_fixed_cells
    runtime_n_fixed_genes = params.runtime_n_fixed_genes

    performance_enabled = params.performance_enabled
    performance_n_runs = params.performance_n_runs
    performance_n_genes = params.performance_n_genes

    bio_enabled = params.bio_enabled
    bio_adj_p_cutoff = params.bio_adj_p_cutoff
    real_data_enabled = negative_enabled || bio_enabled
    if (real_data_enabled) {
        FETCH_PEREZ_CD4()
        PREPARE_PEREZ(FETCH_PEREZ_CD4.out)
    }

    if (performance_enabled) {
        PERFORMANCE(
            performance_n_runs,
            performance_n_genes,
            preprocessing_threshold,
            lfc_threshold
        )
    }

    if (negative_enabled) {
        NEGATIVE(
            PREPARE_PEREZ.out.negative,
            negative_n_runs,
            preprocessing_threshold,
            lfc_threshold,
            negative_p_val
        )
    }
    if (runtime_enabled) {
        RUNTIME(
            runtime_n_runs,
            runtime_n_fixed_cells,
            runtime_n_fixed_genes,
            preprocessing_threshold,
            lfc_threshold
        )
    }

    if (bio_enabled) {
        BIO(
            PREPARE_PEREZ.out.bio,
            preprocessing_threshold,
            lfc_threshold,
            lfc_threshold_lfc05,
            bio_adj_p_cutoff
        )
    }
}
