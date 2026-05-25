include { PREPROCESSING } from '../process/sim/preprocessing'
include { H5_SEURAT } from '../process/sim/h5_seurat'
include { DESEQ2 as DESEQ2_TEST0; DESEQ2 as DESEQ2_TEST05 } from '../process/sim/DESEQ2'
include { LIMMA as LIMMA_TEST0; LIMMA as LIMMA_TEST05 } from '../process/sim/LIMMA'
include { EDGER as EDGER_TEST0; EDGER as EDGER_TEST05 } from '../process/sim/EDGER'
include { WILCOX as WILCOX_TEST0 } from '../process/sim/WILCOX'
include { TTEST as TTEST_TEST0 } from '../process/sim/TTEST'
include { MEMENTO as MEMENTO_TEST0 } from '../process/sim/MEMENTO'
include { NEBULA as NEBULA_TEST0 } from '../process/sim/nebula/main'
include { BOOTSTRAP2 as BOOTSTRAP2_TEST0; BOOTSTRAP2 as BOOTSTRAP2_TEST05 } from '../process/sim/BOOTSTRAP2'
include { BOOTSTRAP3 as BOOTSTRAP3_TEST0 } from '../process/sim/BOOTSTRAP3'
include { BIO_UNIVERSE_GENES } from '../process/bio/universe_genes/main'
include { BIO_FILTER_GENES as FILTER_LFC0_FILTER0; BIO_FILTER_GENES as FILTER_LFC0_FILTER05; BIO_FILTER_GENES as FILTER_LFC05_FILTER05 } from '../process/bio/filter_genes/main'
include { BIO_ENRICH_PATHWAY as REACTOME_LFC0_FILTER0; BIO_ENRICH_PATHWAY as REACTOME_LFC0_FILTER05; BIO_ENRICH_PATHWAY as REACTOME_LFC05_FILTER05 } from '../process/bio/enrich_pathway/main'
include { BIO_ENRICHGO_COMBINED_PLOT as REACTOME_COMBINED_PLOT } from '../process/bio/enrichgo_combined_plot/main'

workflow BIO {
    take:
    prepared_bio_h5ad
    preprocessing_threshold
    lfc0_threshold
    lfc05_threshold
    bio_adj_p_cutoff

    main:

    ch_bio_renamed = prepared_bio_h5ad.map { file -> tuple([scenario: 'bio'], file) }
    PREPROCESSING(ch_bio_renamed, preprocessing_threshold)
    BIO_UNIVERSE_GENES(PREPROCESSING.out)

    H5_SEURAT(PREPROCESSING.out)

    ch_sc_rds = H5_SEURAT.out.rds
    ch_sc_h5ad = H5_SEURAT.out.h5ad

    ch_sc_rds_lfc0 = ch_sc_rds.map { meta, file -> tuple(meta + [scenario: 'bio_lfc0'], file) }
    ch_sc_h5ad_lfc0 = ch_sc_h5ad.map { meta, file -> tuple(meta + [scenario: 'bio_lfc0'], file) }
    ch_sc_rds_lfc05 = ch_sc_rds.map { meta, file -> tuple(meta + [scenario: 'bio_lfc0p5'], file) }
    ch_sc_h5ad_lfc05 = ch_sc_h5ad.map { meta, file -> tuple(meta + [scenario: 'bio_lfc0p5'], file) }

    DESEQ2_TEST0(ch_sc_rds_lfc0, lfc0_threshold)
    LIMMA_TEST0(ch_sc_rds_lfc0, lfc0_threshold)
    EDGER_TEST0(ch_sc_rds_lfc0, lfc0_threshold)
    WILCOX_TEST0(ch_sc_rds_lfc0)
    TTEST_TEST0(ch_sc_rds_lfc0)
    NEBULA_TEST0(ch_sc_rds_lfc0)
    MEMENTO_TEST0(ch_sc_h5ad_lfc0)
    BOOTSTRAP3_TEST0(ch_sc_h5ad_lfc0)
    BOOTSTRAP2_TEST0(ch_sc_h5ad_lfc0, lfc0_threshold)

    DESEQ2_TEST05(ch_sc_rds_lfc05, lfc05_threshold)
    LIMMA_TEST05(ch_sc_rds_lfc05, lfc05_threshold)
    EDGER_TEST05(ch_sc_rds_lfc05, lfc05_threshold)
    BOOTSTRAP2_TEST05(ch_sc_h5ad_lfc05, lfc05_threshold)

    bio_methods_lfc0 = DESEQ2_TEST0.out
        .mix(LIMMA_TEST0.out)
        .mix(EDGER_TEST0.out)
        .mix(WILCOX_TEST0.out)
        .mix(TTEST_TEST0.out)
        .mix(MEMENTO_TEST0.out)
        .mix(NEBULA_TEST0.out)
        .mix(BOOTSTRAP3_TEST0.out)
        .mix(BOOTSTRAP2_TEST0.out)

    bio_methods_lfc05 = DESEQ2_TEST05.out
        .mix(LIMMA_TEST05.out)
        .mix(EDGER_TEST05.out)
        .mix(BOOTSTRAP2_TEST05.out)

    ch_lfc0_filter0 = bio_methods_lfc0.map { meta, file ->
        tuple(meta + [scenario: 'lfc0_filter0'], file)
    }
    ch_lfc0_filter05 = bio_methods_lfc0.map { meta, file ->
        tuple(meta + [scenario: 'lfc0_filter0p5'], file)
    }
    ch_lfc05_filter05 = bio_methods_lfc05.map { meta, file ->
        tuple(meta + [scenario: 'lfc0p5_filter0p5'], file)
    }

    FILTER_LFC0_FILTER0(ch_lfc0_filter0, lfc0_threshold, bio_adj_p_cutoff)
    FILTER_LFC0_FILTER05(ch_lfc0_filter05, lfc05_threshold, bio_adj_p_cutoff)
    FILTER_LFC05_FILTER05(ch_lfc05_filter05, lfc05_threshold, bio_adj_p_cutoff)

    ch_enrich_lfc0_filter0 = FILTER_LFC0_FILTER0.out.filtered
        .combine(BIO_UNIVERSE_GENES.out)
        .map { meta1, filtered, _meta2, universe -> tuple(meta1, filtered, universe) }
    ch_enrich_lfc0_filter05 = FILTER_LFC0_FILTER05.out.filtered
        .combine(BIO_UNIVERSE_GENES.out)
        .map { meta1, filtered, _meta2, universe -> tuple(meta1, filtered, universe) }
    ch_enrich_lfc05_filter05 = FILTER_LFC05_FILTER05.out.filtered
        .combine(BIO_UNIVERSE_GENES.out)
        .map { meta1, filtered, _meta2, universe -> tuple(meta1, filtered, universe) }

    REACTOME_LFC0_FILTER0(ch_enrich_lfc0_filter0)
    REACTOME_LFC0_FILTER05(ch_enrich_lfc0_filter05)
    REACTOME_LFC05_FILTER05(ch_enrich_lfc05_filter05)

    ch_plot_reactome_lfc0_filter0 = REACTOME_LFC0_FILTER0.out.enrich
        .map { _meta, enrich_file -> enrich_file }
        .collect()
    ch_plot_reactome_lfc0_filter05 = REACTOME_LFC0_FILTER05.out.enrich
        .map { _meta, enrich_file -> enrich_file }
        .collect()
    ch_plot_reactome_lfc05_filter05 = REACTOME_LFC05_FILTER05.out.enrich
        .map { _meta, enrich_file -> enrich_file }
        .collect()

    REACTOME_COMBINED_PLOT(
        ch_plot_reactome_lfc0_filter0,
        ch_plot_reactome_lfc0_filter05,
        ch_plot_reactome_lfc05_filter05
    )
}
