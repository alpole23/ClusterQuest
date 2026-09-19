#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

/**
 * Validate pipeline parameters.
 */
def validateParams() {
    def errors = []

    // Boolean params must be real booleans — see Utils.BOOLEAN_PARAMS for why.
    Utils.BOOLEAN_PARAMS.each { name ->
        if (params.containsKey(name)) {
            def value = params[name]
            if (!(value instanceof Boolean)) {
                def kind = value == null ? 'null' : value.getClass().getSimpleName()
                errors << ("params.${name} must be true or false, but is '${value}' (${kind}). " +
                           "A command-line `--${name} ${value}` arrives as a string, and every " +
                           "non-empty string is true — so it would ENABLE ${name} rather than " +
                           "set it. Pass booleans in a params file " +
                           "(-params-file params.json, \"${name}\": false), or omit the flag " +
                           "to keep the configured default.")
            }
        }
    }

    // Validate workflow
    def validWorkflows = ['download', 'bgc_analysis', 'full']
    if (!(params.workflow in validWorkflows)) {
        errors << "Invalid workflow '${params.workflow}'. Valid options: ${validWorkflows.join(', ')}"
    }

    // Validate clustering
    def validClustering = ['none', 'bigscape']
    if (!(params.clustering in validClustering)) {
        errors << "Invalid clustering '${params.clustering}'. Valid options: ${validClustering.join(', ')}"
    }

    // Validate bigscape alignment mode
    def validAlignmentModes = ['auto', 'global', 'glocal']
    if (!(params.bigscape_alignment_mode in validAlignmentModes)) {
        errors << "Invalid bigscape_alignment_mode '${params.bigscape_alignment_mode}'. Valid options: ${validAlignmentModes.join(', ')}"
    }

    // Validate bigscape classify scheme
    def validClassify = ['', 'category', 'class', 'legacy']
    if (!(params.bigscape_classify in validClassify)) {
        errors << "Invalid bigscape_classify '${params.bigscape_classify}'. Valid options: ${validClassify.collect { it ?: '\"\"' }.join(', ')}"
    }

    // Validate bgc_analysis workflow has required input
    if (params.workflow == 'bgc_analysis' && (!params.input_genomes || params.input_genomes == 'null')) {
        errors << "params.input_genomes must be specified for 'bgc_analysis' workflow"
    }

    // Validate reuse_antismash_from if specified
    if (params.reuse_antismash_from) {
        def reuse_dir = "${params.outdir}/antismash_results/${Utils.sanitizeTaxon(params.reuse_antismash_from)}"
        if (!file(reuse_dir).exists()) {
            errors << "reuse_antismash_from directory does not exist: ${reuse_dir}"
        }
    }

    // Validate reuse_gtdbtk_from if specified
    if (params.reuse_gtdbtk_from) {
        def reuse_dir = "${params.outdir}/gtdbtk_results/${Utils.sanitizeTaxon(params.reuse_gtdbtk_from)}/gtdbtk_output"
        if (!file(reuse_dir).exists()) {
            errors << "reuse_gtdbtk_from directory does not exist: ${reuse_dir}"
        }
    }

    // Report errors
    if (errors) {
        log.error "=" * 60
        log.error "PARAMETER VALIDATION FAILED"
        log.error "=" * 60
        errors.each { log.error "  - ${it}" }
        log.error "=" * 60
        error "Please fix the above parameter errors and try again."
    }

    // Warnings
    if (params.run_gtdbtk) {
        log.warn "GTDB-Tk is enabled. This requires ~140 GB disk space and ~56-64 GB RAM."
    }

    if (params.reuse_antismash_from) {
        log.info "antiSMASH result reuse enabled from taxon: ${params.reuse_antismash_from}"
    }

    if (params.reuse_gtdbtk_from) {
        log.info "GTDB-Tk result reuse enabled from taxon: ${params.reuse_gtdbtk_from}"
    }

    // Say which reference databases this run classifies against. storeDir means
    // a database is downloaded once and then never re-checked, so without this
    // the version is invisible — and silently whatever was current on the day of
    // the first run. Never fails the run; upgrading is a deliberate edit.
    log.info "Reference databases:"
    DbVersions.report(params).each { log.info it }
}

// =============================================================================
// SUBWORKFLOWS
// =============================================================================

include { DOWNLOAD_GENOMES } from './subworkflows/download_genomes'
include { BGC_ANALYSIS } from './subworkflows/bgc_analysis'

// =============================================================================
// ENTRY POINT
// =============================================================================

/*
 * Main entry point
 */
workflow {
    validateParams()
    def taxon_dir = Utils.sanitizeTaxon(params.taxon)

    if (params.workflow == "download") {
        DOWNLOAD_GENOMES(params.taxon)

    } else if (params.workflow == "bgc_analysis") {
        // Load from previous download run
        def base_dir = "${params.outdir}/ncbi_genomes/${taxon_dir}"
        def results_dir = "${params.outdir}/main_analysis_results/${taxon_dir}"

        BGC_ANALYSIS(
            params.taxon,
            Channel.fromPath("${params.input_genomes}/*.gbff"),
            Channel.fromPath("${base_dir}/ncbi_dataset/data/assembly_info_table.txt"),
            Channel.fromPath("${base_dir}/name_map.json"),
            Channel.fromPath("${results_dir}/taxonomy_map.json")
        )

    } else if (params.workflow == "full") {
        DOWNLOAD_GENOMES(params.taxon)
        BGC_ANALYSIS(
            params.taxon,
            DOWNLOAD_GENOMES.out.renamed_genomes,
            DOWNLOAD_GENOMES.out.assembly_info,
            DOWNLOAD_GENOMES.out.name_map,
            DOWNLOAD_GENOMES.out.taxonomy_map
        )

    } else {
        error "ERROR: Invalid workflow '${params.workflow}'. Valid options: 'download', 'bgc_analysis', 'full'"
    }
}
