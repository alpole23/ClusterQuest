#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

/**
 * Validate pipeline parameters.
 */
def validateParams() {
    // Declared params, declared types, and no unknown names — nf-schema, nf-core style.
    // This is what stops `--taxn Pantoea` running the default taxon in silence, and
    // `--run_gtdbtk false` enabling what it appears to disable.
    validateParameters()

    def errors = []

    // Types and unknown names are the schema's job (nextflow_schema.json, checked by
    // validateParameters() below); what follows are the cross-param rules a schema
    // cannot express.
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

    // Validate assembly source. Unvalidated until 2026-11-02, so a typo reached the
    // `datasets` CLI inside a task and failed there -- after metadata fetch had started
    // and, at order scale, after a long queue wait. The CLI's own check is
    // case-INSENSITIVE ('genbank' and 'GenBank' both work), so this is too, or it would
    // reject values that work today.
    //
    // `all` is REJECTED although the CLI accepts it. It returns both GenBank and RefSeq,
    // which is most genomes twice -- a GCA_ and a GCF_ accession for one assembly.
    // Measured on Brenneria 2026-11-02: GenBank 60, RefSeq 58, all 119, of which 59 are
    // the same assembly under both accessions. Both twins survive renaming, because the
    // accession is appended when organism names collide, so both are screened, both
    // reach antiSMASH and both enter BiG-SCAPE -- inflating family member counts in a
    // way that reads as biological rather than as double counting. There is no question
    // this pipeline asks that `all` answers better than one source alone.
    def validAssemblySources = ['genbank', 'refseq']
    if (!(params.assembly_source?.toString()?.toLowerCase() in validAssemblySources)) {
        def extra = params.assembly_source?.toString()?.toLowerCase() == 'all'
            ? " 'all' returns each assembly under both accessions, so it double-counts" +
              " genomes; pick one source."
            : ''
        errors << "Invalid assembly_source '${params.assembly_source}'. Valid options: " +
                  "GenBank, RefSeq (case-insensitive).${extra}"
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

    // `publish_mode = 'link'` hard-links published output into outdir, and a hard
    // link cannot cross a filesystem. Point outdir at a second drive, leave the work
    // directory at Nextflow's default ./work, and EVERY publish fails with
    //
    //     Invalid cross-device link
    //
    // after the task has already run -- so the compute is paid for and then thrown
    // away. Hit on 2026-11-02 with outdir on /mnt/bigdisk and work on /. Checked here,
    // where it costs nothing and names both remedies, rather than 30 s in.
    //
    // null means undeterminable (an unreadable parent, an exotic filesystem); only a
    // definite false blocks the run.
    if (params.publish_mode == 'link') {
        if (Utils.sameFileStore(params.outdir, workDir) == false) {
            errors << ("publish_mode 'link' needs outdir and the work directory on ONE " +
                       "filesystem -- a hard link cannot cross devices.\n" +
                       "        outdir: ${params.outdir}\n" +
                       "        work:   ${workDir}\n" +
                       "      Either move the work directory beside outdir, keeping the " +
                       "hard links and their storage saving:\n" +
                       "        -w ${file(params.outdir).toAbsolutePath().parent}/work\n" +
                       "      or publish by copying, which costs a second copy of every " +
                       "published file:\n" +
                       "        add \"publish_mode\": \"copy\" to the params file")
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

    // Run with nothing on the command line and every parameter is a default, which is
    // legitimate -- and is also what someone who meant to configure the run looks like.
    // Say what is about to happen and where the questions are. A notice, not a warning:
    // the defaults are a real configuration, and CI runs this way on purpose.
    //
    // The prompts CANNOT live here. `params` is immutable at runtime -- assigning to it
    // is a no-op that raises nothing -- so this script can ask a question and has no way
    // to apply the answer. See the operational trap in CLAUDE.md.
    //
    // `--foo bar` is a pipeline parameter; `-profile local` is a Nextflow option and
    // does not count as configuring the run, which is why the test is for a DOUBLE dash.
    def cl = workflow.commandLine ?: ''
    if (!(cl =~ /(^|\s)--\S/) && !cl.contains('-params-file')) {
        log.info "No parameters given -- running defaults: " +
                 "workflow=${params.workflow}, taxon=${params.taxon}, " +
                 "outdir=${params.outdir}"
        log.info "To be asked for them instead:  python scripts/configure_run.py"
    }

    // The `local` profile carries the PORTABLE FLOOR -- resource requests sized for
    // the smallest supported machine, 16 GB / 4 cores -- because Nextflow's local
    // executor refuses to ADMIT a task that exceeds the machine, so any absolute in
    // a profile is a hidden minimum spec. A bigger box needs machine.config to be
    // used fully, and there is no way to tell from inside the run whether one was
    // generated and skipped or never generated, so this is a notice either way --
    // and it says which of the two it got.
    //
    // It used to speak up only when machine.config was ABSENT, so a run that WAS
    // sized said nothing at all and the only way to tell sizing had happened was to
    // go and read the wizard's scrollback.
    //
    // What it prints comes from the file's own header, not from the resolved config.
    // `workflow.session.config` does expose the per-label budgets, but
    // `memory = { 6.GB * task.attempt }` arrives there as an UNCALLED closure, and
    // `workflow.session` is not part of the documented metadata object -- reading it
    // to print a number would break on a version bump. The header is a plain file
    // written by the script that did the sizing, and it carries the measurements and
    // the date, so a stale config is visible rather than silent.
    def sized = workflow.configFiles.findAll { it.toString().contains('machine.config') }
    if ((workflow.profile ?: '').contains('local')) {
        if (sized) {
            log.info "Resources sized to this machine by ${sized.join(', ')}"
            sized.each { f ->
                try {
                    def folded = []
                    file(f.toString()).readLines()
                        .findAll { it =~ /^\s*\/\/\s*(Generated by|Detected:|Budgetable|BiG-SCAPE at|partition cap)/ }
                        .collect { it.replaceFirst(/^\s*\/\/\s*/, '') }
                        .collect { it.startsWith('Generated by') ? it.replaceFirst(/\.\s.*$/, '.') : it }
                        .each { ln ->
                            // The header wraps a long budget across two comment lines.
                            // A trailing comma is where it broke, so rejoin there --
                            // one budget, one line.
                            if (folded && folded[-1].endsWith(',')) {
                                folded[-1] = folded[-1] + ' ' + ln
                            }
                            else {
                                folded << ln
                            }
                        }
                    folded.each { log.info "  ${it}" }
                }
                catch (Exception e) {
                    // A header that will not read is no reason to stop a run. The line
                    // above has already named the file that did the sizing.
                }
            }
            log.info "  work directory: ${workDir}"
        }
        else {
            log.info "Resources are the portable floor (sized for 16 GB / 4 cores), not this machine."
            log.info "To size them here:  python scripts/check_machine.py   then add  -c machine.config"
        }
    }

    // Warnings
    if (params.run_gtdbtk) {
        // Corrected twice. The original read "~140 GB disk space and ~56-64 GB RAM".
        // On 2026-11-02 the RAM figure was cut to 4 GB from the RefSeq run's
        // disk_samples.tsv -- and that was wrong, because those were 60-SECOND
        // samples and the peak is a ~30-second spike. Sampling at 5 s and 2 s on
        // 2026-11-03 found ~37 GiB by three independent instruments, so the original
        // 56-64 GB was nearer the truth than the "correction" that replaced it.
        //
        // Disk still gates how much free space you need (~300 GB); RAM gates whether
        // the task is ADMITTED at all. See the GTDBTK_CLASSIFY block in
        // conf/labels.config for the measurements and the saturation curve.
        log.warn "GTDB-Tk is enabled: ~140 GB for the database plus ~150 GB of pplacer " +
                 "work scratch (measured), so plan for ~300 GB free."
        log.warn "It also needs ~37 GB of RAM in a ~30 s spike, so it is allocated " +
                 "40 GB and will not be admitted on a machine under ~44 GB. Set " +
                 "run_gtdbtk = false there. The 4 GB this warning used to quote came " +
                 "from a 60 s sampler that missed the spike entirely."
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

        // These three come from a prior DOWNLOAD_GENOMES run under the SAME
        // outdir and taxon. Pointing --input_genomes at another run's genomes
        // is the normal way to set up an A/B comparison, and then they are
        // simply absent -- which used to kill the run in AGGREGATE_TAXONOMY
        // with FileNotFoundError, after every antiSMASH task had completed.
        // They are metadata: detection and clustering do not need them, so a
        // missing one becomes a placeholder and the steps that need it skip.
        def tax_path = "${results_dir}/taxonomy_map.json"
        def name_path = "${base_dir}/name_map.json"
        def info_path = "${base_dir}/ncbi_dataset/data/assembly_info_table.txt"

        def taxonomy_map_ch = Channel.value(
            file(tax_path).exists() ? file(tax_path) : file('NO_TAXONOMY_MAP'))
        def name_map_ch = Channel.value(
            file(name_path).exists() ? file(name_path) : file('NO_NAME_MAP'))
        def assembly_info_ch = Channel.value(
            file(info_path).exists() ? file(info_path) : file('NO_ASSEMBLY_INFO'))

        if (!file(tax_path).exists()) {
            log.warn "No taxonomy_map.json under ${results_dir} — taxonomy " +
                     "aggregation will be skipped. This is expected when " +
                     "--input_genomes points at another run's genomes."
        }

        BGC_ANALYSIS(
            params.taxon,
            Channel.fromPath("${params.input_genomes}/*.gbff"),
            assembly_info_ch,
            name_map_ch,
            taxonomy_map_ch,
            false   // these genomes never passed through RENAME_GENOMES
        )

    } else if (params.workflow == "full") {
        DOWNLOAD_GENOMES(params.taxon)
        BGC_ANALYSIS(
            params.taxon,
            DOWNLOAD_GENOMES.out.renamed_genomes,
            DOWNLOAD_GENOMES.out.assembly_info,
            DOWNLOAD_GENOMES.out.name_map,
            DOWNLOAD_GENOMES.out.taxonomy_map,
            params.pepm_prescreen   // RENAME_GENOMES already screened
        )

    } else {
        error "ERROR: Invalid workflow '${params.workflow}'. Valid options: 'download', 'bgc_analysis', 'full'"
    }
}
