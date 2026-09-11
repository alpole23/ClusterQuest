process GCF_BIOSYNTHETIC_TREE {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}/gcf_heatmap", mode: params.publish_mode

    input:
    val taxon
    path bigscape_db
    path "antismash_input/*", stageAs: 'antismash_input/*'
    path gtdbtk_summary
    path gtdbtk_db
    path centers_db

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "gcf_species_heatmap.png",                     emit: heatmap_png,          optional: true
    path "gcf_species_heatmap.svg",                     emit: heatmap_svg,          optional: true
    path "gcf_biosynthetic_tree.png",                   emit: gcf_tree_png,         optional: true
    path "gcf_biosynthetic_tree.svg",                   emit: gcf_tree_svg,         optional: true
    path "phosphonate_metadata.json",                   emit: metadata,             optional: true
    path "phosphonate_itol_coupling.txt",               emit: coupling_annotation,  optional: true
    path "phosphonate_coupling_support.tsv",            emit: coupling_support,     optional: true

    script:
    def summary_arg = Utils.optArg('--gtdbtk_summary', gtdbtk_summary)
    // GTDB's own reference phylogeny, shipped in the data package. Orders the heatmap
    // columns by genus without needing a per-run tree, which sharded GTDBTK_CLASSIFY
    // no longer produces.
    def ref_tree = Utils.isValidInput(gtdbtk_db)
        ? "--gtdb_reference_tree ${gtdbtk_db}/release${params.gtdb_release}/pplacer/gtdb_r${params.gtdb_release}_bac120.refpkg/gtdb_r${params.gtdb_release}_bac120_decorated_unrooted.tree"
        : ''
    // Exact centre-to-centre distances when the main database is
    // partitioned; without it the GCF tree substitutes 1.0 for every
    // cross-partition centre pair and its backbone is arbitrary.
    def centers_arg = Utils.optArg('--centers_db', centers_db)
    """
    # Step 1: Generate BGC metadata from BiG-SCAPE database (metadata only — skip slow NJ tree)
    python ${projectDir}/scripts/bgc_pfam_tree.py \\
        --db ${bigscape_db} \\
        --bgc_type phosphonate \\
        --metadata_only \\
        --outdir .

    # Step 2: Generate coupling enzyme class annotations
    python ${projectDir}/scripts/bgc_coupling_annotation.py \\
        --antismash_dir antismash_input \\
        --metadata phosphonate_metadata.json \\
        --outfile phosphonate_itol_coupling.txt \\
        --reference_faa ${projectDir}/assets/reference_sequences/reference_coupling_enzymes.faa \\
        --reference_pepm ${projectDir}/assets/reference_sequences/reference_pepM.faa \\
        --bgc_type phosphonate

    # Step 3: Generate GCF × species heatmap
    python ${projectDir}/scripts/bgc_gcf_heatmap.py \\
        --metadata phosphonate_metadata.json \\
        --coupling_annotation phosphonate_itol_coupling.txt \\
        ${summary_arg} \\
        ${ref_tree} \\
        --outdir .

    # Step 4: Generate GCF biosynthetic NJ tree figure (GCF medoids)
    python ${projectDir}/scripts/bgc_gcf_tree.py \\
        ${centers_arg} \\
        --db ${bigscape_db} \\
        --coupling_annotation phosphonate_itol_coupling.txt \\
        --outdir .
    """
}
