process GCF_BIOSYNTHETIC_TREE {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}/gcf_heatmap", mode: 'copy'

    input:
    val taxon
    path bigscape_db
    path "antismash_input/*", stageAs: 'antismash_input/*'
    path gtdbtk_tree
    path gtdbtk_summary

    output:
    path "gcf_species_heatmap.png",                     emit: heatmap_png,          optional: true
    path "gcf_species_heatmap.svg",                     emit: heatmap_svg,          optional: true
    path "gcf_biosynthetic_tree.png",                   emit: gcf_tree_png,         optional: true
    path "gcf_biosynthetic_tree.svg",                   emit: gcf_tree_svg,         optional: true
    path "all_bgcs_biosynthetic_tree_circular.png",     emit: all_bgcs_tree_png,    optional: true
    path "all_bgcs_biosynthetic_tree_circular.svg",     emit: all_bgcs_tree_svg,    optional: true
    path "phosphonate_metadata.json",                   emit: metadata,             optional: true
    path "phosphonate_itol_coupling.txt",               emit: coupling_annotation,  optional: true
    path "phosphonate_coupling_support.tsv",            emit: coupling_support,     optional: true

    script:
    def tree_arg    = Utils.optArg('--gtdbtk_tree',    gtdbtk_tree)
    def summary_arg = Utils.optArg('--gtdbtk_summary', gtdbtk_summary)
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['bgc_all_bgcs_tree.py', 'bgc_coupling_annotation.py', 'bgc_gcf_heatmap.py', 'bgc_gcf_tree.py', 'bgc_pfam_tree.py', 'utils'])}
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
        ${tree_arg} \\
        ${summary_arg} \\
        --outdir .

    # Step 4: Generate GCF biosynthetic NJ tree figure (GCF medoids)
    python ${projectDir}/scripts/bgc_gcf_tree.py \\
        --db ${bigscape_db} \\
        --coupling_annotation phosphonate_itol_coupling.txt \\
        --outdir .

    # Step 5: Generate all-BGCs circular NJ tree (full distance matrix)
    python ${projectDir}/scripts/bgc_all_bgcs_tree.py \\
        --db ${bigscape_db} \\
        --coupling_annotation phosphonate_itol_coupling.txt \\
        --outdir . \\
        --layout circular
    """
}
