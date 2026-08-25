process COUPLING_ENZYME_TREE {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}/coupling_enzyme_trees", mode: 'copy'

    input:
    val taxon
    path "antismash_input/*", stageAs: 'antismash_input/*'
    path metadata
    path coupling_annotation
    path ref_pepm_faa
    path ref_coupling_faa

    output:
    path "tree_A",                        emit: tree_a,   optional: true
    path "tree_B",                        emit: tree_b,   optional: true
    path "coupling_trees_manifest.json",  emit: manifest, optional: true

    script:
    """
    python ${projectDir}/scripts/bgc_coupling_tree.py \\
        --antismash_dir       antismash_input \\
        --metadata            ${metadata} \\
        --coupling_annotation ${coupling_annotation} \\
        --ref_pepm_faa        ${ref_pepm_faa} \\
        --ref_coupling_faa    ${ref_coupling_faa} \\
        --outdir              . \\
        --hmmbuild            \$(command -v hmmbuild) \\
        --hmmalign            \$(command -v hmmalign) \\
        --hmmsearch           \$(command -v hmmsearch) \\
        --fasttree            \$(command -v FastTree || command -v fasttree) \\
        --tree                ${params.coupling_tree_type}

    # Drop HMM intermediates so the published output matches the documented layout
    rm -rf tree_A/hmm_work tree_B/*/hmm_work
    """
}
