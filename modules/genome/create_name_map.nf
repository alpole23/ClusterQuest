process CREATE_NAME_MAP {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path assembly_info

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "name_map.json", emit: name_map

    script:
    """
    python ${projectDir}/scripts/genome/create_name_map.py ${assembly_info} name_map.json
    """
}
