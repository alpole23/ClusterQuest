process CREATE_NAME_MAP {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(params.taxon)}", mode: 'copy'

    input:
    val taxon
    path assembly_info

    output:
    path "name_map.json", emit: name_map

    script:
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['genome/create_name_map.py'])}
    python ${projectDir}/scripts/genome/create_name_map.py ${assembly_info} name_map.json
    """
}
