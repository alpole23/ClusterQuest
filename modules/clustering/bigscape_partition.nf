/**
 * BiG-SCAPE over one pepM partition.
 *
 * Region GenBanks are staged flat: BiG-SCAPE globs recursively and filters on
 * filenames containing "region"/"cluster", so directory structure is not needed.
 * PARTITION_BGCS fails loudly if two region filenames collide, which is what
 * makes flat staging safe.
 */
process BIGSCAPE_PARTITION {
    tag "${taxon} part ${part_id}"
    label 'process_high'

    input:
    val taxon
    tuple val(part_id), path(gbks)
    path pfam_db

    output:
    path "part_${part_id}.db", emit: db

    script:
    """
    export PFAM_PATH=\$(readlink -f ${pfam_db})/Pfam-A.hmm
    mkdir -p part_input && cp ${gbks} part_input/

    bigscape cluster \\
        -i part_input \\
        -o out_${part_id} \\
        --pfam-path \$PFAM_PATH \\
        --alignment-mode ${params.bigscape_alignment_mode} \\
        --gcf-cutoffs ${params.bigscape_cutoffs} \\
        ${params.bigscape_include_singletons ? '--include-singletons' : ''} \\
        ${params.bigscape_mix ? '--mix' : ''} \\
        --cores ${task.cpus} \\
        ${params.bigscape_classify ? "--classify ${params.bigscape_classify}" : ''}

    # One database per partition, named so MERGE_BIGSCAPE sorts them stably.
    cp \$(find out_${part_id} -name '*.db' | head -1) part_${part_id}.db
    """
}
