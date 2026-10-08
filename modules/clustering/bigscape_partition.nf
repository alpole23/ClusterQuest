/**
 * BiG-SCAPE over one pepM partition.
 *
 * Region GenBanks arrive staged flat, but are re-laid-out into one directory per
 * genome before BiG-SCAPE sees them. BiG-SCAPE records the directory it read each
 * GBK from, and downstream consumers derive the genome from that path — flat
 * input makes every BGC look like it came from `part_input`, which sends
 * bgc_coupling_annotation.py into a fallback that rescans every genome's JSON for
 * every BGC. That took a 2 h task limit and two killed runs to find.
 *
 * **The reference pass happens here, per partition, not after the merge.** BiG-SCAPE
 * computes only the pairs absent from the distance table, so a reference pass over
 * the MERGED database would compute every cross-partition query pair -- ~34 M at a
 * 2-way split of 11,700 BGCs, against the ~199 k actually wanted. That is the n^2
 * partitioning exists to avoid. Each partition database is already complete within
 * itself, so adding references there costs refs x partition size, summing over
 * partitions to exactly the refs x N a monolithic pass would do.
 *
 * It runs here rather than in a process of its own because the input is already
 * staged and already laid out as part_input/<genome>/; a separate process would
 * stage every GBK a second time to rebuild the same tree.
 */
process BIGSCAPE_PARTITION {
    tag "${taxon} part ${part_id}"
    label 'process_high'

    input:
    val taxon
    tuple val(part_id), path(gbks)
    path pfam_db
    path partitions

    // The characterised reference clusters. A run without them stages a placeholder
    // file instead of a directory, the guard below finds no .gbk, and no reference
    // pass runs -- so `refdb` is an optional output rather than an empty database.
    path reference_dir

    // Digest of the reference GenBanks. A val input, not an interpolated path:
    // Nextflow hashes a directory input by its own metadata rather than its
    // contents, so adding a reference would otherwise leave -resume reusing a run
    // without it. Same reasoning as BIGSCAPE_REFERENCES.
    val reference_version

    output:
    path "part_${part_id}.db", emit: db
    path "refpass_${part_id}.db", emit: refdb, optional: true

    script:
    """
    export PFAM_PATH=\$(readlink -f ${pfam_db})/Pfam-A.hmm
    export PYTHONHASHSEED=0  # see BIGSCAPE: load order decides pair orientation

    # Rebuild <genome>/<region>.gbk from the manifest's genome column.
    mkdir -p part_input
    python - <<'LAYOUT'
import csv, os, shutil
genome = {}
with open('${partitions}') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        genome[os.path.basename(row['gbk'])] = row['genome']
missing = []
for f in os.listdir('.'):
    if not f.endswith('.gbk'):
        continue
    g = genome.get(f)
    if g is None:
        missing.append(f)
        continue
    os.makedirs(os.path.join('part_input', g), exist_ok=True)
    shutil.copy2(f, os.path.join('part_input', g, f))
if missing:
    raise SystemExit(f'{len(missing)} staged GBKs absent from the manifest: {missing[:3]}')
LAYOUT

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

    # Reference pass, on a COPY. part_${part_id}.db goes to the merge and must stay
    # dataset-only: a reference inside the GCF cutoff of a query would join its
    # family, and everything downstream reads the database without knowing which
    # records are the dataset. Measured on a 14-genome subset, one leaked reference
    # raised total_bgcs from 18 to 19 and invented a genome named after the
    # reference directory.
    if [ -d "${reference_dir}" ] && ls ${reference_dir}/*.gbk >/dev/null 2>&1; then
        cp part_${part_id}.db refpass_${part_id}.db
        # --db-only-output skips HTML and tree regeneration nothing here reads, and
        # was verified to leave every reference distance identical.
        bigscape cluster \\
            -i part_input \\
            -o refout_${part_id} \\
            --db-path refpass_${part_id}.db \\
            --reference-dir ${reference_dir} \\
            --pfam-path \$PFAM_PATH \\
            --alignment-mode ${params.bigscape_alignment_mode} \\
            --gcf-cutoffs ${params.bigscape_cutoffs} \\
            ${params.bigscape_include_singletons ? '--include-singletons' : ''} \\
            ${params.bigscape_mix ? '--mix' : ''} \\
            ${params.bigscape_classify ? "--classify ${params.bigscape_classify}" : ''} \\
            --db-only-output \\
            --cores ${task.cpus}
    fi
    """
}
