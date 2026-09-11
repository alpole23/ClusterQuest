/**
 * BiG-SCAPE over one pepM partition.
 *
 * Region GenBanks arrive staged flat, but are re-laid-out into one directory per
 * genome before BiG-SCAPE sees them. BiG-SCAPE records the directory it read each
 * GBK from, and downstream consumers derive the genome from that path — flat
 * input makes every BGC look like it came from `part_input`, which sends
 * bgc_coupling_annotation.py into a fallback that rescans every genome's JSON for
 * every BGC. That took a 2 h task limit and two killed runs to find.
 */
process BIGSCAPE_PARTITION {
    tag "${taxon} part ${part_id}"
    label 'process_high'

    input:
    val taxon
    tuple val(part_id), path(gbks)
    path pfam_db
    path partitions

    output:
    path "part_${part_id}.db", emit: db

    script:
    """
    export PFAM_PATH=\$(readlink -f ${pfam_db})/Pfam-A.hmm

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
    """
}
