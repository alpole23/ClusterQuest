/**
 * Rename a batch of genome files, and drop the ones that cannot hold a BGC.
 *
 * Genomes are processed in batches (params.task_batch_size) because the per-genome
 * work is ~0.2 s — far below the scheduler overhead of one job per genome.
 * Files are staged as genome1.gbff, genome2.gbff, ... (their NCBI names all collide
 * on "genomic.gbff"); staging order matches the assembly_ids list.
 *
 * ## Why the pepM screen runs HERE and not after
 *
 * It used to run as its own stage downstream, filtering the channel. That works for
 * scheduling but not for disk: this process published *every* renamed genome, so a
 * run kept 8.7 MB per genome for genomes the screen had already ruled out. Measured
 * on the Erwiniaceae verification run, 24 GB of a 27 GB result directory was renamed
 * genomes, and 2,464 of those 2,771 genomes were rejected by the screen and used for
 * nothing.
 *
 * Deleting them afterwards is not an option. `publishDir` re-publishes from `work/`
 * on `-resume`, which silently undoes any prune, and pruning published output races
 * with tasks still publishing — the two reasons `prune_antismash_results.py` refuses
 * to run during a pipeline run. Screening inside the task that renames means a
 * rejected genome is gone *before* its output is declared, so `publishDir` never
 * sees it and `-resume` cannot bring it back.
 *
 * The screen is unchanged: same script, same references, same thresholds. Only
 * where it runs moved. Genomes supplied through `--input_genomes` never pass
 * through here, so ANTISMASH_ANALYSIS keeps its own PEPM_PRESCREEN for that path.
 *
 * Both outputs are `optional` because a batch in a phosphonate-poor clade can
 * legitimately have zero survivors — at Erwiniaceae's 11% pass rate that is
 * expected, not an error.
 */
process RENAME_GENOMES {
    tag "${assembly_ids.size()} genomes"
    label 'process_medium'
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(params.taxon)}/renamed_genomes",
        mode: params.publish_mode, pattern: '*.gbff'
    publishDir "${params.outdir}/prescreen_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'prescreen_*.tsv'

    input:
    val taxon
    tuple val(assembly_ids), path(genome_files, stageAs: 'genome?.gbff')
    path name_map
    // diamond database of pepM references, or a placeholder when screening is off
    path pepm_db

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "*.gbff", emit: renamed_genome, optional: true
    path "prescreen_${task.index}.tsv", emit: report, optional: true

    script:
    def staged = [genome_files].flatten().collect { it.name }
    def manifest = [assembly_ids, staged].transpose().collect { id, f -> "${id}\t${f}" }.join('\n')
    def screen = params.pepm_prescreen && Utils.isValidInput(pepm_db)
    """
    cat > manifest.tsv <<'MANIFEST_EOF'
${manifest}
MANIFEST_EOF

    # The staged inputs are themselves *.gbff (genome1.gbff, genome2.gbff, ...), so a
    # bare glob after renaming matches both them and the renamed copies -- the first
    # version of this screened all 16 files of an 8-genome batch. Nextflow excludes
    # staged inputs from the OUTPUT glob, which is why only this needed fixing.
    ls -1 *.gbff 2>/dev/null | sort > .staged.txt
    python ${projectDir}/scripts/genome/rename_genome.py manifest.tsv ${name_map}
    ls -1 *.gbff 2>/dev/null | sort > .all.txt
    comm -13 .staged.txt .all.txt > .renamed.txt
    echo "renamed \$(wc -l < .renamed.txt) genomes"

    if ${screen}; then
        # Screen the renamed batch, then delete everything that did not pass, so
        # the rejected genomes are never declared as output and never published.
        python ${projectDir}/scripts/analysis/pepm_prescreen.py \\
            --genomes \$(cat .renamed.txt | tr '\\n' ' ') \\
            --db ${pepm_db} \\
            --out prescreen_${task.index}.tsv \\
            --bitscore ${params.pepm_prescreen_bitscore} \\
            --min_density ${params.pepm_prescreen_min_density}

        python - <<'PRUNE_EOF'
import csv, os
kept = dropped = 0
with open("prescreen_${task.index}.tsv") as fh:
    for row in csv.DictReader(fh, delimiter='\\t'):
        path = row['genome'] + '.gbff'
        if row['pass'] == 'yes':
            kept += 1
        elif os.path.exists(path):
            os.remove(path)
            dropped += 1
print(f'pepM screen: {kept} kept, {dropped} dropped before publishing')
PRUNE_EOF
    fi
    """
}
