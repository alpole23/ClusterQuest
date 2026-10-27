/**
 * Resolve a taxon to an accession list and its metadata, WITHOUT the genomes.
 *
 * This is the first half of what NCBI_DATASETS_DOWNLOAD used to do in one task.
 * That task downloaded every genome in the taxon into one work directory, which
 * set the pipeline's peak disk: 8.7 MB per genome, all of it resident before a
 * single downstream step ran. At the ~150,000 genomes of RefSeq Enterobacterales
 * that is 1.3 TB in one directory, and no amount of screening afterwards helps,
 * because the screen cannot run until the download task has finished.
 *
 * The split is at a seam that already existed: `datasets download --dehydrated`
 * fetches the data report and a manifest of file locations but no payload, and
 * `datasets rehydrate` fetches the payload afterwards. Everything here is the
 * dehydrated half, so the whole stage costs ~3.7 kB per genome -- measured 10.3 MB
 * for 2,771 Erwiniaceae genomes, ~558 MB at Enterobacterales scale.
 *
 * FETCH_RENAME_SCREEN then takes the accession list in batches. It re-applies the
 * same `--assembly-source`, `--assembly-level` and `--exclude-atypical` filters
 * rather than trusting the list alone: `datasets download genome accession`
 * accepts the identical flag set (checked), so the two cannot drift, and a batch
 * download is self-describing if anyone reads it in isolation.
 *
 * `accessions.txt` is sorted, because batch membership is derived from it and an
 * unsorted list would reshuffle every batch between runs and defeat -resume.
 *
 * TWO CHANGES, both forced by GenBank scale (3,350,527 records):
 *
 * 1. `datasets summary --as-json-lines` replaces `download --dehydrated` + unzip.
 *    It yields the identical data report -- `dataformat` consumes it unchanged,
 *    same record count -- while streaming rather than delivering one atomic zip,
 *    and measured 1.9x faster on Bacteroidota (277 s against 514 s). The manifest
 *    the dehydrated mode also produced was never used: FETCH_RENAME_SCREEN
 *    downloads by accession, not by rehydrate.
 *
 * 2. The fetch is SHARDED BY RELEASE DATE above params.metadata_shard_threshold.
 *    A whole-domain request cannot be served in one stream: four consecutive
 *    attempts died with `stream error ... INTERNAL_ERROR; received from peer` at
 *    0.29 / 1.30 / 0.53 / 1.31 GB of ~1.72 GB, one to two hours each, with every
 *    retry restarting from zero. 149,539 records fetch cleanly in ~4.5 min, so
 *    the fix is to ask for less at a time.
 *
 *    Below the threshold this emits ONE shard with no date flags -- the identical
 *    single call it made before, because an extra `datasets` invocation costs
 *    ~0.9 s fixed and sharding a genus would be ~4x slower for nothing.
 *
 * See scripts/genome/plan_metadata_shards.py for why dates rather than sub-taxa.
 */
process NCBI_FETCH_METADATA {
    tag "$taxon"
    label 'process_low'
    label 'retry_on_error'
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode

    input:
    val taxon
    val scripts_version

    output:
    path "ncbi_dataset/data/assembly_info_table.txt", emit: assembly_info
    path "ncbi_dataset/data/assembly_data_report.jsonl", emit: assembly_data_report
    path "ncbi_dataset/data/taxonomy_report.jsonl", emit: taxonomy_report, optional: true
    path "accessions.txt", emit: accessions

    script:
    def level_flag = params.assembly_level ? "--assembly-level ${params.assembly_level}" : ''
    // Selection filter: applied here only. See the note in nextflow.config.
    def ref_flag = params.assembly_reference_only ? '--reference' : ''
    """
    mkdir -p ncbi_dataset/data

    # One line per shard, "<released-after>\\t<released-before>", either side may
    # be empty. A single blank line means one unsharded call -- what every taxon
    # under the threshold gets, identical to this process's previous behaviour.
    python3 ${projectDir}/scripts/genome/plan_metadata_shards.py \\
        --taxon "${taxon}" \\
        --assembly_source ${params.assembly_source} \\
        --threshold ${params.metadata_shard_threshold} \\
        --exclude_atypical > shards.tsv

    : > ncbi_dataset/data/assembly_data_report.jsonl
    n=0
    while IFS=\$'\\t' read -r AFTER BEFORE; do
        n=\$((n + 1))
        DATE_FLAGS=""
        [ "\$AFTER"  != "-" ] && DATE_FLAGS="\$DATE_FLAGS --released-after \$AFTER"
        [ "\$BEFORE" != "-" ] && DATE_FLAGS="\$DATE_FLAGS --released-before \$BEFORE"

        # Each shard retries independently, so a dropped stream costs one shard's
        # minutes rather than the whole taxon's hours. That is the entire point:
        # the previous single call restarted from zero on every retry.
        ok=0
        for attempt in 1 2 3; do
            if datasets summary genome taxon "${taxon}" \\
                    --assembly-source ${params.assembly_source} \\
                    --exclude-atypical \\
                    ${level_flag} \\
                    ${ref_flag} \\
                    \$DATE_FLAGS \\
                    --as-json-lines > shard.jsonl 2> shard.err; then
                cat shard.jsonl >> ncbi_dataset/data/assembly_data_report.jsonl
                ok=1
                break
            fi
            echo "  shard \$n attempt \$attempt failed: \$(tail -1 shard.err)" >&2
            sleep \$((attempt * 20))
        done
        if [ "\$ok" -ne 1 ]; then
            echo "ERROR: shard \$n (\$AFTER .. \$BEFORE) failed three times"
            exit 1
        fi
        rm -f shard.jsonl
    done < shards.tsv
    echo "fetched \$n shard(s) for ${taxon}"

    # Deduplicate BEFORE dataformat. Date shards are disjoint in intent but the
    # bounds are inclusive on both sides, so a genome released exactly on a cut
    # date lands in two shards -- measured: 6 shards over Pantoea gave 1,808
    # records for 1,807 genomes. accessions.txt would absorb that (it is a set),
    # but the duplicate would reach the assembly table and CREATE_NAME_MAP as a
    # real second row. First occurrence wins; shards are concatenated in date
    # order, so that is the earlier range.
    python3 - <<'PY'
import json
seen, kept = set(), 0
src = 'ncbi_dataset/data/assembly_data_report.jsonl'
with open(src) as fh, open(src + '.dedup', 'w') as out:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        acc = json.loads(line).get('accession')
        if acc in seen:
            continue
        seen.add(acc)
        out.write(line + '\\n')
        kept += 1
print(f'{kept} unique records after de-duplicating shard boundaries')
PY
    mv ncbi_dataset/data/assembly_data_report.jsonl.dedup ncbi_dataset/data/assembly_data_report.jsonl

    dataformat tsv genome \\
        --inputfile ncbi_dataset/data/assembly_data_report.jsonl \\
        --fields accession,organism-name,organism-infraspecific-strain,assminfo-biosample-isolation-source,assminfo-biosample-isolate,assminfo-notes,checkm-completeness,checkm-contamination,organism-tax-id \\
        > ncbi_dataset/data/assembly_info_table.txt

    # Sorted: FETCH_RENAME_SCREEN derives batch membership from this file, and an
    # unsorted list would put a genome in a different batch on every run. The set
    # also de-duplicates, which now matters because shards are concatenated --
    # date ranges are disjoint by construction (verified on Pantoea: 310 + 1,497
    # = 1,807, zero shared, zero missing), but a duplicate must never become two
    # batch entries.
    python3 - <<'PY'
import json
seen = set()
with open('ncbi_dataset/data/assembly_data_report.jsonl') as fh:
    for line in fh:
        line = line.strip()
        if line:
            acc = json.loads(line).get('accession')
            if acc:
                seen.add(acc)
with open('accessions.txt', 'w') as out:
    out.write('\\n'.join(sorted(seen)) + '\\n')
print(f'{len(seen)} accessions resolved')
PY

    if [ ! -s accessions.txt ]; then
        echo "ERROR: no accessions resolved for taxon '${taxon}'"
        exit 1
    fi

    # Taxonomy is a separate, small download. A failure here loses the report's
    # taxonomy tree and nothing else, so it must not fail the run -- but it also
    # must not be silent, which an `|| exit 0` in the previous version made it.
    if datasets download taxonomy taxon "${taxon}" --filename taxonomy.zip 2>/dev/null; then
        unzip -o -q taxonomy.zip
    else
        echo "WARNING: taxonomy download failed for '${taxon}'; the report will have no taxonomy tree"
    fi
    """
}
