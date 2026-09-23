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
 */
process NCBI_FETCH_METADATA {
    tag "$taxon"
    label 'process_low'
    label 'retry_on_error'
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode

    input:
    val taxon

    output:
    path "ncbi_dataset/data/assembly_info_table.txt", emit: assembly_info
    path "ncbi_dataset/data/assembly_data_report.jsonl", emit: assembly_data_report
    path "ncbi_dataset/data/taxonomy_report.jsonl", emit: taxonomy_report, optional: true
    path "accessions.txt", emit: accessions

    script:
    def level_flag = params.assembly_level ? "--assembly-level ${params.assembly_level}" : ''
    """
    # Dehydrated: data report and file manifest, no genome payload.
    datasets download genome taxon "${taxon}" \\
        --include gbff \\
        --assembly-source ${params.assembly_source} \\
        --exclude-atypical \\
        ${level_flag} \\
        --filename ncbi_dataset.zip \\
        --dehydrated

    unzip -q ncbi_dataset.zip

    dataformat tsv genome \\
        --inputfile ncbi_dataset/data/assembly_data_report.jsonl \\
        --fields accession,organism-name,organism-infraspecific-strain,assminfo-biosample-isolation-source,assminfo-biosample-isolate,assminfo-notes,checkm-completeness,checkm-contamination,organism-tax-id \\
        > ncbi_dataset/data/assembly_info_table.txt

    # Sorted: FETCH_RENAME_SCREEN derives batch membership from this file, and an
    # unsorted list would put a genome in a different batch on every run.
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
print(f'{len(seen)} accessions resolved for ${taxon}')
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
