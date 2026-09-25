/**
 * Fetch one batch of accessions, rename them, screen them, keep the survivors.
 *
 * The second half of the old NCBI_DATASETS_DOWNLOAD, fused with RENAME_GENOMES.
 * Peak disk was the whole reason: the old task held every genome in the taxon at
 * once. Here a batch is fetched, screened and reduced to its survivors before the
 * next one starts, so transient disk is `maxForks x batch x 8.7 MB` — about
 * 3.5 GB at the defaults, against 1.3 TB for a 150,000-genome taxon.
 *
 * ## maxForks 2 is measured, not guessed
 *
 * Probed against NCBI on 2026-09-23, four batches of eight at each level:
 *
 *   concurrency 1  23.5 MB/s   112 genomes/min   0 failures
 *   concurrency 2  33.5 MB/s   159 genomes/min   0 failures
 *   concurrency 4  29.0 MB/s   138 genomes/min   0 failures
 *   concurrency 8  25.1 MB/s   119 genomes/min   0 failures
 *
 * Throughput PEAKS at 2 and declines above it — the signature of a saturated link
 * (~264 Mbit/s here) rather than server-side throttling. So concurrency above 2
 * buys nothing, costs a little, and is needlessly rude to NCBI. It also means
 * batching this stage is free in wall-clock terms: two fetches already saturate
 * the connection, so splitting one download into hundreds changes throughput not
 * at all.
 *
 * Raise `--download_max_forks` only on a fatter connection, and re-probe first.
 *
 * ## The burst probe's "zero failures" was underpowered
 *
 * A 30-minute soak at concurrency 2 (148 batches, 2,960 genome-fetches, 7.5 GB)
 * found **12 failures, 8.1%** — flat across every 5-minute bucket rather than
 * rising, so a constant transient rate and not progressive throttling. The probe
 * saw 0 of 16 because 16 samples cannot distinguish 0% from 8%, not because the
 * rate was lower. Retries are therefore load-bearing here, not belt-and-braces:
 * `retry_on_error` plus the refetch loop below take the residual to ~0.05%.
 *
 * Sustained throughput was ~98 genome-fetches/min against the burst's 159, with a
 * mild downward drift (108 -> 80 over the half hour) that is within the noise of a
 * cycling accession pool. Budget order-scale downloads at the SUSTAINED figure:
 * 150,000 genomes is ~25 h, not the ~8 h the burst implied.
 *
 * The cause is now known. Captured with full stderr, 3 of 40 batches failed with
 *
 *   Downloading: d.zip  32.8MB invalid zip archive
 *   Error: Internal error (invalid zip archive). Please try again
 *
 * identical in all three: the archive transfers in full (33-43 MB) and then fails
 * NCBI's own validation partway through. Transfer corruption, transient, and the
 * client says to retry — the same failure class as the null-byte files the loop
 * further down was written for. Handled by the inline retry in the script.
 *
 * ## Corruption handling is carried over deliberately
 *
 * NCBI's dehydrated transfers can deliver null-filled files on a timeout. The
 * validate-and-refetch loop below is the old module's, narrowed to a batch: a
 * corrupt file is deleted and refetched, and after the last cycle any survivor
 * that is still null is dropped rather than passed downstream as a genome with no
 * sequence. This logic was earned by a real failure and is not worth rewriting.
 *
 * Both outputs are optional: a batch whose genomes all fail the screen emits
 * nothing, which at Erwiniaceae's 11% pass rate is the common case, not an error.
 */
process FETCH_RENAME_SCREEN {
    tag "${accessions.name}"
    label 'process_medium'
    label 'retry_on_error'
    maxForks params.download_max_forks
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(params.taxon)}/renamed_genomes",
        mode: params.publish_mode, pattern: '*.gbff'
    publishDir "${params.outdir}/prescreen_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'prescreen_*.tsv'

    input:
    val taxon
    path accessions
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
    def level_flag = params.assembly_level ? "--assembly-level ${params.assembly_level}" : ''
    def screen = params.pepm_prescreen && Utils.isValidInput(pepm_db)
    """
    # Same filters as NCBI_FETCH_METADATA. Redundant given the accession list came
    # from that call, and kept anyway so the two cannot silently diverge.
    #
    # Retried inline because this fails ~7.5% of the time with
    #   Error: Internal error (invalid zip archive). Please try again
    # measured 3 of 40 batches, matching the soak's 8.1%. The zip transfers in full
    # (33-43 MB) and then fails NCBI's own validation partway through, so it is
    # transfer corruption and the client itself says to retry. The `retry_on_error`
    # label would also catch it, but at the cost of re-running the whole task for a
    # bad download; three attempts here take 7.5% to ~0.04% without leaving the task.
    FETCHED=no
    for ATTEMPT in 1 2 3; do
        rm -f batch.zip
        if datasets download genome accession \\
                --inputfile ${accessions} \\
                --include gbff \\
                --assembly-source ${params.assembly_source} \\
                --exclude-atypical \\
                ${level_flag} \\
                --filename batch.zip \\
            && unzip -qt batch.zip >/dev/null 2>&1; then
            FETCHED=yes
            break
        fi
        echo "fetch attempt \$ATTEMPT failed (invalid or incomplete archive); retrying"
        sleep \$(( ATTEMPT * 10 ))
    done
    # Falling back to one accession at a time is not belt-and-braces: a whole-batch
    # fetch is all-or-nothing, and CLAUDE.md's rule for every batched process here is
    # that a batch must exit non-zero only when EVERY item failed, or batching turns
    # one bad download into fifty lost genomes. It did exactly that on a P. ananatis
    # run -- one 130 MB archive failed all three inline attempts and all three task
    # retries, killing the run, while the same 51 accessions fetched cleanly in
    # 17-accession chunks minutes later. Risk scales with archive size, so the
    # fallback shrinks the request rather than just repeating it.
    if [ "\$FETCHED" != yes ]; then
        echo "whole-batch fetch failed 3 times; falling back to one accession at a time"
        mkdir -p singles
        SOLO_OK=0
        while read -r ACC; do
            [ -n "\$ACC" ] || continue
            for TRY in 1 2; do
                if datasets download genome accession "\$ACC" --include gbff \\
                        --assembly-source ${params.assembly_source} \\
                        --filename "singles/\$ACC.zip" >/dev/null 2>&1 \\
                   && unzip -qt "singles/\$ACC.zip" >/dev/null 2>&1; then
                    unzip -qo "singles/\$ACC.zip" -x README.md md5sum.txt >/dev/null 2>&1 || true
                    SOLO_OK=\$((SOLO_OK + 1))
                    break
                fi
                rm -f "singles/\$ACC.zip"
                sleep 5
            done
        done < ${accessions}
        rm -rf singles
        if [ "\$SOLO_OK" -eq 0 ]; then
            echo "ERROR: every accession in this batch failed to fetch"
            exit 1
        fi
        echo "recovered \$SOLO_OK of \$(wc -l < ${accessions}) accessions individually"
    fi

    # Only the whole-batch path leaves an archive to open; the fallback already
    # unpacked each accession as it went.
    if [ "\$FETCHED" = yes ]; then
        unzip -q batch.zip
    fi
    sync

    # Null-filled files are what a timed-out transfer leaves behind. Delete and
    # refetch; a valid GenBank file starts with LOCUS.
    for CYCLE in 1 2 3; do
        CORRUPT=0
        for GBFF in ncbi_dataset/data/*/*.gbff; do
            [ -f "\$GBFF" ] || continue
            if [ -z "\$(head -c 100 "\$GBFF" 2>/dev/null | tr -d '\\0' | head -c 5)" ]; then
                echo "corrupt, will refetch: \$GBFF"
                rm -f "\$GBFF"
                CORRUPT=\$((CORRUPT + 1))
            fi
        done
        [ "\$CORRUPT" -eq 0 ] && break
        echo "cycle \$CYCLE: \$CORRUPT corrupt, refetching"
        rm -f batch.zip
        datasets download genome accession --inputfile ${accessions} --include gbff \\
            --assembly-source ${params.assembly_source} --filename batch.zip || true
        unzip -qo batch.zip || true
        sync
    done

    FOUND=\$(find ncbi_dataset/data -name '*.gbff' -type f | wc -l)
    WANTED=\$(wc -l < ${accessions})
    echo "fetched \$FOUND of \$WANTED genomes in this batch"
    if [ "\$FOUND" -eq 0 ]; then
        echo "ERROR: no usable genomes in this batch after 3 fetch cycles"
        exit 1
    fi

    # Pair each .gbff with its accession, which is its parent directory's name.
    python3 - <<'PY' > manifest.tsv
import pathlib
for p in sorted(pathlib.Path('ncbi_dataset/data').glob('*/*.gbff')):
    print(f'{p.parent.name}\\t{p}')
PY

    python ${projectDir}/scripts/genome/rename_genome.py manifest.tsv ${name_map}

    if ${screen}; then
        # Renamed files land in the task directory; the fetched originals stay under
        # ncbi_dataset/, so a bare *.gbff glob here is unambiguous — unlike in
        # RENAME_GENOMES, where the staged inputs share the extension.
        python ${projectDir}/scripts/analysis/pepm_prescreen.py \\
            --genomes *.gbff \\
            --db ${pepm_db} \\
            --out prescreen_${task.index}.tsv \\
            --bitscore ${params.pepm_prescreen_bitscore} \\
            --min_density ${params.pepm_prescreen_min_density}

        python3 - <<'PRUNE_EOF'
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

    # The fetched originals are never an output; dropping them here keeps the task
    # directory at its survivors, so work/ does not accumulate the ~89% of genomes
    # the screen rejected.
    rm -rf ncbi_dataset batch.zip
    """
}
