/**
 * Rebuild genes that a genome's submitted annotation left out.
 *
 * Partial annotation is common and silently destructive. Measured on Erwiniaceae:
 * 327 of 333 BGC regions contain an intergenic gap of 300+ bp, the median region is
 * 78.4% coding where a well-annotated operon is ~90%, and 914 kb inside BGC regions
 * carries no gene call. The worst case is P. ananatis LMG 5342 region 2 -- a
 * lab-confirmed phosphonolipid cluster -- where 9 of 15 genes were never called,
 * including the 2-AEP transaminase that decides the product's headgroup.
 *
 * antiSMASH will not fix this itself: it runs gene finding only on records with ZERO
 * CDS features, so one CDS anywhere on a chromosome disables prodigal entirely. But
 * --genefinding-gff3 merges extra calls into an already-annotated record, which is
 * what these processes feed it.
 *
 * Runs after the pepM pre-screen, so the pool is drawn from the genomes that carry a
 * phosphonate pathway rather than from everything downloaded. That is cheaper and
 * also better: screened genomes are each other's best reference for exactly the genes
 * this is trying to recover.
 */

/**
 * Pool every protein already annotated across the screened genomes.
 *
 * This is the reference set for the homology half of recovery. The premise is that
 * annotation quality is uneven between submissions -- LMG 5342's 2012 deposit is
 * missing aepZ, P. ananatis VY148's 2021 assembly has it -- so a gene absent from one
 * genome is usually present in a sibling.
 *
 * A gather step: every RECOVER_ORFS task waits on this one. It is cheap (concatenating
 * proteins out of a few hundred GenBanks) but it does serialise here.
 */
process BUILD_PROTEIN_POOL {
    tag "$taxon"
    label 'process_low'

    input:
    val taxon
    path genomes

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "protein_pool.faa", emit: pool

    script:
    """
    python ${projectDir}/scripts/genome/build_protein_pool.py \\
        --genomes ${genomes} \\
        --out protein_pool.faa
    """
}

/**
 * Recover missing genes for a batch of genomes, as GFF3 for antiSMASH.
 *
 * Both prodigal and homology, because neither suffices alone. On LMG 5342 prodigal
 * found 8 of the 9 missing genes but truncated aepZ to 162 aa where the gene is 238 --
 * and aepZ is the gene the whole exercise is about. Homology recovers it at 100%
 * identity and names it.
 *
 * The GFF3s are published because they are the evidence for every downstream gene call
 * that was not in the source annotation.
 */
process RECOVER_ORFS {
    tag "${genomes instanceof List ? genomes.size() + ' genomes' : genomes.baseName}"
    label 'process_medium'
    publishDir "${params.outdir}/recovered_orfs/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: '*.gff3'

    input:
    val taxon
    path genomes
    path protein_pool
    val scripts_version

    output:
    path "*.gff3",            emit: gff
    path "recovery_stats.tsv", emit: stats, optional: true

    script:
    """
    for GENOME in ${genomes}; do
        BASE=\$(basename "\$GENOME" .gbff)
        # Every failure is a continue, never an exit: aborting would forfeit the
        # whole batch over one bad genome, and a genome with no recovered genes is
        # a normal outcome rather than an error.
        python ${projectDir}/scripts/genome/recover_orfs.py \\
            --genome "\$GENOME" \\
            --out "\${BASE}.gff3" \\
            --pool ${protein_pool} \\
            --prodigal \$(which prodigal) \\
            --diamond \$(which diamond) \\
            --threads ${task.cpus} \\
            --workdir "_work_\${BASE}" \\
            >> recovery_stats.tsv 2>&1 || {
                echo "  WARN: recovery failed for \$BASE, writing an empty GFF3"
                printf '##gff-version 3\\n' > "\${BASE}.gff3"
            }
        rm -rf "_work_\${BASE}"
    done
    """
}
