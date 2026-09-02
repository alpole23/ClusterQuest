/**
 * All-by-all pepM comparison against gene-neighbourhood similarity.
 *
 * Reproduces Yu et al. PNAS 2013;110(51):20759 Fig. 2B on this run's data, and
 * reports whether pepM identity could partition BiG-SCAPE's all-pairs problem —
 * the constraint that stops a million-genome run clustering in one pass.
 *
 * Needs the Pfam HMM (for PF13714) and the BiG-SCAPE database, which already
 * holds the neighbourhood-similarity axis for every pair.
 */
process PEPM_ALL_BY_ALL {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}/pepm_all_by_all",
        mode: params.publish_mode

    input:
    val taxon
    path bigscape_db
    path pfam_db

    output:
    path "pepm_vs_neighbourhood.tsv", emit: pairs,    optional: true
    path "pepm_all_by_all.json",      emit: summary,  optional: true
    path "pepm_vs_*.png",             emit: figures,  optional: true
    path "pepm_vs_*.svg",             emit: svgs,     optional: true

    script:
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['analysis/pepm_all_by_all.py', 'utils'])}

    python ${projectDir}/scripts/analysis/pepm_all_by_all.py \\
        --db ${bigscape_db} \\
        --pfam ${pfam_db}/Pfam-A.hmm \\
        --outdir . \\
        --hmmfetch \$(which hmmfetch) \\
        --hmmalign \$(which hmmalign)

    # The intermediate alignment and the Pfam index are large and reproducible;
    # only the analysis outputs are worth publishing.
    rm -rf _pfam pepm.sto pepm.faa PF*.hmm
    """
}
