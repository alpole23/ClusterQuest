import java.security.MessageDigest

/**
 * Shared utility functions for the BGC pipeline.
 * These are automatically loaded by Nextflow from the lib/ directory.
 *
 * Usage in modules: Utils.sanitizeTaxon(name)
 * Usage in main.nf: sanitizeTaxon(name) - uses local function definition
 */
class Utils {
    /**
     * Sanitize taxon name for use in file paths and directory names.
     * Replaces special characters with underscores for filesystem safety.
     */
    static String sanitizeTaxon(String name) {
        name.replaceAll('[^a-zA-Z0-9_]', '_').replaceAll('_+', '_').replaceAll('^_|_$', '')
    }

    /**
     * Generate MD5 hash of relevant antiSMASH parameters.
     * Used to detect when parameters change and results need recomputation.
     * Note: clusterhmmer and tigrfam are always enabled, not tracked.
     *
     * @param params The Nextflow params object
     * @return MD5 hash string of the parameter configuration
     */
    static String antismashParamsHash(params) {
        // Collect relevant antiSMASH parameters into a sorted string
        def relevantParams = [
            "minimal=${params.antismash_minimal ?: false}",
            "cb_general=${params.antismash_cb_general ?: false}",
            "cc_mibig=${params.antismash_cc_mibig ?: false}",
            "cb_knownclusters=true",
            "smcog_trees=${params.antismash_smcog_trees ?: false}",
            "hmmdetection_rules=phosphonate"
        ].sort().join(";")

        // Return MD5 hash of the parameter string
        return md5(relevantParams)
    }

    /**
     * Calculate MD5 hash of a string.
     */
    // Memoised across the run: the script block is evaluated once per task, and
    // hashing every .py file for each of thousands of tasks would be wasteful.
    private static final Map<String, String> _scriptsHashCache = [:]

    /**
     * Digest of every Python file under scripts/, so that editing one invalidates
     * the cache of the tasks that run it.
     *
     * Modules invoke their scripts as `python ${projectDir}/scripts/foo.py` — an
     * interpolated path, not a declared `path` input — so Nextflow's task hash does
     * not see the file at all and `-resume` happily reuses output produced by code
     * that no longer exists. Embedding this digest as a comment inside a process's
     * script block puts it into the hashed script text, which fixes that.
     *
     * `paths` are entries under scripts/ — a file, or a package directory hashed
     * recursively. Each process declares only what it actually runs and imports,
     * because a single whole-tree digest would be actively harmful here: it would
     * make RENAME_GENOMES depend on visualisation code, and since RENAME_GENOMES
     * feeds antiSMASH, editing a plotting script would invalidate 1,735 antiSMASH
     * tasks. tests/check_script_deps.py verifies the declared lists still match the
     * scripts' real imports.
     *
     * @param projectDir Pipeline root (the `projectDir` implicit variable)
     * @param paths      Entries under scripts/ (files or package directories)
     * @return 12-character hex digest
     */
    static String scriptsHash(projectDir, List<String> paths) {
        def root = new File("${projectDir}/scripts")
        def key = root.absolutePath + '|' + paths.join(',')
        if (_scriptsHashCache.containsKey(key)) return _scriptsHashCache[key]

        def files = []
        paths.sort().each { rel ->
            def target = new File(root, rel)
            if (target.isDirectory()) {
                target.eachFileRecurse(groovy.io.FileType.FILES) { f ->
                    if (f.name.endsWith('.py')) files << f
                }
            }
            else if (target.isFile()) {
                files << target
            }
        }
        def digest = MessageDigest.getInstance("MD5")
        // sorted so the digest does not depend on filesystem walk order; the relative
        // path is hashed alongside the bytes so renames and deletions register too
        files.unique().sort { root.toPath().relativize(it.toPath()).toString() }.each { f ->
            digest.update(root.toPath().relativize(f.toPath()).toString().getBytes("UTF-8"))
            digest.update(f.bytes)
        }
        def hex = digest.digest().encodeHex().toString().take(12)
        _scriptsHashCache[key] = hex
        return hex
    }

    static String md5(String input) {
        MessageDigest md = MessageDigest.getInstance("MD5")
        byte[] digest = md.digest(input.bytes)
        return digest.collect { String.format("%02x", it) }.join()
    }

    /**
     * Build absolute path for reusing results from a previous taxon run.
     * Handles both absolute and relative outdir paths.
     *
     * @param params The Nextflow params object (needs outdir)
     * @param projectDir The Nextflow projectDir variable
     * @param toolName The tool name (e.g., "antismash", "gtdbtk")
     * @param reuseTaxon The taxon name to reuse results from
     * @param subPath Optional sub-path within the result directory
     * @return Absolute path to the reuse directory
     */
    static String buildReusePath(params, projectDir, String toolName, String reuseTaxon, String subPath = "") {
        def outdir = params.outdir?.toString() ?: "results"
        def outdir_abs = outdir.startsWith('/') ? outdir : "${projectDir}/${outdir}"
        def basePath = "${outdir_abs}/${toolName}_results/${sanitizeTaxon(reuseTaxon)}"
        return subPath ? "${basePath}/${subPath}" : basePath
    }

    /**
     * Check if a file path represents a valid input (not a placeholder).
     * Handles both single files and lists from glob patterns.
     *
     * @param input The file path or list of paths to check
     * @return true if the input is valid, false if it's a placeholder
     */
    static boolean isValidInput(input) {
        if (!input) return false
        if (input instanceof List) {
            return input.size() > 0 && !input[0].name?.startsWith('NO_')
        }
        return !input.name?.startsWith('NO_')
    }

    /**
     * Return the first file from an input that may be a single file or a list.
     *
     * @param input A file path or list of file paths
     * @return The first (or only) file
     */
    static getFirstFile(input) {
        if (input instanceof List) return input[0]
        return input
    }

    /**
     * Build an optional command-line argument for a possibly-placeholder input.
     * Returns "" when the input is a NO_* placeholder, so callers don't have to
     * repeat the sentinel name (and can't get it wrong).
     *
     *   Utils.optArg('--counts', counts_file)  ->  "--counts region_counts.tsv"  or  ""
     *
     * @param flag  The command-line flag, e.g. "--counts"
     * @param input The staged file, list of files, or placeholder
     * @return The flag plus path, or an empty string
     */
    static String optArg(String flag, input) {
        isValidInput(input) ? "${flag} ${getFirstFile(input)}" : ""
    }
}
