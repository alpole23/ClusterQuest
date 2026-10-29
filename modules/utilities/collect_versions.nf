process COLLECT_VERSIONS {
    tag "versions"
    label 'process_low'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode

    input:
    path antismash_done  // Dependency to ensure conda envs exist
    path bigscape_done   // Optional dependency for BiG-SCAPE
    path gtdbtk_done     // Optional dependency for GTDB-Tk

    output:
    path "software_versions.json", emit: versions

    script:
    def nf_version = workflow.nextflow.version
    // Every place a conda environment might be, newest first. This used to be
    // `${workDir}/conda` alone, which stopped being true the day conda.cacheDir
    // moved to ~/.clusterquest/conda so environments would survive a change of
    // work directory. Nothing failed -- find_executable simply found nothing and
    // every tool was recorded as "not installed", in the report header and in
    // software_versions.json, so runs from that period cannot say what they were
    // produced against. Searching both costs one glob over a missing directory.
    def conda_dirs = [
        (System.getenv('NXF_CONDA_CACHEDIR') ?: "${System.getProperty('user.home')}/.clusterquest/conda"),
        "${workDir}/conda",
    ].unique().join(java.io.File.pathSeparator)
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import os
    import re
    from pathlib import Path

    def find_executable(name, search_dirs):
        '''Find an executable in any conda environment, then on PATH.'''
        import shutil
        for search_dir in str(search_dirs).split(os.pathsep):
            if not search_dir:
                continue
            for env_dir in sorted(Path(search_dir).glob("env-*")):
                bin_path = env_dir / "bin" / name
                if bin_path.exists():
                    return str(bin_path)
        # Last resort: the tool is on this task's own PATH. Rare, since each
        # process gets its own environment, but free to check and it is the
        # case for a user running the pipeline inside one combined env.
        return shutil.which(name)
        return None

    def get_version(executable, pattern=None, version_arg="--version"):
        '''Run version command and extract version.'''
        try:
            result = subprocess.run(
                [executable, version_arg],
                capture_output=True, text=True, timeout=30
            )
            output = result.stdout + result.stderr
            if pattern:
                match = re.search(pattern, output)
                if match:
                    return match.group(1)
            # Return first non-empty line
            for line in output.strip().split("\\n"):
                if line.strip():
                    return line.strip()
            return "unknown"
        except Exception as e:
            return f"error: {str(e)}"

    conda_dirs = "${conda_dirs}"
    versions = {}

    # antiSMASH
    antismash_bin = find_executable("antismash", conda_dirs)
    if antismash_bin:
        versions["antismash"] = get_version(antismash_bin, r"antiSMASH\\s+([\\d.]+)")
    else:
        versions["antismash"] = "not installed"

    # BiG-SCAPE - check if it was run (bigscape_done is not a placeholder)
    bigscape_ran = os.path.isdir("${bigscape_done}") or (os.path.isfile("${bigscape_done}") and "NO_BIGSCAPE" not in "${bigscape_done}")
    bigscape_bin = find_executable("bigscape", conda_dirs)
    if bigscape_bin:
        versions["bigscape"] = get_version(bigscape_bin, r"BiG-SCAPE\\s+([\\d.]+)")
    elif bigscape_ran:
        versions["bigscape"] = "used (version unknown)"
    else:
        versions["bigscape"] = "not used"

    # GTDB-Tk - check if it was run (gtdbtk_done is not a placeholder)
    gtdbtk_ran = os.path.isdir("${gtdbtk_done}") or (os.path.isfile("${gtdbtk_done}") and "NO_GTDBTK" not in "${gtdbtk_done}")
    gtdbtk_bin = find_executable("gtdbtk", conda_dirs)
    if gtdbtk_bin:
        versions["gtdbtk"] = get_version(gtdbtk_bin, r"gtdbtk:\\s*v?([\\d.]+)")
    elif gtdbtk_ran:
        versions["gtdbtk"] = "used (version unknown)"
    else:
        versions["gtdbtk"] = "not used"

    # TaxonKit (uses 'version' subcommand, not '--version' flag)
    taxonkit_bin = find_executable("taxonkit", conda_dirs)
    if taxonkit_bin:
        versions["taxonkit"] = get_version(taxonkit_bin, r"taxonkit\\s+v?([\\d.]+)", "version")
    else:
        versions["taxonkit"] = "not installed"

    # Nextflow and pipeline
    versions["nextflow"] = "${nf_version}"
    versions["pipeline_version"] = "1.0.0"

    # Which reference data the results were produced against. Tool versions alone
    # do not identify a run: gtdbtk 2.6.1 against GTDB r226 and against r232 give
    # different taxonomies. storeDir downloads each database once and never
    # re-checks it, so without recording the pins here the version is unknowable
    # after the fact.
    versions["db_gtdb_release"] = "${params.gtdb_release}"
    versions["db_pfam_release"] = "${params.pfam_release}"
    versions["db_taxdump_date"] = "${params.taxdump_date}"

    with open("software_versions.json", "w") as f:
        json.dump(versions, f, indent=2)

    print("Software versions collected:")
    for tool, version in versions.items():
        print(f"  {tool}: {version}")
    """
}
