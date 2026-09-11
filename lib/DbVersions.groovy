/**
 * Reference database version reporting.
 *
 * Every database in this pipeline uses storeDir, which skips its download
 * process whenever the output path already exists. An unpinned URL therefore
 * does not track upstream at all — it freezes whatever was current on the day
 * of the first run and never re-checks. This pipeline classified against a
 * 2026-01-23 NCBI taxdump for seven months without any signal that it had.
 *
 * Pins live in nextflow.config. This class reports them at startup and, when
 * params.check_db_updates is set, says whether upstream has moved past them.
 * It never fails a run and never changes what is downloaded: upgrading a
 * database is a decision with scientific consequences (two runs on different
 * GTDB releases are not directly comparable), so it stays a deliberate edit.
 */
class DbVersions {

    /** Pinned versions, for the log banner and software_versions.json. */
    static Map pinned(params) {
        [
            gtdb_release: params.gtdb_release,
            pfam_release: params.pfam_release,
            taxdump_date: params.taxdump_date,
            // antiSMASH ships its databases keyed to the tool release, so the
            // conda pin on antismash is already the version pin here.
            antismash_db: 'bundled with antismash'
        ]
    }

    /**
     * Latest upstream version for one database, or null when the check cannot
     * be made. Deliberately forgiving: an offline or firewalled machine must
     * still run, so every failure path returns null rather than throwing.
     */
    private static String upstream(String db) {
        try {
            switch (db) {
                case 'gtdb':
                    def txt = fetch('https://data.gtdb.ecogenomic.org/releases/latest/VERSION.txt')
                    def m = txt =~ /v?(\d+)/
                    return m ? m[0][1] : null
                case 'pfam':
                    def txt = fetch('https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt')
                    def m = txt =~ /RELEASE\s+([\d.]+)/
                    return m ? m[0][1] : null
                case 'taxdump':
                    def html = fetch('https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/')
                    def all = (html =~ /taxdmp_(\d{4}-\d{2}-\d{2})\.zip/).collect { it[1] }
                    return all ? all.sort().last() : null
            }
        } catch (Exception ignored) {
            return null
        }
        return null
    }

    private static String fetch(String url) {
        def conn = new URL(url).openConnection()
        conn.setConnectTimeout(5000)
        conn.setReadTimeout(8000)
        return conn.inputStream.getText('UTF-8')
    }

    /**
     * Lines for the startup banner. Each pinned database gets one line, with a
     * note when upstream is ahead. Returns pins alone if checking is disabled.
     */
    static List<String> report(params) {
        def out = []
        def checks = [
            ['GTDB',    'gtdb',    params.gtdb_release as String],
            ['Pfam',    'pfam',    params.pfam_release as String],
            ['taxdump', 'taxdump', params.taxdump_date as String]
        ]
        checks.each { label, key, pin ->
            if (!params.check_db_updates) {
                out << "  ${label.padRight(8)} ${pin}"
                return
            }
            def latest = upstream(key)
            if (latest == null) {
                out << "  ${label.padRight(8)} ${pin}  (upstream check unavailable)"
            } else if (latest != pin) {
                out << "  ${label.padRight(8)} ${pin}  <- upstream now ${latest}; " +
                       "edit params.${key == 'taxdump' ? 'taxdump_date' : key + '_release'} to upgrade"
            } else {
                out << "  ${label.padRight(8)} ${pin}  (current)"
            }
        }
        return out
    }
}
