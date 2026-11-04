#!/usr/bin/env bash
# Pre-demo readiness check. Prints one line per prerequisite.
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB=results_demo/databases
ok=0; bad=0
say() { if [ "$1" = ok ]; then printf '  OK       %s\n' "$2"; ok=$((ok+1));
        else printf '  MISSING  %s\n' "$2"; bad=$((bad+1)); fi; }

echo "ClusterQuest demo readiness"
echo

# antiSMASH: mite is what its prereq check dies on first, so test that specifically.
if [ -d "$DB/antismash/mite" ] && [ -d "$DB/antismash/clusterblast" ]; then
    say ok "antiSMASH databases  ($(du -sh $DB/antismash 2>/dev/null | cut -f1))"
else
    say no "antiSMASH databases  ($DB/antismash)"
fi

# BiG-SCAPE needs Pfam-A.hmm AND the pressed index, not just the download.
if [ -f "$DB/pfam38.2/Pfam-A.hmm" ] && [ -f "$DB/pfam38.2/Pfam-A.hmm.h3p" ]; then
    say ok "Pfam 38.2 (pressed)  ($(du -sh $DB/pfam38.2 2>/dev/null | cut -f1))"
else
    say no "Pfam 38.2 (pressed)  -- run demo/finish_pfam.sh"
fi

if [ -f "$DB/taxonkit/taxdump_2026-08-01/nodes.dmp" ]; then
    say ok "NCBI taxdump 2026-08-01"
else
    say no "NCBI taxdump 2026-08-01"
fi

n=$(ls demo/genomes_offline/*.gbff 2>/dev/null | wc -l)
if [ "$n" -eq 24 ]; then say ok "offline genome set ($n genomes)"
else say no "offline genome set (found $n, expected 24)"; fi

if [ -x "$HOME/miniconda3/envs/nextflow/bin/nextflow" ]; then
    say ok "nextflow + JVM"
else
    say no "nextflow + JVM ($HOME/miniconda3/envs/nextflow)"
fi

echo
if [ "$bad" -eq 0 ]; then echo "Ready. Run:  bash demo/run.sh offline"
else echo "$bad prerequisite(s) missing -- see demo/RUNBOOK.md"; exit 1; fi
