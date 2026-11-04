#!/usr/bin/env bash
# Refuse to start a run when the database volume is not a real mount.
#
# Optional, and specific to a setup that keeps databases on a separate disk. Every
# machine-specific value is a CQ_* override with this workstation's value as the
# default, so it is usable elsewhere rather than being dead weight:
#
#   CQ_MOUNT=/data CQ_LABEL=mydisk CQ_DB_FILES=0 bash scripts/preflight_bigdisk.sh
#
# The WSL `--mount --bare` passthrough and the inner ext4 mount are both lost on
# any WSL shutdown or drive re-plug, and the mountpoint survives as an EMPTY
# DIRECTORY. Nextflow would happily write a fresh work dir onto the internal
# disk and `-resume` would find nothing cached, so this is checked, not assumed.
set -euo pipefail

# Everything specific to one machine is an override, not a constant. The defaults
# are this workstation's external drive; CQ_* lets the same check work elsewhere,
# and CQ_DB_FILES=0 skips the file-count check for a volume whose databases differ.
MNT=${1:-${CQ_MOUNT:-/mnt/bigdisk}}
OUT=${CQ_OUTDIR:-$MNT/clusterquest/results}
WANT_LABEL=${CQ_LABEL:-clusterquest}
WANT_FILES=${CQ_DB_FILES:-140000}
WANT_FREE_GB=${CQ_FREE_GB:-500}

fail() { echo "PREFLIGHT FAIL: $*" >&2; exit 1; }

mountpoint -q "$MNT" || fail "$MNT is not a mountpoint -- re-attach the drive.
  Windows (elevated):  wsl --mount \\\\.\\PHYSICALDRIVE<n> --bare
  WSL:                 sudo mount -L $WANT_LABEL $MNT"

# Mount by label, so a changed device letter cannot silently mount the wrong disk.
if [ -n "$WANT_LABEL" ]; then
    LABEL=$(lsblk -no LABEL "$(findmnt -no SOURCE "$MNT")" 2>/dev/null || echo '')
    [ "$LABEL" = "$WANT_LABEL" ] \
        || fail "$MNT holds label '$LABEL', expected '$WANT_LABEL'"
fi

[ -d "$OUT/databases/gtdbtk" ] || fail "GTDB-Tk database missing under $OUT/databases"

if [ "$WANT_FILES" -gt 0 ]; then
    N=$(find "$OUT/databases" -type f | wc -l)
    [ "$N" -ge "$WANT_FILES" ] \
        || fail "only $N database files, expected >= $WANT_FILES -- incomplete mount?"
else
    N=$(find "$OUT/databases" -type f | wc -l)
fi

AVAIL=$(df -BG --output=avail "$MNT" | tail -1 | tr -dc 0-9)
[ "$AVAIL" -ge "$WANT_FREE_GB" ] \
    || fail "only ${AVAIL} GB free on $MNT, want >= ${WANT_FREE_GB}"

echo "preflight OK: $MNT mounted (label ${WANT_LABEL:-<unchecked>}), $N database files, ${AVAIL} GB free"
