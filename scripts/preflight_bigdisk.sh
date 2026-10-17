#!/usr/bin/env bash
# Refuse to start a run when /mnt/bigdisk is not a real mount.
#
# The WSL `--mount --bare` passthrough and the inner ext4 mount are both lost on
# any WSL shutdown or drive re-plug, and the mountpoint survives as an EMPTY
# DIRECTORY. Nextflow would happily write a fresh work dir onto the internal
# disk and `-resume` would find nothing cached, so this is checked, not assumed.
set -euo pipefail
MNT=${1:-/mnt/bigdisk}
OUT="$MNT/clusterquest/results"
fail() { echo "PREFLIGHT FAIL: $*" >&2; exit 1; }

mountpoint -q "$MNT" || fail "$MNT is not a mountpoint -- re-attach the drive.
  Windows (elevated):  wsl --mount \\\\.\\PHYSICALDRIVE<n> --bare
  WSL:                 sudo mount -L clusterquest $MNT"

# Mount by label, so a changed device letter cannot silently mount the wrong disk.
LABEL=$(lsblk -no LABEL "$(findmnt -no SOURCE "$MNT")")
[ "$LABEL" = clusterquest ] || fail "$MNT holds label '$LABEL', expected 'clusterquest'"

[ -d "$OUT/databases/gtdbtk" ] || fail "GTDB-Tk database missing under $OUT/databases"
N=$(find "$OUT/databases" -type f | wc -l)
[ "$N" -ge 140000 ] || fail "only $N database files, expected ~144,028 -- incomplete mount?"

AVAIL=$(df -BG --output=avail "$MNT" | tail -1 | tr -dc 0-9)
[ "$AVAIL" -ge 500 ] || fail "only ${AVAIL} GB free on $MNT"

echo "preflight OK: $MNT mounted (label clusterquest), $N database files, ${AVAIL} GB free"
