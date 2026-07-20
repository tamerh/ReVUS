#!/usr/bin/env bash
# Fetch the two ClinVar variant_summary archives ReVUS is built from.
# The archives live under the NCBI ClinVar FTP tab_delimited/archive tree.
set -euo pipefail

BASE="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive"
OUT="archive"
mkdir -p "$OUT"

# Freeze (2022-06) and horizon (2026-07) releases. Adjust to reproduce other windows.
FREEZE_YEAR=2022; FREEZE="variant_summary_2022-06.txt.gz"
HORIZON_YEAR=2026; HORIZON="variant_summary_2026-07.txt.gz"

# md5 of the exact release files used to build the shipped dataset. NCBI occasionally re-issues
# archived files; verify against these to guarantee byte-for-byte reproduction.
declare -A MD5=(
    ["$FREEZE"]="e282d5e2a42cb46752311a78a05e521d"
    ["$HORIZON"]="7c6b3e8a910e4054f4e57e2c0551c9c8"
)

for pair in "$FREEZE_YEAR/$FREEZE" "$HORIZON_YEAR/$HORIZON"; do
    name="$(basename "$pair")"
    if [ ! -f "$OUT/$name" ]; then
        echo "downloading $name ..."
        curl -fSL "$BASE/$pair" -o "$OUT/$name"
    fi
    got="$(md5sum "$OUT/$name" | awk '{print $1}')"
    if [ "$got" = "${MD5[$name]}" ]; then
        echo "ok   $name  ($got)"
    else
        echo "WARN $name md5 $got != expected ${MD5[$name]} (NCBI may have re-issued this file)"
    fi
done
echo "done -> $OUT/"
