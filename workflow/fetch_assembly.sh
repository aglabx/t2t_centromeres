#!/usr/bin/env bash
# Download one assembly (GenBank genomic FASTA + assembly report) from NCBI and verify the FASTA
# against NCBI's md5checksums.txt. Refuses on any mismatch.
#
#   fetch_assembly.sh ACCESSION OUTDIR      -> OUTDIR/ACCESSION/{ACCESSION.fna.gz, ACCESSION.assembly_report.txt}
set -euo pipefail
ACC=${1:?accession}; OUT=${2:?outdir}
HERE=$(cd "$(dirname "$0")/.." && pwd)
AS=${ASSEMBLIES:-$HERE/config/assemblies.tsv}   # tests/smoke_test.sh points this at tests/smoke_assemblies.tsv
FTP=$(awk -F'\t' -v a="$ACC" '$1==a{print $7}' "$AS")
[ -n "$FTP" ] || { echo "REFUSED: $ACC not in $AS" >&2; exit 1; }
BASE=$(basename "$FTP")
D="$OUT/$ACC"; mkdir -p "$D"
# every download says why it failed; a stalled transfer (< 1 kB/s for 2 min) is aborted and retried
get () {
    curl -sSfL --retry 5 --retry-delay 10 --retry-all-errors --speed-limit 1000 --speed-time 120 -o "$1" "$2" \
        || { echo "FAILED $ACC download $2" >&2; return 1; }
}
get "$D/md5checksums.txt" "${FTP}md5checksums.txt"
get "$D/$ACC.assembly_report.txt" "${FTP}${BASE}_assembly_report.txt"
if [ ! -s "$D/$ACC.fna.gz" ]; then
    # unique temp name + atomic rename: a concurrent fetch of the same accession cannot clobber it
    part=$(mktemp "$D/.$ACC.fna.gz.XXXXXX")
    get "$part" "${FTP}${BASE}_genomic.fna.gz" || { /bin/rm -f "$part"; exit 1; }
    chmod 644 "$part"; /bin/mv -f "$part" "$D/$ACC.fna.gz"
fi
want=$(awk -v f="./${BASE}_genomic.fna.gz" '$2==f{print $1}' "$D/md5checksums.txt")
have=$(md5sum < "$D/$ACC.fna.gz" | cut -d' ' -f1)
[ -n "$want" ] && [ "$want" = "$have" ] || { echo "REFUSED: $ACC md5 $have != NCBI $want" >&2; exit 1; }
grep -q "^# Assembly name:" "$D/$ACC.assembly_report.txt" || { echo "REFUSED: $ACC assembly report malformed" >&2; exit 1; }
echo "$ACC ok md5=$have"
