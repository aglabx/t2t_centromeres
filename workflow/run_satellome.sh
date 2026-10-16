#!/usr/bin/env bash
# Run satellome (FasTAN engine) on one assembly and accept the catalogue only after content checks.
#
#   run_satellome.sh ASSEMBLY OUTDIR [THREADS]
#
# ASSEMBLY  .fna / .fna.gz / .fa(.gz) / .2bit (2bit needs twoBitToFa on PATH)
# OUTDIR    the catalogue lands in OUTDIR/<accession>/; OUTDIR/<accession>/.complete marks success
#
# Why each step exists:
# - FasTAN truncates the input basename at the first dot, so the input is staged as a dot-free name
#   (GCA_009914755.4 -> GCA_009914755v4); output files carry that prefix.
# - arraysplitter ignores -t for rayon's global pool; RAYON_NUM_THREADS bounds it.
# - satellome's exit code and results.yaml are recorded but not trusted: success means the FasTAN BED
#   and the arraysplitter summary exist, are non-empty and the summary's last row is complete.
set -uo pipefail
IN=${1:?assembly}; OUT=${2:?outdir}; THREADS=${3:-10}
ACC=$(basename "$IN"); ACC=${ACC%.gz}; for s in .2bit .fasta .fna .fa; do ACC=${ACC%$s}; done; ACC=${ACC%_genomic}
SAFE=$(printf '%s' "$ACC" | tr '.' 'v')
D="$OUT/$ACC"
[ -f "$D/.complete" ] && { echo "$ACC already complete"; exit 0; }
/bin/rm -rf "$D"; mkdir -p "$D"
FNA="$D/$SAFE.fna"
case "$IN" in
    *.2bit) twoBitToFa "$IN" "$FNA" ;;
    *.gz)   gzip -cd "$IN" > "$FNA" ;;
    *)      ln -sf "$(readlink -f "$IN")" "$FNA" ;;
esac || { echo "FAIL $ACC decode" >&2; exit 1; }
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-1}
python -m satellome --version > "$D/_satellome.version" 2>&1
python -m satellome -i "$FNA" -o "$D" -p "$SAFE" -t "$THREADS" > "$D/_satellome.log" 2>&1
RC=$?
# arraysplitter can outlive its parent; wait for writers to settle (sizes unchanged 5 s apart, <= 24 h)
sig () { find "$D" -type f -printf '%s ' | md5sum | cut -d' ' -f1; }
prev=""; for _ in $(seq 1 17280); do cur=$(sig); [ "$cur" = "$prev" ] && break; prev=$cur; sleep 5; done
/bin/rm -f "$FNA"
bed=$(ls "$D"/fastan/*.bed 2>/dev/null | head -1)
sum=$(ls "$D"/fastan/*.summary.tsv* 2>/dev/null | head -1)
[ -s "$bed" ] && [ -s "$sum" ] || { echo "FAIL $ACC no_catalogue rc=$RC" >&2; exit 1; }
read -r h l n < <( (case "$sum" in *.gz) gzip -cd "$sum";; *) cat "$sum";; esac) | awk -F'\t' 'NR==1{h=NF}{l=NF} END{print h+0, l+0, NR+0}')
[ "$n" -ge 2 ] && [ "$h" = "$l" ] || { echo "FAIL $ACC truncated_summary cols=$h/$l rows=$n" >&2; exit 1; }
/bin/rm -rf "$D/fasta" "$D/images"
for f in "$D"/*.sat "$D"/fastan/*.tsv "$D"/fastan/*.fasta; do [ -f "$f" ] && gzip -f "$f"; done
# satellome writes gaps.bed and the telomere check after ArraySplitter; without them it stopped early (its
# one-hour ArraySplitter limit under load), and the arrays may be unsorted: refuse, rerun with more threads
ls "$D"/*.gaps.bed > /dev/null 2>&1 && [ -s "$D/reports/telomeres.tsv" ] \
    || { echo "FAIL $ACC incomplete_run rc=$RC (no gaps.bed or reports/telomeres.tsv; rerun alone with more threads)" >&2; exit 1; }
[ -s "$D/results.yaml" ] && y=yes || y=no_results_yaml
printf 'rc=%s results_yaml=%s\n' "$RC" "$y" > "$D/.complete"
echo "$ACC complete rc=$RC results_yaml=$y"
