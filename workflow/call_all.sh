#!/usr/bin/env bash
# Call centromeric arrays and build the GFF (with every published set from config/official_sources.tsv)
# for every assembly in config/assemblies.tsv; then write the panel tables and md5 manifests.
#
#   call_all.sh WORKDIR OFFICIAL_DIR OUTDIR [THREADS] [ACC,ACC,...]
#     WORKDIR       from satellome_all.sh (inputs/<acc>/ + satellome/<acc>/)
#     OFFICIAL_DIR  from fetch_official.sh (raw/ + norm/)
#     OUTDIR        OUTDIR/<acc>/ annotations + OUTDIR/{SUMMARY,ALL_*,STATUS}.tsv + MANIFEST.md5
set -uo pipefail
W=${1:?workdir}; O=${2:?official dir}; OUT=${3:?outdir}; THREADS=${4:-16}; ONLY=${5:-}
HERE=$(cd "$(dirname "$0")/.." && pwd); SRC="$HERE/src/t2t_centromeres"
mkdir -p "$OUT"
FILES="centromeres.gff3 centromeres.tsv centromeric_arrays.tsv centromeric_arrays.bed official_overlap.tsv family_breadth.tsv outliers.tsv exceptions.tsv"
ST="$OUT/STATUS.tsv"; printf 'accession\tasm_name\tcatalogue\tstatus\tsummary\n' > "$ST"
tail -n +2 "$HERE/config/assemblies.tsv" | while IFS=$'\t' read -r acc asm _; do
    [ -n "$ONLY" ] && [[ ",$ONLY," != *",$acc,"* ]] && continue
    cat_dir="$W/satellome/$acc"; rep="$W/inputs/$acc/$acc.assembly_report.txt"; d="$OUT/$acc"
    if [ ! -f "$cat_dir/.complete" ]; then
        printf '%s\t%s\t\tpending\tno satellome catalogue\n' "$acc" "$asm" >> "$ST"; continue
    fi
    if ! res=$(cd "$SRC" && python call_centromeric_arrays.py --acc "$acc" --asm-name "$asm" --catalogue "$cat_dir" \
                --report "$rep" --out "$d" --threads "$THREADS" 2>&1 | tail -1) || [[ "$res" == *REFUSED* ]]; then
        printf '%s\t%s\t%s\trefused\t%s\n' "$acc" "$asm" "$cat_dir" "$res" >> "$ST"; continue
    fi
    printf '%s\t%s\t%s\tcalled\t%s\n' "$acc" "$asm" "$cat_dir" "$res" >> "$ST"
    OFF=()
    while IFS=$'\t' read -r _ src cls file names regex lengths role; do
        L="-"; [ "$lengths" != "-" ] && L="$HERE/$lengths"
        OFF+=(--official-set "$src|$cls|$O/$file|$names|$regex|$L|$role")
    done < <(awk -F'\t' -v a="$acc" 'NR>1 && $1==a' "$HERE/config/official_sources.tsv")
    (cd "$SRC" && python build_gff.py --calls "$d" --report "$rep" --acc "$acc" "${OFF[@]}" \
        --out-gff "$d/centromeres.gff3" --out-tsv "$d/official_overlap.tsv") || echo "[$acc] GFF FAILED" >&2
    /bin/rm -rf "$d/_mmseqs"
    (cd "$d" && md5sum $FILES > MANIFEST.md5)
done
{ head -1 "$(ls "$OUT"/GCA_*/centromeres.tsv | head -1)"; for f in "$OUT"/GCA_*/centromeres.tsv; do tail -n +2 "$f"; done; } > "$OUT/ALL_centromeres.tsv"
{ head -1 "$(ls "$OUT"/GCA_*/official_overlap.tsv | head -1)"; for f in "$OUT"/GCA_*/official_overlap.tsv; do tail -n +2 "$f"; done; } > "$OUT/ALL_official_overlap.tsv"
python -I "$SRC/summarize.py" "$OUT" "$HERE/config/assemblies.tsv" > "$OUT/SUMMARY.tsv"
(cd "$OUT" && md5sum ALL_centromeres.tsv ALL_official_overlap.tsv SUMMARY.tsv > MANIFEST.md5)
echo "EXIT $(date -Is) called=$(awk -F'\t' '$4=="called"' "$ST" | wc -l)"
