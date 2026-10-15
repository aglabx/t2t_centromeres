#!/usr/bin/env bash
# Call centromeric arrays and build the GFF (with every published set from config/official_sources.tsv)
# for every assembly in config/assemblies.tsv; then write the panel tables and md5 manifests.
#
#   call_all.sh WORKDIR OFFICIAL_DIR OUTDIR [THREADS] [ACC,ACC,...]
#     WORKDIR       from satellome_all.sh (inputs/<acc>/ + satellome/<acc>/)
#     OFFICIAL_DIR  from fetch_official.sh (raw/ + norm/)
#     OUTDIR        OUTDIR/<Organism>__<assembly>__<acc>/ annotations (workflow/asm_dir.sh)
#                   + OUTDIR/{SUMMARY,ALL_*,STATUS}.tsv + MANIFEST.md5
#   CALLER_FLAGS="--merge-gap 250000" call_all.sh ...   extra call_centromeric_arrays.py flags (sensitivity
#                 runs); recorded in OUTDIR/CALLER_FLAGS.txt, so such a run never passes for the release
set -uo pipefail
W=${1:?workdir}; O=${2:?official dir}; OUT=${3:?outdir}; THREADS=${4:-16}; ONLY=${5:-}
HERE=$(cd "$(dirname "$0")/.." && pwd); SRC="$HERE/src/t2t_centromeres"
mkdir -p "$OUT"
# the python steps run from $SRC, so every path handed to them must be absolute
W=$(cd "$W" && pwd) && O=$(cd "$O" && pwd) && OUT=$(cd "$OUT" && pwd) || exit 1
read -r -a XF <<< "${CALLER_FLAGS:-}"
if [ ${#XF[@]} -gt 0 ]; then printf '%s\n' "${XF[*]}" > "$OUT/CALLER_FLAGS.txt"; else /bin/rm -f "$OUT/CALLER_FLAGS.txt"; fi
FILES="centromeres.gff3 centromeres.tsv centromeric_arrays.tsv centromeric_arrays.bed official_overlap.tsv family_breadth.tsv outliers.tsv exceptions.tsv satellite_arrays.tsv karyotype.svg"
ST="$OUT/STATUS.tsv"; printf 'accession\tasm_name\tdirectory\tcatalogue\tstatus\tsummary\n' > "$ST"
tail -n +2 "$HERE/config/assemblies.tsv" | while IFS=$'\t' read -r acc asm org _; do
    [ -n "$ONLY" ] && [[ ",$ONLY," != *",$acc,"* ]] && continue
    dn=$(bash "$HERE/workflow/asm_dir.sh" "$acc") || { printf '%s\t%s\t\t\trefused\tno directory name\n' "$acc" "$asm" >> "$ST"; continue; }
    # STATUS.tsv names the catalogue relative to WORKDIR: the table is published, the workdir is local
    cat_rel="satellome/$acc"; cat_dir="$W/$cat_rel"; rep="$W/inputs/$acc/$acc.assembly_report.txt"; d="$OUT/$dn"
    if [ ! -f "$cat_dir/.complete" ]; then
        printf '%s\t%s\t%s\t\tpending\tno satellome catalogue\n' "$acc" "$asm" "$dn" >> "$ST"; continue
    fi
    if ! res=$(cd "$SRC" && python call_centromeric_arrays.py --acc "$acc" --asm-name "$asm" --catalogue "$cat_dir" \
                --report "$rep" --out "$d" --threads "$THREADS" ${XF[@]+"${XF[@]}"} 2>&1 | tail -1) || [[ "$res" == *REFUSED* ]]; then
        printf '%s\t%s\t%s\t%s\trefused\t%s\n' "$acc" "$asm" "$dn" "$cat_rel" "$res" >> "$ST"; continue
    fi
    OFF=()
    while IFS=$'\t' read -r _ src cls file names regex lengths role; do
        L="-"; [ "$lengths" != "-" ] && L="$HERE/$lengths"
        OFF+=(--official-set "$src|$cls|$O/$file|$names|$regex|$L|$role")
    done < <(awk -F'\t' -v a="$acc" 'NR>1 && $1==a' "$HERE/config/official_sources.tsv")
    if ! gff=$(cd "$SRC" && python build_gff.py --calls "$d" --report "$rep" --acc "$acc" "${OFF[@]}" \
            --out-gff "$d/centromeres.gff3" --out-tsv "$d/official_overlap.tsv" 2>&1); then
        printf '%s\t%s\t%s\t%s\tgff_failed\t%s\n' "$acc" "$asm" "$dn" "$cat_rel" "$(printf '%s' "$gff" | tail -1)" >> "$ST"; continue
    fi
    if ! kar=$(cd "$SRC" && python plot_karyotype.py --calls "$d" --report "$rep" --acc "$acc" --asm-name "$asm" \
            --organism "$org" --out "$d/karyotype.svg" 2>&1); then
        printf '%s\t%s\t%s\t%s\tplot_failed\t%s\n' "$acc" "$asm" "$dn" "$cat_rel" "$(printf '%s' "$kar" | tail -1)" >> "$ST"; continue
    fi
    /bin/rm -rf "$d/_mmseqs"
    (cd "$d" && md5sum $FILES > MANIFEST.md5) || { printf '%s\t%s\t%s\t%s\tmanifest_failed\t%s\n' "$acc" "$asm" "$dn" "$cat_rel" "$res" >> "$ST"; continue; }
    printf '%s\t%s\t%s\t%s\tcalled\t%s\n' "$acc" "$asm" "$dn" "$cat_rel" "$res" >> "$ST"
done
{ head -1 "$(ls "$OUT"/*__GCA_*/centromeres.tsv | head -1)"; for f in "$OUT"/*__GCA_*/centromeres.tsv; do tail -n +2 "$f"; done; } > "$OUT/ALL_centromeres.tsv"
{ head -1 "$(ls "$OUT"/*__GCA_*/official_overlap.tsv | head -1)"; for f in "$OUT"/*__GCA_*/official_overlap.tsv; do tail -n +2 "$f"; done; } > "$OUT/ALL_official_overlap.tsv"
python -I "$SRC/summarize.py" "$OUT" "$HERE/config/assemblies.tsv" > "$OUT/SUMMARY.tsv"
(cd "$OUT" && md5sum ALL_centromeres.tsv ALL_official_overlap.tsv SUMMARY.tsv > MANIFEST.md5)
called=$(awk -F'\t' 'NR>1 && $5=="called"' "$ST" | wc -l)
failed=$(awk -F'\t' 'NR>1 && $5!="called" && $5!="pending"' "$ST" | wc -l)
echo "EXIT $(date -Is) called=$called failed=$failed (STATUS.tsv)${XF[@]+ caller_flags: ${XF[*]}}"
[ "$failed" -eq 0 ]
