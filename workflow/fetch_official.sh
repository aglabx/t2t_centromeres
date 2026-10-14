#!/usr/bin/env bash
# Download the published centromere annotations listed in config/official_downloads.tsv and verify
# each against its recorded md5; refuses if any file is missing or differs (the comparison would
# otherwise silently run against a different release). Then extracts the two supplementary-table
# sheets that carry coordinates and normalises everything for build_gff.py.
#
#   fetch_official.sh OFFICIAL_DIR SIAMANG_REPORT... -> OFFICIAL_DIR/raw/*, norm/*, LEDGER.tsv
#   (the NCBI assembly reports of the siamang assemblies fill one name the T2T-Browser alias table lacks)
set -uo pipefail
O=${1:?official dir}; shift; HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$O/raw"
L="$O/LEDGER.tsv"; printf 'name\turl\texpected_md5\tobserved_md5\tstatus\n' > "$L"
fail=0
while IFS=$'\t' read -r name url md5 pub; do
    [ "$name" = name ] && continue
    f="$O/raw/$name"
    [ -s "$f" ] || curl -sfL --retry 3 -o "$f" "$url" || /bin/rm -f "$f"
    got=$([ -s "$f" ] && md5sum < "$f" | cut -d' ' -f1)
    if [ "$got" = "$md5" ]; then st=ok; else st=MISMATCH; fail=1; fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$md5" "${got:-missing}" "$st" >> "$L"
done < "$HERE/config/official_downloads.tsv"
column -t -s$'\t' "$L" | cut -c1-200
[ "$fail" -eq 0 ] || { echo "REFUSED: an official file is missing or differs from the recorded release" >&2; exit 1; }
X="$HERE/src/t2t_centromeres/xlsx_sheet_to_tsv.py"
python -I "$X" "$O/raw/T2T-MFA8_Nature2025_SuppTables.xlsx" --sheet "S39. Coord of active cens" > "$O/raw/MFA8_S39.tsv"
python -I "$X" "$O/raw/T2T-sheep_NatGenet2025_SuppTables.xlsx" --sheet "S8" > "$O/raw/sheep_S8.tsv"
python -I "$HERE/src/t2t_centromeres/prepare_official.py" "$O" "$@"
