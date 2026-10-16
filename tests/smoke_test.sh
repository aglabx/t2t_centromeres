#!/usr/bin/env bash
# Quick end-to-end check on a small input (T2T-MFA8v1.1 chrY deposited on its own, GCA_048601705.1, 17 Mb;
# about 3 min from an empty directory): NCBI download -> satellome -> call -> GFF with the published MFA8 sets -> karyotype, then md5
# against tests/expected/. This accession is a test fixture, not part of the panel: its one sequence (CP168063.1)
# is chrY of GCA_037993035.2, which the panel already contains.
#
#   tests/smoke_test.sh WORKDIR
set -uo pipefail
W=${1:?workdir}; HERE=$(cd "$(dirname "$0")/.." && pwd); ACC=GCA_048601705.1
export ASSEMBLIES="$HERE/tests/smoke_assemblies.tsv" OFFICIAL_SOURCES="$HERE/tests/smoke_official_sources.tsv"
bash "$HERE/workflow/fetch_assembly.sh" "$ACC" "$W/inputs" || exit 1
bash "$HERE/workflow/run_satellome.sh" "$W/inputs/$ACC/$ACC.fna.gz" "$W/satellome" 10 || exit 1
bash "$HERE/workflow/fetch_official.sh" "$W/official" || exit 1
DN=$(bash "$HERE/workflow/asm_dir.sh" "$ACC") || exit 1
M="$HERE/tests/expected/$DN/MANIFEST.md5"
[ -s "$M" ] || { echo "FAIL: $M missing or empty (nothing to compare against)" >&2; exit 1; }
bash "$HERE/workflow/call_all.sh" "$W" "$W/official" "$W/annotations" 16 "$ACC" || exit 1
TD=$(ls -d "$W"/annotations/*/"$DN") && [ "$(printf '%s\n' "$TD" | wc -l)" -eq 1 ] \
    || { echo "FAIL: $DN is not in exactly one tier folder of $W/annotations" >&2; exit 1; }
ok=0; bad=0
while read -r md5 f; do
    got=$(md5sum < "$TD/$f" | cut -d' ' -f1)
    if [ "$got" = "$md5" ]; then ok=$((ok+1)); else bad=$((bad+1)); echo "DIFFERS $DN/$f"; fi
done < "$M"
echo "smoke test: identical=$ok differing=$bad"
[ "$ok" -gt 0 ] && [ "$bad" -eq 0 ]
