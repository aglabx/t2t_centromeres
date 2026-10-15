#!/usr/bin/env bash
# Quick end-to-end check on the smallest assembly (T2T-MFA8v1.1 chrY, GCA_048601705.1, 17 Mb; 1 min 40 s measured):
# NCBI download -> satellome -> call -> GFF with the published MFA8 sets, then md5 against annotations/.
#
#   tests/smoke_test.sh WORKDIR
set -uo pipefail
W=${1:?workdir}; HERE=$(cd "$(dirname "$0")/.." && pwd); ACC=GCA_048601705.1
bash "$HERE/workflow/fetch_assembly.sh" "$ACC" "$W/inputs" || exit 1
bash "$HERE/workflow/run_satellome.sh" "$W/inputs/$ACC/$ACC.fna.gz" "$W/satellome" 10 || exit 1
bash "$HERE/workflow/fetch_official.sh" "$W/official" || exit 1
DN=$(bash "$HERE/workflow/asm_dir.sh" "$ACC") || exit 1
M="$HERE/annotations/$DN/MANIFEST.md5"
[ -s "$M" ] || { echo "FAIL: $M missing or empty (nothing to compare against)" >&2; exit 1; }
bash "$HERE/workflow/call_all.sh" "$W" "$W/official" "$W/annotations" 16 "$ACC" || exit 1
ok=0; bad=0
while read -r md5 f; do
    got=$(md5sum < "$W/annotations/$DN/$f" | cut -d' ' -f1)
    if [ "$got" = "$md5" ]; then ok=$((ok+1)); else bad=$((bad+1)); echo "DIFFERS $DN/$f"; fi
done < "$M"
echo "smoke test: identical=$ok differing=$bad"
[ "$ok" -gt 0 ] && [ "$bad" -eq 0 ]
