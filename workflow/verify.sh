#!/usr/bin/env bash
# Compare a fresh annotation directory with the committed annotations/ file by file (md5).
#
#   verify.sh FRESH_ANNOTATIONS_DIR
#
# Passes only if every file listed in the committed MANIFEST.md5 files exists in FRESH with the same md5,
# no fresh assembly directory is absent from the committed set, and at least one file was compared.
set -uo pipefail
F=${1:?fresh annotations dir}; HERE=$(cd "$(dirname "$0")/.." && pwd); A="$HERE/annotations"
[ -s "$A/MANIFEST.md5" ] || { echo "FAIL: $A/MANIFEST.md5 missing or empty (nothing to verify against)" >&2; exit 1; }
[ -d "$F" ] || { echo "FAIL: $F is not a directory" >&2; exit 1; }
[ -f "$F/CALLER_FLAGS.txt" ] && { echo "FAIL: $F was called with non-default flags: $(cat "$F/CALLER_FLAGS.txt")" >&2; exit 1; }
same=0; diff=0; missing=0; extra=0
while read -r md5 path; do
    f="$F/$path"
    if [ ! -f "$f" ]; then missing=$((missing+1)); echo "MISSING $path"; continue; fi
    if [ "$(md5sum < "$f" | cut -d' ' -f1)" = "$md5" ]; then same=$((same+1)); else diff=$((diff+1)); echo "DIFFERS $path"; fi
done < <(cd "$A" && for m in MANIFEST.md5 GCA_*/MANIFEST.md5; do d=$(dirname "$m"); sed "s|  |  $d/|; s|  \./|  |" "$m"; done)
for d in "$F"/GCA_*/; do
    [ -d "$d" ] || continue
    a=$(basename "$d"); [ -s "$A/$a/MANIFEST.md5" ] || { extra=$((extra+1)); echo "EXTRA $a (not in the committed annotations)"; }
done
echo "identical=$same differing=$diff missing=$missing extra=$extra"
[ "$same" -gt 0 ] && [ "$diff" -eq 0 ] && [ "$missing" -eq 0 ] && [ "$extra" -eq 0 ]
