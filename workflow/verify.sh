#!/usr/bin/env bash
# Compare a fresh annotation directory with the committed annotations/ file by file (md5).
#
#   verify.sh FRESH_ANNOTATIONS_DIR
set -uo pipefail
F=${1:?fresh annotations dir}; HERE=$(cd "$(dirname "$0")/.." && pwd); A="$HERE/annotations"
same=0; diff=0; missing=0
while read -r md5 path; do
    f="$F/$path"
    if [ ! -f "$f" ]; then missing=$((missing+1)); echo "MISSING $path"; continue; fi
    if [ "$(md5sum < "$f" | cut -d' ' -f1)" = "$md5" ]; then same=$((same+1)); else diff=$((diff+1)); echo "DIFFERS $path"; fi
done < <(cd "$A" && for m in MANIFEST.md5 GCA_*/MANIFEST.md5; do d=$(dirname "$m"); sed "s|  |  $d/|; s|  \./|  |" "$m"; done)
echo "identical=$same differing=$diff missing=$missing"
[ "$diff" -eq 0 ] && [ "$missing" -eq 0 ]
