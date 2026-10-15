#!/usr/bin/env bash
# Print the annotation directory of an accession: <Organism>__<assembly name>__<accession>, spaces -> _
# (Homo_sapiens__T2T-CHM13v2.0__GCA_009914755.4), from config/assemblies.tsv.
#
#   asm_dir.sh ACC
set -uo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
n=$(awk -F'\t' -v a="${1:?accession}" 'NR>1 && $1==a{o=$3; gsub(/ /,"_",o); print o "__" $2 "__" $1}' "$HERE/config/assemblies.tsv")
[ -n "$n" ] || { echo "REFUSED: $1 not in config/assemblies.tsv" >&2; exit 1; }
[ "$(printf '%s\n' "$n" | wc -l)" -eq 1 ] || { echo "REFUSED: $1 listed more than once in config/assemblies.tsv" >&2; exit 1; }
case "$n" in */*|*[[:space:]]*) echo "REFUSED: directory name '$n' contains / or whitespace" >&2; exit 1;; esac
echo "$n"
