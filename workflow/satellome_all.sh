#!/usr/bin/env bash
# Fetch every assembly in config/assemblies.tsv and run satellome on it, JOBS at a time.
#
#   satellome_all.sh WORKDIR [JOBS] [THREADS]
#     WORKDIR/inputs/<acc>/     NCBI FASTA + assembly report (fetch_assembly.sh)
#     WORKDIR/satellome/<acc>/  catalogue (run_satellome.sh), .complete on success
set -uo pipefail
W=${1:?workdir}; JOBS=${2:-5}; THREADS=${3:-10}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$W/inputs" "$W/satellome" "$W/logs"
tail -n +2 "${ASSEMBLIES:-$HERE/../config/assemblies.tsv}" | cut -f1 > "$W/logs/accessions.txt"
one () {
    acc=$1
    bash "$HERE/fetch_assembly.sh" "$acc" "$W/inputs" >> "$W/logs/$acc.log" 2>&1 || { echo "$acc FETCH_FAILED"; return; }
    bash "$HERE/run_satellome.sh" "$W/inputs/$acc/$acc.fna.gz" "$W/satellome" "$THREADS" >> "$W/logs/$acc.log" 2>&1 \
        && echo "$acc DONE" || echo "$acc SATELLOME_FAILED"
}
export -f one; export W HERE THREADS
xargs -a "$W/logs/accessions.txt" -P "$JOBS" -I{} bash -c 'one {}' | tee "$W/logs/satellome_all.status"
n=$(wc -l < "$W/logs/accessions.txt"); done_n=$(grep -c ' DONE$' "$W/logs/satellome_all.status")
echo "EXIT $(date -Is) done=$done_n/$n"
grep -v ' DONE$' "$W/logs/satellome_all.status"
[ "$done_n" -eq "$n" ]
