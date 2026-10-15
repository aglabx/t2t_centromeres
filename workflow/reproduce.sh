#!/usr/bin/env bash
# End-to-end reproduction of annotations/ from public inputs.
#
#   reproduce.sh WORKDIR [JOBS] [THREADS]
#
# 1. satellome_all.sh   NCBI GenBank FASTA + assembly report per assembly -> satellome (FasTAN) catalogue
# 2. fetch_official.sh  published centromere annotations, md5-checked against config/official_downloads.tsv
# 3. call_all.sh        our call + GFF with the published sets -> WORKDIR/annotations
# 4. verify.sh          compares WORKDIR/annotations with the committed annotations/ by md5
#
# Needs the environment in environment.yml (satellome 1.17.2 with its managed FasTAN/tanbed,
# arraysplitter, BLAST+). satellome takes ~35-40 min wall per mammalian assembly at 10 threads.
set -uo pipefail
W=${1:?workdir}; JOBS=${2:-5}; THREADS=${3:-10}
HERE=$(cd "$(dirname "$0")" && pwd)
bash "$HERE/satellome_all.sh" "$W" "$JOBS" "$THREADS" \
    || echo "WARNING: not every assembly has a catalogue (see $W/logs/satellome_all.status); they will show as MISSING in verify" >&2
reports=$(ls "$W"/inputs/GCA_028878055.3/*.assembly_report.txt "$W"/inputs/GCA_028878085.3/*.assembly_report.txt)
bash "$HERE/fetch_official.sh" "$W/official" $reports || exit 1
bash "$HERE/call_all.sh" "$W" "$W/official" "$W/annotations" 16   # the thread count the published run used
bash "$HERE/verify.sh" "$W/annotations"
