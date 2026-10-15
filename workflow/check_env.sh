#!/usr/bin/env bash
# Print the tool versions the pipeline will use and compare them with the ones that produced
# annotations/ (config/toolchain.tsv). Exits non-zero on any mismatch.
set -uo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
B=$(python -c "import satellome,os;print(os.path.join(os.path.dirname(satellome.__file__),'bin'))")
sha () { python -c "import json,sys;print(json.load(open(sys.argv[1]))['git_sha'][:7])" "$B/$1.manifest.json" 2>/dev/null || echo unknown; }
declare -A have=(
  [satellome]=$(python -m satellome --version 2>&1 | tail -1 | awk '{print $NF}')
  [arraysplitter]=$(arraysplitter --version 2>&1 | awk '{print $NF}')
  [blastn]=$(blastn -version 2>&1 | head -1 | awk '{print $2}')
  [fastan_git]=$(sha fastan)
  [tanbed_git]=$(sha tanbed)
  [python]=$(python -c "import sys;print('%d.%d' % sys.version_info[:2])")
  [numpy]=$(python -c "import numpy;print(numpy.__version__)" 2>/dev/null || echo missing)
  [zstandard]=$(python -c "import zstandard;print(zstandard.__version__)" 2>/dev/null || echo missing)
  [pyyaml]=$(python -c "import yaml;print(yaml.__version__)" 2>/dev/null || echo missing)
)
bad=0
while IFS=$'\t' read -r tool want note; do
    [ "$tool" = tool ] && continue
    got=${have[$tool]:-missing}
    if [ "$got" = "$want" ]; then s=ok; else s=MISMATCH; bad=1; fi
    printf '%-14s want %-12s have %-12s %s\n' "$tool" "$want" "$got" "$s"
done < "$HERE/config/toolchain.tsv"
[ "$bad" -eq 0 ] || { echo "toolchain differs from the one that produced annotations/" >&2; exit 1; }
