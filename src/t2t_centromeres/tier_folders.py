"""Move every called assembly folder into the folder of its completeness tier and name the tiered path in STATUS.tsv.

Usage: python -I tier_folders.py OUTDIR      (after tier.py has written OUTDIR/TIERS.tsv)

OUTDIR/<Organism>__<assembly>__<acc>/ -> OUTDIR/<tier>/<Organism>__<assembly>__<acc>/, tier one of real_t2t,
gap_free, gapped. A copy of the same folder left in any tier folder by an earlier run is removed first, so an
assembly whose tier changed does not stay in two places. Rows of STATUS.tsv that were not called keep their flat
directory name and nothing of theirs is moved. STATUS.tsv is rewritten line by line with only its directory field
changed (no CSV quoting). Refuses (exit 1) when a called assembly has no tier, a tier is not one of the three, or a
called folder is missing.
"""
import csv
import os
import shutil
import sys

TIERS = ("real_t2t", "gap_free", "gapped")
out = sys.argv[1]
st_path = os.path.join(out, "STATUS.tsv")
with open(st_path) as fh:
    lines = fh.read().split("\n")
if lines[-1] == "":
    lines.pop()
head = lines[0].split("\t")
i_acc, i_dir, i_st = head.index("accession"), head.index("directory"), head.index("status")
with open(os.path.join(out, "TIERS.tsv")) as fh:
    tier = {r["accession"]: r["completeness_tier"] for r in csv.DictReader(fh, delimiter="\t")}

# check every row before moving anything, so a refusal leaves OUTDIR as it was
rows, moves = [line.split("\t") for line in lines[1:]], []
for f in rows:
    if f[i_st] != "called":
        continue
    acc, dn = f[i_acc], f[i_dir]
    t = tier.get(acc)
    if t not in TIERS:
        sys.exit(f"tier_folders REFUSED: {acc} has tier {t!r}, expected one of {TIERS}")
    if "/" in dn:
        sys.exit(f"tier_folders REFUSED: {acc} directory {dn!r} is already tiered (STATUS.tsv rewritten twice?)")
    if not os.path.isdir(os.path.join(out, dn)):
        sys.exit(f"tier_folders REFUSED: {acc} folder {os.path.join(out, dn)} missing")
    moves.append((f, dn, str(t)))
for f, dn, t in moves:
    for old in TIERS:
        stale = os.path.join(out, old, dn)
        if os.path.lexists(stale):
            shutil.rmtree(stale)
    os.makedirs(os.path.join(out, t), exist_ok=True)
    os.rename(os.path.join(out, dn), os.path.join(out, t, dn))
    f[i_dir] = f"{t}/{dn}"
new = [lines[0]] + ["\t".join(f) for f in rows]
with open(st_path, "w") as fh:
    fh.write("\n".join(new) + "\n")
