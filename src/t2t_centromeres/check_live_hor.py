"""External check: does our centromeric locus, and
the largest array inside it, overlap the published live/active alpha HOR of the same chromosome?
Also scores the null method (largest array of the chromosome, same catalogue). The published
track is a check only, never an input to the call.

  CHM13:  --names ucsc     --live 'hor_.*L\\)'      (chm13 censat v2.x)
  apes:   --names seqname  --live '^active_hor'    (NHGRI ape CenSat, cur.20231122 coordinates)

A chromosome is compared only when the track's assembly has the same length for it as the
assembly we called on (--track-fai); otherwise it is reported as `not_comparable`, not scored.

Usage: check_chm13_live_hor.py --calls CENTROMERES_TSV --track BED --report ASSEMBLY_REPORT
           --catalogue RUN_DIR --names ucsc|seqname --live REGEX [--track-fai FAI] --out TSV
"""
import argparse
import re
from collections import defaultdict

from catalogue import read_bed

ap = argparse.ArgumentParser()
for a in ("--calls", "--track", "--report", "--catalogue", "--out"):
    ap.add_argument(a, required=True)
ap.add_argument("--names", choices=["ucsc", "seqname"], required=True)
ap.add_argument("--live", required=True)
ap.add_argument("--track-fai", default=None)
args = ap.parse_args()
live_re = re.compile(args.live)

track_name, length = {}, {}
with open(args.report) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if f[1] == "assembled-molecule" and f[3] == "Chromosome":
            t = f[9] if args.names == "ucsc" else f[0]
            for key in (f[0], f[4], f[6]):
                track_name[key] = t
            length[t] = int(f[8])

track_len = {}
if args.track_fai:
    with open(args.track_fai) as fh:
        for ln in fh:
            f = ln.split("\t")
            track_len[f[0]] = int(f[1])

live = defaultdict(list)
with open(args.track) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) > 3 and live_re.search(f[3]):
            live[f[0]].append((int(f[1]), int(f[2]), f[3]))


def hits(chrom, s, e):
    return [n for a, b, n in live[chrom] if a < e and s < b]


largest = {}
for seq, s, e, p, _ in read_bed(args.catalogue):
    if seq in track_name and (seq not in largest or e - s > largest[seq][1] - largest[seq][0]):
        largest[seq] = (s, e, p)

rows = []
score = defaultdict(int)
with open(args.calls) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    c = {k: i for i, k in enumerate(head)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        u = track_name[f[c["chromosome"]]]
        if args.track_fai and track_len.get(u) != length[u]:
            status = "not_comparable_length" if u in track_len else "not_in_track_assembly"
        elif not live[u]:
            status = "no_live_hor_in_track"
        else:
            status = "compared"
        loc = hits(u, int(f[c["locus_start"]]), int(f[c["locus_end"]]))
        big = hits(u, int(f[c["largest_start"]]), int(f[c["largest_end"]]))
        ns, ne, npd = largest[f[c["seq_name"]]]
        null = hits(u, ns, ne)
        score[status] += 1
        if status == "compared":
            score["locus"] += bool(loc)
            score["largest"] += bool(big)
            score["null"] += bool(null)
        rows.append([u, status, f[c["locus_start"]], f[c["locus_end"]], int(bool(loc)), int(bool(big)),
                     int(bool(null)), ns, ne, npd, ";".join(loc) or "-",
                     ";".join(f"{a}-{b}" for a, b, _ in live[u]) or "-"])
with open(args.out, "w") as out:
    out.write("chrom\tstatus\tlocus_start\tlocus_end\tlocus_on_live_hor\tlargest_array_on_live_hor\t"
              "null_largest_array_on_live_hor\tnull_start\tnull_end\tnull_period\tlive_hor_hit\t"
              "live_hor_intervals\n")
    for r in rows:
        out.write("\t".join(str(x) for x in r) + "\n")
print(f"chromosomes={len(rows)} compared={score['compared']} locus_on_live_hor={score['locus']} "
      f"largest_array_in_locus_on_live_hor={score['largest']} null_largest_array_on_live_hor={score['null']} "
      f"not_comparable={len(rows) - score['compared']}")
