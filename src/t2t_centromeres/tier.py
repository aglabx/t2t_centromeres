"""Completeness tier of every assembly, from its sequence and report, never from its name or NCBI level.

  real_t2t   no gap, no unplaced or unlocalized sequence, telomeric repeat at both ends of every chromosome
  gap_free   no gap, no unplaced or unlocalized sequence, telomeric repeat missing at one or more chromosome ends
  gapped     at least one gap (run of N) or one unplaced / unlocalized sequence

Gaps: satellome `<prefix>.gaps.bed` (every non-header line is one run of N). Unplaced: assembly-report rows whose
role is not `assembled-molecule`, organelles excluded (a mitochondrion deposited as an unlocalized scaffold is not a
piece of a nuclear chromosome). Chromosomes: `assembled-molecule` rows of type Chromosome or Linkage Group (some fish
reports use the latter for whole chromosomes), as in the caller.
Telomeres: satellome `reports/telomeres.tsv` (telomere-check: TTAGGG / CCCTAA in the terminal 10 kb, density
>= 0.15), status PRESENT at both ends. A chromosome is looked up by its GenBank accession or, for a catalogue made
from the RefSeq copy of the assembly, its RefSeq accession. Refuses (exit 1) when a file or a chromosome row is
missing, or a chromosome has rows under both names.

Usage: python -I tier.py WORKDIR STATUS.tsv > TIERS.tsv
  WORKDIR/inputs/<acc>/<acc>.assembly_report.txt and WORKDIR/satellome/<acc>/ as written by satellome_all.sh;
  STATUS.tsv of call_all.sh (every row with status "called", in its order)
"""
import csv
import glob
import os
import sys

work, manifest = sys.argv[1], sys.argv[2]
ORGANELLE = {"Mitochondrion", "Chloroplast", "Plastid", "Apicoplast", "Kinetoplast"}


def table(path):
    """Rows of a tab-separated file whose header may be a '# '-prefixed comment line (NCBI reports)."""
    with open(path, newline="") as fh:
        lines = [ln.rstrip("\r\n") for ln in fh]
    head = [i for i, ln in enumerate(lines) if not ln.startswith("#") or ln.startswith("# Sequence-Name")]
    if not head:
        sys.exit(f"tier REFUSED: no header in {path}")
    body = [lines[i] for i in head]
    if body[0].startswith("# "):
        body[0] = body[0][2:]
    return list(csv.DictReader(body, delimiter="\t"))


cols = ["accession", "chromosomes", "gaps", "gaps_on_chromosomes", "unplaced_sequences", "both_telomeres",
        "one_telomere", "no_telomere", "completeness_tier"]
print("\t".join(cols))
errors = []
with open(manifest) as fh:
    accs = [r["accession"] for r in csv.DictReader(fh, delimiter="\t") if r.get("status", "called") == "called"]
for acc in accs:
    rep_path = os.path.join(work, "inputs", acc, f"{acc}.assembly_report.txt")
    sat = os.path.join(work, "satellome", acc)
    gaps_beds = glob.glob(os.path.join(sat, "*.gaps.bed"))
    tel_path = os.path.join(sat, "reports", "telomeres.tsv")
    missing = [p for p in (rep_path, tel_path) if not os.path.isfile(p)]
    if missing or len(gaps_beds) != 1:
        errors.append(f"{acc}: missing {missing or ''} gaps.bed found {len(gaps_beds)}")
        continue
    rep = table(rep_path)
    chroms = [r for r in rep if r["Sequence-Role"] == "assembled-molecule"
              and r["Assigned-Molecule-Location/Type"] in ("Chromosome", "Linkage Group")]
    unplaced = [r for r in rep if r["Sequence-Role"] != "assembled-molecule"
                and r["Assigned-Molecule-Location/Type"] not in ORGANELLE]
    with open(gaps_beds[0]) as fh:
        gaps = [ln.split("\t")[0] for ln in fh if ln.strip() and not ln.startswith("#")]
    # the catalogue names a chromosome by its GenBank or, when the sequence came from the RefSeq copy, its RefSeq
    # accession (as the caller accepts)
    names = {c["GenBank-Accn"]: [n for n in (c["GenBank-Accn"], c["RefSeq-Accn"]) if n and n != "na"] for c in chroms}
    chrom_ids = {n for ns in names.values() for n in ns}
    tel = {r["chr"]: r for r in table(tel_path)}
    both = one = none = 0
    for c in chroms:
        rows = [tel[n] for n in names[c["GenBank-Accn"]] if n in tel]
        if len(rows) != 1:
            errors.append(f"{acc}: chromosome {c['GenBank-Accn']} has "
                          + ("no telomere row" if not rows else "a telomere row under both its GenBank and RefSeq names"))
            continue
        t = rows[0]
        n = (t["left_status"] == "PRESENT") + (t["right_status"] == "PRESENT")
        both, one, none = both + (n == 2), one + (n == 1), none + (n == 0)
    if not chroms:
        errors.append(f"{acc}: no chromosome in the assembly report")
        continue
    if gaps or unplaced:
        tier = "gapped"
    elif both == len(chroms):
        tier = "real_t2t"
    else:
        tier = "gap_free"
    print("\t".join(str(x) for x in (acc, len(chroms), len(gaps), sum(g in chrom_ids for g in gaps), len(unplaced),
                                     both, one, none, tier)))
if errors:
    sys.exit("tier REFUSED:\n  " + "\n  ".join(errors))
