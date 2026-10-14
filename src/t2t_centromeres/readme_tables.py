"""Markdown tables for README.md from annotations/SUMMARY.tsv (so no number in the README is typed by hand).
Usage: python -I readme_tables.py annotations/SUMMARY.tsv annotations/ALL_centromeres.tsv annotations/ALL_official_overlap.tsv > tables.md
"""
import csv
import sys
from collections import OrderedDict

summary, allcen = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(summary), delimiter="\t"))
per = OrderedDict()
for r in rows:
    per.setdefault(r["accession"], r)
fam = {}
for r in csv.DictReader(open(allcen), delimiter="\t"):
    fam.setdefault(r["accession"], r["family"])

print("| organism | accession | assembly | chromosomes with a call | centromeric family (monomer bp) |")
print("|---|---|---|---:|---|")
for acc, r in sorted(per.items(), key=lambda x: (x[1]["organism"], x[1]["asm_name"])):
    cof = " + co-family" if "+" in fam.get(acc, "") else ""
    print(f"| *{r['organism']}* | {acc} | {r['asm_name']} | {r['called']} / {r['chromosomes']} | "
          f"{r['family_base_period']}{cof} |")
print()
print("| assembly | published set | chromosomes compared | locus overlaps it | largest array overlaps it | "
      "published inside our locus | partial | ours inside published | no overlap |")
print("|---|---|---:|---:|---:|---:|---:|---:|---:|")
for r in rows:
    if r["official_set"] in ("", "none_published") or r["official_chroms_comparable"] in ("", "0"):
        continue
    print(f"| {r['asm_name']} | {r['official_set'].split(':')[0]} ({r['official_set'].split(':')[1]}) | "
          f"{r['official_chroms_comparable']} | {r['locus_on_official']} | {r['largest_array_on_official']} | "
          f"{r['official_within_ours']} | {r['partial_overlap']} | {r['ours_within_official']} | {r['disjoint']} |")
print()
print("Chromosomes where our locus does not overlap the published interval at all:")
print()
allov = sys.argv[3] if len(sys.argv) > 3 else None
if allov:
    names = {r["accession"]: r["asm_name"] for r in rows}
    dis = OrderedDict()
    for r in csv.DictReader(open(allov), delimiter="\t"):
        if r["relation"] == "disjoint":
            c, m = r["chromosome"], r.get("molecule", "")
            label = c if (not m or c.split("_")[0].lower().removeprefix("chr") == m.lower()) else f"{c} (molecule {m})"
            dis.setdefault((names.get(r["accession"], r["accession"]), r["set"].split(":")[0]), []).append(label)
    for (asm, src), chs in dis.items():
        print(f"- {asm} vs {src}: {', '.join(chs)}")
