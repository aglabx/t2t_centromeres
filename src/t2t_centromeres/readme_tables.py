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

print("| organism | assembly | accession | chromosomes with a call | centromeric family (monomer bp) | karyotype |")
print("|---|---|---|---:|---|---|")
for acc, r in sorted(per.items(), key=lambda x: (x[1]["organism"], x[1]["asm_name"])):
    cof = " + co-family" if "+" in fam.get(acc, "") else ""
    d = f"annotations/{r['directory']}"
    print(f"| *{r['organism']}* | [{r['asm_name']}]({d}) | {acc} | {r['called']} / {r['chromosomes']} | "
          f"{r['family_base_period']}{cof} | [svg]({d}/karyotype.svg) |")
print()
print("| assembly | published set | chromosomes compared | agrees | contested: no overlap | contested: no call |")
print("|---|---|---:|---:|---:|---:|")
for r in rows:
    if r["official_set"] in ("", "none_published") or r["official_chroms_comparable"] in ("", "0"):
        continue
    print(f"| {r['asm_name']} | {r['official_set'].split(':')[0]} ({r['official_set'].split(':')[1]}) | "
          f"{r['official_chroms_comparable']} | {r['agrees']} | {r['contested_locus']} | {r['contested_no_call']} |")
print()
print("**Contested** — every chromosome where our call contradicts a published annotation "
      "(`agreement` in `official_overlap.tsv`, `published_agreement=contested` in the GFF):")
print()
allov = sys.argv[3] if len(sys.argv) > 3 else None
if allov:
    names = {r["accession"]: r["asm_name"] for r in rows}
    kinds = {"contested_locus": "no overlap", "contested_no_call": "published centromere, no call"}
    dis = OrderedDict()
    for r in csv.DictReader(open(allov), delimiter="\t"):
        if r["agreement"].startswith("contested"):
            c, m = r["chromosome"], r.get("molecule", "")
            num = c.split("_")[0].lower().removeprefix("chr")
            label = c if (not m or (num.lstrip("0") or num) == (m.lower().lstrip("0") or m.lower())) else f"{c} (molecule {m})"
            dis.setdefault((names.get(r["accession"], r["accession"]), r["set"].split(":")[0], kinds[r["agreement"]]), []).append(label)
    for (asm, src, kind), chs in dis.items():
        print(f"- {asm} vs {src}, {kind}: {', '.join(chs)}")
