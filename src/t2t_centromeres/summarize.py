"""One row per assembly: what was called, which satellite superfamily, and how the call meets every
official set (relation counts; locus and largest array on the official interval, counted over the
chromosomes where the set has an interval and is comparable).
Usage: python -I summarize.py CALLS_DIR ASSEMBLIES.tsv > SUMMARY.tsv   (config/assemblies.tsv)
"""
import csv
import os
import sys
from collections import defaultdict

calls, manifest = sys.argv[1], sys.argv[2]
org = {}
with open(manifest) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("tier", "t2t") == "t2t":
            org[r["accession"]] = r["organism"]
with open(os.path.join(calls, "STATUS.tsv")) as fh:
    status = list(csv.DictReader(fh, delimiter="\t"))

cols = ["accession", "organism", "asm_name", "directory", "status", "chromosomes", "called", "no_family_array",
        "family_base_period", "family_n_chrom", "family_n_seq_families", "near_tie", "official_set",
        "official_chroms_comparable", "official_within_ours", "partial_overlap", "ours_within_official",
        "disjoint", "locus_on_official", "largest_array_on_official", "not_comparable",
        "agrees", "contested_locus", "contested_active_array", "contested_no_call", "note"]
print("\t".join(cols))
for st in status:
    acc = st["accession"]
    row = {k: "" for k in cols}
    row.update({"accession": acc, "organism": org.get(acc, ""), "asm_name": st["asm_name"], "directory": st["directory"],
                "status": st["status"]})
    d = os.path.join(calls, st["directory"])
    if st["status"] != "called":
        row["note"] = st["summary"][:160]
        print("\t".join(str(row[k]) for k in cols))
        continue
    kv = dict(x.split("=", 1) for x in st["summary"].split()[1:] if "=" in x)
    row.update({"chromosomes": kv.get("chromosomes"), "called": kv.get("called"),
                "no_family_array": kv.get("no_family_array"), "family_base_period": kv.get("base_period"),
                "family_n_chrom": kv.get("n_chrom")})
    with open(os.path.join(d, "family_breadth.tsv")) as fh:
        for b in csv.DictReader(fh, delimiter="\t"):
            if b["chosen"] == "1":
                row["family_n_seq_families"] = b["n_seq_families"]
    with open(os.path.join(d, "exceptions.tsv")) as fh:
        ties = [e for e in csv.DictReader(fh, delimiter="\t") if e["kind"].startswith("breadth_")]
    row["near_tie"] = ";".join(f"{e['kind']}:{e['magnitude']}" for e in ties)
    ov = os.path.join(d, "official_overlap.tsv")
    if not os.path.exists(ov):
        sys.exit(f"{acc}: status 'called' but {ov} is missing")
    sets = defaultdict(list)
    with open(ov) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sets[r["set"]].append(r)
    if not sets or list(sets) == ["none"]:
        row["official_set"] = "none_published"
        print("\t".join(str(row[k]) for k in cols))
        continue
    for name, rs in sets.items():
        out = dict(row)
        rel, agr = defaultdict(int), defaultdict(int)
        for r in rs:
            rel[r["relation"]] += 1
            if not r.get("agreement"):
                sys.exit(f"{acc}: official_overlap.tsv row without agreement ({r['chromosome']}, {name})")
            agr[r["agreement"]] += 1
        comp = [r for r in rs if r["relation"] in ("official_within_ours", "partial_overlap",
                                                   "ours_within_official", "disjoint")]
        out.update({"official_set": name, "official_chroms_comparable": len(comp),
                    "official_within_ours": rel["official_within_ours"], "partial_overlap": rel["partial_overlap"],
                    "ours_within_official": rel["ours_within_official"], "disjoint": rel["disjoint"],
                    "locus_on_official": sum(1 for r in comp if int(r["overlap_bp"] or 0) > 0),
                    "largest_array_on_official": sum(1 for r in comp if int(r["largest_array_overlap_bp"] or 0) > 0),
                    "not_comparable": rel["official_not_comparable"], "agrees": agr["agrees"],
                    "contested_locus": agr["contested_locus"], "contested_active_array": agr["contested_active_array"],
                    "contested_no_call": agr["contested_no_call"]})
        print("\t".join(str(out[k]) for k in cols))
