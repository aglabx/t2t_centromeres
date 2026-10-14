"""GFF3 of our centromere call for one assembly, with the published (official) centromere
annotations in the same file and, per chromosome and per official set, where and how they intersect.

Features (1-based, closed; seqid = the name the assembly's 2bit uses):
  centromere      satellome_cen   our locus, one per chromosome
  satellite_DNA   satellome_cen   each array of the locus, Parent = locus; largest=true = candidate active array
  <SO type>       <set source>    each official interval; type by class: active arrays -> satellite_DNA,
                                  spans/regions -> centromere, genes -> gene
  region          intersection    each stretch where our locus and an official interval overlap

Per chromosome and official set, the relation of our locus to the set's span on that chromosome
(min start .. max end of its intervals):
  official_within_ours | ours_within_official | partial_overlap | disjoint |
  no_official_on_chromosome | official_not_comparable (track assembly length differs) | no_call
Written to the GFF (attributes of our locus, one block per set) and to --out-tsv (one row per
chromosome x set; a row with set `none` when the assembly has no official annotation).

--official-set "SOURCE|CLASS|BED|NAMES|REGEX|FAI|ROLE"  (repeatable)
   NAMES: ucsc | seqname | accession | molecule ("chr7"/"Chr07"/"7" -> assigned molecule 7)
   REGEX selects rows by the BED name column; FAI ('-' for none) gates on equal chromosome length;
   ROLE primary | extra (primary is summarised first; both are written in full)
Usage: build_gff.py --calls DIR --report REPORT --acc ACC --out-gff F --out-tsv F [--twobit 2BIT]
           [--official-set SPEC ...]
"""
import argparse
import csv
import re
import shutil
import subprocess
from collections import defaultdict
from urllib.parse import quote

ap = argparse.ArgumentParser()
ap.add_argument("--calls", required=True)
ap.add_argument("--report", required=True)
ap.add_argument("--acc", required=True)
ap.add_argument("--out-gff", required=True)
ap.add_argument("--out-tsv", required=True)
ap.add_argument("--twobit", default=None)
ap.add_argument("--official-set", action="append", default=[])
ap.add_argument("--evidence", default="t2t_centromeres", help="value of the evidence attribute on our features")
args = ap.parse_args()

SO_TYPE = {"active_hor": "satellite_DNA", "presumed_active_array": "satellite_DNA",
           "centromere_region": "centromere", "alpha_satellite_region": "centromere",
           "gene_in_centromeric_region": "gene"}


def esc(v):
    return quote(str(v), safe=" ()/:|.-_+*")


def attrs(d):
    return ";".join(f"{k}={esc(v)}" for k, v in d.items() if v is not None and v != "")


def mol_key(name):
    m = re.sub(r"^[Cc]hr", "", name)
    return m.lstrip("0") or m


# ---- chromosomes ----
chroms = []
with open(args.report) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if f[1] == "assembled-molecule" and f[3] == "Chromosome":
            chroms.append({"name": f[0], "molecule": f[2], "genbank": f[4], "refseq": f[6],
                           "length": int(f[8]), "ucsc": f[9]})
by_name = {c["name"]: c for c in chroms}
twobit_names = None
if args.twobit:
    tool = shutil.which("twoBitInfo")
    if tool is None:
        raise SystemExit("twoBitInfo (UCSC kent tools) is not on PATH")
    out = subprocess.run([tool, args.twobit, "stdout"], check=True,
                         capture_output=True, text=True).stdout
    twobit_names = {ln.split("\t")[0] for ln in out.splitlines()}
for c in chroms:
    if twobit_names is None or c["genbank"] in twobit_names:
        c["seqid"] = c["genbank"]
    elif c["refseq"] in twobit_names:
        c["seqid"] = c["refseq"]
    else:
        raise SystemExit(f"{args.acc}: chromosome {c['name']} is in neither naming of the 2bit")

# ---- our calls ----
with open(f"{args.calls}/centromeres.tsv") as fh:
    calls = {r["chromosome"]: r for r in csv.DictReader(fh, delimiter="\t")}
arrays = defaultdict(list)
with open(f"{args.calls}/centromeric_arrays.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        arrays[r["chromosome"]].append(r)

# ---- official sets ----
sets = []
for spec in args.official_set:
    source, cls, bed, names, regex, fai, role = spec.split("|")
    key = {}
    for c in chroms:
        if names == "ucsc":
            key[c["ucsc"]] = c["name"]
        elif names == "seqname":
            key[c["name"]] = c["name"]
        elif names == "accession":
            key[c["genbank"]] = c["name"]
            key[c["refseq"]] = c["name"]
        elif names == "molecule":
            key[c["molecule"]] = c["name"]
        else:
            raise SystemExit(f"unknown names mode {names}")
    keyf = (lambda n: key.get(mol_key(n))) if names == "molecule" else key.get
    not_comp = set()
    if fai != "-":
        with open(fai) as fh:
            tlen = {ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in fh}
        for tname, cname in key.items():
            if tname in tlen and tlen[tname] != by_name[cname]["length"]:
                not_comp.add(cname)
    live = re.compile(regex)
    iv, unmatched = defaultdict(list), 0
    with open(bed) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4 or not live.search(f[3]):
                continue
            cname = keyf(f[0])
            if cname is None:
                unmatched += 1
                continue
            if cname not in not_comp:
                iv[cname].append((int(f[1]), int(f[2]), f[3]))
    sets.append({"source": source, "class": cls, "role": role, "iv": iv, "not_comp": not_comp,
                 "unmatched": unmatched})
sets.sort(key=lambda s: s["role"] != "primary")


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def relate(call, s, nm):
    offs = sorted(s["iv"].get(nm, []))
    r = {"set": f"{s['source']}:{s['class']}", "official_start": "", "official_end": "",
         "official_n": len(offs), "official_bp": sum(b - a for a, b, _ in offs), "overlap_bp": 0,
         "official_frac_in_our_locus": "", "our_locus_frac_official": "",
         "largest_array_overlap_bp": "", "gap_bp": ""}
    if nm in s["not_comp"]:
        r["relation"] = "official_not_comparable"
        return r, offs
    if not offs:
        r["relation"] = "no_official_on_chromosome"
        return r, offs
    os_, oe = offs[0][0], max(b for _, b, _ in offs)
    r.update({"official_start": os_, "official_end": oe})
    if not call:
        r["relation"] = "no_call"
        return r, offs
    ls, le = int(call["locus_start"]), int(call["locus_end"])
    bs, be = int(call["largest_start"]), int(call["largest_end"])
    o_bp = sum(ov(ls, le, a, b) for a, b, _ in offs)
    if s["class"] == "gene_in_centromeric_region":
        # genes are points inside the authors' region: relate our locus to the region's span
        # (first..last gene), count genes inside our locus; bp overlap with genes is meaningless
        r["n_official_inside_our_locus"] = sum(1 for a, b, _ in offs if a < le and b > ls)
        r["official_bp"] = oe - os_
        o_bp = ov(ls, le, os_, oe)
        offs_for_largest = [(os_, oe, "")]
    else:
        offs_for_largest = offs
    if o_bp == 0:
        rel = "disjoint"
        r["gap_bp"] = min(abs(ls - oe), abs(os_ - le))
    elif ls <= os_ and oe <= le:
        rel = "official_within_ours"
    elif os_ <= ls and le <= oe:
        rel = "ours_within_official"
    else:
        rel = "partial_overlap"
    r.update({"relation": rel, "overlap_bp": o_bp,
              "official_frac_in_our_locus": f"{o_bp / r['official_bp']:.4f}" if r["official_bp"] else "",
              "our_locus_frac_official": f"{o_bp / (le - ls):.4f}",
              "largest_array_overlap_bp": sum(ov(bs, be, a, b) for a, b, _ in offs_for_largest)})
    return r, offs


gff = ["##gff-version 3"] + [f"##sequence-region {c['seqid']} 1 {c['length']}" for c in chroms]
tsv = []
for c in chroms:
    nm, sid = c["name"], c["seqid"]
    call = calls.get(nm)
    rels = [relate(call, s, nm) for s in sets]
    base = {"chromosome": nm, "molecule": c["molecule"], "seqid": sid,
            "our_locus_start": call["locus_start"] if call else "",
            "our_locus_end": call["locus_end"] if call else ""}
    if not sets:
        tsv.append({**base, "set": "none", "relation": "no_official_annotation" if call else "no_call"})
    for r, _ in rels:
        tsv.append({**base, **r})
    if call:
        ls, le = int(call["locus_start"]), int(call["locus_end"])
        bs, be = int(call["largest_start"]), int(call["largest_end"])
        cid = f"cen.{nm}"
        a = {"ID": cid, "Name": f"CEN_{c['molecule']}", "chromosome": nm, "molecule": c["molecule"],
             "family": call["family"], "family_bp": call["family_bp"], "n_arrays": call["n_arrays"],
             "rel_pos_mid": call["rel_pos_mid"], "second_locus_bp": call["second_locus_bp"],
             "colocated_families": call.get("colocated_families", ""),
             "evidence": args.evidence}
        if not sets:
            a["official"] = "none_published"
        for (r, _), s in zip(rels, sets):
            p = s["class"]
            a[f"{p}.source"] = s["source"]
            a[f"{p}.relation"] = r["relation"]
            a[f"{p}.overlap_bp"] = r["overlap_bp"]
            a[f"{p}.official_frac_in_locus"] = r["official_frac_in_our_locus"]
            a[f"{p}.locus_frac_official"] = r["our_locus_frac_official"]
            a[f"{p}.largest_array_overlap_bp"] = r["largest_array_overlap_bp"]
        gff.append("\t".join([sid, "satellome_cen", "centromere", str(ls + 1), str(le), ".", ".", ".", attrs(a)]))
        for i, x in enumerate(sorted(arrays[nm], key=lambda r: int(r["start"]))):
            s0, e0 = int(x["start"]), int(x["end"])
            gff.append("\t".join([sid, "satellome_cen", "satellite_DNA", str(s0 + 1), str(e0), ".", ".", ".", attrs({
                "ID": f"arr.{nm}.{i + 1}", "Parent": cid, "family": x["family"], "membership": x["membership"],
                "fastan_period": x["fastan_period"], "mono_period": x["mono_period"],
                "hor_period": x["hor_period"], "mono_consensus": x["mono_consensus"],
                "largest": "true" if (s0, e0) == (bs, be) else None})]))
    for (r, offs), s in zip(rels, sets):
        so = SO_TYPE.get(s["class"], "sequence_feature")
        for j, (a0, b0, lab) in enumerate(offs):
            oid = f"off.{s['class']}.{nm}.{j + 1}"
            if call:
                ob = ov(ls, le, a0, b0)
                r2 = "inside_our_locus" if ob == b0 - a0 else ("partial" if ob else "outside_our_locus")
            else:
                ob, r2 = 0, "no_call"
            gff.append("\t".join([sid, s["source"], so, str(a0 + 1), str(b0), ".", ".", ".", attrs({
                "ID": oid, "Name": lab, "official_class": s["class"], "relation_to_our_locus": r2,
                "overlap_bp": ob})]))
            if ob:
                gff.append("\t".join([sid, "intersection", "region", str(max(a0, ls) + 1), str(min(b0, le)),
                                      ".", ".", ".", attrs({"ID": f"ix.{s['class']}.{nm}.{j + 1}",
                                                            "ours": f"cen.{nm}", "official": oid,
                                                            "official_class": s["class"], "relation": r2,
                                                            "overlap_bp": ob})]))

with open(args.out_gff, "w") as fh:
    fh.write("\n".join(gff) + "\n")
cols = ["chromosome", "molecule", "seqid", "our_locus_start", "our_locus_end", "set", "relation",
        "official_start", "official_end", "official_n", "n_official_inside_our_locus", "official_bp", "overlap_bp",
        "official_frac_in_our_locus", "our_locus_frac_official", "largest_array_overlap_bp", "gap_bp"]
with open(args.out_tsv, "w") as fh:
    fh.write("\t".join(["accession"] + cols) + "\n")
    for r in tsv:
        fh.write("\t".join([args.acc] + [str(r.get(k, "")) for k in cols]) + "\n")
tally = defaultdict(int)
for r in tsv:
    tally[(r["set"], r["relation"])] += 1
msg = " ".join(f"{s}:{rel}={n}" for (s, rel), n in sorted(tally.items()))
um = " ".join(f"{s['source']}:{s['class']}:unmatched_rows={s['unmatched']}" for s in sets if s["unmatched"])
print(f"[{args.acc}] {msg} {um}".rstrip())
