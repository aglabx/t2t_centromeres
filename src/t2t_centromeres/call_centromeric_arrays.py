"""Centromeric satellite arrays per chromosome from one satellome FasTAN catalogue.

  1. arrays >= --min-array from the FasTAN BED, joined to the arraysplitter summary
     (monomer / HOR period and consensus)
  2. sequence families: monomer consensus vs a doubled + rc-doubled target set (rotation- and
     strand-aware), BLASTN word 7, single linkage at identity >= --ident and query cover >= --cov;
     consensus shorter than --short-len grouped by canonical rotation
  3. superfamilies: sequence families whose base monomer period is equal to or 2x that of the
     superfamily founder, its largest family (+- --period-tol); no chaining member to member
  4. centromeric superfamily: among those with a base period (a consensus >= --short-len), present
     (>= --breadth-min-bp) on the most chromosomes, ties by bp;
     co-families = broad superfamilies co-located with it (--cofamily-breadth, --cofamily-coloc)
  5. per chromosome: its arrays (plus monomer-period-class arrays, k = 1..4) merged across gaps
     < --merge-gap; the cluster with the most bp is the locus. Cap mode (the superfamily also sits in
     both terminal 10 % of most chromosomes): an interior cluster is preferred.

Writes centromeres.tsv, centromeric_arrays.tsv, centromeric_arrays.bed, family_breadth.tsv,
satellite_arrays.tsv (every array used, with its superfamily), outliers.tsv, exceptions.tsv into --out. Refuses (non-zero exit) when the catalogue's sequence
names cannot be tied to the assembly's chromosomes (NCBI assembly report).

Usage: call_centromeric_arrays.py --acc ACC --catalogue RUN_DIR --report ASSEMBLY_REPORT \
           --asm-name NAME --out DIR [--threads N]
"""
import argparse
import os
import shutil
import sys
from collections import defaultdict

from catalogue import array_id, read_bed, read_summary
from families import UF, canonical, rc, run

ap = argparse.ArgumentParser()
ap.add_argument("--acc", required=True)
ap.add_argument("--asm-name", required=True)
ap.add_argument("--catalogue", required=True)
ap.add_argument("--report", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--threads", type=int, default=16)
ap.add_argument("--min-array", type=int, default=10_000)
ap.add_argument("--ident", type=float, default=0.70)
ap.add_argument("--cov", type=float, default=0.80)
ap.add_argument("--short-len", type=int, default=30)
ap.add_argument("--chrom-min-bp", type=int, default=100_000,
                help="bp of a superfamily that counts for co-location, both-ends (cap mode) and co-located families")
ap.add_argument("--breadth-min-bp", type=int, default=10_000,
                help="bp of a superfamily on a chromosome for it to count toward breadth (default: one array of --min-array)")
ap.add_argument("--merge-gap", type=int, default=1_000_000)
ap.add_argument("--ambiguous-frac", type=float, default=0.50)
ap.add_argument("--end-frac", type=float, default=0.01)
ap.add_argument("--cofamily-breadth", type=float, default=0.75,
                help="a runner-up superfamily on >= this fraction of the chosen one's chromosomes ...")
ap.add_argument("--cofamily-coloc", type=float, default=0.80,
                help="... joins the centromeric class when its largest cluster lies within --merge-gap of the chosen one's on >= this fraction of shared chromosomes")
ap.add_argument("--period-tol", type=float, default=0.03,
                help="an array joins the centromeric family's locus when its monomer period is k x the family base period (k=1..4) within this tolerance")
ap.add_argument("--kmer", type=int, default=0, help="mmseqs -k (only with --aligner mmseqs); 0 = mmseqs default")
ap.add_argument("--aligner", choices=["mmseqs", "blastn"], default="blastn",
                help="blastn (default) or mmseqs on the same doubled-target design; the mmseqs prefilter misses diverged alpha (CHM13 chr18/11/6/X active arrays come out as one-array families)")
ap.add_argument("--max-targets", type=int, default=200,
                help="blastn -max_target_seqs; enough for union-find connectivity, keeps dense satellite families tractable")
ap.add_argument("--keep-hits", action="store_true", help="keep the MMseqs hits table as family_hits.tsv")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)


def die(msg):
    sys.stderr.write(f"[{args.acc}] REFUSED: {msg}\n")
    raise SystemExit(2)


# ---- chromosomes from the NCBI assembly report ----
chroms = {}      # catalogue-side name -> chromosome record (both GenBank and RefSeq keys)
chrom_list = []
with open(args.report) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if f[1] != "assembled-molecule" or f[3] != "Chromosome":
            continue
        rec = {"name": f[0], "genbank": f[4], "refseq": f[6], "length": int(f[8])}
        chrom_list.append(rec)
        for key in (f[4], f[6]):
            if key and key != "na":
                chroms[key] = rec
if not chrom_list:
    die("assembly report lists no assembled-molecule chromosomes")

# ---- arrays ----
bed = read_bed(args.catalogue)
bed_seqs = {r[0] for r in bed}
mapped = bed_seqs & set(chroms)
if not mapped:
    die(f"no catalogue sequence is a chromosome of this assembly "
        f"(catalogue e.g. {sorted(bed_seqs)[:3]}, report e.g. {sorted(chroms)[:3]})")
for s in mapped:
    over = [r for r in bed if r[0] == s and r[2] > chroms[s]["length"]]
    if over:
        die(f"{s}: array ends at {over[0][2]} beyond chromosome length {chroms[s]['length']}")

summ = read_summary(args.catalogue, min_len=args.min_array)
arrays, outliers = [], []
for seq, s, e, p, _ in bed:
    if e - s < args.min_array or seq not in chroms:
        continue
    m = summ.get(array_id(seq, s, e, p))
    cons = (m["mono_consensus"].upper() if m else "")
    if not m:
        outliers.append(("array_without_summary_row", chroms[seq]["name"], seq, s, e, p, ""))
        continue
    if not cons or set(cons) - set("ACGTN"):
        outliers.append(("array_without_usable_consensus", chroms[seq]["name"], seq, s, e, p, cons[:40]))
        continue
    arrays.append({"id": f"a{len(arrays)}", "seq": seq, "chrom": chroms[seq]["name"], "start": s,
                   "end": e, "len": e - s, "period": p, "mono_period": m["mono_period"],
                   "hor_period": m["hor_period"], "cons": cons, "hor_cons": m["hor_consensus"]})
if not arrays:
    die("no array >= min_array with a usable consensus on any chromosome")
byid = {a["id"]: a for a in arrays}

# ---- sequence families (rotation/strand-aware doubled-target search, single linkage) ----
uf = UF()
for a in arrays:
    uf.find(a["id"])
long_ids = [a["id"] for a in arrays if len(a["cons"]) >= args.short_len]
short_ids = [a["id"] for a in arrays if len(a["cons"]) < args.short_len]
tmp = os.path.join(args.out, "_mmseqs")
if long_ids:
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    qfa, tfa, hits = (os.path.join(tmp, x) for x in ("query.fasta", "target.fasta", "hits.tsv"))
    with open(qfa, "w") as q, open(tfa, "w") as t:
        for mid in long_ids:
            s = byid[mid]["cons"]
            q.write(f">{mid}\n{s}\n")
            t.write(f">{mid}__f\n{s + s}\n")
            t.write(f">{mid}__r\n{rc(s) + rc(s)}\n")
    if args.aligner == "mmseqs":
        run(["mmseqs", "createdb", qfa, os.path.join(tmp, "qdb")])
        run(["mmseqs", "createdb", tfa, os.path.join(tmp, "tdb")])
        run(["mmseqs", "search", os.path.join(tmp, "qdb"), os.path.join(tmp, "tdb"),
             os.path.join(tmp, "res"), os.path.join(tmp, "mmtmp"),
             "--min-seq-id", "0.5", "-c", "0.5", "--cov-mode", "2", "-s", "7.5",
             "--search-type", "3", "-a", "--threads", str(args.threads), "--max-seqs", "2000",
             "-k", str(args.kmer)])
        run(["mmseqs", "convertalis", os.path.join(tmp, "qdb"), os.path.join(tmp, "tdb"),
             os.path.join(tmp, "res"), hits, "--format-output", "query,target,fident,qcov"])
    else:
        # Same doubled + rc-doubled target and the same fident/qcov columns as the mmseqs path.
        # Only the plus strand of the target is searched: the rc copy is already in the DB.
        run(["makeblastdb", "-in", tfa, "-dbtype", "nucl", "-out", os.path.join(tmp, "tdb")])
        raw = os.path.join(tmp, "blast.tsv")
        run(["blastn", "-task", "blastn", "-word_size", "7", "-strand", "plus", "-evalue", "1e-5",
             "-query", qfa, "-db", os.path.join(tmp, "tdb"), "-num_threads", str(args.threads),
             "-max_target_seqs", str(args.max_targets), "-max_hsps", "1", "-outfmt", "6 qseqid sseqid pident qstart qend qlen",
             "-out", raw])
        with open(raw) as fh, open(hits, "w") as out:
            for ln in fh:
                q, t, pid, qs, qe, ql = ln.rstrip("\n").split("\t")
                out.write(f"{q}\t{t}\t{float(pid) / 100:.3f}\t{(int(qe) - int(qs) + 1) / int(ql):.3f}\n")
    if args.keep_hits:
        shutil.copy(hits, os.path.join(args.out, "family_hits.tsv"))
        with open(os.path.join(args.out, "array_ids.tsv"), "w") as fh:
            for a in arrays:
                fh.write(f"{a['id']}\t{a['chrom']}\t{a['start']}\t{a['len']}\t{a['mono_period']}\n")
    with open(hits) as fh:
        for ln in fh:
            q, tgt, fid, qc = ln.rstrip("\n").split("\t")
            tb = tgt.rsplit("__", 1)[0]
            if q != tb and float(fid) >= args.ident and float(qc) >= args.cov:
                uf.union(q, tb)
    shutil.rmtree(tmp)
by_canon = defaultdict(list)
for mid in short_ids:
    by_canon[canonical(byid[mid]["cons"])].append(mid)
for grp in by_canon.values():
    for x in grp[1:]:
        uf.union(grp[0], x)

fam_members = defaultdict(list)
for a in arrays:
    fam_members[uf.find(a["id"])].append(a)
# stable family names: F001.. by total bp
fam_order = sorted(fam_members, key=lambda r: -sum(a["len"] for a in fam_members[r]))
fam_name = {r: f"F{i + 1:03d}" for i, r in enumerate(fam_order)}
for r, mem in fam_members.items():
    for a in mem:
        a["family"] = fam_name[r]

def _mono(a):
    try:
        return int(a["mono_period"])
    except ValueError:
        return 0


def _base(mem):
    """bp-weighted median monomer period of arrays with a consensus >= short_len; 0 if none."""
    long_ = sorted((a for a in mem if _mono(a) >= args.short_len), key=_mono)
    half, acc = sum(a["len"] for a in long_) / 2, 0
    for a in long_:
        acc += a["len"]
        if acc >= half:
            return _mono(a)
    return 0


def _k_of(p, q, kmax):
    """k if p = k x q within period_tol x k (p >= q), else 0."""
    ratio = p / q
    k = round(ratio)
    return k if 1 <= k <= kmax and abs(ratio - k) <= args.period_tol * k else 0


# ---- superfamilies: sequence families whose base monomer periods are equal or 2x (+-period_tol) ----
# Breadth on raw sequence families is unstable: in bonobo mat alpha split into 18 families, the largest
# on 18 chromosomes - the same breadth as the 32-bp subterminal satellite, which then won on bp.
# Each family is compared with a superfamily's founder (its largest family by bp), never with another
# member: single linkage on period drifts. bonobo pri: subterminal families of 62..141 bp chained in
# 3 % steps up to ~86 bp, 2 x 86 ~ 171, and a 70-bp satellite joined alpha (loci moved to 0.4 Mb).
fam_base = {fam_name[r]: _base(fam_members[r]) for r in fam_order}
suf = UF()
for fn in fam_base:
    suf.find(fn)
founders = []
for fn in (fam_name[r] for r in fam_order):          # largest family first
    b = fam_base[fn]
    if b < args.short_len:
        continue
    for f0 in founders:
        b0 = fam_base[f0]
        if _k_of(max(b, b0), min(b, b0), 2):
            suf.union(f0, fn)
            break
    else:
        founders.append(fn)
sup_members = defaultdict(list)
for a in arrays:
    sup_members[suf.find(a["family"])].append(a)
sup_order = sorted(sup_members, key=lambda r: -sum(a["len"] for a in sup_members[r]))
sup_name = {r: f"S{i + 1:03d}" for i, r in enumerate(sup_order)}
for r, mem in sup_members.items():
    for a in mem:
        a["sfam"] = sup_name[r]

# ---- breadth ----
sup_chrom_bp = defaultdict(lambda: defaultdict(int))
for a in arrays:
    sup_chrom_bp[a["sfam"]][a["chrom"]] += a["len"]
chrom_len = {c["name"]: c["length"] for c in chrom_list}


def _both_ends(sf):
    """chromosomes carrying >= chrom_min_bp of the superfamily in BOTH terminal 10 %."""
    pe, qe = defaultdict(int), defaultdict(int)
    for a in sup_members_by_name[sf]:
        L = chrom_len[a["chrom"]]
        if a["end"] <= 0.10 * L:
            pe[a["chrom"]] += a["len"]
        if a["start"] >= 0.90 * L:
            qe[a["chrom"]] += a["len"]
    return sum(1 for c in pe if pe[c] >= args.chrom_min_bp and qe.get(c, 0) >= args.chrom_min_bp)


sup_members_by_name = {sup_name[r]: m for r, m in sup_members.items()}
breadth = []
for r in sup_order:
    sn = sup_name[r]
    mem = sup_members[r]
    rep = max(mem, key=lambda a: a["len"])
    # Breadth = chromosomes carrying the superfamily at all. At >= 100 kb bonobo mat alpha (present on
    # 24 chromosomes, centromeres as small as 93 kb) tied the 32-bp subterminal satellite 18 vs 18.
    n_ch = sum(1 for bp in sup_chrom_bp[sn].values() if bp >= args.breadth_min_bp)
    n_big = sum(1 for bp in sup_chrom_bp[sn].values() if bp >= args.chrom_min_bp)
    fams = sorted({a["family"] for a in mem})
    breadth.append({"family": sn, "n_arrays": len(mem), "total_bp": sum(a["len"] for a in mem),
                    "n_chrom": n_ch, "n_chrom_big": n_big, "rep": rep, "n_families": len(fams),
                    "families": ",".join(fams[:20]) + ("..." if len(fams) > 20 else ""),
                    "base_period": _base(mem), "both_ends": _both_ends(sn)})
breadth.sort(key=lambda b: (-b["n_chrom"], -b["total_bp"]))
# A superfamily without a base period (every consensus < --short-len: microsatellites, the telomeric
# repeat) is never the centromeric one. Under breadth = presence the telomere, >= 10 kb on all 28 sheep
# chromosomes, outranked SatI (26) and moved every sheep locus to a chromosome end.
eligible = [b for b in breadth if b["base_period"] > 0]
if not eligible:
    die(f"no superfamily has a monomer consensus >= --short-len {args.short_len}; none can be centromeric")
cen = eligible[0]
cen_fam = cen["family"]


def _largest_cluster(sf, chrom):
    arr = sorted((a for a in sup_members_by_name[sf] if a["chrom"] == chrom), key=lambda a: a["start"])
    cl = []
    for a in arr:
        if cl and a["start"] - cl[-1][1] < args.merge_gap:
            cl[-1] = [cl[-1][0], max(cl[-1][1], a["end"]), cl[-1][2] + a["len"]]
        else:
            cl.append([a["start"], a["end"], a["len"]])
    return max(cl, key=lambda c: c[2]) if cl else None


# ---- co-families: broad superfamilies sitting in the same locus as the chosen one ----
# goat: SatI (816) and SatII (703) are both on ~every chromosome at the same end; breadth alone flipped
# between assemblies of one individual (goat1.0 816, goat2.0 703 by 30 vs 29 chromosomes).
cofams, coloc_note = [], {}
for b in breadth:
    if b is cen or b["n_chrom"] < args.cofamily_breadth * cen["n_chrom"] or b["base_period"] == 0:
        continue
    shared = [c for c in sup_chrom_bp[cen_fam] if sup_chrom_bp[cen_fam][c] >= args.chrom_min_bp
              and sup_chrom_bp[b["family"]].get(c, 0) >= args.chrom_min_bp]
    near = 0
    for c in shared:
        x, y = _largest_cluster(cen_fam, c), _largest_cluster(b["family"], c)
        if x and y and max(x[0], y[0]) - min(x[1], y[1]) < args.merge_gap:
            near += 1
    coloc_note[b["family"]] = f"{near}/{len(shared)}"
    if shared and near / len(shared) >= args.cofamily_coloc:
        cofams.append(b["family"])
cen_label = "+".join([cen_fam] + cofams)
# Cap mode: the chosen superfamily also forms subterminal caps (both terminal 10 % on > half its
# chromosomes; siamang alpha: 23 of 26). Then a terminal cluster is never the centromere if an
# interior one (midpoint within 10-90 %) exists. Off for acrocentric genomes (caps absent).
cap_mode = cen["both_ends"] > 0.5 * cen["n_chrom"]

# ---- period class: arrays whose monomer period is the family's, k = 1..4 ----
# One consensus per array is a fragile family key for alpha: on CHM13 the largest active arrays of
# chr2/6/7/X miss the 0.70/0.80 cut against the rest of alpha (0.69 at 0.94 cov; 0.74 at 0.57 cov).
for a in arrays:
    a["membership"] = "family" if a["sfam"] == cen_fam else ("cofamily" if a["sfam"] in cofams else "")
fam_long = [a for a in arrays if a["sfam"] == cen_fam and _mono(a) >= args.short_len]
base_period = 0
if fam_long:
    acc_bp, half = 0, sum(a["len"] for a in fam_long) / 2
    for a in sorted(fam_long, key=_mono):
        acc_bp += a["len"]
        if acc_bp >= half:
            base_period = _mono(a)
            break
    for a in arrays:
        m = _mono(a)
        if a["membership"] or m < args.short_len:
            continue
        ratio = m / base_period
        k = round(ratio)
        if 1 <= k <= 4 and abs(ratio - k) <= args.period_tol * k:
            a["membership"] = f"period_class_k{k}"

# ---- loci ----
exceptions = []
if cap_mode:
    exceptions.append(("cap_mode", "*", cen_fam,
                       f"{cen['both_ends']} of {cen['n_chrom']} chromosomes carry it in both terminal 10%",
                       "interior clusters preferred over terminal ones on every chromosome"))
for f in cofams:
    exceptions.append(("cofamily_joined_locus", "*", f, f"co-located on {coloc_note[f]} shared chromosomes",
                       "arrays of this superfamily are part of the centromeric locus"))
for b in breadth:
    if b is cen or b["family"] in cofams:
        continue
    if b["n_chrom"] >= cen["n_chrom"] - 1:
        if b["base_period"] == 0:
            kind, applied = "short_consensus_family_as_broad", "not eligible: no consensus >= --short-len"
        else:
            kind = "breadth_tie_broken_by_bp" if b["n_chrom"] == cen["n_chrom"] else "breadth_near_tie"
            applied = "kept the broader family; see colocated_families"
        exceptions.append((kind, "*", b["family"],
                           f"{b['n_chrom']} chromosomes, {b['total_bp']} bp, mono {b['rep']['mono_period']} vs "
                           f"{cen['n_chrom']} chromosomes, {cen['total_bp']} bp, mono {cen['rep']['mono_period']}",
                           applied))
rows, loci_arrays, chrom_rows_missing = [], [], []
for ch in chrom_list:
    mine = sorted((a for a in arrays if a["chrom"] == ch["name"] and a["membership"]),
                  key=lambda a: a["start"])
    if not mine:
        outliers.append(("chromosome_without_centromeric_family_array", ch["name"],
                         ch["genbank"], "", "", "", f"length {ch['length']}"))
        continue
    clusters = []
    for a in mine:
        if clusters and a["start"] - clusters[-1]["end"] < args.merge_gap:
            c = clusters[-1]
            c["end"] = max(c["end"], a["end"])
            c["bp"] += a["len"]
            c["arrays"].append(a)
        else:
            clusters.append({"start": a["start"], "end": a["end"], "bp": a["len"], "arrays": [a]})
    clusters.sort(key=lambda c: -c["bp"])
    if cap_mode:
        Lc = ch["length"]
        inner = [c for c in clusters if 0.10 * Lc <= (c["start"] + c["end"]) / 2 <= 0.90 * Lc]
        if inner and inner[0] is not clusters[0]:
            exceptions.append(("subterminal_cap_skipped", ch["name"], cen_fam,
                               f"terminal {clusters[0]['bp']} bp at {clusters[0]['start']}-{clusters[0]['end']}; "
                               f"interior {inner[0]['bp']} bp chosen", "cap mode"))
            clusters = inner + [c for c in clusters if c not in inner]
        elif not inner:
            exceptions.append(("cap_mode_no_interior_cluster", ch["name"], cen_fam,
                               f"all {len(clusters)} clusters terminal", "kept the largest terminal cluster"))
    loc = clusters[0]
    second = clusters[1]["bp"] if len(clusters) > 1 else 0
    big = max(loc["arrays"], key=lambda a: a["len"])
    L = ch["length"]
    mid = (loc["start"] + loc["end"]) / 2 / L
    near = defaultdict(int)
    for a in arrays:
        if (a["chrom"] == ch["name"] and not a["membership"]
                and a["start"] < loc["end"] + args.merge_gap and a["end"] > loc["start"] - args.merge_gap):
            near[a["sfam"]] += a["len"]
    mono_of = {b["family"]: b["rep"]["mono_period"] for b in breadth}
    coloc = ";".join(f"{f}:{bp}:m{mono_of[f]}" for f, bp in sorted(near.items(), key=lambda x: -x[1])
                     if bp >= args.chrom_min_bp)
    rows.append([args.acc, args.asm_name, ch["name"], big["seq"], ch["genbank"], ch["refseq"], L,
                 cen_label, loc["start"], loc["end"], loc["end"] - loc["start"], loc["bp"],
                 len(loc["arrays"]), f"{mid:.4f}", second, big["start"], big["end"], big["len"],
                 big["period"], big["mono_period"], big["hor_period"], big["cons"], coloc or "-"])
    for a in loc["arrays"]:
        loci_arrays.append(a)
    by_period = sum(a["len"] for a in loc["arrays"] if a["membership"].startswith("period"))
    if by_period:
        exceptions.append(("locus_includes_period_class_arrays", ch["name"], cen_fam,
                           f"{by_period} of {loc['bp']} bp", "joined by monomer period, not by sequence"))
    if second >= args.ambiguous_frac * loc["bp"]:
        exceptions.append(("second_locus_comparable", ch["name"], cen_fam,
                           f"second {second} bp vs first {loc['bp']} bp", "kept the larger"))
    if loc["start"] < args.end_frac * L or loc["end"] > (1 - args.end_frac) * L:
        exceptions.append(("locus_at_chromosome_end", ch["name"], cen_fam,
                           f"{loc['start']}-{loc['end']} of {L}", "informative: telocentric"))
    cls_bp = sum(a["len"] for a in arrays if a["chrom"] == ch["name"] and a["membership"])
    other = defaultdict(int)
    for a in arrays:
        if a["chrom"] == ch["name"] and not a["membership"]:
            other[a["sfam"]] += a["len"]
    if other:
        top = max(other, key=other.get)
        if other[top] > cls_bp:
            exceptions.append(("other_family_larger_on_chromosome", ch["name"], top,
                               f"{other[top]} bp vs {cls_bp} bp", "centromeric family kept"))

# ---- write ----
def write(name, header, data):
    with open(os.path.join(args.out, name), "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in data:
            fh.write("\t".join(str(x) for x in r) + "\n")


write("centromeres.tsv",
      ["accession", "asm_name", "chromosome", "seq_name", "genbank_acc", "refseq_acc", "chrom_len",
       "family", "locus_start", "locus_end", "locus_span", "family_bp", "n_arrays", "rel_pos_mid",
       "second_locus_bp", "largest_start", "largest_end", "largest_len", "fastan_period",
       "mono_period", "hor_period", "mono_consensus", "colocated_families"], rows)
write("centromeric_arrays.tsv",
      ["chromosome", "seq_name", "start", "end", "length", "family", "seq_family", "membership",
       "fastan_period", "mono_period", "hor_period", "mono_consensus", "hor_consensus"],
      [[a["chrom"], a["seq"], a["start"], a["end"], a["len"], a["sfam"], a["family"], a["membership"], a["period"],
        a["mono_period"], a["hor_period"], a["cons"], a["hor_cons"]] for a in loci_arrays])
with open(os.path.join(args.out, "centromeric_arrays.bed"), "w") as fh:
    for a in loci_arrays:
        fh.write(f"{a['seq']}\t{a['start']}\t{a['end']}\t{a['chrom']}|{a['sfam']}|p{a['period']}"
                 f"|m{a['mono_period']}\t0\t.\n")
write("family_breadth.tsv",
      ["family", "n_arrays", "total_bp", "n_chrom_breadth", "n_chrom_ge_chrom_min_bp", "chosen", "base_period", "n_seq_families",
       "chroms_both_ends", "rep_seq", "rep_start", "rep_len", "rep_fastan_period", "rep_mono_period",
       "rep_mono_consensus", "seq_families"],
      [[b["family"], b["n_arrays"], b["total_bp"], b["n_chrom"], b["n_chrom_big"], (1 if b["family"] == cen_fam else (2 if b["family"] in cofams else 0)),
        b["base_period"], b["n_families"], b["both_ends"], b["rep"]["seq"], b["rep"]["start"],
        b["rep"]["len"], b["rep"]["period"], b["rep"]["mono_period"], b["rep"]["cons"], b["families"]]
       for b in breadth])
chrom_order = {c["name"]: i for i, c in enumerate(chrom_list)}
write("satellite_arrays.tsv",
      ["chromosome", "seq_name", "start", "end", "length", "family", "seq_family", "membership",
       "fastan_period", "mono_period", "hor_period"],
      [[a["chrom"], a["seq"], a["start"], a["end"], a["len"], a["sfam"], a["family"], a["membership"] or "-",
        a["period"], a["mono_period"], a["hor_period"]]
       for a in sorted(arrays, key=lambda a: (chrom_order[a["chrom"]], a["start"], a["end"]))])
write("outliers.tsv", ["kind", "chromosome", "seq_or_acc", "start", "end", "period", "note"], outliers)
write("exceptions.tsv", ["kind", "chromosome", "family", "magnitude", "default_applied"], exceptions)
print(f"[{args.acc}] chromosomes={len(chrom_list)} called={len(rows)} "
      f"no_family_array={sum(1 for o in outliers if o[0].startswith('chromosome'))} "
      f"family={cen_label} n_chrom={cen['n_chrom']} base_period={base_period} "
      f"exceptions={len(exceptions)} arrays_unusable={sum(1 for o in outliers if o[0].startswith('array'))}")
