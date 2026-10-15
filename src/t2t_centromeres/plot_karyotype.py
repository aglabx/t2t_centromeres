"""Karyotype figure of one assembly, two panels.

  a  whole chromosomes to one scale: every satellite array >= 10 kb coloured by superfamily, our
     centromeric locus boxed, its largest array marked, each published centromere interval bracketed
     beside it, one dot per published set: agrees / contested (no overlap, other active array, no call)
  b  the centromeric region of every chromosome, each window to its own scale (width printed above it):
     the arrays inside the locus, the locus, its largest array and the published interval side by side;
     a published interval more than 3 Mb away is drawn as an arrow with its distance

Hover over any element for coordinates. Plain SVG written directly (no plotting library), so the file is
byte-identical between runs and sits in MANIFEST.md5 with the tables.

Usage: plot_karyotype.py --calls DIR --report ASSEMBLY_REPORT --acc ACC --asm-name NAME --organism ORG --out SVG
  DIR = call_centromeric_arrays.py + build_gff.py output (centromeres.tsv, satellite_arrays.tsv,
        family_breadth.tsv, official_overlap.tsv, exceptions.tsv)
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

ap = argparse.ArgumentParser()
for a in ("--calls", "--report", "--acc", "--asm-name", "--organism", "--out"):
    ap.add_argument(a, required=True)
args = ap.parse_args()


def rows(name):
    path = os.path.join(args.calls, name)
    if not os.path.isfile(path):
        sys.exit(f"[{args.acc}] plot REFUSED: {path} missing")
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def mb(x):
    return f"{int(x) / 1e6:.2f}"


def size(x):
    x = int(x)
    return f"{x / 1e6:.2f} Mb" if x >= 1_000_000 else f"{x / 1e3:.0f} kb"


def short_size(x):
    return f"{x / 1e6:.1f} Mb" if x >= 1_000_000 else f"{x / 1e3:.0f} kb"


def f1(x):
    return f"{x:.1f}"


# ---- inputs ----
chroms, molecule = [], {}
with open(args.report, encoding="utf-8") as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if f[1] == "assembled-molecule" and f[3] == "Chromosome":
            chroms.append((f[0], int(f[8])))
            molecule[f[0]] = f[2]
if not chroms:
    sys.exit(f"[{args.acc}] plot REFUSED: assembly report lists no assembled-molecule chromosomes")
known = {c for c, _ in chroms}
chrom_len = dict(chroms)

breadth = rows("family_breadth.tsv")
arrays = rows("satellite_arrays.tsv")
loci = {r["chromosome"]: r for r in rows("centromeres.tsv")}
overlap = rows("official_overlap.tsv")
exc_by_chrom = defaultdict(list)
for r in rows("exceptions.tsv"):
    exc_by_chrom[r["chromosome"]].append(r["kind"])
for name, rs in (("satellite_arrays.tsv", arrays), ("centromeres.tsv", loci.values()),
                 ("official_overlap.tsv", overlap)):
    bad = sorted({r["chromosome"] for r in rs} - known)
    if bad:
        sys.exit(f"[{args.acc}] plot REFUSED: {name} names chromosomes not in the assembly report: {bad[:5]}")
# build_gff writes set "none" (relation no_official_annotation / no_call) when the assembly has no published set
KNOWN = {"official_within_ours", "ours_within_official", "partial_overlap", "disjoint", "no_call",
         "no_official_on_chromosome", "official_not_comparable", "no_official_annotation"}
odd = sorted({r["relation"] for r in overlap} - KNOWN)
if odd:
    sys.exit(f"[{args.acc}] plot REFUSED: official_overlap.tsv has unknown relations {odd}")
AGREE = {"agrees", "contested_locus", "contested_active_array", "contested_no_call", "not_compared"}
odd = sorted({r.get("agreement") or "<empty>" for r in overlap} - AGREE)
if odd:
    sys.exit(f"[{args.acc}] plot REFUSED: official_overlap.tsv has unknown or missing agreement {odd}")
overlap = [r for r in overlap if r["set"] != "none"]
sets = []
for r in overlap:
    if r["set"] not in sets:
        sets.append(r["set"])
rel = {(r["chromosome"], r["set"]): r for r in overlap}


def published(c, s):
    """(start, end, relation) of set s on chromosome c, or None when there is nothing to draw"""
    o = rel.get((c, s))
    if not o or not o["official_start"] or o["relation"] in ("no_official_on_chromosome", "official_not_comparable"):
        return None
    return int(o["official_start"]), int(o["official_end"]), o["relation"]


def label(c):
    """the report's name, plus the assigned molecule where the name says otherwise (siamang chr19_hap1 = 12)"""
    m = molecule.get(c, "")
    num = c.split("_")[0].lower().removeprefix("chr")
    same = not m or m == "na" or (num.lstrip("0") or num) == (m.lower().lstrip("0") or m.lower())
    return c if same else f"{c} (mol {m})"


# ---- colours ----
C_FAMILY, C_COFAM, C_PERIOD, C_OTHER = "#1d4ed8", "#0f766e", "#93c5fd", "#cbd5e1"
OTHERS = ["#f97316", "#a855f7", "#84cc16", "#b45309"]
C_LOCUS, C_PUB = "#0f172a", "#c026d3"
C_BODY, C_BODY_LINE, C_GRID, C_TEXT, C_MUTED = "#f1f5f9", "#94a3b8", "#e2e8f0", "#0f172a", "#64748b"
C_AGREE, C_CONTEST, C_NOCALL, C_NONE = "#16a34a", "#dc2626", "#f59e0b", "#94a3b8"
DASH = ' stroke-dasharray="2.5 1.5"'
NEAR = 3_000_000          # a published interval farther than this from our locus is an arrow in panel b


def period_label(b):
    p = int(b["base_period"])
    m = b["rep_mono_consensus"]
    return f"{p:,} bp" if p else f"({m if len(m) <= 12 else m[:10] + '…'})n"


chosen = [b for b in breadth if b["chosen"] == "1"]
if len(chosen) != 1:
    sys.exit(f"[{args.acc}] plot REFUSED: family_breadth.tsv has {len(chosen)} chosen superfamilies")
chosen = chosen[0]
cofams = [b for b in breadth if b["chosen"] == "2"]
rest = sorted((b for b in breadth if b["chosen"] == "0"), key=lambda b: (-int(b["total_bp"]), b["family"]))
other_colour = {b["family"]: OTHERS[i] for i, b in enumerate(rest[:len(OTHERS)])}


def colour(a):
    if a["membership"] == "family":
        return C_FAMILY
    if a["membership"] == "cofamily":
        return C_COFAM
    if a["membership"].startswith("period_class"):
        return C_PERIOD
    return other_colour.get(a["family"], C_OTHER)


sat_items = [("swatch", C_FAMILY, f"{chosen['family']} centromeric superfamily, {period_label(chosen)}, {size(chosen['total_bp'])}")]
sat_items += [("swatch", C_COFAM, f"{b['family']} co-family, {period_label(b)}, {size(b['total_bp'])}") for b in cofams]
if any(a["membership"].startswith("period_class") for a in arrays):
    sat_items.append(("swatch", C_PERIOD, f"joined by monomer period (k × {chosen['base_period']} bp)"))
sat_items += [("swatch", other_colour[b["family"]], f"{b['family']} {period_label(b)}, {size(b['total_bp'])}")
              for b in rest[:len(OTHERS)]]
sat_items.append(("swatch", C_OTHER, "other satellite arrays"))
mark_items = [("locus", C_LOCUS, "our centromeric locus"), ("largest", C_LOCUS, "its largest array (candidate active array)")]
for i, s in enumerate(sets):
    mark_items.append(("pub", i, f"published{' (dashed)' if i else ''}: {s.replace(':', ' ')}"))
if sets:
    mark_items += [("agree", "agrees", "agrees with the published annotation"),
                   ("agree", "contested_locus", "contested: no overlap"),
                   ("agree", "contested_active_array", "contested: published active array is not our largest array"),
                   ("agree", "contested_no_call", "contested: published centromere, no call"),
                   ("agree", "not_compared", "nothing published to compare")]
mark_items.append(("nocall", "#dc2626", "no array of the centromeric family"))

# ---- header text ----
n = len(chroms)
n_called = sum(1 for c, _ in chroms if c in loci)
sub = [f"{args.acc} · centromeric locus called on {n_called} of {n} chromosomes"]
for s in sets:
    ag = defaultdict(int)
    for r in overlap:
        if r["set"] == s:
            ag[r["agreement"]] += 1
    comp = sum(v for k, v in ag.items() if k != "not_compared")
    con = comp - ag["agrees"]
    detail = ", ".join(f"{t} {ag[k]}" for k, t in (("contested_locus", "no overlap"), ("contested_active_array",
                       "other active array"), ("contested_no_call", "no call")) if ag[k])
    sub.append(f"published {s.replace(':', ' ')}: " + (f"agrees on {ag['agrees']} of {comp} chromosomes"
               + (f", contested on {con} ({detail})" if con else "") if comp else "no interval on these chromosomes"))
if not sets:
    sub.append("no published centromere intervals for this assembly")

# ---- geometry ----
pitch = 44 if n <= 24 else (40 if n <= 28 else 36)
bar_w, left, right = 14, 80, 40
width = max(left + n * pitch + right, 1000)
max_len = max(L for _, L in chroms)
max_h = 480
scale = max_h / max_len
cw = 16                      # panel b bar width
hb = 230                     # panel b bar height
label_rot = any(len(label(c)) > 5 for c, _ in chroms)
lab_h = 96 if label_rot else 26


def wrap(items, x0):
    out_rows, row, x = [], [], x0
    for it in items:
        w = 24 + 6.7 * len(it[2]) + 20
        if row and x + w > width - right:
            out_rows.append(row)
            row, x = [], x0
        row.append((x, it))
        x += w
    out_rows.append(row)
    return out_rows


lx = left + 92
sat_rows, mark_rows = wrap(sat_items, lx), wrap(mark_items, lx)
y_sub = 62
y_leg = y_sub + 18 * (len(sub) - 1) + 34
y_marks = y_leg + 20 * len(sat_rows) + 6
y_a = y_marks + 20 * len(mark_rows) + 22
top_a = y_a + (40 if sets else 22)
y_b = top_a + max_h + lab_h + 24
top_b = y_b + (60 if sets else 46)
height = int(top_b + hb + 22 + lab_h + 48)

out = []
w_ = out.append
w_(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
   f'font-family="Helvetica Neue, Helvetica, Arial, sans-serif">')
w_(f'<title>{esc(args.organism)} {esc(args.asm_name)} {esc(args.acc)}: centromeric satellite arrays</title>')
w_(f'<rect width="{width}" height="{height}" fill="#ffffff"/>')

# ---- header ----
w_(f'<text x="{left - 40}" y="38" font-size="22" fill="{C_TEXT}"><tspan font-style="italic" font-weight="bold">'
   f'{esc(args.organism)}</tspan><tspan dx="10" font-weight="bold">{esc(args.asm_name)}</tspan></text>')
for k, line in enumerate(sub):
    w_(f'<text x="{left - 40}" y="{y_sub + 18 * k}" font-size="13" fill="#334155">{esc(line)}</text>')


def dot(x, y, agreement, title=""):
    t = f"<title>{esc(title)}</title>" if title else ""
    if agreement == "agrees":
        return f'<circle cx="{f1(x)}" cy="{f1(y)}" r="4.5" fill="{C_AGREE}">{t}</circle>'
    if agreement == "contested_locus":
        return f'<circle cx="{f1(x)}" cy="{f1(y)}" r="4.5" fill="{C_CONTEST}">{t}</circle>'
    if agreement == "contested_active_array":
        return f'<circle cx="{f1(x)}" cy="{f1(y)}" r="3.6" fill="#ffffff" stroke="{C_CONTEST}" stroke-width="2">{t}</circle>'
    if agreement == "contested_no_call":
        return f'<circle cx="{f1(x)}" cy="{f1(y)}" r="4.5" fill="{C_NOCALL}">{t}</circle>'
    return f'<circle cx="{f1(x)}" cy="{f1(y)}" r="4" fill="#ffffff" stroke="{C_NONE}" stroke-width="1.3">{t}</circle>'


def agreement_dots(i, c, y):
    for si, s in enumerate(sets):
        o = rel.get((c, s))
        agreement = o["agreement"] if o else "not_compared"
        detail = f"{o['relation']}, largest array overlap {o['largest_array_overlap_bp'] or 0} bp" if o else "nothing published"
        w_(dot(col_x(i) + bar_w / 2 + (si - (len(sets) - 1) / 2) * 11, y, agreement,
               f"{label(c)} vs {s}: {agreement} ({detail})"))


def legend_rows(rws, y0, title):
    w_(f'<text x="{left - 40}" y="{y0}" font-size="11" font-weight="bold" fill="{C_MUTED}" letter-spacing="0.6">{title}</text>')
    for ri, row in enumerate(rws):
        y = y0 + 20 * ri
        for x, (kind, col, text) in row:
            if kind == "swatch":
                w_(f'<rect x="{f1(x)}" y="{y - 9}" width="16" height="10" rx="2" fill="{col}"/>')
            elif kind == "locus":
                w_(f'<rect x="{f1(x + 2)}" y="{y - 11}" width="12" height="14" rx="2" fill="none" stroke="{col}" stroke-width="1.5"/>')
            elif kind == "largest":
                w_(f'<path d="M{f1(x + 14)},{y - 4} l-8,-5 v10 z" fill="{col}"/>')
            elif kind == "pub":
                w_(f'<path d="M{f1(x + 5)},{y - 11} h5 V{y + 3} h-5" fill="none" stroke="{C_PUB}" stroke-width="2"{DASH if col else ""}/>')
            elif kind == "agree":
                w_(dot(x + 8, y - 4, col))
            else:
                w_(f'<rect x="{f1(x + 4)}" y="{y - 11}" width="8" height="14" rx="4" fill="{C_BODY}" stroke="{col}" stroke-width="1"{DASH}/>')
            w_(f'<text x="{f1(x + 22)}" y="{y}" font-size="12" fill="#334155">{esc(text)}</text>')


legend_rows(sat_rows, y_leg, "SATELLITES")
legend_rows(mark_rows, y_marks, "CALLS")

arrays_by_chrom = defaultdict(list)
for a in arrays:
    arrays_by_chrom[a["chromosome"]].append(a)


def col_x(i):
    return left + i * pitch + (pitch - bar_w) / 2


def chrom_label(i, c, y):
    cx = col_x(i) + bar_w / 2
    lc = C_TEXT if c in loci else "#dc2626"
    if label_rot:
        w_(f'<text x="{f1(cx + 4)}" y="{y}" font-size="11" fill="{lc}" text-anchor="end" '
           f'transform="rotate(-60 {f1(cx + 4)} {y})">{esc(label(c))}</text>')
    else:
        w_(f'<text x="{f1(cx)}" y="{y}" font-size="12" fill="{lc}" text-anchor="middle">{esc(label(c))}</text>')


def array_rect(a, x, w, y0, sc, lo=0, hi=None):
    s, e = int(a["start"]), int(a["end"])
    if hi is not None:
        s, e = max(s, lo), min(e, hi)
        if e <= s:
            return
    hh = max(1.0, (e - s) * sc)
    w_(f'<rect x="{f1(x)}" y="{f1(y0 + (s - lo) * sc)}" width="{w}" height="{f1(hh)}" fill="{colour(a)}">'
       f'<title>{esc(a["chromosome"])} {mb(a["start"])}-{mb(a["end"])} Mb ({size(a["length"])}), {esc(a["family"])}, '
       f'monomer {esc(a["mono_period"])} bp, HOR {esc(a["hor_period"])} bp</title></rect>')


def locus_title(c, r):
    return (f'{esc(label(c))} locus {mb(r["locus_start"])}-{mb(r["locus_end"])} Mb ({size(r["locus_span"])}): '
            f'{size(r["family_bp"])} of {esc(r["family"])} in {r["n_arrays"]} arrays; largest array '
            f'{mb(r["largest_start"])}-{mb(r["largest_end"])} Mb ({size(r["largest_len"])}), monomer {esc(r["mono_period"])} bp')


# ---- panel a: whole chromosomes ----
w_(f'<text x="{left - 40}" y="{y_a}" font-size="14" fill="{C_TEXT}"><tspan font-weight="bold">a</tspan>'
   f'<tspan dx="10">Whole chromosomes, one scale</tspan></text>')
step = next(s for s in (5, 10, 20, 25, 50, 100, 200) if max_len / 1e6 / s <= 8) * 1_000_000
t = 0
while t <= max_len:
    y = top_a + t * scale
    w_(f'<line x1="{left - 8}" y1="{f1(y)}" x2="{width - right + 4}" y2="{f1(y)}" stroke="{C_GRID}" stroke-width="1"/>')
    w_(f'<text x="{left - 14}" y="{f1(y + 4)}" font-size="10" fill="{C_MUTED}" text-anchor="end">{t // 1_000_000}</text>')
    t += step
w_(f'<text x="{left - 14}" y="{top_a - 10}" font-size="10" fill="{C_MUTED}" text-anchor="end">Mb</text>')

for i, (c, L) in enumerate(chroms):
    cx = col_x(i)
    called = c in loci
    tip = [f"{label(c)}: {L:,} bp"] + ([] if called else ["no array of the centromeric family: no call"])
    tip += [f"exception: {k}" for k in sorted(set(exc_by_chrom.get(c, [])))]
    w_(f'<rect x="{f1(cx)}" y="{top_a}" width="{bar_w}" height="{f1(L * scale)}" rx="7" fill="{C_BODY}" '
       f'stroke="{"#dc2626" if not called else C_BODY_LINE}" stroke-width="1"{"" if called else DASH}>'
       f'<title>{esc(chr(10).join(tip))}</title></rect>')
    for a in sorted(arrays_by_chrom.get(c, []), key=lambda a: (int(a["length"]), int(a["start"]))):
        array_rect(a, cx + 1, bar_w - 2, top_a, scale)
    if called:
        r = loci[c]
        s, e = int(r["locus_start"]), int(r["locus_end"])
        hh = max(5.0, (e - s) * scale)
        y0 = top_a + (s + e) / 2 * scale - hh / 2
        w_(f'<rect x="{f1(cx - 3)}" y="{f1(y0)}" width="{bar_w + 6}" height="{f1(hh)}" rx="2" fill="none" '
           f'stroke="{C_LOCUS}" stroke-width="1.5"><title>{locus_title(c, r)}</title></rect>')
        ym = top_a + (int(r["largest_start"]) + int(r["largest_end"])) / 2 * scale
        w_(f'<path d="M{f1(cx - 5)},{f1(ym)} l-7,-4.5 v9 z" fill="{C_LOCUS}"/>')
    for si, s in enumerate(sets):
        p = published(c, s)
        if p:
            x0 = cx + bar_w + 5 + 5 * si
            y1, y2 = top_a + p[0] * scale, max(top_a + p[0] * scale + 3, top_a + p[1] * scale)
            w_(f'<path d="M{f1(x0)},{f1(y1)} h4 V{f1(y2)} h-4" fill="none" stroke="{C_PUB}" stroke-width="2"{DASH if si else ""}>'
               f'<title>{esc(s)}: {mb(p[0])}-{mb(p[1])} Mb; {p[2]}</title></path>')
    agreement_dots(i, c, top_a - 16)
    chrom_label(i, c, top_a + max_h + 16)

# ---- panel b: centromeric region close-ups ----
w_(f'<text x="{left - 40}" y="{y_b}" font-size="14" fill="{C_TEXT}"><tspan font-weight="bold">b</tspan>'
   f'<tspan dx="10">Centromeric region of each chromosome, each window to its own scale (width shown above it)</tspan></text>')
for i, (c, L) in enumerate(chroms):
    bx = col_x(i) + bar_w / 2 - cw / 2
    agreement_dots(i, c, top_b - 40)
    if c not in loci:
        w_(f'<rect x="{f1(bx)}" y="{top_b}" width="{cw}" height="{hb}" rx="3" fill="#ffffff" stroke="#dc2626" '
           f'stroke-width="1"{DASH}><title>{esc(label(c))}: no call</title></rect>')
        w_(f'<text x="{f1(bx + cw / 2)}" y="{f1(top_b + hb / 2)}" font-size="10" fill="#dc2626" text-anchor="middle" '
           f'transform="rotate(-90 {f1(bx + cw / 2)} {f1(top_b + hb / 2)})">no call</text>')
        chrom_label(i, c, top_b + hb + 30)
        continue
    r = loci[c]
    s0, e0 = int(r["locus_start"]), int(r["locus_end"])
    lo, hi, far = s0, e0, []
    for si, s in enumerate(sets):
        p = published(c, s)
        if not p:
            continue
        if p[0] - e0 > NEAR or s0 - p[1] > NEAR:
            far.append((si, s, p))
        else:
            lo, hi = min(lo, p[0]), max(hi, p[1])
    pad = max(0.08 * (hi - lo), 20_000)
    lo, hi = max(0, int(lo - pad)), min(L, int(hi + pad))
    sc = hb / (hi - lo)
    w_(f'<text x="{f1(bx + cw / 2)}" y="{top_b - 22}" font-size="9.5" fill="{C_MUTED}" text-anchor="middle">{short_size(hi - lo)}</text>')
    w_(f'<rect x="{f1(bx)}" y="{top_b}" width="{cw}" height="{hb}" rx="2" fill="{C_BODY}" stroke="{C_BODY_LINE}" stroke-width="1">'
       f'<title>{esc(label(c))} {mb(lo)}-{mb(hi)} Mb</title></rect>')
    for a in sorted(arrays_by_chrom.get(c, []), key=lambda a: (int(a["length"]), int(a["start"]))):
        if int(a["end"]) > lo and int(a["start"]) < hi:
            array_rect(a, bx + 1, cw - 2, top_b, sc, lo, hi)
    y1, y2 = top_b + (s0 - lo) * sc, top_b + (e0 - lo) * sc
    w_(f'<path d="M{f1(bx - 2)},{f1(y1)} h-4 V{f1(max(y2, y1 + 2))} h4" fill="none" stroke="{C_LOCUS}" stroke-width="1.8">'
       f'<title>{locus_title(c, r)}</title></path>')
    ym = top_b + ((int(r["largest_start"]) + int(r["largest_end"])) / 2 - lo) * sc
    w_(f'<path d="M{f1(bx - 8)},{f1(ym)} l-7,-4.5 v9 z" fill="{C_LOCUS}"/>')
    for si, s in enumerate(sets):
        p = published(c, s)
        if not p or (si, s, p) in far:
            continue
        x0 = bx + cw + 3 + 5 * si
        p1, p2 = top_b + (max(p[0], lo) - lo) * sc, top_b + (min(p[1], hi) - lo) * sc
        w_(f'<path d="M{f1(x0)},{f1(p1)} h4 V{f1(max(p2, p1 + 2))} h-4" fill="none" stroke="{C_PUB}" stroke-width="2"{DASH if si else ""}>'
           f'<title>{esc(s)}: {mb(p[0])}-{mb(p[1])} Mb; {p[2]}</title></path>')
    for k, (si, s, p) in enumerate(far):
        up = p[1] < s0
        dist = (s0 - p[1]) if up else (p[0] - e0)
        ty = top_b - 8 if up else top_b + hb + 14
        w_(f'<text x="{f1(bx + cw / 2)}" y="{f1(ty + (-10 if up else 10) * k)}" font-size="9.5" font-weight="bold" '
           f'fill="{C_PUB}" text-anchor="middle"><title>{esc(s)}: {mb(p[0])}-{mb(p[1])} Mb, {short_size(dist)} '
           f'{"before" if up else "after"} our locus; {p[2]}</title>{"↑" if up else "↓"}{dist / 1e6:.0f}</text>')
    chrom_label(i, c, top_b + hb + 30)

w_(f'<text x="{left - 40}" y="{height - 28}" font-size="10" fill="#94a3b8">Satellite arrays ≥ 10 kb from the '
   f'satellome (FasTAN) catalogue, coloured by sequence superfamily. Published intervals are compared, never used as '
   f'input.</text>')
w_(f'<text x="{left - 40}" y="{height - 14}" font-size="10" fill="#94a3b8">Arrows in b: a published interval more '
   f'than 3 Mb from our locus, distance in Mb. Hover for coordinates. github.com/aglabx/t2t_centromeres</text>')
w_('</svg>')
with open(args.out, "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
print(f"[{args.acc}] karyotype chromosomes={n} called={n_called} arrays={len(arrays)} sets={len(sets)}")
