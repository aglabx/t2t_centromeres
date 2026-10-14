"""Normalise the fetched published centromere annotations to BED4 (chrom, start0, end, label)
with chromosome names that build_gff.py can join without guessing.
Downloaded files are data only; this script reads them, it does not execute them.

  siamang  CenSat v2.0 on v2.1   names -> GenBank accession through the T2T-Browser alias table
           (NCBI's Sequence-Name carries the pre-v2.1 chr12/chr19 labels; never join on it)
  MFA8     Supp. Table 39 presumed-active arrays; release centromere.bed (alpha-sat span)
  MMU8     release alpha_satellite.bed (RepeatMasker alpha-sat interval per chromosome)
  sheep    Supp. Table 8 genes "in the centromeric regions" (gene intervals, 1-based in the table)
Usage: python -I prepare_official.py OFFICIAL_DIR [NCBI_REPORT ...]   (reads OFFICIAL_DIR/raw, writes OFFICIAL_DIR/norm)
  NCBI reports fill CenSat names the alias table lacks (2026-10-14: chr24_hap2 only), by NCBI
  Sequence-Name, never for chr12_* / chr19_* (the labels NCBI carries swapped); each fill is listed.
"""
import os
import sys

D0 = sys.argv[1]
D = os.path.join(D0, "raw")
OUT = os.path.join(D0, "norm")
os.makedirs(OUT, exist_ok=True)
report = []


def write(name, rows):
    with open(os.path.join(OUT, name), "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    report.append(f"{name}\t{len(rows)}")


# siamang
alias = {}
with open(os.path.join(D, "mSymSyn1_v2.1.alias")) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) > 2 and f[2]:
            alias[f[0]] = f[2]
filled = {}
for rep in sys.argv[2:]:
    with open(rep) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            name = f[0]
            if (f[1] == "assembled-molecule" and name not in alias and f[4] != "na"
                    and not name.startswith(("chr12_", "chr19_"))):
                alias[name] = f[4]
                filled[name] = f[4]
report.append(f"siamang_alias_filled_from_ncbi\t{len(filled)}\t{','.join(f'{k}={v}' for k, v in filled.items())}")
rows, unmapped = [], set()
with open(os.path.join(D, "SymSyn.cenSatv2.0.bed")) as fh:
    for ln in fh:
        if ln.startswith(("track", "#")):
            continue
        f = ln.rstrip("\n").split("\t")
        if f[0] in alias:
            rows.append((alias[f[0]], f[1], f[2], f[3]))
        else:
            unmapped.add(f[0])
write("siamang_cenSat_v2.1.by_genbank.bed", rows)
report.append(f"siamang_unmapped_names\t{len(unmapped)}\t{','.join(sorted(unmapped)[:5])}")

# MFA8
rows = []
with open(os.path.join(D, "MFA8_S39.tsv")) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if f[0].startswith("chr") and f[1].isdigit():
            rows.append((f[0], f[1], f[2], "presumed_active_array(S39)"))
write("MFA8_S39_presumed_active.bed", rows)
rows = []
with open(os.path.join(D, "T2T-MFA8v1.1.centromere.bed")) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        rows.append((f[0], f[1], f[2], "centromere_alpha_sat_span"))
write("MFA8_centromere_alpha_span.bed", rows)

# MMU8
rows = []
with open(os.path.join(D, "T2T-MMU8v2.0.alpha_satellite.bed")) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        rows.append((f[0], f[1], f[2], "alpha_satellite_RepeatMasker"))
write("MMU8_alpha_satellite.bed", rows)

# sheep
rows = []
with open(os.path.join(D, "sheep_S8.tsv")) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if f[0].startswith("Chr") and f[1].isdigit():
            rows.append((f[0], int(f[1]) - 1, f[2], f"gene_in_centromeric_region({f[3]})"))
write("sheep_S8_centromeric_genes.bed", rows)

print("\n".join(report))
