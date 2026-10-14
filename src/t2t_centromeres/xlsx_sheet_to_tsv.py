"""Dump one sheet of an .xlsx to TSV with the standard library only (zipfile + ElementTree).
`--list` prints the sheet names. Cells are written as stored (shared strings resolved, numbers raw).
Usage: python -I xlsx_sheet_to_tsv.py BOOK.xlsx --list
       python -I xlsx_sheet_to_tsv.py BOOK.xlsx --sheet NAME > out.tsv
"""
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
book = sys.argv[1]
z = zipfile.ZipFile(book)
wb = ET.fromstring(z.read("xl/workbook.xml"))
rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
target = {r.get("Id"): r.get("Target") for r in rels}
sheets = [(s.get("name"), target[s.get(f"{{{NS['r']}}}id")]) for s in wb.find("m:sheets", NS)]
if "--list" in sys.argv:
    for n, t in sheets:
        print(f"{n}\t{t}")
    raise SystemExit(0)
name = sys.argv[sys.argv.index("--sheet") + 1]
path = dict(sheets).get(name)
if path is None:
    raise SystemExit(f"no sheet {name!r}; have {[n for n, _ in sheets]}")
path = "xl/" + path.lstrip("/").removeprefix("xl/")
shared = []
if "xl/sharedStrings.xml" in z.namelist():
    for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
        shared.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))


def col_index(ref):
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n - 1


for row in ET.fromstring(z.read(path)).find("m:sheetData", NS):
    cells = {}
    for c in row:
        v = c.find("m:v", NS)
        if c.get("t") == "s" and v is not None:
            val = shared[int(v.text)]
        elif c.get("t") == "inlineStr":
            val = "".join(t.text or "" for t in c.iter(f"{{{NS['m']}}}t"))
        else:
            val = v.text if v is not None else ""
        cells[col_index(c.get("r"))] = val.replace("\t", " ").replace("\n", " ")
    if cells:
        print("\t".join(cells.get(i, "") for i in range(max(cells) + 1)))
