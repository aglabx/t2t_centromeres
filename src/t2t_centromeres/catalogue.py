"""Read one satellome output directory (FasTAN layout) without guessing.

The analysis tree holds catalogues in three encodings: plain files, gzip, and the
`satellome compact` columnar container (.satz). Every reader here fails loudly when the
file it needs is missing or unreadable; there is no fallback to a neighbouring file.
"""
import glob
import gzip
import io
import os
import shutil
import subprocess


class CatalogueError(Exception):
    pass


def _one(pattern, what):
    hits = sorted(glob.glob(pattern))
    if len(hits) != 1:
        raise CatalogueError(f"{what}: expected exactly one match for {pattern}, got {hits}")
    return hits[0]


def _text_lines(path):
    if path.endswith(".satz"):
        # catalogues compacted by `satellome compact`; needs the satellome package
        from satellome.compact.columnar import restore_table
        buf = io.BytesIO()
        restore_table(path, buf)
        return buf.getvalue().decode().splitlines()
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return fh.read().splitlines()
    with open(path) as fh:
        return fh.read().splitlines()


def bed_path(run_dir):
    return _one(os.path.join(run_dir, "fastan", "*.bed*"), "FasTAN BED")


def summary_path(run_dir):
    return _one(os.path.join(run_dir, "fastan", "*.summary.tsv*"), "arraysplitter summary")


def read_bed(run_dir):
    """FasTAN/tanbed BED: seq, start, end, period, score. Returns list of tuples."""
    rows = []
    for ln in _text_lines(bed_path(run_dir)):
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 4:
            raise CatalogueError(f"malformed BED row in {run_dir}: {ln!r}")
        rows.append((f[0], int(f[1]), int(f[2]), int(f[3]), f[4] if len(f) > 4 else ""))
    if not rows:
        raise CatalogueError(f"empty BED in {run_dir}")
    return rows


def read_summary(run_dir, min_len=0):
    """arraysplitter summary keyed by array_id; only arrays >= min_len kept."""
    lines = _text_lines(summary_path(run_dir))
    head = lines[0].split("\t")
    col = {c: i for i, c in enumerate(head)}
    need = ["array_id", "array_length", "orientation", "hor_period", "hor_consensus",
            "mono_period", "mono_consensus"]
    missing = [c for c in need if c not in col]
    if missing:
        raise CatalogueError(f"summary in {run_dir} lacks columns {missing}")
    out = {}
    for ln in lines[1:]:
        f = ln.split("\t")
        if int(f[col["array_length"]]) < min_len:
            continue
        out[f[col["array_id"]]] = {k: f[col[k]] for k in need}
    return out


def array_id(seq, start, end, period):
    """summary.tsv keys arrays as <seq>_<start>_<end>_<length>_<period>."""
    return f"{seq}_{start}_{end}_{end - start}_{period}"


def twobit_sizes(path):
    tool = shutil.which("twoBitInfo")
    if tool is None:
        raise CatalogueError("twoBitInfo (UCSC kent tools) is not on PATH")
    out = subprocess.run([tool, path, "stdout"], check=True,
                         capture_output=True, text=True).stdout
    sizes = {}
    for ln in out.splitlines():
        name, size = ln.split("\t")
        sizes[name] = int(size)
    if not sizes:
        raise CatalogueError(f"twoBitInfo returned nothing for {path}")
    return sizes
