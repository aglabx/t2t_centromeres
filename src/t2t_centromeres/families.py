"""Small helpers shared by the caller: reverse complement, rotation/strand-canonical form of a
monomer (Booth's least rotation), union-find, and a subprocess runner that fails loudly."""
import subprocess
import sys

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def least_rotation(s):
    """Booth's algorithm: lexicographically minimal rotation, O(n)."""
    s2 = s + s
    n = len(s2)
    f = [-1] * n
    k = 0
    for j in range(1, n):
        sj = s2[j]
        i = f[j - k - 1]
        while i != -1 and sj != s2[k + i + 1]:
            if sj < s2[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if sj != s2[k + i + 1]:
            if sj < s2[k]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return s2[k:k + len(s)]


def canonical(s):
    s = s.upper()
    return min(least_rotation(s), least_rotation(rc(s)))


class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.stderr.write(p.stderr[-2000:])
        raise SystemExit(f"command failed ({p.returncode}): {' '.join(cmd)}")
