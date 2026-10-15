#!/usr/bin/env python3
"""Compare two Roary runs as partitions of the same proteins.

Roary is not deterministic (two runs on identical input share ~98.9% of family member sets and
renumber `group_N`), so pangenomes are compared by member sets, never by bytes or names.
Locus tags are normalised by stripping a prefix (the reproduction's Prokka ran with
`--locustag iso_<id>`; v3 uses `<id>`).

  python dev/compare_pangenome.py A/clustered_proteins B/clustered_proteins --strip-a iso_ --out cmp.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from s2g2p.features.gpa import read_clustered_proteins  # noqa: E402


def partition(path, strip):
    fams = read_clustered_proteins(path)
    return {frozenset(m[len(strip):] if strip and m.startswith(strip) else m for m in mem) for mem in fams.values()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--strip-a", default="")
    ap.add_argument("--strip-b", default="")
    ap.add_argument("--out")
    x = ap.parse_args()
    A, B = partition(x.a, x.strip_a), partition(x.b, x.strip_b)
    pa, pb = set().union(*A), set().union(*B)
    shared = A & B
    row = {"families_a": len(A), "families_b": len(B), "identical_member_sets": len(shared),
           "frac_of_a": round(len(shared) / len(A), 4), "frac_of_b": round(len(shared) / len(B), 4),
           "proteins_a": len(pa), "proteins_b": len(pb), "proteins_both": len(pa & pb)}
    for k, v in row.items():
        print(f"{k}\t{v}")
    if x.out:
        with open(x.out, "w") as fh:
            fh.write("\t".join(row) + "\n" + "\t".join(map(str, row.values())) + "\n")


if __name__ == "__main__":
    main()
