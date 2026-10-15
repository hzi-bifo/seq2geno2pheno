#!/usr/bin/env python3
"""Find the Roary stage where two runs on the same proteins diverge.

Roary (run with -z) keeps its intermediates. Each stage is compared after normalising locus
tags (--strip-a/--strip-b remove a prefix, e.g. the reproduction's `iso_`):

  _combined_files       every protein (id -> sequence)
  _clustered.clstr      CD-HIT's clusters, as a partition of protein ids
  _clustered            CD-HIT's representatives
  _blast_results        all-vs-all hits between representatives, as (query, subject) pairs
  _uninflated_mcl_groups, _inflated_mcl_groups   MCL's groups, as partitions
  clustered_proteins    the final families, as a partition

  python dev/compare_roary_stages.py A/roary B/roary --strip-a iso_ --out stages.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def norm(strip):
    return (lambda s: s[len(strip):] if strip and s.startswith(strip) else s)


def fasta(path, n):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = n(line[1:].split()[0]), []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def clstr(path, n):
    groups, cur = [], None
    rx = re.compile(r">(\S+?)\.\.\.")
    with open(path) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                cur = set()
                groups.append(cur)
            else:
                m = rx.search(line)
                if m:
                    cur.add(n(m.group(1)))
    return {frozenset(g) for g in groups}


def lines_partition(path, n, sep=None, skip_label=False):
    out = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if skip_label and ":" in line.split()[0]:
                line = line.split(":", 1)[1]
            out.add(frozenset(n(x) for x in line.split(sep)))
    return out


def blast_pairs(path, n):
    with open(path) as fh:
        return {(n(f[0]), n(f[1])) for f in (l.split("\t") for l in fh) if len(f) > 1}


def compare(a, b):
    both = len(a & b)
    return {"a": len(a), "b": len(b), "identical": both,
            "frac_a": round(both / len(a), 4) if a else 1.0, "frac_b": round(both / len(b), 4) if b else 1.0}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", type=Path)
    ap.add_argument("b", type=Path)
    ap.add_argument("--strip-a", default="")
    ap.add_argument("--strip-b", default="")
    ap.add_argument("--out")
    x = ap.parse_args()
    na, nb = norm(x.strip_a), norm(x.strip_b)
    rows = []

    fa, fb = fasta(x.a / "_combined_files", na), fasta(x.b / "_combined_files", nb)
    rows.append(("_combined_files (id, sequence)", compare(set(fa.items()), set(fb.items()))))
    del fa, fb
    rows.append(("_clustered.clstr (CD-HIT clusters)", compare(clstr(x.a / "_clustered.clstr", na), clstr(x.b / "_clustered.clstr", nb))))
    rows.append(("_clustered (representatives)", compare(set(fasta(x.a / "_clustered", na)), set(fasta(x.b / "_clustered", nb)))))
    rows.append(("_blast_results (query, subject)", compare(blast_pairs(x.a / "_blast_results", na), blast_pairs(x.b / "_blast_results", nb))))
    for f in ("_uninflated_mcl_groups", "_inflated_mcl_groups"):
        if (x.a / f).exists() and (x.b / f).exists():
            rows.append((f, compare(lines_partition(x.a / f, na), lines_partition(x.b / f, nb))))
    rows.append(("clustered_proteins (families)", compare(lines_partition(x.a / "clustered_proteins", na, skip_label=True),
                                                          lines_partition(x.b / "clustered_proteins", nb, skip_label=True))))

    head = ["stage", "a", "b", "identical", "frac_a", "frac_b"]
    text = "\t".join(head) + "\n" + "".join("\t".join([s] + [str(r[k]) for k in head[1:]]) + "\n" for s, r in rows)
    print(text, end="")
    if x.out:
        Path(x.out).write_text(text)


if __name__ == "__main__":
    main()
