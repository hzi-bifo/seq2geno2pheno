#!/usr/bin/env python3
"""Check the s2g2p ports against a finished Seq2Geno denovo run (the published pipeline).

Each check recomputes one published intermediate from that run's own upstream files and
compares:

  kmers     count_canonical(spades/<s>/contigs.fasta)    vs  KMC dump merge_6mers_<id>.txt
  gpa       gpa_from_roary(roary/gene_presence_absence)  vs  gpa.mat
  compress  compress(gpa.mat), compress(indel.mat)       vs  *_NONRDNT and *_GROUPS
  families  select_families(legacy=True)                 vs  the families with an indels/<f>.vcf
  vcf2indel vcf2indel(indels/<f>.vcf)                    vs  indels/<f>_indels.txt
  indel     indel_matrix(...)                            vs  indel.mat

Usage (a SLURM job; the run's tables are tens of MB):
  python dev/validate_reproduction.py --s2g-dir <software_output/<species>> \
      --kmc-dir <software_output/cano6mer> --out validation.tsv
"""

from __future__ import annotations

import argparse
import random
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from s2g2p.features.compress import compress  # noqa: E402
from s2g2p.features.gpa import gpa_from_roary, read_rename_map  # noqa: E402
from s2g2p.features.indel import indel_matrix, select_families, vcf2indel  # noqa: E402
from s2g2p.features.kmer import count_canonical, read_kmc_dump  # noqa: E402

ROWS = []


def report(check, ok, detail=""):
    ROWS.append({"check": check, "result": "identical" if ok else "DIFFERENT", "detail": detail})
    print(f"[{'OK ' if ok else 'BAD'}] {check}: {detail}", flush=True)


def read_table(path):
    return pd.read_csv(path, sep="\t", index_col=0, dtype=str, keep_default_na=False, na_values=[""])


def check_kmers(den, kmc_dir, strains, n):
    rng = random.Random(0)
    pick = strains if n <= 0 else rng.sample(strains, min(n, len(strains)))
    bad = []
    for s in pick:
        mine = count_canonical(den / "spades" / s / "contigs.fasta", 6)
        theirs = read_kmc_dump(kmc_dir / f"merge_6mers_{s.replace('iso_', '')}.txt", 6)
        if not np.array_equal(mine, theirs):
            bad.append((s, int(np.abs(mine - theirs).sum())))
    report("kmers", not bad, f"{len(pick) - len(bad)}/{len(pick)} genomes identical" + (f"; e.g. {bad[:3]}" if bad else ""))


def check_gpa(den, strains):
    rename = read_rename_map(den / "roary" / "roary_abricate.txt")
    mine = gpa_from_roary(den / "roary" / "gene_presence_absence.csv", strains, rename)
    theirs = read_table(den / "gpa.mat")
    same_cols = list(mine.columns) == list(theirs.columns)
    same_rows = list(mine.index) == list(theirs.index)
    same_vals = same_cols and same_rows and (mine.astype(str).to_numpy() == theirs.to_numpy()).all()
    report("gpa", same_vals, f"{theirs.shape}; columns {'same' if same_cols else 'differ'}, rows {'same' if same_rows else 'differ'}")
    return rename


def check_compress(den):
    for name in ("gpa.mat", "indel.mat"):
        src = read_table(den / name)
        out, groups = compress(src)
        ref = read_table(den / f"{name}_NONRDNT")
        ok = list(out.columns) == list(ref.columns) and (out.astype(str).to_numpy() == ref.to_numpy()).all()
        ref_groups = {}
        with open(den / f"{name}_GROUPS") as fh:
            for line in fh:
                k, _, v = line.rstrip("\n").partition(":\t")
                ref_groups[k] = v.split("\t")
        report(f"compress {name}", ok and groups == ref_groups,
               f"{src.shape[1]} -> {out.shape[1]} (published {ref.shape[1]}); groups {'same' if groups == ref_groups else 'differ'}")


def check_indels(den, rename):
    sel = select_families(den / "roary", min_isolates=2, legacy=True)
    vcfs = {p.stem for p in (den / "indels").glob("*.vcf")}
    report("families aligned", set(sel.groups) == vcfs,
           f"{len(sel.aligned)} families in {len(sel.groups)} alignments, {len(vcfs)} published vcf")

    calls, diff = {}, []
    for key in sorted(vcfs):
        mine = vcf2indel(den / "indels" / f"{key}.vcf")
        ref = pd.read_csv(den / "indels" / f"{key}_indels.txt", sep="\t", index_col=0).iloc[:, 0]
        if {str(k): int(v) for k, v in ref.items()} != mine:
            diff.append(key)
        calls[key] = mine
    report("vcf2indel", not diff, f"{len(vcfs) - len(diff)}/{len(vcfs)} families identical" + (f"; e.g. {diff[:3]}" if diff else ""))

    gpa_raw = gpa_from_roary(den / "roary" / "gene_presence_absence.csv",
                             list(read_table(den / "gpa.mat").index))
    mine = indel_matrix(gpa_raw, sel.candidates, calls, rename=rename)
    ref = read_table(den / "indel.mat")
    same_cols = list(mine.columns) == list(ref.columns)
    ok = same_cols and list(mine.index) == list(ref.index) and (mine.astype(str).to_numpy() == ref.to_numpy()).all()
    detail = f"{mine.shape[1]} columns (published {ref.shape[1]})"
    if not same_cols:
        a, b = set(mine.columns), set(ref.columns)
        detail += f"; only ours {len(a - b)} e.g. {sorted(a - b)[:3]}, only published {len(b - a)} e.g. {sorted(b - a)[:3]}"
        if a == b:
            detail += "; same set, different order"
    report("indel.mat", ok, detail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s2g-dir", required=True, type=Path)
    ap.add_argument("--kmc-dir", type=Path)
    ap.add_argument("--n-kmer-genomes", type=int, default=0, help="0 = all")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    den = a.s2g_dir / "results" / "denovo"
    strains = [l.split("\t")[0] for l in open(a.s2g_dir / "assemble_list") if l.strip()]
    t0 = time.time()
    rename = check_gpa(den, strains)
    check_compress(den)
    check_indels(den, rename)
    if a.kmc_dir:
        check_kmers(den, a.kmc_dir, strains, a.n_kmer_genomes)
    pd.DataFrame(ROWS).to_csv(a.out, sep="\t", index=False)
    print(f"done in {time.time() - t0:.0f}s -> {a.out}")
    sys.exit(0 if all(r["result"] == "identical" for r in ROWS) else 1)


if __name__ == "__main__":
    main()
