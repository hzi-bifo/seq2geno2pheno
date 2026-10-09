"""Indel features (Seq2Geno's ``indel_*`` rules), ported to Python 3.

Per gene family: align the members' nucleotide sequences (mafft), call variants against the
alignment (jvarkit ``msa2vcf``), and mark a sample 0 when one of its copies carries an indel of
at least 8 bp relative to the majority, 1 when it has the family without one, and 0 when it
lacks the family. A family becomes a column only if some sample carries an indel.

Which families are aligned, and which become columns, is where the published rules have
quirks. They are reproduced under ``legacy=True`` (the ``hu2024``, ``noleak`` and ``fast``
presets) because they change the columns, and each was confirmed against the reproduction's
own ``indel.mat`` on seven species:

* **Aligned** are the first N lines of Roary's ``clustered_proteins``, N being the number of
  families in more than ``min_isolates`` isolates (``head -n $core_genes_wc``) -- not those N
  families. Rare families that sort early are aligned; some common ones are not.
* **Alignment files** are named after the family with every non-word character replaced by
  ``_`` (gene_clusters2multi_fasta.py). Families whose names collide that way (``blaR1-2`` and
  ``blaR1_2``) share one file, and their sequences are aligned *together*.
* **Columns** are taken from every family in Roary's table order (the rule reads the ``.Rtab``,
  not the core list it declares as input). A family gets a column when a file of exactly its
  name exists: an aligned family whose name has a non-word character never does, a family whose
  sanitised twin was aligned does.

With ``legacy=False`` the core families themselves are aligned, each in its own file, and each
core family can become a column.

The external tools are taken from ``PATH``; the workflow's containers pin their versions
(mafft 7.310, jvarkit-msa2vcf 1.0 in the published environment).
"""

from __future__ import annotations

import collections
import csv
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pandas as pd

from ..io import read_fasta
from .gpa import read_clustered_proteins, read_roary_gpa

MIN_INDEL_LENGTH = 8
_NONWORD = re.compile(r"\W")


def file_key(family: str) -> str:
    """The name a family's alignment file gets in Seq2Geno (gene_clusters2multi_fasta.py)."""
    return _NONWORD.sub("_", family)


def unique_key(family: str) -> str:
    """A collision-free file name for ``legacy=False``."""
    return f"{file_key(family)}__{hashlib.md5(family.encode()).hexdigest()[:8]}"


@dataclass
class FamilySelection:
    groups: "collections.OrderedDict[str, list[str]]"   # alignment file stem -> member families
    candidates: list[tuple[str, str]]                   # (family, stem) that may become a column
    renamed: list[str] = field(default_factory=list)    # aligned families whose file is not their name

    @property
    def aligned(self) -> list[str]:
        return [f for fams in self.groups.values() for f in fams]


def table_order(roary_dir: Path) -> list[str]:
    """Families in the order of ``gene_presence_absence.Rtab`` (what the rule read), falling back
    to the CSV, which Roary writes in the same order."""
    rtab = roary_dir / "gene_presence_absence.Rtab"
    if rtab.is_file():
        with open(rtab) as fh:
            next(fh)
            return [l.split("\t", 1)[0].replace('"', "") for l in fh if l.strip()]
    return [str(f) for f in read_roary_gpa(roary_dir / "gene_presence_absence.csv").index]


def select_families(roary_dir: str | os.PathLike, min_isolates: int = 2, legacy: bool = True) -> FamilySelection:
    """Aligned families and indel-column candidates (``filter.awk`` with
    ``threshold=min_isolates``: strictly more isolates than the threshold)."""
    roary_dir = Path(roary_dir)
    gpa = read_roary_gpa(roary_dir / "gene_presence_absence.csv")
    n_iso = pd.to_numeric(gpa["No. isolates"], errors="raise")
    core = [str(f) for f in gpa.index[n_iso > min_isolates]]
    groups: "collections.OrderedDict[str, list[str]]" = collections.OrderedDict()
    if not legacy:
        for f in core:
            groups[unique_key(f)] = [f]
        return FamilySelection(groups, [(f, unique_key(f)) for f in core])
    head = list(read_clustered_proteins(roary_dir / "clustered_proteins"))[:len(core)]
    for f in head:
        groups.setdefault(file_key(f), []).append(f)
    candidates = [(f, f) for f in table_order(roary_dir) if f in groups]
    renamed = [f for f in head if file_key(f) != f]
    return FamilySelection(groups, candidates, renamed)


def write_family_fastas(clustered_proteins: str | os.PathLike, groups: Mapping[str, Sequence[str]],
                        ffn_files: Iterable[str | os.PathLike], out_dir: str | os.PathLike) -> list[Path]:
    """One multi-FASTA per alignment group, ``<stem>.fasta``, records in the order the ``.ffn``
    files and their records are read (gene_clusters2multi_fasta.py; it was given the files as
    a shell glob, so pass them sorted by path)."""
    fam2mem = read_clustered_proteins(clustered_proteins)
    fam2stem = {f: stem for stem, fams in groups.items() for f in fams}
    prot2stem: dict[str, str] = {}
    for fam, members in fam2mem.items():
        if fam in fam2stem:
            for p in members:
                prot2stem.setdefault(p, fam2stem[fam])
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seqs: dict[str, list[str]] = collections.defaultdict(list)
    for f in ffn_files:
        with open(f) as fh:
            header, chunk = None, []
            for line in fh:
                if line.startswith(">"):
                    if header is not None:
                        pid = header.split(" ")[0][1:].rstrip("\n")
                        if pid in prot2stem:
                            seqs[prot2stem[pid]].append(header + "".join(chunk))
                    header, chunk = line, []
                else:
                    chunk.append(line)
            if header is not None:
                pid = header.split(" ")[0][1:].rstrip("\n")
                if pid in prot2stem:
                    seqs[prot2stem[pid]].append(header + "".join(chunk))
    paths = []
    for stem in groups:
        if stem in seqs:
            p = out_dir / f"{stem}.fasta"
            p.write_text("".join(seqs[stem]))
            paths.append(p)
    return paths


def write_family_index(sel: FamilySelection, path: str | os.PathLike) -> None:
    """``families.tsv`` beside the alignments: which file holds which family, and which
    families are column candidates. Placement reads it to find a column's alignment."""
    cand = {f for f, _ in sel.candidates}
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["family", "stem", "candidate"])
        for stem, fams in sel.groups.items():
            for f in fams:
                w.writerow([f, stem, int(f in cand)])
        for f, stem in sel.candidates:
            if f not in {x for xs in sel.groups.values() for x in xs}:
                w.writerow([f, stem, 1])


def read_family_index(path: str | os.PathLike) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def sample_of(record_name: str) -> str:
    """Sample id from a Prokka locus tag ``<sample>_<number>`` (vcf2indel.py)."""
    return "_".join(record_name.split(" ")[0].split("_")[:-1])


def vcf2indel(vcf_path: str | os.PathLike) -> dict[str, int]:
    """Per-sample call of one family's msa2vcf output: 0 = carries an indel >= 8 bp,
    1 = does not. A line-by-line port of Seq2Geno's ``vcf2indel.py`` (Python 2)."""
    with open(vcf_path) as fh:
        line = ""
        for line in fh:
            if not line.startswith("##"):
                break
        header = line.strip("#").rstrip("\n").split("\t")
        names = [sample_of(c) for c in header[9:]]
        has = {n: 1 for n in names}
        half = len(names) / 2
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            ref, alts = cols[3], cols[4].split(",")
            gts = [c.split("/")[0] for c in cols[9:]]
            # "." is a gap. When most samples have the gap, the gap is the majority allele and
            # the samples carrying sequence are the ones with the indel.
            dot_majority = sum(1 for g in gts if g == ".") > half
            for i, alt in enumerate(alts, start=1):
                if alt == "" or abs(len(alt) - len(ref)) < MIN_INDEL_LENGTH:
                    continue
                if ref == "." or alt == ".":
                    hits = [n for n, g in zip(names, gts) if (g != "." if dot_majority else g == ".")]
                elif dot_majority:
                    hits = [n for n, g in zip(names, gts) if g == str(i)]
                else:
                    hits = [n for n, g in zip(names, gts) if g == str(i) or g == "."]
                for h in hits:
                    has[h] = 0
    return has


def indel_matrix(gpa: pd.DataFrame, candidates: Sequence[tuple[str, str]],
                 calls: Mapping[str, Mapping[str, int]],
                 rename: Mapping[str, str] | None = None) -> pd.DataFrame:
    """samples x ``<family>_indel`` (generate_indel_features.py). ``gpa`` has raw Roary family
    names as columns; ``candidates`` are (family, alignment stem); ``calls[stem]`` is that
    alignment's vcf2indel result."""
    cols = {}
    for fam, stem in candidates:
        if stem not in calls:
            continue
        c = calls[stem]
        if not any(v == 0 for v in c.values()):
            continue
        col = gpa[fam].astype(int).copy()
        for s, v in c.items():
            if s in col.index:
                col.at[s] = v
        cols[f"{rename.get(fam, fam) if rename else fam}_indel"] = col
    out = pd.DataFrame(cols, index=gpa.index)
    out.index.name = "sample"
    return out.astype(int) if len(cols) else out


def _run(cmd: list[str], cwd: str, stdin=None, stdout=None) -> None:
    subprocess.run(cmd, cwd=cwd, stdin=stdin, stdout=stdout, stderr=subprocess.PIPE, check=True)


def align_and_call(fasta: str | os.PathLike, workdir: str | os.PathLike,
                   mafft: str = "mafft", msa2vcf: str = "msa2vcf") -> tuple[Path, dict[str, int]]:
    """mafft (default options, as the published rule) then msa2vcf then vcf2indel."""
    fasta = Path(fasta)
    workdir = Path(workdir)
    aln = workdir / (fasta.stem + ".aln")
    vcf = workdir / (fasta.stem + ".vcf")
    with open(aln, "w") as out:
        _run([mafft, str(fasta.resolve())], cwd=str(workdir), stdout=out)
    with open(aln) as fi, open(vcf, "w") as fo:
        _run([msa2vcf], cwd=str(workdir), stdin=fi, stdout=fo)
    return aln, vcf2indel(vcf)


def build_indel_table(roary_dir: str | os.PathLike, ffn_files: Sequence[str | os.PathLike],
                      gpa: pd.DataFrame, out_dir: str | os.PathLike, threads: int = 1,
                      min_isolates: int = 2, legacy: bool = True,
                      mafft: str = "mafft", msa2vcf: str = "msa2vcf") -> tuple[pd.DataFrame, dict]:
    """The whole indel stage for one pangenome. Alignments stay in ``out_dir/aln`` with their
    index ``families.tsv``, because held-out genomes are later added to them (s2g2p.place)."""
    roary_dir, out_dir = Path(roary_dir), Path(out_dir)
    sel = select_families(roary_dir, min_isolates=min_isolates, legacy=legacy)
    aln_dir = out_dir / "aln"
    fastas = write_family_fastas(roary_dir / "clustered_proteins", sel.groups, sorted(map(str, ffn_files)), aln_dir)
    write_family_index(sel, aln_dir / "families.tsv")
    calls: dict[str, dict[str, int]] = {}

    def one(p: Path):
        return p.stem, align_and_call(p, aln_dir, mafft=mafft, msa2vcf=msa2vcf)[1]

    with ThreadPoolExecutor(max(1, threads)) as ex:
        for stem, c in ex.map(one, fastas):
            calls[stem] = c
    table = indel_matrix(gpa, sel.candidates, calls)
    stats = {"alignments": len(fastas), "families_aligned": len(sel.aligned),
             "indel_candidates": len(sel.candidates), "indel_columns": table.shape[1],
             "aligned_under_another_name": len(sel.renamed)}
    return table, stats


def add_and_call(train_aln: str | os.PathLike, records: Sequence[tuple[str, str]], sample: str,
                 mafft: str = "mafft", msa2vcf: str = "msa2vcf") -> int:
    """One held-out sample's call for one family: its copies are added to the training
    alignment (``mafft --add``), and only that sample is ever in the alignment besides the
    training genomes."""
    tmp = tempfile.mkdtemp(prefix="s2g2p_add_")
    try:
        with open(os.path.join(tmp, "add.fa"), "w") as fh:
            for h, s in records:
                fh.write(f">{h}\n{s}\n")
        with open(os.path.join(tmp, "comb.aln"), "w") as out:
            _run([mafft, "--add", "add.fa", "--thread", "1", str(Path(train_aln).resolve())], cwd=tmp, stdout=out)
        with open(os.path.join(tmp, "comb.aln")) as fi, open(os.path.join(tmp, "v.vcf"), "w") as fo:
            _run([msa2vcf], cwd=tmp, stdin=fi, stdout=fo)
        calls = vcf2indel(os.path.join(tmp, "v.vcf"))
        if sample not in calls:
            raise RuntimeError(f"{sample} absent from the msa2vcf output of {train_aln}")
        return calls[sample]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def read_ffn(path: str | os.PathLike) -> dict[str, tuple[str, str]]:
    """locus tag -> (full header, sequence)."""
    return {h.split(" ")[0]: (h, s) for h, s in read_fasta(path)}
