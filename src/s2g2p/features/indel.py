"""Indel features (Seq2Geno's ``indel_*`` rules), ported to Python 3.

Per gene family: align the members' nucleotide sequences (mafft), call variants against the
alignment (jvarkit ``msa2vcf``), and mark a sample 0 when one of its copies carries an indel of
at least 8 bp relative to the majority, 1 when it has the family without one, and 0 when it
lacks the family. A family becomes a column only if some sample carries an indel.

Two behaviours of the published rules are reproduced only under ``legacy=True`` (the
``hu2024`` and ``noleak`` presets), because they change which families get a column:

* The families aligned are the first N lines of Roary's ``clustered_proteins``, where N is the
  number of families in more than ``min_isolates`` isolates (``head -n $core_genes_wc``), not
  those N families themselves. A family gets an indel column only if it is in both sets.
* File names replace every non-word character of a family name with ``_``; the core list
  replaced only the first. A family whose name contains such a character therefore never found
  its indel file. Here it is skipped and counted, never silently.

The external tools are taken from ``PATH``; the workflow's containers pin their versions
(mafft 7.310, jvarkit-msa2vcf 1.0 in the published environment).
"""

from __future__ import annotations

import collections
import os
import re
import shutil
import subprocess
import tempfile
import warnings
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
    """The name a family's alignment file gets (gene_clusters2multi_fasta.py)."""
    return _NONWORD.sub("_", family)


@dataclass
class FamilySelection:
    aligned: list[str]                     # families whose alignment is computed
    candidates: list[str]                  # families that can become an indel column, in order
    skipped_nonword: list[str] = field(default_factory=list)


def select_families(roary_csv: str | os.PathLike, clustered_proteins: str | os.PathLike,
                    min_isolates: int = 2, legacy: bool = True) -> FamilySelection:
    """Which families are aligned and which can yield an indel column (``filter.awk`` with
    ``threshold=min_isolates``: strictly more isolates than the threshold)."""
    gpa = read_roary_gpa(roary_csv)
    n_iso = pd.to_numeric(gpa["No. isolates"], errors="raise")
    core = [str(f) for f in gpa.index[n_iso > min_isolates]]
    if not legacy:
        return FamilySelection(aligned=core, candidates=core)
    head = list(read_clustered_proteins(clustered_proteins))[:len(core)]
    aligned_keys = {file_key(f) for f in head}
    candidates, skipped = [], []
    for f in core:
        if file_key(f) != f:
            skipped.append(f)
        elif f in aligned_keys:
            candidates.append(f)
    if skipped:
        warnings.warn(f"{len(skipped)} families with non-word characters in their names get no "
                      f"indel column (legacy rule), e.g. {skipped[:3]}")
    return FamilySelection(aligned=head, candidates=candidates, skipped_nonword=skipped)


def write_family_fastas(clustered_proteins: str | os.PathLike, families: Sequence[str],
                        ffn_files: Iterable[str | os.PathLike], out_dir: str | os.PathLike) -> list[Path]:
    """One multi-FASTA of member nucleotide sequences per family, ``<file_key>.fasta``, records
    in the order the ``.ffn`` files and their records are read (gene_clusters2multi_fasta.py)."""
    fam2mem = read_clustered_proteins(clustered_proteins)
    wanted = set(families)
    prot2fam: dict[str, str] = {}
    for fam, members in fam2mem.items():
        if fam in wanted:
            for p in members:
                prot2fam.setdefault(p, fam)
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
                        if pid in prot2fam:
                            seqs[prot2fam[pid]].append(header + "".join(chunk))
                    header, chunk = line, []
                else:
                    chunk.append(line)
            if header is not None:
                pid = header.split(" ")[0][1:].rstrip("\n")
                if pid in prot2fam:
                    seqs[prot2fam[pid]].append(header + "".join(chunk))
    paths = []
    for fam in families:
        if fam in seqs:
            p = out_dir / f"{file_key(fam)}.fasta"
            p.write_text("".join(seqs[fam]))
            paths.append(p)
    return paths


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


def indel_matrix(gpa: pd.DataFrame, candidates: Sequence[str], calls: Mapping[str, Mapping[str, int]],
                 rename: Mapping[str, str] | None = None) -> pd.DataFrame:
    """samples x ``<family>_indel`` (generate_indel_features.py). ``gpa`` is the presence table
    with raw Roary family names as columns; ``calls[family]`` the vcf2indel result."""
    cols = {}
    for fam in candidates:
        if fam not in calls:
            continue
        c = calls[fam]
        if not any(v == 0 for s, v in c.items()):
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
    """The whole indel stage for one pangenome. Alignments are kept in ``out_dir/aln`` because
    held-out genomes are later added to them (see :mod:`s2g2p.place`)."""
    roary_dir, out_dir = Path(roary_dir), Path(out_dir)
    sel = select_families(roary_dir / "gene_presence_absence.csv", roary_dir / "clustered_proteins",
                          min_isolates=min_isolates, legacy=legacy)
    fasta_dir = out_dir / "aln"
    fastas = write_family_fastas(roary_dir / "clustered_proteins", sel.aligned, ffn_files, fasta_dir)
    key2fam = {file_key(f): f for f in sel.aligned}
    calls: dict[str, dict[str, int]] = {}

    def one(p: Path):
        return key2fam[p.stem], align_and_call(p, fasta_dir, mafft=mafft, msa2vcf=msa2vcf)[1]

    with ThreadPoolExecutor(max(1, threads)) as ex:
        for fam, c in ex.map(one, fastas):
            calls[fam] = c
    table = indel_matrix(gpa, sel.candidates, calls)
    stats = {"families_aligned": len(fastas), "indel_candidates": len(sel.candidates),
             "indel_columns": table.shape[1], "skipped_nonword": len(sel.skipped_nonword)}
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
