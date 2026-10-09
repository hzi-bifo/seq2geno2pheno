"""Placing genomes into a pangenome they did not help build.

The leak-free presets build the pangenome, the family alignments and the indel calls from
training genomes only; every other genome (an outer test fold, or a genome to predict) is then
placed into what the training genomes defined, one genome at a time:

* Gene presence: each protein joins the family of its best training hit under Roary's own edge
  criterion -- blastp ``-evalue 1E-6 -max_target_seqs 2000``, kept if identity > 95% and the
  alignment covers > 95% of both sequences (Roary 3.8.2, ``Bio/Roary/External/Blastp.pm``). A
  protein identical to a training protein takes that protein's family without a search.
* Indels: the genome's copies of each indel family are added to the training alignment and
  called as in :func:`s2g2p.features.indel.add_and_call`.
"""

from __future__ import annotations

import collections
import gzip
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd

from .features.gpa import protein_to_family, read_clustered_proteins
from .features.indel import add_and_call, file_key, read_ffn
from .io import read_fasta

BLAST_EVALUE = "1E-6"
BLAST_MAX_TARGETS = "2000"
MIN_PIDENT = 95.0
MIN_COVER = 0.95
BLAST_FMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"


def training_sequence_families(clustered_proteins: str | os.PathLike,
                               train_faa: Mapping[str, str | os.PathLike]) -> tuple[dict[str, str], dict]:
    """Distinct training protein sequence -> family. A sequence found in several families takes
    the one most of its copies are in (ties: the smaller name)."""
    prot2fam = protein_to_family(read_clustered_proteins(clustered_proteins))
    votes: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    found = total = 0
    for _, faa in train_faa.items():
        for name, seq in read_fasta(faa):
            total += 1
            pid = name.split(" ")[0]
            if pid in prot2fam:
                votes[seq][prot2fam[pid]] += 1
                found += 1
    if prot2fam and found < 0.95 * len(prot2fam):
        raise ValueError(f"only {found} of {len(prot2fam)} clustered proteins were found in the "
                         "training .faa files; are these the genomes Roary ran on?")
    seq2fam = {s: sorted(c.items(), key=lambda x: (-x[1], x[0]))[0][0] for s, c in votes.items()}
    return seq2fam, {"train_proteins": total, "train_clustered": found}


def write_sequence_families(seq2fam: Mapping[str, str], path: str | os.PathLike) -> None:
    """The training side of placement, saved with a model so predicting needs no training
    annotations: one ``family<TAB>sequence`` row per distinct training protein."""
    with gzip.open(path, "wt") as fh:
        fh.write("family\tsequence\n")
        for seq, fam in seq2fam.items():
            fh.write(f"{fam}\t{seq}\n")


def read_sequence_families(path: str | os.PathLike) -> dict[str, str]:
    with gzip.open(path, "rt") as fh:
        next(fh)
        return {seq: fam for fam, seq in (l.rstrip("\n").split("\t") for l in fh)}


def best_hits(hits_tsv: str | os.PathLike) -> dict[str, int]:
    """query id -> index of its best subject under Roary's criterion (by bitscore, then
    identity, then the earlier subject). Ids are ``q<i>`` / ``s<i>``."""
    best: dict[str, tuple[tuple, int]] = {}
    with open(hits_tsv) as fh:
        for line in fh:
            a = line.rstrip("\n").split("\t")
            pid, length, slen, qlen = float(a[2]), int(a[3]), int(a[12]), int(a[13])
            if not (pid > MIN_PIDENT and length / slen > MIN_COVER and length / qlen > MIN_COVER):
                continue
            si = int(a[1][1:])
            key = (float(a[11]), pid, -si)
            if a[0] not in best or key > best[a[0]][0]:
                best[a[0]] = (key, si)
    return {q: si for q, (_, si) in best.items()}


def place_gpa(seq2fam: Mapping[str, str], query_faa: Mapping[str, str | os.PathLike],
              columns: Sequence[str], threads: int = 1, blastp: str = "blastp",
              makeblastdb: str = "makeblastdb", workdir: str | os.PathLike | None = None,
              rename: Mapping[str, str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """(presence table over ``columns``, per-protein assignment, stats) for the query genomes.

    ``columns`` are the training table's columns, so the placed rows line up with it. A family
    a query protein joins that is not a column is an error: it means the training table and the
    pangenome disagree.
    """
    q_seqs: "collections.OrderedDict[str, list[tuple[str, str]]]" = collections.OrderedDict()
    for sample, faa in query_faa.items():
        for name, seq in read_fasta(faa):
            q_seqs.setdefault(seq, []).append((sample, name.split(" ")[0]))
    assign = {q: seq2fam[q] for q in q_seqs if q in seq2fam}
    rest = [q for q in q_seqs if q not in seq2fam]
    if rest:
        db = list(seq2fam)
        tmp = tempfile.mkdtemp(prefix="s2g2p_place_", dir=workdir)
        with open(os.path.join(tmp, "db.faa"), "w") as fh:
            for i, s in enumerate(db):
                fh.write(f">s{i}\n{s}\n")
        with open(os.path.join(tmp, "q.faa"), "w") as fh:
            for i, s in enumerate(rest):
                fh.write(f">q{i}\n{s}\n")
        subprocess.run([makeblastdb, "-in", "db.faa", "-dbtype", "prot", "-out", "db"], cwd=tmp,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        with open(os.path.join(tmp, "hits.tsv"), "w") as out:
            subprocess.run([blastp, "-query", "q.faa", "-db", "db", "-evalue", BLAST_EVALUE,
                            "-num_threads", str(threads), "-outfmt", BLAST_FMT,
                            "-max_target_seqs", BLAST_MAX_TARGETS], cwd=tmp, check=True, stdout=out,
                           stderr=subprocess.PIPE)
        for qid, si in best_hits(os.path.join(tmp, "hits.tsv")).items():
            assign[rest[int(qid[1:])]] = seq2fam[db[si]]
    colset = set(columns)
    present: dict[str, set] = collections.defaultdict(set)
    rows, n_prot, n_hit = [], collections.Counter(), collections.Counter()
    unknown = set()
    for q, owners in q_seqs.items():
        fam = assign.get(q)
        for sample, locus in owners:
            n_prot[sample] += 1
            if fam is None:
                continue
            n_hit[sample] += 1
            name = rename.get(fam, fam) if rename else fam
            (present[sample] if name in colset else unknown).add(name)
            rows.append((sample, locus, fam))
    if unknown:
        raise ValueError(f"{len(unknown)} assigned families are not training columns, e.g. {sorted(unknown)[:3]}")
    samples = list(query_faa)
    table = pd.DataFrame([[1 if c in present[s] else 0 for c in columns] for s in samples],
                         index=samples, columns=list(columns), dtype=int)
    table.index.name = "sample"
    assignment = pd.DataFrame(rows, columns=["sample", "locus", "family"])
    frac = [n_hit[s] / n_prot[s] for s in samples if n_prot[s]]
    stats = {"query_genomes": len(samples), "query_proteins": sum(n_prot.values()),
             "query_distinct_proteins": len(q_seqs), "query_exact": len(q_seqs) - len(rest),
             "query_blasted": len(rest),
             "assigned_frac_min": round(min(frac), 4) if frac else "",
             "assigned_frac_mean": round(sum(frac) / len(frac), 4) if frac else ""}
    return table, assignment, stats


def place_indels(indel_columns: Sequence[str], aln_dir: str | os.PathLike, assignment: pd.DataFrame,
                 query_ffn: Mapping[str, str | os.PathLike], threads: int = 1,
                 mafft: str = "mafft", msa2vcf: str = "msa2vcf") -> tuple[pd.DataFrame, dict]:
    """Indel calls of the query genomes for the training table's indel columns
    (``<family>_indel``). A genome without the family is 0, as the training table encodes it;
    genomes with identical copies of a family are aligned once."""
    by_fam: dict[str, dict[str, list[str]]] = collections.defaultdict(lambda: collections.defaultdict(list))
    for s, locus, fam in assignment.itertuples(index=False):
        by_fam[fam][s].append(locus)
    ffn = {s: read_ffn(p) for s, p in query_ffn.items()}
    jobs, groups = [], {}
    for col in indel_columns:
        fam = col[: -len("_indel")]
        aln = Path(aln_dir) / f"{file_key(fam)}.aln"
        if not aln.exists():
            raise FileNotFoundError(f"no training alignment {aln} for indel column {col}")
        by_content: "collections.OrderedDict[tuple, list]" = collections.OrderedDict()
        for s in query_ffn:
            recs = [ffn[s][l] for l in by_fam.get(fam, {}).get(s, []) if l in ffn[s]]
            if recs:
                by_content.setdefault(tuple(sorted(q for _, q in recs)), []).append((s, recs))
        for members in by_content.values():
            rep, recs = members[0]
            groups[(col, rep)] = [m[0] for m in members]
            jobs.append((col, str(aln), rep, recs))

    def one(job):
        col, aln, rep, recs = job
        return col, rep, add_and_call(aln, recs, rep, mafft=mafft, msa2vcf=msa2vcf)

    calls = {}
    with ThreadPoolExecutor(max(1, threads)) as ex:
        for col, rep, v in ex.map(one, jobs):
            for s in groups[(col, rep)]:
                calls[(col, s)] = v
    samples = list(query_ffn)
    table = pd.DataFrame([[calls.get((c, s), 0) for c in indel_columns] for s in samples],
                         index=samples, columns=list(indel_columns), dtype=int)
    table.index.name = "sample"
    return table, {"indel_columns": len(indel_columns), "indel_jobs": len(jobs)}
