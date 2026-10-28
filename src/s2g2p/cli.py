"""``s2g2p`` command line.

Feature subcommands (``kmers``, ``gpa``, ``indels``, ``compress``, ``place``) are what the
Nextflow workflow calls inside its containers; they can also be run by hand on existing Prokka
and Roary output. Model subcommands (``folds``, ``evaluate``, ``train``, ``predict``) need only
this package.
"""

from __future__ import annotations

import argparse
import datetime
import json
import platform
import sys
from pathlib import Path

import pandas as pd

from . import __version__, io, presets


def _log(msg: str) -> None:
    print(f"[s2g2p {datetime.datetime.now():%H:%M:%S}] {msg}", file=sys.stderr, flush=True)


def _samples_arg(args, sheet: pd.DataFrame | None = None) -> list[str] | None:
    if getattr(args, "samples_file", None):
        return [l.strip() for l in open(args.samples_file) if l.strip()]
    if sheet is not None:
        return list(sheet.index)
    return None


def _annotation(annot_dir: str, sample: str, ext: str) -> str:
    p = Path(annot_dir) / sample / f"{sample}.{ext}"
    if not p.is_file():
        raise FileNotFoundError(f"no {ext} for {sample}: {p}")
    return str(p)


# ----------------------------------------------------------------------------------- features
def cmd_kmers(args):
    from .features.kmer import assemblies_from_sheet, kmer_matrix

    sheet = io.read_samplesheet(args.samplesheet)
    asm = assemblies_from_sheet(sheet)
    keep = _samples_arg(args)
    if keep is not None:
        asm = {s: asm[s] for s in keep}
    df = kmer_matrix(asm, k=args.k, n_jobs=args.threads)
    io.write_matrix(df, args.out)
    _log(f"{df.shape[0]} samples x {df.shape[1]} canonical {args.k}-mers -> {args.out}")


def cmd_gpa(args):
    from .features.gpa import gpa_from_roary, read_rename_map

    rename = read_rename_map(args.rename) if args.rename else None
    df = gpa_from_roary(args.roary_csv, _samples_arg(args), rename)
    io.write_matrix(df, args.out)
    _log(f"{df.shape[0]} samples x {df.shape[1]} families -> {args.out}")


def cmd_compress(args):
    from .features.compress import compress, write_groups

    df = io.read_matrix(args.input)
    out, groups = compress(df)
    io.write_matrix(out, args.out)
    if args.groups:
        write_groups(groups, args.groups)
    _log(f"{df.shape[1]} -> {out.shape[1]} columns -> {args.out}")


def cmd_indels(args):
    from .features.indel import build_indel_table

    gpa = io.read_matrix(args.gpa, dtype=int)
    samples = list(gpa.index)
    ffn = [_annotation(args.annotations, s, "ffn") for s in samples]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    table, stats = build_indel_table(args.roary_dir, ffn, gpa, out_dir, threads=args.threads,
                                     min_isolates=args.min_isolates, legacy=not args.no_legacy_rules,
                                     mafft=args.mafft, msa2vcf=args.msa2vcf)
    io.write_matrix(table, out_dir / "indel.tsv")
    io.write_manifest(out_dir, [dict(file="indel.tsv", samples=table.shape[0], columns=table.shape[1], **stats)],
                      ["file", "samples", "columns"] + list(stats))
    _log(f"indels: {stats}")


def cmd_place(args):
    from .features.gpa import read_rename_map
    from .place import place_gpa, place_indels, read_sequence_families, training_sequence_families

    pg = Path(args.pangenome)
    train_gpa = io.read_matrix(pg / "gpa.tsv", dtype=int)
    train_indel = io.read_matrix(pg / "indel.tsv", dtype=int) if (pg / "indel.tsv").is_file() else None
    query = [l.strip() for l in open(args.query) if l.strip()]
    overlap = set(query) & set(train_gpa.index)
    if overlap:
        raise SystemExit(f"{len(overlap)} query genomes are in the training pangenome, e.g. {sorted(overlap)[:3]}")
    query_faa = {s: _annotation(args.annotations, s, "faa") for s in query}
    rename = read_rename_map(args.rename) if args.rename else None
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    saved = pg / "train_sequences.tsv.gz"
    if saved.is_file():
        seq2fam, st1 = read_sequence_families(saved), {"train_sequences": str(saved)}
    else:
        train_faa = {s: _annotation(args.annotations, s, "faa") for s in train_gpa.index}
        seq2fam, st1 = training_sequence_families(pg / "roary" / "clustered_proteins", train_faa)
    q_gpa, assignment, st2 = place_gpa(seq2fam, query_faa, list(train_gpa.columns), threads=args.threads,
                                       blastp=args.blastp, makeblastdb=args.makeblastdb, workdir=out_dir,
                                       rename=rename)
    assignment.to_csv(out_dir / "assignment.tsv.gz", sep="\t", index=False)
    io.write_matrix(pd.concat([train_gpa, q_gpa]), out_dir / "gpa.tsv")
    rows = [dict(file="gpa.tsv", train=len(train_gpa), placed=len(q_gpa), columns=train_gpa.shape[1], **st1, **st2)]
    if train_indel is not None:
        q_ind, st3 = place_indels(list(train_indel.columns), pg / "aln", assignment,
                                  {s: _annotation(args.annotations, s, "ffn") for s in query}, q_gpa,
                                  threads=args.threads, mafft=args.mafft, msa2vcf=args.msa2vcf)
        io.write_matrix(pd.concat([train_indel, q_ind]), out_dir / "indel.tsv")
        rows.append(dict(file="indel.tsv", train=len(train_indel), placed=len(q_ind),
                         columns=train_indel.shape[1], **st3))
    cols = list(dict.fromkeys(k for r in rows for k in r))
    io.write_manifest(out_dir, rows, cols)
    _log(f"placed {len(query)} genomes: {st2}")


def cmd_sequences(args):
    from .place import training_sequence_families, write_sequence_families

    pg = Path(args.pangenome)
    train = list(io.read_matrix(pg / "gpa.tsv").index)
    seq2fam, st = training_sequence_families(pg / "roary" / "clustered_proteins",
                                             {s: _annotation(args.annotations, s, "faa") for s in train})
    write_sequence_families(seq2fam, args.out)
    _log(f"{len(seq2fam)} distinct training proteins ({st}) -> {args.out}")


# -------------------------------------------------------------------------------------- model
def _load_blocks(args, fold: int | None, blocks: list[str]) -> dict[str, pd.DataFrame]:
    out = {}
    for b in blocks:
        tmpl = getattr(args, b)
        if not tmpl:
            raise SystemExit(f"the preset uses the {b} block; pass --{b}")
        path = io.expand_fold_path(tmpl, fold) if fold is not None else tmpl
        out[b] = io.read_matrix(path, dtype=float)
    return out


def _phenotypes(sheet: pd.DataFrame, wanted: str) -> list[str]:
    cols = io.phenotype_columns(sheet)
    return cols if wanted == "all" else [wanted]


def _folds_for(args, labels: pd.Series, phenotype: str) -> list[list[str]]:
    from .folds import restrict, stratified_folds, validate

    if args.folds:
        folds = restrict(io.read_folds(args.folds.replace("{phenotype}", phenotype)), list(labels.index))
    else:
        folds = stratified_folds(labels, args.n_folds, args.seed)
    validate(folds, list(labels.index))
    return folds


def _groups(sheet: pd.DataFrame) -> pd.Series | None:
    return sheet["group"] if "group" in sheet.columns and (sheet["group"] != "").all() else None


def cmd_evaluate(args):
    from .ml.evaluate import evaluate

    cfg = presets.load(args.preset)
    sheet = io.read_samplesheet(args.samplesheet)
    cache: dict = {}

    def blocks_for(fold):
        key = tuple(io.expand_fold_path(getattr(args, b) or "", fold) for b in cfg["ml"]["blocks"])
        if key not in cache:
            cache.clear()
            cache[key] = _load_blocks(args, fold, cfg["ml"]["blocks"])
        return cache[key]

    for ph in _phenotypes(sheet, args.phenotype):
        labels = io.labels_for(sheet, ph)
        folds = _folds_for(args, labels, ph)
        out = Path(args.out_dir) / ph
        out.mkdir(parents=True, exist_ok=True)
        _log(f"{ph}: {len(labels)} labelled genomes, {len(folds)} folds, preset {cfg['name']}")
        t0 = datetime.datetime.now()
        res = evaluate(blocks_for, labels, folds, cfg["ml"], groups=_groups(sheet), n_jobs=args.threads, log=_log)
        files = {"candidates.tsv": res.candidates, "folds.tsv": res.folds,
                 "predictions.tsv": res.predictions, "summary.tsv": res.summary}
        for name, df in files.items():
            df.to_csv(out / name, sep="\t", index=False)
        io.write_folds(folds, out / "folds_used.tsv")
        run = {"s2g2p_version": __version__, "python": platform.python_version(),
               "sklearn_version": __import__("sklearn").__version__, "preset": cfg,
               "phenotype": ph, "inputs": {b: getattr(args, b) for b in cfg["ml"]["blocks"]},
               "folds": args.folds or f"stratified n={args.n_folds} seed={args.seed}",
               "dropped_grid_keys": res.dropped_grid_keys,
               "started": t0.isoformat(timespec="seconds"),
               "finished": datetime.datetime.now().isoformat(timespec="seconds")}
        (out / "run.json").write_text(json.dumps(run, indent=2, default=str))
        f1 = res.summary.set_index("metric").loc["f1_macro"]
        io.write_manifest(out, [{"file": n, "rows": len(d)} for n, d in files.items()] +
                          [{"file": "folds_used.tsv", "rows": sum(map(len, folds))}, {"file": "run.json", "rows": 1}],
                          ["file", "rows"])
        _log(f"{ph}: f1_macro {f1['mean']:.4f} +- {f1['std']:.4f} (support-weighted {f1['weighted_mean']:.4f})")


def cmd_train(args):
    from .ml import model

    cfg = presets.load(args.preset)
    sheet = io.read_samplesheet(args.samplesheet)
    labels = io.labels_for(sheet, args.phenotype)
    folds = _folds_for(args, labels, args.phenotype)
    bundle = model.train(_load_blocks(args, None, cfg["ml"]["blocks"]), labels, folds, cfg,
                         groups=_groups(sheet), n_jobs=args.threads, log=_log)
    model.save(bundle, args.out)
    _log(f"model: {bundle['feature_set']} {bundle['classifier']} cv f1_macro {bundle['cv_f1_macro']:.4f} -> {args.out}")


def cmd_predict(args):
    from .ml import model

    bundle = model.load(args.model)
    blocks = _load_blocks(args, None, bundle["blocks"])
    samples = _samples_arg(args)
    out = model.predict(bundle, blocks, samples)
    out.to_csv(args.out, sep="\t", index=False)
    _log(f"{len(out)} predictions -> {args.out}")


def cmd_folds(args):
    """Outer folds per phenotype by one of s2g2p.folds.METHODS; --out-dir writes
    <phenotype>_<method>_cv.json per phenotype and manifest.tsv."""
    from . import folds as F

    sheet = io.read_samplesheet(args.samplesheet)
    phenotypes = _phenotypes(sheet, args.phenotype)
    if bool(args.out) == bool(args.out_dir):
        raise SystemExit("pass one of --out-dir, or --out with a single --phenotype")
    if args.out and len(phenotypes) != 1:
        raise SystemExit("--out takes a single --phenotype; use --out-dir for several")
    if args.method == "phylogeny" and not args.tree:
        raise SystemExit("--method phylogeny needs --tree (Newick)")
    if args.method == "kma" and not args.kma_clusters:
        raise SystemExit("--method kma needs --kma-clusters (kma_clustering's log; '{phenotype}' is replaced)")
    leaves = F.newick_leaves(Path(args.tree).read_text()) if args.method == "phylogeny" else None
    rows = []
    for ph in phenotypes:
        labels = io.labels_for(sheet, ph)
        samples = list(labels.index)
        info: dict = {}
        if args.method == "random":
            folds = F.stratified_folds(labels, args.n_folds, args.seed)
            info = {"source": f"stratified seed={args.seed}"}
        elif args.method == "phylogeny":
            folds = F.tree_folds(leaves, samples, args.n_folds)
            info = {"source": args.tree}
        else:
            path = args.kma_clusters.replace("{phenotype}", ph)
            clusters = F.read_kma_clusters(path)
            folds = F.kma_folds(clusters, samples, args.n_folds, args.kma_seed)
            info = {"source": path, "clusters": len(clusters),
                    "unclustered": len(F.unclustered(clusters, samples))}
        F.validate(folds, samples)
        if args.out:
            target = Path(args.out)
        else:
            Path(args.out_dir).mkdir(parents=True, exist_ok=True)
            target = Path(args.out_dir) / f"{ph}_{args.method}_cv.json"
        io.write_folds(folds, target)
        rows.append({"file": target.name, "phenotype": ph, "method": args.method, "genomes": len(samples),
                     "folds": len(folds), "fold_sizes": ",".join(str(len(f)) for f in folds),
                     "resistant_per_fold": ",".join(str(int(labels[f].sum())) for f in folds), **info})
        _log(f"{ph}: {args.method} folds over {len(samples)} genomes, sizes {rows[-1]['fold_sizes']} -> {target}")
    if args.out_dir:
        cols = ["file", "phenotype", "method", "genomes", "folds", "fold_sizes", "resistant_per_fold", "source"]
        cols += [c for c in ("clusters", "unclustered") if any(c in r for r in rows)]
        io.write_manifest(args.out_dir, rows, cols)


def cmd_kma_input(args):
    from .folds import write_kma_input

    sheet = io.read_samplesheet(args.samplesheet)
    labels = io.labels_for(sheet, args.phenotype)
    asm = (io.read_samplesheet(args.assemblies) if args.assemblies else sheet)["assembly"]
    write_kma_input([(s, asm[s]) for s in labels.index], args.out)
    _log(f"{args.phenotype}: {len(labels)} genomes, contigs joined by Ns -> {args.out}")


def cmd_split(args):
    """Per phenotype and outer fold: the held-out genomes, and the pangenome genomes (the whole
    cohort minus them -- genomes unlabelled for this phenotype stay in, as NOLEAK.md decided)."""
    sheet = io.read_samplesheet(args.samplesheet)
    rows = []
    for ph in _phenotypes(sheet, args.phenotype):
        labels = io.labels_for(sheet, ph)
        folds = _folds_for(args, labels, ph)
        d = Path(args.out_dir) / ph
        d.mkdir(parents=True, exist_ok=True)
        io.write_folds(folds, d / "folds_used.tsv")
        for k, test in enumerate(folds):
            held = set(test)
            train = [s for s in sheet.index if s not in held]
            (d / f"fold_{k}.test.txt").write_text("".join(s + "\n" for s in test))
            (d / f"fold_{k}.train.txt").write_text("".join(s + "\n" for s in train))
            rows.append({"phenotype": ph, "fold": k, "test": len(test), "train": len(train),
                         "files": f"{ph}/fold_{k}.train.txt,{ph}/fold_{k}.test.txt"})
    io.write_manifest(args.out_dir, rows, ["phenotype", "fold", "test", "train", "files"])
    _log(f"{len(rows)} (phenotype, fold) splits -> {args.out_dir}")


def cmd_preset_get(args):
    v = presets.load(args.name)
    for k in args.key.split("."):
        v = v[k]
    print(json.dumps(v) if isinstance(v, (dict, list)) else v)


def cmd_presets(args):
    if args.name:
        print(presets.path_of(args.name).read_text())
    else:
        for n in presets.available():
            print(presets.describe(n))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="s2g2p", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"s2g2p {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("presets", help="list presets, or print one")
    p.add_argument("name", nargs="?")
    p.set_defaults(fn=cmd_presets)

    p = sub.add_parser("preset-get", help="print one value of a preset (for the workflow)")
    p.add_argument("name")
    p.add_argument("key", help="dotted key, e.g. features.roary_args")
    p.set_defaults(fn=cmd_preset_get)

    p = sub.add_parser("kmers", help="canonical k-mer counts from the sample sheet's assemblies")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--samples-file", help="only these samples (one per line)")
    p.add_argument("--k", type=int, default=6)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_kmers)

    p = sub.add_parser("gpa", help="gene presence/absence from Roary's gene_presence_absence.csv")
    p.add_argument("--roary-csv", required=True)
    p.add_argument("--samples-file", help="rows to emit (default: every isolate column)")
    p.add_argument("--rename", help="Seq2Geno roary_abricate.txt (family -> column name)")
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_gpa)

    p = sub.add_parser("compress", help="featCompress: one column per distinct pattern")
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--groups", help="also write the group membership file")
    p.set_defaults(fn=cmd_compress)

    tools = argparse.ArgumentParser(add_help=False)
    tools.add_argument("--mafft", default="mafft")
    tools.add_argument("--msa2vcf", default="msa2vcf")
    tools.add_argument("--threads", type=int, default=1)

    p = sub.add_parser("indels", parents=[tools], help="indel table from a Roary run")
    p.add_argument("--roary-dir", required=True)
    p.add_argument("--annotations", required=True, help="directory of <sample>/<sample>.ffn (Prokka layout)")
    p.add_argument("--gpa", required=True, help="the presence table of the same run (gpa.tsv)")
    p.add_argument("--min-isolates", type=int, default=2)
    p.add_argument("--no-legacy-rules", action="store_true",
                   help="align every family in > min-isolates isolates (not Seq2Geno's head -n subset)")
    p.add_argument("--out-dir", required=True)
    p.set_defaults(fn=cmd_indels)

    p = sub.add_parser("place", parents=[tools], help="place genomes into a training pangenome")
    p.add_argument("--pangenome", required=True, help="dir with roary/, gpa.tsv, indel.tsv, aln/")
    p.add_argument("--annotations", required=True, help="directory of <sample>/<sample>.{faa,ffn}")
    p.add_argument("--query", required=True, help="file of genome ids to place, one per line")
    p.add_argument("--rename", help="Seq2Geno roary_abricate.txt, if the training table used one")
    p.add_argument("--blastp", default="blastp")
    p.add_argument("--makeblastdb", default="makeblastdb")
    p.add_argument("--out-dir", required=True)
    p.set_defaults(fn=cmd_place)

    p = sub.add_parser("sequences", help="save a pangenome's training protein -> family table")
    p.add_argument("--pangenome", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_sequences)

    model_in = argparse.ArgumentParser(add_help=False)
    model_in.add_argument("--kmer", help="k-mer table (TSV)")
    model_in.add_argument("--gpa", help="presence table; '{fold}' is replaced by the outer fold")
    model_in.add_argument("--indel", help="indel table; '{fold}' is replaced by the outer fold")
    model_in.add_argument("--threads", type=int, default=1)

    folds_in = argparse.ArgumentParser(add_help=False)
    folds_in.add_argument("--folds", help="fold file (JSON list of lists, or TSV sample/fold); "
                                          "'{phenotype}' is replaced. Default: stratified random")
    folds_in.add_argument("--n-folds", type=int, default=10)
    folds_in.add_argument("--seed", type=int, default=0)

    p = sub.add_parser("evaluate", parents=[model_in, folds_in], help="nested cross-validation")
    p.add_argument("--preset", default="noleak")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--phenotype", default="all")
    p.add_argument("--out-dir", required=True)
    p.set_defaults(fn=cmd_evaluate)

    p = sub.add_parser("train", parents=[model_in, folds_in], help="fit a final model on all labelled genomes")
    p.add_argument("--preset", default="noleak")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--phenotype", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_train)

    p = sub.add_parser("predict", parents=[model_in], help="apply a trained model")
    p.add_argument("--model", required=True)
    p.add_argument("--samples-file", help="only these samples")
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_predict)

    p = sub.add_parser("split", parents=[folds_in], help="per-fold train/test genome lists for the workflow")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--phenotype", default="all")
    p.add_argument("--out-dir", required=True)
    p.set_defaults(fn=cmd_split)

    p = sub.add_parser("folds", help="outer folds: random, phylogeny-aware or homology-aware (KMA)")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--method", choices=["random", "phylogeny", "kma"], default="random")
    p.add_argument("--phenotype", default="all")
    p.add_argument("--n-folds", type=int, default=10)
    p.add_argument("--seed", type=int, default=0, help="random: the StratifiedKFold seed")
    p.add_argument("--tree", help="phylogeny: Newick tree whose leaves are the sample ids")
    p.add_argument("--kma-clusters", help="kma: kma_clustering's log; '{phenotype}' is replaced")
    p.add_argument("--kma-seed", type=int, default=42, help="kma: the cluster allocation seed (published: 42)")
    p.add_argument("--out-dir", help="<phenotype>_<method>_cv.json per phenotype, and manifest.tsv")
    p.add_argument("--out", help="one fold file (.json, else TSV); needs a single --phenotype")
    p.set_defaults(fn=cmd_folds)

    p = sub.add_parser("kma-input", help="one record per genome, contigs joined by Ns (kma_clustering input)")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--phenotype", required=True, help="genomes labelled for it, in sample-sheet order")
    p.add_argument("--assemblies", help="sample/assembly table overriding the sheet's paths (staged copies)")
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_kma_input)

    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
