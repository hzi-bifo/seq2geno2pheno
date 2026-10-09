import random

import pandas as pd
import pytest

from s2g2p.features.indel import file_key, indel_matrix, select_families, vcf2indel, write_family_fastas


def vcf2indel_reference(vcf):
    """Seq2Geno denovo/lib/vcf2indel.py (Python 2) with prints and file outputs removed and
    py2's integer `/` written as `//`; returns sample_has_indel."""
    with open(vcf) as vcf_open:
        for l in vcf_open:
            if not l.startswith("##"):
                break
        columns = l.strip("#").split("\t")
        sample_names = pd.Series(["_".join(i.split(" ")[0].split("_")[:-1]) for i in columns[9:]])
        sample_has_indel = dict([(i, 1) for i in sample_names])
        for l in vcf_open:
            columns = pd.Series(l.strip().split("\t"))
            _, pos, _, ref, alt, _, _, _, _ = columns[:9]
            samples = pd.Series([i.split("/")[0] for i in columns[9:]])
            alts = alt.split(",")
            for i, alt in zip(range(1, len(alts) + 1), alts):
                if not alt == "":
                    if abs(len(alt) - len(ref)) >= 8:
                        if ref == ".":
                            if len(sample_names.loc[samples == "."]) > sample_names.shape[0] // 2:
                                hits = sample_names.loc[samples != "."]
                            else:
                                hits = sample_names.loc[samples == "."]
                        elif alt == ".":
                            if len(sample_names.loc[samples == "."]) > sample_names.shape[0] // 2:
                                hits = sample_names.loc[samples != "."]
                            else:
                                hits = sample_names.loc[samples == "."]
                        elif len(sample_names.loc[samples == "."]) > sample_names.shape[0] // 2:
                            hits = sample_names.loc[(str(i) == samples) & (samples != ".")]
                        else:
                            hits = sample_names.loc[(str(i) == samples) | (samples == ".")]
                        for hit in hits:
                            sample_has_indel[hit] = 0
    return sample_has_indel


def random_vcf(path, rng, n_samples):
    names = [f"S{i}_{rng.randint(1, 99999):05d}" for i in range(n_samples)]
    lines = ["##fileformat=VCFv4.2\n",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n"]
    for p in range(rng.randint(1, 12)):
        ref = rng.choice([".", "A", "ACGTACGTACGTA", "AC"])
        alts = [rng.choice([".", "T", "TTTTTTTTTTTT", "ACGTACGTACGTACGTAC", "G"]) for _ in range(rng.randint(1, 3))]
        gts = [rng.choice(["."] + [str(i) for i in range(len(alts) + 1)]) for _ in names]
        lines.append(f"chr\t{p + 1}\t.\t{ref}\t{','.join(alts)}\t.\t.\t.\tGT\t" + "\t".join(gts) + "\n")
    path.write_text("".join(lines))


@pytest.mark.parametrize("seed", range(40))
def test_vcf2indel_matches_published(tmp_path, seed):
    rng = random.Random(seed)
    f = tmp_path / "x.vcf"
    random_vcf(f, rng, rng.randint(1, 9))
    assert vcf2indel(f) == vcf2indel_reference(f)


def roary_csv(path, fams):
    """fams: list of (name, n_isolates, {sample: locus})"""
    meta = ["Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
            "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment",
            "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc",
            "Avg group size nuc"]
    samples = sorted({s for _, _, m in fams for s in m})
    rows = [",".join(f'"{c}"' for c in meta + samples)]
    for name, n, m in fams:
        vals = [name, "", "", str(n), str(n), "1", "1", "1", "", "", "", "1", "1", "1"]
        vals += [m.get(s, "") for s in samples]
        rows.append(",".join(f'"{v}"' for v in vals))
    path.write_text("\n".join(rows) + "\n")


def test_select_families_legacy_rules(tmp_path):
    fams = [("famA", 3, {}), ("famB", 1, {}), ("fam-C", 4, {}), ("famD", 5, {}), ("famE", 3, {})]
    roary_csv(tmp_path / "gpa.csv", fams)
    # clustered_proteins in another order: Roary sorts it by its own rule
    (tmp_path / "cp").write_text("famD: a\nfamB: b\nfam-C: c\nfamA: d\nfamE: e\n")
    sel = select_families(tmp_path / "gpa.csv", tmp_path / "cp", min_isolates=2, legacy=True)
    # core (> 2 isolates): famA, fam-C, famD, famE -> 4 -> head -n 4 of clustered_proteins
    assert sel.aligned == ["famD", "famB", "fam-C", "famA"]
    assert sel.candidates == ["famA", "famD"]          # famE not aligned, fam-C has a non-word char
    assert sel.skipped_nonword == ["fam-C"]
    native = select_families(tmp_path / "gpa.csv", tmp_path / "cp", min_isolates=2, legacy=False)
    assert native.candidates == ["famA", "fam-C", "famD", "famE"]
    assert file_key("fam-C") == "fam_C"


def test_family_fastas(tmp_path):
    (tmp_path / "cp").write_text("famA: s1_00001\ts2_00001\nfamB: s1_00002\n")
    for s in ("s1", "s2"):
        (tmp_path / f"{s}.ffn").write_text(f">{s}_00001 prod\nACGT\nAC\n>{s}_00002 other\nGGG\n")
    paths = write_family_fastas(tmp_path / "cp", ["famA"], [tmp_path / "s1.ffn", tmp_path / "s2.ffn"], tmp_path / "out")
    assert [p.name for p in paths] == ["famA.fasta"]
    assert paths[0].read_text() == ">s1_00001 prod\nACGT\nAC\n>s2_00001 prod\nACGT\nAC\n"


def generate_indel_features_reference(gene_list, gpa, calls):
    """Seq2Geno denovo/lib/generate_indel_features.py without file I/O."""
    gpa = gpa.loc[:, gene_list].copy()
    indel_genes = []
    for l in gene_list:
        if l in calls:
            indels = pd.Series(calls[l])
            if (indels == 0).any():
                indel_genes.append(l)
                gpa.loc[indels.index, l] = indels
    out = gpa.loc[:, indel_genes]
    out.columns = ["%s_indel" % l for l in out.columns]
    return out


def test_indel_matrix_matches_published():
    gpa = pd.DataFrame({"f1": [1, 1, 0], "f2": [1, 1, 1], "f3": [1, 0, 1]}, index=["a", "b", "c"])
    calls = {"f1": {"a": 0, "b": 1}, "f2": {"a": 1, "b": 1, "c": 1}, "f3": {"a": 1, "c": 0}}
    ref = generate_indel_features_reference(["f1", "f2", "f3"], gpa, calls)
    out = indel_matrix(gpa, ["f1", "f2", "f3"], calls)
    assert list(out.columns) == ["f1_indel", "f3_indel"]
    assert (out.to_numpy() == ref.to_numpy()).all()
