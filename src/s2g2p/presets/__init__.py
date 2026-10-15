"""Named configurations. Each preset is one YAML file in this directory, read by both the
Python package (the ``ml`` section) and the Nextflow workflow (the ``features`` section), so a
preset cannot mean one thing to the features and another to the model."""

from __future__ import annotations

from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
TRANSFORMS = {"standard_global", "standard", "l1", "dedup", "none"}
FEATURE_SETS = {"concat", "separate"}
TIE_BREAKS = {"first", "lowest_std"}
BLOCKS = ("kmer", "gpa", "indel")


def available() -> list[str]:
    return sorted(p.stem for p in HERE.glob("*.yml"))


def path_of(name_or_path: str) -> Path:
    p = Path(name_or_path)
    if p.suffix in (".yml", ".yaml") and p.is_file():
        return p
    q = HERE / f"{name_or_path}.yml"
    if not q.is_file():
        raise KeyError(f"unknown preset {name_or_path!r}; available: {', '.join(available())}")
    return q


def load(name_or_path: str) -> dict:
    with open(path_of(name_or_path)) as fh:
        cfg = yaml.safe_load(fh)
    validate(cfg)
    return cfg


def validate(cfg: dict) -> None:
    ml = cfg.get("ml") or {}
    for key in ("blocks", "transforms", "feature_sets", "classifiers", "grid", "inner_cv", "tie_break"):
        if key not in ml:
            raise ValueError(f"preset {cfg.get('name')}: ml.{key} missing")
    bad = [b for b in ml["blocks"] if b not in BLOCKS]
    if bad:
        raise ValueError(f"preset {cfg.get('name')}: unknown blocks {bad}")
    for b in ml["blocks"]:
        t = ml["transforms"].get(b)
        if t not in TRANSFORMS:
            raise ValueError(f"preset {cfg.get('name')}: transform of {b} is {t!r}; one of {sorted(TRANSFORMS)}")
    if ml["feature_sets"] not in FEATURE_SETS:
        raise ValueError(f"preset {cfg.get('name')}: feature_sets must be one of {sorted(FEATURE_SETS)}")
    if ml["tie_break"] not in TIE_BREAKS:
        raise ValueError(f"preset {cfg.get('name')}: tie_break must be one of {sorted(TIE_BREAKS)}")
    kind = str(ml["inner_cv"]).split(":")[0]
    if kind not in ("other_folds", "group_kfold", "stratified"):
        raise ValueError(f"preset {cfg.get('name')}: inner_cv {ml['inner_cv']!r} not understood")
    if cfg.get("predict") and any(ml["transforms"][b] == "standard_global" for b in ml["blocks"]):
        raise ValueError(f"preset {cfg.get('name')}: predict needs no cohort-fitted transform")


def describe(name: str) -> str:
    cfg = load(name)
    return f"{cfg['name']:<8} {cfg.get('status', ''):<14} {' '.join(str(cfg.get('summary', '')).split())}"


__all__ = ["available", "load", "validate", "path_of", "describe", "HERE", "BLOCKS"]
