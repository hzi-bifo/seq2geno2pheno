"""One version, everywhere it is written (docs/RELEASING.md)."""

import re
from pathlib import Path

import s2g2p

ROOT = Path(__file__).resolve().parents[1]
# PEP 440 release or pre-release: 3.0.0, 3.0.0.dev0, 3.0.0rc1, 3.1.0a2 ...
VERSION = re.compile(r"^\d+\.\d+\.\d+((a|b|rc)\d+|\.dev\d+)?$")


def test_package_version_is_pep440_semver():
    assert VERSION.match(s2g2p.__version__), s2g2p.__version__


def test_workflow_manifest_matches():
    text = (ROOT / "nextflow.config").read_text()
    m = re.search(r"^\s*version\s*=\s*'([^']+)'", text, re.M)
    assert m and m.group(1) == s2g2p.__version__


def test_citation_matches():
    text = (ROOT / "CITATION.cff").read_text()
    m = re.search(r"^version:\s*(\S+)", text, re.M)
    assert m and m.group(1) == s2g2p.__version__


def test_changelog_top_entry_matches():
    text = (ROOT / "CHANGELOG.md").read_text()
    m = re.search(r"^## \[?(\d[^\] ]*)\]?", text, re.M)
    release = re.sub(r"(\.dev\d+|(a|b|rc)\d+)$", "", s2g2p.__version__)
    assert m and m.group(1) == release, "CHANGELOG.md's newest entry must be this version"
