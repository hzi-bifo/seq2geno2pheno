# Versioning and releases

## The version number

`MAJOR.MINOR.PATCH` ([Semantic Versioning](https://semver.org)), written in PEP 440 form so the
same string works for the Python package and the workflow: `3.0.0`, pre-releases `3.0.0rc1`,
development `3.0.0.dev0`. Git tags carry a `v`: `v3.0.0`.

It is written in four places, and `tests/test_version.py` fails when they disagree:

| File | Field |
|---|---|
| `src/s2g2p/__init__.py` | `__version__` (the source; `pyproject.toml` reads it) |
| `nextflow.config` | `manifest.version` |
| `CITATION.cff` | `version` |
| `CHANGELOG.md` | the newest `## <version>` heading (without the pre-release suffix) |

Every evaluation records the version and the whole preset in `run.json`.

## What each part means here

Results are what users cite, so a version bump is decided by what happens to **results**, not
only to the interface.

- **MAJOR**: a sample sheet, fold file, preset file, command line or output that worked before
  now fails or means something else; or a saved model (`--mode train`) can no longer be loaded.
- **MINOR**: new presets, options or outputs; or the numbers an existing **non-reproduction**
  preset gives on the same input change, e.g. a placement fix in `noleak`. The CHANGELOG says
  which presets and why.
- **PATCH**: fixes that change no preset's numbers: crashes, messages, documentation,
  installation.

**Reproduction presets are frozen.** `hu2024` (status `reproduction`) must keep giving the
published numbers. If the published method is found to have been reproduced wrongly, the fix is
a MINOR release with the old behaviour available under the same preset only if it was itself
published; otherwise a new preset name. Its tool versions (`envs/*.lock.txt`) change only with a
MAJOR release.

## Releasing

1. All of `dev/run-tests.sbatch` passes; the release's validation is recorded in
   `docs/validation/`.
2. Set the version in the four places above (drop `.devN`), date the CHANGELOG entry.
3. Commit, tag `vX.Y.Z` (annotated), push the tag, make the GitHub release from it. Zenodo
   archives each release and mints its DOI.
4. Bump to the next `.dev0` on `main`.

v2 (`v2.00001`, 2021) is kept at its tag and on the `legacy-v2` branch.
