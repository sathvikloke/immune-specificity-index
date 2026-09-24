#!/usr/bin/env python3
"""Write ENVIRONMENT.md from the RUNNING interpreter, never from memory.

    python scripts/23_export_environment.py            # print, do not write
    python scripts/23_export_environment.py --write    # write ENVIRONMENT.md

WHY THIS EXISTS
===============
A9. The NSCLC ISI is 0.318240180054566 on macOS/arm64 and 0.2964 on
Linux/x86_64 with the same seeds, the same pinned versions and inputs verified
byte-identical. A reader who reproduces a third number has no way to tell
whether they have hit A9 or their own environment drift. This file gives them
something to diff against.

WHY IT READS EVERYTHING LIVE
============================
Handoff #20 asserted that the macOS/HPC4 contrast was Accelerate-versus-
OpenBLAS. Measured, this Mac is **openblas** (numpy's build record says 0.3.21;
the library loaded at run time is 0.3.29). That claim was plausible,
repeated across sessions, and wrong, and it survived only because nobody
printed the value. Every row below is therefore read from the interpreter at
run time. Nothing here is a literal, so nothing here can go stale silently.

NOT A SUBSTITUTE FOR THE VERSION CHECK
======================================
This script REPORTS. `20_check_versions.py` DECIDES: it compares the four
packages whose arithmetic reaches the estimand against requirements-lock.txt
and exits non-zero on drift. A report nobody diffs is not a check.
"""

from __future__ import annotations

import argparse
import importlib
import os
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Reported for context. The FIVE that 20_check_versions.py actually enforces
# are numpy, pandas, scikit-learn, scipy and statsmodels -- the ones whose
# arithmetic reaches the estimand (statsmodels added 2026-09-17, F6.10).
# matplotlib and pyarrow render and load; they cannot move the number.
PACKAGES = ("numpy", "pandas", "scipy", "sklearn", "statsmodels", "matplotlib",
            "pyarrow")

THREAD_VARS = ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
               "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONHASHSEED")


def blas_info() -> dict:
    """numpy's own record of what it was linked against.

    Wrapped because the `mode="dicts"` form is newer than some numpy builds and
    the shape of the dict is not guaranteed across versions. A missing key must
    print as `<unknown>` rather than raise -- an environment report that
    crashes on an unexpected environment is useless precisely when it matters.
    """
    try:
        import numpy
        return numpy.__config__.show(mode="dicts")["Build Dependencies"]["blas"]
    except Exception as exc:  # noqa: BLE001 - reporting tool, must not die
        return {"name": f"<unavailable: {type(exc).__name__}>"}


def scipy_build_blas() -> dict:
    try:
        import scipy
        return scipy.show_config(mode="dicts")["Build Dependencies"]["blas"]
    except Exception as exc:  # noqa: BLE001 - reporting tool, must not die
        return {"name": f"<unavailable: {type(exc).__name__}>"}


def blas_runtime() -> list[dict]:
    """The BLAS libraries actually LOADED in this process (ledger D5, 2026-09-16).

    numpy's record above is what numpy was BUILT against, which is not
    necessarily what runs. Measured 2026-09-16: on the machine that produced the
    frozen results numpy records OpenBLAS 0.3.21 while the only OpenBLAS in the
    process is conda's 0.3.29, shared by numpy and scipy; on HPC4 the numpy and
    scipy wheels each load their own (0.3.27 and 0.3.28, job 129294). File names
    only -- the directories are local paths.
    """
    try:
        import numpy  # noqa: F401 - load both BLAS consumers before asking
        import scipy.linalg  # noqa: F401
        from threadpoolctl import threadpool_info
        return [{"api": d.get("internal_api"), "version": d.get("version"),
                 "file": Path(str(d.get("filepath", ""))).name,
                 "architecture": d.get("architecture")}
                for d in threadpool_info() if d.get("user_api") == "blas"]
    except Exception as exc:  # noqa: BLE001
        return [{"api": f"<unavailable: {type(exc).__name__}>", "version": "", "file": "",
                 "architecture": ""}]


def git_commit() -> str:
    try:
        r = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                           capture_output=True, text=True, check=True)
        return r.stdout.strip()
    except Exception:  # noqa: BLE001
        return "<not a git checkout>"


def render() -> str:
    b = blas_info()
    sb = scipy_build_blas()
    loaded = "; ".join(f"{r['api']} {r['version']} (`{r['file']}`, {r['architecture']})"
                       for r in blas_runtime()) or "<none reported>"
    pkgs = {}
    for m in PACKAGES:
        try:
            pkgs[m] = importlib.import_module(m).__version__
        except Exception as exc:  # noqa: BLE001
            pkgs[m] = f"<not importable: {type(exc).__name__}>"

    rows_pkg = "\n".join(f"| {k} | {v} |" for k, v in pkgs.items())
    rows_thr = "\n".join(f"| `{v}` | {os.environ.get(v, '(unset)')} |"
                         for v in THREAD_VARS)

    return f"""# Exact environment, read from the running interpreter

**Generated by `scripts/23_export_environment.py`. Every value below was read
from the running interpreter, not written from memory.** Regenerate with
`python scripts/23_export_environment.py --write`.

This file exists because of A9: the NSCLC ISI is 0.318240180054566 on this
platform and 0.2964 on Linux/x86_64 with the same pinned versions and
byte-identical inputs. A reader who reproduces a different number needs to be
able to tell whether they are seeing A9 or their own setup drift. Comparing
against this file is how.

**The table below describes the machine this file was last generated on, which
is not necessarily the machine that produced the frozen results.** The two were
the same until 2026-09-20, when the laptop was upgraded from **Darwin 25.5.0**
-- the value this file carried from its first version on 2026-09-04 until that
day -- to Darwin 27.0.0. The frozen results were re-verified on the new OS
immediately after the upgrade and nothing moved: provenance 220/220, the
stable-sort partition still hashing to `9f31e33c2caaf1e3`, and the number
checker green. The title said "the environment that produced the frozen
results" until 2026-09-20; it was renamed when that stopped being true, rather
than regenerated into a quiet contradiction.

**Which macOS produced the frozen results, answered 2026-09-21.** Until
2026-09-20 nothing in this project recorded the macOS *marketing* version --
only `platform.release()`, the Darwin string -- so the `product version` row
below is new, and the manuscript's "macOS 15 / arm64" had never been checked.
It was wrong. The laptop's own install record (`softwareupdate --history`,
cross-checked against `/Library/Receipts/InstallHistory.plist`) shows macOS
26.5.2 installed on 2026-07-19 and no further OS install until macOS 27.0 on
2026-09-20, and the frozen runs were launched on this laptop on 2026-08-19/20.
The frozen results were therefore produced on **macOS 26.5.2 (Darwin 25.5.0) /
arm64**, which agrees with the Darwin string this file carried from 2026-09-04.
The two version numbers do NOT coincide on that release -- Darwin 25.5.0 is
macOS 26.5.2 -- so neither can be inferred from the other, and this file now
records both.

## Platform

| | |
|---|---|
| system | {platform.system()} {platform.release()} |
| product version | {platform.mac_ver()[0] or "(not macOS)"} |
| machine | {platform.machine()} |
| processor | {platform.processor() or "(not reported)"} |
| python | {sys.version.split()[0]} ({platform.python_implementation()}) |
| build | {" ".join(platform.python_build())} |
| compiler | {platform.python_compiler()} |
| git commit | `{git_commit()}` |

## BLAS -- the part A9 turns on

| | |
|---|---|
| name | **{b.get("name", "<unknown>")}** |
| version numpy was built against | {b.get("version", "<unknown>")} |
| version scipy was built against | {sb.get("version", "<unknown>")} |
| **loaded at run time** | **{loaded}** |
| detection | {b.get("detection method", "<unknown>")} |
| openblas config | {b.get("openblas configuration", "<unknown>")} |
| pc file | {b.get("pc file directory", "<unknown>")} |

**Handoff #20 asserted this contrast was Accelerate-vs-OpenBLAS. It is not.**
This machine is OpenBLAS, confirmed three independent ways. That wrong claim
survived several sessions because nobody printed it; this file prints it.

## Packages

| package | version |
|---|---|
{rows_pkg}

`scripts/20_check_versions.py` checks the five whose arithmetic reaches the
estimand (numpy, pandas, scikit-learn, scipy, statsmodels) against
`requirements-lock.txt` and exits non-zero on drift. The others are listed for
context only and are **NOT** checked -- if another ever reaches the estimand,
add it to that script's watch list rather than relying on this file.

## Threading

**No BLAS thread-count variable is set for macOS runs, deliberately.** Thread
count changes floating-point reduction order, and A9 *is* a numerical
reproducibility investigation, so the environment is kept identical to the one
that produced the frozen results. The cluster sbatch scripts DO set them;
cluster numbers are diagnostic, never manuscript numbers.

Observed at generation time:

| variable | value |
|---|---|
{rows_thr}

`PYTHONHASHSEED=0` is exported by `reproduce.sh`, because str hashing is
randomised per process and gene-pool ordering must be deterministic. If it
reads `(unset)` above, this file was generated outside `reproduce.sh` -- that
does not affect any other value, but it is why the row is printed rather than
assumed.

## Seed audit

Every source of randomness in the pipeline, and where its seed is set.
**Re-derived 2026-09-05 by reading every call site, not by grep.** The previous
version of this table was written from a grep and listed four sources; there are
**ten** `np.random.default_rng` call sites in `src/aacr27/`. [Corrected
2026-09-22: this sentence said "twelve". By the method it names there were ten
at the commit that wrote it (`cd69375`), ten at its parent, eight at the
pre-registration snapshot, and ten now; the only count that reaches twelve is of
the generators' USES (`rng.permutation(`, `rng.integers(` and the like), which is
a different quantity. The same re-read found one call site this table never
listed -- the random-patient split's, now its own row -- and one row whose
function is never called.] All ten are
deterministic, but six of them were undocumented and four take their seed from a
DEFAULT ARGUMENT rather than from `AuditConfig`, which is a weaker guarantee than
this table used to imply: they agree with the configured seed because
`seed: int = 0` happens to match `AuditConfig.seed = 0`, not because anything
wires them together. Changing `AuditConfig.seed` would silently NOT change them.

| source | call site | seeded by | notes |
|---|---|---|---|
| null gene-set draws | `signatures.random_gene_sets` | `AuditConfig.seed` (0) | per-signature nulls deliberately SHARE the seed, so two signatures of identical size and expression profile draw identical sets. Intentional for the family-level statistic; it means the 16 nulls are not independent (Limitations 7). |
| *(none -- never called)* | `stats.bootstrap_corr` | default `seed=0` | **DEAD CODE.** Defined, and seeded by its own default, but nothing in `src/`, `scripts/` or `tests/` calls it -- measured by grep 2026-09-22, and equally true at `cd69375`, the commit that wrote this table. No draw in the delivered analysis comes from it. Until 2026-09-22 this row was titled "patient-clustered bootstrap" and described "one patient set resampled per draw and shared across signatures"; that is the pooled-ISI bootstrap's behaviour (next row), which a function taking one x/y pair cannot have. |
| pooled-ISI bootstrap | `experiment._pooled_isi_bootstrap` | `config.seed`, explicit | the patient-clustered interval: one patient set resampled per draw and SHARED across signatures, which preserves between-signature correlation. Named rather than numbered: this row read `experiment.py:944` until 2026-09-08, by which point 944 was a blank line and the function had moved to 946. |
| Δr split-robustness bootstrap | `stats.paired_delta` | `config.seed`, explicit in `experiment.run_audit`'s `paired_delta` call | Previously undocumented. |
| site-to-fold partition | `splits._assign_sites_greedy` | `AuditConfig.fold_seed` | a derived property, not a field: `AuditConfig.fold_seed` returns `self.seed if self.split_seed is None else self.split_seed`, so the partition follows `seed` **unless** `split_seed` is set. Varied ONLY by `12_partition_variance.py`, which holds `seed` fixed so the spread it reports is attributable to the partition and nothing else. The exact path (`preserved_site_split(exact=True)`, a cvxpy QP) consumes no randomness, and its two greedy fallbacks pass NO seed -- but nothing calls `exact=True` (measured 2026-09-22), so the delivered partition always takes the seeded greedy path. **The seed is not the whole story here — see A9 below.** |
| random-patient split | `splits.random_patient_split` | `AuditConfig.fold_seed`, explicit in `experiment.run_audit` | the site-agnostic comparator. `run_audit` passes no `stratify`, so the path that runs is the `default_rng` shuffle of patient order; with `stratify` it would instead pass `random_state=` to `StratifiedGroupKFold`. **Missing from this table until 2026-09-22:** its stratified branch was listed in the row above as the site-to-fold partition's "stratified path", which it is not, and the branch that actually runs was not listed at all. Its fold assignment was platform-dependent for a reason no seed controls -- A12 in `14-SCIENCE-AUDIT.md` -- until 2026-09-24, when the group order was pinned (`splits.stable_group_kfold`); the frozen runs used the library's. |
| split-half reliability | `signatures.split_half_reliability` | explicit in `13_scorer_sensitivity.py`'s `selfcheck_split_half`; default `seed=0` | the scorer-agnostic reliability used by A7's reconstruction. Previously undocumented. |
| Control C permutation calibration | `decomposition.permutation_calibration` | **DEFAULT `seed=0` only** — `11_close_science_gaps.py`'s `permutation_calibration` call passes no seed | deterministic, but not wired to `AuditConfig`. Previously undocumented. This row cited `11_close_science_gaps.py:177` until 2026-09-08, which is a docstring about the fingerprint guard; the call is in that script's `main()`. |
| ancestry Cramér's V calibration | `ancestry.cramers_v_calibrated` | **DEFAULT `seed=0` only** — neither call in `06_run_ancestry.py` passes a seed | the identifiability result in Supplementary Table S8. Previously undocumented. |
| ancestry pooling | `ancestry.pooled_over_signatures` | **DEFAULT `seed=0` only** | Previously undocumented. |
| rotation null | *(no RNG)* | — | `stats.rotation_null` **consumes no randomness at all**: it reuses the already-drawn `null_r` matrix and computes a permutation p from it. The previous table listed it as "`AuditConfig.seed` (0), 1,000 draws", which was true of the draws it reads but wrong about the function, which is a pure transformation of its input. |
| str hashing / gene-pool order | — | `PYTHONHASHSEED=0` | exported by `reproduce.sh`. NOT set by default in an interactive shell -- see the table above. |

Every call site above is named by FUNCTION, not by line. Until 2026-09-22 all but
one were `file.py:NNN` coordinates, and by then two had drifted onto blank or
unrelated lines (`signatures.py:353`, `experiment.py:448`) while the one row
named on 2026-09-08 was still right -- a pointer is correct only until someone
edits above it.

**What is NOT seeded, and does not need to be:** the ridge solve and
`scipy.stats.pearsonr` are deterministic given their inputs. That is precisely
why A9 is interesting -- with every seed pinned and the inputs hashed
identical, the answer still moves across platforms.

<!-- END GENERATED -- everything below this line is hand-authored.

     `test_environment_md_matches_a_fresh_render_on_this_platform` compares the
     committed file to a fresh `render()` only DOWN TO this marker. The sections
     below record things no machine can introspect: the second platform, and
     which non-stable sorts were load-bearing and why. Before this marker
     existed, adding either meant choosing between an accurate document and a
     green test, and the test would have lost. Keep new hand-authored material
     below the marker; keep anything derivable from the running interpreter
     above it, in `render()`, where it cannot go stale. -->
"""


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true",
                    help="write ENVIRONMENT.md instead of printing to stdout")
    args = ap.parse_args()
    text = render()
    if args.write:
        out = ROOT / "ENVIRONMENT.md"
        text = _splice(text, out)
        out.write_text(text)
        print(f"written to {out}")
    else:
        sys.stdout.write(text)
    return 0


MARKER = "<!-- END GENERATED"


def _splice(fresh: str, out) -> str:
    """Keep everything the existing file has BELOW the marker.

    MEASURED BUG, 2026-09-05: `--write` used to do `out.write_text(render())`,
    which silently deleted every hand-authored section below the marker. It
    destroyed the entire "Sort stability" section -- the A9/A10 record, the
    six-row table of every non-stable ordering and the guidance for anyone
    extending the pipeline -- and it was recovered only because the file was in
    git. The marker itself survived, because `render()` emits it, so a check
    that looked only for the marker reported success. That is precisely the
    shape of failure this project keeps meeting: the guard was present and the
    thing it guarded was gone.

    The contract the marker promises is "everything below this line is
    hand-authored", and a writer that honours it must therefore READ the file
    before overwriting it. If the file does not exist yet, or has no marker,
    the fresh render is written whole -- there is nothing to preserve.
    """
    if not out.exists():
        return fresh
    existing = out.read_text()
    if MARKER not in existing or MARKER not in fresh:
        return fresh
    return fresh[:fresh.index(MARKER)] + existing[existing.index(MARKER):]


if __name__ == "__main__":
    sys.exit(main())
