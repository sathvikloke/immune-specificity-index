#!/usr/bin/env python3
"""Compare the installed package versions against `requirements-lock.txt`.

    python scripts/20_check_versions.py          # warn on a mismatch, exit 0
    python scripts/20_check_versions.py --strict # exit 1 on a mismatch

WHY THIS EXISTS
===============
The repository pins numpy 2.1.3, pandas 2.2.3, scikit-learn 1.6.1 and scipy
1.15.3 in `requirements-lock.txt`, and every diagnostic script PRINTS the
versions it is running under. Until now nothing CHECKED them. Printing a version
is not verifying it: a reader running numpy 2.3 sees the number scroll past in
the environment banner and gets different results with no signal that anything
is wrong.

That gap matters more here than in most projects. A9 -- the open question about
why the NSCLC ISI is 0.318240180054566 on macOS/arm64 and 0.2964 on HPC4 --
exists precisely because a numerical result moved without the code, the seeds or
the inputs moving. An unnoticed version drift would look exactly like it, and
would waste a session's budget being investigated as a platform difference.

WHY IT WARNS RATHER THAN FAILS BY DEFAULT
=========================================
Aborting would be the wrong default. Someone reading the code, regenerating a
figure, or running the tests on a newer numpy should not be blocked; they should
be told, loudly, that the numbers they produce are not comparable to the frozen
ones. `--strict` is for the reproduction path, where a mismatch means the run
cannot claim to reproduce anything and should stop.

WHAT IT DOES NOT DO
===================
It checks the five packages whose arithmetic reaches the estimand (statsmodels
joined numpy, pandas, scikit-learn and scipy on 2026-09-17). The lockfile pins
the whole environment; a full audit of every transitive pin would fail on any
platform with a different wheel set and would therefore be ignored. The five
checked here are the ones a version change would silently move a number through.

WHY `pyarrow` IS NOT ON THE LIST -- decided 2026-09-05, deliberately
=====================================================================
`pyarrow` reads every parquet input in the project, so the question was raised
of whether it belongs in the enforced set. It does not, and the reason is the
criterion above rather than an oversight.

The criterion is "packages whose ARITHMETIC reaches the estimand". Parquet is a
lossless binary container: a float64 written to it is the same float64 when read
back, and `pyarrow` performs no arithmetic on the values on the way through. A
version change therefore cannot move a reported number the way a change to
numpy's reductions or scipy's `pearsonr` can -- which is precisely the class of
thing A9 is about.

Input fidelity IS guarded, just not here. It is guarded where it belongs:
`21_provenance_manifest.py` hashes 90 frozen artefacts, and `14_repro_probe.py`
stage 1 hashes the inputs themselves. A `pyarrow` bug that silently altered a
decoded value would change those hashes and fail both, which is a stronger
check than a version equality test -- it detects the harm rather than a proxy
for it.

`matplotlib` is excluded for the same reason and a weaker one still: it touches
no reported number at all, only rendering. Both are REPORTED in
`ENVIRONMENT.md`, so a reader diffing environments still sees them; they are
simply not gates.
"""

from __future__ import annotations

import argparse
import importlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "requirements-lock.txt"

# (distribution name in the lockfile, module to import). They differ for
# scikit-learn, which is the one people get wrong.
WATCHED = [
    ("numpy", "numpy"),
    ("pandas", "pandas"),
    ("scikit-learn", "sklearn"),
    ("scipy", "scipy"),
    # Added 2026-09-17 (session 50, ledger F6.10). statsmodels fits the OLS in
    # globalaxis, decomposition, outcome and experiment -- the index's own path
    # and the secondary tables -- and was pinned in the lockfile but never
    # compared against it. 0.14.4 on macOS and on HPC4's environment alike.
    ("statsmodels", "statsmodels"),
]

# The platform the frozen results were produced on, so a mismatch report says
# what the numbers are being compared against rather than just "expected X".
# BLAS: the library LOADED when the results were produced, 0.3.29 (conda); numpy's
# build record says 0.3.21, which is what this line said until 2026-09-16 (D5).
FROZEN_ON = ("macOS 15 / arm64, Python 3.13.9, BLAS openblas 0.3.29 at run time "
             "(numpy built against 0.3.21)")


def pinned() -> dict[str, str]:
    """Read the pins from the lockfile. Missing lockfile is fatal, not skipped."""
    if not LOCK.exists():
        raise SystemExit(
            f"LOCKFILE MISSING: {LOCK.relative_to(ROOT.parent)} -- refusing to "
            "check versions against remembered constants.")
    out: dict[str, str] = {}
    for line in LOCK.read_text().splitlines():
        m = re.match(r"^([A-Za-z0-9._-]+)==([^\s;#]+)", line.strip())
        if m:
            out[m.group(1).lower()] = m.group(2)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 on any mismatch (use on the reproduction path)")
    args = ap.parse_args()

    pins = pinned()
    print(f"Python  {sys.version.split()[0]}")
    try:
        import numpy as _np
        cfg = _np.show_config(mode="dicts")
        blas = cfg.get("Build Dependencies", {}).get("blas", {})
        if blas:
            print(f"BLAS    {blas.get('name', '?')} {blas.get('version', '?')} "
                  "(numpy's build record)")
        import scipy.linalg  # noqa: F401 -- load scipy's BLAS too before asking
        from threadpoolctl import threadpool_info
        for d in threadpool_info():
            if d.get("user_api") == "blas":
                print(f"BLAS    {d.get('internal_api')} {d.get('version')} loaded at run time "
                      f"({Path(str(d.get('filepath', ''))).name})")
    except Exception:  # noqa: BLE001 -- purely informational
        pass
    print(f"frozen results were produced on: {FROZEN_ON}\n")

    mismatches, unpinned, missing = [], [], []
    for dist, module in WATCHED:
        want = pins.get(dist)
        try:
            got = importlib.import_module(module).__version__
        except Exception as exc:  # noqa: BLE001
            missing.append(f"{dist}: cannot import {module} ({exc})")
            continue
        if want is None:
            unpinned.append(f"{dist}: installed {got}, NOT PINNED in the lockfile")
            print(f"  ?  {dist:14s} {got}  (no pin found)")
        elif got != want:
            mismatches.append(f"{dist}: installed {got}, lockfile pins {want}")
            print(f"  !! {dist:14s} {got}  != pinned {want}")
        else:
            print(f"  ok {dist:14s} {got}")

    problems = mismatches + unpinned + missing
    if not problems:
        print("\nOK: every watched package matches the lockfile.")
        return 0

    print("\n" + "=" * 72)
    for p in problems:
        print(f"  {p}")
    print(
        "\nNumbers produced under these versions are NOT comparable to the\n"
        "frozen results, and a difference must not be attributed to A9 (the\n"
        "open macOS-vs-HPC4 platform question) until the versions match.")
    print("=" * 72)
    if args.strict:
        print("\nFAILED (--strict).")
        return 1
    print("\nContinuing anyway (pass --strict to make this fatal).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
