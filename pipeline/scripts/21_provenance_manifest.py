#!/usr/bin/env python3
"""SHA-256 every FROZEN artefact, so an edit to one cannot pass unnoticed.

    python scripts/21_provenance_manifest.py            # verify against the manifest
    python scripts/21_provenance_manifest.py --write    # (re)generate it

THE GAP THIS CLOSES
===================
`pipeline/results/README.md` classifies each directory as FROZEN, SUPERSEDED,
DIAGNOSTIC or DEAD, and the project's first guardrail is "do not change any
frozen result". Nothing enforced it. `19_check_numbers.py` checks that the
DOCUMENTS agree with the artefacts -- it reads the artefacts as ground truth, so
if a frozen CSV were edited, every document would be re-verified against the
edited value and the checker would go green. That is the one failure mode the
guardrails most want to prevent and the one nothing detected.

This is the complement: it checks that the ARTEFACTS THEMSELVES have not moved.
Run together, `19` says the documents match the artefacts and `21` says the
artefacts are the ones that were frozen.

WHAT IS COVERED, AND WHY NOT EVERYTHING
=======================================
Only files that back a quoted number: the frozen and superseded run directories'
`summary.json` / `*.csv` / `*.npz`, plus the A9 probe `.npz` files that
`19_check_numbers.py` now reads as an authority. Deliberately excluded:

  * `.log` files -- they are appended to and their timestamps vary; the run
    directories they describe are what carry the numbers.
  * `results/figures/` -- regenerating figures is an expected operation
    (workstream 3), and a changed PNG is not a changed result.
  * `_cache_global_axis_*` -- caches, documented as safe to delete.
  * `results/README.md` -- prose, and it is meant to be edited.

A file that is in the manifest and has vanished is reported as MISSING, which is
a failure. A file that exists and is not in the manifest is reported as NEW,
which is NOT a failure -- new runs are expected -- but it is printed, so a new
artefact cannot quietly accumulate without someone deciding whether it is frozen.

RUNNING THIS INSIDE THE PUBLIC DEPOSIT
======================================
Found 2026-09-07 by running this script inside a staged snapshot, which nobody
had done: it printed **FAILED: 0 changed, 4 missing** to any reader who tried to
verify the deposit. The four are `bisect_{macos,hpc4}.npz` and
`alpha_margin{,_fast}_macos.npz` -- `.npz` files the snapshot excludes on
purpose, because they are megabytes and answer an internal platform question.
Nothing was wrong with the deposit. The integrity checker was reporting
corruption that did not exist, which is worse than checking nothing: a check
that prints a claim it cannot verify is a liability.

The manifest therefore records, per artefact, whether the public snapshot stages
it. That flag is computed at `--write` time by IMPORTING `22_build_public_snapshot.py`
and asking its own rule sets -- never by restating them here. A second definition
of "what is staged" would be free to disagree with the first, which is the defect
this whole class keeps producing.

At verify time an absent artefact is tolerated ONLY when it is flagged
`deposited: false` AND this tree is a public deposit. The deposit is detected by
the absence of `scripts/22_build_public_snapshot.py`, which excludes itself by
design (its docstring necessarily quotes every string it forbids). On the working
tree that file is present, so every absence is still a failure there. A manifest
written before this flag existed has no `deposited` key; that is treated as
"unknown" and fails on absence, i.e. the old behaviour.

HOW TO USE IT WHEN A CHANGE IS LEGITIMATE
=========================================
If a frozen artefact is *supposed* to change -- it almost never is -- rerun with
`--write` and commit the manifest diff **in its own commit, with the reason**.
The point is not that the bytes can never change; it is that changing them must
be a deliberate, reviewable act rather than a side effect.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
MANIFEST = RESULTS / "PROVENANCE.json"


def _show(path: Path) -> str:
    """Path for display. Falls back to the absolute form when the manifest
    is not under ROOT -- which happens whenever a caller points RESULTS at a
    temporary tree, as the test that proves this script can FAIL does.
    `relative_to` raised ValueError there and crashed the check itself."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)

# Run directories whose contents back a quoted number. `results/README.md` is
# the authority on the classification; this list mirrors it.
TRACKED_DIRS = [
    "nsclc_v3", "pancancer_v3",            # FROZEN -- the headlines
    "nsclc_v2", "pancancer_v2", "pancancer", "nsclc_real", "nsclc_fixed",
    "partition_variance_nsclc",            # FROZEN -- A5
    "scorer_sensitivity_local",            # FROZEN -- A7
    "scorer_sensitivity_hpc4",             # DIAGNOSTIC, but quoted as provenance
    "ancestry",                            # FROZEN -- the supplement
    # Added 2026-09-06. All three are SENSITIVITY runs, not frozen primaries --
    # but each is now read as an authority by 19_check_numbers.py, and the note
    # on bisect_hpc4.npz above applies unchanged: an authority that is not
    # hashed is an authority that can drift, and this script exists precisely
    # because the number checker cannot detect an edited artefact.
    "nsclc_v3_stablesort",                 # A9/A10 -- Limitations 8's second number
    "scorer_sensitivity_ssgsea_stable",    # A10 -- the pinned ssGSEA arm
    "partition_variance_pancancer",        # the pan-cancer partition-variance run
]

# Standalone files read as an authority by 19_check_numbers.py, or quoted.
TRACKED_FILES = [
    "bisect_macos.npz",
    # Added 2026-09-06 with its macOS counterpart's reasoning: 19_check_numbers.py
    # now reads `pooled_z_excess` out of this file as the authority behind three
    # full-precision quotes in 14-SCIENCE-AUDIT.md, and the checker cannot detect
    # an edited artefact by construction -- which is the whole reason this script
    # exists. An authority that is not hashed is an authority that can drift.
    "bisect_hpc4.npz",
    "alpha_margin_macos.npz",
    "alpha_margin_fast_macos.npz",
    # Added 2026-09-05 for the same reason as bisect_hpc4.npz above:
    # 19_check_numbers.py reads five keys out of each of these as the authority
    # behind the four full-precision synthetic-ISI values and the tied-site and
    # moved-patient counts quoted in 09-PAPER-DRAFT.md and 14-SCIENCE-AUDIT.md.
    # One file per machine, each written by a real run on that machine.
    "split_determinism_darwin_arm64.json",
    "split_determinism_linux_x86_64.json",
    # Added 2026-09-05, same reasoning again: `19_check_numbers.py` reads eleven
    # keys out of each of these as the authority behind the A10 numbers quoted in
    # 09-PAPER-DRAFT.md and 14-SCIENCE-AUDIT.md, AND asserts a relationship
    # BETWEEN the two files -- stable depth hashes must agree across platforms,
    # quicksort hashes must not. That assertion is only as trustworthy as the
    # files, so both are hashed. One file per machine, each written by a real run
    # on that machine (HPC4's by cluster job 121782).
    "sort_audit_darwin_arm64.json",
    "sort_audit_linux_x86_64.json",
    # Added 2026-09-07. This one is hashed for a DIFFERENT reason than the rest:
    # no document quotes it. It is the 16 post-filtering panel sizes figure 2
    # plots against, deposited so that the public snapshot -- which excludes
    # `pipeline/data/` on purpose -- can render figure 2 at all. Before it
    # existed, `10_make_figures.py` exited 1 there for every reader reproducing
    # from the deposit. It is hashed because it is the only figure input that is
    # a COPY of a quantity living somewhere else: `10_make_figures.py` compares
    # the two whenever data/ is present, but a reader of the snapshot has no
    # data/ and therefore no way to notice drift. Hashing is their check.
    "gene_set_sizes.csv",
]

TRACKED_SUFFIXES = {".json", ".csv", ".npz", ".npy", ".tsv"}


BUILDER = ROOT / "scripts" / "22_build_public_snapshot.py"


def _is_public_deposit() -> bool:
    """True when this tree looks like the published snapshot rather than the
    working repository. The marker is the ABSENCE of the snapshot builder, which
    excludes itself by design and is therefore never present in a deposit."""
    return not BUILDER.is_file()


def _snapshot_rules() -> tuple[set[str], set[str]]:
    """The builder's OWN (suffix, name) keep-sets, loaded from source once.

    Deliberately not reimplemented here: two definitions of "what is staged"
    would be free to disagree, and every defect in this class so far has lived
    in exactly that gap."""
    spec = importlib.util.spec_from_file_location("_snapshot_builder", BUILDER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.RESULT_KEEP_SUFFIX, mod.RESULT_KEEP_NAMES


def _deposited(rel: str, rules: tuple[set[str], set[str]]) -> bool:
    """Does the public snapshot stage `results/<rel>`?"""
    keep_suffix, keep_names = rules
    p = Path(rel)
    return p.suffix in keep_suffix or p.name in keep_names


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collect() -> dict[str, dict[str, object]]:
    """Path (relative to results/) -> {sha256, bytes}, sorted for a stable diff."""
    out: dict[str, dict[str, object]] = {}
    for name in TRACKED_DIRS:
        d = RESULTS / name
        if not d.is_dir():
            continue
        for p in sorted(d.rglob("*")):
            if p.is_file() and p.suffix in TRACKED_SUFFIXES:
                out[str(p.relative_to(RESULTS))] = {
                    "sha256": _sha256(p), "bytes": p.stat().st_size}
    for name in TRACKED_FILES:
        p = RESULTS / name
        if p.is_file():
            out[str(p.relative_to(RESULTS))] = {
                "sha256": _sha256(p), "bytes": p.stat().st_size}
    return dict(sorted(out.items()))


def verify() -> int:
    if not MANIFEST.exists():
        print(f"NO MANIFEST at {_show(MANIFEST)}.")
        print("Generate it with:  python scripts/21_provenance_manifest.py --write")
        return 1

    recorded = json.loads(MANIFEST.read_text())["files"]
    current = collect()

    deposit = _is_public_deposit()
    changed, missing, new, undeposited = [], [], [], []
    for rel, meta in recorded.items():
        cur = current.get(rel)
        if cur is None:
            # Absent. Tolerated ONLY inside the public deposit, and only for an
            # artefact the manifest says the snapshot does not stage. A manifest
            # written before that flag existed reports None -> still a failure.
            if deposit and meta.get("deposited") is False:
                undeposited.append(rel)
            else:
                missing.append(rel)
        elif cur["sha256"] != meta["sha256"]:
            changed.append((rel, meta, cur))
    for rel in current:
        if rel not in recorded:
            new.append(rel)

    print(f"Provenance: {len(recorded)} recorded artefact(s), "
          f"{len(current)} found on disk.")

    for rel in new:
        print(f"  NEW      {rel}  ({current[rel]['bytes']} b) -- not in the "
              "manifest. Not a failure; decide whether it is frozen, then "
              "--write.")

    for rel in undeposited:
        print(f"  NOT DEPOSITED  {rel}  -- excluded from the public snapshot on "
              "purpose,\n                 so its absence here is expected, not "
              "corruption.")

    if not changed and not missing:
        if undeposited:
            print(f"\nOK: every DEPOSITED artefact still hashes to its frozen "
                  f"value\n({len(recorded) - len(undeposited)} of "
                  f"{len(recorded)} verified; {len(undeposited)} not staged "
                  "into this snapshot).")
        else:
            print("\nOK: every recorded artefact still hashes to its frozen value.")
        return 0

    for rel in missing:
        print(f"  MISSING  {rel}  -- recorded, absent on disk.")
    for rel, meta, cur in changed:
        print(f"  CHANGED  {rel}")
        print(f"    recorded {meta['sha256'][:16]}  ({meta['bytes']} b)")
        print(f"    on disk  {cur['sha256'][:16]}  ({cur['bytes']} b)")

    print(f"\nFAILED: {len(changed)} changed, {len(missing)} missing.")
    print("A frozen artefact changed. If that was deliberate, rerun with "
          "--write and\ncommit the manifest diff on its own, with the reason.")
    return 1


def write() -> int:
    files = collect()
    # The flag is written here, on the working tree, because the builder it is
    # derived from is excluded from the snapshot and cannot be consulted there.
    rules = _snapshot_rules()
    for rel, meta in files.items():
        meta["deposited"] = _deposited(rel, rules)
    MANIFEST.write_text(json.dumps(
        {"note": "SHA-256 of every frozen/quoted artefact. Regenerate ONLY "
                 "when a change is deliberate; see scripts/21_provenance_manifest.py.",
         "files": files}, indent=1) + "\n")
    total = sum(int(m["bytes"]) for m in files.values())
    n_dep = sum(1 for m in files.values() if m["deposited"])
    print(f"wrote {_show(MANIFEST)}: {len(files)} artefact(s), "
          f"{total:,} bytes total")
    print(f"  {n_dep} staged into the public snapshot, "
          f"{len(files) - n_dep} excluded from it on purpose")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true",
                    help="(re)generate the manifest instead of verifying it")
    args = ap.parse_args()
    return write() if args.write else verify()


if __name__ == "__main__":
    raise SystemExit(main())
