#!/usr/bin/env python3
"""NSCLC sensitivity of the index to four analysis choices (ledger E1, E3, E4, E5).

    python scripts/29_nsclc_sensitivity.py                      # every variant
    python scripts/29_nsclc_sensitivity.py --variant global_axis

Each variant changes ONE thing in the pre-registered NSCLC configuration and
recomputes the index; nothing else moves (seed 0, the same 1,000 null sets
unless the variant is about them, the same bootstrap draws). Every variant runs
the ISI-only configuration of `12_partition_variance.py --isi-only` (no
covariate baselines, no decomposition, no ComBat, no legacy null, no outcome
arm), which leaves the index untouched; the `reference` variant proves that on
this machine by returning the corrected-ordering primary
(`results/nsclc_v3_stablesort/`) bit for bit, and it is what every other
variant is compared against.

  reference        nothing changed
  global_axis      E1: the global expression axis is computed over both
                   cancer types together and residualised globally, instead of
                   within type (AuditConfig.axis_within_type=False)
  unmatched_null   E3: random gene sets matched on size only, not on mean
                   expression. A DECOMPOSITION of the registered estimand --
                   what expression matching buys -- never a robustness check.
  null_boot_400    E4: the bootstrap recomputes the null mean from 400 of the
  null_boot_1000       1,000 draws per signature (registered: 200), then all
  k3, k10          E5: 3 or 10 cross-validation folds instead of 5

Writes results/nsclc_sensitivity/<variant>/ (the audit's usual files) and
results/nsclc_sensitivity/summary.json. Needs data/interim.
"""

from __future__ import annotations

import argparse
import functools
import json
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import experiment, globalaxis  # noqa: E402

OUT = ROOT / "results" / "nsclc_sensitivity"
REFERENCE = ROOT / "results" / "nsclc_v3_stablesort" / "summary.json"
ISI_ONLY = {"baselines": (), "purity_col": "__isi_only_no_such_column__",
            "run_combat": False, "run_legacy_venet_null": False}
VARIANTS = {
    "reference": {},
    "global_axis": {"axis_within_type": False},
    "unmatched_null": {"match_null_expression": False},
    "null_boot_400": {"_max_null_for_boot": 400},
    "null_boot_1000": {"_max_null_for_boot": 1000},
    "k3": {"n_folds": 3},
    "k10": {"n_folds": 10},
}


def _load():
    import importlib.util

    spec = importlib.util.spec_from_file_location("run05", ROOT / "scripts" / "05_run_nsclc.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.load_cohort()


def run_variant(name: str, cohort, X, expr, sigs) -> dict:
    change = dict(VARIANTS[name])
    boot_n = change.pop("_max_null_for_boot", None)
    original = experiment._pooled_isi_bootstrap
    if boot_n is not None:
        experiment._pooled_isi_bootstrap = functools.partial(original, max_null_for_boot=boot_n)
    try:
        within = change.get("axis_within_type", True)
        axis = globalaxis.compute_global_axis(
            expr, method="pc1",
            cancer_type=cohort.set_index("patient_id")["cancer_type"] if within else None)
        cfg = experiment.AuditConfig(seed=0, n_boot=1000, n_null_sets=1000,
                                     **{**ISI_ONLY, **change})
        t0 = time.time()
        res = experiment.run_audit(cohort, X, [c for c in cohort.columns if c.startswith("SIG_")],
                                   cfg, expression=expr, signature_set=sigs,
                                   global_axis=axis, survival=None, verbose=False)
        res.save(OUT / name)
    finally:
        experiment._pooled_isi_bootstrap = original
    est = res.primary_excess
    per = res.per_signature if hasattr(res, "per_signature") else None
    beats = None
    ie = OUT / name / "immune_excess.csv"
    if ie.exists():
        import pandas as pd
        beats = int(pd.read_csv(ie)["beats_null"].sum())
    return {"variant": name, "change": VARIANTS[name], "isi": est.value, "lo": est.lo,
            "hi": est.hi, "signatures_beating_null": beats,
            "seconds": round(time.time() - t0, 1), "has_per_signature": per is not None}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--variant", choices=sorted(VARIANTS), action="append")
    args = ap.parse_args(argv)
    names = args.variant or list(VARIANTS)
    if "reference" not in names:
        names = ["reference"] + names
    OUT.mkdir(parents=True, exist_ok=True)
    cohort, X, expr, sigs = _load()
    ref = json.loads(REFERENCE.read_text())["primary_excess"]
    rows = []
    for name in names:
        row = run_variant(name, cohort, X, expr, sigs)
        row["shift_from_reference"] = row["isi"] - ref["value"]
        rows.append(row)
        print(f"  {name:15s} ISI {row['isi']:.4f} [{row['lo']:.4f}, {row['hi']:.4f}]  "
              f"shift {row['shift_from_reference']:+.4f}  beats {row['signatures_beating_null']}/16  "
              f"({row['seconds']:.0f}s)", flush=True)
    gate = rows[0]
    identical = (gate["isi"], gate["lo"], gate["hi"]) == (ref["value"], ref["ci_lo"], ref["ci_hi"])
    summary = {"reference_source": str(REFERENCE.relative_to(ROOT)),
               "reference_reproduced_bitwise": identical, "variants": rows}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(f"\nreference reproduced bit for bit: {identical}")
    if not identical:
        print("FAILED: the ISI-only reference does not reproduce results/nsclc_v3_stablesort; "
              "no variant here can be compared against it.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
