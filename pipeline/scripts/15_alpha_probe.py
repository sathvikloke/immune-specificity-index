#!/usr/bin/env python3
"""Are the per-fold ridge alphas platform-stable? (the 0.3182 vs 0.2964 hunt)

    python scripts/15_alpha_probe.py

WHY THIS EXISTS, AND WHY 14_repro_probe.py WAS NOT ENOUGH
========================================================
Job 117403 established the fact that has to be explained: with BIT-IDENTICAL
inputs (14_repro_probe.py stages 1-6 hash the same), the same pinned versions
(numpy 2.1.3, pandas 2.2.3, sklearn 1.6.1, scipy 1.15.3), the same seeds and the
STORED cohort columns -- no rescoring anywhere -- the NSCLC ISI is

    macOS / arm64    0.318240180054566
    HPC4  / x86_64   0.2964

`14_repro_probe.py` stage 7 fits RidgeCV ONCE on the full cohort over the 16
observed signatures and hashes the selected alphas. They matched, and that was
briefly over-read as clearing the ridge. It does not clear it: that is 16
selections on the easiest possible fit, and it is not the path the audit takes.

THE PATH THE AUDIT ACTUALLY TAKES
=================================
`experiment._select_alpha_per_fold` is called once PER SIGNATURE, and returns one
scalar per fold, each chosen by `RidgeCV(alphas=ALPHAS)` on that fold's training
rows only. So the real run makes 16 x 5 = 80 scalar selections off a 24-point
grid, and each one is then held FIXED across that signature's observed target and
all 1,000 null targets (`cross_val_predict_multi(..., fixed_alpha=fold_alpha)`).

That makes these 80 numbers the highest-leverage discrete branch in the pipeline.
`ALPHAS = np.logspace(-2, 5, 24)`, so neighbours differ by a factor of ~1.96: one
flip does not perturb a prediction, it roughly halves or doubles the shrinkage for
one signature-fold, for the observed AND the null side at once. A handful of flips
across 80 selections is entirely capable of moving a pooled estimate by 0.022,
whereas ordinary libm noise would move it by ~1e-12. The size of the shift is
itself the argument that a discrete branch flipped.

WHY THIS IS CHEAP
=================
`fold_alpha` is selected on `obs_resid` -- the residualised OBSERVED score -- and
on the fold geometry. It does not depend on the null sets at all. Residualisation
is column-wise, so column 0's residual does not depend on the other columns
either. Therefore the 80 alphas at n_null=2 are IDENTICAL to the 80 alphas at
n_null=1000, and this probe runs in seconds instead of six minutes.

Run on both machines and diff. If the alphas differ, the reproducibility defect is
located. If they are identical, the cause is downstream of alpha selection and the
next suspect is the ridge solve itself.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import experiment, globalaxis, signatures  # noqa: E402

INTERIM = ROOT / "data" / "interim"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"


def main() -> int:
    print(f"python   {sys.version.split()[0]}")
    print(f"numpy    {np.__version__}")
    import sklearn
    print(f"sklearn  {sklearn.__version__}")
    print(f"platform {sys.platform}")
    print()

    cohort = pd.read_parquet(INTERIM / "cohort_nsclc.parquet")
    expr = pd.read_parquet(INTERIM / "expr_nsclc.parquet")
    X = np.load(INTERIM / "X_nsclc.npy")
    if not expr.index.equals(pd.Index(cohort["patient_id"])):
        expr = expr.reindex(cohort["patient_id"])

    sigs = signatures.SignatureSet.from_gmt(GMT)
    sigs = signatures.SignatureSet(
        name=sigs.name, sets={f"SIG_{k}": v for k, v in sigs.sets.items()})
    axis = globalaxis.compute_global_axis(
        expr, method="pc1",
        cancer_type=cohort.set_index("patient_id")["cancer_type"])

    sig_cols = [c for c in cohort.columns if c.startswith("SIG_")]
    print(f"cohort {cohort.shape}, {len(sig_cols)} signature columns, X {X.shape}")

    # ---- record every _select_alpha_per_fold return, in call order ----------
    calls: list[tuple[str, np.ndarray]] = []
    original = experiment._select_alpha_per_fold

    def recording(X_, y_, split_, alphas=None):
        out = (original(X_, y_, split_) if alphas is None
               else original(X_, y_, split_, alphas))
        # Tag by the target vector so the diff survives any change in call order.
        tag = hashlib.sha256(
            np.ascontiguousarray(np.asarray(y_, dtype=float)).tobytes()
        ).hexdigest()[:12]
        calls.append((tag, np.asarray(out, dtype=float).copy()))
        return out

    experiment._select_alpha_per_fold = recording
    try:
        # n_null=2 ONLY to make the run cheap; the alphas do not depend on it.
        cfg = experiment.AuditConfig(
            n_folds=5, seed=0, n_boot=50, n_null_sets=2,
            run_combat=False, run_legacy_venet_null=False,
        )
        experiment.run_audit(
            cohort, X, sig_cols, cfg, expression=expr, signature_set=sigs,
            global_axis=axis, survival=None, verbose=False)
    finally:
        experiment._select_alpha_per_fold = original

    print(f"\n{len(calls)} alpha-selection calls recorded "
          f"({len(calls)} x {len(calls[0][1]) if calls else 0} folds)\n")

    grid = np.logspace(-2, 5, 24)
    for i, (tag, a) in enumerate(calls):
        idx = [int(np.argmin(np.abs(grid - v))) if np.isfinite(v) else -1 for v in a]
        pretty = "  ".join("nan" if not np.isfinite(v) else f"{v:11.4f}" for v in a)
        print(f"  call {i:3d}  y={tag}  grid_idx={idx}  {pretty}")

    flat = np.concatenate([a for _, a in calls]) if calls else np.array([])
    digest = hashlib.sha256(
        np.ascontiguousarray(flat).tobytes()).hexdigest()[:16]
    print(f"\nALPHA DIGEST (all {flat.size} values): {digest}")
    print("grid indices are what matter: a differing index IS a discrete flip.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
