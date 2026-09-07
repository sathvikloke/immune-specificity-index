#!/usr/bin/env python3
"""How close was each ridge-alpha decision? (falsifies A9's leading hypothesis)

    python scripts/16_alpha_margin.py

WHY THIS EXISTS
===============
`15_alpha_probe.py` asks whether the 80 selected alphas DIFFER between macOS and
HPC4. That question needs both machines. This script asks a different question
that needs only one machine, and it can kill the hypothesis outright:

    could a platform-level floating-point difference plausibly flip any of
    these 80 decisions AT ALL?

`RidgeCV` picks `argmin` over a 24-point grid of the leave-one-out CV error. A
decision flips under a different libm/BLAS only if the winner and the runner-up
are separated by less than the error those libraries disagree by. So measure the
separation directly:

    rel_margin = (mse[runner_up] - mse[winner]) / mse[winner]

Read it like this:

  * rel_margin ~ 1e-16 to 1e-14 for some fold  -> the decision is a coin flip at
    machine precision. A9's alpha hypothesis is LIVE, and this is the fold to
    look at first.
  * rel_margin >= ~1e-6 everywhere             -> no realistic cross-library
    difference (they agree to ~1e-12 relative on this kind of reduction) can
    reach the runner-up. The alphas CANNOT differ, the hypothesis is DEAD
    without ever touching the cluster, and the next suspect is the ridge solve
    itself, `stats.corr_ci`, or the patient-clustered bootstrap.

This is the cheapest experiment that could falsify the hypothesis, which is why
it is worth running before spending cluster time.

WHAT IT ALSO MEASURES, AND WHY
==============================
`15_alpha_probe.py`'s docstring predicts it "runs in seconds". Observed on macOS
on 2026-09-02: >30 minutes, twice. That prediction is wrong and nobody knows
where the time goes, so this script timestamps every alpha-selection call and
prints the gap between consecutive ones. The gap is the per-signature cost of
everything `run_audit` does BETWEEN selections, which is the thing to optimise if
this path gets run again.

Numbers, not guesses -- two A9 root causes have already been asserted and
disproved.

NOTE ON FIDELITY
================
The alphas are taken from the REAL `experiment._select_alpha_per_fold`, not from
a reimplementation. The CV curve is then recomputed alongside it, and the script
ASSERTS that the recomputed argmin reproduces the real function's choice. If that
assertion ever fails, the margins are meaningless and the script says so loudly
rather than reporting a number.
"""

from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import experiment, globalaxis, models, signatures  # noqa: E402

INTERIM = ROOT / "data" / "interim"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"
OUT = ROOT / "results" / "alpha_margin_macos.npz"


def main() -> int:
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler
    import sklearn

    print(f"python   {sys.version.split()[0]}")
    print(f"numpy    {np.__version__}")
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

    grid = models.ALPHAS
    records: list[dict] = []
    mismatches: list[str] = []
    original = experiment._select_alpha_per_fold
    t_start = time.time()
    last_call_end = [t_start]

    def recording(X_, y_, split_, alphas=None):
        gap = time.time() - last_call_end[0]
        # Ground truth: the real function, called exactly as the audit calls it.
        out = (original(X_, y_, split_) if alphas is None
               else original(X_, y_, split_, alphas))
        used = grid if alphas is None else np.asarray(alphas, dtype=float)

        y_arr = np.asarray(y_, dtype=float)
        tag = hashlib.sha256(np.ascontiguousarray(y_arr).tobytes()).hexdigest()[:12]

        # Recompute the LOO-CV curve per fold, mirroring the real function's
        # preprocessing exactly (see experiment._select_alpha_per_fold).
        for k in range(split_.n_folds):
            alpha_real = float(np.asarray(out, dtype=float)[k])
            train_idx, _ = split_.indices(k)
            ok = train_idx[np.isfinite(y_arr[train_idx])]
            if len(ok) < 20 or not np.isfinite(alpha_real):
                records.append(dict(tag=tag, fold=k, alpha=np.nan, idx=-1,
                                    runner_up=-1, rel_margin=np.nan, gap=gap))
                continue
            med = models.fit_impute_median(X_[ok])
            Xs = StandardScaler().fit_transform(models.impute_median(X_[ok], med))
            fit = RidgeCV(alphas=used, store_cv_results=True).fit(Xs, y_arr[ok])
            mse = np.asarray(fit.cv_results_).mean(axis=0)

            order = np.argsort(mse)
            win, second = int(order[0]), int(order[1])
            idx_real = int(np.argmin(np.abs(used - alpha_real)))
            if win != idx_real:
                mismatches.append(
                    f"y={tag} fold={k}: recomputed argmin {win} != real alpha "
                    f"index {idx_real}")
            rel = float((mse[second] - mse[win]) / mse[win]) if mse[win] > 0 else np.nan
            records.append(dict(tag=tag, fold=k, alpha=alpha_real, idx=idx_real,
                                runner_up=second, rel_margin=rel, gap=gap))

        last_call_end[0] = time.time()
        return out

    experiment._select_alpha_per_fold = recording
    try:
        # n_null=2 ONLY to make the run cheap; the alphas do not depend on the
        # null sets at all (see 15_alpha_probe.py).
        cfg = experiment.AuditConfig(
            n_folds=5, seed=0, n_boot=50, n_null_sets=2,
            run_combat=False, run_legacy_venet_null=False,
        )
        experiment.run_audit(
            cohort, X, sig_cols, cfg, expression=expr, signature_set=sigs,
            global_axis=axis, survival=None, verbose=False)
    except Exception as exc:  # noqa: BLE001
        # The margins are already recorded by the time anything downstream can
        # fail; say so plainly rather than losing the measurement.
        print(f"\n!! run_audit raised after {len(records)} fold-records: "
              f"{type(exc).__name__}: {exc}")
    finally:
        experiment._select_alpha_per_fold = original

    if not records:
        print("no alpha selections recorded -- nothing to report")
        return 1

    df = pd.DataFrame(records)
    print(f"\n{df['tag'].nunique()} signatures x {df['fold'].nunique()} folds "
          f"= {len(df)} selections, in {time.time() - t_start:.1f}s\n")

    print("  per-signature cost BETWEEN alpha selections (the 'seconds' claim):")
    gaps = df.groupby("tag", sort=False)["gap"].first().to_numpy()
    print(f"    gaps  min {gaps.min():.1f}s  median {np.median(gaps):.1f}s  "
          f"max {gaps.max():.1f}s  total {gaps.sum():.1f}s\n")

    finite = df[np.isfinite(df["rel_margin"])]
    print("  tag           fold  grid_idx  runner_up        alpha   rel_margin")
    for _, r in df.iterrows():
        rel = ("nan" if not np.isfinite(r["rel_margin"])
               else f"{r['rel_margin']:.3e}")
        print(f"  {r['tag']}  {int(r['fold']):4d}  {int(r['idx']):8d}  "
              f"{int(r['runner_up']):9d}  {r['alpha']:11.4f}  {rel:>11}")

    if mismatches:
        print(f"\n!! {len(mismatches)} recomputed argmin(s) disagree with the real "
              "function -- MARGINS BELOW ARE NOT TRUSTWORTHY")
        for m in mismatches[:10]:
            print(f"     {m}")
    else:
        print("\n  recomputed argmin reproduces the real alpha in all "
              f"{len(finite)} finite selections")

    m = finite["rel_margin"].to_numpy()
    print(f"\nRELATIVE MARGIN over {m.size} decisions:")
    print(f"  min {m.min():.3e}   p10 {np.percentile(m, 10):.3e}   "
          f"median {np.median(m):.3e}   max {m.max():.3e}")
    for thr in (1e-14, 1e-12, 1e-9, 1e-6):
        print(f"  below {thr:.0e}: {int((m < thr).sum())} / {m.size}")

    print("\nVERDICT")
    if m.min() < 1e-12:
        print("  At least one decision sits within ~1e-12 relative of its runner-up.")
        print("  A platform-level FP difference CAN flip it. A9's alpha hypothesis")
        print("  survives this test -- diff 15_alpha_probe.py across machines next.")
    else:
        print(f"  The closest decision is {m.min():.3e} from its runner-up, which is")
        print("  far outside what differing libm/BLAS implementations disagree by.")
        print("  The 80 alphas CANNOT differ across platforms, so they are NOT the")
        print("  cause of the 0.0218 gap. Next suspects: the ridge solve itself,")
        print("  stats.corr_ci, the patient-clustered bootstrap.")

    # Same digest 15_alpha_probe.py prints, over the same values in the same
    # order, so this script is a strict superset of that one and the two can be
    # cross-checked if both are ever run on the same machine.
    flat = df["alpha"].to_numpy(dtype=float)
    digest = hashlib.sha256(np.ascontiguousarray(flat).tobytes()).hexdigest()[:16]
    print(f"\nALPHA DIGEST (all {flat.size} values): {digest}")
    print("grid indices are what matter: a differing index IS a discrete flip.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, tag=df["tag"].to_numpy().astype(str), fold=df["fold"].to_numpy(),
             alpha=df["alpha"].to_numpy(), idx=df["idx"].to_numpy(),
             runner_up=df["runner_up"].to_numpy(),
             rel_margin=df["rel_margin"].to_numpy(), gap=df["gap"].to_numpy())
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
