#!/usr/bin/env python3
"""A9 alpha-margin measurement WITHOUT run_audit. Minutes, not hours.

    python scripts/17_alpha_margin_fast.py

WHY THIS EXISTS
===============
`16_alpha_margin.py` asks the right question -- could a platform-level
floating-point difference plausibly flip any of the 80 ridge-alpha decisions at
all? -- but it answers it by monkeypatching `experiment._select_alpha_per_fold`
and then running the FULL `run_audit`. Measured on macOS 2026-09-02: setup ~16
min, then >3h23m inside `run_audit` without reaching a verdict.

That cost is entirely incidental. Read `experiment.py:1099`:

    fold_alpha = _select_alpha_per_fold(X, obs_resid.to_numpy(), split)

The alpha selection depends on exactly three things:

  * `X`                -- the image embedding, loaded straight off disk
  * `obs_resid`        -- the observed signature score, residualised on the
                          global axis within cancer type
  * `split`            -- the preserved-site fold partition

It does NOT depend on the null gene sets, Cronbach's alpha, split-half
reliability, the four baselines, or the patient-clustered bootstrap -- and those
are the entire multi-hour cost of `run_audit`. So build the three inputs
directly and call the real function. Same 80 numbers, none of the waiting.

WHAT IT REPORTS
===============
`RidgeCV` picks `argmin` over the 24-point grid of leave-one-out CV error, so a
decision can only flip under a different libm/BLAS if the winner and runner-up
are closer than the error those libraries disagree by:

    rel_margin = (mse[runner_up] - mse[winner]) / mse[winner]

  * min rel_margin <= ~1e-12  -> some decision is a coin flip at machine
    precision. A9's alpha hypothesis is LIVE; diff the grid indices across
    machines next.
  * min rel_margin >= ~1e-6   -> no realistic cross-library difference reaches
    the runner-up. The 80 alphas CANNOT differ, the hypothesis is DEAD without
    the cluster, and the next suspects are the ridge solve itself,
    `stats.corr_ci`, and the patient-clustered bootstrap.

    This branch is the one that fired. Measured: min 1.034e-04, median
    5.099e-03, 0/80 below 1e-06.

    NOTE ON WHAT THE BLAS CONTRAST ACTUALLY IS. Earlier versions of this text
    named the next suspect as "Accelerate vs OpenBLAS/MKL". That was wrong and
    is not a live hypothesis: this Mac was measured three independent ways and
    runs **openblas 0.3.21**, not Accelerate, and `20_check_versions.py` now
    asserts it on every run. If HPC4's numpy also ships the OpenBLAS wheel, the
    contrast is one library's arm64 NEON kernels against its own x86_64 AVX
    kernels -- a much narrower question. `18_bisect_platform.py` prints the BLAS
    name on both machines for exactly this reason; read those two lines before
    theorising.

Sanity anchor measured 2026-09-02: on a synthetic 200x30 ridge problem the
best-vs-runner-up relative margin was 3.65e-05. Margins are not automatically
tiny; this test can genuinely come back either way.

FIDELITY: THE ONE DIFFERENCE FROM THE REAL PATH, NOW MEASURED TO BE A NO-OP
===========================================================================
RESOLVED 2026-09-03. This script's full 80-row run printed
`ALPHA DIGEST fc299c9c0f47466e`, which is **the digest `16_alpha_margin.py`
printed over the same 80 alphas** -- a run through `run_audit` itself, taking
19,123 s against this script's ~1 s of setup. Two independent code paths agree
bit-for-bit on every alpha. The difference described below is therefore not a
risk that was argued away; it is a no-op that was measured. The description is
kept because it explains WHY the two paths could in principle have differed.


The alphas come from the REAL `experiment._select_alpha_per_fold`, not a
reimplementation, and the script ASSERTS that its recomputed `argmin`
reproduces that function's choice. If the assertion fails it says so loudly
instead of reporting a meaningless number.

The one deliberate difference: `run_audit` residualises the observed score in a
single `residualise_matrix` call STACKED with that signature's null columns
(experiment.py:1085-1093), whereas this script residualises all 16 observed
columns stacked together. Both are one per-group lstsq against the same design
matrix, and the design matrix -- hence its SVD -- does not depend on the
right-hand side, so the columns are mathematically independent. Different
right-hand-side blocking CAN differ in the last bits. That cannot move a margin
of order 1e-6, and if it moved a margin of order 1e-16 that would itself BE the
coin-flip finding. Stated here rather than buried.

Fixes `16_alpha_margin.py`'s known weakness: this prints one line per signature
as it is recorded, so a live run is distinguishable from a hung one.
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

from aacr27 import experiment, globalaxis, models, splits  # noqa: E402

INTERIM = ROOT / "data" / "interim"
OUT = ROOT / "results" / "alpha_margin_fast_macos.npz"
# The cache stores the patient ids alongside the values. An earlier version
# stored a bare .npy, which silently dropped the axis's patient_id index and
# turned every residual into NaN -- the run burned 791 s and reported "no finite
# margins". A cache that cannot be checked against the cohort is a bug waiting to
# happen, so this one is checked on every load.
AXIS_CACHE = ROOT / "results" / "_cache_global_axis_nsclc.npz"
AXIS_CACHE_LEGACY = ROOT / "results" / "_cache_global_axis_nsclc.npy"


def main() -> int:
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler
    import sklearn

    print(f"python   {sys.version.split()[0]}")
    print(f"numpy    {np.__version__}")
    print(f"sklearn  {sklearn.__version__}")
    print(f"platform {sys.platform}", flush=True)
    print()

    t0 = time.time()
    frame = pd.read_parquet(INTERIM / "cohort_nsclc.parquet")
    X = np.load(INTERIM / "X_nsclc.npy")
    sig_cols = [c for c in frame.columns if c.startswith("SIG_")]
    print(f"cohort {frame.shape}, {len(sig_cols)} signature columns, X {X.shape} "
          f"({time.time() - t0:.1f}s)", flush=True)

    # The global axis is the expensive input (PC1 over 41,046 genes). It depends
    # only on the frozen expression matrix, so cache it -- a rerun of this probe
    # should not pay 16 minutes twice.
    pids = frame["patient_id"].to_numpy()
    if AXIS_CACHE.exists():
        z = np.load(AXIS_CACHE, allow_pickle=False)
        if not np.array_equal(z["patient_id"].astype(str), pids.astype(str)):
            raise SystemExit(
                f"{AXIS_CACHE.name} was built for a different patient order; "
                "delete it and rerun")
        axis_raw = z["axis"]
        print(f"global axis loaded from cache {AXIS_CACHE.name}, patient order "
              f"verified ({time.time() - t0:.1f}s)", flush=True)
    elif AXIS_CACHE_LEGACY.exists():
        # The legacy .npy holds the same values in cohort row order -- the axis
        # was computed on an `expr` already reindexed to frame["patient_id"] --
        # but carries no ids of its own, so the most that can be checked is the
        # length. Upgrade it to the checkable format on the way past.
        axis_raw = np.load(AXIS_CACHE_LEGACY)
        if axis_raw.shape != (len(frame),):
            raise SystemExit(
                f"{AXIS_CACHE_LEGACY.name} has shape {axis_raw.shape}, expected "
                f"({len(frame)},); delete it and rerun")
        np.savez(AXIS_CACHE, axis=axis_raw, patient_id=pids.astype(str))
        print(f"global axis loaded from legacy {AXIS_CACHE_LEGACY.name} "
              f"(length checked) and upgraded to {AXIS_CACHE.name} "
              f"({time.time() - t0:.1f}s)", flush=True)
    else:
        t_axis = time.time()
        expr = pd.read_parquet(INTERIM / "expr_nsclc.parquet")
        if not expr.index.equals(pd.Index(frame["patient_id"])):
            expr = expr.reindex(frame["patient_id"])
        print(f"expression {expr.shape} loaded ({time.time() - t_axis:.1f}s); "
              "computing PC1...", flush=True)
        axis = globalaxis.compute_global_axis(
            expr, method="pc1",
            cancer_type=frame.set_index("patient_id")["cancer_type"])
        # KEEP THE INDEX. `axis.values` is a Series indexed by patient_id, which
        # is what experiment.py:696-698 reindexes against. Casting it to a bare
        # array here and reindexing afterwards is exactly the bug described at
        # AXIS_CACHE above.
        axis_ser = axis.values if hasattr(axis, "values") else pd.Series(axis)
        axis_raw = (pd.Series(axis_ser).reindex(pids).to_numpy(dtype=float)
                    if isinstance(getattr(axis_ser, "index", None), pd.Index)
                    else np.asarray(axis_ser, dtype=float))
        AXIS_CACHE.parent.mkdir(parents=True, exist_ok=True)
        np.savez(AXIS_CACHE, axis=axis_raw, patient_id=pids.astype(str))
        del expr
        print(f"global axis computed in {time.time() - t_axis:.1f}s, "
              f"cached to {AXIS_CACHE.name}", flush=True)

    n_finite_axis = int(np.isfinite(axis_raw).sum())
    print(f"global axis: {n_finite_axis}/{len(frame)} finite", flush=True)
    if n_finite_axis < len(frame) // 2:
        raise SystemExit(
            "the global axis is mostly NaN -- residualising against it would "
            "make every alpha NaN and the margins meaningless")

    # Rebuild axis_values / within / split EXACTLY as experiment.py:696-718 does.
    # `axis_raw` is already in cohort row order (checked above), so it is aligned
    # positionally rather than reindexed by label.
    axis_values = pd.Series(axis_raw, index=frame.index)
    within = frame["cancer_type"]
    split = splits.preserved_site_split(
        frame["patient_id"], frame["tss"].astype("string"), n_folds=5, seed=0,
    )
    print(f"split: {split.n_folds} preserved-site folds, seed 0", flush=True)

    # experiment.py:1085-1093, minus the null columns (see FIDELITY above).
    obs_raw = np.column_stack(
        [pd.to_numeric(frame[c], errors="coerce").to_numpy(dtype=float)
         for c in sig_cols])
    obs_resid = globalaxis.residualise_matrix(
        obs_raw,
        pd.to_numeric(axis_values, errors="coerce").to_numpy(dtype=float),
        groups=within.to_numpy(),
    )
    # FAIL FAST. `_select_alpha_per_fold` silently returns NaN for any fold with
    # fewer than 20 finite targets, so an all-NaN `y` produces 80 NaN alphas, no
    # margins, and no error -- 791 s to discover a one-line alignment bug. Check
    # here, where it costs nothing.
    finite_frac = float(np.isfinite(obs_resid).mean())
    print(f"residualised {obs_resid.shape} observed scores, "
          f"{100 * finite_frac:.1f}% finite ({time.time() - t0:.1f}s)\n",
          flush=True)
    if finite_frac < 0.5:
        raise SystemExit(
            f"only {100 * finite_frac:.1f}% of residualised scores are finite; "
            "every alpha would be NaN and every margin meaningless")

    grid = models.ALPHAS
    records: list[dict] = []
    mismatches: list[str] = []

    print("  signature                       fold  idx  runner  "
          "      alpha   rel_margin", flush=True)
    for s, sig in enumerate(sig_cols):
        y = obs_resid[:, s]
        # Ground truth: the real function, called exactly as run_audit calls it.
        alphas_real = np.asarray(
            experiment._select_alpha_per_fold(X, y, split), dtype=float)

        for k in range(split.n_folds):
            alpha_real = float(alphas_real[k])
            train_idx, _ = split.indices(k)
            ok = train_idx[np.isfinite(y[train_idx])]
            if len(ok) < 20 or not np.isfinite(alpha_real):
                records.append(dict(sig=sig, fold=k, alpha=np.nan, idx=-1,
                                    runner_up=-1, rel_margin=np.nan))
                print(f"  {sig[:30]:30s}  {k:4d}    -       -            nan"
                      "          nan", flush=True)
                continue
            # Mirror _select_alpha_per_fold's preprocessing exactly.
            med = models.fit_impute_median(X[ok])
            Xs = StandardScaler().fit_transform(models.impute_median(X[ok], med))
            fit = RidgeCV(alphas=grid, store_cv_results=True).fit(Xs, y[ok])
            mse = np.asarray(fit.cv_results_).mean(axis=0)

            order = np.argsort(mse)
            win, second = int(order[0]), int(order[1])
            idx_real = int(np.argmin(np.abs(grid - alpha_real)))
            if win != idx_real:
                mismatches.append(
                    f"{sig} fold={k}: recomputed argmin {win} != real alpha "
                    f"index {idx_real}")
            rel = (float((mse[second] - mse[win]) / mse[win])
                   if mse[win] > 0 else np.nan)
            records.append(dict(sig=sig, fold=k, alpha=alpha_real, idx=idx_real,
                                runner_up=second, rel_margin=rel))
            print(f"  {sig[:30]:30s}  {k:4d}  {idx_real:3d}  {second:6d}  "
                  f"{alpha_real:11.4f}  {rel:11.3e}", flush=True)

    df = pd.DataFrame(records)
    print(f"\n{df['sig'].nunique()} signatures x {split.n_folds} folds = "
          f"{len(df)} selections, in {time.time() - t0:.1f}s total", flush=True)

    if mismatches:
        print(f"\n!! {len(mismatches)} recomputed argmin(s) disagree with the real "
              "function -- THE MARGINS ARE NOT TRUSTWORTHY")
        for m in mismatches[:10]:
            print(f"     {m}")
    else:
        print("\n  recomputed argmin reproduces the real alpha in every "
              "finite selection")

    finite = df[np.isfinite(df["rel_margin"])]
    m = finite["rel_margin"].to_numpy()
    if m.size == 0:
        print("no finite margins -- nothing to report")
        return 1

    print(f"\nRELATIVE MARGIN over {m.size} decisions:")
    print(f"  min {m.min():.3e}   p10 {np.percentile(m, 10):.3e}   "
          f"median {np.median(m):.3e}   max {m.max():.3e}")
    for thr in (1e-14, 1e-12, 1e-9, 1e-6):
        print(f"  below {thr:.0e}: {int((m < thr).sum())} / {m.size}")

    print("\nVERDICT")
    if m.min() < 1e-12:
        worst = finite.loc[finite["rel_margin"].idxmin()]
        print("  At least one decision sits within ~1e-12 relative of its runner-up.")
        print("  A platform-level FP difference CAN flip it. A9's alpha hypothesis")
        print("  SURVIVES this test -- diff the grid indices across machines next.")
        print(f"  Closest: {worst['sig']} fold {int(worst['fold'])}, "
              f"grid idx {int(worst['idx'])} vs runner-up {int(worst['runner_up'])}, "
              f"margin {worst['rel_margin']:.3e}")
    else:
        print(f"  The closest decision is {m.min():.3e} from its runner-up, which is")
        print("  far outside what differing libm/BLAS implementations disagree by")
        print("  (~1e-12 relative on this kind of reduction).")
        print("  The 80 alphas CANNOT differ across platforms, so they are NOT the")
        print("  cause of the 0.0218 NSCLC gap. A9's leading hypothesis is DEAD.")
        print("  Next suspects, in order: the ridge solve itself, stats.corr_ci,")
        print("  the patient-clustered bootstrap. Discriminate them with")
        print("  18_bisect_platform.py, which dumps raw float64 per stage.")
        print("  NOT 'Accelerate vs OpenBLAS' -- this Mac was measured three ways")
        print("  and runs openblas 0.3.21. Read the BLAS line 18 prints on each")
        print("  machine before theorising about which library is which.")

    # Same digest 15_alpha_probe.py / 16_alpha_margin.py print, over the same
    # values in the same order, so the three can be cross-checked. NOTE: the
    # signature ITERATION ORDER must match for the digest to be comparable.
    flat = df["alpha"].to_numpy(dtype=float)
    digest = hashlib.sha256(np.ascontiguousarray(flat).tobytes()).hexdigest()[:16]
    print(f"\nALPHA DIGEST (all {flat.size} values): {digest}")
    print("grid indices are what matter: a differing index IS a discrete flip.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, sig=df["sig"].to_numpy().astype(str), fold=df["fold"].to_numpy(),
             alpha=df["alpha"].to_numpy(), idx=df["idx"].to_numpy(),
             runner_up=df["runner_up"].to_numpy(),
             rel_margin=df["rel_margin"].to_numpy())
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
