#!/usr/bin/env python3
"""A9 stage-by-stage bisection of the 0.0218 macOS-vs-HPC4 ISI gap.

    # on each machine
    python scripts/18_bisect_platform.py --tag macos
    python scripts/18_bisect_platform.py --tag hpc4

    # then, with both .npz files on one machine
    python scripts/18_bisect_platform.py --compare \\
        results/bisect_macos.npz results/bisect_hpc4.npz

WHY THIS EXISTS, AND WHAT IT REPLACED
=====================================
The NSCLC ISI is 0.318240180054566 on macOS 26.5.2 / arm64 / Python 3.13.9 and
0.2964 on Einstein HPC4 / x86_64 / Python 3.12.14, from the same code, the same
seeds and the same pinned versions (numpy 2.1.3, pandas 2.2.3, sklearn 1.6.1,
scipy 1.15.3). `14_repro_probe.py` showed stages 1-6 -- input hashes, gene-axis
order, mean expression, qcut bins, sampled null sets, mean-z scores -- are
IDENTICAL on both machines. So the divergence enters somewhere after scoring.

The leading suspect was the ridge alpha: `_select_alpha_per_fold` makes 16
signatures x 5 folds = 80 argmin decisions over a 24-point grid whose adjacent
points differ by ~1.96x, and one flip roughly halves or doubles the shrinkage
for a signature-fold. `16_alpha_margin.py` MEASURED that on macOS 2026-09-02,
19,123 s, all 80 selections:

    min 1.034e-04   p10 7.917e-04   median 5.099e-03   max 1.125e-02
    below 1e-06: 0 / 80

The closest decision is 1e-4 from its runner-up -- roughly eleven orders of
magnitude further than differing libm/BLAS implementations disagree by.

FROM THAT, SIX SESSIONS CONCLUDED "the 80 alphas CANNOT differ across platforms".
THAT CONCLUSION WAS WRONG, and this paragraph is kept as the record of it.

The margins above are correctly measured. The inference from them is not, for two
reasons that only became visible on 2026-09-06 when this script was finally run
on both machines. First, both runs that "independently confirmed" the invariance
were on the SAME machine, so what they established was that the fast path is
faithful -- not that the result is platform-invariant. Second, a wide margin only
implies invariance if nothing between the inputs and the decision amplifies; the
amplification here is enormous. Measured: 11 of the 80 alphas differ, each by
exactly one grid step, and the digests are fc299c9c0f47466e (macOS) against
946033557a8ad08b (HPC4).

AND THE ALPHA IS NOT THE ROOT CAUSE EITHER. Eight of the sixteen signatures have
IDENTICAL alphas on both machines, and their out-of-fold predictions still
diverge -- by up to 3.023e-01, 177% of those predictions' own standard deviation,
correlating at only r ~= 0.96 -- from inputs that agree to 9.77e-15. So the ridge
SOLVE differs at identical penalties on identical inputs, and the alpha flips are
a symptom of the same instability rather than its origin. `--compare` now runs
that control automatically instead of asserting the amplification story; see
`_amplification_check`.

The structural fact that most likely explains it, measured rather than assumed:
X is (944, 768) and the split is 5 preserved-site folds, so every training set
is about 755 samples against 768 features -- fewer rows than columns. What is NOT
yet established is which part of the solve turns that into an O(0.1) difference.
That is the open question; do not write down an answer to it that has not been
run.

THE THREE REMAINING SUSPECTS
============================
Reading `experiment._residualised_null_scores` (lines 1099-1140 at `8fe1187`),
everything after the alpha is:

  S4  models.cross_val_predict_multi(X, stacked, split, patients,
                                     fixed_alpha=fold_alpha)
      The ridge SOLVE -- the only stage running a large dense factorisation, so
      the only one where a BLAS difference has room to accumulate.

      MEASURED 2026-09-02, and it corrects the standing assumption: numpy on
      this Mac reports `BLAS openblas 0.3.21`, NOT Accelerate. So if HPC4's
      numpy also ships the OpenBLAS wheel, the contrast is OpenBLAS's arm64 NEON
      kernels against its x86_64 AVX kernels -- same library, different
      hand-written kernels and different blocking -- not two separate libraries.
      This script PRINTS the BLAS name on both machines so the assumption is
      checked rather than carried. Read those two lines before theorising.

  S5  stats.corr_ci(...).value  -- observed r, and 1000 null r values.
      Wraps scipy.stats.pearsonr. Amplifies whatever S4 hands it; a divergence
      that FIRST appears here rather than in S4 would implicate scipy, not BLAS.

  S6/S7  Fisher z, disattenuation, and _pooled_isi_bootstrap.
      The bootstrap resamples patients with a seeded Generator. PCG64 is
      bit-exact across platforms, so the DRAWS cannot differ -- but the
      per-draw arithmetic can, and the ISI point estimate is a pooled mean over
      quantities the earlier stages produced.

HOW TO READ THE OUTPUT
======================
`--compare` prints, per stage, the max absolute difference between the two
machines' dumps. The stage where that difference FIRST leaves the ~1e-15
round-off floor is where the divergence enters.

It is tempting to add "and everything after it is downstream amplification, not
a separate bug" -- this docstring did say exactly that, and it is false here.
`--compare` therefore ends with an AMPLIFICATION CHECK that tests the claim
against the signatures whose alphas did NOT flip. Read that block; it is the
part that distinguishes one defect from two.

  * S4 diverges  -> the ridge solve. Expect ~1e-13 there growing to ~1e-2 by
    S7, i.e. a conditioning/amplification story, not a coding error.
  * S4 clean, S5 diverges -> scipy's pearsonr, not the BLAS.
  * S4-S6 clean, S7 diverges -> the bootstrap. S0 below isolates that case.

S0 is a STANDALONE bootstrap probe on deterministic synthetic input, so the
bootstrap can be diffed WITHOUT depending on any upstream stage. If S0 differs,
the bootstrap is implicated on its own; if S0 matches and S7 differs, the
bootstrap is only propagating what S4-S6 handed it.

COST, AND WHY --n-null IS SMALL BY DEFAULT
==========================================
This is a bisection, not an estimate. A platform divergence in the ridge solve
shows up in the FIRST signature with 8 null columns exactly as clearly as with
1000, and 1000 is what makes `run_audit` a multi-hour job. `--n-null 8` keeps
the whole thing to 16 alpha selections' worth of RidgeCV plus trivial extras.

The alpha selection is still the dominant cost -- 5 RidgeCV fits per signature,
nominally ~1 s each but 50-100x slower on a contended Mac. Use `--signatures N`
to bisect on the first N signatures only; one is enough to localise a stage.

MEASURED, not estimated: `--n-null 4 --signatures 1` cost 540.5 s per signature
on macOS at load average ~215 with two other analysis jobs running, plus 3.5 s
of setup. So all 16 signatures is ~2.4 h under that contention and far less on
an idle machine or on HPC4.

MODELLED ON 17_alpha_margin_fast.py: minimal inputs built directly rather than
via `run_audit`, the REAL library functions at every stage, fail-fast assertions
on input shape (an all-NaN target once produced 80 silent NaN alphas and cost
791 s to diagnose), and one progress line per signature so a live run is
distinguishable from a hung one.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import experiment, globalaxis, models, signatures as sig_mod, splits, stats  # noqa: E402

INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"
AXIS_CACHE = RESULTS / "_cache_global_axis_nsclc.npz"
AXIS_CACHE_LEGACY = RESULTS / "_cache_global_axis_nsclc.npy"

# The stages, in pipeline order. `--compare` walks this list and reports the
# first one to leave the round-off floor.
STAGES = [
    ("S0_boot_synth", "standalone bootstrap on synthetic input"),
    ("S2_obs_resid", "global-axis-residualised observed score"),
    ("S2_null_resid", "global-axis-residualised null scores"),
    ("S3_fold_alpha", "per-fold ridge alpha -- 11 of 80 differ across platforms"),
    ("S4_pred_obs", "out-of-fold ridge predictions, observed"),
    ("S4_pred_null", "out-of-fold ridge predictions, null"),
    ("S5_r_obs", "stats.corr_ci observed r"),
    ("S5_r_null", "stats.corr_ci null r"),
    ("S6_z_excess", "Fisher-z excess per signature (the per-signature ISI term)"),
]

# Anything at or below this is float64 round-off from a different summation
# order, not a divergence. Two independent BLAS implementations agree on a
# well-conditioned dense solve to about here.
FLOOR = 1e-12


def _load_axis(frame: pd.DataFrame, t0: float) -> np.ndarray:
    """The PC1 global axis, cached. Same contract as 17_alpha_margin_fast.py.

    The cache carries the patient ids and is CHECKED against the cohort on every
    load: an earlier bare-.npy cache silently dropped the axis's patient_id
    index, turned every residual into NaN, and produced 80 NaN alphas with no
    error.
    """
    pids = frame["patient_id"].to_numpy()
    if AXIS_CACHE.exists():
        z = np.load(AXIS_CACHE, allow_pickle=False)
        if not np.array_equal(z["patient_id"].astype(str), pids.astype(str)):
            raise SystemExit(f"{AXIS_CACHE.name} was built for a different "
                             "patient order; delete it and rerun")
        print(f"global axis from cache, patient order verified "
              f"({time.time() - t0:.1f}s)", flush=True)
        return z["axis"]
    if AXIS_CACHE_LEGACY.exists():
        axis_raw = np.load(AXIS_CACHE_LEGACY)
        if axis_raw.shape != (len(frame),):
            raise SystemExit(f"{AXIS_CACHE_LEGACY.name} has shape "
                             f"{axis_raw.shape}, expected ({len(frame)},)")
        np.savez(AXIS_CACHE, axis=axis_raw, patient_id=pids.astype(str))
        print(f"global axis from legacy .npy (length checked), upgraded to "
              f".npz ({time.time() - t0:.1f}s)", flush=True)
        return axis_raw

    t_axis = time.time()
    expr = pd.read_parquet(INTERIM / "expr_nsclc.parquet")
    if not expr.index.equals(pd.Index(frame["patient_id"])):
        expr = expr.reindex(frame["patient_id"])
    axis = globalaxis.compute_global_axis(
        expr, method="pc1",
        cancer_type=frame.set_index("patient_id")["cancer_type"])
    axis_ser = axis.values if hasattr(axis, "values") else pd.Series(axis)
    # KEEP THE INDEX -- see the cache note above.
    axis_raw = (pd.Series(axis_ser).reindex(pids).to_numpy(dtype=float)
                if isinstance(getattr(axis_ser, "index", None), pd.Index)
                else np.asarray(axis_ser, dtype=float))
    np.savez(AXIS_CACHE, axis=axis_raw, patient_id=pids.astype(str))
    del expr
    print(f"global axis computed in {time.time() - t_axis:.1f}s and cached",
          flush=True)
    return axis_raw


def _stage0_bootstrap_synthetic() -> np.ndarray:
    """Isolate the patient-clustered bootstrap from every upstream stage.

    Feeds `_pooled_isi_bootstrap` deterministic synthetic inputs built from a
    seeded Generator, so its output depends on NOTHING this pipeline computed.
    If S0 differs across platforms the bootstrap is implicated on its own; if S0
    matches and S7 differs, the bootstrap is only propagating S4-S6.

    PCG64 is bit-exact across platforms by specification, so the resampled
    INDICES are guaranteed identical and any difference here is arithmetic.
    """
    rng = np.random.default_rng(20260902)
    n, n_sig, b = 300, 4, 32
    boot_inputs = []
    point_values = []
    for k in range(n_sig):
        y = rng.standard_normal(n)
        p = 0.4 * y + rng.standard_normal(n)
        Yn = rng.standard_normal((n, b))
        Pn = 0.05 * Yn + rng.standard_normal((n, b))
        # Keys verified against `experiment._immune_specific_excess` (lines
        # 826-831 at `8fe1187`), NOT guessed.
        boot_inputs.append({
            "signature": f"SYNTH_{k}",
            "y": y, "p": p, "Yn": Yn, "Pn": Pn,
            # `rel_null` is indexed COLUMN-WISE in `experiment._pooled_isi_bootstrap`
            # (line 988 at `8fe1187`), so it is a
            # per-null-draw array, not a scalar. Passing a scalar raises
            # IndexError deep inside the bootstrap.
            "rel_obs": 0.9, "rel_null": np.full(b, 0.8),
        })
        point_values.append(
            stats.fisher_z(stats.corr_ci(y, p).value)
            - float(np.mean(stats.fisher_z(
                [stats.corr_ci(Yn[:, j], Pn[:, j]).value for j in range(b)])))
        )

    cfg = experiment.AuditConfig(n_boot=200, seed=0)
    notes: list[str] = []
    try:
        est = experiment._pooled_isi_bootstrap(
            np.asarray(point_values, dtype=float), boot_inputs, cfg, notes)
    except Exception as exc:  # noqa: BLE001 -- an internal contract, not an API
        # `boot_inputs` is an internal shape and a future refactor may change it.
        # Say so loudly rather than dumping a fabricated stage.
        print(f"  S0 SKIPPED: synthetic boot_inputs rejected "
              f"({type(exc).__name__}: {exc}). The bootstrap cannot be isolated "
              "on its own; read the S6 pooled value instead.", flush=True)
        return np.array([np.nan, np.nan, np.nan])
    return np.array([est.value, est.lo, est.hi], dtype=float)


def run(tag: str, n_null: int, n_sig_limit: int | None) -> int:
    import sklearn
    import scipy

    print(f"python   {sys.version.split()[0]}")
    print(f"numpy    {np.__version__}")
    print(f"scipy    {scipy.__version__}")
    print(f"sklearn  {sklearn.__version__}")
    print(f"platform {sys.platform} / {__import__('platform').machine()}")
    # The BLAS name is the whole point of the S4 suspicion, so print it -- but
    # numpy's config API has changed shape across versions, so never let a
    # failure to introspect it kill the run.
    try:
        blas = np.__config__.show(mode="dicts")["Build Dependencies"]["blas"]
        print(f"BLAS     {blas.get('name')} {blas.get('version', '')}")
    except Exception as exc:  # noqa: BLE001
        print(f"BLAS     (could not introspect: {type(exc).__name__})")
    print(flush=True)

    t0 = time.time()
    frame = pd.read_parquet(INTERIM / "cohort_nsclc.parquet")
    X = np.load(INTERIM / "X_nsclc.npy")
    sig_cols = [c for c in frame.columns if c.startswith("SIG_")]
    if n_sig_limit:
        sig_cols = sig_cols[:n_sig_limit]
    print(f"cohort {frame.shape}, {len(sig_cols)} signatures, X {X.shape} "
          f"({time.time() - t0:.1f}s)", flush=True)

    axis_raw = _load_axis(frame, t0)
    n_finite = int(np.isfinite(axis_raw).sum())
    print(f"global axis: {n_finite}/{len(frame)} finite", flush=True)
    if n_finite < len(frame) // 2:
        raise SystemExit("the global axis is mostly NaN -- every downstream "
                         "stage would be meaningless")

    axis_values = pd.to_numeric(pd.Series(axis_raw, index=frame.index),
                                errors="coerce").to_numpy(dtype=float)
    within = frame["cancer_type"].to_numpy()
    split = splits.preserved_site_split(
        frame["patient_id"], frame["tss"].astype("string"), n_folds=5, seed=0)
    patients = frame["patient_id"].to_numpy()
    print(f"split: {split.n_folds} preserved-site folds, seed 0", flush=True)

    # The null sets need the expression matrix. `14_repro_probe.py` already
    # showed the SAMPLED SETS hash identically on both machines, so this stage
    # is carried for completeness rather than as a suspect.
    t_expr = time.time()
    expr = pd.read_parquet(INTERIM / "expr_nsclc.parquet")
    if not expr.index.equals(pd.Index(frame["patient_id"])):
        expr = expr.reindex(frame["patient_id"])
    print(f"expression {expr.shape} ({time.time() - t_expr:.1f}s)", flush=True)
    cfg = experiment.AuditConfig(n_null_sets=n_null, seed=0)
    # The SIG_ prefixing convention, from `05_run_nsclc.py`'s `load_cohort`
    # (lines 48-53 at `8fe1187`): gene-SET names
    # in the GMT are bare, score COLUMNS in the cohort frame carry `SIG_`.
    raw_sigs = sig_mod.SignatureSet.from_gmt(GMT)
    sigs = sig_mod.SignatureSet(
        name=raw_sigs.name, sets={f"SIG_{k}": v for k, v in raw_sigs.sets.items()})
    filtered = sigs.filter_to(set(expr.columns))
    missing = [c for c in sig_cols if c not in filtered.sets]
    if missing:
        raise SystemExit(f"no filtered gene set for {missing} -- the null sets "
                         "are size- and expression-matched to the real "
                         "signature, so its membership is required")
    gene_means = expr.mean(axis=0)

    print("\n  S0: standalone bootstrap probe on synthetic input", flush=True)
    s0 = _stage0_bootstrap_synthetic()
    print(f"     {s0}", flush=True)

    out: dict[str, list] = {k: [] for k, _ in STAGES if k != "S0_boot_synth"}
    print(f"\n  {'signature':<34}{'r_obs':>10}{'null_r_mean':>13}"
          f"{'z_excess':>11}{'secs':>8}", flush=True)

    for col in sig_cols:
        t_sig = time.time()
        # Null sets, matched on size and mean expression, exactly as
        # `experiment._residualised_null_scores` built them at `8fe1187` (lines
        # 1068-1078).
        single = sig_mod.SignatureSet(name=col, sets={col: filtered.sets[col]})
        family = sig_mod.random_gene_sets(
            sorted(expr.columns), single, n_per_signature=n_null,
            seed=cfg.seed, match_expression=gene_means,
        )[col]
        scores = sig_mod.score_mean_z(expr, family)
        scores = scores.reindex(frame["patient_id"].to_numpy())
        scores.index = frame.index

        obs_raw = pd.to_numeric(frame[col], errors="coerce").to_numpy(dtype=float)
        stacked_scores = np.column_stack([obs_raw, scores.to_numpy(dtype=float)])
        stacked_resid = globalaxis.residualise_matrix(
            stacked_scores, axis_values, groups=within)
        obs_resid, null_resid = stacked_resid[:, 0], stacked_resid[:, 1:]

        # FAIL FAST. `_select_alpha_per_fold` returns NaN for any fold with <20
        # finite targets, silently; an all-NaN y once produced 80 NaN alphas.
        finite_frac = float(np.isfinite(obs_resid).mean())
        if finite_frac < 0.5:
            raise SystemExit(f"{col}: only {100 * finite_frac:.1f}% of "
                             "residualised scores are finite")

        fold_alpha = experiment._select_alpha_per_fold(X, obs_resid, split)
        if not np.all(np.isfinite(np.asarray(fold_alpha, dtype=float))):
            raise SystemExit(f"{col}: non-finite fold alpha {fold_alpha}")

        stacked = np.column_stack([obs_resid, null_resid])
        P = models.cross_val_predict_multi(X, stacked, split, patients,
                                           fixed_alpha=fold_alpha)
        pred_obs, pred_null = P[:, 0], P[:, 1:]

        r_obs = stats.corr_ci(obs_resid, pred_obs).value
        r_null = np.asarray([stats.corr_ci(null_resid[:, j], pred_null[:, j]).value
                             for j in range(null_resid.shape[1])], dtype=float)
        # The per-signature ISI term, UNDISATTENUATED. Disattenuation divides by
        # sqrt(reliability), which is a separate multiplicative stage; keeping it
        # out means a divergence here is unambiguously the ridge/corr chain.
        z_excess = float(stats.fisher_z(r_obs)
                         - np.nanmean(stats.fisher_z(r_null)))

        out["S2_obs_resid"].append(obs_resid)
        out["S2_null_resid"].append(null_resid)
        out["S3_fold_alpha"].append(np.asarray(fold_alpha, dtype=float))
        out["S4_pred_obs"].append(pred_obs)
        out["S4_pred_null"].append(pred_null)
        out["S5_r_obs"].append(np.array([r_obs]))
        out["S5_r_null"].append(r_null)
        out["S6_z_excess"].append(np.array([z_excess]))

        print(f"  {col[:32]:<34}{r_obs:>10.6f}{np.nanmean(r_null):>13.6f}"
              f"{z_excess:>11.6f}{time.time() - t_sig:>8.1f}", flush=True)

    payload = {"S0_boot_synth": s0}
    for key in out:
        payload[key] = np.stack([np.atleast_2d(a) if a.ndim > 1 else a
                                 for a in out[key]]) if out[key] else np.array([])
    payload["signatures"] = np.array(sig_cols, dtype=object).astype(str)
    payload["pooled_z_excess"] = np.array(
        [float(np.mean(np.concatenate(out["S6_z_excess"])))])

    dest = RESULTS / f"bisect_{tag}.npz"
    np.savez(dest, **payload)
    print(f"\n  pooled undisattenuated z-excess over {len(sig_cols)} "
          f"signatures: {payload['pooled_z_excess'][0]:.12f}")
    print(f"  total {time.time() - t0:.1f}s")
    print(f"\nwrote {dest.relative_to(ROOT)}")
    print("Run the SAME command on the other machine, then:")
    print("  python scripts/18_bisect_platform.py --compare "
          "results/bisect_macos.npz results/bisect_hpc4.npz")
    return 0


def compare(path_a: Path, path_b: Path) -> int:
    a, b = np.load(path_a, allow_pickle=True), np.load(path_b, allow_pickle=True)
    print(f"A  {path_a}")
    print(f"B  {path_b}\n")
    if "signatures" in a and "signatures" in b:
        if not np.array_equal(a["signatures"], b["signatures"]):
            print("!! the two dumps cover different signatures -- not comparable")
            return 1

    print(f"  {'stage':<18}{'max |A-B|':>14}{'max rel':>13}   description")
    first_diverged = None
    for key, desc in STAGES:
        if key not in a.files or key not in b.files:
            print(f"  {key:<18}{'MISSING':>14}{'':>13}   {desc}")
            continue
        xa, xb = np.asarray(a[key], float), np.asarray(b[key], float)
        if xa.shape != xb.shape:
            print(f"  {key:<18}{'SHAPE':>14}{'':>13}   {xa.shape} vs {xb.shape}")
            continue
        ok = np.isfinite(xa) & np.isfinite(xb)
        if not ok.any():
            print(f"  {key:<18}{'all-NaN':>14}{'':>13}   {desc}")
            continue
        d = np.abs(xa[ok] - xb[ok])
        scale = np.maximum(np.abs(xa[ok]), 1e-300)
        amax, rmax = float(d.max()), float((d / scale).max())
        flag = "" if amax <= FLOOR else "  <-- DIVERGES"
        if amax > FLOOR and first_diverged is None:
            first_diverged = (key, desc, amax)
        print(f"  {key:<18}{amax:>14.3e}{rmax:>13.3e}   {desc}{flag}")

    for key in ("pooled_z_excess",):
        if key in a.files and key in b.files:
            print(f"\n  {key}: A={float(a[key][0]):.12f}  B={float(b[key][0]):.12f}  "
                  f"delta={float(a[key][0]) - float(b[key][0]):+.12f}")

    print()
    if first_diverged is None:
        print(f"VERDICT: every stage agrees to within {FLOOR:.0e}. The gap does "
              "NOT enter\n  anywhere measured here. Widen --n-null, or the "
              "divergence is in the\n  reliability/disattenuation stage this "
              "script deliberately excludes.")
    else:
        key, desc, amax = first_diverged
        print(f"VERDICT: the divergence FIRST appears at {key} "
              f"(max |A-B| = {amax:.3e}),\n  which is: {desc}.")
        _amplification_check(a, b)
    return 0


def _amplification_check(a, b) -> None:
    """Test -- not assert -- that later stages are amplification of S3.

    This function exists because the sentence it replaces was wrong. Until
    2026-09-06 this script closed with "Every later stage is downstream
    amplification of this one, not a separate defect. Attack this stage." That
    is an inference, and the data needed to check it were already in the two
    dumps nobody had compared.

    The check is a natural control that costs nothing: some signatures have
    IDENTICAL fold alphas on both machines. If S3 were the only defect, those
    signatures' downstream stages would agree to round-off. If they diverge
    anyway, there is a second, independent divergence and "attack S3" sends the
    next session to the wrong stage.
    """
    need = ("S3_fold_alpha", "S4_pred_obs", "S4_pred_null", "S5_r_null")
    if not all(k in a.files and k in b.files for k in need):
        print("  (amplification check skipped: dumps predate S3/S4/S5 keys)")
        return

    aa, ba = np.asarray(a["S3_fold_alpha"], float), np.asarray(b["S3_fold_alpha"], float)
    if aa.shape != ba.shape:
        print("  (amplification check skipped: S3 shapes differ)")
        return

    same = np.where((aa == ba).all(axis=1))[0]
    n_sig = aa.shape[0]
    n_flip = int((aa != ba).sum())
    print(f"\n  AMPLIFICATION CHECK: {n_flip} of {aa.size} fold alphas differ; "
          f"{len(same)} of {n_sig}\n  signature(s) have IDENTICAL alphas on both "
          "machines and act as a control.")
    if len(same) == 0:
        print("  No control signatures -- every signature flipped, so this run "
              "cannot\n  separate 'amplification of S3' from a second defect.")
        return

    po = np.abs(np.asarray(a["S4_pred_obs"], float)[same]
                - np.asarray(b["S4_pred_obs"], float)[same]).max()
    pn = np.abs(np.asarray(a["S4_pred_null"], float)[same]
                - np.asarray(b["S4_pred_null"], float)[same]).max()
    worst = max(po, pn)
    if worst <= FLOOR:
        print(f"  Control signatures agree downstream to {worst:.3e}. The "
              "amplification\n  reading HOLDS: S3 really is the single origin. "
              "Attack that stage.")
        return

    scale = float(np.asarray(a["S4_pred_obs"], float)[same].std())
    print(f"  Control signatures DIVERGE downstream anyway: max |A-B| in the "
          f"out-of-fold\n  predictions is {worst:.3e}, which is "
          f"{100 * worst / max(scale, 1e-300):.0f}% of their own sd ({scale:.4f}).")
    print("  So S3 is the FIRST divergence but NOT the only one: the ridge solve "
          "differs\n  at identical penalties on identical inputs. Attack "
          "models.cross_val_predict_multi,\n  not the alpha grid.")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--tag", default="local",
                   help="machine tag; writes results/bisect_<tag>.npz")
    p.add_argument("--n-null", type=int, default=8,
                   help="null draws per signature. A platform divergence shows "
                        "at 8 as clearly as at 1000 (default: 8)")
    p.add_argument("--signatures", type=int, default=None,
                   help="bisect on the first N signatures only; one is enough "
                        "to localise a stage")
    p.add_argument("--compare", nargs=2, metavar=("A.npz", "B.npz"),
                   help="diff two dumps and name the first diverging stage")
    args = p.parse_args()
    if args.compare:
        return compare(Path(args.compare[0]), Path(args.compare[1]))
    return run(args.tag, args.n_null, args.signatures)


if __name__ == "__main__":
    raise SystemExit(main())
