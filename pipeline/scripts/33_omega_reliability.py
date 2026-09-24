#!/usr/bin/env python3
"""McDonald's omega beside Cronbach's alpha, and the NSCLC index under omega (ledger E8).

    python scripts/33_omega_reliability.py      # writes results/omega_reliability/

WHY
===
The index divides each correlation by the square root of its target's
reliability, and the registered reliability is Cronbach's alpha of the
axis-residualised composite (`signatures.alpha_from_score_and_item_variances`).
Alpha assumes every gene loads equally on one factor (tau-equivalence). Where
that fails, alpha understates reliability, and it does so by different amounts
for curated and random sets. McDonald's omega relaxes the assumption.

    omega = (sum_j lambda_j)^2 / ((sum_j lambda_j)^2 + sum_j psi_j)

lambda and psi are the loadings and unique variances of a one-factor model
(sklearn FactorAnalysis, LAPACK SVD, so no randomness), fitted to the SAME
residualised z-scored items the registered alpha uses. The sum is signed, as
the composite is the unweighted mean.

WHAT IT COMPUTES (NSCLC, from stored artefacts plus expression)
===============================================================
* alpha and omega for the 16 curated sets, and for the first 25 null draws of
  each (the same draws the audit uses: same pool, seed 0, expression-matched).
  Gate: this alpha for the curated sets equals the registered `alpha_observed`
  in results/nsclc_v3_stablesort/immune_excess.csv to 1e-10.
* The disattenuated excess reconstructed with omega, and, as the like-for-like
  control, with alpha. Both use the stored null correlations
  (null_draws.npz), the stored refitted r, the mean null reliability of the 25
  draws, and `globalaxis.excess_over_null`'s analytic interval. Pooled across
  signatures with the across-signature standard error, as the split-half
  reconstruction does (script 13). It is not the registered bootstrap.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import FactorAnalysis

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import globalaxis, signatures  # noqa: E402

_spec = importlib.util.spec_from_file_location("sens13", ROOT / "scripts" / "13_scorer_sensitivity.py")
m13 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m13)

RUN = ROOT / "results" / "nsclc_v3_stablesort"
CRIT = 1.959963984540054


def omega(items: np.ndarray) -> float:
    """McDonald's omega-total of the unweighted composite of `items` (n x k)."""
    if items.shape[1] < 2:
        return float("nan")
    fa = FactorAnalysis(n_components=1, svd_method="lapack").fit(items)
    lam = fa.components_[0]
    common = float(lam.sum()) ** 2
    return common / (common + float(fa.noise_variance_.sum()))


def pooled(vals: np.ndarray) -> dict:
    vals = vals[np.isfinite(vals)]
    se = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else float("nan")
    m = float(np.mean(vals)) if len(vals) else float("nan")
    return {"isi": m, "lo": m - CRIT * se, "hi": m + CRIT * se, "n_signatures": int(len(vals))}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-null-reliability", type=int, default=25)
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "omega_reliability")
    args = ap.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    cohort, _, expr, sigs = m13.load_nsclc()
    filtered = sigs.filter_to(set(expr.columns))
    axis = globalaxis.compute_global_axis(
        expr, method="pc1", cancer_type=cohort.set_index("patient_id")["cancer_type"])
    axis_values = pd.Series(axis.values).reindex(expr.index)
    within = pd.Series(cohort["cancer_type"].to_numpy(), index=expr.index)
    # The registered items: z-scored over the expression's own samples, then
    # residualised on the axis within type (experiment.run_audit).
    z = (expr - expr.mean(axis=0)) / expr.std(axis=0).replace(0, np.nan)
    z = z.dropna(axis=1, how="all").fillna(0.0)
    R = globalaxis.residualise_matrix(
        z.to_numpy(dtype=float), axis_values.to_numpy(dtype=float),
        groups=within.to_numpy())
    col = {g: j for j, g in enumerate(z.columns)}
    var = np.nanvar(R, axis=0, ddof=1)

    def reliabilities(genes):
        idx = [col[g] for g in genes if g in col]
        items = R[:, idx]
        a = signatures.alpha_from_score_and_item_variances(items.mean(axis=1), var[idx])
        return float(a), omega(items)

    ie = pd.read_csv(RUN / "immune_excess.csv", float_precision="round_trip").set_index("signature")
    draws = np.load(RUN / "null_draws.npz")
    rows, worst = [], 0.0
    for name, genes in filtered.sets.items():
        a_obs, w_obs = reliabilities(genes)
        worst = max(worst, abs(a_obs - float(ie.loc[name, "alpha_observed"])))
        single = signatures.SignatureSet(name=name, sets={name: genes})
        fam = signatures.random_gene_sets(
            sorted(expr.columns), single, n_per_signature=1000, seed=0,
            match_expression=expr.mean(axis=0))[name]
        null = [reliabilities(g) for g in list(fam.sets.values())[: args.n_null_reliability]]
        a_null = float(np.nanmean([x[0] for x in null]))
        w_null = float(np.nanmean([x[1] for x in null]))
        r_obs, n = float(ie.loc[name, "r_residual_refit"]), int(ie.loc[name, "n"])
        nulls = np.asarray(draws[name], dtype=float)
        e_a = globalaxis.excess_over_null(r_obs, nulls, n=n, reliability_observed=a_obs,
                                          reliability_null=a_null)
        e_w = globalaxis.excess_over_null(r_obs, nulls, n=n, reliability_observed=w_obs,
                                          reliability_null=w_null)
        rows.append({"signature": name, "k": len(genes),
                     "alpha_observed": a_obs, "omega_observed": w_obs,
                     "alpha_null_mean": a_null, "omega_null_mean": w_null,
                     "excess_alpha": e_a.value, "excess_alpha_lo": e_a.lo,
                     "excess_omega": e_w.value, "excess_omega_lo": e_w.lo,
                     "excess_omega_hi": e_w.hi,
                     "beats_null_omega": bool(np.isfinite(e_w.lo) and e_w.lo > 0)})
        print(f"  {name.replace('SIG_HALLMARK_', ''):38s} alpha {a_obs:.3f}/{a_null:.3f} "
              f"omega {w_obs:.3f}/{w_null:.3f}  excess alpha {e_a.value:+.3f} "
              f"omega {e_w.value:+.3f}", flush=True)
    if worst > 1e-10:
        raise SystemExit(f"GATE FAILED: alpha differs from the registered alpha_observed by {worst:.3g}")
    d = pd.DataFrame(rows)
    d.to_csv(args.outdir / "per_signature.csv", index=False)
    summary = {
        "cohort": "nsclc", "source_run": str(RUN.relative_to(ROOT)),
        "n_null_reliability": args.n_null_reliability,
        "alpha_gate_max_abs_diff": worst,
        "median_alpha_observed": float(d["alpha_observed"].median()),
        "median_omega_observed": float(d["omega_observed"].median()),
        "median_alpha_null": float(d["alpha_null_mean"].median()),
        "median_omega_null": float(d["omega_null_mean"].median()),
        "median_gap_alpha": float((d["alpha_observed"] - d["alpha_null_mean"]).median()),
        "median_gap_omega": float((d["omega_observed"] - d["omega_null_mean"]).median()),
        "reconstructed_alpha": pooled(d["excess_alpha"].to_numpy(dtype=float)),
        "reconstructed_omega": pooled(d["excess_omega"].to_numpy(dtype=float)),
        "beats_null_omega": int(d["beats_null_omega"].sum()),
        "registered_isi": json.loads((RUN / "summary.json").read_text())["primary_excess"],
        "method": "excess_over_null with mean null reliability of the first draws; across-"
                  "signature SE; NOT the registered patient-clustered bootstrap",
        "seconds": round(time.time() - t0, 1),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "registered_isi"}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
