#!/usr/bin/env python3
"""A third scorer for the scorer-sensitivity comparison: PLAGE (ledger E6).

    python scripts/32_third_scorer.py                 # NSCLC, n_null 1000; ~30 min on 8 CPUs
    python scripts/32_third_scorer.py --n-null 200    # a quick look; NOT the main analysis

WHY
===
Limitation 9's scorer-dependence rests on two scorers. Under mean-z the index is
positive with 16/16 signatures. Under ssGSEA the registered estimand is
undefined, because its reliability estimator assumes a mean composite, and the
split-half reconstruction reverses sign (`results/scorer_sensitivity_local/`,
A7; pinned-ordering ssGSEA arm in `results/scorer_sensitivity_ssgsea_stable/`).
Two points cannot say whether that is a mean-versus-rank split or something
about ssGSEA alone. PLAGE is a third kind: a data-driven weighted composite,
the first singular vector of the gene-standardised set (`signatures.score_plage`).

WHAT IT RUNS (NSCLC, current code, pinned orderings)
====================================================
It uses `13_scorer_sensitivity.py`'s own functions, imported unchanged: the
rescoring, the arm runner, the batched split-half reliability and its self-check.
  mean_z__uncorrected   the like-for-like control, under the same code as PLAGE
  plage__uncorrected    the scorer comparison proper
  plage__disattenuated  the registered estimand under PLAGE, as configured. Its
                        reliability estimator assumes a mean composite, so this
                        is expected to be undefined, as it is under ssGSEA.
  reconstruction        the disattenuated excess with split-half reliability,
                        for mean_z and plage, exactly as script 13 builds it.
The ssGSEA values are read from the stored pinned-ordering run, not recomputed.

Writes results/scorer_sensitivity_plage/{summary.json,per_signature.csv}.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import globalaxis, signatures  # noqa: E402

_spec = importlib.util.spec_from_file_location("sens13", ROOT / "scripts" / "13_scorer_sensitivity.py")
m13 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m13)

ARMS = (("mean_z", False), ("plage", False), ("plage", True))
SSGSEA_STABLE = ROOT / "results" / "scorer_sensitivity_ssgsea_stable" / "summary.json"


def reconstruct(table, rel):
    """Script 13's reconstruction, for one scorer (same function, same inputs)."""
    draws = table.attrs.get("null_draws", {})
    rows = []
    for _, row in table.iterrows():
        name = row["signature"]
        nulls = np.asarray(draws.get(name, []), dtype=float)
        rel_o, rel_n = rel["observed"].get(name, np.nan), rel["null_mean"].get(name, np.nan)
        if not len(nulls) or not np.isfinite(rel_o) or not np.isfinite(rel_n):
            continue
        est = globalaxis.excess_over_null(
            float(row["r_residual_refit"]), nulls, n=int(row["n"]),
            reliability_observed=rel_o, reliability_null=rel_n)
        rows.append({"signature": name, "excess_z": est.value, "excess_lo": est.lo,
                     "excess_hi": est.hi, "beats_null": bool(np.isfinite(est.lo) and est.lo > 0),
                     "reliability_observed": rel_o, "reliability_null": rel_n})
    d = pd.DataFrame(rows)
    vals = d["excess_z"].to_numpy(dtype=float) if len(d) else np.array([])
    vals = vals[np.isfinite(vals)]
    se = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else np.nan
    return {
        "isi": float(np.mean(vals)) if len(vals) else None,
        "lo": float(np.mean(vals) - m13.CRIT * se) if len(vals) > 1 else None,
        "hi": float(np.mean(vals) + m13.CRIT * se) if len(vals) > 1 else None,
        "n_signatures": int(len(vals)),
        "beats_null": int(d["beats_null"].sum()) if len(d) else 0,
        "method": "excess-over-null-z-disattenuated-with-split-half-reliability; CI is "
                  "the across-signature SE, NOT the registered patient-clustered bootstrap",
    }, d


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-null", type=int, default=1000)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--n-splits", type=int, default=20)
    ap.add_argument("--n-null-reliability", type=int, default=25)
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "scorer_sensitivity_plage")
    args = ap.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    t_all = time.time()

    cohort, X, expr, sigs = m13.load_nsclc()
    filtered = sigs.filter_to(set(expr.columns))
    axis = globalaxis.compute_global_axis(
        expr, method="pc1", cancer_type=cohort.set_index("patient_id")["cancer_type"])
    axis_values = pd.Series(axis.values).reindex(expr.index)
    axis_values.index = expr.index
    within = pd.Series(cohort["cancer_type"].to_numpy(), index=expr.index)
    ref = m13.selfcheck_split_half(expr, filtered.sets[sorted(filtered.sets)[0]], axis_values, within)
    print(f"nsclc: {len(cohort)} patients; split-half self-check passed ({ref:.4f}); "
          f"n_null {args.n_null}", flush=True)

    arms, tables, rel, recon, pieces = {}, {}, {}, {}, []
    for scorer, dis in ARMS:
        key = f"{scorer}__{'disattenuated' if dis else 'uncorrected'}"
        frame, sig_cols = m13.rescore_cohort(cohort, expr, sigs, scorer)
        res, dt = m13.run_arm(frame, X, sig_cols, expr, sigs, axis,
                              scorer=scorer, disattenuate=dis, args=args)
        pooled, table = m13.summarise(res)
        pooled["seconds"] = round(dt, 1)
        pooled["beats_null"] = int(table["beats_null"].sum()) if len(table) else 0
        arms[key], tables[key] = pooled, table
        print(f"  {key:22s} {m13.fmt(pooled)}  beats {pooled['beats_null']}  ({dt:.0f}s)", flush=True)
        if len(table):
            have = [c for c in ("excess_z", "excess_lo", "excess_hi", "beats_null",
                                "r_residual_refit", "null_mean_r") if c in table.columns]
            pieces.append(table[["signature"] + have].rename(
                columns={c: f"{c}__{key}" for c in have}))

    for scorer in ("mean_z", "plage"):
        t0 = time.time()
        fn = getattr(signatures, f"score_{scorer}")
        obs = m13.split_half_reliabilities(expr, filtered.sets, fn, axis_values, within,
                                           n_splits=args.n_splits, seed=0)
        sampled = {}
        for name, genes in filtered.sets.items():
            single = signatures.SignatureSet(name=name, sets={name: genes})
            fam = signatures.random_gene_sets(
                sorted(expr.columns), single, n_per_signature=args.n_null, seed=0,
                match_expression=expr.mean(axis=0))[name]
            sampled.update(dict(list(fam.sets.items())[: args.n_null_reliability]))
        per_draw = m13.split_half_reliabilities(expr, sampled, fn, axis_values, within,
                                                n_splits=args.n_splits, seed=0)
        null = {name: float(np.nanmean([v for k, v in per_draw.items()
                                        if k.startswith(f"{name}__random_")] or [np.nan]))
                for name in filtered.sets}
        rel[scorer] = {"observed": {k: float(v) for k, v in obs.items()}, "null_mean": null}
        recon[scorer], d = reconstruct(tables[f"{scorer}__uncorrected"], rel[scorer])
        if len(d):
            pieces.append(d.rename(columns={c: f"{c}__{scorer}__reconstructed"
                                            for c in d.columns if c != "signature"}))
        print(f"  {scorer} reconstructed  isi {recon[scorer]['isi']}  beats "
              f"{recon[scorer]['beats_null']}  reliability observed "
              f"{np.nanmean(list(obs.values())):.4f} null {np.nanmean(list(null.values())):.4f}"
              f"  ({time.time() - t0:.0f}s)", flush=True)

    merged = pieces[0]
    for piece in pieces[1:]:
        merged = merged.merge(piece, on="signature", how="outer")
    merged.to_csv(args.outdir / "per_signature.csv", index=False)
    ss = json.loads(SSGSEA_STABLE.read_text()) if SSGSEA_STABLE.exists() else {}
    summary = {
        "cohort": "nsclc", "n_patients": int(len(cohort)), "n_null_sets": args.n_null,
        "n_null_sets_matches_main_analysis": args.n_null == 1000, "n_boot": args.n_boot,
        "arms": arms, "reconstructed": recon,
        "split_half_reliability": {s: {"observed_mean": float(np.nanmean(list(r["observed"].values()))),
                                       "null_mean": float(np.nanmean(list(r["null_mean"].values())))}
                                   for s, r in rel.items()},
        "ssgsea_pinned_from": str(SSGSEA_STABLE.relative_to(ROOT)),
        "ssgsea_pinned": {"uncorrected": (ss.get("arms") or {}).get("ssgsea__uncorrected"),
                          "reconstructed": (ss.get("reconstructed") or {}).get("ssgsea")},
        "runtime_s": round(time.time() - t_all, 1),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"wrote {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
