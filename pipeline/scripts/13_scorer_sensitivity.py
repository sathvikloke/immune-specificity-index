#!/usr/bin/env python3
"""Does the ISI headline survive a RANK-BASED scorer? (science-audit item A7)

    python scripts/13_scorer_sensitivity.py --n-null 1000

    # just the ssGSEA question; skips the two mean_z arms (3.2 h on a loaded Mac)
    python scripts/13_scorer_sensitivity.py --arms ssgsea__uncorrected \
                                                   ssgsea__disattenuated

WHY THIS HAD NEVER BEEN RUN
===========================
`AuditConfig.scorer` has always accepted "ssgsea", and every number in this
repository was produced with "mean_z". The reason was arithmetic, not judgement:
`signatures.score_ssgsea` cost 6.26 ms per patient-set, and one NSCLC ISI scores
944 patients x 16 signatures x (1 + n_null) gene sets, so the sensitivity was a
24-hour job. The scorer now costs 0.0127 ms per patient-set (494x; the per-sample
rank order and |rank|^alpha weights were being recomputed once per gene set even
though neither depends on the gene set), which puts the whole comparison inside
half an hour and makes n_null_sets=1000 -- the same null as the main analysis --
affordable rather than a compromise.

WHAT THIS COMPARISON CAN AND CANNOT SAY
=======================================
The pre-registered ISI disattenuates for reliability, and the audit estimates
that reliability with `signatures.alpha_from_score_and_item_variances`: the
general-form Cronbach alpha for a composite that is the MEAN of its k z-scored
genes. A mean-z score is exactly that composite. An ssGSEA score is NOT -- it is
a rank-walk statistic over the whole transcriptome, whose SD on NSCLC is ~20x
smaller than the mean-z score's. Feeding it to that formula compares two
unrelated scales:

    HALLMARK_INTERFERON_GAMMA_RESPONSE, k=199, residualised, NSCLC
      mean_z   sum(item var) 170.19   k^2 var(score) 8512.49   alpha 0.985
      ssGSEA   sum(item var) 170.19   k^2 var(score)   18.20   alpha 0.000

so `globalaxis.disattenuate` returns NaN for every draw and the registered
estimand is undefined under ssGSEA. That is a property of the reliability
ESTIMATOR, not of ssGSEA, and reporting it as "ssGSEA kills the ISI" would be
wrong. This script therefore reports three things, in decreasing authority:

  1. UNCORRECTED excess under both scorers. Like for like, the correction held
     absent on both sides, so any difference is attributable to the scorer. This
     is the arm that answers "does the scorer change the answer". It is NOT the
     primary endpoint -- the uncorrected contrast had a 100% false-positive rate
     in simulation (see `globalaxis.excess_over_null`) -- so it cannot be quoted
     as an ISI, only as a scorer comparison.
  2. The registered disattenuated ISI under each scorer, exactly as configured,
     including its failure under ssGSEA.
  3. A RECONSTRUCTION of the disattenuated ISI under each scorer using
     split-half reliability, which is scorer-agnostic and is what
     `signatures.split_half_reliability`'s own docstring says to use for ssGSEA.
     Its interval is the analytic one from `excess_over_null`, not the registered
     patient-clustered bootstrap, so it is a secondary estimate.

n_null_sets
===========
Chosen by measurement, not by guess. The script times one real ssGSEA scoring
call before starting and prints the extrapolation. Measured on this machine:
a 1000-draw family over 944 patients scores in 12.0 s, and a full audit at
n_null=1000 runs in about 8 minutes per scorer per arm -- four arms, well inside
two hours -- so the default is 1000 and the null is IDENTICAL to the main
analysis. If you lower it, the printed banner says so and the comparison stops
being equivalent to the headline: at B draws the smallest attainable one-sided
permutation p is 1/(B+1), and BH-FDR over 16 signatures needs p <= 0.05/16 =
0.0031, so anything below B=319 cannot reject at all.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import experiment, globalaxis, signatures  # noqa: E402

INTERIM = ROOT / "data" / "interim"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"

SCORERS = ("mean_z", "ssgsea")
# The four arms, as `<scorer>__<correction>` keys. `--arms` selects a subset.
# The mean_z arms cost 3.2 h on a loaded 14-core Mac and their answer is already
# known (+0.3182, 16/16 — the control that exonerated this harness), so a rerun
# aimed at the ssGSEA question should not be forced to recompute them.
ARMS = tuple(f"{s}__{c}" for s in SCORERS for c in ("uncorrected", "disattenuated"))
CRIT = 1.959963984540054


def plan_arms(requested) -> dict:
    """Resolve `--arms` into the work plan. PURE: no I/O, no cohort, no compute.

    Extracted from `main()` on 2026-09-05 because the partial-run path was the
    oldest untested behaviour in the project. The test that existed said so in
    its own docstring -- "reaching that code requires the real NSCLC cohort and
    hours of compute. It has been smoke-tested by hand" -- which is another way
    of saying the flag that exists to make a partial run affordable could only
    be checked by paying for a full one. Everything the flag actually decides is
    a set operation over `ARMS`, so none of it needs the cohort.

    Returns the canonical selection plus the three derived questions `main()`
    asks: which scorers must run at all, which can be reconstructed, and whether
    the cross-scorer agreement section is computable. `selected` is always in
    ARMS order, never the user's, so a user's argument order cannot leak into
    `arms_requested` in summary.json.
    """
    wanted = set(requested)
    unknown = sorted(wanted - set(ARMS))
    if unknown:
        raise ValueError(f"unknown arm(s): {unknown}; choices are {list(ARMS)}")
    selected = [a for a in ARMS if a in wanted]
    # A scorer runs if any of its arms was selected. Reconstruction and
    # split-half reliability additionally need that scorer's UNCORRECTED arm,
    # because they are built from the null draws that arm stores in
    # `table.attrs["null_draws"]`; without it there is nothing to reconstruct.
    recon = [s for s in SCORERS if f"{s}__uncorrected" in selected]
    return {
        "selected": selected,
        "skipped": [a for a in ARMS if a not in wanted],
        "scorers_run": [s for s in SCORERS
                        if any(a.startswith(f"{s}__") for a in selected)],
        "recon_scorers": recon,
        "complete": selected == list(ARMS),
        "agreement_possible": len(recon) >= 2,
    }


def arm_status_rows(plan: dict, results) -> list[dict]:
    """One row per arm, in canonical order, each with an explicit status. PURE.

    The point is that a missing arm is never a `KeyError` and never silently
    absent: an arm the user EXCLUDED and an arm that was requested and then
    failed to produce a result are different facts, and a reader of a partial
    summary.json is entitled to tell them apart. Iterating `ARMS` rather than
    `results` is what guarantees every arm appears.
    """
    rows = []
    for arm in ARMS:
        if arm not in plan["selected"]:
            status, value = "NOT RUN (excluded by --arms)", None
        elif arm not in results:
            status, value = "REQUESTED BUT MISSING", None
        else:
            status, value = "ok", results[arm]
        rows.append({"arm": arm, "status": status, "value": value})
    return rows


def load_nsclc():
    cohort = pd.read_parquet(INTERIM / "cohort_nsclc.parquet")
    X = np.load(INTERIM / "X_nsclc.npy")
    expr = pd.read_parquet(INTERIM / "expr_nsclc.parquet")
    if not expr.index.equals(pd.Index(cohort["patient_id"])):
        expr = expr.reindex(cohort["patient_id"])
    sigs = signatures.SignatureSet.from_gmt(GMT)
    sigs = signatures.SignatureSet(
        name=sigs.name, sets={f"SIG_{k}": v for k, v in sigs.sets.items()}, source=str(GMT))
    return cohort, X, expr, sigs


# ------------------------------------------------------------- scoring ---

def rescore_cohort(cohort, expr, sigs, scorer):
    """`cohort` with its SIG_ columns REPLACED by scores from `scorer`.

    The cohort parquet ships mean-z columns. The audit scores the NULL with
    `AuditConfig.scorer`, so leaving those columns in place while switching the
    null to ssGSEA would put a mean-z observation against an ssGSEA null -- the
    exact scale mismatch `experiment._score_sets` warns about, and the two
    scorers' SDs differ by ~20x here. Rescoring both sides is not optional.

    `_score_sets` is reached into deliberately: it IS the dispatch the null side
    uses, so borrowing it is what guarantees the two sides cannot drift apart.
    """
    filtered = sigs.filter_to(set(expr.columns))
    scored = experiment._score_sets(
        expr, filtered, experiment.AuditConfig(scorer=scorer))
    frame = cohort.copy()
    scored = scored.reindex(frame["patient_id"].to_numpy())
    scored.index = frame.index
    for c in scored.columns:
        frame[c] = scored[c].to_numpy()
    return frame, list(scored.columns)


# ------------------------------------------------ scorer-agnostic reliability ---

def split_half_reliabilities(
    expr, gene_lists, scorer, axis_values, within, *, n_splits=20, seed=0
):
    """Reliability of each gene set AS SCORED, every half in ONE scorer call.

    Same estimator as `signatures.split_half_reliability` -- same RNG stream,
    same halving rule, same Spearman-Brown step-up r_full = 2r/(1+r), same
    clipping -- and `selfcheck_split_half` below asserts the two agree. What
    differs is that all 2*n_splits halves of every gene set are handed to the
    scorer at once. That is the only affordable shape for ssGSEA, which pays a
    fixed ~9 s per call to rank the whole 41k-gene matrix regardless of how many
    sets it is asked for; called per split per gene set, the fixed cost is 50x
    the scoring it was preparing for.

    Split-half is used here because Cronbach's alpha -- in either of its closed
    forms in `signatures` -- assumes the score is the mean of its k z-scored
    genes. That is true of mean-z and false of ssGSEA, where it returns 0 and
    takes the whole disattenuated estimand with it.
    """
    halves: dict[str, list[str]] = {}
    plan: dict[str, list[tuple[str, str]]] = {}
    for name, genes in gene_lists.items():
        cols = [g for g in genes if g in expr.columns]
        if len(cols) < 4:
            continue
        rng = np.random.default_rng(seed)
        pairs = []
        for s in range(n_splits):
            perm = rng.permutation(cols)
            a, b = list(perm[: len(perm) // 2]), list(perm[len(perm) // 2:])
            ka, kb = f"{name}\x00{s}\x00a", f"{name}\x00{s}\x00b"
            halves[ka], halves[kb] = a, b
            pairs.append((ka, kb))
        plan[name] = pairs
    if not plan:
        return {}

    scored = scorer(expr, signatures.SignatureSet(name="halves", sets=halves))
    resid = globalaxis.residualise_matrix(
        scored.to_numpy(dtype=float),
        pd.to_numeric(axis_values, errors="coerce").to_numpy(dtype=float),
        groups=None if within is None else np.asarray(within),
    )
    at = {c: j for j, c in enumerate(scored.columns)}

    out = {}
    for name, pairs in plan.items():
        rs = []
        for ka, kb in pairs:
            if ka not in at or kb not in at:
                continue
            sa, sb = resid[:, at[ka]], resid[:, at[kb]]
            ok = np.isfinite(sa) & np.isfinite(sb)
            if ok.sum() < 10:
                continue
            r = float(np.corrcoef(sa[ok], sb[ok])[0, 1])
            if np.isfinite(r):
                rs.append(r)
        if not rs:
            out[name] = float("nan")
            continue
        r_half = float(np.median(rs))
        r_full = 2 * r_half / (1 + r_half) if r_half > -1 else float("nan")
        out[name] = float(np.clip(r_full, 0.0, 1.0))
    return out


def selfcheck_split_half(expr, genes, axis_values, within, *, n_splits=6, seed=0):
    """Assert the batched split-half equals the library's one-set-at-a-time form.

    Run on mean-z and on ONE gene set, because that is cheap and the batching is
    scorer-independent. If this drifts, every reconstructed ISI below is wrong,
    so it fails loudly rather than warning.
    """
    def residualiser(s):
        return globalaxis.residualise(s, axis_values, within=within)

    reference = signatures.split_half_reliability(
        expr, genes, scorer=signatures.score_mean_z, residualiser=residualiser,
        n_splits=n_splits, seed=seed)["reliability"]
    batched = split_half_reliabilities(
        expr, {"x": genes}, signatures.score_mean_z, axis_values, within,
        n_splits=n_splits, seed=seed)["x"]
    if not np.isclose(reference, batched, atol=1e-10):
        raise AssertionError(
            f"batched split-half {batched:.6f} != library "
            f"split_half_reliability {reference:.6f}; the reconstruction below "
            "would be estimating a different quantity"
        )
    return float(batched)


# ----------------------------------------------------------------- runs ---

def run_arm(frame, X, sig_cols, expr, sigs, axis, *, scorer, disattenuate, args):
    # `run_combat=False` and `run_legacy_venet_null=False` differ from the frozen
    # nsclc_v3 config, which had both True. NEITHER TOUCHES THE ISI, so the
    # mean_z arms remain a valid control for the stored +0.3182:
    #   * run_combat gates only `site_prediction_control_oof_combat`, a negative
    #     control on site predictability, plus the note it appends. Nothing it
    #     produces feeds the residualised null or the disattenuated excess.
    #   * run_legacy_venet_null gates the descriptive raw-signature null, which
    #     the ISI does not consume -- it computes its own residualised null.
    # This is measured, not assumed: 12_partition_variance.py runs this same
    # config at split_seed=0 and reproduces nsclc_v3's primary to full double
    # precision (0.318240180054566, CI 0.26138392827249185-0.37633421961128716).
    # Both are off here purely for runtime -- the Venet null alone is ~a third of
    # each arm, and this script pays for four arms.
    cfg = experiment.AuditConfig(
        n_folds=5, seed=0, n_boot=args.n_boot, n_null_sets=args.n_null,
        run_combat=False, run_legacy_venet_null=False,
        scorer=scorer, disattenuate=disattenuate,
    )
    t0 = time.time()
    res = experiment.run_audit(
        frame, X, sig_cols, cfg, expression=expr, signature_set=sigs,
        global_axis=axis, survival=None, verbose=False)
    return res, time.time() - t0


def summarise(res):
    """Pooled estimate and per-signature table from one audit, or empty markers.

    The table is returned by reference, not copied: `immune_excess.attrs` carries
    the raw null draws, which the reconstruction below needs.
    """
    est = res.primary_excess
    ie = res.immune_excess
    pooled = {
        "isi": None if est is None else float(est.value),
        "lo": None if est is None else float(est.lo),
        "hi": None if est is None else float(est.hi),
        "n_signatures": None if est is None else int(est.n),
        "method": None if est is None else est.method,
    }
    return pooled, (pd.DataFrame() if ie is None else ie)


def fmt(pooled):
    if pooled["isi"] is None or not np.isfinite(pooled["isi"]):
        return "NOT COMPUTED"
    lo, hi = pooled["lo"], pooled["hi"]
    if lo is None or not np.isfinite(lo):
        return f"{pooled['isi']:+.4f} [CI not reportable]"
    return f"{pooled['isi']:+.4f} [{lo:+.4f}, {hi:+.4f}]"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-null", type=int, default=1000,
                    help="random gene sets per signature; 1000 matches the main analysis")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--n-splits", type=int, default=20,
                    help="split-half repeats for the scorer-agnostic reliability")
    ap.add_argument("--n-null-reliability", type=int, default=25,
                    help="null draws whose split-half reliability is measured; the "
                         "MEAN of these is used as the null reliability, so this "
                         "need only be enough to estimate a mean")
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--arms", nargs="+", choices=ARMS, default=list(ARMS),
                    metavar="ARM",
                    help="which of the four <scorer>__<correction> arms to run; "
                         f"default all. Choices: {', '.join(ARMS)}. A partial run "
                         "writes a partial summary.json — the skipped arms are "
                         "ABSENT, not null, and `arms_requested` records what was "
                         "asked for.")
    args = ap.parse_args()

    # Single source of truth for what a partial run does; `plan_arms` is pure
    # and is what the tests exercise, so this path cannot drift away from them.
    plan = plan_arms(args.arms)
    selected = plan["selected"]
    scorers_run = plan["scorers_run"]
    recon_scorers = plan["recon_scorers"]

    outdir = args.outdir or (ROOT / "results" / "scorer_sensitivity")
    outdir.mkdir(parents=True, exist_ok=True)

    if not plan["complete"]:
        print(f"!! PARTIAL RUN: {len(selected)}/{len(ARMS)} arms — "
              f"{', '.join(selected)}", flush=True)
        print(f"   skipping: {', '.join(plan['skipped'])}", flush=True)
        if not plan["agreement_possible"]:
            print("   cross-scorer agreement (report section 4) needs both "
                  "scorers' uncorrected arms and will be empty", flush=True)

    cohort, X, expr, sigs = load_nsclc()
    filtered = sigs.filter_to(set(expr.columns))
    print(f"nsclc: {len(cohort)} patients, {cohort['tss'].nunique()} sites, "
          f"{len(filtered)} signatures, {expr.shape[1]} genes", flush=True)

    # MEASURE, THEN COMMIT. One real 200-draw ssGSEA family, extrapolated, so the
    # n_null_sets on the banner is a measured choice and not an aspiration.
    probe_single = signatures.SignatureSet(
        name="probe", sets={"p": next(iter(filtered.sets.values()))})
    probe_family = signatures.random_gene_sets(
        sorted(expr.columns), probe_single, n_per_signature=200, seed=0,
        match_expression=expr.mean(axis=0))["p"]
    t0 = time.time()
    signatures.score_ssgsea(expr, probe_family)
    t_probe = time.time() - t0
    ms_per_set = 1000 * t_probe / (200 * len(expr))
    print(f"ssGSEA scoring: {t_probe:.1f}s for 200 sets x {len(expr)} patients "
          f"= {ms_per_set:.4f} ms per patient-set "
          f"(the pre-optimisation loop cost 6.26)", flush=True)
    print(f"  extrapolated ssGSEA scoring for the whole run at n_null={args.n_null}: "
          f"{len(filtered) * args.n_null * len(expr) * ms_per_set / 1000:.0f}s "
          f"per arm, x2 ssGSEA arms", flush=True)
    print(f"n_null_sets = {args.n_null}", flush=True)
    if args.n_null < 1000:
        print("  *** NOT EQUIVALENT TO THE MAIN ANALYSIS, which uses 1000 draws. "
              "A reduced null widens every null mean's standard error and, below "
              "319 draws, cannot clear BH-FDR at 16 signatures at all. ***",
              flush=True)

    axis = globalaxis.compute_global_axis(
        expr, method="pc1", cancer_type=cohort.set_index("patient_id")["cancer_type"])
    axis_values = pd.Series(axis.values).reindex(expr.index)
    axis_values.index = expr.index
    within = pd.Series(cohort["cancer_type"].to_numpy(), index=expr.index)

    ref = selfcheck_split_half(
        expr, filtered.sets[sorted(filtered.sets)[0]], axis_values, within)
    print(f"split-half self-check passed (batched == library, {ref:.4f})", flush=True)

    results, timings, tables = {}, {}, {}
    reliab: dict[str, dict] = {}

    for scorer in scorers_run:
        t0 = time.time()
        frame, sig_cols = rescore_cohort(cohort, expr, sigs, scorer)
        print(f"\n[{scorer}] observed signatures rescored in {time.time() - t0:.1f}s; "
              f"score sd {frame[sig_cols].std().mean():.4f}", flush=True)
        if scorer == "mean_z":
            # The parquet's stored columns should already BE these. If they are
            # not, every frozen result rests on a scoring path this script cannot
            # reproduce, which matters far more than the sensitivity analysis.
            drift = float(np.nanmax(np.abs(
                frame[sig_cols].to_numpy() - cohort[sig_cols].to_numpy())))
            print(f"  max drift vs the cohort's stored mean-z columns: {drift:.2e}",
                  flush=True)

        for dis in (True, False):
            key = f"{scorer}__{'disattenuated' if dis else 'uncorrected'}"
            if key not in selected:
                print(f"  {key:32s} SKIPPED (not in --arms)", flush=True)
                continue
            res, dt = run_arm(frame, X, sig_cols, expr, sigs, axis,
                              scorer=scorer, disattenuate=dis, args=args)
            pooled, table = summarise(res)
            results[key], timings[key], tables[key] = pooled, dt, table
            beats = int(table["beats_null"].sum()) if len(table) else 0
            finite = int(np.isfinite(table["excess_z"]).sum()) if len(table) else 0
            print(f"  {key:32s} {dt:6.0f}s  pooled {fmt(pooled):34s} "
                  f"finite {finite}/{len(filtered)}  beats null {beats}", flush=True)
            for note in res.notes:
                if "DISATTENUAT" in note.upper() or "NOT COMPUTED" in note.upper():
                    print(f"      note: {note[:150]}", flush=True)

        # Split-half reliability feeds ONLY the reconstruction, which needs this
        # scorer's uncorrected arm. Computing it anyway would spend ~20 split-half
        # repeats plus 25 null draws on a number nothing downstream reads.
        if scorer not in recon_scorers:
            print(f"  [{scorer}] split-half reliability SKIPPED — its uncorrected "
                  "arm was not selected, so there is nothing to reconstruct",
                  flush=True)
            continue

        # Scorer-agnostic reliability, on the observed sets and on a sample of
        # null draws from the SAME families the audit used — same pool ordering,
        # same seed, same expression matching, so draw i here IS draw i there.
        t0 = time.time()
        score_fn = getattr(signatures, f"score_{scorer}")
        obs_rel = {k: float(v) for k, v in split_half_reliabilities(
            expr, filtered.sets, score_fn, axis_values, within,
            n_splits=args.n_splits, seed=0).items()}

        # ONE call for every signature's null draws together. Per signature it
        # would be 16 ssGSEA calls each paying the fixed whole-matrix rank cost,
        # which is most of the bill; the draws are named "<signature>__random_NNN"
        # so they regroup unambiguously afterwards.
        sampled: dict[str, list[str]] = {}
        for name, genes in filtered.sets.items():
            single = signatures.SignatureSet(name=name, sets={name: genes})
            fam = signatures.random_gene_sets(
                sorted(expr.columns), single, n_per_signature=args.n_null, seed=0,
                match_expression=expr.mean(axis=0))[name]
            sampled.update(dict(list(fam.sets.items())[: args.n_null_reliability]))
        per_draw = split_half_reliabilities(
            expr, sampled, score_fn, axis_values, within,
            n_splits=args.n_splits, seed=0)
        null_rel = {}
        for name in filtered.sets:
            vals = [v for k, v in per_draw.items() if k.startswith(f"{name}__random_")]
            null_rel[name] = float(np.nanmean(vals)) if vals else float("nan")

        reliab[scorer] = {"observed": obs_rel, "null_mean": null_rel}
        o_bar = float(np.nanmean(list(obs_rel.values())))
        n_bar = float(np.nanmean(list(null_rel.values())))
        print(f"  split-half reliability ({args.n_splits} splits, "
              f"{args.n_null_reliability} null draws/signature): "
              f"observed {o_bar:.4f}  null {n_bar:.4f}  gap {o_bar - n_bar:+.4f}"
              f"  ({time.time() - t0:.0f}s)", flush=True)

    # ------------------------------------------------------ reconstruction ---
    # The registered estimand, recomputed from the audit's own stored null draws
    # with a reliability the scorer actually supports. `excess_over_null` is the
    # same function the audit calls; only the reliability input changes. The
    # interval is its analytic one, NOT the patient-clustered bootstrap, because
    # the bootstrap needs per-patient predictions the audit does not return.
    recon = {}
    for scorer in recon_scorers:
        table = tables[f"{scorer}__uncorrected"]
        draws = table.attrs.get("null_draws", {})
        rows = []
        for _, row in table.iterrows():
            name = row["signature"]
            nulls = np.asarray(draws.get(name, []), dtype=float)
            rel_o = reliab[scorer]["observed"].get(name, np.nan)
            rel_n = reliab[scorer]["null_mean"].get(name, np.nan)
            if not len(nulls) or not np.isfinite(rel_o) or not np.isfinite(rel_n):
                continue
            est = globalaxis.excess_over_null(
                float(row["r_residual_refit"]), nulls, n=int(row["n"]),
                reliability_observed=rel_o, reliability_null=rel_n)
            rows.append({"signature": name, "excess_z": est.value,
                         "excess_lo": est.lo, "excess_hi": est.hi,
                         "beats_null": bool(np.isfinite(est.lo) and est.lo > 0),
                         "reliability_observed": rel_o, "reliability_null": rel_n})
        d = pd.DataFrame(rows)
        vals = d["excess_z"].to_numpy(dtype=float) if len(d) else np.array([])
        vals = vals[np.isfinite(vals)]
        se = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else np.nan
        recon[scorer] = {
            "isi": float(np.mean(vals)) if len(vals) else None,
            "lo": float(np.mean(vals) - CRIT * se) if len(vals) > 1 else None,
            "hi": float(np.mean(vals) + CRIT * se) if len(vals) > 1 else None,
            "n_signatures": int(len(vals)),
            "beats_null": int(d["beats_null"].sum()) if len(d) else 0,
            "method": "excess-over-null-z-disattenuated-with-split-half-reliability; "
                      "CI is the across-signature SE, NOT the registered "
                      "patient-clustered bootstrap",
            "per_signature": d,
        }

    # -------------------------------------------------------- comparisons ---
    keep_cols = ["excess_z", "excess_lo", "excess_hi", "beats_null",
                 "r_residual_refit", "null_mean_r", "null_sd_r",
                 "alpha_observed", "alpha_null_mean"]
    pieces = []
    for key, table in tables.items():
        if not len(table):
            continue
        have = [c for c in keep_cols if c in table.columns]
        pieces.append(table[["signature"] + have].rename(
            columns={c: f"{c}__{key}" for c in have}))
    for scorer in recon_scorers:
        d = recon[scorer]["per_signature"]
        if len(d):
            pieces.append(d.rename(columns={
                c: f"{c}__{scorer}__reconstructed" for c in d.columns
                if c != "signature"}))
    # Drop `.attrs` before merging. `AuditResult` frames carry the embedding
    # matrix under attrs["_X"] and the 16x1000 null correlations under
    # attrs["null_draws"], and pandas' `__finalize__` decides whether to
    # propagate attrs with `all(obj.attrs == attrs for obj in ...)`. On a
    # numpy array that `==` is elementwise, so `all()` raises
    #   ValueError: The truth value of an array with more than one element is
    #   ambiguous
    # and the whole run dies AFTER every arm has been computed -- which is
    # exactly what happened on HPC4 job 116796: four arms, 34 minutes, zero
    # files written. None of these attrs belong in a per-signature summary
    # table, so the fix is to drop them rather than to make them comparable.
    pieces = [p.copy() for p in pieces]
    for p in pieces:
        p.attrs = {}
    merged = pieces[0]
    for piece in pieces[1:]:
        merged = merged.merge(piece, on="signature", how="outer")
        merged.attrs = {}
    merged = merged.sort_values("signature").reset_index(drop=True)
    merged.to_csv(outdir / "per_signature.csv", index=False)

    def agreement(col_a, col_b, label):
        if col_a not in merged.columns or col_b not in merged.columns:
            return {"label": label, "n": 0, "pearson": None,
                    "spearman": None, "mean_abs_diff": None}
        a = pd.to_numeric(merged[col_a], errors="coerce").to_numpy(dtype=float)
        b = pd.to_numeric(merged[col_b], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 3:
            return {"label": label, "n": int(ok.sum()), "pearson": None,
                    "spearman": None, "mean_abs_diff": None}
        return {
            "label": label, "n": int(ok.sum()),
            "pearson": float(np.corrcoef(a[ok], b[ok])[0, 1]),
            "spearman": float(sps.spearmanr(a[ok], b[ok]).statistic),
            "mean_abs_diff": float(np.mean(np.abs(a[ok] - b[ok]))),
        }

    agree = [
        agreement("excess_z__mean_z__uncorrected", "excess_z__ssgsea__uncorrected",
                  "uncorrected excess, mean_z vs ssgsea"),
        agreement("excess_z__mean_z__reconstructed", "excess_z__ssgsea__reconstructed",
                  "reconstructed disattenuated excess, mean_z vs ssgsea"),
        agreement("r_residual_refit__mean_z__uncorrected",
                  "r_residual_refit__ssgsea__uncorrected",
                  "residualised observed r, mean_z vs ssgsea"),
    ]

    # ------------------------------------------------------------- report ---
    print("\n" + "=" * 74)
    print(f"SCORER SENSITIVITY — NSCLC, n_null_sets={args.n_null}, "
          f"n_boot={args.n_boot}")
    print("=" * 74)
    n_sig = len(filtered)

    print("\n  1. UNCORRECTED excess (like for like; NOT the ISI, "
          "not to be quoted as one)")
    for scorer in SCORERS:
        if f"{scorer}__uncorrected" not in results:
            print(f"     {scorer:8s} NOT RUN (excluded by --arms)")
            continue
        p, t = results[f"{scorer}__uncorrected"], tables[f"{scorer}__uncorrected"]
        print(f"     {scorer:8s} {fmt(p):36s} beats null "
              f"{int(t['beats_null'].sum())}/{n_sig}")

    print("\n  2. REGISTERED disattenuated ISI, exactly as configured")
    for scorer in SCORERS:
        if f"{scorer}__disattenuated" not in results:
            print(f"     {scorer:8s} NOT RUN (excluded by --arms)")
            continue
        p, t = results[f"{scorer}__disattenuated"], tables[f"{scorer}__disattenuated"]
        finite = int(np.isfinite(t["excess_z"]).sum()) if len(t) else 0
        print(f"     {scorer:8s} {fmt(p):36s} beats null "
              f"{int(t['beats_null'].sum())}/{n_sig}  (defined for {finite}/{n_sig})")
        if len(t):
            print(f"              audit alpha: observed {np.nanmean(t['alpha_observed']):.4f}"
                  f"  null {np.nanmean(t['alpha_null_mean']):.4f}")

    print("\n  3. RECONSTRUCTED disattenuated ISI (split-half reliability, "
          "analytic CI)")
    for scorer in SCORERS:
        if scorer not in recon:
            print(f"     {scorer:8s} NOT RUN (its uncorrected arm was excluded "
                  "by --arms)")
            continue
        r = recon[scorer]
        band = ("NOT COMPUTED" if r["isi"] is None
                else f"{r['isi']:+.4f} [{r['lo']:+.4f}, {r['hi']:+.4f}]")
        print(f"     {scorer:8s} {band:36s} beats null "
              f"{r['beats_null']}/{n_sig}  (defined for {r['n_signatures']}/{n_sig})")
        print(f"              split-half reliability: observed "
              f"{np.nanmean(list(reliab[scorer]['observed'].values())):.4f}"
              f"  null {np.nanmean(list(reliab[scorer]['null_mean'].values())):.4f}")

    print("\n  4. Per-signature agreement between the scorers")
    for a in agree:
        if a["pearson"] is None:
            print(f"     {a['label']:52s} n={a['n']} — too few finite pairs")
        else:
            print(f"     {a['label']:52s} n={a['n']}  r={a['pearson']:+.3f}  "
                  f"rho={a['spearman']:+.3f}  mean|diff|={a['mean_abs_diff']:.4f}")

    print("\n  runtimes (s): " + "  ".join(f"{k}={v:.0f}" for k, v in timings.items()))

    payload = {
        "cohort": "nsclc",
        "n_patients": int(len(cohort)),
        "n_signatures": int(n_sig),
        "n_null_sets": int(args.n_null),
        "n_null_sets_matches_main_analysis": bool(args.n_null == 1000),
        "n_boot": int(args.n_boot),
        "ssgsea_ms_per_patient_set": float(1000 * t_probe / (200 * len(expr))),
        # A partial run records what was ASKED FOR alongside what it produced, so
        # a reader can tell "this arm was excluded" from "this arm failed".
        "arms_requested": list(selected),
        "arms_complete": bool(plan["complete"]),
        "arms": {k: {kk: vv for kk, vv in v.items()} for k, v in results.items()},
        # Every arm appears here with an explicit status, so a partial
        # summary.json distinguishes "excluded by --arms" from "requested and
        # missing" without the reader inferring either from an absent key.
        "arm_status": [{"arm": r["arm"], "status": r["status"]}
                       for r in arm_status_rows(plan, results)],
        "reconstructed": {
            s: {k: v for k, v in recon[s].items() if k != "per_signature"}
            for s in recon_scorers
        },
        "split_half_reliability": {
            s: {"observed_mean": float(np.nanmean(list(reliab[s]["observed"].values()))),
                "null_mean": float(np.nanmean(list(reliab[s]["null_mean"].values()))),
                "per_signature_observed": reliab[s]["observed"],
                "per_signature_null_mean": reliab[s]["null_mean"]}
            for s in recon_scorers
        },
        "beats_null": {k: int(t["beats_null"].sum()) if len(t) else 0
                       for k, t in tables.items()},
        "agreement": agree,
        "runtime_s": timings,
        "caveat": (
            "The registered disattenuated ISI is UNDEFINED under scorer='ssgsea': "
            "the audit estimates reliability with the general-form Cronbach alpha "
            "for a mean-of-z-scores composite, which an ssGSEA rank-walk score is "
            "not, so alpha collapses to 0 and every draw is dropped. This is a "
            "property of the reliability estimator, not of ssGSEA. Arms 1 and 3 "
            "are the comparable ones."
        ),
    }
    (outdir / "summary.json").write_text(json.dumps(payload, indent=2, default=str))
    print(f"\nWritten to {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
