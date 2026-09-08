#!/usr/bin/env python3
"""Close the pre-registered analyses that were computed-and-discarded or skipped.

    python scripts/11_close_science_gaps.py --results results/pancancer_v2

Found by the 2026-08-19 science audit (`../14-SCIENCE-AUDIT.md`). Each of these
was committed to in the pre-registration or the brief and is absent from the
result files:

  A1  BH-FDR at q=0.05 across the signature family. `beats_null` was an
      UNCORRECTED 95% interval. This is the one that can change a headline: the
      outcome arm's five winners include hypoxia at excess_lo = +0.0069.
  A2  `stats.rotation_null` — computed on every run and stored in
      `DataFrame.attrs`, which `to_csv` silently drops.
  A3  `stats.effective_tests` (Li & Ji M_eff), promised descriptively.
  A4  Control C, label-side site variance — brief §4C, never called by a runner.
  A6  The small-panel bracket. "Small clinical panels are most exposed" is
      extrapolated from k=36..200; nothing below 36 was measured.

Everything except A4 and A6 is computed post-hoc from the saved result files, so
no model is refitted and no number can drift.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import decomposition, signatures as sig_mod, stats  # noqa: E402

GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"


def join_scores_to_meta(meta, scored):
    """Join per-patient signature scores onto the clinical table.

    Returns a frame carrying `patient_id` as a COLUMN, whatever the inputs'
    index names are.

    The defect this pins, found 2026-09-06 by running A4 for the first time
    since the cohort-scope fix shipped: `DataFrame.join` preserves the LEFT
    frame's index name only when the right frame's index is named too.
    `scored` is built by reassigning a bare `pd.Index`, which is unnamed, so
    the joined index came back anonymous and `reset_index()` produced a column
    called `index`. `frame["patient_id"]` then raised KeyError, A4 died, A4b
    cascaded off A4's unset `cols`, and because each stage swallows its own
    exception the script still exited 0 and wrote a `science_gaps.json` full of
    nulls that read as a completed artefact.

    It was invisible for a session because the cohort-scope fix was verified
    against pan-cancer, where the patient filter is a measured no-op -- but the
    KeyError is upstream of the filter and cohort-independent, so the one case
    checked was the one case whose *filter* could not misbehave while the
    *join* was already broken.

    `rename_axis` restores the key unconditionally rather than depending on
    what the join chose to keep.
    """
    joined = meta.join(scored, how="inner").rename_axis("patient_id").reset_index()
    ids = joined["patient_id"]
    if ids.isna().any() or ids.duplicated().any():
        raise AssertionError(
            f"the score join produced {int(ids.isna().sum())} null and "
            f"{int(ids.duplicated().sum())} duplicate patient keys; a "
            "per-patient control cannot be computed from that")
    return joined


def repo_relative(p) -> str:
    """Render a path relative to the pipeline root, for anything written to JSON.

    `--cohort` and `--expr` default to `ROOT / "data/interim/..."`, i.e. ABSOLUTE
    paths containing the author's home directory. Writing `str(args.cohort)` put
    `/Users/<user>/...` into `science_gaps.json`, where it is two things at once:
    a personal-path leak that `22_build_public_snapshot.py` refuses to publish
    (it caught this on 2026-09-07, the first rebuild after the field was added),
    and a value that is meaningless on any other machine. `results_dir` beside it
    was relative only by accident of how the script happened to be invoked.

    Falls back to the basename when the path is outside the repo, so a caller
    pointing at data elsewhere still cannot write their home directory into a
    published artefact.
    """
    p = Path(p)
    try:
        return str(p.resolve().relative_to(ROOT))
    except ValueError:
        return p.name


def fingerprint_patients(patients) -> str:
    """SHA-256 of the sorted unique patient ids.

    Deliberately identical in construction to
    `aacr27.experiment.AuditRun.cohort_fingerprint`, which is what writes the
    digest this module compares against. If the two ever diverge the check
    below turns into a permanent false alarm, so
    `test_science_gaps_fingerprint_matches_the_writer` pins them together by
    running BOTH on the same input rather than trusting that they agree.
    """
    import hashlib

    return hashlib.sha256(
        "\n".join(sorted({str(p) for p in patients})).encode()).hexdigest()


def cohort_fingerprint_status(summary: dict | None, run_patients) -> dict:
    """Compare a run's recorded cohort against the one its predictions show.

    PURE: no I/O, so the mismatch branch is reachable in a test. That matters
    more than usual here -- as of 2026-09-07 NO results directory on disk
    carries the fingerprint (the writer was added on 2026-09-06 and no run has
    been executed since), so every real invocation takes the "absent" branch.
    A check exercised only where it is a no-op is not a check, which is the
    precise defect that let A11 through in the first place.

    Three outcomes:

    "absent"    -- no fingerprint recorded. Every run frozen before 2026-09-06
                   is in this state, and it means "predates the writer", NOT
                   "cohort unknown". Never an error.
    "match"     -- summary.json and predictions.csv.gz describe the same
                   patients.
    "mismatch"  -- they do not. The directory is internally inconsistent: this
                   is A11's defect class, and A11 itself was caught only
                   because three unrelated values happened to agree with
                   another cohort's to seventeen significant figures.

    `summary["cohort"]` is NOT assumed to be the fingerprint. Several other
    scripts (12_partition_variance.py, 13_scorer_sensitivity.py) already write
    `"cohort": "nsclc"` -- a bare STRING -- into their own summaries, so the
    key is overloaded. Only a mapping carrying `patient_set_sha256` counts;
    anything else reads as absent.
    """
    block = (summary or {}).get("cohort")
    if not isinstance(block, dict) or "patient_set_sha256" not in block:
        return {"status": "absent", "recorded": None,
                "observed_n": len(set(map(str, run_patients)))}

    observed = fingerprint_patients(run_patients)
    recorded = str(block.get("patient_set_sha256"))
    observed_n = len(set(map(str, run_patients)))
    recorded_n = block.get("n_patients")
    ok = observed == recorded and recorded_n == observed_n
    return {
        "status": "match" if ok else "mismatch",
        "recorded": recorded, "observed": observed,
        "recorded_n": recorded_n, "observed_n": observed_n,
    }


def read_and_check_cohort(res):
    """Read a results directory's cohort and check it against its own summary.

    Does the I/O that `cohort_fingerprint_status` deliberately avoids: recovers
    the run's patients from `predictions.csv.gz`, loads `summary.json` if it is
    there, compares the two, and refuses to continue on a mismatch.

    Extracted from `main()` so the *integration* can be exercised on a fixture
    directory. The pure function's branches were unit-tested from the day it
    was written, but this wiring -- which file is read, which column holds the
    patient id, whether the digest is compared against the right set, and
    whether a mismatch actually stops the run -- had never executed on a
    non-absent case, because no results directory on disk carries a
    fingerprint. Four tested branches behind untested plumbing is not a guard;
    it is a guard-shaped object, which is the A11 defect class exactly.

    Returns `(preds, run_patients, status)`. Raises SystemExit on a mismatch:
    every downstream stage is keyed to `run_patients`, so if the summary and
    the predictions disagree there is no defensible answer to "whose data is
    this", and writing gap numbers into the directory would repeat A11.
    """
    preds = pd.read_csv(res / "predictions.csv.gz")
    # The run's own patient set, straight from its predictions. Before the
    # writer existed this was the ONLY cohort record a results directory
    # carried -- config.json and summary.json both describe the estimator's
    # settings but never say which patients it was fitted on.
    run_patients = set(preds["patient"].astype(str))

    summary_path = res / "summary.json"
    fp = cohort_fingerprint_status(
        json.loads(summary_path.read_text()) if summary_path.exists() else None,
        run_patients)
    if fp["status"] == "mismatch":
        raise SystemExit(
            f"\nFATAL: {res.name}/summary.json records a cohort of "
            f"{fp['recorded_n']} patients (sha256 {str(fp['recorded'])[:16]}...) "
            f"but predictions.csv.gz contains {fp['observed_n']} "
            f"(sha256 {fp['observed'][:16]}...). The directory describes two "
            "different cohorts; refusing to write gap numbers into it.")
    if fp["status"] == "match":
        print(f"  cohort fingerprint: MATCHES summary.json "
              f"({fp['observed_n']} patients)")
    else:
        print("  cohort fingerprint: not recorded (run predates 2026-09-06) "
              "— cohort taken from predictions.csv.gz")
    return preds, run_patients, fp


def cohort_scope_status(available, run_patients, scope: str) -> dict:
    """Resolve which patients a cohort-derived stage may use. PURE: no I/O.

    `available` is whatever patients the cohort/expression tables offer;
    `run_patients` is the set the run in `--results` was actually fitted on,
    recovered from its predictions. Under scope "run" the stage is restricted
    to their intersection; under "file" it keeps everything, which is the
    pre-2026-09-06 behaviour that let a `--results results/nsclc_v3_stablesort`
    run report 202-site pan-TCGA numbers under an NSCLC heading.

    Returned separately from the filtering itself so a test can assert on the
    decision without building a cohort, and so the numbers reach science_gaps
    .json -- the mismatch hid precisely because nothing recorded the cohort.
    """
    if scope not in ("run", "file"):
        raise ValueError(f"unknown cohort scope {scope!r}; expected 'run'/'file'")
    avail = {str(p) for p in available}
    run = {str(p) for p in run_patients}
    selected = (avail & run) if scope == "run" else avail
    return {
        "scope": scope,
        "selected": selected,
        "n_available": len(avail),
        "n_run": len(run),
        "n_selected": len(selected),
        # True only when the stage describes exactly the run's patients. Under
        # scope "file" a superset is the normal case, and it is still a mismatch.
        "matches_run": selected == run,
        "n_missing_from_cohort": len(run - avail),
    }


def two_sided_p(value: pd.Series, lo: pd.Series, hi: pd.Series) -> np.ndarray:
    """p from an estimate and its 95% interval, via the implied standard error."""
    se = (hi - lo) / (2 * sps.norm.ppf(0.975))
    with np.errstate(invalid="ignore", divide="ignore"):
        return 2 * sps.norm.sf(np.abs(value / se))


def add_bh(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    """A1 — Benjamini-Hochberg across the family, as pre-registered."""
    out = frame.copy()
    p = two_sided_p(out["excess_z"], out["excess_lo"], out["excess_hi"])
    out["p_two_sided"] = p
    ok = np.isfinite(p)
    q = np.full(len(out), np.nan)
    rej = np.zeros(len(out), dtype=bool)
    if ok.any():
        rej_ok, q_ok = stats.bh_fdr(p[ok])
        q[ok], rej[ok] = q_ok, rej_ok
    out["q_bh"] = q
    # An estimate can be significant by BH yet point the WRONG way; "beats null"
    # must mean significantly ABOVE it.
    out["beats_null_bh"] = rej & (out["excess_z"] > 0).to_numpy()
    n_before = int(out["beats_null"].sum()) if "beats_null" in out else -1
    n_after = int(out["beats_null_bh"].sum())
    print(f"  {label}: uncorrected {n_before}/{len(out)} -> BH q<0.05 {n_after}/{len(out)}")
    if n_before != n_after:
        print(f"    ** BH CHANGED THE COUNT ({n_before} -> {n_after}) **")
        lost = out[out.get("beats_null", False) & ~out["beats_null_bh"]]
        for _, r in lost.iterrows():
            print(f"       dropped: {r['signature']} "
                  f"excess={r['excess_z']:+.4f} q={r['q_bh']:.4f}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=ROOT / "results" / "pancancer_v2")
    ap.add_argument("--cohort", type=Path, default=ROOT / "data/interim/meta.parquet")
    ap.add_argument("--expr", type=Path, default=ROOT / "data/interim/expression_hugo.parquet")
    # `--cohort`/`--expr` are the WHOLE pan-TCGA tables and do not know which
    # cohort `--results` came from. Before 2026-09-06 the cohort-derived stages
    # (A4, A4b, A6) silently scored all 7,168 pan-TCGA patients no matter what
    # `--results` pointed at, so a `--results results/nsclc_v3_stablesort` run
    # wrote pan-cancer numbers into an NSCLC directory under an NSCLC heading.
    # Nothing in the output recorded the cohort, which is exactly how it hid.
    # "run" restricts to the patients the run actually used, recovered from
    # predictions.csv.gz. "file" is the pre-fix behaviour, kept because
    # pancancer_v3/science_gaps.json's A6 block was produced under it.
    ap.add_argument("--cohort-scope", choices=("run", "file"), default="run",
                    help="which patients the cohort-derived stages (A4, A4b, A6) "
                         "use: 'run' = only those in --results (default, correct); "
                         "'file' = every patient in --cohort/--expr (pre-2026-09-06 "
                         "behaviour; needed to reproduce pancancer_v3's A6 block)")
    args = ap.parse_args()
    res = args.results
    report: dict[str, object] = {"results_dir": str(res)}
    failures: list[tuple[str, str]] = []

    print(f"=== closing science gaps for {res.name} ===\n")

    # ---------------------------------------------------------------- A1 ---
    print("A1  BH-FDR across the family (pre-registration section 4)")
    ie = pd.read_csv(res / "immune_excess.csv")
    ie = add_bh(ie, "ISI")
    ie.to_csv(res / "immune_excess.csv", index=False)
    report["isi_beats_null_bh"] = int(ie["beats_null_bh"].sum())

    oa_path = res / "outcome_arm.csv"
    if oa_path.exists():
        oa = pd.read_csv(oa_path)
        oa = add_bh(oa, "outcome arm")
        oa.to_csv(oa_path, index=False)
        report["outcome_beats_null_uncorrected"] = int(oa["beats_null"].sum())
        report["outcome_beats_null_bh"] = int(oa["beats_null_bh"].sum())

    # ---------------------------------------------------------------- A2 ---
    print("\nA2  rotation null (family-level p, respects between-signature correlation)")
    draws_path = res / "null_draws.npz"
    if draws_path.exists():
        z = np.load(draws_path)
        names = list(z.files)
        width = min(len(z[n]) for n in names)
        R = np.vstack([np.asarray(z[n][:width], dtype=float) for n in names])
        observed = float(np.nanmean(ie["r_residual_refit"]))
        rot = stats.rotation_null(R, observed)
        print(f"  {len(names)} signatures x {width} draws")
        print(f"  observed family mean r = {rot['observed']:.4f}; null family mean "
              f"{rot['null_family_mean']:.4f} [{rot['ci_lo']:.4f}, {rot['ci_hi']:.4f}]")
        print(f"  p = {rot['p']:.4g}  (B = {rot['B']})")
        report["rotation_null"] = rot
    else:
        print(f"  SKIPPED: {draws_path.name} not present. Null draws are only "
              "persisted by runs after 2026-08-19; re-run to enable.")
        report["rotation_null"] = None

    # ---------------------------------------------------------------- A3 ---
    print("\nA3  effective number of tests (Li & Ji) — DESCRIPTIVE ONLY")
    # Make the fingerprint LOAD-BEARING. Writing it was the source-level fix for
    # A11; until something reads it, a directory could still carry a summary
    # describing one cohort beside predictions from another and nothing would
    # object. The read, the comparison and the refusal all live in
    # `read_and_check_cohort` so that a fixture test can drive them.
    preds, run_patients, fp = read_and_check_cohort(res)
    report["run_n_patients"] = len(run_patients)
    report["cohort_fingerprint"] = fp
    wide = (preds[preds["model"] == "ridge_embedding"]
            .pivot_table(index="patient", columns="target", values="y_true"))
    corr = wide.corr().to_numpy()
    corr = np.nan_to_num(corr, nan=0.0)
    meff = stats.effective_tests(corr)
    print(f"  {corr.shape[0]} signatures -> M_eff = {meff:.2f}")
    print("  Reported as a description of panel redundancy. NOT used to divide "
          "alpha — combining that with BH is a category error.")
    report["effective_tests"] = float(meff)
    report["n_signatures"] = int(corr.shape[0])

    # ---------------------------------------------------------------- A4 ---
    print("\nA4  Control C — label-side site variance (brief section 4C, RNA only)")
    try:
        meta = pd.read_parquet(args.cohort)
        sig_cols = [c for c in pd.read_csv(res / "immune_excess.csv")["signature"]]
        expr = pd.read_parquet(args.expr)
        expr.index = pd.Index(expr.index).astype(str).str[:12]
        expr = expr.groupby(level=0).mean()
        sets = sig_mod.SignatureSet.from_gmt(GMT)
        sets = sig_mod.SignatureSet(name=sets.name,
                                    sets={f"SIG_{k}": v for k, v in sets.sets.items()})
        scored = sig_mod.score_mean_z(expr, sets.filter_to(set(expr.columns)))
        frame = meta.drop_duplicates("patient_id").set_index("patient_id")
        frame = join_scores_to_meta(frame, scored)
        frame = frame.rename(columns={"cancer type abbreviation": "cancer_type"})
        # Restrict to the patients the run actually used. Filtering HERE rather
        # than before score_mean_z is load-bearing: the scorer z-scores across
        # whatever patients it is handed, so filtering earlier would move the
        # baseline. Measured 2026-09-06: meta INNER JOIN expr is exactly 7,168
        # patients, which is precisely pancancer_v3's run, so for pan-cancer
        # this filter is a no-op and the frozen A4/A4b numbers still reproduce.
        st = cohort_scope_status(frame["patient_id"], run_patients,
                                 args.cohort_scope)
        frame = frame[frame["patient_id"].astype(str).isin(st["selected"])]
        report["cohort_path"] = repo_relative(args.cohort)
        report["cohort_scope"] = st["scope"]
        report["cohort_n_patients_in_file"] = st["n_available"]
        report["cohort_n_patients_used"] = st["n_selected"]
        report["cohort_matches_run"] = st["matches_run"]
        print(f"  cohort scope '{st['scope']}': {st['n_selected']} of "
              f"{st['n_available']} joined patients used; the run itself used "
              f"{st['n_run']}")
        if st["n_selected"] == 0:
            raise ValueError(
                f"no overlap between {args.cohort} and the {st['n_run']} "
                f"patients in {res.name}/predictions.csv.gz -- wrong --cohort?")
        if not st["matches_run"]:
            # Loud, but not fatal: a legitimately smaller overlap (a patient
            # lacking expression) is a result, not a crash. A SILENT mismatch
            # is what produced pan-cancer numbers in an NSCLC directory.
            print(f"  !! COHORT MISMATCH: the stages below describe "
                  f"{st['n_selected']} patients, the run used {st['n_run']} "
                  f"({st['n_missing_from_cohort']} absent from "
                  f"{args.cohort.name}). These are DIFFERENT patient sets.")
        cols = [c for c in sig_cols if c in frame.columns]
        lsv = decomposition.label_site_variance(frame, cols, plate_col="plate")
        lsv.to_csv(res / "label_site_variance.csv", index=False)
        print(f"  {len(lsv)} signatures, {int(lsv['n_sites'].iloc[0])} sites")
        print(f"  median R^2 of SITE on the LABEL (image-free): "
              f"{lsv['r2_site_only'].median():.4f}")
        if "r2_site_given_type" in lsv:
            print(f"  median R^2 of site GIVEN cancer type:      "
                  f"{lsv['r2_site_given_type'].median():.4f}")
        if "r2_plate_within_site" in lsv:
            print(f"  median R^2 of plate WITHIN site (technical floor): "
                  f"{lsv['r2_plate_within_site'].median():.4f}")
        report["label_site_variance_median_r2"] = float(lsv["r2_site_only"].median())
    except Exception as exc:  # noqa: BLE001 — a failed control is a result, not a crash
        print(f"  FAILED: {exc}")
        report["label_site_variance_median_r2"] = None
        failures.append(("A4", f"{type(exc).__name__}: {exc}"))

    # ---------------------------------------------------- A4b calibration ---
    print("\nA4b Calibrating Control C against permuted site labels")
    print("  Dummy-coded site inflates R^2 by roughly p/(n-1) BY CONSTRUCTION;")
    print("  with 200+ sites that is several points of free R^2, and it lands on")
    print("  the covariate rungs, so an uncalibrated number overstates site.")
    try:
        probe = cols[0]
        cal = decomposition.permutation_calibration(
            frame, probe, site_col="tss", within_col="cancer_type", n_perm=200)
        print(f"  {probe.replace('SIG_HALLMARK_', '')}: observed R^2 "
              f"{cal['observed_r2']:.4f}, permuted median {cal['perm_median_r2']:.4f}"
              f" -> CALIBRATED {cal['calibrated_r2']:.4f}  (p={cal['perm_p']:.4g})")
        report["control_c_calibrated_r2"] = cal["calibrated_r2"]
        report["control_c_perm_median"] = cal["perm_median_r2"]
    except Exception as exc:  # noqa: BLE001
        print(f"  FAILED: {exc}")
        failures.append(("A4b", f"{type(exc).__name__}: {exc}"))

    # ---------------------------------------------------------------- A6 ---
    print("\nA6  Small-panel bracket — is 'small panels are most exposed' MEASURED?")
    print("  The claim is extrapolated from k=36..200. Subsampling the Hallmark")
    print("  sets down to k=10 and 20 tests it inside the measured range.")
    try:
        av = pd.read_csv(res / "immune_excess.csv")
        # A6 z-scores over PATIENTS, so the patient set changes every alpha.
        # Unlike A4 there is no join to narrow it: `expr` still holds all 9,781
        # patients with expression, which is not even pancancer_v3's 7,168. Under
        # scope 'run' this now matches the run; under 'file' it reproduces the
        # stored pancancer_v3 A6 block, which was computed on all 9,781.
        a6 = cohort_scope_status(expr.index, run_patients, args.cohort_scope)
        a6_expr = expr[expr.index.astype(str).isin(a6["selected"])]
        report["a6_n_patients"] = a6["n_selected"]
        print(f"  alpha computed over {len(a6_expr)} patients "
              f"(scope '{args.cohort_scope}')")
        z = (a6_expr - a6_expr.mean(axis=0)) / a6_expr.std(axis=0).replace(0, np.nan)
        z = z.dropna(axis=1, how="all").fillna(0.0)
        rng = np.random.default_rng(0)
        pool = sorted(z.columns)
        rows = []
        for k in (10, 20, 40, 80, 160):
            gaps = []
            for name in av["signature"]:
                genes = [g for g in sets.sets.get(name, []) if g in z.columns]
                if len(genes) < k:
                    continue
                sub = list(rng.choice(genes, size=k, replace=False))
                rand = list(rng.choice(pool, size=k, replace=False))
                a_c = sig_mod.cronbach_alpha(z, sub)
                a_r = sig_mod.cronbach_alpha(z, rand)
                if np.isfinite(a_c) and np.isfinite(a_r):
                    gaps.append(a_c - a_r)
            if gaps:
                rows.append({"k": k, "n_sets": len(gaps),
                             "median_gap": float(np.median(gaps)),
                             "mean_gap": float(np.mean(gaps))})
        bracket = pd.DataFrame(rows)
        bracket.to_csv(res / "small_panel_bracket.csv", index=False)
        print(bracket.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
        if len(bracket) > 2:
            rho, pv = sps.spearmanr(bracket["k"], bracket["median_gap"])
            print(f"  Spearman(k, gap) = {rho:+.3f}, p = {pv:.4f}  "
                  "(RAW-score alpha, so this is a lower bound on the residualised gap)")
            report["small_panel_spearman"] = [float(rho), float(pv)]
        report["small_panel_bracket"] = rows
    except Exception as exc:  # noqa: BLE001
        print(f"  FAILED: {exc}")
        failures.append(("A6", f"{type(exc).__name__}: {exc}"))

    # A stage that DIED must not look like a stage that ran. Each handler above
    # keeps swallowing its exception -- a control that genuinely cannot be
    # computed is a result -- but the outcome is now recorded in the artefact
    # and in the exit code. Before this, three of five stages could fail and the
    # script still returned 0 while writing a science_gaps.json whose nulls were
    # indistinguishable from "not applicable to this cohort".
    report["stage_failures"] = [{"stage": s, "error": e} for s, e in failures]
    (res / "science_gaps.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWritten to {res}/science_gaps.json")
    if failures:
        print(f"\n!! {len(failures)} STAGE(S) FAILED: "
              + ", ".join(s for s, _ in failures))
        for stage, err in failures:
            print(f"   {stage}: {err}")
        print("   science_gaps.json records these under 'stage_failures'.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
