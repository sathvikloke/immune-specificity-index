#!/usr/bin/env python3
"""The ancestry arm, pan-TCGA. SUPPLEMENTARY — see aacr27.ancestry for the scope.

    python scripts/06_run_ancestry.py --outdir results/ancestry

Question: a model trained on a cohort that is 81% EUR, does it predict TME
signatures as well for everyone else? One model, preserved-site CV over the whole
pan-cancer cohort, then correlations computed WITHIN each ancestry group and
contrasted against EUR in Fisher z.

This is a non-inferiority analysis and is reported as a supplementary table. Every
contrast is printed next to its minimum detectable effect, because with AMR at
n=181 the honest statement will usually be "we could not detect a difference of at
least X" rather than "there is no difference".
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import ancestry as anc_mod  # noqa: E402
from aacr27 import data, models, signatures, splits, stats  # noqa: E402

INTERIM = ROOT / "data" / "interim"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"
TYPE_COL = "cancer type abbreviation"


def build_cohort():
    meta = pd.read_parquet(INTERIM / "meta.parquet")
    X = np.load(INTERIM / "X.npy")
    meta, X = data.collapse_to_patient(meta, X)

    ancestry, _ = data.load_ancestry(ROOT / "data/raw/ancestry/UCSF_Ancestry_Calls.csv")
    cohort = meta.merge(ancestry[["patient_id", "ancestry_broad"]], on="patient_id", how="left")
    cohort = cohort.rename(columns={TYPE_COL: "cancer_type"})

    keep = cohort["ancestry_broad"].isin(anc_mod.ANALYSIS_GROUPS).to_numpy()
    cohort, X = cohort.loc[keep].reset_index(drop=True), X[keep]

    expr = pd.read_parquet(INTERIM / "expression_hugo.parquet")
    # Expression is keyed on the 15-character SAMPLE barcode; the cohort carries
    # the 12-character patient id. Collapse expression to patient (mean over
    # samples) before joining, rather than assuming a one-to-one match.
    expr.index = pd.Index(expr.index).astype(str).str[:12]
    expr = expr.groupby(level=0).mean()

    shared = cohort["patient_id"].isin(expr.index).to_numpy()
    cohort, X = cohort.loc[shared].reset_index(drop=True), X[shared]
    expr = expr.reindex(cohort["patient_id"])

    sigs = signatures.SignatureSet.from_gmt(GMT)
    scored = signatures.score_mean_z(expr, sigs.filter_to(set(expr.columns)))
    scored.index = cohort.index
    for c in scored.columns:
        cohort[f"SIG_{c}"] = scored[c].to_numpy()

    return cohort, X, [f"SIG_{c}" for c in scored.columns]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "ancestry")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cohort, X, sig_cols = build_cohort()
    print(f"{len(cohort)} patients with embeddings + ancestry + expression · "
          f"{X.shape[1]}-dim · {len(sig_cols)} signatures")

    comp = anc_mod.composition_report(
        cohort["ancestry_broad"], cohort["cancer_type"], cohort["tss"]
    )
    print("\nPATIENT-LEVEL counts (NOT slide-level):")
    for g in anc_mod.ANALYSIS_GROUPS:
        print(f"  {g:6s} {comp['counts'].get(g, 0)}")
    # Calibrated against a permutation null, because the analytic correction is
    # unreliable on a table this sparse (4 groups x ~600 sites over ~6,600
    # patients is under 3 patients per cell).
    print("\nCONFOUNDING  (Cramer's V, permutation-calibrated)")
    sa = anc_mod.cramers_v_calibrated(cohort["ancestry_broad"], cohort["tss"], n_perm=200)
    ta = anc_mod.cramers_v_calibrated(cohort["ancestry_broad"], cohort["cancer_type"],
                                      n_perm=200)
    for label, st_ in (("site", sa), ("type", ta)):
        print(f"  ancestry x {label:5s} V={st_['v']:.3f} (uncorrected {st_['v_uncorrected']:.3f}) "
              f"vs permutation null median {st_['v_perm_median']:.3f} / p95 "
              f"{st_['v_perm_p95']:.3f}  ->  calibrated {st_['v_calibrated']:+.3f}, "
              f"p={st_['perm_p']:.4f}  ({st_['n_cols']} levels)")

    patients = cohort["patient_id"]
    sites = cohort["tss"].astype("string")
    split = splits.preserved_site_split(patients, sites, n_folds=args.folds, seed=args.seed)
    print(f"\npreserved-site CV: {args.folds} folds over {sites.nunique()} sites")

    Y = cohort[sig_cols].apply(pd.to_numeric, errors="coerce").to_numpy()
    print("fitting one multi-target ridge head over all signatures ...", flush=True)
    P = models.cross_val_predict_multi(X, Y, split, patients.to_numpy(), target_names=sig_cols)

    ident = anc_mod.site_identifiability(
        cohort["ancestry_broad"], cohort["tss"], cohort["cancer_type"]
    )
    print("\nIDENTIFIABILITY — can this contrast be separated from site at all?")
    for _, r in ident.iterrows():
        print(f"  {r['group']:6s} {r['n_sites_supporting']:3d} of {r['n_sites_total']:3d} sites "
              f"have >=10 of both ({r['frac_sites_supporting']:.1%}); "
              f"{r['n_type_site_cells_supporting']} (type,site) cells")

    per_group_rows, crude_rows, strat_rows, site_rows = [], [], [], []
    for j, sig in enumerate(sig_cols):
        pg = anc_mod.group_correlations(Y[:, j], P[:, j], cohort["ancestry_broad"])
        pg.insert(0, "signature", sig)
        per_group_rows.append(pg)

        cr = anc_mod.contrast_vs_reference(pg)
        if len(cr):
            cr.insert(0, "signature", sig)
            crude_rows.append(cr)

        st = anc_mod.stratified_contrast(
            Y[:, j], P[:, j], cohort["ancestry_broad"], cohort["cancer_type"]
        )
        if len(st):
            st.insert(0, "signature", sig)
            strat_rows.append(st)

        # THE DECISIVE SENSITIVITY: same contrast, but stratified by SITE rather
        # than cancer type. If a difference survives here it is not a site
        # effect. Cells are tiny, so min_cell_n drops to 10 and the result is
        # expected to be inconclusive — which is itself the answer.
        sw = anc_mod.stratified_contrast(
            Y[:, j], P[:, j], cohort["ancestry_broad"], cohort["tss"], min_cell_n=10
        )
        if len(sw):
            sw.insert(0, "signature", sig)
            site_rows.append(sw)

    per_group = pd.concat(per_group_rows, ignore_index=True)
    crude = pd.concat(crude_rows, ignore_index=True)
    strat = pd.concat(strat_rows, ignore_index=True)
    site_strat = pd.concat(site_rows, ignore_index=True) if site_rows else pd.DataFrame()

    # Multiplicity: BH-FDR within this family only. It is a FOURTH declared
    # family and must be counted as such in the paper's multiplicity statement.
    ok = strat["delta_z"].notna() & strat["se"].notna()
    p = np.full(len(strat), np.nan)
    from scipy import stats as sps
    p[ok.to_numpy()] = 2 * sps.norm.sf(
        np.abs(strat.loc[ok, "delta_z"] / strat.loc[ok, "se"]))
    strat["p"] = p
    rej = np.zeros(len(strat), dtype=bool)
    q = np.full(len(strat), np.nan)
    if ok.any():
        rej_ok, q_ok = stats.bh_fdr(p[ok.to_numpy()])
        rej[ok.to_numpy()], q[ok.to_numpy()] = rej_ok, q_ok
    strat["q_bh"], strat["significant_bh"] = q, rej

    per_group.to_csv(args.outdir / "per_group_correlations.csv", index=False)
    crude.to_csv(args.outdir / "contrasts_crude.csv", index=False)
    strat.to_csv(args.outdir / "contrasts_type_stratified.csv", index=False)
    comp["type_share"].to_csv(args.outdir / "cancer_type_composition.csv")
    ident.to_csv(args.outdir / "site_identifiability.csv", index=False)
    if len(site_strat):
        site_strat.to_csv(args.outdir / "contrasts_site_stratified.csv", index=False)

    print("\n" + "=" * 78)
    print("PER-GROUP PREDICTABILITY (median r across the 16 signatures)")
    print("=" * 78)
    med = per_group.groupby("group").agg(median_r=("r", "median"), n=("n", "max"))
    print(med.reindex(list(anc_mod.ANALYSIS_GROUPS)).to_string(
        float_format=lambda v: f"{v:.3f}"))

    print("\n" + "=" * 78)
    print("CONTRAST vs EUR, in Fisher z   (negative = worse than EUR)")
    print("=" * 78)
    summary = {}
    for g in ("AFR", "ASIAN", "AMR"):
        pc = anc_mod.pooled_over_signatures(crude, group=g)
        ps = anc_mod.pooled_over_signatures(strat, group=g)
        sub = strat[strat["group"] == g]
        mde = float(sub["mde_80pct"].median()) if len(sub) else float("nan")
        n_sig = int(sub["significant_bh"].sum())
        het = int(sub["heterogeneous"].sum())
        print(f"\n  {g}  (n={comp['counts'].get(g, 0)})")
        print(f"    crude               mean delta_z = {pc.value:+.4f} [{pc.lo:+.4f}, {pc.hi:+.4f}]")
        print(f"    TYPE-STRATIFIED     mean delta_z = {ps.value:+.4f} [{ps.lo:+.4f}, {ps.hi:+.4f}]"
              f"   ({int(sub['n_strata_used'].median())} of 31 types)")
        pw = (anc_mod.pooled_over_signatures(site_strat, group=g)
              if len(site_strat) else None)
        if pw is not None and np.isfinite(pw.value):
            print(f"    SITE-STRATIFIED     mean delta_z = {pw.value:+.4f} "
                  f"[{pw.lo:+.4f}, {pw.hi:+.4f}]   <- decisive sensitivity")
        print(f"    median MDE (80% power, per signature) = {mde:.4f}")
        print(f"    signatures with q<0.05: {n_sig}/16   heterogeneous strata (I2): {het}/16")
        row_id = ident[ident["group"] == g]
        n_supp = int(row_id["n_sites_supporting"].iloc[0]) if len(row_id) else 0
        if n_sig and n_supp < 10:
            verdict = (f"difference detected, BUT only {n_supp} sites carry both groups — "
                       "NOT separable from site")
        elif n_sig:
            verdict = "DIFFERENCE DETECTED"
        else:
            verdict = f"no difference detected; could not rule out |delta_z| < {mde:.3f}"
        print(f"    -> {verdict}")
        summary[g] = {
            "n": comp["counts"].get(g, 0),
            "crude_mean_delta_z": pc.value, "crude_lo": pc.lo, "crude_hi": pc.hi,
            "stratified_mean_delta_z": ps.value, "stratified_lo": ps.lo, "stratified_hi": ps.hi,
            "site_stratified_mean_delta_z": None if pw is None else pw.value,
            "site_stratified_lo": None if pw is None else pw.lo,
            "site_stratified_hi": None if pw is None else pw.hi,
            "n_sites_supporting": n_supp,
            "n_types_used": int(sub["n_strata_used"].median()) if len(sub) else 0,
            "median_mde_80pct": mde, "n_signatures_q005": n_sig, "n_heterogeneous": het,
            "verdict": verdict,
        }

    payload = {
        "n_patients": int(len(cohort)),
        "counts_patient_level": comp["counts"],
        "site_cramers_v": sa, "type_cramers_v": ta,
        "contrasts": summary,
        "scope": ("SUPPLEMENTARY non-inferiority analysis. Fourth multiplicity "
                  "family. Counts are patient-level, not slide-level."),
    }
    (args.outdir / "summary.json").write_text(json.dumps(payload, indent=2, default=str))
    print(f"\nWritten to {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
