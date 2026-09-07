#!/usr/bin/env python3
"""Run the audit.

Two modes:

    python scripts/02_run_audit.py --demo
        Synthetic cohort with a planted site confound. Proves the pipeline works
        end to end in about a minute, no data required.

    python scripts/02_run_audit.py --real --outdir results/run01
        The real thing, once scripts/01_fetch_data.py has been run and the
        column-name guesses in data.py have been verified.

PRE-REGISTER BEFORE THE REAL RUN
--------------------------------
Write the primary endpoint, the analysis plan and the multiplicity plan down and
timestamp them (an OSF registration or a signed git tag both work) BEFORE the
first real run. This costs nothing and pre-empts the most common reviewer
objection. The endpoint is recorded in AuditConfig and written to config.json in
the output directory, so the run itself is self-documenting.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import data, experiment, globalaxis, outcome, signatures  # noqa: E402


def build_demo_cohort(n_patients: int = 800, n_sites: int = 24, dim: int = 96,
                      seed: int = 0):
    """Synthetic cohort: three signatures with different amounts of site artefact."""
    rng = np.random.default_rng(seed)
    site_ids = [f"{i:02d}" for i in range(n_sites)]
    site_of = rng.choice(site_ids, size=n_patients)
    patients = [f"TCGA-{site_of[i]}-{i:04d}" for i in range(n_patients)]

    biology = rng.normal(size=n_patients)
    site_fingerprint = {s: rng.normal(size=dim) for s in site_ids}
    site_signal = np.vstack([site_fingerprint[s] for s in site_of])
    site_scalar = np.array([rng.normal() for _ in site_ids])
    site_map = dict(zip(site_ids, site_scalar))
    site_vec = np.array([site_map[s] for s in site_of])
    purity = np.clip(0.55 + 0.15 * rng.normal(size=n_patients), 0.05, 0.99)

    load_bio, load_pur = rng.normal(size=dim), rng.normal(size=dim)
    X = (
        np.outer(biology, load_bio)
        + site_signal
        + np.outer(purity, load_pur)
        + rng.normal(size=(n_patients, dim))
    )

    frame = pd.DataFrame({
        "patient_id": patients,
        "tss": site_of,
        "cancer_type": rng.choice(["LUAD", "LUSC", "BRCA", "COAD"], size=n_patients),
        "purity": purity,
        "stage": rng.choice(["I", "II", "III", "IV"], size=n_patients),
        # Mostly real biology.
        "SIG_immune_clean": 1.5 * biology + 0.2 * purity + rng.normal(scale=.5, size=n_patients),
        # Mostly site artefact — should collapse under preserved-site CV.
        "SIG_immune_siteconf": 0.3 * biology + 2.5 * site_vec + rng.normal(scale=.5, size=n_patients),
        # Mostly purity — should be absorbed by the purity baseline.
        "SIG_immune_purity": 0.2 * biology + 4.0 * purity + rng.normal(scale=.4, size=n_patients),
    })
    return frame, X


def run_demo(outdir: Path, n_boot: int) -> int:
    print("Synthetic demo — three signatures with known, different confounding.\n")
    frame, X = build_demo_cohort()
    sig_cols = [c for c in frame.columns if c.startswith("SIG_")]

    cfg = experiment.AuditConfig(
        n_folds=5, n_boot=n_boot, min_patients_per_signature=100, run_combat=True
    )
    result = experiment.run_audit(frame, X, sig_cols, cfg, verbose=True)

    print("\n" + result.headline())
    print("\nPer-signature detail:")
    cols = ["signature", "n", "r_random_patient", "r_preserved_site", "delta_r",
            "r_site_only", "r_purity_only", "r_covariates"]
    print(result.per_signature[[c for c in cols if c in result.per_signature.columns]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    if result.decomposition is not None and len(result.decomposition):
        print("\nPurity decomposition (incremental R^2 of image over covariates):")
        dcols = ["signature", "incremental_image_r2", "r_image", "r_purity", "r_partial"]
        print(result.decomposition[[c for c in dcols if c in result.decomposition.columns]]
              .to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    result.save(outdir)
    print(f"\nWritten to {outdir}")
    print("\nExpected pattern if the pipeline is behaving:")
    print("  SIG_immune_siteconf  largest delta_r, and r_site_only close to r_preserved_site")
    print("  SIG_immune_purity    r_purity_only high, low incremental image R^2")
    print("  SIG_immune_clean     small delta_r, high r under both schemes")
    return 0


def run_real(outdir: Path, args) -> int:
    raw = ROOT / "data" / "raw"
    emb_files = list((raw / "provgigapath").rglob("*.parquet"))
    if not emb_files:
        print("No embeddings found. Run: python scripts/01_fetch_data.py --embeddings")
        return 1

    print(f"Loading embeddings from {emb_files[0]}")
    meta, X, prov_emb = data.load_embeddings(emb_files[0])
    print(f"  {len(meta)} slides, embedding dim {X.shape[1]}")
    for note in prov_emb.notes:
        print(f"  note: {note}")

    meta, X = data.collapse_to_patient(meta, X)
    print(f"  collapsed to {len(meta)} patients")

    provenance = [prov_emb]
    ancestry = purity = sigs = None

    anc_files = list((raw / "ancestry").glob("*"))
    anc_files = [f for f in anc_files if f.suffix in {".xlsx", ".xls", ".csv"}]
    if anc_files:
        ancestry, prov_anc = data.load_ancestry(anc_files[0])
        provenance.append(prov_anc)
        print(f"  ancestry: {len(ancestry)} patients")
        for note in prov_anc.notes:
            print(f"    note: {note}")
    else:
        print("  !! no ancestry file — the ancestry arm will be skipped")

    pur_files = list((raw / "purity").glob("*.csv"))
    if pur_files:
        purity = pd.read_csv(pur_files[0])
        print(f"  purity: {len(purity)} rows")
    else:
        print("  !! no purity file — the purity decomposition will be skipped")

    sig_files = list((raw / "signatures").glob("*"))
    if sig_files:
        path = sig_files[0]
        sigs = (signatures.SignatureSet.from_gmt(path) if path.suffix == ".gmt"
                else signatures.SignatureSet.from_json(path))
        print(f"  signatures: {len(sigs)} sets")
    else:
        print("  !! no signature file — using the CYT smoke-test signature only")
        sigs = signatures.example_signature_set()

    cohort = data.assemble_cohort(meta, ancestry=ancestry, purity=purity)
    print("\n" + data.coverage_report(cohort))

    # Score the signatures rather than requiring the caller to have done it.
    # The previous version built `sigs`, never used it, and always exited here.
    expr = axis = surv = None
    expr_files = [f for f in (raw / "expression").glob("*")
                  if f.suffix in {".parquet", ".gz", ".tsv", ".csv"}]
    if expr_files:
        expr, prov_expr = data.load_expression(expr_files[0])
        provenance.append(prov_expr)
        print(f"  expression: {expr.shape[0]} samples x {expr.shape[1]} genes")

        scored = signatures.score_mean_z(expr, sigs.filter_to(set(expr.columns)))
        scored.columns = [f"SIG_{c}" for c in scored.columns]
        cohort = cohort.merge(
            scored, left_on="patient_id", right_index=True, how="left"
        )
        print(f"  scored {scored.shape[1]} signatures")

        axis = globalaxis.compute_global_axis(
            expr, method="pc1",
            cancer_type=cohort.set_index("patient_id")["cancer_type"]
            if "cancer_type" in cohort.columns else None,
        )
        print(f"  global axis: {axis}")
    else:
        print("  !! no expression matrix — the PRIMARY endpoint cannot be computed.")
        print("     Fetch it: python scripts/01_fetch_data.py --expression")

    try:
        surv = outcome.from_frame(cohort, endpoint="PFI")
        print(f"  survival: PFI, {surv.n_events} events")
    except KeyError:
        print("  !! no PFI columns — the outcome arm will be skipped")

    sig_cols = [c for c in cohort.columns if c.startswith("SIG_")]
    if not sig_cols:
        print("\nNo signature score columns and no expression matrix to build them from.")
        data.write_provenance(provenance, outdir / "provenance.json")
        return 1

    cfg = experiment.AuditConfig(n_folds=args.folds, n_boot=args.n_boot)
    result = experiment.run_audit(
        cohort, X, sig_cols, cfg,
        expression=expr, signature_set=sigs, global_axis=axis, survival=surv,
        verbose=True,
    )
    print("\n" + result.headline())
    result.save(outdir)
    data.write_provenance(provenance, outdir / "provenance.json")
    print(f"\nWritten to {outdir}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--demo", action="store_true", help="synthetic smoke test")
    mode.add_argument("--real", action="store_true", help="run on fetched data")
    parser.add_argument("--outdir", type=Path, default=None)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--n-boot", type=int, default=2000)
    args = parser.parse_args()

    outdir = args.outdir or (ROOT / "results" / ("demo" if args.demo else "run01"))
    outdir.mkdir(parents=True, exist_ok=True)

    return run_demo(outdir, args.n_boot) if args.demo else run_real(outdir, args)


if __name__ == "__main__":
    sys.exit(main())
