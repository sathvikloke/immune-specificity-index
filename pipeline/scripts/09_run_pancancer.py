#!/usr/bin/env python3
"""The ISI pan-TCGA. Same estimand as the NSCLC arm, 31 cancer types.

    python scripts/09_run_pancancer.py --outdir results/pancancer

WHY RUN THIS AT ALL, given NSCLC already answered the question
==============================================================
Three reasons, in order of how much they matter:

1. **Generalisation across 31 tumour types is a partial substitute for an
   external cohort.** It is not a true external validation — same consortium,
   same slide pipeline, same source sites — and must never be described as one.
   But "the ISI holds in 31 diseases" is a materially stronger claim than "it
   holds in lung", and it is available today rather than after a GPU rental.

2. **n crosses 1,000.** Cohort size n >= 1,000 appeared in 12.9% of AACR 2026
   computational ORAL abstracts against 7.5% of posters (landscape audit §4.3).
   The NSCLC arm is n=944 — just under, for no scientific reason.

3. **It tests the one thing NSCLC cannot.** With a single disease, "the image
   reads composition" and "the image reads immune biology" are hard to separate
   from tissue of origin. Pan-cancer, the axis is computed WITHIN cancer type, so
   tissue of origin is removed by construction and what remains is within-disease
   variation.

WHAT WILL LIKELY DIFFER FROM NSCLC, stated before running
=========================================================
- The global axis explains less within type than pooled, so residualisation
  should cost more here than the ~0 it cost in NSCLC.
- Signature predictability varies a lot by tumour type, so expect between-type
  heterogeneity; the pooled ISI is a mean over a heterogeneous family and should
  be read alongside the per-signature table, not instead of it.
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

INTERIM = ROOT / "data" / "interim"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"
TYPE_COL = "cancer type abbreviation"


def build():
    meta = pd.read_parquet(INTERIM / "meta.parquet")
    X = np.load(INTERIM / "X.npy")
    meta, X = data.collapse_to_patient(meta, X)
    cohort = meta.rename(columns={TYPE_COL: "cancer_type"})

    purity, _ = data.load_purity(ROOT / "data/raw/purity/TCGA_ABSOLUTE_purity.csv")
    cohort = cohort.merge(purity, on="patient_id", how="left")

    expr = pd.read_parquet(INTERIM / "expression_hugo.parquet")
    expr.index = pd.Index(expr.index).astype(str).str[:12]
    expr = expr.groupby(level=0).mean()

    keep = cohort["patient_id"].isin(expr.index).to_numpy()
    cohort, X = cohort.loc[keep].reset_index(drop=True), X[keep]
    expr = expr.reindex(cohort["patient_id"])

    sigs = signatures.SignatureSet.from_gmt(GMT)
    sigs = signatures.SignatureSet(
        name=sigs.name, sets={f"SIG_{k}": v for k, v in sigs.sets.items()}, source=str(GMT)
    )
    scored = signatures.score_mean_z(expr, sigs.filter_to(set(expr.columns)))
    scored.index = cohort.index
    for c in scored.columns:
        cohort[c] = scored[c].to_numpy()
    return cohort, X, expr, sigs, list(scored.columns)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "pancancer")
    ap.add_argument("--n-null", type=int, default=1000)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    cohort, X, expr, sigs, sig_cols = build()
    print(f"{len(cohort)} patients · {cohort['cancer_type'].nunique()} cancer types · "
          f"{cohort['tss'].nunique()} sites · {X.shape[1]}-dim · {len(sig_cols)} signatures",
          flush=True)

    axis = globalaxis.compute_global_axis(
        expr, method="pc1", cancer_type=cohort.set_index("patient_id")["cancer_type"]
    )
    print(f"  axis: {axis}", flush=True)

    surv = None
    try:
        surv = outcome.from_frame(cohort, endpoint="PFI")
        print(f"  survival: PFI, {surv.n_events} events", flush=True)
    except KeyError:
        print("  !! no PFI columns; outcome arm skipped", flush=True)

    cfg = experiment.AuditConfig(
        n_folds=args.folds, seed=args.seed, n_boot=args.n_boot,
        n_null_sets=args.n_null, run_combat=False,   # descriptive only, and slow at this n
        stage_col=None,
    )
    result = experiment.run_audit(
        cohort, X, sig_cols, cfg,
        expression=expr, signature_set=sigs, global_axis=axis, survival=surv,
        verbose=True,
    )
    print("\n" + result.headline())
    result.save(args.outdir)
    print(f"\nWritten to {args.outdir}  ({result.runtime_s / 60:.1f} min)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
