#!/usr/bin/env python3
"""Run the audit on the NSCLC (LUAD + LUSC) discovery cohort.

This is the canonical runner for the primary analysis. It lives in the repo
rather than a scratch directory because the exact invocation IS part of the
method: the cohort definition, the SIG_ prefixing convention, and the axis
construction all have to be reproducible from the manuscript.

    python scripts/05_run_nsclc.py --outdir results/nsclc_v2

Inputs are the interim artefacts written by scripts/01_fetch_data.py:
    data/interim/cohort_nsclc.parquet   944 patients, TCGA-CDR merged
    data/interim/X_nsclc.npy            944 x 768 Prov-GigaPath slide embeddings
    data/interim/expr_nsclc.parquet     944 x 41,046 log2(TPM+1), HUGO symbols

NAMING CONVENTION, AND WHY IT MATTERS
-------------------------------------
Signature COLUMNS in the cohort frame are prefixed `SIG_`; gene-SET names in
the SignatureSet must carry the identical prefix, because
`_residualised_null_scores` matches them by exact string equality. A mismatch
there does not crash — it silently skips every signature and returns an empty
primary. The prefix is applied to both sides here, in one place.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import experiment, globalaxis, outcome, signatures  # noqa: E402

INTERIM = ROOT / "data" / "interim"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"


def load_cohort() -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame, signatures.SignatureSet]:
    cohort = pd.read_parquet(INTERIM / "cohort_nsclc.parquet")
    X = np.load(INTERIM / "X_nsclc.npy")
    expr = pd.read_parquet(INTERIM / "expr_nsclc.parquet")

    sigs = signatures.SignatureSet.from_gmt(GMT)
    # Prefix the gene-SET names to match the SIG_ prefixed score COLUMNS.
    sigs = signatures.SignatureSet(
        name=sigs.name,
        sets={f"SIG_{k}": v for k, v in sigs.sets.items()},
        source=str(GMT),
    )

    if len(cohort) != len(X):
        raise ValueError(f"cohort {len(cohort)} rows vs X {len(X)} rows")
    if not expr.index.equals(pd.Index(cohort["patient_id"])):
        # Reindex rather than assume; the audit depends on row alignment.
        expr = expr.reindex(cohort["patient_id"])
    return cohort, X, expr, sigs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "nsclc_v2")
    ap.add_argument("--n-null", type=int, default=1000)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-outcome", action="store_true")
    ap.add_argument("--no-combat", action="store_true")
    args = ap.parse_args()

    cohort, X, expr, sigs = load_cohort()
    sig_cols = [c for c in cohort.columns if c.startswith("SIG_")]
    print(f"{len(cohort)} patients · {X.shape[1]}-dim embeddings · "
          f"{expr.shape[1]} genes · {len(sig_cols)} signatures")

    axis = globalaxis.compute_global_axis(
        expr, method="pc1", cancer_type=cohort.set_index("patient_id")["cancer_type"]
    )
    print(f"  axis: {axis}")

    surv = None
    if not args.no_outcome:
        surv = outcome.from_frame(cohort, endpoint="PFI")
        print(f"  survival: PFI, {surv.n_events} events")

    cfg = experiment.AuditConfig(
        n_folds=args.folds, seed=args.seed, n_boot=args.n_boot,
        n_null_sets=args.n_null, run_combat=not args.no_combat,
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
