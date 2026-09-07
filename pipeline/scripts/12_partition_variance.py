#!/usr/bin/env python3
"""How much of the ISI's uncertainty is missing because we fixed one partition?

    python scripts/12_partition_variance.py --cohort nsclc --partitions 6

THE PROBLEM THIS MEASURES
=========================
Preserved-site cross-validation assigns whole tissue source sites to folds. There
are astronomically many valid assignments, and the pipeline uses exactly one — the
deterministic greedy solution at `seed=0`. The reported interval is a
patient-clustered bootstrap *conditional on that one partition*, so it answers
"how much would this move with different patients" and is silent on "how much
would this move with a different fold assignment".

`splits.repeated_preserved_site_splits` already recorded the scale of the omission
on NSCLC: across 15 valid partitions, Delta_r had sd 0.019 against a bootstrap
half-width of 0.048 — roughly **27% of the honest interval missing**. That number
was measured, written into a docstring, and never applied, because
`report_partition_variance` defaults to False.

WHAT THIS SCRIPT DOES
=====================
Re-runs the PRIMARY endpoint under several site-to-fold assignments, varying
`split_seed` ONLY. The null gene sets, the bootstrap draws and every other random
choice are held fixed by `seed`, so the spread that comes out is attributable to
the partition and nothing else. Conflating them — by varying `seed` — would have
been the easy mistake.

It then reports the combined interval:

    SE_total = sqrt(SE_bootstrap^2 + sd_partition^2)

which treats partition choice and patient sampling as independent sources, the
standard two-component variance decomposition. That is an approximation: the two
are not perfectly independent, because a partition that happens to isolate an
unusual site also changes which patients dominate a fold. It is nonetheless much
closer to honest than ignoring one source entirely.

COST
====
The legacy Venet null, ComBat and the outcome arm are all switched off — only the
ISI is needed. NSCLC runs in a few minutes per partition; pan-cancer is roughly
an order of magnitude more and is worth doing once the NSCLC magnitude is known.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import data, experiment, globalaxis, signatures  # noqa: E402

INTERIM = ROOT / "data" / "interim"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"


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


def load_pancancer():
    meta = pd.read_parquet(INTERIM / "meta.parquet")
    X = np.load(INTERIM / "X.npy")
    meta, X = data.collapse_to_patient(meta, X)
    cohort = meta.rename(columns={"cancer type abbreviation": "cancer_type"})
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
        name=sigs.name, sets={f"SIG_{k}": v for k, v in sigs.sets.items()}, source=str(GMT))
    return cohort, X, expr, sigs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["nsclc", "pancancer"], default="nsclc")
    ap.add_argument("--partitions", type=int, default=6,
                    help="number of site-to-fold assignments; split_seed 0..N-1")
    ap.add_argument("--n-null", type=int, default=1000)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()
    outdir = args.outdir or (ROOT / "results" / f"partition_variance_{args.cohort}")
    outdir.mkdir(parents=True, exist_ok=True)

    cohort, X, expr, sigs = (load_nsclc() if args.cohort == "nsclc" else load_pancancer())
    sig_cols = [c for c in cohort.columns if c.startswith("SIG_")]
    if not sig_cols:
        scored = signatures.score_mean_z(expr, sigs.filter_to(set(expr.columns)))
        scored.index = cohort.index
        for c in scored.columns:
            cohort[c] = scored[c].to_numpy()
        sig_cols = list(scored.columns)

    print(f"{args.cohort}: {len(cohort)} patients, {cohort['tss'].nunique()} sites, "
          f"{len(sig_cols)} signatures, {args.partitions} partitions", flush=True)

    axis = globalaxis.compute_global_axis(
        expr, method="pc1", cancer_type=cohort.set_index("patient_id")["cancer_type"])

    rows, skipped = [], []
    for sp in range(args.partitions):
        t0 = time.time()
        cfg = experiment.AuditConfig(
            n_folds=5, seed=0, split_seed=sp,          # <- ONLY the partition moves
            n_boot=args.n_boot, n_null_sets=args.n_null,
            run_combat=False, run_legacy_venet_null=False,
            stage_col=None if args.cohort == "pancancer" else "stage",
        )
        try:
            res = experiment.run_audit(
                cohort, X, sig_cols, cfg, expression=expr, signature_set=sigs,
                global_axis=axis, survival=None, verbose=False)
        except (RuntimeError, ValueError) as exc:
            # A degenerate site-to-fold assignment is refused by design. Record
            # it rather than silently retrying — how often the geometry fails is
            # itself informative.
            print(f"  split_seed={sp}: PARTITION REFUSED ({exc})", flush=True)
            skipped.append({"split_seed": sp, "reason": str(exc)[:200]})
            continue
        est = res.primary_excess
        if est is None or not np.isfinite(est.value):
            skipped.append({"split_seed": sp, "reason": "primary not computed"})
            continue
        rows.append({"split_seed": sp, "isi": est.value, "lo": est.lo, "hi": est.hi,
                     "n_signatures": est.n, "runtime_s": time.time() - t0})
        print(f"  split_seed={sp}: ISI {est.value:.4f} [{est.lo:.4f}, {est.hi:.4f}]"
              f"  ({time.time() - t0:.0f}s)", flush=True)

    if len(rows) < 3:
        print(f"\nOnly {len(rows)} usable partitions; a spread needs at least 3.")
        return 1

    d = pd.DataFrame(rows)
    d.to_csv(outdir / "per_partition.csv", index=False)

    crit = 1.959963984540054
    isi = d["isi"].to_numpy()
    boot_se = float(np.mean((d["hi"] - d["lo"]) / (2 * crit)))
    part_sd = float(np.std(isi, ddof=1))
    total_se = float(np.sqrt(boot_se**2 + part_sd**2))
    centre = float(np.mean(isi))

    print("\n" + "=" * 66)
    print(f"PARTITION VARIANCE — {args.cohort}, {len(d)} usable partitions")
    print("=" * 66)
    print(f"  ISI across partitions: mean {centre:.4f}, sd {part_sd:.4f}, "
          f"range {isi.min():.4f}-{isi.max():.4f}")
    print(f"  bootstrap SE (mean over partitions): {boot_se:.4f}")
    print(f"  partition sd:                       {part_sd:.4f}")
    print(f"  combined SE:                        {total_se:.4f}")
    print(f"\n  reported (one partition):  {centre:.4f} "
          f"[{centre - crit * boot_se:.4f}, {centre + crit * boot_se:.4f}]"
          f"  width {2 * crit * boot_se:.4f}")
    print(f"  HONEST (both sources):     {centre:.4f} "
          f"[{centre - crit * total_se:.4f}, {centre + crit * total_se:.4f}]"
          f"  width {2 * crit * total_se:.4f}")
    inflation = total_se / boot_se - 1
    print(f"\n  the reported interval is {inflation:.1%} too narrow")
    print(f"  partition variance is {part_sd**2 / total_se**2:.1%} of total variance")
    if skipped:
        print(f"\n  {len(skipped)} partition(s) refused by the degeneracy guard: "
              f"{[s['split_seed'] for s in skipped]}")
    still_positive = centre - crit * total_se > 0
    print(f"\n  CI still excludes zero: {still_positive}")

    (outdir / "summary.json").write_text(json.dumps({
        "cohort": args.cohort, "n_partitions_used": len(d),
        "n_partitions_refused": len(skipped), "refused": skipped,
        "isi_mean": centre, "isi_sd_across_partitions": part_sd,
        "bootstrap_se": boot_se, "combined_se": total_se,
        "reported_ci": [centre - crit * boot_se, centre + crit * boot_se],
        "honest_ci": [centre - crit * total_se, centre + crit * total_se],
        "interval_understated_by": inflation,
        "ci_still_excludes_zero": bool(still_positive),
    }, indent=2))
    print(f"\nWritten to {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
