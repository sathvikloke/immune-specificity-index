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
ISI is needed. Measured on HPC4 (8 CPUs, one partition per job):
  NSCLC                              266-365 s    (job 129201)
  pan-cancer, default mode           7,399-14,135 s for the 11 of 24 that
                                     finished inside the 4-hour cap (129207)
  pan-cancer, --isi-only             8,297-12,270 s (4 shards, 129667)
  pan-cancer, --isi-only + the       5,360-6,630 s on six shards; 13,403 and
    shared-mask fast path            14,326 s for two that shared a node with
                                     four other shards (129874, 2026-09-16/17)
Node placement moves a partition's cost by more than 2x at fixed arithmetic;
peak memory is 22-23 GiB per pan-cancer partition in every mode (sacct).

WHAT A "REFUSED" PARTITION ACTUALLY WAS (measured 2026-09-16, HPC4 job 129294)
=============================================================================
Every pan-cancer partition this script has "refused" -- frozen split_seed=5 on
macOS, stable-sort seeds 6, 11 and 12 on HPC4 -- carried the message "SVD did
not converge". That is numpy.linalg.LinAlgError, which SUBCLASSES ValueError, so
the `except (RuntimeError, ValueError)` below caught a LAPACK convergence failure
in the same clause as the designed guard in splits.py ("degenerate site-disjoint
partition"). The failing call is the COVARIATE BASELINE -- RidgeCV on one-hot
site + type + purity, fitted under the RANDOM-patient split -- whose training
design has 29-34 all-zero dummy columns and rank 590-595 of 650. scipy's gesdd
fails on it deterministically; gesvd and numpy's gesdd both converge on the same
matrix. The ISI never reads that baseline. So those partitions were discarded
for a failure in a computation whose output this script throws away.

Seed 15 failed at a SECOND site with the same message: the variance
decomposition's statsmodels OLS (site dummies again), whose pinv calls numpy's
SVD. The decomposition regresses on the preserved-site out-of-fold predictions,
which is why it varies by partition. The ISI never reads it either.

`--isi-only` removes both computations. The refusal record now says which KIND
of failure it was, so a LinAlgError can no longer pass for geometry.
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


# run_audit reads `purity_col` in exactly two places: the covariate baselines and
# the gate on the variance decomposition (`if config.purity_col in
# frame.columns`). Emptying the baselines and pointing purity_col at a name no
# cohort carries switches both off without touching the library.
_NOT_A_COLUMN = "__isi_only_no_such_column__"
_ISI_ONLY = {"baselines": (), "purity_col": _NOT_A_COLUMN}


def _failure_kind(exc: BaseException) -> str:
    """Name the failure: only the first kind says anything about the partition."""
    if isinstance(exc, RuntimeError) and str(exc).startswith("degenerate site-disjoint partition"):
        return "degeneracy_guard"
    if isinstance(exc, np.linalg.LinAlgError):
        return "linalg_error"
    return f"other_{type(exc).__name__}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["nsclc", "pancancer"], default="nsclc")
    ap.add_argument("--partitions", type=int, default=6,
                    help="number of site-to-fold assignments; split_seed "
                         "START..START+N-1, where START is --start-seed")
    # Added 2026-09-16 so a long sweep can be split across Slurm jobs that each
    # fit the 4-hour cap. The "~1-2 h" that sized the first such array was an
    # inference from a different workload and was wrong: a pan-cancer partition
    # measured 2.1-3.6 h on HPC4, and five of sixteen exceeded 4 h.
    ap.add_argument("--start-seed", type=int, default=0,
                    help="first split_seed (default 0). Shards of one sweep "
                         "must not overlap: seeds are the identity of a "
                         "partition, not an index into this run.")
    ap.add_argument("--shard", action="store_true",
                    help="this run is a PIECE of a larger sweep: write "
                         "per_partition.csv and skip the across-partition "
                         "summary, which is meaningless on a shard. Without "
                         "it the exit code and output are exactly as before.")
    ap.add_argument("--isi-only", action="store_true",
                    help="skip the two secondary analyses the ISI never reads: "
                         "the covariate baselines and the variance decomposition. "
                         "Their SVDs are what failed on every 'refused' pan-cancer "
                         "partition. Off by default so every earlier invocation "
                         "is unchanged.")
    ap.add_argument("--n-null", type=int, default=1000)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()
    outdir = args.outdir or (ROOT / "results" / f"partition_variance_{args.cohort}")
    outdir.mkdir(parents=True, exist_ok=True)

    cohort, X, expr, sigs = (load_nsclc() if args.cohort == "nsclc" else load_pancancer())
    assert _NOT_A_COLUMN not in cohort.columns, "the --isi-only sentinel is a real column"
    sig_cols = [c for c in cohort.columns if c.startswith("SIG_")]
    if not sig_cols:
        scored = signatures.score_mean_z(expr, sigs.filter_to(set(expr.columns)))
        scored.index = cohort.index
        for c in scored.columns:
            cohort[c] = scored[c].to_numpy()
        sig_cols = list(scored.columns)

    print(f"{args.cohort}: {len(cohort)} patients, {cohort['tss'].nunique()} sites, "
          f"{len(sig_cols)} signatures, {args.partitions} partitions "
          f"(split_seed {args.start_seed}..{args.start_seed + args.partitions - 1})"
          f"{' [SHARD]' if args.shard else ''}"
          f"{' [ISI ONLY]' if args.isi_only else ''}", flush=True)

    axis = globalaxis.compute_global_axis(
        expr, method="pc1", cancer_type=cohort.set_index("patient_id")["cancer_type"])

    rows, skipped = [], []
    for sp in range(args.start_seed, args.start_seed + args.partitions):
        t0 = time.time()
        cfg = experiment.AuditConfig(
            n_folds=5, seed=0, split_seed=sp,          # <- ONLY the partition moves
            n_boot=args.n_boot, n_null_sets=args.n_null,
            run_combat=False, run_legacy_venet_null=False,
            stage_col=None if args.cohort == "pancancer" else "stage",
            **(_ISI_ONLY if args.isi_only else {}),
        )
        try:
            res = experiment.run_audit(
                cohort, X, sig_cols, cfg, expression=expr, signature_set=sigs,
                global_axis=axis, survival=None, verbose=False)
        except (RuntimeError, ValueError) as exc:
            # Record it rather than silently retrying. But say WHICH failure it
            # was: numpy.linalg.LinAlgError subclasses ValueError, so this clause
            # catches a LAPACK convergence failure as readily as the designed
            # degeneracy guard, and only the latter is a statement about the
            # partition's geometry. See the docstring.
            kind = _failure_kind(exc)
            print(f"  split_seed={sp}: PARTITION NOT COMPUTED [{kind}] ({exc})", flush=True)
            skipped.append({"split_seed": sp, "kind": kind,
                            "exception": type(exc).__name__, "reason": str(exc)[:200]})
            continue
        est = res.primary_excess
        if est is None or not np.isfinite(est.value):
            skipped.append({"split_seed": sp, "kind": "primary_not_computed",
                            "exception": "", "reason": "primary not computed"})
            continue
        rows.append({"split_seed": sp, "isi": est.value, "lo": est.lo, "hi": est.hi,
                     "n_signatures": est.n, "runtime_s": time.time() - t0})
        print(f"  split_seed={sp}: ISI {est.value:.4f} [{est.lo:.4f}, {est.hi:.4f}]"
              f"  ({time.time() - t0:.0f}s)", flush=True)

    # Persist BEFORE any early return: a shard's whole purpose is its rows, and
    # a sweep that threw away the data because one piece could not summarise it
    # would have to be re-run from the start.
    # The refusal record is written whether or not any partition succeeded: a
    # one-partition shard that fails has NOTHING else to show for itself, and
    # before 2026-09-16 it wrote no file at all.
    if rows:
        d = pd.DataFrame(rows)
        d.to_csv(outdir / "per_partition.csv", index=False)
    if skipped:
        pd.DataFrame(skipped).to_csv(outdir / "refused_partitions.csv", index=False)

    if args.shard:
        print(f"\nSHARD: wrote {len(rows)} partition(s) "
              f"(split_seed {args.start_seed}..{args.start_seed + args.partitions - 1}) "
              f"to {outdir / 'per_partition.csv'}. The across-partition summary is "
              "NOT computed here -- merge the shards and compute it once.")
        if skipped:
            print(f"SHARD: {len(skipped)} partition(s) not computed, recorded in "
                  f"{outdir / 'refused_partitions.csv'}")
        # Exit 0 when every requested seed is ACCOUNTED FOR -- computed or
        # recorded as not computed. A recorded failure is data; exiting 1 made
        # sacct report it as FAILED, indistinguishable from a crash (job 129207).
        return 0 if len(rows) + len(skipped) == args.partitions else 1

    if len(rows) < 3:
        print(f"\nOnly {len(rows)} usable partitions; a spread needs at least 3.")
        return 1

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
        print(f"\n  {len(skipped)} partition(s) not computed: "
              f"{[(s['split_seed'], s['kind']) for s in skipped]}")
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
