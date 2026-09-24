#!/usr/bin/env python3
"""Merge a sharded partition sweep, and derive what the manuscript compares.

    # pan-cancer: shards from three passes -> one table, computed ONCE
    python scripts/27_partition_sweep.py merge --cohort pancancer \\
        --source default=SHARDS_A --source isi_only=SHARDS_B \\
        --source fast_isi_only=SHARDS_C --expect 0-23 \\
        --not-computed-kind 0=timeout ... \\
        --out results/partition_variance_pancancer_stable24

    # NSCLC ran unsharded: derive the same quantities beside its summary
    python scripts/27_partition_sweep.py derive \\
        results/partition_variance_nsclc_stable24/per_partition.csv \\
        --out results/partition_variance_nsclc_stable24/derived.json

    # the cross-cohort scaling the Results section states
    python scripts/27_partition_sweep.py compare \\
        --nsclc results/partition_variance_nsclc_stable24/derived.json \\
        --pancancer results/partition_variance_pancancer_stable24/summary.json \\
        --out results/partition_variance_stable24_scaling.json

WHY A SEPARATE SCRIPT
=====================
`12_partition_variance.py --shard` deliberately writes rows and no summary: a
standard deviation over one shard is meaningless, and summarising each shard
and averaging would be wrong. The pan-cancer sweep under the corrected ordering
arrived in three passes (default mode, array 129207; `--isi-only`, array
129667; fast path plus `--isi-only`, array 129874), because the 4-hour cap
killed some partitions and a LAPACK failure in two secondary analyses killed
others. This script pools them and computes the spread exactly once.

WHAT MAKES POOLING LEGITIMATE (checked, not assumed)
====================================================
The three passes compute the ISI by the same arithmetic: `--isi-only` removes
two analyses the ISI never reads (job 129666: bitwise equal on NSCLC), and the
shared-mask fast path is bitwise equal to the per-column reference on real
pan-cancer folds (job 129737). A seed that appears in more than one source is
a direct check of that claim, so this script REQUIRES such rows to be bitwise
equal (isi, lo, hi) and refuses to merge otherwise. The first source listed
wins; the duplicate is recorded, not discarded silently.

WHAT IS COMPUTED
================
The same quantities 12_partition_variance.py prints for an unsharded sweep --
mean, sd (ddof=1), bootstrap SE as the mean of (hi - lo) / (2 * 1.959963984540054),
combined SE, the reported and honest intervals, how much too narrow the reported
one is -- plus three that the manuscript needs and that script never stored:

  * the partition share of total variance, sd^2 / (sd^2 + SE_boot^2);
  * a chi-square 95% interval for the true sd at k - 1 degrees of freedom,
    [sd * sqrt((k-1) / q_0.975), sd * sqrt((k-1) / q_0.025)];
  * that interval propagated into the share (the share is monotone in sd).

A seed in the expected range with neither a row nor a refusal record is
`no_output` unless `--not-computed-kind SEED=KIND` names what happened (the
shard script cannot record a timeout: the scheduler kills it first).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
CRIT = 1.959963984540054
# pandas' default CSV float parser is NOT round-trip exact (measured 2026-09-16
# on both partition-variance tables: some values come back a ULP away), so
# every read here asks for the exact one.
EXACT = {"float_precision": "round_trip"}
ROW_COLUMNS = ["split_seed", "isi", "lo", "hi", "n_signatures", "runtime_s"]

# Cohort sizes that set the sqrt(n) scaling predictions. Patients: the frozen
# primary runs' analysed cohorts (n=944 and n=7,168). Sites: every tissue
# source site in each cohort, counted from the barcodes in each run's
# predictions.csv.gz -- the quantity site-to-fold assignment ranges over (see
# 14-SCIENCE-AUDIT.md, A5 correction of 2026-09-07, for why Control C's usable
# site counts are the wrong quantity here).
COHORT_SIZE = {"nsclc": {"patients": 944, "sites": 68},
               "pancancer": {"patients": 7168, "sites": 619}}


def _parse_range(text: str) -> list[int]:
    lo, _, hi = text.partition("-")
    return list(range(int(lo), int(hi or lo) + 1))


def _read_source(label: str, root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rows and refusal records from one source: a sweep dir or a shard tree."""
    files = []
    if (root / "per_partition.csv").is_file():
        files.append(root / "per_partition.csv")
    files += sorted(root.glob("seed_*/per_partition.csv"))
    refusals = []
    if (root / "refused_partitions.csv").is_file():
        refusals.append(root / "refused_partitions.csv")
    refusals += sorted(root.glob("seed_*/refused_partitions.csv"))
    if not files and not refusals:
        raise SystemExit(f"FATAL: source {label}={root} holds no per_partition.csv "
                         "and no refused_partitions.csv -- wrong path?")
    rows = [pd.read_csv(f, **EXACT).assign(mode=label, source=str(f)) for f in files]
    rows = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=ROW_COLUMNS + ["mode", "source"])
    missing = set(ROW_COLUMNS) - set(rows.columns)
    if missing:
        raise SystemExit(f"FATAL: {label} rows lack columns {sorted(missing)}")
    refs = [pd.read_csv(f, **EXACT).assign(mode=label, source=str(f)) for f in refusals]
    refs = pd.concat(refs, ignore_index=True) if refs else pd.DataFrame(
        columns=["split_seed", "kind", "mode", "source"])
    return rows, refs


def derive(rows: pd.DataFrame) -> dict:
    """The summary quantities, computed once over the whole table."""
    if rows["split_seed"].duplicated().any():
        raise ValueError("split_seed must be unique before deriving a spread")
    k = len(rows)
    if k < 3:
        raise ValueError(f"only {k} partitions; a spread needs at least 3")
    # The SAME operations, in the same order, as 12_partition_variance.py's
    # unsharded summary: pandas' std and mean differ from numpy's by a ULP, and a
    # derived number that disagrees with its source in the last digit is a
    # number the checker cannot match.
    isi = rows["isi"].to_numpy()
    boot_se = float(np.mean((rows["hi"] - rows["lo"]) / (2 * CRIT)))
    sd = float(np.std(isi, ddof=1))
    total_se = float(np.sqrt(boot_se**2 + sd**2))
    centre = float(np.mean(isi))
    df = k - 1
    sd_lo = sd * math.sqrt(df / stats.chi2.ppf(0.975, df))
    sd_hi = sd * math.sqrt(df / stats.chi2.ppf(0.025, df))

    def share(s: float) -> float:
        return s**2 / (s**2 + boot_se**2)

    return {
        "n_partitions_used": k,
        "split_seeds": [int(s) for s in sorted(rows["split_seed"])],
        "isi_mean": centre,
        "isi_min": float(isi.min()), "isi_max": float(isi.max()),
        "isi_sd_across_partitions": sd,
        "isi_sd_chi2_95ci": [sd_lo, sd_hi],
        "bootstrap_se": boot_se,
        "combined_se": total_se,
        "reported_ci": [centre - CRIT * boot_se, centre + CRIT * boot_se],
        "honest_ci": [centre - CRIT * total_se, centre + CRIT * total_se],
        "interval_understated_by": total_se / boot_se - 1,
        "partition_share_of_total_variance": share(sd),
        "partition_share_from_sd_chi2_95ci": [share(sd_lo), share(sd_hi)],
        "ci_still_excludes_zero": bool(centre - CRIT * total_se > 0),
    }


def _portable(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return f"<outside the repository>/{path.name}"


def merge(args) -> int:
    expected = _parse_range(args.expect)
    kinds = {}
    for spec in args.not_computed_kind:
        seed, _, kind = spec.partition("=")
        kinds[int(seed)] = kind
    sources = []
    for spec in args.source:
        label, _, path = spec.partition("=")
        if not path:
            raise SystemExit(f"FATAL: --source must be LABEL=DIR, got {spec!r}")
        sources.append((label, Path(path)))
    origins = {}
    for spec in args.origin:
        label, _, text = spec.partition("=")
        origins[label] = text
    unknown = set(origins) - {label for label, _ in sources}
    if unknown:
        raise SystemExit(f"FATAL: --origin names no --source: {sorted(unknown)}")

    all_rows, all_refs = [], []
    for label, root in sources:
        r, f = _read_source(label, root)
        all_rows.append(r)
        all_refs.append(f)
    rows = pd.concat(all_rows, ignore_index=True)
    refs = pd.concat(all_refs, ignore_index=True)

    duplicates = []
    for seed, grp in rows.groupby("split_seed", sort=True):
        if len(grp) < 2:
            continue
        first = grp.iloc[0]
        equal = all(
            all(float(o[c]) == float(first[c]) for c in ("isi", "lo", "hi"))
            for _, o in grp.iloc[1:].iterrows())
        duplicates.append({"split_seed": int(seed), "modes": list(grp["mode"]),
                           "isi": [repr(float(v)) for v in grp["isi"]],
                           "bitwise_equal": equal})
        if not equal:
            raise SystemExit(
                f"FATAL: split_seed={seed} differs between {list(grp['mode'])}: "
                f"{list(grp['isi'])}. Pooling rests on these being bitwise equal; "
                "do not merge until the difference is explained.")
    rows = rows.drop_duplicates("split_seed", keep="first")
    stray = sorted(set(rows["split_seed"]) - set(expected))
    if stray:
        raise SystemExit(f"FATAL: rows for seeds outside --expect: {stray}")

    computed = set(int(s) for s in rows["split_seed"])
    not_computed = []
    for seed in expected:
        if seed in computed:
            continue
        rec = refs[refs["split_seed"] == seed]
        # Records written before 2026-09-16 carry a `reason` and no `kind`.
        rec_kinds = (rec["kind"].astype(str) if "kind" in rec.columns
                     else pd.Series(["unrecorded"] * len(rec), index=rec.index))
        detail = "; ".join(f"{m}:{k}" for m, k in zip(rec["mode"], rec_kinds))
        recorded = rec_kinds.iloc[-1] if len(rec) else "no_output"
        not_computed.append({
            "split_seed": seed,
            "kind": kinds.get(seed, recorded),
            "refusal_records": detail})
    for seed in sorted(kinds):
        if seed in computed:
            raise SystemExit(f"FATAL: --not-computed-kind names seed {seed}, which has a row")

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    for name in ("per_partition.csv", "not_computed.csv", "summary.json"):
        if (out / name).exists() and not args.overwrite_new:
            raise SystemExit(f"FATAL: {out / name} exists; refusing to overwrite")
    table = rows.sort_values("split_seed")[ROW_COLUMNS + ["mode"]]
    table.to_csv(out / "per_partition.csv", index=False)
    pd.DataFrame(not_computed, columns=["split_seed", "kind", "refusal_records"]).to_csv(
        out / "not_computed.csv", index=False)

    summary = {"cohort": args.cohort, "n_partitions_requested": len(expected)}
    summary.update(derive(table))
    summary.update({
        "n_partitions_not_computed": len(not_computed),
        "not_computed": not_computed,
        "rows_by_mode": {m: int(n) for m, n in table["mode"].value_counts().sort_index().items()},
        "duplicate_seeds": duplicates,
        # A deposited summary must not carry this machine's paths: a source
        # inside the repository is recorded relative to it, one outside by its
        # last component only, and --origin says where the shards came from.
        "sources": {label: {"read_from": _portable(root), "origin": origins.get(label)}
                    for label, root in sources},
        "pooling_basis": (
            "Rows from different modes are pooled because each mode computes the ISI "
            "by the same arithmetic: --isi-only was bitwise equal on NSCLC (HPC4 job "
            "129666) and the shared-mask fast path was bitwise equal on real pan-cancer "
            "folds (HPC4 job 129737). Every seed present in more than one source is "
            "required to be bitwise equal; see duplicate_seeds."),
    })
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    _print(summary)
    print(f"\nWritten to {out}")
    return 0


def _print(s: dict) -> None:
    print(f"{s.get('cohort', '?')}: {s['n_partitions_used']} partitions "
          f"(seeds {s['split_seeds']})")
    print(f"  ISI mean {s['isi_mean']:.4f}, sd {s['isi_sd_across_partitions']:.6f} "
          f"[chi-square 95%: {s['isi_sd_chi2_95ci'][0]:.6f}, {s['isi_sd_chi2_95ci'][1]:.6f}]")
    print(f"  bootstrap SE {s['bootstrap_se']:.6f}, combined {s['combined_se']:.6f}, "
          f"reported interval {s['interval_understated_by']:.2%} too narrow")
    print(f"  partition share {s['partition_share_of_total_variance']:.2%} "
          f"[{s['partition_share_from_sd_chi2_95ci'][0]:.2%}, "
          f"{s['partition_share_from_sd_chi2_95ci'][1]:.2%}]")
    for rec in s.get("not_computed", []):
        print(f"  not computed: split_seed={rec['split_seed']} [{rec['kind']}]")


def derive_cmd(args) -> int:
    rows = pd.read_csv(args.per_partition, **EXACT)
    s = {"cohort": args.cohort, "source": str(args.per_partition)}
    s.update(derive(rows))
    if args.out.exists():
        raise SystemExit(f"FATAL: {args.out} exists; refusing to overwrite")
    args.out.write_text(json.dumps(s, indent=2) + "\n")
    _print(s)
    print(f"\nWritten to {args.out}")
    return 0


def compare(args) -> int:
    a = json.loads(args.nsclc.read_text())
    b = json.loads(args.pancancer.read_text())
    n, p = COHORT_SIZE["nsclc"], COHORT_SIZE["pancancer"]
    pred_patients = math.sqrt(p["patients"] / n["patients"])
    pred_sites = math.sqrt(p["sites"] / n["sites"])
    boot_ratio = a["bootstrap_se"] / b["bootstrap_se"]
    part_ratio = a["isi_sd_across_partitions"] / b["isi_sd_across_partitions"]
    # The ratio's own interval: the two sd estimates are independent chi-square
    # variates, so (sd_a/sd_b)^2 / (sigma_a/sigma_b)^2 ~ F(k_a-1, k_b-1).
    ka, kb = a["n_partitions_used"], b["n_partitions_used"]
    f_hi = stats.f.ppf(0.975, ka - 1, kb - 1)
    f_lo = stats.f.ppf(0.025, ka - 1, kb - 1)
    ratio_ci = [part_ratio / math.sqrt(f_hi), part_ratio / math.sqrt(f_lo)]
    out = {
        "nsclc_source": str(args.nsclc), "pancancer_source": str(args.pancancer),
        "k": {"nsclc": ka, "pancancer": kb},
        "cohort_size": COHORT_SIZE,
        "sqrt_patients_prediction": pred_patients,
        "sqrt_sites_prediction": pred_sites,
        "bootstrap_component_ratio": boot_ratio,
        "bootstrap_observed_over_sqrt_patients": boot_ratio / pred_patients,
        "partition_component_ratio": part_ratio,
        "partition_component_ratio_f_95ci": ratio_ci,
        "partition_observed_over_sqrt_patients": part_ratio / pred_patients,
        "partition_observed_over_sqrt_sites": part_ratio / pred_sites,
        "sqrt_patients_inside_ratio_ci": bool(ratio_ci[0] <= pred_patients <= ratio_ci[1]),
        "sqrt_sites_inside_ratio_ci": bool(ratio_ci[0] <= pred_sites <= ratio_ci[1]),
        "one_inside_ratio_ci": bool(ratio_ci[0] <= 1.0 <= ratio_ci[1]),
        "sd_chi2_95ci_overlap": bool(
            max(a["isi_sd_chi2_95ci"][0], b["isi_sd_chi2_95ci"][0])
            <= min(a["isi_sd_chi2_95ci"][1], b["isi_sd_chi2_95ci"][1])),
    }
    if args.out.exists():
        raise SystemExit(f"FATAL: {args.out} exists; refusing to overwrite")
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    for key in ("sqrt_patients_prediction", "sqrt_sites_prediction",
                "bootstrap_component_ratio", "partition_component_ratio"):
        print(f"  {key}: {out[key]:.4f}")
    print(f"  partition ratio 95% (F): [{ratio_ci[0]:.4f}, {ratio_ci[1]:.4f}]; "
          f"sqrt(patients) inside: {out['sqrt_patients_inside_ratio_ci']}, "
          f"sqrt(sites) inside: {out['sqrt_sites_inside_ratio_ci']}, "
          f"1 inside: {out['one_inside_ratio_ci']}")
    print(f"  sd chi-square intervals overlap: {out['sd_chi2_95ci_overlap']}")
    print(f"\nWritten to {args.out}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("merge", help="pool shard rows and compute the spread once")
    m.add_argument("--cohort", choices=sorted(COHORT_SIZE), required=True)
    m.add_argument("--source", action="append", required=True,
                   help="LABEL=DIR; DIR holds per_partition.csv and/or seed_*/ shards. "
                        "Listed in priority order: the first source wins a duplicate seed.")
    m.add_argument("--expect", required=True, help="the requested split_seed range, e.g. 0-23")
    m.add_argument("--not-computed-kind", action="append", default=[],
                   help="SEED=KIND for a seed that left no record (e.g. 4=timeout)")
    m.add_argument("--origin", action="append", default=[],
                   help="LABEL=TEXT: where that source's shards came from (job, "
                        "machine directory); recorded in summary.json")
    m.add_argument("--out", type=Path, required=True)
    m.add_argument("--overwrite-new", action="store_true",
                   help="allow replacing this script's own earlier output in --out")
    m.set_defaults(fn=merge)

    d = sub.add_parser("derive", help="derived quantities for an unsharded sweep")
    d.add_argument("per_partition", type=Path)
    d.add_argument("--cohort", choices=sorted(COHORT_SIZE), required=True)
    d.add_argument("--out", type=Path, required=True)
    d.set_defaults(fn=derive_cmd)

    c = sub.add_parser("compare", help="cross-cohort scaling of the two components")
    c.add_argument("--nsclc", type=Path, required=True)
    c.add_argument("--pancancer", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    c.set_defaults(fn=compare)

    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
