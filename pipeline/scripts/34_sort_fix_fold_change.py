#!/usr/bin/env python3
"""How much of each cohort's partition did the A9 sort fix change? (ledger F1.1, session 50)

    python scripts/34_sort_fix_fold_change.py     # writes results/sort_fix_fold_change.json

WHY
===
Limitation 8 used to explain why the corrected ordering moved the pan-cancer
index by fewer partition standard deviations than the NSCLC index with an
ARGUMENT: pan-TCGA assigns 619 sites rather than 68, "so any single tie
resolved differently perturbs a correspondingly smaller share of the
partition". Nothing had measured that share. Both cohorts' frozen runs and
their corrected-ordering re-runs store the preserved-site fold of every
patient in `predictions.csv.gz`, so the share can be read off directly.

It turned out to point the other way: the fix moved MORE of the pan-cancer
partition, not less. This script records the measurement and the yardstick
ratios Limitation 8 quotes, so the prose has an authority to be checked against.

WHAT IT READS (all committed or stored artefacts; nothing is recomputed)
  results/{nsclc_v3,nsclc_v3_stablesort,pancancer_v3,pancancer_v3_stablesort}/
      predictions.csv.gz   (scheme == "preserved_site": one fold per patient)
      summary.json         (primary_excess)
  results/partition_variance_{nsclc,pancancer}{,_stable24}/summary.json
  results/scorer_sensitivity_hpc4/summary.json   (the Linux shipped-sort ISI)
  results/split_determinism_darwin_arm64.json    (the 68 / 51 site counts)

The site is the TCGA tissue-source-site code, characters 5-6 of the patient
barcode, which is how `tss` is defined; the NSCLC site and tied-site counts are
asserted against the split-determinism record so a different derivation fails
loudly rather than producing plausible numbers.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = RES / "sort_fix_fold_change.json"

COHORTS = {
    "nsclc": ("nsclc_v3", "nsclc_v3_stablesort", "partition_variance_nsclc",
              "partition_variance_nsclc_stable24"),
    "pancancer": ("pancancer_v3", "pancancer_v3_stablesort",
                  "partition_variance_pancancer", "partition_variance_pancancer_stable24"),
}


def _folds(run: str) -> pd.Series:
    d = pd.read_csv(RES / run / "predictions.csv.gz", usecols=["patient", "fold", "scheme"])
    d = d[d["scheme"] == "preserved_site"]
    per = d.groupby("patient")["fold"].nunique()
    if (per != 1).any():
        raise SystemExit(f"{run}: a patient carries more than one preserved-site fold")
    return d.drop_duplicates("patient").set_index("patient")["fold"].sort_index()


def _isi(run: str) -> float:
    return float(json.loads((RES / run / "summary.json").read_text())["primary_excess"]["value"])


def _sd(run: str) -> float:
    return float(json.loads((RES / run / "summary.json").read_text())["isi_sd_across_partitions"])


def cohort(frozen: str, fixed: str, pv_small: str, pv24: str) -> dict:
    fa, fb = _folds(frozen), _folds(fixed)
    if not fa.index.equals(fb.index):
        raise SystemExit(f"{frozen} and {fixed} do not hold the same patients")
    site = pd.Series(fa.index.str.slice(5, 7), index=fa.index)
    per_site = site.value_counts()
    size_counts = per_site.value_counts()
    tied_sizes = size_counts[size_counts > 1].index
    tied_sites = per_site.index[per_site.isin(tied_sizes)]
    in_tied = site.isin(tied_sites)
    changed = fa != fb
    n = len(fa)
    isi_frozen, isi_fixed = _isi(frozen), _isi(fixed)
    move = abs(isi_fixed - isi_frozen)
    sd_small, sd24 = _sd(pv_small), _sd(pv24)
    # The fix's move against what a change of partition actually does: every
    # pairwise |ISI difference| among the 24 stable-sort partitions.
    per = pd.read_csv(RES / pv24 / "per_partition.csv").drop_duplicates("split_seed")
    if len(per) != 24:
        raise SystemExit(f"{pv24}: expected 24 partitions, found {len(per)}")
    v = per["isi"].to_numpy()
    diffs = [abs(a - b) for a, b in itertools.combinations(v, 2)]
    return {
        "frozen_run": frozen,
        "corrected_run": fixed,
        "n_patients": n,
        "n_sites": int(len(per_site)),
        "n_sites_tied": int(len(tied_sites)),
        "n_patients_in_tied_sites": int(in_tied.sum()),
        "pct_patients_in_tied_sites": 100.0 * float(in_tied.mean()),
        "n_patients_fold_changed": int(changed.sum()),
        "pct_patients_fold_changed": 100.0 * float(changed.mean()),
        "n_changed_in_untied_sites": int((changed & ~in_tied).sum()),
        "fold_sizes_frozen": [int(v) for v in fa.value_counts().sort_index()],
        "fold_sizes_corrected": [int(v) for v in fb.value_counts().sort_index()],
        "adjusted_rand_index": float(adjusted_rand_score(fa, fb)),
        "isi_frozen": isi_frozen,
        "isi_corrected": isi_fixed,
        "isi_move": move,
        "partition_sd_small_sweep": sd_small,
        "partition_sd_24": sd24,
        "move_in_small_sweep_sds": move / sd_small,
        "move_in_24_partition_sds": move / sd24,
        "n_partition_pairs_24": len(diffs),
        "n_pairs_differing_by_at_least_the_move": int(sum(d >= move for d in diffs)),
        "median_pairwise_abs_difference_24": float(pd.Series(diffs).median()),
        "range_across_24": float(v.max() - v.min()),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    out = {name: cohort(*runs) for name, runs in COHORTS.items()}

    # Guards: the derivation must reproduce what is already on record.
    sd_rec = json.loads((RES / "split_determinism_darwin_arm64.json").read_text())
    ns = out["nsclc"]
    if (ns["n_sites"], ns["n_sites_tied"]) != (sd_rec["n_sites"], sd_rec["n_sites_tied"]):
        raise SystemExit(f"site derivation disagrees with split_determinism: "
                         f"{ns['n_sites']}/{ns['n_sites_tied']} vs "
                         f"{sd_rec['n_sites']}/{sd_rec['n_sites_tied']}")
    if ns["n_patients_fold_changed"] != sd_rec["n_patients_fold_changed"]:
        raise SystemExit("NSCLC fold-change count disagrees with split_determinism's macOS record")
    for name, c in out.items():
        if c["n_changed_in_untied_sites"] != 0:
            raise SystemExit(f"{name}: a patient in an untied site changed fold -- the fix "
                             "is supposed to reorder tied sites only")
        if c["fold_sizes_frozen"] != c["fold_sizes_corrected"]:
            raise SystemExit(f"{name}: fold sizes differ between the two partitions")

    linux = json.loads((RES / "scorer_sensitivity_hpc4" / "summary.json").read_text())
    linux_isi = float(linux["arms"]["mean_z__disattenuated"]["isi"])
    gap = ns["isi_frozen"] - linux_isi
    pc = out["pancancer"]
    out["nsclc_platform_gap"] = {
        "source": "frozen macOS NSCLC ISI minus the Linux shipped-sort ISI "
                  "(scorer_sensitivity_hpc4, mean_z disattenuated arm)",
        "linux_shipped_isi": linux_isi,
        "gap": gap,
        "gap_in_small_sweep_sds": gap / ns["partition_sd_small_sweep"],
        "gap_in_24_partition_sds": gap / ns["partition_sd_24"],
    }
    out["contrast"] = {
        "nsclc_over_pancancer_in_small_sweep_sds":
            ns["move_in_small_sweep_sds"] / pc["move_in_small_sweep_sds"],
        "nsclc_over_pancancer_in_24_partition_sds":
            ns["move_in_24_partition_sds"] / pc["move_in_24_partition_sds"],
        "raw_move_ratio": ns["isi_move"] / pc["isi_move"],
    }
    out["reading"] = (
        "The fix reorders equally sized sites only: no patient in an untied site "
        "changed fold and every fold kept its size, in both cohorts. Pan-TCGA has "
        "more of its patients in tied sites, and the fix moved more of its partition "
        "(ARI near 0 means a nearly unrelated partition), yet the index moved by fewer "
        "partition sds. The share of the partition changed does not explain the contrast."
    )
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
    print(f"\nwrote {args.out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
