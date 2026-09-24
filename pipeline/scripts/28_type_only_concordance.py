#!/usr/bin/env python3
"""How prognostic is cancer type ALONE on pooled pan-TCGA PFI?

    python scripts/28_type_only_concordance.py        # writes results/type_only_concordance.json

WHY THIS EXISTS (ledger E13, 2026-09-16)
========================================
The outcome arm is stratified by cancer type because, pooled across 31 types,
cancer type alone reaches C = 0.676 on progression-free interval. That number
was measured once, on 2026-08-19 (commit 39664aa), and then lived only as a
LITERAL: in a comment and in the note `experiment.run_audit` writes into every
summary.json, NSCLC's included. A checker reading that note was checking the
manuscript against a string the code prints regardless of the data. This
script recomputes it from the cohort and stores the result, so the number has
a source.

WHAT "CANCER TYPE ALONE" MEANS HERE
===================================
Harrell's C (`outcome.concordance_index`, unstratified) of a score that is a
function of cancer type only. The recorded 0.676 is reproduced by the type's
PFI event proportion used as the risk score; two other type-only scores are
reported beside it so the number is not read as more exact than it is. All
three are in-sample: the score is fitted on the patients it is evaluated on.

The cohort is the pan-cancer analysis cohort: slides collapsed to patients,
restricted to patients with any TOIL sample (7,168). Needs data/interim.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import data, outcome  # noqa: E402

INTERIM = ROOT / "data" / "interim"
OUT = ROOT / "results" / "type_only_concordance.json"


def main() -> int:
    meta = pd.read_parquet(INTERIM / "meta.parquet")
    patients, _ = data.collapse_to_patient(meta, np.zeros((len(meta), 1)))
    ids = pq.read_table(INTERIM / "expression_hugo.parquet",
                        columns=["__index_level_0__"]).column(0).to_pylist()
    with_expr = set(pd.Index(ids).astype(str).str[:12])
    cohort = (patients[patients["patient_id"].isin(with_expr)]
              .rename(columns={"cancer type abbreviation": "cancer_type"})
              .reset_index(drop=True))
    surv = outcome.from_frame(cohort, endpoint="PFI")
    ctype = cohort["cancer_type"]
    valid = np.asarray(surv.valid)
    frame = pd.DataFrame({"type": ctype, "time": surv.time, "event": surv.event})[valid]
    by_type = frame.groupby("type")

    scores = {
        "event_proportion": ctype.map(by_type["event"].mean()),
        "events_per_person_time": ctype.map(by_type["event"].sum() / by_type["time"].sum()),
    }
    from statsmodels.duration.hazard_regression import PHReg

    X = pd.get_dummies(ctype[valid], drop_first=True, dtype=float).to_numpy()
    fit = PHReg(surv.time[valid], X, status=surv.event[valid]).fit(disp=False)
    lp = np.full(len(cohort), np.nan)
    lp[valid] = X @ fit.params
    scores["cox_type_dummies"] = pd.Series(lp)

    result = {
        "cohort": "pancancer", "endpoint": "PFI", "stratified": False, "in_sample": True,
        "n_patients": int(len(cohort)), "n_valid": int(valid.sum()),
        "n_events": int(surv.n_events), "n_types": int(ctype.nunique()),
        "c_index": {k: outcome.concordance_index(v.to_numpy(dtype=float), surv)
                    for k, v in scores.items()},
        "reported_as": "event_proportion",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    for k, v in result["c_index"].items():
        print(f"  {k:24s} C = {v:.6f}")
    print(f"{result['n_patients']} patients, {result['n_valid']} with PFI, "
          f"{result['n_events']} events, {result['n_types']} types -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
