#!/usr/bin/env python3
"""Power of the outcome arm to see a given C-index advantage (ledger E10).

    python scripts/30_outcome_power.py            # NSCLC; writes results/outcome_power_nsclc.json

WHY
===
NSCLC's outcome arm found 0 of 32 signature-by-adjustment comparisons beating
their random-set null, and the manuscript attaches a minimum detectable effect
(80% power) of about 0.046 in C-index. A single MDE answers "what could we have
seen at 80%"; it does not answer "could we have seen 0.03?", which is the size
of the pan-cancer advantage. This computes the power curve from the frozen
table, with the SAME normal approximation `outcome.outcome_null` uses for the
MDE, so the two numbers cannot disagree:

    se      = (excess_hi - excess_lo) / (2 * z_0.975)          per comparison
    delta_z = arctanh(2 * delta_C)      (the inverse of the MDE's tanh(z) / 2)
    power   = Phi(delta_z / se - z_0.975) + Phi(-delta_z / se - z_0.975)

Reads only `results/<cohort>_v3/outcome_arm.csv`, so it runs from the deposit.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from pathlib import Path

from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DELTAS = (0.02, 0.03, 0.04, 0.05, 0.06)
Z = stats.norm.ppf(0.975)


def power(delta_c: float, se: float) -> float:
    dz = math.atanh(2 * delta_c)
    return float(stats.norm.cdf(dz / se - Z) + stats.norm.cdf(-dz / se - Z))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cohort", default="nsclc", choices=["nsclc", "pancancer"])
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    src = ROOT / "results" / f"{args.cohort}_v3" / "outcome_arm.csv"
    with open(src, newline="") as fh:
        rows = list(csv.DictReader(fh))
    ses = [(float(r["excess_hi"]) - float(r["excess_lo"])) / (2 * Z) for r in rows]
    # The MDE by outcome_null's own formula. Frozen tables written before that
    # column existed (NSCLC's) lack it; where it IS stored, it must agree.
    mde = [math.tanh((Z + stats.norm.ppf(0.80)) * s) / 2 for s in ses]
    for r, m in zip(rows, mde):
        if r.get("mde_80pct_c_index") not in (None, ""):
            assert abs(float(r["mde_80pct_c_index"]) - m) < 1e-12, (r["signature"], m)
    # At delta_C = MDE the power is 0.80 plus the negligible opposite tail.
    for s, m in zip(ses, mde):
        assert abs(power(m, s) - 0.80) < 1e-3, (s, m, power(m, s))
    curve = {}
    for d in DELTAS:
        p = [power(d, s) for s in ses]
        curve[f"{d:.2f}"] = {"median": statistics.median(p), "min": min(p), "max": max(p)}
    out = args.out or ROOT / "results" / f"outcome_power_{args.cohort}.json"
    result = {"cohort": args.cohort, "source": str(src.relative_to(ROOT)),
              "n_comparisons": len(rows),
              "n_beating_null": sum(r["beats_null"] == "True" for r in rows),
              "mde_80pct_c_index_median": statistics.median(mde),
              "power_by_c_index_advantage": curve}
    out.write_text(json.dumps(result, indent=2) + "\n")
    for d, v in curve.items():
        print(f"  C-index advantage {d}: power median {v['median']:.3f} "
              f"(range {v['min']:.3f}-{v['max']:.3f})")
    print(f"{len(rows)} comparisons -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
