#!/usr/bin/env python3
"""Is the pan-cancer outcome arm's five-of-five composition beyond chance? (ledger E11)

    python scripts/31_outcome_category_test.py     # writes results/outcome_category_test.json

WHY
===
Pan-cancer, five signatures beat their random-set null on progression-free
interval (G2M checkpoint, E2F targets, angiogenesis, epithelial-mesenchymal
transition, hypoxia; `results/pancancer_v3/outcome_arm.csv`, axis-residualised).
All five sit in the six signatures that MSigDB's process categories file under
proliferation, development or pathway (`10_make_figures.py`'s CAT, the table the
manuscript and the checker already use). The pre-specified contrast is immune
against everything else, reported in Methods. This asks the set-level question
the Results prose implies: how often would five winners drawn from the 16 land
entirely inside those six?

Two exact answers, both by enumeration (no sampling, so no seed):
  count   P(all k winners inside the group) under a uniformly random choice of
          k of 16 signatures: C(g, k) / C(16, k).
  excess  the difference in mean axis-residualised outcome excess between the
          group and the rest, against all C(16, g) relabellings (one-sided,
          ">=" the observed difference).

READ WITH TWO CAVEATS, both recorded in the output:
  * the three-category grouping is named AFTER seeing which signatures won, so
    neither number tests a pre-specified hypothesis; they size the pattern;
  * relabelling assumes exchangeable signatures, and these are not: G2M and E2F
    share genes and move together, so the effective number of independent
    winners is below five and the true tail is heavier than enumerated.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "pancancer_v3" / "outcome_arm.csv"
OUT = ROOT / "results" / "outcome_category_test.json"
GROUP = ("proliferation", "development", "pathway")


def _categories() -> dict[str, str]:
    spec = importlib.util.spec_from_file_location("figs10", ROOT / "scripts" / "10_make_figures.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return dict(mod.CAT)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=OUT)
    out = ap.parse_args(argv).out
    d = pd.read_csv(SRC, float_precision="round_trip")
    d = d[d["adjustment"] == "axis_residualised"].copy()
    d["short"] = d["signature"].str.replace("SIG_HALLMARK_", "", regex=False)
    cat = _categories()
    missing = sorted(set(d["short"]) - set(cat))
    if missing:
        raise SystemExit(f"no category for {missing}")
    d["cat"] = d["short"].map(cat)
    n = len(d)
    in_group = d["cat"].isin(GROUP).to_numpy()
    g = int(in_group.sum())
    winners = d.loc[d["beats_null"].astype(bool), "short"].tolist()
    k = len(winners)
    inside = int(d.loc[d["beats_null"].astype(bool), "cat"].isin(GROUP).sum())

    p_count = math.comb(g, inside) * math.comb(n - g, k - inside) / math.comb(n, k) \
        if inside == k else None
    # One-sided tail for "at least this many winners inside the group".
    p_count_tail = sum(math.comb(g, j) * math.comb(n - g, k - j)
                       for j in range(inside, min(g, k) + 1)) / math.comb(n, k)

    x = d["excess_z"].to_numpy(dtype=float)
    total = float(x.sum())
    observed = float(x[in_group].mean() - x[~in_group].mean())
    at_least, count = 0, 0
    for idx in itertools.combinations(range(n), g):
        s = float(x[list(idx)].sum())
        diff = s / g - (total - s) / (n - g)
        count += 1
        at_least += diff >= observed - 1e-12
    result = {
        "source": str(SRC.relative_to(ROOT)),
        "adjustment": "axis_residualised",
        "group_categories": list(GROUP),
        "n_signatures": n,
        "group_size": g,
        "group_members": sorted(d.loc[in_group, "short"]),
        "winners": sorted(winners),
        "winners_inside_group": inside,
        "p_all_winners_inside_group": p_count,
        "p_at_least_this_many_inside": p_count_tail,
        "mean_excess_difference_group_minus_rest": observed,
        "relabellings_enumerated": count,
        "p_excess_difference_one_sided": at_least / count,
        "caveats": [
            "the grouping was named after the winners were known: descriptive, not a test "
            "of a pre-specified hypothesis (the pre-specified contrast is immune vs other)",
            "relabelling assumes exchangeable signatures; G2M and E2F share genes, so the "
            "true tail is heavier than enumerated",
        ],
    }
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{inside} of {k} winners inside the {g}-signature group "
          f"(P all inside = {p_count}); excess difference {observed:.4f}, "
          f"one-sided P = {at_least}/{count} = {at_least / count:.5f}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
