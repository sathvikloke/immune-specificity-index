#!/usr/bin/env python3
"""Does the A10 fix cost ssGSEA runtime?  Item 4 of the science list.

WHAT IS WRONG WITH THE EXISTING NUMBER
--------------------------------------
`results/README.md` records ssGSEA scoring at 0.0347 ms per patient-set for the
pinned run.  The pre-fix counterpart is 0.0759.  Both are single 6.5-second
micro-benchmarks taken inside much longer runs on contended machines, so the
2.19x ratio confounds the algorithm with whatever else the machine was doing --
and it points the WRONG WAY for the change actually made.  `kind="stable"`
selects a mergesort; the default is an introsort.  A stable sort is normally
SLOWER on ties, not 2x faster, so a ratio that says the fix made scoring faster
is evidence about the machines, not about the sort.

WHAT THIS MEASURES INSTEAD
--------------------------
The one line the A10 fix changed, on the real array it runs on:

    depth[i, np.argsort(-ranks[i])]                  <- shipped before
    depth[i, np.argsort(-ranks[i], kind="stable")]   <- shipped after

Ties are the whole point of A10 -- log expression has a floor, every gene
resting on it shares a rank, and 35,394,448 of 38,747,424 rank entries are tied
-- and tie density is exactly what separates these two algorithms.  So the
benchmark runs on the REAL NSCLC rank table, not a synthetic array.

HOW IT AVOIDS THE DEFECT IT IS REPLACING
----------------------------------------
  * both orderings in ONE process, on the SAME array, so no cross-run drift;
  * A/B INTERLEAVED, so a machine that slows down mid-benchmark penalises both
    equally instead of whichever ran second;
  * many repetitions, reporting the MEDIAN and the full spread rather than one
    sample -- a single timing is what produced the number being replaced;
  * a stated load context, because this project has been bitten repeatedly by
    timings quoted without one.
"""
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "results"
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 7
N_ROWS = int(sys.argv[2]) if len(sys.argv) > 2 else 200   # rows per repetition


def load_ranks() -> np.ndarray:
    """The real rank table the A10 line consumes."""
    expr = pd.read_parquet(REPO / "data/interim/expr_nsclc.parquet")
    print(f"expr: {expr.shape[0]} samples x {expr.shape[1]} genes", flush=True)
    # signatures.py ranks genes WITHIN each sample, descending.
    ranks = expr.rank(axis=1).to_numpy(dtype=np.float64)
    return ranks


def tie_profile(ranks: np.ndarray) -> dict:
    tied_rows = 0
    tied_elems = 0
    for i in range(ranks.shape[0]):
        u = np.unique(ranks[i])
        t = ranks.shape[1] - u.size
        if t:
            tied_rows += 1
            tied_elems += t
    return {"tied_rows": tied_rows, "rows": int(ranks.shape[0]),
            "tied_elements": int(tied_elems),
            "elements": int(ranks.size)}


def time_sort(ranks: np.ndarray, rows: np.ndarray, kind: str | None) -> float:
    depth = np.empty(ranks.shape[1], dtype=np.int32)
    countdown = np.arange(ranks.shape[1], 0, -1, dtype=np.int32)
    t0 = time.perf_counter()
    for i in rows:
        if kind is None:
            depth[np.argsort(-ranks[i])] = countdown
        else:
            depth[np.argsort(-ranks[i], kind=kind)] = countdown
    return time.perf_counter() - t0


def main() -> int:
    print(f"python {platform.python_version()}  numpy {np.__version__}  "
          f"{platform.system()}/{platform.machine()}", flush=True)
    try:
        batt = subprocess.run(["pmset", "-g", "batt"], capture_output=True,
                              text=True, timeout=10).stdout.strip().splitlines()[-1]
    except Exception:
        batt = "(unavailable)"
    print(f"power: {batt}", flush=True)

    ranks = load_ranks()
    prof = tie_profile(ranks)
    print(f"ties: {prof['tied_rows']}/{prof['rows']} rows carry ties; "
          f"{prof['tied_elements']:,} of {prof['elements']:,} entries tied", flush=True)

    rng = np.random.default_rng(0)
    rows = rng.choice(ranks.shape[0], size=min(N_ROWS, ranks.shape[0]), replace=False)
    print(f"timing {len(rows)} rows per repetition, {REPS} repetitions, "
          "A/B interleaved\n", flush=True)

    # one untimed pass so page faults and caches are not charged to whichever
    # ordering happens to run first
    time_sort(ranks, rows[:20], None)
    time_sort(ranks, rows[:20], "stable")

    quick, stable = [], []
    for r in range(REPS):
        q = time_sort(ranks, rows, None)
        s = time_sort(ranks, rows, "stable")
        quick.append(q)
        stable.append(s)
        print(f"  rep {r + 1}/{REPS}: quicksort {q:.4f}s   stable {s:.4f}s   "
              f"ratio {s / q:.4f}", flush=True)

    q = np.array(quick)
    s = np.array(stable)
    per_row_q = np.median(q) / len(rows) * 1000
    per_row_s = np.median(s) / len(rows) * 1000
    out = {
        "platform": f"{platform.system()}/{platform.machine()}",
        "numpy": np.__version__,
        "power": batt,
        "reps": REPS,
        "rows_per_rep": int(len(rows)),
        "ties": prof,
        "quicksort_s": list(map(float, q)),
        "stable_s": list(map(float, s)),
        "quicksort_median_s": float(np.median(q)),
        "stable_median_s": float(np.median(s)),
        "ratio_stable_over_quicksort_median": float(np.median(s) / np.median(q)),
        "ratio_min": float(np.min(s) / np.max(q)),
        "ratio_max": float(np.max(s) / np.min(q)),
        "ms_per_row_quicksort": float(per_row_q),
        "ms_per_row_stable": float(per_row_s),
    }
    print("\n" + "=" * 68)
    print(f"quicksort  median {np.median(q):.4f}s  ({per_row_q:.4f} ms/row)")
    print(f"stable     median {np.median(s):.4f}s  ({per_row_s:.4f} ms/row)")
    print(f"ratio stable/quicksort: median {out['ratio_stable_over_quicksort_median']:.4f}"
          f"   worst-case bracket [{out['ratio_min']:.4f}, {out['ratio_max']:.4f}]")
    print("=" * 68)
    Path(ROOT / "ssgsea_sort_bench.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {ROOT / 'ssgsea_sort_bench.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
