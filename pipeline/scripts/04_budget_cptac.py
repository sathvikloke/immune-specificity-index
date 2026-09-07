#!/usr/bin/env python3
"""Cost and time budget for extracting Prov-GigaPath features from CPTAC lung.

WHY THIS EXISTS
===============
Week-1 verification established that the public embeddings parquet contains
**zero CPTAC patients**. So external validation is not a matter of downloading
someone else's features — it means running Prov-GigaPath ourselves. This script
prices that, from measured anchors rather than assertions, so the scope decision
is made deliberately in August rather than discovered in October.

VERIFIED INPUTS (fetched from TCIA, 2026-08-12)
----------------------------------------------
CPTAC-LUAD  243 subjects  1,137 SVS  431.5 GB histopathology
CPTAC-LSCC  212 subjects  1,081 SVS  414.0 GB histopathology
TOTAL       455 subjects  2,218 SVS  845.5 GB

TWO DENOMINATORS, RECONCILED 2026-09-04 -- READ THIS BEFORE QUOTING A COST
--------------------------------------------------------------------------
The numbers above are the TCIA **collection** totals. They are correct and they
are NOT what will be extracted. The TCIA **cohort manifest**
(`data/raw/cptac/tcia-luad-lusc-cohort.csv`) is a different, smaller set, and
17-CPTAC-PLAN.md works from it:

    TCIA collection pages            455 subjects   2,218 slides   <- this file
    cohort manifest CSV              441 cases      2,138 slides
    after the pre-registered filter  439 cases      1,138 slides
      (tumor_tissue AND Tumor_Segment_Acceptable == "Yes")
    after dropping the 4 `11LU`      435 cases      1,134 slides   <- the real job

Both sets were labelled "verified" in different documents and the difference was
never stated, so cost figures derived from each were being quoted side by side.

**The working set is 1,134 slides, which is 51% of 2,218.** Every per-slide
figure this script prints for "every slide" is therefore about 2x the real
extraction cost: option C's $52-105 corresponds to roughly **$27-54** on the
actual worklist, which is where that second pair of numbers came from. The
conclusion is unchanged and was never close -- the GPU is not the expensive part
either way -- but quote one denominator and say which.

This script is deliberately NOT rewritten to the manifest denominator: it prices
the collection as fetched from TCIA, which is a real and separately useful
figure. `scripts/20_cptac_manifest.py` (planned) is what emits the frozen
worklist, and any cost quoted in the paper or a budget request should be derived
from that artefact, not from here.

Note the ratio: 2,218 / 455 = 4.9 slides per subject. That is the fingerprint of
the CPTAC "qualification workflow" — multiple top/bottom slides per case, not one
diagnostic slide. It is also the single biggest cost lever, because you do not
need every slide to validate on every patient.

THROUGHPUT ANCHORS (from the landscape audit, measured not modelled)
--------------------------------------------------------------------
  SurGen:  1,020 WSIs, UNI (ViT-L/16, ~300M), 1x V100 32GB = 58 GPU-h = 17.5 slides/h
  THUNDER: ~4,000 WSIs, Virchow2 (632M),      1x V100      =  7.8 slides/h

Prov-GigaPath's tile encoder is ViT-g/14 (~1.1B params), so roughly 3x UNI's
per-tile cost; an RTX 4090 in fp16 is roughly 2.5x a V100 for this workload.
Net: about 14-15 slides/h on a 4090, which is the central estimate below. The
range spans the two anchors scaled the same way.

    python scripts/04_budget_cptac.py
    python scripts/04_budget_cptac.py --slides-per-patient 1 --gpu 4090
"""

from __future__ import annotations

import argparse
import sys

# --- Verified TCIA figures, 2026-08-12 ---
COLLECTIONS = {
    "CPTAC-LUAD": {"subjects": 243, "slides": 1137, "gb": 431.5},
    "CPTAC-LSCC": {"subjects": 212, "slides": 1081, "gb": 414.0},
}

# --- Rental prices (USD/hour), landscape audit section 7.6 ---
GPUS = {
    "4090":  (0.34, 0.69, "RTX 4090, community/interruptible tiers"),
    "l40s":  (0.79, 0.79, "L40S"),
    "a100":  (1.07, 2.06, "A100 80GB"),
    "h100":  (1.50, 3.00, "H100, neocloud"),
    "h100-hyperscaler": (6.88, 12.29, "H100 at AWS/GCP/Azure — DO NOT USE"),
}

# Slides/hour for Prov-GigaPath tile encoding, scaled from the measured anchors.
THROUGHPUT = {"low": 9.0, "central": 14.6, "high": 21.0}

STORAGE_USD_PER_GB_MONTH = 0.10
DOWNLOAD_MB_PER_S = 50.0          # conservative for TCIA


def budget(slides: int, gb: float, subjects: int, gpu: str,
           streaming: bool, batch_slides: int) -> dict:
    lo_rate, hi_rate, _ = GPUS[gpu]

    hours = {k: slides / v for k, v in THROUGHPUT.items()}
    cost = {k: (h * lo_rate, h * hi_rate) for k, h in hours.items()}

    download_h = (gb * 1024) / DOWNLOAD_MB_PER_S / 3600

    if streaming:
        peak_gb = gb * (batch_slides / slides) if slides else gb
        # Disk is held only for the run, so charge a fraction of a month.
        storage_usd = peak_gb * STORAGE_USD_PER_GB_MONTH * (hours["central"] / 720)
    else:
        peak_gb = gb
        storage_usd = gb * STORAGE_USD_PER_GB_MONTH

    return {
        "subjects": subjects, "slides": slides, "gb": gb,
        "gpu_hours": hours, "gpu_cost": cost,
        "download_hours": download_h,
        "peak_disk_gb": peak_gb, "storage_usd": storage_usd,
        "wall_clock_days_1gpu": hours["central"] / 24,
    }


def show(name: str, b: dict) -> None:
    lo, hi = b["gpu_cost"]["central"]
    print(f"\n  {name}")
    print(f"    {b['subjects']} subjects · {b['slides']} slides · {b['gb']:.0f} GB")
    print(f"    GPU hours   {b['gpu_hours']['high']:.0f} – {b['gpu_hours']['low']:.0f} "
          f"(central {b['gpu_hours']['central']:.0f})")
    print(f"    GPU cost    ${lo:.0f} – ${hi:.0f}")
    print(f"    Download    {b['download_hours']:.1f} h at {DOWNLOAD_MB_PER_S:.0f} MB/s")
    print(f"    Peak disk   {b['peak_disk_gb']:.0f} GB  (storage ${b['storage_usd']:.2f})")
    print(f"    Wall clock  {b['wall_clock_days_1gpu']:.1f} days on ONE GPU "
          f"({b['wall_clock_days_1gpu'] * 24 / 4:.0f} h on four)")
    print(f"    TOTAL       ${lo + b['storage_usd']:.0f} – ${hi + b['storage_usd']:.0f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gpu", default="4090", choices=list(GPUS))
    ap.add_argument("--slides-per-patient", type=int, default=0,
                    help="0 = all slides; 1 or 2 = cap per patient (the big lever)")
    ap.add_argument("--no-streaming", action="store_true",
                    help="hold the whole collection on disk instead of batching")
    ap.add_argument("--batch-slides", type=int, default=50)
    args = ap.parse_args()

    subjects = sum(c["subjects"] for c in COLLECTIONS.values())
    all_slides = sum(c["slides"] for c in COLLECTIONS.values())
    all_gb = sum(c["gb"] for c in COLLECTIONS.values())
    per_slide_gb = all_gb / all_slides

    print("=" * 74)
    print("CPTAC LUNG — Prov-GigaPath FEATURE EXTRACTION BUDGET")
    print("=" * 74)
    print(f"  GPU: {GPUS[args.gpu][2]}  (${GPUS[args.gpu][0]:.2f}–${GPUS[args.gpu][1]:.2f}/h)")
    print(f"  Verified: {subjects} subjects, {all_slides} slides, {all_gb:.0f} GB "
          f"({all_slides / subjects:.1f} slides/subject)")

    print("\n" + "-" * 74)
    print("SCOPE OPTIONS")
    print("-" * 74)

    options = [
        ("A. One slide per patient", subjects),
        ("B. Two slides per patient", min(subjects * 2, all_slides)),
        ("C. Every slide", all_slides),
    ]
    if args.slides_per_patient:
        options = [(f"Requested: {args.slides_per_patient}/patient",
                    min(subjects * args.slides_per_patient, all_slides))]

    for name, n in options:
        show(name, budget(n, n * per_slide_gb, subjects, args.gpu,
                          not args.no_streaming, args.batch_slides))

    print("\n" + "-" * 74)
    print("READ THIS BEFORE CHOOSING")
    print("-" * 74)
    print("""
  The GPU is NOT the expensive part. Every option above costs well under the
  $300-600 extension envelope. The decision is methodological, not financial.

  ONE SLIDE PER PATIENT (A) is cheapest but introduces a researcher degree of
  freedom: CPTAC slides are qualification slides (top/bottom of the block), so
  *which* one you keep is a choice that can be made after seeing results. If you
  take this option, pre-register the rule (e.g. lowest filename ordinal) and
  state it.

  EVERY SLIDE (C) removes that choice and lets collapse_to_patient() average
  them, which is what the pipeline already does for TCGA. It costs about 5x the
  GPU time and is still cheap. For consistency with the TCGA arm, C is the
  defensible default.

  STORAGE IS THE REAL TRAP, and only if you get it wrong. Holding all 845 GB on
  a persistent volume costs ~$85/month for nothing. Stream in batches of ~50
  slides, extract, delete, repeat — peak disk stays near 20 GB.

  DOWNLOAD IS THE WALL-CLOCK BOTTLENECK, not compute, unless you overlap them.
  Prefetch batch n+1 while encoding batch n or the GPU idles for hours.

  ONE FRICTION: prov-gigapath is gated:auto on HuggingFace — a request that is
  auto-approved. Not the MahmoodLab situation (which states gmail will be
  denied), but do it in week 1 rather than the night you rent the GPU.

  WHAT THIS BUYS: a genuine external cohort, no TCGA patient overlap, and full
  control of tiling/level/pooling on BOTH sides — which resolves the week-1
  magnification ambiguity by construction, because you set it yourself.
""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
