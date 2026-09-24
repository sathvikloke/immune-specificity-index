#!/usr/bin/env python3
"""Triage the CPTAC lung slide inventory before committing to extraction.

    python scripts/07_triage_cptac.py

WHY THIS EXISTS
===============
Two earlier planning decisions were made about CPTAC without anyone opening the
manifest. The first dropped the cohort outright, on the grounds that its slides
are "qualification-workflow specimen-segment" slides rather than diagnostic ones.
The second reinstated it and re-budgeted the extraction, treating the
4.9-slides-per-subject figure purely as a cost lever. Both were reasoning from
the cohort's description; neither checked what is actually in it.

(Those decisions live in the project's planning documents, which sit ABOVE this
directory and are deliberately not part of the published pipeline. They are
paraphrased here rather than cited by path so this rationale still reads when the
pipeline is distributed on its own.)

This script does. It reads the TCIA CPTAC lung cohort table and answers three
questions the budget cannot:

  1. How many slides are TUMOUR? (a third are not)
  2. How many passed CPTAC's own segment QC?
  3. What tumour-nuclei content were the accepted segments selected for?

MEASURED, 2026-08-18 (2,138 slides, 441 cases)
----------------------------------------------
    normal_tissue                 765  (35.8%)  -- not tumour at all
    tumour, segment NOT acceptable 235
    tumour, segment acceptable   1,138  from 439 cases
        LUAD  615 slides / 229 cases
        LSCC  523 slides / 210 cases
    accepted segments: median 70% tumour nuclei (IQR 60-75), necrosis <= 20%

THE OBJECTION THAT ACTUALLY MATTERS, AND IT IS NOT SLIDE TYPE
=============================================================
CPTAC segments are SELECTED for high tumour-nuclei content. This project's
target is TME signatures — immune and stromal content — so the selection acts
directly on the compartment being predicted. Selecting for 60-75% tumour nuclei
systematically depletes stroma and immune infiltrate, which:

  - compresses the range of the signature scores, and
  - attenuates any correlation with the image by construction.

That is selection on the outcome variable, not merely a domain shift. It is
still usable as a validation cohort, because the direction of the bias is
predictable (attenuation, so a surviving ISI is a conservative result), but it
must be declared in advance and the range restriction must be reported next to
the result rather than discovered afterwards.

WHAT THIS SCRIPT CANNOT SETTLE
==============================
Whether the accepted slides are FFPE or frozen sections. The CPTAC GBM protocol
describes an ADJACENT segment being formalin-fixed, paraffin-embedded and H&E
stained for quality assessment, which would make them comparable to TCGA's DX
slides; the "top and bottom" sectioning language points at frozen OCT blocks,
which would not be. The public manifest does not record it. Resolve by
inspecting a handful of downloaded slides before committing the full run.

For reference the DISCOVERY side is 100% DX (diagnostic FFPE): all 10,170 slides
in the embeddings parquet. This docstring said 99.48% (10,117) until 2026-09-17;
the 53 it missed carry lettered suffixes (`-DXA` to `-DXU`), which a pattern
requiring a digit after `DX` does not match.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "raw" / "cptac" / "tcia-luad-lusc-cohort.csv"
SOURCE = ("https://raw.githubusercontent.com/GeorgeBatch/"
          "TCIA-CPTAC-lung-histology-download/master/tcia-luad-lusc-cohort.csv")

# Measured from the TCIA collection pages, 2026-08-12.
GB_TOTAL = 845.5
SLIDES_TOTAL_TCIA = 2218


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=MANIFEST)
    args = ap.parse_args()

    if not args.manifest.exists():
        print(f"Manifest not found at {args.manifest}\n  fetch it with:\n"
              f"  curl -sL {SOURCE} -o {args.manifest}")
        return 1

    d = pd.read_csv(args.manifest)
    print("=" * 74)
    print("CPTAC LUNG SLIDE TRIAGE")
    print("=" * 74)
    print(f"  manifest: {len(d)} rows, {d.Case_ID.nunique()} cases, "
          f"{d.Slide_ID.nunique()} slides")

    normal = d[d.Specimen_Type == "normal_tissue"]
    tumour = d[d.Specimen_Type == "tumor_tissue"]
    ok = tumour[tumour.Tumor_Segment_Acceptable == "Yes"]
    rejected = tumour[tumour.Tumor_Segment_Acceptable == "No"]

    print("\n  ATTRITION")
    print(f"    all slides                      {len(d):5d}")
    print(f"    - normal tissue                 {len(normal):5d}  ({len(normal)/len(d):.1%})")
    print(f"    - tumour, segment NOT acceptable {len(rejected):4d}")
    print(f"    = USABLE                        {len(ok):5d}  from {ok.Case_ID.nunique()} cases")
    print(f"      {len(ok)/max(ok.Case_ID.nunique(),1):.1f} slides per case")

    print("\n  BY COHORT (usable only)")
    by = ok.groupby("Tumor").agg(slides=("Slide_ID", "nunique"), cases=("Case_ID", "nunique"))
    for name, row in by.iterrows():
        print(f"    {name:6s} {row['slides']:4d} slides / {row['cases']:3d} cases")

    nuc = pd.to_numeric(ok.Tumor_Percent_Tumor_Nuclei, errors="coerce").dropna()
    nec = pd.to_numeric(ok.Tumor_Percent_Necrosis, errors="coerce").dropna()
    print("\n  SELECTION ON THE TARGET COMPARTMENT  <-- the real objection")
    print(f"    tumour nuclei %  median {nuc.median():.0f}  "
          f"IQR {nuc.quantile(.25):.0f}-{nuc.quantile(.75):.0f}  "
          f"min {nuc.min():.0f}  max {nuc.max():.0f}")
    print(f"    necrosis %       median {nec.median():.0f}  max {nec.max():.0f}")
    print("    These segments were CHOSEN for high tumour content. The target of")
    print("    this project is stromal/immune signature score, so the selection")
    print("    acts on the outcome and will attenuate the correlation. A surviving")
    print("    ISI is therefore conservative; a null one is uninterpretable.")

    frac = len(ok) / SLIDES_TOTAL_TCIA
    print("\n  REVISED BUDGET (usable slides only)")
    print(f"    {len(ok)} of {SLIDES_TOTAL_TCIA} TCIA slides = {frac:.0%} of the "
          f"{GB_TOTAL:.0f} GB collection ~ {GB_TOTAL*frac:.0f} GB")
    print(f"    at 14.6 slides/h that is {len(ok)/14.6:.0f} GPU-hours "
          f"(${len(ok)/14.6*0.34:.0f}-${len(ok)/14.6*0.69:.0f} on a 4090)")
    print("    Filtering first roughly halves both the download and the compute.")

    print("\n  UNRESOLVED: FFPE or frozen? Not in the manifest. Inspect a few")
    print("  downloaded slides before committing. Discovery side is 100% DX.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
