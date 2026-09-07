#!/usr/bin/env python3
"""Fetch every input the audit needs. Run this first.

Everything here is open-access and ungated. Nothing requires dbGaP, an EGA data
access committee, an institutional email, or a signing official.

    python scripts/01_fetch_data.py --all
    python scripts/01_fetch_data.py --embeddings --ancestry

Some artefacts cannot be fetched programmatically (supplementary tables behind
publisher UI, HuggingFace datasets needing the `datasets` client). For those
this script prints exact instructions and the destination path rather than
pretending to download them.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"


EMBEDDINGS_INSTRUCTIONS = f"""
SLIDE-LEVEL PROV-GIGAPATH EMBEDDINGS (primary substrate)
--------------------------------------------------------
Repo:    seandavis/tcga_provgigapath_embeddings   (HuggingFace, dataset)
Licence: CC-BY-4.0, ungated
Size:    ~466 MB, single parquet, ~11,948 TCGA slides
Encoder: Prov-GigaPath (Apache-2.0) — no non-commercial restriction

  pip install huggingface_hub
  huggingface-cli download seandavis/tcga_provgigapath_embeddings \\
      --repo-type dataset --local-dir {RAW / 'provgigapath'}

VERIFY ON THE DATASET CARD before writing methods text (red-team could not
confirm these): tiling parameters, magnification, and whether vectors are
mean-pooled tile embeddings or a CLS token.

DELIBERATELY NOT USED:
  MahmoodLab/UNI2-h-features   gated; gate text says gmail will be denied;
                               UNI2-h weights are CC-BY-NC
  W8Yi/tcga-wsi-uni2h-features ungated but derived from UNI2-h, so it inherits
                               the non-commercial restriction
Both are avoided because of the industry affiliation. If that stops mattering,
`data.load_embeddings` accepts any parquet with a barcode and a vector column.
"""

ANCESTRY_INSTRUCTIONS = f"""
TCGA GENETIC ANCESTRY CALLS
---------------------------
Carrot-Zhang J, Chambwe N, Damrauer JS, et al. "Comprehensive Analysis of
Genetic Ancestry and Its Molecular Correlates in Cancer."
Cancer Cell 2020;37(5):639-654.e6.  doi:10.1016/j.ccell.2020.04.012

10,678 patients, 33 cancer types. This is what makes the ancestry arm powered
pan-cancer (AFR ~717, EAS ~535, AMR ~249) where it is hopeless in any single
tumour type (NSCLC has EAS n=17).

Get Table S1 from either:
  1. GDC PanCanAtlas publication page (preferred, stable):
     https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020
  2. The article's supplementary material:
     https://www.cell.com/cancer-cell/fulltext/S1535-6108(20)30211-7
  3. Open-access mirror:
     https://escholarship.org/uc/item/8q63d5vs

Save as: {RAW / 'ancestry' / 'carrot_zhang_2020_TableS1.xlsx'}

On first load, `data.load_ancestry` guesses the patient and ancestry column
names and records its guesses in the provenance file. CHECK THEM, then pass the
names explicitly.
"""

PURITY_INSTRUCTIONS = f"""
TUMOUR PURITY (ABSOLUTE)
------------------------
Needed for the purity decomposition, which is one of the two scientific cores.
Aran D, Sirota M, Butte AJ. "Systematic pan-cancer analysis of tumour purity."
Nat Commun 2015;6:8971  — provides ESTIMATE, ABSOLUTE, LUMP, IHC and CPE.

Source options:
  - PanCanAtlas: https://gdc.cancer.gov/about-data/publications/pancanatlas
    (look for the ABSOLUTE purity/ploidy file)
  - The Aran et al. supplementary table (CPE consensus column)

Save as: {RAW / 'purity' / 'tcga_purity.csv'}
Required columns: a TCGA patient barcode, and at least one purity estimate.
Prefer CPE (consensus); fall back to ABSOLUTE. Record which you used — the
red-team flagged that H&E-predicted ABSOLUTE purity reaches Spearman 0.418-0.655
(Oner et al., Patterns 2021), which is the number the decomposition argues with.
"""

SIGNATURES_INSTRUCTIONS = f"""
TME GENE SIGNATURES
-------------------
NOT bundled with this package, deliberately: signature definitions must come
from their sources so provenance is auditable.

Options:
  1. HistoTME's curated set (the direct precedent this audit targets):
     https://github.com/spatkar94/HistoTME
     Its ground-truth file carries the signature columns; the underlying
     definitions trace to Bagaev et al., Cancer Cell 2021 (doi:10.1016/j.ccell.2021.04.014).
  2. MSigDB Hallmark / C7 immunologic sets (GMT format):
     https://www.gsea-msigdb.org/gsea/msigdb
  3. Any GMT you trust.

Save as: {RAW / 'signatures' / 'tme_signatures.gmt'}  (or .json)

Loaders: `signatures.SignatureSet.from_gmt()` / `.from_json()`.
A single-signature smoke test is available via `signatures.example_signature_set()`
(CYT = GZMA + PRF1, Rooney et al. Cell 2015) so you can exercise the pipeline
before the full set is in place.
"""

EXPRESSION_INSTRUCTIONS = f"""
TCGA BULK RNA-SEQ
-----------------
Needed to compute signature ground truth AND to run the Venet null (random gene
sets must be scored through the identical pipeline, which requires the full
expression matrix, not just precomputed signature scores).

Open tier, no dbGaP. Options:
  - GDC Data Portal / API, "STAR - Counts" for TCGA projects:
    https://portal.gdc.cancer.gov/
    A pan-cancer pull is ~45 GB at ~4.2 MB per file.
  - Recount3 or the UCSC Xena pan-cancer TOIL matrix are far smaller and
    already harmonised; prefer one of these if bandwidth is a constraint.
    https://xenabrowser.net/datapages/

Save as: {RAW / 'expression' / 'tcga_expression.parquet'}
Shape: samples x genes, log-transformed, HUGO symbols as columns.

WARNING that bites people: if you z-score or median-scale WITHIN cohort, then
TCGA-derived and CPTAC-derived ground truth sit on different scales and your
"external validation" target has been silently redefined by cohort composition.
Decide the scaling policy once, write it down, and apply it identically.
"""


def ensure_dirs() -> None:
    for sub in ("provgigapath", "ancestry", "purity", "signatures", "expression"):
        (RAW / sub).mkdir(parents=True, exist_ok=True)


def fetch_embeddings() -> bool:
    """Attempt the HF download; fall back to printing instructions."""
    dest = RAW / "provgigapath"
    try:
        from huggingface_hub import snapshot_download
    except ImportError:
        print(EMBEDDINGS_INSTRUCTIONS)
        print("!! huggingface_hub not installed:  pip install huggingface_hub")
        return False

    print(f"Downloading seandavis/tcga_provgigapath_embeddings -> {dest}")
    try:
        snapshot_download(
            repo_id="seandavis/tcga_provgigapath_embeddings",
            repo_type="dataset",
            local_dir=str(dest),
        )
    except Exception as exc:  # noqa: BLE001 - report anything and continue
        print(f"!! download failed: {exc}")
        print(EMBEDDINGS_INSTRUCTIONS)
        return False

    files = list(dest.rglob("*.parquet"))
    print(f"   done: {len(files)} parquet file(s)")
    for f in files:
        print(f"     {f.relative_to(RAW)}  {f.stat().st_size / 1e6:.1f} MB")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="fetch/describe everything")
    parser.add_argument("--embeddings", action="store_true")
    parser.add_argument("--ancestry", action="store_true")
    parser.add_argument("--purity", action="store_true")
    parser.add_argument("--signatures", action="store_true")
    parser.add_argument("--expression", action="store_true")
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        return 1

    ensure_dirs()
    want = lambda flag: args.all or flag  # noqa: E731

    if want(args.embeddings):
        fetch_embeddings()
    if want(args.ancestry):
        print(ANCESTRY_INSTRUCTIONS)
    if want(args.purity):
        print(PURITY_INSTRUCTIONS)
    if want(args.signatures):
        print(SIGNATURES_INSTRUCTIONS)
    if want(args.expression):
        print(EXPRESSION_INSTRUCTIONS)

    print(f"\nRaw data root: {RAW}")
    print("Next:  python scripts/02_run_audit.py --demo    (synthetic smoke test)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
