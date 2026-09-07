#!/usr/bin/env python3
"""Fetch TME gene signatures under a licence an industry-affiliated author can use.

THE PROBLEM THIS SOLVES
=======================
The ~29 TME signatures used by HistoTME trace to Bagaev et al., Cancer Cell 2021
(BostonGene). HistoTME's README points at BostonGene/MFP to compute them. That
repository's licence is academic/non-profit ONLY, verbatim:

    "use the Software solely for academic and non-profit purposes"

    "Any use, reproduction, or distribution of the Software or Derivative
     Software for direct or indirect commercial (including strategic) gain,
     purpose, or advantage, INCLUDING FOR ANY RESEARCH AND/OR DEVELOPMENT
     PURPOSE BY A FOR-PROFIT ENTITY OR ON BEHALF OF A FOR-PROFIT ENTITY,
     requires a separately executed written license agreement."

and its definitions section explicitly sweeps data in, not just code:

    '"Source Form" shall mean software source code, documentation source,
     configuration files, AND DATA.'

So a for-profit-affiliated researcher cannot use the MFP gene lists without a
separate agreement (askusepermission@bostongene.com).

WHAT IS ACTUALLY AVAILABLE
==========================
1. MSigDB — Creative Commons Attribution 4.0. Commercial use permitted with
   attribution. Verified downloadable with no login (HTTP 200). THIS IS THE
   DEFAULT AND IT UNBLOCKS THE PROJECT TODAY.
   Caveat: gene sets with a "KEGG_" or "KEGG_MEDICUS_" prefix carry additional
   Kanehisa Laboratories restrictions and are EXCLUDED by this script.

2. HistoTME itself is BSD 3-Clause (permissive, no commercial restriction) and
   ships `example_data/pantcga_tme_signatures.csv` — 8,024 TCGA samples x 29
   signature SCORES. Usable as a cross-check target, but it contains scores, not
   gene lists, so it cannot support the random-set null on its own.

3. Bagaev et al., Cancer Cell 2021 (doi:10.1016/j.ccell.2021.04.014) publishes
   the signature definitions in its supplement, under journal terms rather than
   BostonGene's software licence. Whether that is a usable alternative route to
   the same lists is a LEGAL question, not a technical one. Do not assume it;
   ask your employer's counsel.

    python scripts/03_fetch_signatures.py --hallmark
    python scripts/03_fetch_signatures.py --hallmark --c7
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data" / "raw" / "signatures"

MSIGDB_BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"
MSIGDB_VERSION = "2024.1.Hs"

COLLECTIONS = {
    "hallmark": ("h.all", "Hallmark — 50 well-curated sets, the default"),
    "c7": ("c7.all", "C7 immunologic signatures — large, immune-specific"),
    "c2cp": ("c2.cp", "C2 canonical pathways (KEGG subsets are filtered out)"),
}

# Sets carrying restrictions beyond CC-BY. Excluded on principle.
RESTRICTED_PREFIXES = ("KEGG_", "KEGG_MEDICUS_", "BIOCARTA_")

# The subset of Hallmark that maps onto TME / immune biology. Sizes are from
# v2024.1.Hs and are worth noting: 36-200 genes, i.e. MUCH larger than the
# BostonGene panels. Larger k pushes you into the DEGENERACY regime rather than
# the attenuation regime — `globalaxis.null_degeneracy` reports which one you are
# in, and both corrections are implemented. Pre-register the choice.
TME_RELEVANT = (
    "HALLMARK_ALLOGRAFT_REJECTION",            # 200
    "HALLMARK_ANGIOGENESIS",                   #  36
    "HALLMARK_COAGULATION",                    # 138
    "HALLMARK_COMPLEMENT",                     # 200
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",  # 200
    "HALLMARK_HYPOXIA",                        # 200
    "HALLMARK_INFLAMMATORY_RESPONSE",          # 200
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",      #  97
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",      # 200
    "HALLMARK_TGF_BETA_SIGNALING",             #  54
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    # Proliferation — include deliberately. These are the global-axis proxies,
    # and having them in the panel makes the axis argument concrete rather than
    # theoretical (cf. Venet's meta-PCNA control).
    "HALLMARK_G2M_CHECKPOINT",                 # 200
    "HALLMARK_E2F_TARGETS",                    # 200
    "HALLMARK_MYC_TARGETS_V1",                 # 200
)

ATTRIBUTION = (
    "MSigDB is CC-BY 4.0, (c) 2004-2025 Broad Institute, MIT, and Regents of the "
    "University of California. Cite: Liberzon et al., Cell Syst 2015;1:417-425 "
    "(Hallmark) and Subramanian et al., PNAS 2005;102:15545-15550 (GSEA)."
)


def fetch(collection: str) -> Path | None:
    stem, _ = COLLECTIONS[collection]
    fname = f"{stem}.v{MSIGDB_VERSION}.symbols.gmt"
    url = f"{MSIGDB_BASE}/{MSIGDB_VERSION}/{fname}"
    DEST.mkdir(parents=True, exist_ok=True)
    out = DEST / fname

    print(f"  {url}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read()
    except Exception as exc:  # noqa: BLE001
        print(f"  !! failed: {exc}")
        return None

    if len(data) < 2000:
        print(f"  !! suspiciously small ({len(data)} bytes) — not written")
        return None

    out.write_bytes(data)
    print(f"  -> {out.name}  ({len(data) / 1000:.1f} KB)")
    return out


def filter_gmt(path: Path, *, tme_only: bool) -> Path:
    """Drop restricted sets, and optionally keep only the TME-relevant subset."""
    kept, dropped_restricted, dropped_other = [], 0, 0
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        name = line.split("\t", 1)[0]
        if name.startswith(RESTRICTED_PREFIXES):
            dropped_restricted += 1
            continue
        if tme_only and name not in TME_RELEVANT:
            dropped_other += 1
            continue
        kept.append(line)

    suffix = ".tme.gmt" if tme_only else ".permissive.gmt"
    out = path.with_suffix(suffix)
    out.write_text("\n".join(kept) + "\n")

    print(f"  -> {out.name}: {len(kept)} sets kept")
    if dropped_restricted:
        print(f"     {dropped_restricted} dropped (KEGG/BioCarta extra restrictions)")
    if dropped_other:
        print(f"     {dropped_other} dropped (not in the TME-relevant subset)")

    sizes = sorted(len(line.split("\t")) - 2 for line in kept)
    if sizes:
        median = sizes[len(sizes) // 2]
        print(f"     set sizes: min={sizes[0]} median={median} max={sizes[-1]}")
        print(f"     -> at k~{median}, expect a {'DEGENERATE' if median > 50 else 'ATTENUATED'} null; "
              "confirm with globalaxis.null_degeneracy()")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hallmark", action="store_true", help="Hallmark (recommended default)")
    ap.add_argument("--c7", action="store_true", help="C7 immunologic signatures")
    ap.add_argument("--c2cp", action="store_true", help="C2 canonical pathways")
    ap.add_argument("--all-sets", action="store_true",
                    help="keep every permissive set, not just the TME subset")
    args = ap.parse_args()

    wanted = [c for c in ("hallmark", "c7", "c2cp") if getattr(args, c.replace("2cp", "2cp"))]
    if not wanted:
        wanted = ["hallmark"]
        print("No collection specified; defaulting to --hallmark\n")

    print(f"MSigDB {MSIGDB_VERSION} — CC-BY 4.0, commercial use permitted with attribution.\n")
    ok = []
    for c in wanted:
        print(f"[{c}] {COLLECTIONS[c][1]}")
        p = fetch(c)
        if p:
            ok.append(filter_gmt(p, tme_only=(c == "hallmark" and not args.all_sets)))
        print()

    if not ok:
        print("Nothing fetched. Check connectivity, or download manually from")
        print("https://www.gsea-msigdb.org/gsea/msigdb/")
        return 1

    (DEST / "ATTRIBUTION.txt").write_text(ATTRIBUTION + "\n")
    print(ATTRIBUTION)
    print(f"\nWritten to {DEST}")
    print("\nLoad with:")
    print("  from aacr27.signatures import SignatureSet")
    print(f"  sigs = SignatureSet.from_gmt('{ok[0]}')")
    print("\nNOTE: these are NOT the BostonGene/HistoTME 29 signatures. They are a")
    print("permissively-licensed substitute. State that explicitly in the methods —")
    print("it is a strength (fully reproducible without a licence gate), not a")
    print("weakness, but it must not be misrepresented as replicating HistoTME.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
