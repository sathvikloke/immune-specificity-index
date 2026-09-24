#!/usr/bin/env python3
"""Fetch every raw input the audit needs, and prove each one is the file the results used.

    python scripts/01_fetch_data.py --all          # fetch what is missing, verify everything
    python scripts/01_fetch_data.py --verify       # verify only; download nothing

Everything here is open access. Nothing requires dbGaP, a data access
committee, an institutional email or a signing official.

WHAT THIS WRITES (rewritten 2026-09-17, ledger C4)
==================================================
Until 2026-09-17 this script downloaded only the embeddings -- under the
dataset's own file name, which nothing downstream reads -- and printed
instructions naming files (`carrot_zhang_2020_TableS1.xlsx`, `tcga_purity.csv`,
`tme_signatures.gmt`, `tcga_expression.parquet`) that no script opens, and its
purity note preferred CPE where `data.load_purity` prefers ABSOLUTE. A reader
following it could not reach `00_build_interim.py`. Every entry in RAW_FILES is
now the exact path the pipeline reads, the public source it came from, and the
sha256 of the copy the frozen results were computed from. Each source below
was re-downloaded on 2026-09-16/17 and reproduced that hash, except where noted.

    data/raw/provgigapath/embeddings.parquet
        HuggingFace dataset seandavis/tcga_provgigapath_embeddings, revision
        073115403c2fc5134ee8d1332c603edba591dddb, file
        provgigapath_embeddings_with_metadata.parquet (CC-BY-4.0, ungated).
    data/raw/expression/tcga_RSEM_gene_tpm.gz
        UCSC Xena TOIL hub, log2(TPM + 0.001).
    data/raw/expression/ensembl_to_hugo.csv
        HGNC's complete set as served on 2026-08-17, reduced to the rows with an
        Ensembl id: `ensembl_gene_id,symbol`, in HGNC's order. HGNC (CC0) does
        not archive that snapshot -- its 2026-08-04 and 2026-08-07 monthly
        archives map 41,037 of TOIL's genes where this map maps 41,046 -- so the
        map itself is committed at pipeline/resources/ensembl_to_hugo.csv and
        copied from there. `derive_gene_map()` rebuilds it from any HGNC file and
        reports whether the result matches.
    data/raw/ancestry/UCSF_Ancestry_Calls.csv
        GDC publication page CCG-AIM-2020 (Carrot-Zhang et al. 2020).
    data/raw/purity/TCGA_ABSOLUTE_purity.tsv
        GDC PanCanAtlas ABSOLUTE, `TCGA_mastercalls.abs_tables_JSedit.fixed.txt`.
    data/raw/purity/TCGA_ABSOLUTE_purity.csv
        The same file with tabs replaced by commas (what `data.load_purity` reads).
    data/raw/signatures/h.all.v2024.1.Hs.symbols.gmt, and its .tme.gmt subset
        MSigDB 2024.1.Hs Hallmark, via `03_fetch_signatures.py`, which also
        writes the 16-set TME subset the analysis uses.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RESOURCES = ROOT / "resources"

HF_REPO = "seandavis/tcga_provgigapath_embeddings"
HF_REVISION = "073115403c2fc5134ee8d1332c603edba591dddb"
HF_FILE = "provgigapath_embeddings_with_metadata.parquet"
GDC = "https://api.gdc.cancer.gov/data/"

# path under data/raw -> (sha256, how to get it)
RAW_FILES: dict[str, tuple[str, str]] = {
    "provgigapath/embeddings.parquet": (
        "234379b0d85cce26817358d7c63f707042ee8fd8f7145c18940a6893fb5deed6", "huggingface"),
    "expression/tcga_RSEM_gene_tpm.gz": (
        "2ac7215f35fbe2cdc671c03c2b40b934ac1a0c908757f05f6b8267e5ba0a1b6d",
        "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz"),
    "expression/ensembl_to_hugo.csv": (
        "0315430e6dddde71c79fcc1f38a0b053507828f2bf6e47e53c6172cc9a7ad54b", "resources"),
    "ancestry/UCSF_Ancestry_Calls.csv": (
        "8e7b1f92c3d957ac9ab04e40fe731120ffa05328352126afe2ae89641f65fdb9",
        GDC + "fdfa536a-c3c8-405d-99d9-bc9375b5084c"),
    "purity/TCGA_ABSOLUTE_purity.tsv": (
        "f430a975433d82e0098d7405619d4f12a0c765fcd97e7d63cc9b1de7f2d763cd",
        GDC + "4f277128-f793-4354-a13d-30cc7fe9f6b5"),
    "purity/TCGA_ABSOLUTE_purity.csv": (
        "7a96dac3a253b7d0888cff5dc068256699ff2ffce8126f6a75a63fb4537becd5", "tabs-to-commas"),
    "signatures/h.all.v2024.1.Hs.symbols.gmt": (
        "ee2463540042078bfa3f67828e1e223bb354446d9fbb4d22845866835ba5c772",
        "03_fetch_signatures.py"),
    "signatures/h.all.v2024.1.Hs.symbols.tme.gmt": (
        "098e63374eea867925911f5136d0f147ff1e73257a18de91931a5be7f6e55c94",
        "03_fetch_signatures.py"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _download(url: str, dest: Path) -> None:
    tmp = dest.with_name(dest.name + ".part")
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, length=1 << 20)
    tmp.replace(dest)


def _fetch_embeddings(dest: Path) -> None:
    from huggingface_hub import hf_hub_download

    got = hf_hub_download(repo_id=HF_REPO, repo_type="dataset", revision=HF_REVISION,
                          filename=HF_FILE)
    shutil.copyfile(got, dest)


def derive_gene_map(hgnc_complete_set: Path, dest: Path) -> str:
    """Rebuild the Ensembl -> HUGO map from an HGNC complete set; return its sha256."""
    import pandas as pd

    h = pd.read_csv(hgnc_complete_set, sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id"])
    h.dropna(subset=["ensembl_gene_id"])[["ensembl_gene_id", "symbol"]].to_csv(dest, index=False)
    return sha256(dest)


def obtain(rel: str, how: str) -> None:
    dest = RAW / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    if how == "huggingface":
        _fetch_embeddings(dest)
    elif how == "resources":
        shutil.copyfile(RESOURCES / dest.name, dest)
    elif how == "tabs-to-commas":
        src = dest.with_suffix(".tsv")
        if not src.exists():
            raise FileNotFoundError(f"{src} is needed first")
        dest.write_bytes(src.read_bytes().replace(b"\t", b","))
    elif how.startswith("https://"):
        _download(how, dest)
    else:
        raise FileNotFoundError(f"run `python scripts/{how}` to create {rel}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--all", action="store_true", help="fetch what is missing, then verify")
    ap.add_argument("--verify", action="store_true", help="verify only; download nothing")
    args = ap.parse_args(argv)
    if not (args.all or args.verify):
        ap.print_help()
        return 1

    bad = 0
    for rel, (want, how) in RAW_FILES.items():   # insertion order: .tsv before .csv
        path = RAW / rel
        if not path.exists() and args.all:
            try:
                print(f"  fetching {rel} <- {how}", flush=True)
                obtain(rel, how)
            except Exception as exc:  # noqa: BLE001 - report every missing input, then fail
                print(f"  COULD NOT FETCH {rel}: {type(exc).__name__}: {exc}")
        if not path.exists():
            print(f"  MISSING        {rel}  (source: {how})")
            bad += 1
            continue
        got = sha256(path)
        if got == want:
            print(f"  VERIFIED       {rel}")
        else:
            print(f"  HASH DIFFERS   {rel}: {got[:16]} (expected {want[:16]}) -- "
                  "not the file the frozen results were computed from")
            bad += 1
    print(f"\nRaw data root: {RAW}")
    if bad:
        print(f"FAILED: {bad} input(s) missing or different.")
        return 1
    print("OK: every raw input is the file the frozen results were computed from.")
    print("Next:  python scripts/00_build_interim.py --out data/interim --verify")
    return 0


if __name__ == "__main__":
    sys.exit(main())
