#!/usr/bin/env python3
"""Build the six `data/interim/` files the cohort scripts read, from `data/raw/`.

    python scripts/00_build_interim.py --out data/interim            # everything
    python scripts/00_build_interim.py --out DIR --cohort nsclc      # the NSCLC three only
    python scripts/00_build_interim.py --out DIR --tumour-only       # A16 sensitivity inputs
    python scripts/00_build_interim.py --out DIR --verify            # compare with recorded hashes

WHY THIS EXISTS (14-SCIENCE-AUDIT.md, A13)
==========================================
Until 2026-09-16 nothing in this repository wrote `data/interim/`: the six files
were made on 2026-08-17 by code that was never committed, so a reader with the
raw sources could not start the full reproduction path. This script is that
missing step, reconstructed from the files themselves and verified against
their recorded sha256 prefixes (RECORDED below; results under VERIFIED). Where a choice
below looks odd, it is because the frozen files were made that way, and the
frozen results were computed from them.

WHAT EACH FILE IS
=================
meta.parquet (10,170 slides x 40)
    `data.load_embeddings(embeddings.parquet, barcode_col="sample",
    embedding_col="embedding")`: the embeddings table's own clinical columns
    (TCGA-CDR, bundled with the Prov-GigaPath release), barcode fallback from
    `filename` for 124 null `sample` values, primary-tumour rows only. The
    parquet repeats a slide once per clinical sample of its patient (11,948
    rows, 10,542 slide files); the primary filter leaves one row per slide.
X.npy (10,170 x 768)
    The LAST of the 14 slide-encoder layers, in `meta` row order.
expression_hugo.parquet (10,535 samples x 41,046 genes)
    UCSC Xena TOIL `tcga_RSEM_gene_tpm`, transposed to samples x genes, Ensembl
    ids mapped to HUGO symbols by `ensembl_to_hugo.csv` (version suffix
    stripped), first occurrence of a symbol kept, unmapped ids dropped. Rows keep
    TOIL's column order and its 15-character sample ids, ALL sample types.
    Values are TOIL's own scale, log2(TPM + 0.001) -- its floor is -9.9658 =
    log2(0.001). They are NOT log2(TPM + 1), whatever older docstrings say.
X_nsclc.npy, cohort_nsclc.parquet, expr_nsclc.parquet (944 patients)
    LUAD + LUSC slides from `meta`, collapsed to patients (mean embedding),
    `cancer type abbreviation` renamed `cancer_type`, ancestry (UCSF calls) and
    ABSOLUTE purity left-joined, restricted to patients with ANY TOIL sample,
    then 16 `SIG_` mean-z scores (the TME Hallmark GMT) computed on
    `expr_nsclc` as read back from disk -- numpy's summation order depends on
    memory layout, and the frozen scores were computed on the re-read frame.
    `expr_nsclc` takes each patient's FIRST TOIL sample in TOIL column order.

A16 -- THE EXPRESSION IS NOT ALWAYS THE TUMOUR (found 2026-09-16)
=================================================================
Neither cohort restricts expression to primary tumour (sample code 01):
  * NSCLC keeps the first TOIL sample per patient. For 51 of 944 patients that
    is not the tumour: 48 adjacent-normal (11) where a tumour sample exists,
    2 recurrences (02), and 1 patient with only a normal sample.
  * Pan-cancer (09_run_pancancer.py and others) averages every TOIL sample of
    a patient after truncating the id: 578 of 7,168 patients are averages
    (520 tumour + normal) and 40 have no tumour sample at all (28 normal-only).
This script reproduces the frozen files as they are. `--tumour-only` builds
the corrected inputs (code 01 only; patients without one are dropped) for the
sensitivity analysis, into a DIFFERENT directory.

VERIFIED (2026-09-16)
=====================
HPC4, Linux/x86_64, pyarrow 24.0.0 (job 129998): meta.parquet, X.npy,
expression_hugo.parquet, X_nsclc.npy and expr_nsclc.parquet byte-identical;
cohort_nsclc.parquet identical except its 16 SIG_ scores, which differ by at
most 1.1e-15 (Linux rounds the mean-z reductions differently). macOS arm64,
pyarrow 24.0.0: all five files it can build within the laptop's memory
(everything except expression_hugo.parquet) byte-identical, cohort_nsclc
included.

MEMORY
======
The full TOIL matrix is 60,498 x 10,535 float64 (4.7 GiB) before the transpose
copy. Build pan-cancer where that fits (HPC4). `--cohort nsclc` reads only the
NSCLC patients' TOIL columns (about 1,036) and runs on a laptop.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import data, signatures  # noqa: E402

NSCLC_TYPES = ("LUAD", "LUSC")
GMT_NAME = "h.all.v2024.1.Hs.symbols.tme.gmt"

# sha256 prefixes of the files the frozen results were computed from, recorded
# in 14-SCIENCE-AUDIT.md (A13) and identical on macOS and HPC4. Byte identity
# additionally depends on the writer (pyarrow 24.0.0 wrote the parquet files).
RECORDED = {
    "meta.parquet": "86305b9f7c97950f",
    "X.npy": "6048e461aab97f87",
    "expression_hugo.parquet": "2a36f4a3d7f36a43",
    "cohort_nsclc.parquet": "0b396a1e2780b0b0",
    "X_nsclc.npy": "7a1408eac2f64f7c",
    "expr_nsclc.parquet": "8532e31cf0176acd",
}

# The one file whose BYTES depend on the platform. Linux/x86_64 with pyarrow
# 24.0.0 wrote exactly these bytes on two separate HPC4 runs (jobs 129998 and
# 130188, 2026-09-16/17), and job 129998 showed their content equals the frozen
# file except the 16 SIG_ columns, by at most 1.1e-15. A match here is reported
# as that known build, never as the frozen file itself.
RECORDED_PLATFORM = {
    "cohort_nsclc.parquet": {
        "7fe02f65d5083280": "the recorded Linux/x86_64 build (HPC4 jobs 129998, 130188); "
                            "SIG_ columns within 1.1e-15 of the frozen file",
    },
}


def sha16(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()[:16]


def _write_parquet(frame: pd.DataFrame, path: Path) -> None:
    # pandas >= 2.1 stores DataFrame.attrs in the file. The frozen files carry
    # none (assemble_cohort's join_coverage would otherwise change the bytes).
    frame = frame.copy()
    frame.attrs = {}
    frame.to_parquet(path)


def toil_header(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def read_toil(path: Path, samples: list[str] | None = None) -> pd.DataFrame:
    """TOIL genes x samples. `samples` restricts the columns read (header order kept).

    pandas' DEFAULT float parser, deliberately: it is what made the frozen files,
    and A15 showed it is not always exact. On TOIL's four-decimal values the
    rebuilt expression files are byte-identical on both platforms, so it is kept.
    """
    usecols = None
    if samples is not None:
        want = set(samples)
        usecols = [h for i, h in enumerate(toil_header(path)) if i == 0 or h in want]
    return pd.read_csv(path, sep="\t", index_col=0, usecols=usecols)


def to_hugo(toil: pd.DataFrame, gene_map: Path) -> pd.DataFrame:
    """Samples x HUGO symbols, first occurrence of a symbol kept, unmapped ids dropped."""
    mp = pd.read_csv(gene_map)
    lookup = dict(zip(mp.iloc[:, 0].astype(str), mp.iloc[:, 1].astype(str)))
    frame = toil.T
    frame.columns = [lookup.get(str(g).split(".")[0], str(g)) for g in toil.index]
    frame = frame.loc[:, ~frame.columns.duplicated()]
    frame = frame.loc[:, [not str(c).startswith("ENSG") for c in frame.columns]]
    frame.index = pd.Index(frame.index.to_numpy())
    return frame


def _sample_mix(ids: pd.Index, patients: set[str]) -> dict:
    codes = pd.Series(ids.str[13:15], index=ids.str[:12])
    codes = codes[codes.index.isin(patients)]
    per = codes.groupby(level=0).agg(lambda c: tuple(sorted(c)))
    return {"patients": int(len(per)),
            "more_than_one_sample": int((per.map(len) > 1).sum()),
            "no_primary_tumour_sample": int((~per.map(lambda t: "01" in t)).sum())}


def build(raw: Path, out: Path, *, cohort: str, tumour_only: bool) -> dict[str, Path]:
    t0 = time.time()
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    toil_path = raw / "expression" / "tcga_RSEM_gene_tpm.gz"
    gene_map = raw / "expression" / "ensembl_to_hugo.csv"
    header = toil_header(toil_path)[1:]
    if tumour_only:
        header = [h for h in header if h[13:15] == "01"]
    toil_patients = {h[:12] for h in header}

    meta, X, prov = data.load_embeddings(
        raw / "provgigapath" / "embeddings.parquet",
        barcode_col="sample", embedding_col="embedding")
    print(f"[{time.time() - t0:6.0f}s] meta {meta.shape}, X {X.shape}; "
          + "; ".join(prov.notes[1:4]), flush=True)
    if cohort in ("all", "pancancer"):
        _write_parquet(meta, out / "meta.parquet")
        np.save(out / "X.npy", X)
        written.update({"meta.parquet": out / "meta.parquet", "X.npy": out / "X.npy"})

    if cohort in ("all", "pancancer"):
        expr = to_hugo(read_toil(toil_path, header if tumour_only else None), gene_map)
        print(f"[{time.time() - t0:6.0f}s] expression {expr.shape}", flush=True)
        pan_patients = set(data.collapse_to_patient(meta, np.zeros((len(meta), 1)))[0]["patient_id"])
        mix = _sample_mix(expr.index, pan_patients)
        print(f"  A16 pan-cancer: {mix}", flush=True)
        _write_parquet(expr, out / "expression_hugo.parquet")
        written["expression_hugo.parquet"] = out / "expression_hugo.parquet"
        del expr

    if cohort in ("all", "nsclc"):
        keep = meta["cancer type abbreviation"].isin(NSCLC_TYPES).to_numpy()
        pts, Xn = data.collapse_to_patient(meta[keep].reset_index(drop=True), X[keep])
        pts = pts.rename(columns={"cancer type abbreviation": "cancer_type"})
        ancestry, _ = data.load_ancestry(raw / "ancestry" / "UCSF_Ancestry_Calls.csv")
        purity, _ = data.load_purity(raw / "purity" / "TCGA_ABSOLUTE_purity.csv")
        co = data.assemble_cohort(pts, ancestry=ancestry, purity=purity)
        print(f"  NSCLC slide patients {len(co)}; join coverage "
              f"{co.attrs.get('join_coverage')}", flush=True)
        has = co["patient_id"].isin(toil_patients).to_numpy()
        co, Xn = co.loc[has].reset_index(drop=True), Xn[has]
        samples = [h for h in header if h[:12] in set(co["patient_id"])]
        e = to_hugo(read_toil(toil_path, samples), gene_map)
        # TOIL column order, then the FIRST sample per patient (see A16).
        order = {h: i for i, h in enumerate(header)}
        e = e.loc[sorted(e.index, key=order.get)]
        print(f"  A16 NSCLC: {_sample_mix(e.index, set(co['patient_id']))}", flush=True)
        e.index = e.index.str[:12]
        e = e[~e.index.duplicated(keep="first")].reindex(pd.Index(co["patient_id"].to_numpy()))
        _write_parquet(e, out / "expr_nsclc.parquet")
        np.save(out / "X_nsclc.npy", Xn)
        e = pd.read_parquet(out / "expr_nsclc.parquet")      # score the re-read frame
        gmt = raw / "signatures" / GMT_NAME
        sigs = signatures.SignatureSet.from_gmt(gmt)
        sigs = signatures.SignatureSet(
            name=sigs.name, sets={f"SIG_{k}": v for k, v in sigs.sets.items()}, source=str(gmt))
        scores = signatures.score_mean_z(e, sigs.filter_to(set(e.columns)))
        for c in scores.columns:
            co[c] = scores[c].to_numpy()
        _write_parquet(co, out / "cohort_nsclc.parquet")
        print(f"[{time.time() - t0:6.0f}s] NSCLC {len(co)} patients, {len(scores.columns)} "
              "signatures", flush=True)
        written.update({n: out / n for n in
                        ("cohort_nsclc.parquet", "X_nsclc.npy", "expr_nsclc.parquet")})
    return written


# Largest absolute difference in a float column that `verify` accepts as
# platform arithmetic. Measured 2026-09-16 on HPC4 (job 129998): the rebuilt
# cohort_nsclc.parquet matched the macOS-built file in every column except the
# 16 SIG_ mean-z scores, which differed by at most 1.1e-15 -- the reductions in
# `score_mean_z` round differently on Linux/x86_64. Everything else, and the
# other five files, were byte-identical.
FLOAT_TOLERANCE = 1e-12


def _frame_difference(a: pd.DataFrame, b: pd.DataFrame) -> tuple[bool, float, list[str]]:
    """(non-float content identical, max |diff| over float columns, differing columns)."""
    if list(a.columns) != list(b.columns) or not a.index.equals(b.index):
        return False, float("inf"), ["<shape, columns or index>"]
    worst, cols, exact = 0.0, [], True
    for c in a.columns:
        x, y = a[c], b[c]
        if x.dtype.kind == "f" and y.dtype.kind == "f":
            xa, ya = x.to_numpy(), y.to_numpy()
            if not np.array_equal(xa, ya, equal_nan=True):
                if not np.array_equal(np.isnan(xa), np.isnan(ya)):
                    return False, float("inf"), [c]
                worst = max(worst, float(np.nanmax(np.abs(xa - ya))))
                cols.append(c)
        elif not x.equals(y):
            exact = False
            cols.append(c)
    return exact, worst, cols


def verify(files: dict[str, Path], reference: Path | None) -> int:
    """Hash against RECORDED; where bytes differ, compare CONTENT with `reference`.

    Content that differs only in float columns, by at most FLOAT_TOLERANCE, is
    reported with its size and accepted; anything else counts as a failure.
    Without a reference, bytes matching a RECORDED_PLATFORM build are accepted
    and named as that build.
    """
    bad = 0
    for name, path in sorted(files.items()):
        got = sha16(path)
        if got == RECORDED[name]:
            print(f"  BYTE-IDENTICAL {name} {got}")
            continue
        known = RECORDED_PLATFORM.get(name, {}).get(got)
        if known is not None and reference is None:
            print(f"  KNOWN BUILD    {name} {got}: {known}")
            continue
        verdict = "no reference to compare content against"
        if reference is not None and (reference / name).exists():
            ref = reference / name
            if name.endswith(".npy"):
                A, B = np.load(path), np.load(ref)
                exact, worst = A.shape == B.shape, (float(np.max(np.abs(A - B)))
                                                    if A.shape == B.shape else float("inf"))
                cols = [] if worst == 0 else ["<array>"]
            else:
                exact, worst, cols = _frame_difference(pd.read_parquet(path), pd.read_parquet(ref))
            if exact and worst == 0:
                verdict = "CONTENT-EQUAL"
            elif exact and worst <= FLOAT_TOLERANCE:
                verdict = (f"NUMERICALLY EQUAL: {len(cols)} float column(s) differ by at "
                           f"most {worst:.2e}, all else identical")
            else:
                verdict = f"CONTENT DIFFERS in {cols[:5]} (max float |diff| {worst:.3g})"
                bad += 1
        else:
            bad += 1
        print(f"  BYTES DIFFER   {name} {got} (recorded {RECORDED[name]}): {verdict}")
    return bad


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", type=Path, default=ROOT / "data" / "raw")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cohort", choices=["all", "nsclc", "pancancer"], default="all")
    ap.add_argument("--tumour-only", action="store_true",
                    help="A16 sensitivity: primary-tumour (01) expression only; "
                         "patients without one are dropped. Never write this into "
                         "data/interim/.")
    ap.add_argument("--verify", action="store_true",
                    help="hash the outputs against RECORDED; compare content with "
                         "--reference where bytes differ")
    ap.add_argument("--reference", type=Path, default=None)
    ap.add_argument("--overwrite", action="store_true",
                    help="allow writing over existing files in --out")
    args = ap.parse_args(argv)

    if args.tumour_only and args.out.resolve() == (ROOT / "data" / "interim").resolve():
        raise SystemExit("FATAL: --tumour-only must not write into data/interim/")
    names = {"all": list(RECORDED),
             "pancancer": ["meta.parquet", "X.npy", "expression_hugo.parquet"],
             "nsclc": ["cohort_nsclc.parquet", "X_nsclc.npy", "expr_nsclc.parquet"]}[args.cohort]
    clash = [n for n in names if (args.out / n).exists()]
    if clash and not args.overwrite:
        raise SystemExit(f"FATAL: {args.out} already holds {clash}; refusing to overwrite")

    files = build(args.raw, args.out, cohort=args.cohort, tumour_only=args.tumour_only)
    if args.verify and not args.tumour_only:
        bad = verify(files, args.reference)
        print("OK: every output reproduces the frozen interim file." if not bad
              else f"NOT VERIFIED: {bad} output(s) differ from the recorded bytes and "
                   "could not be shown equal in content.")
        # 3, not 1: the files WERE built; the caller decides whether an
        # unverified file is fatal (a crash exits 1 through the traceback).
        return 3 if bad else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
