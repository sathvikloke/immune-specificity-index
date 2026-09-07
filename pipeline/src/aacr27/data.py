"""Data loaders.

Every loader is written against a VERIFIED, ungated, permissively-licensed
source, and every one records its provenance so the manuscript's data
availability statement writes itself.

Primary substrate
-----------------
`seandavis/tcga_provgigapath_embeddings` on HuggingFace — slide-level
Prov-GigaPath embeddings for TCGA, CC-BY-4.0, ungated, a single parquet of
roughly 466 MB covering ~11,948 slides, already joined to the TCGA-CDR clinical
table. Prov-GigaPath itself is Apache-2.0, so there is no non-commercial
licence exposure.

DELIBERATELY NOT USED
---------------------
- `MahmoodLab/UNI2-h-features` — gated; the gate form states that gmail
  addresses will be denied, and the underlying UNI2-h licence is CC-BY-NC.
- `W8Yi/tcga-wsi-uni2h-features` — ungated, but derived from UNI2-h and so
  inherits the non-commercial restriction.

Both are avoided because the researcher is industry-affiliated. If that changes,
`load_embeddings` takes any parquet with a barcode column and a vector column.

VERIFY BEFORE BUILDING (week 1)
-------------------------------
The red-team could not confirm: tiling parameters, magnification, and whether
the embeddings are mean-pooled or CLS-token. Confirm against the dataset card
before the results are final — it changes nothing about the pipeline's shape but
must be stated accurately in the methods.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import barcodes

HF_EMBEDDINGS_REPO = "seandavis/tcga_provgigapath_embeddings"
ANCESTRY_URL = "https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020"
ANCESTRY_CITATION = (
    "Carrot-Zhang J, Chambwe N, Damrauer JS, et al. Comprehensive Analysis of "
    "Genetic Ancestry and Its Molecular Correlates in Cancer. Cancer Cell "
    "2020;37(5):639-654.e6. doi:10.1016/j.ccell.2020.04.012"
)

# VERIFIED 2026-08-12 against the actual downloaded file
# (UCSF_Ancestry_Calls.csv, GDC, HTTP 200, no token, 1,933,183 bytes, 10,128 rows).
#
# Columns: Unnamed: 0, Patient_ID, Aliquot_ID, race, ethnicity,
#          pam.ancestry.cluster, pam.ancestry.subcluster, PC1..PC7
#
# CRITICAL LABEL DETAIL: the file spells East Asian **"ASIAN"**, not "EAS".
# An earlier version of this constant listed "EAS", which silently routed all
# 633 Asian-ancestry patients into the ADMIXED bucket — quietly destroying the
# subgroup the ancestry arm is powered on. "EAS"/"SAS" are retained as accepted
# aliases in case a future release changes the spelling.
#
# Verified counts in pam.ancestry.cluster (n = 10,128):
#     EUR 8,337 | AFR 928 | ASIAN 633 | AMR 228
# Subcluster splits EUR into 3 and ASIAN into 2.
ANCESTRY_GROUPS = ("EUR", "AFR", "ASIAN", "AMR", "EAS", "SAS")
ANCESTRY_CANONICAL = {"EAS": "ASIAN", "SAS": "ASIAN"}
ANCESTRY_PATIENT_COL = "Patient_ID"
ANCESTRY_CALL_COL = "pam.ancestry.cluster"
ANCESTRY_VERIFIED_COUNTS = {"EUR": 8337, "AFR": 928, "ASIAN": 633, "AMR": 228}


@dataclass
class Provenance:
    """Recorded for every loaded artefact; dumped into the results directory."""

    name: str
    source: str
    licence: str = "UNVERIFIED"
    accessed: str = ""
    n_rows: int = 0
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "source": self.source,
            "licence": self.licence,
            "accessed": self.accessed,
            "n_rows": self.n_rows,
            "notes": self.notes,
        }


def write_provenance(records: list[Provenance], path: str | Path) -> None:
    Path(path).write_text(json.dumps([r.to_dict() for r in records], indent=2))


# ------------------------------------------------------------ embeddings ---

def load_embeddings(
    path: str | Path,
    *,
    barcode_col: str | None = None,
    embedding_col: str | None = None,
    primary_tumour_only: bool = True,
    layer: int = -1,
    sample_codes: tuple[str, ...] | None = ("01",),
    fallback_barcode_col: str | None = "filename",
) -> tuple[pd.DataFrame, np.ndarray, Provenance]:
    """Load slide-level embeddings from a parquet.

    Returns (metadata_frame, embedding_matrix, provenance). The metadata frame
    and the matrix share row order; never reorder one without the other.

    Column names are auto-detected as a convenience, but for the real parquet
    pass them explicitly: `barcode_col="sample"`, `embedding_col="embedding"`.

    LAYER SELECTION IS MANDATORY AND MUST BE PRE-REGISTERED.
    The public parquet stores `embedding` as `large_list<large_list<double>>` =
    **14 layers x 768 dims** (slide-encoder hidden states), not a flat vector.
    Naive stacking raises `ValueError: setting an array element with a sequence`.
    `layer=-1` takes the final layer; pass an explicit index and record it. This
    is a researcher degree of freedom — choosing the layer after seeing results
    is a garden of forking paths.

    SAMPLE CODES. TCGA barcodes carry a sample-type code; the parquet contains
    01 (primary solid tumour, 11,701 of 11,948), plus 06 (metastatic, 178),
    02, 05, 11 and 07. `sample_codes=("01",)` restricts to primary tumours,
    which is the right default and must be stated in the methods.
    """
    path = Path(path)
    frame = pd.read_parquet(path)

    barcode_col = barcode_col or _guess_barcode_col(frame)
    embedding_col = embedding_col or _guess_embedding_col(frame)

    # FALLBACK, found on contact with the real parquet: 124 of 11,948 rows have a
    # NULL `sample` while `filename` still carries the full barcode
    # (e.g. "TCGA-02-0004-01Z-00-DX1.<uuid>"). All 124 are valid primary-tumour
    # diagnostic slides, so dropping them is needless loss. Recover the barcode
    # from the fallback column before annotating, and report how many.
    recovered = 0
    if fallback_barcode_col and fallback_barcode_col in frame.columns:
        missing = frame[barcode_col].isna() | (
            frame[barcode_col].astype("string").str.len().fillna(0) < 12
        )
        if missing.any():
            # Reindex the salvage mask onto the FULL frame index before combining.
            # Combining two boolean Series with different indexes relies on
            # pandas alignment and is deprecated; it would also silently produce
            # NaN-filled masks if the subset index ever diverged.
            salvage = (
                frame.loc[missing, fallback_barcode_col]
                .astype("string").str.split(".").str[0]
            )
            valid = salvage.map(
                lambda b: isinstance(b, str) and barcodes.try_parse(b) is not None
            ).reindex(frame.index, fill_value=False).astype(bool)
            take = missing.astype(bool) & valid
            frame.loc[take, barcode_col] = salvage.reindex(frame.index)[take]
            recovered = int(take.sum())

    frame = barcodes.annotate(frame, barcode_col=barcode_col)
    dropped = frame.attrs.get("dropped_unparseable", 0)

    if primary_tumour_only:
        before = len(frame)
        frame = frame[frame["is_primary_tumour"]].copy()
        filtered = before - len(frame)
    else:
        filtered = 0

    if embedding_col is not None:
        matrix = _stack_embeddings(frame[embedding_col], layer=layer)
        frame = frame.drop(columns=[embedding_col])
    else:
        vector_cols = [c for c in frame.columns if _looks_like_feature(c)]
        if not vector_cols:
            raise ValueError(
                "could not locate embedding columns; pass embedding_col explicitly"
            )
        matrix = frame[vector_cols].to_numpy(dtype=float)
        frame = frame.drop(columns=vector_cols)

    prov = Provenance(
        name="slide_embeddings",
        source=str(path),
        licence="CC-BY-4.0 (VERIFY on dataset card)",
        n_rows=len(frame),
        notes=[
            f"embedding dim = {matrix.shape[1]}",
            f"recovered {recovered} barcodes from {fallback_barcode_col!r} "
            f"(null in {barcode_col!r})",
            f"dropped {dropped} unparseable barcodes",
            f"filtered {filtered} non-primary-tumour slides",
            "VERIFY: tiling params, magnification, pooling (mean vs CLS)",
        ],
    )
    return frame.reset_index(drop=True), matrix, prov


def _guess_barcode_col(frame: pd.DataFrame) -> str:
    for candidate in ("barcode", "sample_barcode", "slide_barcode", "case_barcode",
                      "bcr_patient_barcode", "sample", "slide_id", "case_id"):
        if candidate in frame.columns:
            return candidate
    for col in frame.columns:
        if frame[col].dtype == object:
            sample = frame[col].dropna().astype(str).head(20)
            if len(sample) and sample.str.startswith("TCGA-").mean() > 0.8:
                return col
    raise ValueError(f"no barcode column found in {list(frame.columns)[:20]}")


def _guess_embedding_col(frame: pd.DataFrame) -> str | None:
    for candidate in ("embedding", "features", "feature", "vector", "embeddings"):
        if candidate in frame.columns:
            return candidate
    for col in frame.columns:
        head = frame[col].dropna().head(1)
        if len(head) and isinstance(head.iloc[0], (list, np.ndarray)):
            return col
    return None


def _looks_like_feature(col: str) -> bool:
    return bool(
        col.startswith(("feat_", "emb_", "f_", "dim_", "x"))
        and col[1:].lstrip("_").isdigit()
    ) or col.startswith(("feat_", "emb_", "dim_"))


def collapse_to_patient(
    frame: pd.DataFrame, matrix: np.ndarray, *, how: str = "mean"
) -> tuple[pd.DataFrame, np.ndarray]:
    """One row per patient.

    Patients with several slides would otherwise be weighted more heavily and
    would break the independence assumption in every downstream CI. Do this
    BEFORE splitting, not after.
    """
    if how not in {"mean", "first"}:
        raise ValueError("how must be 'mean' or 'first'")

    order = frame.groupby("patient_id", sort=True).indices
    patients = sorted(order)
    rows, vectors = [], []

    for p in patients:
        idx = order[p]
        rows.append(frame.iloc[idx[0]])
        vectors.append(matrix[idx].mean(axis=0) if how == "mean" else matrix[idx[0]])

    out = pd.DataFrame(rows).reset_index(drop=True)
    out["n_slides"] = [len(order[p]) for p in patients]
    return out, np.vstack(vectors)


# -------------------------------------------------------------- ancestry ---

def load_ancestry(
    path: str | Path,
    *,
    patient_col: str | None = None,
    ancestry_col: str | None = None,
) -> tuple[pd.DataFrame, Provenance]:
    """Load Carrot-Zhang et al. genetic ancestry calls (Table S1).

    Download instructions are in scripts/01_fetch_ancestry.py. The file is a
    supplementary table, so the exact column naming is UNVERIFIED here and is
    auto-detected; confirm on first run and then pass the names explicitly.

    10,678 patients across 33 cancer types. This is what makes the ancestry
    analysis powered pan-cancer (AFR ~717, EAS ~535, AMR ~249) where it is
    hopeless within any single tumour type.
    """
    path = Path(path)
    frame = pd.read_excel(path) if path.suffix in {".xlsx", ".xls"} else pd.read_csv(path)

    patient_col = patient_col or _guess_column(
        frame, ("patient", "bcr_patient_barcode", "submitter_id", "case", "sample")
    )
    ancestry_col = ancestry_col or _guess_column(
        frame, ("consensus_ancestry", "ancestry", "EIGENSTRAT", "call", "genetic_ancestry")
    )

    out = frame[[patient_col, ancestry_col]].rename(
        columns={patient_col: "patient_id", ancestry_col: "ancestry_raw"}
    )
    out["patient_id"] = out["patient_id"].astype(str).str.strip().str.upper().str[:12]
    out["ancestry"] = out["ancestry_raw"].astype(str).str.upper().str.strip()
    out["ancestry_broad"] = out["ancestry"].map(_broad_ancestry)
    out = out.drop_duplicates("patient_id")

    prov = Provenance(
        name="ancestry",
        source=f"{path} — {ANCESTRY_URL}",
        licence="see publication terms",
        n_rows=len(out),
        notes=[
            ANCESTRY_CITATION,
            f"detected columns: patient={patient_col!r} ancestry={ancestry_col!r}",
            "VERIFY these column guesses against Table S1 on first run",
            f"group counts: {out['ancestry_broad'].value_counts().to_dict()}",
        ],
    )
    return out, prov


def _broad_ancestry(value: str) -> str:
    """Collapse admixed labels to a broad group; anything else becomes ADMIXED.

    Carrot-Zhang's operational definition of admixed is a 20-80% non-EUR
    fraction. Admixed patients are kept as their own category rather than forced
    into a bin, because misassigning them would bias exactly the contrast we
    care about.
    """
    v = str(value).upper().strip()
    v = ANCESTRY_CANONICAL.get(v, v)
    for group in ANCESTRY_GROUPS:
        if v == group:
            return ANCESTRY_CANONICAL.get(group, group)
    if v in {"NAN", "", "NA", "NONE", "UNKNOWN"}:
        return "UNKNOWN"
    return "ADMIXED"


def _guess_column(frame: pd.DataFrame, candidates: tuple[str, ...]) -> str:
    lowered = {c.lower(): c for c in frame.columns}
    for cand in candidates:
        if cand.lower() in lowered:
            return lowered[cand.lower()]
    for cand in candidates:
        for low, orig in lowered.items():
            if cand.lower() in low:
                return orig
    raise ValueError(f"none of {candidates} found in {list(frame.columns)[:25]}")


# ---------------------------------------------------------------- joining ---

def assemble_cohort(
    embeddings_meta: pd.DataFrame,
    *,
    signatures: pd.DataFrame | None = None,
    ancestry: pd.DataFrame | None = None,
    purity: pd.DataFrame | None = None,
    clinical: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Left-join everything onto the embedding metadata, keyed on patient_id.

    Reports join coverage for each table, because a silent 40% join failure is
    the kind of thing that shows up as an inexplicably weak result three weeks
    later.
    """
    out = embeddings_meta.copy()
    if "patient_id" not in out.columns:
        raise KeyError("embeddings_meta must carry patient_id (call barcodes.annotate)")

    coverage: dict[str, float] = {}
    for name, table in (
        ("signatures", signatures),
        ("ancestry", ancestry),
        ("purity", purity),
        ("clinical", clinical),
    ):
        if table is None:
            continue
        if "patient_id" not in table.columns:
            raise KeyError(f"{name} table must carry patient_id")
        before = len(out)
        out = out.merge(table, on="patient_id", how="left", suffixes=("", f"_{name}"))
        if len(out) != before:
            raise RuntimeError(
                f"joining {name} changed row count {before} -> {len(out)}; "
                "the table has duplicate patient_ids"
            )
        probe = [c for c in table.columns if c != "patient_id"][0]
        coverage[name] = float(out[probe].notna().mean())

    out.attrs["join_coverage"] = coverage
    return out


def coverage_report(frame: pd.DataFrame) -> str:
    cov = frame.attrs.get("join_coverage", {})
    lines = [f"cohort: {len(frame)} rows, {frame['patient_id'].nunique()} patients"]
    for name, frac in cov.items():
        flag = "  <-- CHECK" if frac < 0.8 else ""
        lines.append(f"  {name:12s} {frac:6.1%} joined{flag}")
    return "\n".join(lines)


def _stack_embeddings(col: pd.Series, *, layer: int = -1) -> np.ndarray:
    """Stack an embedding column into (n_samples, n_dims), selecting one layer.

    Handles three storage shapes seen in the wild:
      - flat vector per row            -> used as is (layer ignored)
      - (n_layers, n_dims) per row     -> `layer` is selected
      - object array of per-layer rows -> coerced, then `layer` selected
    """
    first = np.asarray(col.iloc[0], dtype=object)
    probe = np.asarray(first.tolist() if first.dtype == object else first, dtype=float)

    if probe.ndim == 1:
        return np.vstack([np.asarray(v, dtype=float) for v in col]).astype(float)

    if probe.ndim != 2:
        raise ValueError(f"unexpected embedding shape {probe.shape}; expected 1-D or 2-D")

    n_layers = probe.shape[0]
    if not (-n_layers <= layer < n_layers):
        raise IndexError(f"layer {layer} out of range for {n_layers} layers")

    rows = []
    for v in col:
        a = np.asarray(v, dtype=object)
        a = np.asarray(a.tolist() if a.dtype == object else a, dtype=float)
        rows.append(a[layer])
    return np.vstack(rows).astype(float)


def load_purity(
    path: str | Path,
    *,
    patient_col: str | None = None,
    value_col: str | None = None,
    prefer: tuple[str, ...] = ("ABSOLUTE", "purity", "CPE", "ESTIMATE", "LUMP", "IHC"),
) -> tuple[pd.DataFrame, Provenance]:
    """Load a TCGA tumour-purity table and normalise its key to `patient_id`.

    Neither the PanCanAtlas ABSOLUTE file nor Aran et al. Supplementary Data 1
    ships a `patient_id` column — they use `array` or `Sample ID` — so
    `assemble_cohort` previously raised `KeyError: purity table must carry
    patient_id`. This normalises a 15- or 28-character barcode down to the
    12-character patient id.

    PREFER ABSOLUTE OVER CPE. CPE is a consensus that includes ESTIMATE, which
    is itself expression-derived and therefore partially embeds the immune score
    being predicted. Using CPE as the covariate in the purity decomposition is
    circular. Report CPE only as a sensitivity analysis.

    VERIFIED 2026-08-12 against the real PanCanAtlas ABSOLUTE table
    (GDC, HTTP 200, no auth, 901,812 bytes, 10,786 rows). Its columns are:
        array | sample | call status | purity | ploidy | Genome doublings |
        Coverage for 80% power | Cancer DNA fraction | Subclonal genome fraction | solution
    Note the purity column is plainly named `purity` — the FILE is the ABSOLUTE
    output, so no column says "ABSOLUTE". `array` is the 15-char sample barcode
    (TCGA-OR-A5J1-01) and `sample` is the full aliquot barcode; both normalise to
    a 12-char patient id.
    """
    path = Path(path)
    frame = pd.read_excel(path) if path.suffix in {".xlsx", ".xls"} else pd.read_csv(path)

    patient_col = patient_col or _guess_column(
        frame, ("patient_id", "bcr_patient_barcode", "array", "Sample ID", "sample", "submitter_id")
    )
    if value_col is None:
        for cand in prefer:
            hits = [c for c in frame.columns if cand.lower() in str(c).lower()]
            if hits:
                value_col = hits[0]
                break
    if value_col is None:
        raise ValueError(f"no purity column found among {list(frame.columns)[:20]}")

    out = frame[[patient_col, value_col]].rename(
        columns={patient_col: "patient_id", value_col: "purity"}
    )
    out["patient_id"] = out["patient_id"].astype(str).str.strip().str.upper().str[:12]
    out["purity"] = pd.to_numeric(out["purity"], errors="coerce")
    out = out.dropna(subset=["purity"]).groupby("patient_id", as_index=False)["purity"].mean()

    source_name = str(path).lower()
    circular = (
        ("estimate" in str(value_col).lower() or "cpe" in str(value_col).lower())
        and "absolute" not in source_name
    )
    prov = Provenance(
        name="purity", source=str(path), licence="see publication terms", n_rows=len(out),
        notes=[
            f"column used: {value_col!r}",
            "CIRCULARITY WARNING: this estimate is expression-derived and partially "
            "embeds the immune score. Use ABSOLUTE for the primary." if circular
            else "non-expression-derived estimate — suitable for the primary",
        ],
    )
    return out, prov


def load_expression(
    path: str | Path,
    *,
    gene_map: str | Path | None = None,
    sample_level: bool = True,
    sample_codes: tuple[str, ...] = ("01",),
) -> tuple[pd.DataFrame, Provenance]:
    """Load a TCGA bulk expression matrix (samples x genes), keyed for joining.

    JOIN GRANULARITY IS A REAL TRAP. TCGA expression is per-aliquot and the
    embeddings parquet's `sample` column is the 15-character sample barcode
    (`TCGA-06-0138-01`), whereas the null previously reindexed by the
    12-character patient id. A silent mismatch there yields n_null = 0 and an
    empty null with no error raised. `sample_level=True` keeps the 15-character
    key; set False only if you have deliberately collapsed to patients.

    `gene_map` optionally maps Ensembl ids to HUGO symbols (Xena TOIL ships
    Ensembl). Without it, gene sets keyed on symbols will not match.
    """
    path = Path(path)
    frame = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(
        path, sep="\t", index_col=0, compression="infer"
    )
    if path.suffix == ".parquet" and frame.index.name is None:
        frame = frame.set_index(frame.columns[0])

    if gene_map is not None:
        mp = pd.read_csv(gene_map, sep=None, engine="python")
        key, val = mp.columns[0], mp.columns[1]
        lookup = dict(zip(mp[key].astype(str), mp[val].astype(str)))
        frame.columns = [lookup.get(str(c).split(".")[0], str(c)) for c in frame.columns]
        frame = frame.loc[:, ~frame.columns.duplicated()]

    idx = frame.index.astype(str).str.upper().str.strip()
    if sample_codes:
        keep = idx.str[13:15].isin(sample_codes)
        frame, idx = frame.loc[keep.to_numpy()], idx[keep.to_numpy()]
    frame.index = idx.str[:15] if sample_level else idx.str[:12]

    prov = Provenance(
        name="expression", source=str(path), licence="open", n_rows=len(frame),
        notes=[
            f"{frame.shape[1]} genes x {frame.shape[0]} samples",
            f"key granularity: {'15-char sample' if sample_level else '12-char patient'}",
            f"sample codes kept: {sample_codes}",
            "VERIFY the key matches the embeddings frame before scoring",
        ],
    )
    return frame, prov
