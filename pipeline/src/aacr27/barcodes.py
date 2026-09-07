"""TCGA barcode parsing.

The tissue source site (TSS) is the second field of a TCGA barcode and is the
confounder at the centre of this study. Getting this parse right matters more
than almost anything else in the pipeline: every split, every negative control
and every variance-decomposition term keys off it.

    TCGA-02-0001-01C-01D-0182-01
         ^^        ^^
         |         └── sample type ("01" = primary solid tumour)
         └── tissue source site

Reference: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

# Sample-type codes we treat as primary tumour. 01 = Primary Solid Tumor,
# 03 = Primary Blood Derived Cancer (heme, kept out by default).
PRIMARY_TUMOUR_CODES = frozenset({"01"})
NORMAL_CODES = frozenset({"10", "11", "12", "13", "14"})

_BARCODE_RE = re.compile(
    r"^TCGA-(?P<tss>[0-9A-Z]{2})-(?P<participant>[0-9A-Z]{4})"
    r"(?:-(?P<sample>[0-9]{2})(?P<vial>[A-Z])?)?"
    r"(?:-(?P<portion>[0-9]{2})(?P<analyte>[A-Z])?)?"
    r"(?:-(?P<plate>[0-9A-Z]{4}))?"
    r"(?:-(?P<centre>[0-9]{2}))?"
)


@dataclass(frozen=True)
class Barcode:
    """A parsed TCGA barcode. Fields beyond `patient_id` may be None."""

    full: str
    tss: str
    participant: str
    sample: str | None = None
    vial: str | None = None
    portion: str | None = None
    analyte: str | None = None
    plate: str | None = None
    centre: str | None = None

    @property
    def patient_id(self) -> str:
        """TCGA-XX-YYYY. The unit of analysis — never split on anything finer."""
        return f"TCGA-{self.tss}-{self.participant}"

    @property
    def is_primary_tumour(self) -> bool:
        return self.sample in PRIMARY_TUMOUR_CODES

    @property
    def is_normal(self) -> bool:
        return self.sample in NORMAL_CODES


def parse(barcode: str) -> Barcode:
    """Parse a TCGA barcode of any depth. Raises ValueError if malformed."""
    if not isinstance(barcode, str):
        raise ValueError(f"barcode must be a string, got {type(barcode)!r}")
    match = _BARCODE_RE.match(barcode.strip().upper())
    if match is None:
        raise ValueError(f"not a TCGA barcode: {barcode!r}")
    return Barcode(full=barcode.strip().upper(), **match.groupdict())


def try_parse(barcode: str) -> Barcode | None:
    """Non-raising variant, for use over columns of mixed provenance."""
    try:
        return parse(barcode)
    except ValueError:
        return None


def annotate(
    frame: pd.DataFrame,
    barcode_col: str = "barcode",
    *,
    drop_unparseable: bool = True,
) -> pd.DataFrame:
    """Add patient_id / tss / sample_type columns derived from a barcode column.

    Anything unparseable is dropped by default and the count reported, rather
    than silently coerced — a silent coercion here would corrupt the site
    analysis, which is the whole point of the study.
    """
    if barcode_col not in frame.columns:
        raise KeyError(f"{barcode_col!r} not in frame; have {list(frame.columns)[:12]}")

    parsed = frame[barcode_col].map(try_parse)
    bad = int(parsed.isna().sum())
    if bad and not drop_unparseable:
        raise ValueError(f"{bad} unparseable barcodes and drop_unparseable=False")

    out = frame.loc[parsed.notna()].copy()
    good = parsed.dropna()
    out["patient_id"] = [b.patient_id for b in good]
    out["tss"] = [b.tss for b in good]
    out["sample_type"] = [b.sample for b in good]
    out["is_primary_tumour"] = [b.is_primary_tumour for b in good]
    # Plate (barcode positions 22-25) is pure processing batch with no plausible
    # biological reading, which makes plate-within-site the technical-only lower
    # bound on site variance. See decomposition.label_site_variance.
    out["plate"] = [b.plate for b in good]

    if bad:
        out.attrs["dropped_unparseable"] = bad
    return out


def site_summary(frame: pd.DataFrame, tss_col: str = "tss") -> pd.DataFrame:
    """Per-site patient counts — read this before choosing a fold count.

    Preserved-site CV assigns whole sites to folds, so the number of usable
    folds is bounded by how evenly patients are distributed across sites. A
    cohort where one site holds 40% of patients cannot support 5 balanced
    site-disjoint folds.
    """
    if "patient_id" not in frame.columns:
        raise KeyError("call annotate() first — patient_id required")
    counts = (
        frame.groupby(tss_col)["patient_id"]
        .nunique()
        # `kind="stable"` is LOAD-BEARING even though nothing calls this
        # function today. Site sizes tie heavily -- 51 tied counts in the
        # pan-TCGA cohort -- and pandas' default sort is quicksort, which is
        # not stable, so tied sites came back in an order that depends on the
        # platform's introsort. That is precisely the A9/A10 defect, and
        # `26_sort_audit.py` reported this call as its last DANGEROUS case
        # until 2026-09-06. Fixed here rather than left for whoever calls it
        # next, because a latent instance of a known defect class is still an
        # instance. `cum_frac` below is read in fold-count decisions, so the
        # row order is not cosmetic.
        .sort_values(ascending=False, kind="stable")
        .rename("n_patients")
        .to_frame()
    )
    total = counts["n_patients"].sum()
    counts["frac"] = counts["n_patients"] / total
    counts["cum_frac"] = counts["frac"].cumsum()
    return counts
