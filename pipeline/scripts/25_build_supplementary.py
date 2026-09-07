#!/usr/bin/env python3
"""Render the promised supplementary tables from the frozen artefacts.

    python3 scripts/25_build_supplementary.py            # render
    python3 scripts/25_build_supplementary.py --check     # verify, render nothing

WHY THIS IS A SCRIPT AND NOT SEVEN HAND-MADE FILES
==================================================
`18-SUPPLEMENTARY-INVENTORY.md` records seven promised supplementary items with
their numbers frozen on disk and **zero** formatted deliverables. Typing those
tables out by hand would create seven new places for a number to drift away from
`results/`, and this project has already been burned by exactly that failure --
which is the whole reason `19_check_numbers.py` exists. Rendering them from the
artefacts creates none.

`--check` goes further: it re-reads the rendered files and compares a named value
in each against `19_check_numbers.py`'s `authority()`, so the supplement cannot
silently disagree with the manuscript. It is designed to be able to fail --
corrupt a rendered TSV and it goes red.

FORMAT
======
One TSV per item plus `MANIFEST.md`. TSV rather than XLSX deliberately: *Cancer
Research*'s supplementary packaging rules are **unverified** -- `aacrjournals.org`
returns HTTP 403 to an automated request, and that block is recorded in
`16-YOUR-TASKS.md`. A plain delimited table converts to whatever the journal
turns out to want; a styled workbook does not convert back. Revisit once the
Instructions for Authors are in the repo.

PROVENANCE
==========
Every table names its source artefact in `MANIFEST.md`, and every source is
hashed in `results/PROVENANCE.json`, so a reader can verify the chain from a
rendered cell back to the run that produced it.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
OUT = RESULTS / "supplementary"
sys.path.insert(0, str(REPO / "src"))


def _authority() -> dict:
    """Import 19_check_numbers.py by path -- its name is not an identifier."""
    spec = importlib.util.spec_from_file_location(
        "_chk19", Path(__file__).parent / "19_check_numbers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.authority()


def _add_bh(frame: pd.DataFrame) -> pd.DataFrame:
    """Add p_two_sided / q_bh / beats_null_bh using the pipeline's OWN function.

    Imported rather than reimplemented on purpose: `11_close_science_gaps.py`
    `add_bh` is what produced pan-cancer's stored columns, so deriving NSCLC's
    with the same code guarantees one definition of p across the table. A local
    copy could drift from it silently, which is the whole failure mode S3's
    caption was worried about.
    """
    spec = importlib.util.spec_from_file_location(
        "_gaps11", Path(__file__).parent / "11_close_science_gaps.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.add_bh(frame.copy(), "ISI")


def _csv(rel: str) -> pd.DataFrame:
    p = RESULTS / rel
    if not p.exists():
        raise SystemExit(f"MISSING SOURCE {p.relative_to(REPO)} -- refusing to "
                         "render a supplement from a file that is not there.")
    return pd.read_csv(p)


def _json(rel: str) -> dict:
    p = RESULTS / rel
    if not p.exists():
        raise SystemExit(f"MISSING SOURCE {p.relative_to(REPO)}")
    return json.loads(p.read_text())


# ---------------------------------------------------------------- the tables
# Each builder returns (dataframe, one-line description, source list).

def s1_ancestry():
    j = _json("ancestry/summary.json")
    rows = [{"quantity": "patients analysed", "value": j["n_patients"]}]
    for grp, val in j["counts_patient_level"].items():
        rows.append({"quantity": f"patients, {grp}", "value": val})
    rows += [
        {"quantity": "Cramer's V, ancestry vs tissue source site",
         "value": j["site_cramers_v"]["v"]},
        {"quantity": "Cramer's V, permutation 95th percentile",
         "value": j["site_cramers_v"]["v_perm_p95"]},
    ]
    return (pd.DataFrame(rows),
            "Ancestry arm: composition and the site-confounding check. The "
            "manuscript calls this a supplementary table by name, so it is S1.",
            ["ancestry/summary.json"])


def s2_small_panel():
    d = _csv("pancancer_v3/small_panel_bracket.csv")
    return (d,
            "Curated-versus-random reliability gap by panel size, pan-TCGA. The "
            "manuscript quotes only the endpoints (0.084 at k=160, 0.453 at "
            "k=10); this is the full bracket.",
            ["pancancer_v3/small_panel_bracket.csv"])


def s3_per_signature_excess():
    # The asymmetry this used to report is GONE, and the reason it survived so
    # long is worth stating. NSCLC's stored immune_excess.csv predates the
    # BH-FDR columns, so the table left them empty and the caption said they
    # "cannot be recovered without re-running the cohort". That was wrong.
    # Measured 2026-09-06: pan-cancer's STORED p_two_sided is exactly the normal
    # approximation from its own excess_z/lo/hi -- recomputing it reproduces the
    # stored column to 3.9e-121, i.e. bit-for-bit. There was never a second
    # definition of p to collide with, so the identical function applied to
    # NSCLC's stored excess_z/lo/hi yields a p in the SAME definition, with no
    # re-run and no imputation.
    #
    # Derived HERE rather than written back: results/nsclc_v3/ is frozen and
    # hashed in PROVENANCE.json, so the supplement computes the columns on the
    # fly and the frozen artefact is never touched.
    base = ["signature", "excess_z", "excess_lo", "excess_hi", "beats_null"]
    optional = ["p_two_sided", "q_bh", "beats_null_bh"]
    frames, derived = [], []
    for cohort, label in (("pancancer_v3", "pan-TCGA"), ("nsclc_v3", "NSCLC")):
        d = _csv(f"{cohort}/immune_excess.csv")
        missing = [c for c in optional if c not in d.columns]
        if missing:
            d = _add_bh(d)
            still = [c for c in optional if c not in d.columns]
            assert not still, f"{label}: {still} still absent after derivation"
            derived.append(f"{label} ({', '.join(missing)})")
        d = d[base + optional].copy()
        d.insert(0, "cohort", label)
        frames.append(d)
    note = ("Per-signature immune-specificity excess. The manuscript quotes "
            "only the RANGE (0.148-0.440 NSCLC); all 16 rows per cohort "
            "belong in the supplement.")
    if derived:
        note += (" MULTIPLICITY COLUMNS DERIVED for " + "; ".join(derived)
                 + ": that run predates the columns, so they are recomputed "
                 "from the stored excess and its interval by the same normal "
                 "approximation that produced the other cohort's stored "
                 "values. No cell is imputed and no frozen file is modified.")
    return (pd.concat(frames, ignore_index=True), note,
            ["pancancer_v3/immune_excess.csv", "nsclc_v3/immune_excess.csv"])


def s4_reliability():
    frames = []
    for cohort, label in (("pancancer_v3", "pan-TCGA"), ("nsclc_v3", "NSCLC")):
        d = _csv(f"{cohort}/immune_excess.csv")[
            ["signature", "alpha_observed", "alpha_null_mean",
             "alpha_observed_raw", "alpha_null_mean_raw",
             "reliability_gap"]].copy()
        d.insert(0, "cohort", label)
        frames.append(d)
    return (pd.concat(frames, ignore_index=True),
            "Cronbach's alpha and the curated-versus-random reliability gap, "
            "per signature, residualised and raw. NOTE ON CONVENTION: the "
            "manuscript quotes alpha levels as MEDIANS across signatures and "
            "reliability gaps as MEANS. That mix was recorded only in a "
            "handoff document until now; it is stated here where a reader "
            "meets the numbers.",
            ["pancancer_v3/immune_excess.csv", "nsclc_v3/immune_excess.csv"])


def s5_negative_controls():
    frames = []
    for cohort, label in (("pancancer_v3", "pan-TCGA"), ("nsclc_v3", "NSCLC")):
        d = _csv(f"{cohort}/site_control.csv").copy()
        d.insert(0, "cohort", label)
        d.insert(1, "control", "site AUROC from the embedding")
        frames.append(d)
    combat = _csv("nsclc_v3/site_control_combat.csv").copy()
    combat.insert(0, "cohort", "NSCLC")
    combat.insert(1, "control", "site AUROC after ComBat")
    frames.append(combat)
    return (pd.concat(frames, ignore_index=True),
            "Negative controls: per-site AUROC for recovering tissue source "
            "site from the embedding, plus the ComBat arm. Figure 4 shows two "
            "of these; the ComBat estimability claim had no table at all.",
            ["pancancer_v3/site_control.csv", "nsclc_v3/site_control.csv",
             "nsclc_v3/site_control_combat.csv"])


def s6_label_site_variance():
    d = _csv("pancancer_v3/label_site_variance.csv")
    return (d,
            "Control C, label-side site variance: each signature's ground "
            "truth regressed on site FROM RNA ALONE, no image. Prose-only in "
            "the manuscript until now.",
            ["pancancer_v3/label_site_variance.csv"])


def s7_decomposition():
    frames = []
    for cohort, label in (("pancancer_v3", "pan-TCGA"), ("nsclc_v3", "NSCLC")):
        d = _csv(f"{cohort}/decomposition.csv").copy()
        d.insert(0, "cohort", label)
        frames.append(d)
    return (pd.concat(frames, ignore_index=True),
            "Nested-model decomposition: partial correlation and incremental "
            "image R-squared over purity, type and site.",
            ["pancancer_v3/decomposition.csv", "nsclc_v3/decomposition.csv"])


def s8_partition_variance():
    frames = []
    for cohort, label in (("nsclc", "NSCLC"), ("pancancer", "pan-TCGA")):
        p = RESULTS / f"partition_variance_{cohort}" / "per_partition.csv"
        if not p.exists():
            continue
        d = pd.read_csv(p).copy()
        d.insert(0, "cohort", label)
        frames.append(d)
    if not frames:
        raise SystemExit("no partition-variance artefacts at all")
    note = ("Per-partition ISI under varying site-to-fold assignments "
            "(split_seed only; every other random choice held fixed).")
    if len(frames) == 1:
        note += (" PAN-CANCER IS ABSENT: it had not finished when this was "
                 "rendered. Re-run this script once it has.")
    return (pd.concat(frames, ignore_index=True), note,
            [f"partition_variance_{c}/per_partition.csv"
             for c in ("nsclc", "pancancer")])


def s9_platform():
    rows = []
    for tag, fname in (("macOS arm64", "split_determinism_darwin_arm64.json"),
                       ("Linux x86_64", "split_determinism_linux_x86_64.json")):
        p = RESULTS / fname
        if not p.exists():
            continue
        j = json.loads(p.read_text())
        rows.append({
            "platform": tag,
            "python": j["platform"]["python"],
            "numpy": j["platform"]["numpy"],
            "sites": j["n_sites"],
            "sites_tied_on_size": j["n_sites_tied"],
            "hash_sizes": j["hash_sizes"],
            "hash_fold_shipped": j["hash_fold_shipped"],
            "hash_fold_stable": j["hash_fold_stable"],
            "synthetic_isi_shipped": j["synthetic_isi_shipped"],
            "synthetic_isi_stable": j["synthetic_isi_stable"],
        })
    if not rows:
        raise SystemExit("no split-determinism artefacts")
    return (pd.DataFrame(rows),
            "A9, cross-platform reproducibility. Inputs hash identically on "
            "both machines; the SHIPPED fold assignment does not, because "
            "np.argsort's default sort is not stable and most sites are tied "
            "on size. Under a stable sort both machines agree exactly. "
            "Supports Limitations 8 directly.",
            ["split_determinism_darwin_arm64.json",
             "split_determinism_linux_x86_64.json"])


TABLES = [
    ("S1_ancestry", s1_ancestry),
    ("S2_small_panel_bracket", s2_small_panel),
    ("S3_per_signature_excess", s3_per_signature_excess),
    ("S4_reliability_alpha", s4_reliability),
    ("S5_negative_controls", s5_negative_controls),
    ("S6_label_site_variance", s6_label_site_variance),
    ("S7_decomposition", s7_decomposition),
    ("S8_partition_variance", s8_partition_variance),
    ("S9_platform_reproducibility", s9_platform),
]

# Journal-style captions. MANIFEST.md's one-liners are internal notes -- they
# say things like "prose-only in the manuscript until now", which is provenance
# for us and noise for a reviewer. A journal caption has to be a titled,
# numbered sentence that STANDS ALONE: a reader who sees only the table and its
# caption, with no manuscript in front of them, must be able to say what the
# table shows, what the columns mean, and where the numbers came from.
#
# Written here rather than in a separate document so they regenerate with the
# tables and cannot drift out of step with them. Where a table has a defect a
# reviewer would otherwise discover for themselves -- S3's missing NSCLC
# multiplicity columns, S8's absent pan-cancer arm -- the caption states it.
# Concealing a gap in a supplement is how a reviewer stops trusting the rest.
CAPTIONS = {
    "S1_ancestry": (
        "Genetic-ancestry composition of the pan-TCGA cohort and the "
        "site-confounding check.",
        "Counts and proportions by inferred ancestry group, with Cramer's V "
        "for the association between ancestry and tissue source site together "
        "with its permutation 95th percentile. The comparison of V against its "
        "permutation null is the identifiability result: ancestry and "
        "collection site are not separable in this cohort, so no "
        "ancestry-stratified claim is made anywhere in the manuscript."),
    "S2_small_panel_bracket": (
        "Curated-versus-random reliability gap as a function of gene-panel "
        "size, pan-TCGA.",
        "Signatures were subsampled to fixed panel sizes and the "
        "curated-minus-random difference in Cronbach's alpha recomputed at "
        "each. The manuscript quotes only the endpoints; the full bracket is "
        "given here. The gap grows monotonically as panels shrink "
        "(Spearman rho = -1.00), so the smallest panels are the most "
        "misleading, not the least."),
    "S3_per_signature_excess": (
        "Immune-specificity excess for each of the 16 Hallmark signatures, "
        "both cohorts.",
        "Residualized image-signature correlation, the mean of its size- and "
        "expression-matched random-set null, the disattenuated excess in "
        "Fisher z with its patient-clustered bootstrap interval, and whether "
        "the excess exceeds the null. The manuscript quotes only the range. "
        "NOTE ON THE MULTIPLICITY COLUMNS: pan-TCGA's p_two_sided, q_bh and "
        "beats_null_bh are read from its stored results; NSCLC's run predates "
        "those columns, so they are DERIVED here from its stored excess and "
        "bootstrap interval. This is not an imputation and the two cohorts do "
        "not carry two definitions of p: pan-TCGA's stored p is itself the "
        "normal approximation from its interval, which recomputing reproduces "
        "to within 4e-121, so the same function applied to NSCLC yields the "
        "same quantity. NSCLC's frozen file is read, never modified. Both "
        "cohorts give 16 of 16 signatures beating their null after "
        "Benjamini-Hochberg, so the correction changes no conclusion."),
    "S4_reliability_alpha": (
        "Per-signature reliability, residualized and raw, both cohorts.",
        "Cronbach's alpha for each curated signature, the mean alpha of its "
        "matched random sets, and the difference. Reported on both the "
        "residualized scores the index uses and the raw scores. NOTE ON "
        "CONVENTION: the manuscript quotes alpha levels as MEDIANS across "
        "signatures and reliability gaps as MEANS; both are given per "
        "signature here so either can be recomputed."),
    "S5_negative_controls": (
        "Negative controls: recoverability of tissue source site from the "
        "image embedding, with and without batch correction.",
        "One-versus-rest AUROC for predicting collection site from the "
        "embedding, per site, in both cohorts, plus the ComBat-corrected NSCLC "
        "arm. Site is almost perfectly recoverable (median AUROC 0.998 "
        "pan-TCGA, 0.992 NSCLC), which is why preserved-site cross-validation "
        "is used throughout. The ComBat arm is tabulated so that the "
        "post-correction AUROCs can be seen rather than described, and it is "
        "NOT offered as a corrected confounding estimate: per-site centring "
        "imposes a within-site zero-sum constraint that drives the site "
        "classifier below chance (median 0.005), so the corrected number is "
        "uninterpretable in either direction and the manuscript quotes the raw "
        "figure instead. Separately, and this is the deployment-relevant "
        "point rather than a property of this table, ComBat is estimable for "
        "0% of held-out samples under site-disjoint folds by construction, "
        "which is the position of a model meeting a new hospital."),
    "S6_label_site_variance": (
        "Control C: site-attributable variance on the LABEL side, from RNA "
        "alone.",
        "Each signature's ground-truth score regressed on tissue source site "
        "with no image involved, so the quantity measured is how much of the "
        "target is site structure before any prediction is attempted. "
        "Permutation-calibrated. Referred to in the manuscript's prose but not "
        "previously tabulated."),
    "S7_decomposition": (
        "Nested-model decomposition of the image signal over purity, cancer "
        "type and site.",
        "Partial correlations and incremental image R-squared as covariates "
        "are added in a fixed order, both cohorts. This is the basis for the "
        "statement that pan-cancer the embedding adds nothing over covariates."),
    "S8_partition_variance": (
        "Sensitivity of the index to the choice of site-to-fold partition.",
        "The index recomputed under repeated site-to-fold assignments, varying "
        "only split_seed and holding every other random choice fixed, so the "
        "spread is attributable to the partition alone. The between-partition "
        "standard deviation is the quantity the manuscript compares against "
        "the bootstrap standard error. Both cohorts are present as of "
        "2026-09-05. Six assignments were requested per cohort; NSCLC yielded "
        "six and pan-TCGA five, because one pan-TCGA assignment was refused by "
        "the degeneracy guard when the residualization SVD failed to converge. "
        "That refusal is reported rather than replaced: the estimand is not "
        "computable for every valid partition, and substituting another seed "
        "would have concealed it. Runtime per partition is included because "
        "the pan-TCGA figures (3,364-11,404 s) are the first measured for that "
        "cohort and correct an estimate that three handoffs carried."),
    "S9_platform_reproducibility": (
        "Cross-platform reproducibility of the fold assignment, and the "
        "non-stable sort responsible.",
        "Hashes of the inputs and of the resulting site-to-fold partition, "
        "computed independently on macOS/arm64 and Linux/x86_64. The inputs "
        "hash identically on both machines; the partition produced by the "
        "originally shipped code does not, because numpy's default sort is not "
        "stable and most tissue source sites share a patient count with "
        "another site. Under a stable sort both machines agree exactly. This "
        "table supports Limitation 8 directly, and the same defect was "
        "subsequently found on the single-sample GSEA rank table."),
}

# (table, column, row-selector, authority key, tolerance). Checked BY RE-READING
# the rendered file, not the source -- otherwise this would only prove pandas
# can round-trip a CSV.
CHECKS = [
    ("S1_ancestry", "value", ("quantity", "patients analysed"),
     "ancestry_n", 0.5),
    ("S1_ancestry", "value",
     ("quantity", "Cramer's V, ancestry vs tissue source site"),
     "ancestry_cramers_v", 5e-4),
    ("S3_per_signature_excess", "excess_z", None, "pancancer_isi", 5e-3),
    ("S9_platform_reproducibility", "synthetic_isi_shipped",
     ("platform", "macOS arm64"), "splitdet_macos_isi_shipped", 1e-12),
    ("S9_platform_reproducibility", "synthetic_isi_stable",
     ("platform", "Linux x86_64"), "splitdet_hpc4_isi_stable", 1e-12),
]


def run_checks(table: dict) -> int:
    bad = 0
    for name, col, sel, key, tol in CHECKS:
        p = OUT / f"{name}.tsv"
        if not p.exists():
            print(f"  FAIL {name}: not rendered")
            bad += 1
            continue
        d = pd.read_csv(p, sep="\t")
        if sel is None:
            # The mean over the pan-TCGA rows must reproduce the headline ISI.
            got = float(d[d["cohort"] == "pan-TCGA"][col].mean())
            what = f"mean {col} over pan-TCGA rows"
        else:
            scol, sval = sel
            hit = d[d[scol].astype(str) == sval]
            if hit.empty:
                print(f"  FAIL {name}: no row where {scol} == {sval!r}")
                bad += 1
                continue
            got = float(hit.iloc[0][col])
            what = f"{col} where {scol}={sval!r}"
        want = table.get(key)
        if want is None:
            print(f"  FAIL {name}: authority has no key {key!r}")
            bad += 1
        elif abs(got - want) > tol:
            print(f"  FAIL {name}: {what} = {got!r}, authority {key} = "
                  f"{want!r} (tol {tol})")
            bad += 1
        else:
            print(f"  ok   {name}: {what} agrees with {key}")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="verify the rendered tables against authority(); "
                         "render nothing, exit 1 on disagreement")
    args = ap.parse_args()

    if args.check:
        print(f"checking {len(CHECKS)} values in {OUT.relative_to(REPO)} "
              "against 19_check_numbers.py's authority()\n")
        bad = run_checks(_authority())
        print()
        if bad:
            print(f"FAILED: {bad} disagreement(s) between the supplement and "
                  "the frozen results.")
            return 1
        print("OK: every checked supplementary value agrees with the "
              "manuscript's authority.")
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    manifest = ["# Supplementary tables — rendered, do not hand-edit", "",
                "Generated by `pipeline/scripts/25_build_supplementary.py` from "
                "the frozen artefacts in `pipeline/results/`. Re-run it rather "
                "than editing a cell; `--check` verifies these files against "
                "the same authority table the manuscript is checked against.",
                "",
                "Format is TSV because *Cancer Research*'s supplementary "
                "packaging rules are still unverified (HTTP 403 to an "
                "automated request). See `16-YOUR-TASKS.md`.", ""]
    # A table without a caption is not submittable, so this is an assertion and
    # not a warning: adding a table to TABLES and forgetting CAPTIONS should
    # stop the build, not produce a silently incomplete package.
    uncaptioned = [n for n, _ in TABLES if n not in CAPTIONS]
    assert not uncaptioned, f"tables with no caption in CAPTIONS: {uncaptioned}"

    captions = [
        "# Supplementary table captions",
        "",
        "Generated by `pipeline/scripts/25_build_supplementary.py` alongside "
        "the tables themselves, so a caption cannot drift out of step with the "
        "table it describes. Each is written to stand alone: a reader with the "
        "table and this caption, and no manuscript, should be able to say what "
        "is shown and where it came from.",
        "",
        "Paste these into the submission package in this order. Numbering "
        "follows the S-numbers used in `09-PAPER-DRAFT.md`.",
        "",
    ]
    for name, fn in TABLES:
        df, desc, sources = fn()
        path = OUT / f"{name}.tsv"
        df.to_csv(path, sep="\t", index=False)
        print(f"  {name:32s} {len(df):>4d} rows x {len(df.columns):>2d} cols")
        num = name.split("_")[0]
        title, body = CAPTIONS[name]
        manifest += [f"## {name}", "", desc, "",
                     f"- rows: {len(df)}, columns: {len(df.columns)}",
                     "- source(s): "
                     + ", ".join(f"`results/{s}`" for s in sources),
                     f"- caption: see `CAPTIONS.md` ({num})", ""]
        captions += [f"**Supplementary Table {num}.** {title}", "", body, "",
                     f"*{len(df)} rows x {len(df.columns)} columns. "
                     + "Source: "
                     + ", ".join(f"`{s}`" for s in sources) + ".*", "", "---", ""]
    (OUT / "MANIFEST.md").write_text("\n".join(manifest))
    (OUT / "CAPTIONS.md").write_text("\n".join(captions))
    print(f"\nwrote {len(TABLES)} tables + MANIFEST.md + CAPTIONS.md to "
          f"{OUT.relative_to(REPO)}")
    print("run with --check to verify them against the frozen results")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
