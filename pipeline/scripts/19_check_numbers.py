#!/usr/bin/env python3
"""Fail if any document quotes a number the frozen results do not support.

    python scripts/19_check_numbers.py             # check, exit 1 on mismatch
    python scripts/19_check_numbers.py --list      # print the authority table
    python scripts/19_check_numbers.py --verbose   # show every claim, not just misses
    python scripts/19_check_numbers.py --self-test # prove the check can still FAIL

WHY THIS EXISTS
===============
This project has repeatedly shipped documents whose numbers were REMEMBERED
rather than re-measured. The test count and the abstract character count have
both drifted. Every handoff since #17 has carried a rule about it, and a rule is
not a fix: it depends on a human noticing, every time, forever.

So this checks mechanically. It parses the frozen artefacts into a table of
authoritative values, then scans the markdown for claims tied to those values
and reports every disagreement with a file and line number.

WHY --self-test EXISTS
======================
Version 1 of this script reported "checked 62 claims ... OK" while checking
nothing at all: its regexes matched only the CORRECT literal, so they passed by
construction. That was found by hand, by deliberately injecting a wrong number.
A green check nobody has ever seen go red is not evidence.

`--self-test` makes that permanent. It copies the documents to a temp tree,
substitutes a deliberately wrong value for each family of claim, and asserts the
checker reports a MISMATCH for it. A pattern that has gone vacuous fails the
self-test instead of passing silently. Run it whenever you add a claim.

WHAT IT DELIBERATELY DOES NOT DO
================================
It does not try to parse free prose for arbitrary numeric assertions. A
general-purpose parser would be wrong often enough that its output would be
ignored, which is worse than no checker. Instead it holds a CURATED list of
(regex -> authoritative key) pairs, and it PRINTS ITS OWN COVERAGE so the
limitation is visible rather than implied: "checked N claims across M files".

A claim not on the list is not checked. That is honest, and adding one is a
two-line edit to the tables below.

WHAT IS DELIBERATELY EXEMPT
===========================
Three kinds of number are historical records, and "correcting" them would
falsify the record. They are skipped by file or by rule, and the reasons are
stated at EXEMPT below rather than left as tribal knowledge:

  * `05-PRE-REGISTRATION.md`'s test count (67) is frozen at what was true when
    the protocol was tagged.
  * `02-PROJECT-DECISION.md` records what was true when it was written.
  * `results/scorer_sensitivity_hpc4/` numbers are diagnostic and not quotable,
    so a document must NOT match them.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
RESULTS = ROOT / "results"
BUILDER = ROOT / "scripts" / "22_build_public_snapshot.py"


def _is_deposit() -> bool:
    """True when this tree is the published snapshot, not the working repo.

    Same marker `21_provenance_manifest.py` uses: the snapshot builder excludes
    itself by design -- its docstring quotes every string it forbids -- so its
    ABSENCE is what identifies a deposit."""
    return not BUILDER.is_file()


def _missing_artefact_message(missing: list[str]) -> str:
    """Refuse, but name the RIGHT reason.

    This check refuses either way: it cannot verify a full-precision claim whose
    authority is not on disk, and quietly dropping the authority would leave the
    claim looking checked when it is not. What changes here is only the
    diagnosis, and only because the wrong one is expensive.

    `bisect_macos.npz` and `bisect_hpc4.npz` are 2.18 MB EACH and are excluded
    from the public snapshot on purpose. Measured 2026-09-07 inside a freshly
    staged tree: this refusal fired, reading "FROZEN ARTEFACT MISSING", which
    tells a reader of an intact deposit that it has been tampered with. That is
    the same liability `21_provenance_manifest.py` carried until it learned to
    say NOT DEPOSITED -- a check that prints a claim it cannot support is worse
    than no check, because the next reader cannot tell it from real corruption.
    """
    head = ("FROZEN ARTEFACT MISSING -- refusing to check against remembered "
            "constants:\n  " + "\n  ".join(missing))
    if not _is_deposit():
        return head
    return (head + "\n\n"
            "This tree is a PUBLIC DEPOSIT, not the working repository, and "
            "the files above\nare excluded from it deliberately -- they are "
            "multi-megabyte internal platform\nprobes, not results the paper "
            "reports. NOTHING IS WRONG WITH THIS DEPOSIT.\n"
            "This script verifies the working repository's documents against "
            "its own frozen\nartefacts; it is not one of the commands that "
            "runs from the deposit alone.\nSee `pipeline/README.md`, "
            "'What runs from the deposit alone'.")

# Files that record history and must not be retro-corrected. Listed with the
# reason, because an unexplained exemption becomes a place to hide drift.
#
# EVERY KEY HERE MUST ALSO APPEAR IN `DOCS`, and a test enforces that. Until
# 2026-09-07 neither key was in DOCS, so this table was unreachable: the files
# were unscanned because nothing listed them, not because anything exempted
# them. The net effect was identical, which is exactly why it survived -- the
# docstring section above described a mechanism that had never once executed,
# and handoffs carried "in the checker's EXEMPT table" forward as a live fact.
# An exemption that cannot fire is not an exemption; it is a comment.
EXEMPT = {
    "05-PRE-REGISTRATION.md":
        "test count frozen at 67, deliberately, at the protocol's git tag",
    "02-PROJECT-DECISION.md":
        "historical record; its numbers state what was true when written",
}


def _dig(obj, key):
    """First occurrence of `key` anywhere in a nested JSON structure."""
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for v in obj.values():
            got = _dig(v, key)
            if got is not None:
                return got
    elif isinstance(obj, list):
        for v in obj:
            got = _dig(v, key)
            if got is not None:
                return got
    return None


def _rows(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _median(values) -> float:
    return float(statistics.median(values))


def _mde_c_index(rows: list[dict]) -> float:
    """Median minimum detectable effect in C-index, 80% power, alpha 0.05.

    Derived rather than read, because `nsclc_v3/outcome_arm.csv` predates the
    `mde_80pct_c_index` column and does not carry it. The derivation is the one
    in `10_make_figures.py`: recover the standard error from the reported CI
    half-width, scale it for 80% power, then convert Fisher z to C-index via
    Somers' D (C = (tanh z)/2 + 1/2, so a difference of z is (tanh z)/2 in C).

    This is NOT asserted -- it is cross-checked. For pan-cancer, where the
    column DOES exist, `authority()` asserts the derived value reproduces the
    stored one. If the derivation ever stops matching, the run aborts.
    """
    import math

    from scipy import stats as sps

    crit = float(sps.norm.ppf(0.975))
    power = float(sps.norm.ppf(0.80))
    ax = [r for r in rows if r["adjustment"] == "axis_residualised"]
    if not ax:
        raise SystemExit("outcome_arm.csv has no axis_residualised rows")
    se = [(float(r["excess_hi"]) - float(r["excess_lo"])) / (2 * crit) for r in ax]
    z = _median([(crit + power) * s for s in se])
    return math.tanh(z) / 2


def _hallmark_categories() -> dict[str, str] | None:
    """MSigDB process categories, read from `10_make_figures.py`'s CAT literal.

    SINGLE SOURCE OF TRUTH, deliberately. The immune-versus-other contrast in
    the manuscript depends entirely on which signatures count as immune, and a
    second copy of that mapping here would be a second place for it to drift --
    exactly the failure this file exists to prevent.

    The literal is extracted with `ast`, NOT by importing the module. Importing
    would pull in matplotlib for a dictionary, and matplotlib needs a writable
    config directory that is not guaranteed here. Parsing the source cannot run
    anything and cannot fail for environmental reasons.
    """
    import ast

    src = ROOT / "scripts" / "10_make_figures.py"
    if not src.exists():
        return None
    try:
        tree = ast.parse(src.read_text())
    except SyntaxError:
        return None
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == "CAT":
                try:
                    val = ast.literal_eval(node.value)
                except ValueError:
                    return None
                return val if isinstance(val, dict) else None
    return None


def _immune_contrast(rows: list[dict]) -> dict[str, float] | None:
    """The immune-versus-other outcome contrast, recomputed from the artefact.

    The manuscript reports "(Mann-Whitney p = 0.016; Fisher p = 0.093)" and the
    Figure 3 caption reports "p = 0.016 -> 0.008" for the counterfactual in
    which IL6/JAK/STAT3 is counted as immune. None of those four numbers lives
    in any frozen artefact -- they were computed once and typed in -- so until
    now the paper's only inferential p-value for its headline biological claim
    was checked by nothing.

    TAIL CONVENTIONS ARE NOT ASSUMED, THEY ARE RECORDED. Measured 2026-09-05:
    the reported 0.016 reproduces as a ONE-SIDED Mann-Whitney (0.015609) while
    the reported Fisher 0.093 reproduces as a TWO-SIDED Fisher (0.093407). Both
    tails of both tests are therefore returned, and the CLAIMS entries below
    point at whichever one the document actually quotes. That keeps the checker
    honest about a mixed convention instead of silently blessing one.
    """
    cat = _hallmark_categories()
    if not cat:
        return None
    try:
        from scipy import stats as sps
    except ImportError:
        return None

    ax = [r for r in rows if r["adjustment"] == "axis_residualised"]
    if not ax:
        return None

    def short(s: str) -> str:
        return s.replace("SIG_HALLMARK_", "")

    out: dict[str, float] = {}
    base = {k for k, v in cat.items() if v == "immune"}
    for tag, immune in (("", base), ("_il6", base | {"IL6_JAK_STAT3_SIGNALING"})):
        a = [float(r["excess_z"]) for r in ax if short(r["signature"]) in immune]
        b = [float(r["excess_z"]) for r in ax if short(r["signature"]) not in immune]
        if not a or not b:
            return None
        out[f"immune_mwu_p_1s{tag}"] = float(
            sps.mannwhitneyu(a, b, alternative="less").pvalue)
        out[f"immune_mwu_p_2s{tag}"] = float(
            sps.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        ka = sum(r["beats_null"] == "True" for r in ax
                 if short(r["signature"]) in immune)
        kb = sum(r["beats_null"] == "True" for r in ax
                 if short(r["signature"]) not in immune)
        out[f"immune_fisher_p_2s{tag}"] = float(sps.fisher_exact(
            [[ka, len(a) - ka], [kb, len(b) - kb]], alternative="two-sided")[1])
        out[f"immune_beat{tag}"] = float(ka)
        out[f"immune_n{tag}"] = float(len(a))
        out[f"other_beat{tag}"] = float(kb)
        out[f"other_n{tag}"] = float(len(b))
    return out


def _test_count() -> int | None:
    """The number of tests pytest actually collects, right now.

    Read rather than remembered -- the count has drifted before. Returns None
    if pytest cannot be run, which is reported rather than silently skipped.
    """
    try:
        out = subprocess.run(
            [sys.executable, "-m", "pytest", "--collect-only", "-q"],
            cwd=ROOT, capture_output=True, text=True, timeout=600)
    except (OSError, subprocess.TimeoutExpired):
        return None
    m = re.search(r"(\d+)\s+tests?\s+collected", out.stdout)
    return int(m.group(1)) if m else None


def _site_counts() -> dict[str, float] | None:
    """Distinct tissue source sites per cohort, recomputed from the interim data.

    WHY THIS IS LIVE-MEASURED RATHER THAN READ FROM AN ARTEFACT. The manuscript
    and the poster both state "619 tissue source sites", and NO frozen artefact
    contains that number: `site_control.csv` holds only the sites that were
    *evaluable* in the negative control (98 pan-cancer, 15 NSCLC), which is a
    different quantity. Until now 619 was an unsourced count in a manuscript --
    exactly the failure mode this script exists to prevent.

    It is recomputed the same way `12_partition_variance.py:119` prints it --
    `cohort['tss'].nunique()` after the same patient filter -- but reading the
    expression index with `columns=[]` instead of the full matrix, which takes
    1.5 s rather than loading gigabytes. Confirmed 2026-09-04: 619 sites over
    7,168 patients pan-cancer, 68 over 944 NSCLC.

    Returns None if the interim data is absent, which is reported rather than
    silently skipped.
    """
    interim = ROOT / "data" / "interim"
    try:
        import numpy as np
        import pandas as pd

        sys.path.insert(0, str(ROOT / "src"))
        from aacr27 import data as _data

        out: dict[str, float] = {}
        nsclc = interim / "cohort_nsclc.parquet"
        if nsclc.exists():
            out["nsclc_n_sites"] = float(
                pd.read_parquet(nsclc)["tss"].nunique())
        meta_p, x_p, expr_p = (interim / "meta.parquet", interim / "X.npy",
                               interim / "expression_hugo.parquet")
        if meta_p.exists() and x_p.exists() and expr_p.exists():
            meta = pd.read_parquet(meta_p)
            X = np.asarray(np.load(x_p, mmap_mode="r"))
            meta, _ = _data.collapse_to_patient(meta, X)
            idx = pd.Index(pd.read_parquet(expr_p, columns=[]).index)
            idx = idx.astype(str).str[:12].unique()
            keep = meta["patient_id"].isin(idx)
            out["pancancer_n_sites"] = float(meta.loc[keep, "tss"].nunique())
        return out or None
    except Exception:
        return None


def _abstract_chars() -> tuple[int, int] | None:
    """(counted, limit) from 08_count_abstract.py, run now rather than recalled."""
    script = ROOT / "scripts" / "08_count_abstract.py"
    if not script.exists():
        return None
    try:
        out = subprocess.run([sys.executable, str(script)],
                             cwd=ROOT, capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return None
    m = re.search(r"TOTAL\s+([\d,]+)\s*/\s*([\d,]+)", out.stdout)
    if not m:
        return None
    return (int(m.group(1).replace(",", "")), int(m.group(2).replace(",", "")))


def _snapshot_size() -> tuple[int, int] | None:
    """(staged file count, staged bytes) from a LIVE rebuild, not from memory.

    Three documents quote this pair and nothing verified it, so it drifted
    every time anyone edited a staged file: 169 -> 170 -> 172 -> 173 files and
    1,447,396 -> 1,493,616 -> 1,756,296 -> 1,764,887 bytes inside two sessions,
    twice leaving a stale pair in the manuscript. It is derived by RUNNING
    `22_build_public_snapshot.py` into a temp directory and reading the numbers
    it prints, rather than by re-implementing its staging rules here -- a second
    definition of "what is staged" would be free to disagree with the first, and
    this project has already been bitten by exactly that (S3's two definitions
    of p). ~0.2 s, so it sits behind the `slow` gate with the other subprocesses.
    """
    script = ROOT / "scripts" / "22_build_public_snapshot.py"
    if not script.exists():
        return None
    with tempfile.TemporaryDirectory() as td:
        try:
            out = subprocess.run(
                [sys.executable, str(script), "--out", str(Path(td) / "snap")],
                cwd=ROOT, capture_output=True, text=True, timeout=300)
        except (OSError, subprocess.TimeoutExpired):
            return None
        if out.returncode != 0:
            # A leaking tree is a different failure, reported by that script.
            # Returning None here drops the two claims rather than checking
            # them against a number produced by a run that refused to finish.
            return None
        files = re.search(r"Audited ([\d,]+) staged file", out.stdout)
        size = re.search(r"Staged tree: ([\d,]+) bytes", out.stdout)
        if not (files and size):
            return None
        return (int(files.group(1).replace(",", "")),
                int(size.group(1).replace(",", "")))


def authority(*, slow: bool = True) -> dict[str, float]:
    """The authoritative values, read from the frozen artefacts on every run.

    NOTHING here is hard-coded from memory -- that is the whole point. If a
    frozen file is missing the run fails loudly rather than checking against a
    remembered constant.
    """
    table: dict[str, float] = {}
    missing: list[str] = []

    for cohort, prefix in (("nsclc_v3", "nsclc"), ("pancancer_v3", "pancancer")):
        d = RESULTS / cohort
        p = d / "summary.json"
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
            continue
        j = json.loads(p.read_text())

        pe = _dig(j, "primary_excess")
        if pe is None:
            missing.append(f"{p.relative_to(REPO)}:primary_excess")
            continue
        table[f"{prefix}_isi"] = float(pe["value"])
        table[f"{prefix}_ci_lo"] = float(pe["ci_lo"])
        table[f"{prefix}_ci_hi"] = float(pe["ci_hi"])

        # The secondary endpoint is quoted as a triple in the paper's Results
        # ("Delta-MAE = 0.0088 [0.0078, 0.0098]"). Before this key existed those
        # two triples were REPORTED AS UNMATCHED rather than checked.
        dm = _dig(j, "secondary_delta_mae")
        if dm is not None:
            table[f"{prefix}_dmae"] = float(dm["value"])
            table[f"{prefix}_dmae_lo"] = float(dm["ci_lo"])
            table[f"{prefix}_dmae_hi"] = float(dm["ci_hi"])
            # `n` on the secondary endpoint is the cohort size: 944 / 7168.
            if dm.get("n") is not None:
                table[f"{prefix}_n"] = float(dm["n"])

        # The number of cancer types is recorded only in the outcome-arm note,
        # e.g. "stratified by cancer_type (31 levels)". Parsed from the frozen
        # artefact, not asserted.
        notes = " ".join(str(x) for x in (_dig(j, "notes") or []))
        m = re.search(r"cancer_type\s*\((\d+)\s+levels?\)", notes)
        if m:
            table[f"{prefix}_n_types"] = float(m.group(1))

        for name, fname, fn in (
            # 16/16 signatures beat their null.
            ("sigs_beat", "immune_excess.csv",
             lambda r: sum(x["beats_null"] == "True" for x in r)),
            ("n_sigs", "immune_excess.csv", len),
            # 10/32 and 0/32 in the outcome arm -- 16 signatures x 2 adjustments.
            ("outcome_beat", "outcome_arm.csv",
             lambda r: sum(x["beats_null"] == "True" for x in r)),
            ("n_outcome_tests", "outcome_arm.csv", len),
            ("n_events", "outcome_arm.csv",
             lambda r: int(r[0]["n_events"])),
            # Headline site AUROC is the MEDIAN over evaluable sites; the paper
            # says so explicitly ("recoverable ... at median AUROC 0.998").
            ("site_auroc", "site_control.csv",
             lambda r: _median([float(x["auroc"]) for x in r])),
            ("n_sites_evaluable", "site_control.csv", len),
        ):
            f = d / fname
            if not f.exists():
                missing.append(str(f.relative_to(REPO)))
                continue
            table[f"{prefix}_{name}"] = float(fn(_rows(f)))

        oa = d / "outcome_arm.csv"
        if oa.exists():
            rows = _rows(oa)
            derived = _mde_c_index(rows)
            table[f"{prefix}_mde"] = derived
            # Cross-check the derivation wherever the stored column exists.
            # Pan-cancer carries it; NSCLC predates it. If the derivation stops
            # reproducing the stored value, every MDE claim below is untrustworthy
            # and we must NOT quietly keep checking against it.
            stored = [float(r["mde_80pct_c_index"]) for r in rows
                      if r.get("mde_80pct_c_index") not in (None, "")
                      and r["adjustment"] == "axis_residualised"]
            if stored and abs(_median(stored) - derived) > 1e-9:
                raise SystemExit(
                    f"MDE DERIVATION BROKEN for {cohort}: derived "
                    f"{derived!r} != stored {_median(stored)!r}. Refusing to "
                    "check MDE claims against a derivation that no longer "
                    "reproduces the artefact.")

            # The REAL effect on the C-index scale, for the signatures that
            # actually beat their null. 07-ABSTRACT-DRAFT.md and
            # 00-PROJECT-BRIEF.md both compared "the real pan-cancer effect"
            # against the 0.046 MDE, but said "~0.03" -- a number that matched
            # NO authority key (mde 0.0167, median_delta_r 0.0199,
            # incr_r2 0.0362). --coverage flagged it as unchecked and the
            # source turned out to be this: observed - null_mean over the
            # beats_null rows. The prose now quotes the measured median
            # instead of a rounded recollection.
            beat = [float(r["observed"]) - float(r["null_mean"]) for r in rows
                    if str(r.get("beats_null", "")).strip().lower()
                    in ("true", "1")]
            if beat:
                table[f"{prefix}_outcome_delta_c_median"] = _median(beat)
                table[f"{prefix}_outcome_delta_c_max"] = max(beat)

    # The A9/A10 stable-sort sensitivity run. Read as an AUTHORITY, not as a
    # frozen result: `nsclc_v3/` remains the pre-registered primary and this
    # sits beside it. The shift is DERIVED from the two runs rather than stored,
    # so it cannot drift away from the pair it describes -- the failure mode
    # that let a stale "+0.180" survive on the poster.
    p = RESULTS / "nsclc_v3_stablesort" / "summary.json"
    if p.exists():
        j = json.loads(p.read_text())["primary_excess"]
        table["sortaudit_isi"] = float(j["value"])
        table["sortaudit_ci_lo"] = float(j["ci_lo"])
        table["sortaudit_ci_hi"] = float(j["ci_hi"])
        frozen = RESULTS / "nsclc_v3" / "summary.json"
        if frozen.exists():
            f = json.loads(frozen.read_text())["primary_excess"]
            table["sortaudit_shift"] = float(f["value"]) - float(j["value"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # S3's multiplicity columns for the FROZEN NSCLC cohort. DERIVED here, from
    # the frozen run's own excess and interval, rather than read from a stored
    # column -- the stored file predates the columns and is never modified.
    # Deriving keeps the authority tied to the frozen artefact: if excess_z or
    # its interval ever changed, these move with it and the documents go red.
    ie = RESULTS / "nsclc_v3" / "immune_excess.csv"
    if ie.exists():
        import numpy as _np
        import pandas as _pd

        sys.path.insert(0, str(ROOT / "src"))
        from aacr27 import stats as _st
        _d = _pd.read_csv(ie)
        _sp = __import__("scipy.stats", fromlist=["stats"])
        _se = (_d["excess_hi"] - _d["excess_lo"]) / (2 * _sp.norm.ppf(0.975))
        with _np.errstate(invalid="ignore", divide="ignore"):
            _pv = 2 * _sp.norm.sf(_np.abs(_d["excess_z"] / _se))
        _ok = _np.isfinite(_pv)
        _q = _np.full(len(_pv), _np.nan)
        _rej, _qv = _st.bh_fdr(_pv[_ok])
        _q[_ok] = _qv
        table["nsclc_bh_min_p"] = float(_np.nanmin(_pv))
        table["nsclc_bh_max_q"] = float(_np.nanmax(_q))
        # Asserted, not stored: an authority key nothing reads is the
        # models.py:341 defect. The 16/16 count itself is checked by the
        # FRACTIONS mechanism wherever a document writes it.
        assert int(_np.nansum(_q < 0.05)) == len(_d), (
            "frozen NSCLC no longer gives 16/16 after BH")
    else:
        missing.append(str(ie.relative_to(REPO)))

    # A5, both cohorts. The NSCLC keys are UNPREFIXED for history: they were the
    # only partition-variance numbers in the project for five handoffs and the
    # documents' patterns anchor on them. The pan-cancer run landed 2026-09-05
    # and takes the `pancancer_` prefix the rest of this table uses; NSCLC also
    # gets prefixed aliases so a new pattern never has to guess which cohort an
    # unprefixed key means.
    for cohort, prefix in (("nsclc", "nsclc"), ("pancancer", "pancancer")):
        p = RESULTS / f"partition_variance_{cohort}" / "summary.json"
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
            continue
        j = json.loads(p.read_text())
        vals = {
            "partition_sd": float(j["isi_sd_across_partitions"]),
            "bootstrap_se": float(j["bootstrap_se"]),
            "combined_se": float(j["combined_se"]),
            "isi_mean": float(j["isi_mean"]),
            "honest_ci_lo": float(j["honest_ci"][0]),
            "honest_ci_hi": float(j["honest_ci"][1]),
            "reported_ci_lo": float(j["reported_ci"][0]),
            "reported_ci_hi": float(j["reported_ci"][1]),
            "interval_understated": float(j["interval_understated_by"]),
            # Documents quote this as a PERCENT ("4.5% too narrow"), the JSON
            # stores a fraction. Without this key the percent pattern silently
            # matched nothing -- a gap the coverage count would not have shown.
            "interval_understated_pct":
                100.0 * float(j["interval_understated_by"]),
            "n_partitions_used": float(j["n_partitions_used"]),
            "n_partitions_refused": float(j["n_partitions_refused"]),
        }
        # Share of total variance attributable to partition choice. Derived
        # here rather than read, because `12_partition_variance.py` prints it
        # but does not persist it -- and the paper quotes it for both cohorts.
        sd, se = vals["partition_sd"], vals["bootstrap_se"]
        vals["partition_var_pct"] = 100.0 * sd**2 / (sd**2 + se**2)
        for k, v in vals.items():
            table[f"{prefix}_{k}"] = v
            if prefix == "nsclc":
                table[k] = v  # legacy unprefixed aliases, NSCLC only

    # A7, the scorer-sensitivity arms. macOS only: `scorer_sensitivity_hpc4/`
    # is diagnostic and deliberately NOT an authority (its mean-z arms return
    # 0.2964, which is A9, not a scorer effect).
    #
    # `ssgsea__disattenuated` is `null` in the JSON by design -- the registered
    # estimand is UNDEFINED under ssGSEA -- so every read is guarded rather than
    # assumed present. A None slipping into the table would compare as a crash,
    # not as a mismatch.
    p = RESULTS / "scorer_sensitivity_local" / "summary.json"
    if p.exists():
        j = json.loads(p.read_text())
        for arm, short in (("mean_z__disattenuated", "a7_meanz_dis"),
                           ("mean_z__uncorrected", "a7_meanz_unc"),
                           ("ssgsea__uncorrected", "a7_ssgsea_unc")):
            a = j["arms"][arm]
            if a.get("isi") is None:
                continue
            table[short] = float(a["isi"])
            table[short + "_lo"] = float(a["lo"])
            table[short + "_hi"] = float(a["hi"])
        for scorer, short in (("mean_z", "a7_recon_meanz"),
                              ("ssgsea", "a7_recon_ssgsea")):
            r = j["reconstructed"][scorer]
            # The ssGSEA reconstruction is NEGATIVE, and the two rules that
            # read it need opposite conventions, so both are stored rather
            # than one being made to serve both:
            #   * LONG_DECIMAL captures UNSIGNED digits (its regex starts at a
            #     digit), so the magnitude must be present or a correct
            #     full-precision quote is reported as unknown.
            #   * TRIPLES compares a SIGNED point estimate, so the signed value
            #     must be present or a sign flip in the manuscript passes.
            # Storing only the magnitude would have silently accepted a
            # reconstruction whose sign had been reversed -- which, for this
            # particular result, is the whole finding.
            table[short] = abs(float(r["isi"]))
            table[short + "_lo"] = abs(float(r["lo"]))
            table[short + "_hi"] = abs(float(r["hi"]))
            table[short + "_signed"] = float(r["isi"])
            table[short + "_signed_lo"] = float(r["lo"])
            table[short + "_signed_hi"] = float(r["hi"])
        rel = j["split_half_reliability"]
        table["a7_meanz_rel_obs"] = float(rel["mean_z"]["observed_mean"])
        table["a7_meanz_rel_null"] = float(rel["mean_z"]["null_mean"])
        table["a7_ssgsea_rel_obs"] = float(rel["ssgsea"]["observed_mean"])
        table["a7_ssgsea_rel_null"] = float(rel["ssgsea"]["null_mean"])
        for row in j["agreement"]:
            if row["label"].startswith("uncorrected"):
                table["a7_agree_unc_pearson"] = float(row["pearson"])
                table["a7_agree_unc_spearman"] = float(row["spearman"])
            elif row["label"].startswith("reconstructed"):
                table["a7_agree_recon_pearson"] = float(row["pearson"])
            elif row["label"].startswith("residualised"):
                table["a7_agree_resid_pearson"] = float(row["pearson"])
                table["a7_agree_resid_absdiff"] = float(row["mean_abs_diff"])
        table["a7_beats_ssgsea_unc"] = float(j["beats_null"]["ssgsea__uncorrected"])
        table["a7_n_sigs"] = float(j["n_signatures"])
        if "a7_ssgsea_unc" in table and "a7_meanz_unc" in table:
            # "the effect survives at roughly N% of its magnitude" -- DERIVED
            # from the pair so it cannot drift away from the two arms it
            # compares, which is the failure the sortaudit block also avoids.
            table["a7_frac_of_meanz_pct"] = (
                100.0 * table["a7_ssgsea_unc"] / table["a7_meanz_unc"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # A10's ssGSEA re-run under the PINNED (stable) tie ordering. This sits
    # beside `scorer_sensitivity_local/` exactly as `nsclc_v3_stablesort/` sits
    # beside `nsclc_v3/`: the pre-fix arm is what was run and reported, this is
    # what the fixed code produces, and both are quoted.
    #
    # Only the two ssGSEA arms were requested, so `arms_complete` is False and
    # the mean-z keys are absent by design -- that scorer never touches the rank
    # table, so re-running it would measure nothing. `ssgsea__disattenuated` is
    # `null` here for the SAME reason it is null pre-fix (the reliability
    # estimator collapses, not the sort), so its absence is not evidence about
    # the fix and is guarded rather than read.
    #
    # The shift is DERIVED from the pair rather than stored, so it cannot drift
    # away from the two numbers it describes.
    p = RESULTS / "scorer_sensitivity_ssgsea_stable" / "summary.json"
    if p.exists():
        j = json.loads(p.read_text())
        a = j["arms"]["ssgsea__uncorrected"]
        if a.get("isi") is not None:
            table["a10ss_ssgsea_unc"] = float(a["isi"])
            table["a10ss_ssgsea_unc_lo"] = float(a["lo"])
            table["a10ss_ssgsea_unc_hi"] = float(a["hi"])
            if "a7_ssgsea_unc" in table:
                table["a10ss_shift"] = table["a7_ssgsea_unc"] - float(a["isi"])
                # The shift expressed against the pre-fix interval's own half
                # width, as a PERCENT. A raw shift of 0.006 means nothing until
                # it is set against the precision of the thing it moves.
                half = (table["a7_ssgsea_unc_hi"]
                        - table["a7_ssgsea_unc_lo"]) / 2.0
                table["a10ss_shift_pct_of_halfwidth"] = (
                    100.0 * table["a10ss_shift"] / half)
                if "a7_meanz_unc" in table:
                    table["a10ss_frac_of_meanz_pct"] = (
                        100.0 * float(a["isi"]) / table["a7_meanz_unc"])
        table["a10ss_beats_ssgsea_unc"] = float(
            j["beats_null"]["ssgsea__uncorrected"])
        r = j["reconstructed"]["ssgsea"]
        table["a10ss_recon_signed"] = float(r["isi"])
        table["a10ss_rel_obs"] = float(
            j["split_half_reliability"]["ssgsea"]["observed_mean"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # A9's macOS bisection half. Not a headline -- it is one side of a
    # two-machine experiment and answers nothing alone -- but its pooled value
    # is quoted in 14-SCIENCE-AUDIT.md and results/README.md to full precision,
    # and the full-precision rule correctly refused to accept a 12-digit number
    # backed by no artefact. Reading it here makes those quotes VERIFIED rather
    # than exempted, which is the whole point of that rule.
    p = RESULTS / "bisect_macos.npz"
    if p.exists():
        import numpy as _np
        with _np.load(p, allow_pickle=True) as _z:
            table["bisect_macos_pooled"] = float(_z["pooled_z_excess"].ravel()[0])
    else:
        missing.append(str(p.relative_to(REPO)))

    # The HPC4 half, retrieved from cluster job 121439 on 2026-09-06. Same
    # reasoning as the macOS half above: 14-SCIENCE-AUDIT.md quotes the pooled
    # value and the platform delta to full precision, and the full-precision
    # rule refused all three the moment they were pasted in -- which is the rule
    # working. The delta is DERIVED here rather than typed, so it cannot drift
    # away from the two values it is the difference of.
    p = RESULTS / "bisect_hpc4.npz"
    if p.exists():
        import numpy as _np
        with _np.load(p, allow_pickle=True) as _z:
            table["bisect_hpc4_pooled"] = float(_z["pooled_z_excess"].ravel()[0])
        if "bisect_macos_pooled" in table:
            table["bisect_platform_delta"] = (
                table["bisect_macos_pooled"] - table["bisect_hpc4_pooled"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # The A2/A3 gap closures live in their own artefact, not in summary.json.
    p = RESULTS / "pancancer_v3" / "science_gaps.json"
    if p.exists():
        j = json.loads(p.read_text())
        rn = j.get("rotation_null", {})
        table["rotation_p"] = float(rn["p"])
        table["rotation_observed"] = float(rn["observed"])
        table["rotation_null_mean"] = float(rn["null_family_mean"])
        table["rotation_null_lo"] = float(rn["ci_lo"])
        table["rotation_null_hi"] = float(rn["ci_hi"])
        table["rotation_B"] = float(rn["B"])
        table["m_eff"] = float(j["effective_tests"])
        table["control_c_calibrated_r2"] = float(j["control_c_calibrated_r2"])
        table["label_site_variance_r2"] = float(j["label_site_variance_median_r2"])
        # The NSCLC counterpart, added 2026-09-06. It exists at all only because
        # the A11 cohort-scope fix was repaired: before that, this file's twin
        # in `nsclc_v3_stablesort/` held THIS directory's pan-cancer values, and
        # the two agreeing to seventeen significant figures is how A11 was
        # found. Holding both as authorities means a future recurrence shows up
        # as two keys with one value rather than as a coincidence someone has to
        # notice by eye.
        q = RESULTS / "nsclc_v3_stablesort" / "science_gaps.json"
        if q.exists():
            jn = json.loads(q.read_text())
            if jn.get("label_site_variance_median_r2") is not None:
                table["nsclc_label_site_variance_r2"] = float(
                    jn["label_site_variance_median_r2"])
                table["nsclc_control_c_calibrated_r2"] = float(
                    jn["control_c_calibrated_r2"])
                table["nsclc_control_c_perm_median"] = float(
                    jn["control_c_perm_median"])
                table["nsclc_gaps_n_patients_used"] = float(
                    jn["cohort_n_patients_used"])
                table["nsclc_gaps_n_patients_in_file"] = float(
                    jn["cohort_n_patients_in_file"])
        # The manuscript quotes both of these as PERCENTAGES ("a median 44% of
        # the label", "gives 4.1%"), and the fraction-valued keys above could
        # never match a percent literal. `label_site_variance_r2` in particular
        # was read into this table and then used by NO pattern at all -- an
        # authority value with no claim attached, which looks like coverage and
        # is not. That is precisely what `--coverage` now exists to surface.
        table["control_c_calibrated_pct"] = 100.0 * float(j["control_c_calibrated_r2"])
        table["label_site_variance_pct"] = (
            100.0 * float(j["label_site_variance_median_r2"]))
        table["small_panel_spearman"] = float(j["small_panel_spearman"][0])
        # The small-panel bracket table. Quoted in prose as "0.084 at k=160,
        # 0.141 at 80, 0.214 at 40, 0.290 at 20 and 0.453 at k=10" and, until
        # now, entirely unchecked -- five numbers carrying the paper's claim
        # that the correction matters most for clinical-sized panels.
        for row in j["small_panel_bracket"]:
            table[f"small_panel_k{int(row['k'])}"] = float(row["median_gap"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # Median image-signature correlation and median split-robustness Delta r.
    # These are quoted in the Results opening ("median ... 0.630 pan-cancer and
    # 0.400 in NSCLC", "median Delta r = +0.020 and +0.023") and had no source
    # in this table. They are medians over per_signature.csv, recomputed here
    # rather than remembered -- the same rule as everything else in authority().
    for cohort, prefix in (("nsclc_v3", "nsclc"), ("pancancer_v3", "pancancer")):
        p = RESULTS / cohort / "per_signature.csv"
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
            continue
        rows = _rows(p)
        table[f"{prefix}_median_r"] = _median(
            float(r["r_preserved_site"]) for r in rows)
        table[f"{prefix}_median_delta_r"] = _median(
            float(r["delta_r"]) for r in rows)

    # ---- added 2026-09-05: the reliability, decomposition and label-side -----
    # numbers. These are the three families handoff #24 listed as still
    # unchecked, and between them they carry the manuscript's entire
    # "reliability, not composition, is the correction that matters" section
    # and its negative-control section.
    #
    # AGGREGATION IS MEASURED, NOT ASSUMED, and it is not uniform. Verified
    # 2026-09-05 by computing both: the alpha LEVELS the paper quotes are
    # MEDIANS over the 16 signatures (median alpha_null_mean_raw = 0.974044 ->
    # "median alpha 0.974"), while the reliability GAPS it quotes are MEANS
    # (mean reliability_gap = 0.184291 -> "+0.184"; the median is 0.152609 and
    # would have been reported as drift). Guessing one convention for both
    # would have produced four false mismatches on correct prose.
    for cohort, prefix in (("nsclc_v3", "nsclc"), ("pancancer_v3", "pancancer")):
        p = RESULTS / cohort / "immune_excess.csv"
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
            continue
        rows = _rows(p)
        table[f"{prefix}_alpha_obs_raw"] = _median(
            float(r["alpha_observed_raw"]) for r in rows)
        table[f"{prefix}_alpha_obs_resid"] = _median(
            float(r["alpha_observed"]) for r in rows)
        table[f"{prefix}_alpha_null_raw"] = _median(
            float(r["alpha_null_mean_raw"]) for r in rows)
        table[f"{prefix}_alpha_null_resid"] = _median(
            float(r["alpha_null_mean"]) for r in rows)
        table[f"{prefix}_null_mean_r"] = _median(
            float(r["null_mean_r"]) for r in rows)
        # The per-signature excess RANGE, quoted twice in the manuscript
        # ("per-signature excess 0.165-0.390 and 0.148-0.440") and by nothing
        # else. It is the spread behind "all 16 signatures exceeded their null".
        ez = [float(r["excess_z"]) for r in rows]
        table[f"{prefix}_excess_min"] = min(ez)
        table[f"{prefix}_excess_max"] = max(ez)
        gap_resid = [float(r["reliability_gap"]) for r in rows]
        gap_raw = [float(r["alpha_observed_raw"]) - float(r["alpha_null_mean_raw"])
                   for r in rows]
        table[f"{prefix}_rel_gap_resid"] = sum(gap_resid) / len(gap_resid)
        table[f"{prefix}_rel_gap_raw"] = sum(gap_raw) / len(gap_raw)
        # The "13-fold and 5.3-fold" claim is the RATIO OF THE MEANS, which is
        # what the prose's own two numbers divide to. The ratio of the medians
        # is 11.4, so this is a real distinction and not a rounding choice.
        table[f"{prefix}_rel_gap_fold"] = (table[f"{prefix}_rel_gap_resid"]
                                           / table[f"{prefix}_rel_gap_raw"])
        for r in rows:
            if "ANGIOGENESIS" in r["signature"]:
                table[f"{prefix}_rel_gap_angio"] = float(r["reliability_gap"])

        p = RESULTS / cohort / "three_way.csv"
        if p.exists():
            tw = _rows(p)
            table[f"{prefix}_r_resid_refit"] = _median(
                float(r["r_residual_refit"]) for r in tw)
            table[f"{prefix}_r_axis"] = _median(float(r["r_axis"]) for r in tw)

        p = RESULTS / cohort / "decomposition.csv"
        if p.exists():
            de = _rows(p)
            table[f"{prefix}_r_partial"] = _median(
                float(r["r_partial"]) for r in de)
            table[f"{prefix}_r_unadjusted"] = _median(
                float(r["r_image"]) for r in de)
            table[f"{prefix}_incr_r2"] = _median(
                float(r["incremental_image_r2"]) for r in de)

        p = RESULTS / cohort / "per_signature.csv"
        if p.exists():
            ps = _rows(p)
            table[f"{prefix}_r_covariates"] = _median(
                float(r["r_covariates"]) for r in ps)
            table[f"{prefix}_emb_over_cov"] = _median(
                float(r["embedding_over_covariates"]) for r in ps)

        oa = RESULTS / cohort / "outcome_arm.csv"
        if oa.exists() and prefix == "pancancer":
            got = _immune_contrast(_rows(oa))
            if got:
                table.update(got)

    # The label-side site-variance control. `label_site_variance_median_r2` is
    # already read from science_gaps.json above; these two are the numbers the
    # same sentence quotes and nothing sourced them.
    p = RESULTS / "pancancer_v3" / "label_site_variance.csv"
    if p.exists():
        lsv = _rows(p)
        table["label_site_given_type_pct"] = 100.0 * _median(
            float(r["r2_site_given_type"]) for r in lsv)
        table["label_plate_within_site"] = _median(
            float(r["r2_plate_within_site"]) for r in lsv)
    else:
        missing.append(str(p.relative_to(REPO)))

    # The ancestry arm. Six numbers in 09-PAPER-DRAFT.md, none of them checked
    # by any earlier version of this file, and the arm is reported as a
    # supplementary table so nothing else would catch a drifted count.
    p = RESULTS / "ancestry" / "summary.json"
    if p.exists():
        j = json.loads(p.read_text())
        table["ancestry_n"] = float(j["n_patients"])
        for grp, val in j["counts_patient_level"].items():
            table[f"ancestry_{grp.lower()}"] = float(val)
        table["ancestry_cramers_v"] = float(j["site_cramers_v"]["v"])
        table["ancestry_perm_p95"] = float(j["site_cramers_v"]["v_perm_p95"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # A9's two-platform record. ONE FILE PER MACHINE, each written by an actual
    # run of 24_split_determinism.py --write on that machine, so the
    # cross-platform table in 14-SCIENCE-AUDIT.md and 09-PAPER-DRAFT.md is
    # quoted from artefacts rather than from a transcript. Without these keys
    # the full-precision rule correctly refuses the four synthetic-ISI values --
    # which is how this gap was found.
    for tag, fname in (("macos", "split_determinism_darwin_arm64.json"),
                       ("hpc4", "split_determinism_linux_x86_64.json")):
        p = RESULTS / fname
        if p.exists():
            j = json.loads(p.read_text())
            table[f"splitdet_{tag}_isi_shipped"] = float(j["synthetic_isi_shipped"])
            table[f"splitdet_{tag}_isi_stable"] = float(j["synthetic_isi_stable"])
            table[f"splitdet_{tag}_moved"] = float(j["n_patients_fold_changed"])
            table[f"splitdet_{tag}_sites"] = float(j["n_sites"])
            table[f"splitdet_{tag}_sites_tied"] = float(j["n_sites_tied"])
        else:
            missing.append(str(p.relative_to(REPO)))

    # A10's two-platform record, on exactly the pattern above: one file per
    # machine, each written by a real `26_sort_audit.py --write` on that machine.
    # Both hashes are recorded deliberately. `hash_depth_stable` is claimed to
    # AGREE across platforms and `hash_depth_quicksort` to DIFFER, so keeping
    # both makes the A10 table falsifiable in either direction rather than only
    # confirmable -- a check that can only ever agree is the vacuous kind this
    # project has shipped before.
    for tag, fname in (("macos", "sort_audit_darwin_arm64.json"),
                       ("hpc4", "sort_audit_linux_x86_64.json")):
        p = RESULTS / fname
        if p.exists():
            j = json.loads(p.read_text())
            table[f"sortaudit_{tag}_samples"] = float(j["n_samples"])
            table[f"sortaudit_{tag}_genes"] = float(j["n_genes"])
            table[f"sortaudit_{tag}_tied_rows"] = float(j["tied_rows"])
            table[f"sortaudit_{tag}_tied_elements"] = float(j["tied_elements"])
            table[f"sortaudit_{tag}_cells_diff"] = float(j["ssgsea_cells_differing"])
            table[f"sortaudit_{tag}_cells_total"] = float(j["ssgsea_cells_total"])
            table[f"sortaudit_{tag}_sigs_affected"] = float(j["ssgsea_signatures_affected"])
            table[f"sortaudit_{tag}_max_diff"] = float(j["ssgsea_max_abs_diff"])
            table[f"sortaudit_{tag}_score_sd"] = float(j["ssgsea_mean_score_sd"])
            table[f"sortaudit_{tag}_frac_sd"] = float(j["ssgsea_max_diff_as_frac_of_sd"])
            table[f"sortaudit_{tag}_tied_sites"] = float(j["splits_tied_site_sizes"])
            table[f"sortaudit_{tag}_hash_stable"] = j["hash_depth_stable"]
            table[f"sortaudit_{tag}_hash_quicksort"] = j["hash_depth_quicksort"]
        else:
            missing.append(str(p.relative_to(REPO)))

    # The A10 cross-platform claim, asserted directly against the two artefacts
    # rather than left to prose. The claim has TWO halves and both are checked,
    # because only checking the agreement half would pass just as happily if the
    # two files were copies of each other.
    if "sortaudit_macos_hash_stable" in table and "sortaudit_hpc4_hash_stable" in table:
        if table["sortaudit_macos_hash_stable"] != table["sortaudit_hpc4_hash_stable"]:
            raise SystemExit(
                "A10 BROKEN: the STABLE ssGSEA depth hash was supposed to be "
                "platform-invariant but the two artefacts disagree:\n"
                f"  macOS {table['sortaudit_macos_hash_stable']}\n"
                f"  HPC4  {table['sortaudit_hpc4_hash_stable']}")
        if table["sortaudit_macos_hash_quicksort"] == table["sortaudit_hpc4_hash_quicksort"]:
            raise SystemExit(
                "A10 UNSUPPORTED: the QUICKSORT depth hashes agree across "
                "platforms, so the two artefacts do not demonstrate a "
                "platform-dependent tie order and the A10 write-up overstates "
                f"what was measured (both {table['sortaudit_macos_hash_quicksort']}). "
                "Check that the two files really came from different machines.")
    # The hashes are strings; drop them before the numeric comparison stage.
    for k in [k for k in table if k.startswith("sortaudit_") and "_hash_" in k]:
        del table[k]

    if missing:
        raise SystemExit(_missing_artefact_message(missing))

    # The number of entries in the manuscript's reference list, COUNTED from
    # the list itself. Unlike everything else here the authority is a document
    # rather than an artefact, and that is the point: the failure this catches
    # is the prose disagreeing with the list beneath it. The list dropped from
    # thirteen entries to twelve on 2026-09-04 and the prose still said "all
    # thirteen re-verified" -- invisible to every pattern in this file, because
    # they all match digits and the claim was a word.
    p = REPO / "09-PAPER-DRAFT.md"
    if p.exists():
        body = p.read_text(errors="replace")
        _, _, refs = body.partition("\n## References\n")
        if refs:
            nums = [int(m.group(1)) for m in
                    re.finditer(r"^(\d{1,2})\.\s+[A-Z]", refs, re.MULTILINE)]
            if nums:
                table["n_references"] = float(max(nums))
                # A gap or a duplicate in the numbering is its own defect --
                # citation-sequence order means 1..N with no holes.
                if sorted(nums) != list(range(1, max(nums) + 1)):
                    raise SystemExit(
                        "REFERENCE NUMBERING IS NOT 1..N: found "
                        f"{sorted(nums)}. Citation-sequence order requires a "
                        "contiguous run; refusing to check counts against it.")

    if slow:
        n = _test_count()
        if n is not None:
            table["test_count"] = float(n)
        ab = _abstract_chars()
        if ab is not None:
            table["abstract_chars"], table["abstract_limit"] = float(ab[0]), float(ab[1])
            # DERIVED, not stored: the two counts were checked and the
            # subtraction between them was not, so 07-ABSTRACT-DRAFT.md read
            # "2,515 / 2,600 with 83 characters of headroom" -- 2,515 correct,
            # 2,600 correct, and the difference between them wrong by two, on
            # the line that says how much room is left before a hard deadline.
            # Found 2026-09-07 by running the counter, not by reading.
            table["abstract_headroom"] = float(ab[1] - ab[0])
        sc = _site_counts()
        if sc is not None:
            table.update(sc)
        snap = _snapshot_size()
        if snap is not None:
            table["snapshot_files"], table["snapshot_bytes"] = (
                float(snap[0]), float(snap[1]))
    return table


def every_stored_run() -> dict[str, float]:
    """Every value any stored run produced, frozen or superseded.

    Used ONLY to whitelist full-precision numbers quoted as provenance --
    `results/README.md` says what each superseded run produced, and those are
    real measurements, not drift. A number matching no stored run at all is the
    thing worth failing on.
    """
    known: dict[str, float] = {}
    for p in sorted(RESULTS.glob("*/summary.json")):
        run = p.parent.name
        try:
            j = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        for block in ("primary_excess", "secondary_delta_mae"):
            pe = _dig(j, block)
            if isinstance(pe, dict):
                for field in ("value", "ci_lo", "ci_hi"):
                    v = pe.get(field)
                    if isinstance(v, (int, float)):
                        known[f"{run}.{block}.{field}"] = float(v)
        for field in ("isi_mean", "isi_sd_across_partitions", "bootstrap_se",
                      "combined_se", "interval_understated_by"):
            v = j.get(field)
            if isinstance(v, (int, float)):
                known[f"{run}.{field}"] = float(v)
        for field in ("reported_ci", "honest_ci"):
            v = j.get(field)
            if isinstance(v, list) and len(v) == 2:
                known[f"{run}.{field}[0]"] = float(v[0])
                known[f"{run}.{field}[1]"] = float(v[1])
        # Scorer-sensitivity runs store their results under `arms` and
        # `reconstructed` rather than `primary_excess`, so none of the loops
        # above saw them. That made every HPC4 arm value an unknown
        # full-precision number the moment 14-SCIENCE-AUDIT.md tabulated the
        # macOS/HPC4 contrast. These are real measurements from a stored
        # DIAGNOSTIC run -- provenance, not headline claims -- which is exactly
        # what this whitelist is for. Reliabilities are included because the
        # same contrast quotes them.
        for block in ("arms", "reconstructed"):
            b = j.get(block)
            if not isinstance(b, dict):
                continue
            for name, vals in b.items():
                if not isinstance(vals, dict):
                    continue
                for field in ("isi", "lo", "hi"):
                    v = vals.get(field)
                    if isinstance(v, (int, float)):
                        known[f"{run}.{block}.{name}.{field}"] = float(v)
                        known[f"{run}.{block}.{name}.{field}|abs"] = abs(float(v))
        shr = j.get("split_half_reliability")
        if isinstance(shr, dict):
            for scorer, vals in shr.items():
                if not isinstance(vals, dict):
                    continue
                for field in ("observed_mean", "null_mean"):
                    v = vals.get(field)
                    if isinstance(v, (int, float)):
                        known[f"{run}.reliability.{scorer}.{field}"] = float(v)
    p = RESULTS / "pancancer_v3" / "science_gaps.json"
    if p.exists():
        j = json.loads(p.read_text())
        for k, v in (j.get("rotation_null") or {}).items():
            if isinstance(v, (int, float)):
                known[f"science_gaps.rotation_null.{k}"] = float(v)
        for k, v in j.items():
            if isinstance(v, (int, float)):
                known[f"science_gaps.{k}"] = float(v)
    return known


# THE CHECK MUST BE ABLE TO FAIL.
#
# The obvious design -- a regex matching the CORRECT literal, then comparing it
# to the authoritative value -- is vacuous: such a regex can only ever match the
# right answer, so it passes by construction. The first version of this script
# did exactly that and reported "62 claims OK" while checking nothing. Every
# pattern below therefore anchors on CONTEXT and CAPTURES whatever number the
# document actually wrote, so a wrong number matches the pattern and then fails
# the comparison. `--self-test` proves that is still true.
#
# (key, regex with one capture group, tolerance, description, line-guard or None)
#
# The line-guard is a second regex the LINE must match before the pattern is
# applied. It exists so a pattern can be loose enough to catch a wholesale
# replacement ("p = 0.5") without firing on unrelated numbers elsewhere.
#
# A guard written as "PARA:<regex>" is matched against the whole PARAGRAPH the
# line sits in (blank-line delimited) instead of the line alone. This exists
# because line-scoped guards have silently died four times in this file's
# history, every time for the same reason: the prose was rewrapped and the
# guard word ended up on a different line from the number, so the guard
# excluded the only line carrying the value and the pattern checked nothing.
# Session 25 found one such guard completely dead; two more were caught by
# --self-test; #27 shipped a fourth. A paragraph is the honest scope for a
# guard like /NSCLC/ -- "this paragraph is about NSCLC" is a claim about the
# text that rewrapping cannot invalidate, whereas "this LINE mentions NSCLC"
# is an accident of where the wrap fell. Line scope remains the default and is
# unchanged, so no existing guard changes behaviour; paragraph scope is opt-in
# per pattern. A paragraph guard is still a guard: --self-test proves it can
# fail, and it is NOT a licence to write a guard so broad it never excludes.
#
# Tolerance is absolute and set by the digits the prose quotes: a document
# writing "0.3182" cannot be held to 1e-15.
CLAIMS: list[tuple[str, str, float, str, str | None]] = [
    # Labelled scalars: the label is the anchor, the number is captured.
    # The A9/A10 stable-sort sensitivity run. These are quoted only in
    # Limitations 8 and in `results/README.md`, both of which name the cohort in
    # the surrounding paragraph, so the guards are paragraph-scoped on the
    # phrase that is unique to that discussion rather than on "NSCLC" -- the
    # limitation sits in a section that names both cohorts repeatedly.
    ("sortaudit_isi", r"ISI of\s+\**(\d\.\d+)", 5e-5,
     "NSCLC ISI under the corrected sorts", "PARA:corrected ordering"),
    ("sortaudit_shift", r"moves by\s+\**(\d\.\d+)", 5e-5,
     "how far the corrected sorts move the NSCLC point estimate",
     "PARA:corrected ordering"),
    # S3's DERIVED NSCLC multiplicity values. These are the project's first
    # claims written in scientific notation, which is why `_tol_for` had to
    # learn to scale the half-width by the exponent -- under the old decimal
    # rule a claim of 1.669e-41 carried a 5e-8 tolerance and could not fail.
    ("nsclc_bh_min_p", r"min p = (\d\.\d+e-\d+)", 0.0,
     "smallest per-signature p in the frozen NSCLC cohort",
     "PARA:schema asymmetry is RESOLVED"),
    ("nsclc_bh_max_q", r"max q = (\d\.\d+e-\d+)", 0.0,
     "largest BH q in the frozen NSCLC cohort",
     "PARA:schema asymmetry is RESOLVED"),
    # A5, partition variance. Both cohorts now report the same five quantities,
    # so every pattern here is PARAGRAPH-guarded on its cohort. Line guards were
    # tried first and could not work: the manuscript wraps at ~80 columns and
    # "NSCLC" repeatedly landed on a different line from the number it governs.
    # Rewrapping the prose to force them together made the sentences read like
    # a checker's output rather than a manuscript ("an NSCLC within-partition
    # bootstrap standard error"), which is the wrong trade -- the tool exists to
    # serve the paper.
    #
    # NOTE the capture is `(\d+\.\d+)`, not `(\d\.\d+)`. The single-digit form
    # was live here until 2026-09-05 and could not capture a two-digit percent
    # at all: on "is 16.0% too narrow" the `\s+` before the capture cannot be
    # satisfied by the "1", so the pattern simply did not match and the number
    # went unchecked in silence. It passed for five handoffs only because the
    # sole value it ever saw was NSCLC's 4.47%. The pan-cancer counterpart is
    # 16.0%, which is exactly the case it was blind to.
    ("partition_sd", r"partition[- ](?:choice )?sd[^\d\n]{0,12}(\d+\.\d+)", 5e-5,
     "NSCLC partition-choice sd", "PARA:NSCLC!!pan-cancer"),
    ("pancancer_partition_sd",
     r"partition[- ](?:choice )?sd[^\d\n]{0,12}(\d+\.\d+)", 5e-5,
     "pan-cancer partition-choice sd", "PARA:pan-cancer!!NSCLC"),
    # "bootstrap SE" in the checklist and on the poster, "bootstrap standard
    # error" in the manuscript. The narrow form matched only the former, so the
    # paper's two bootstrap SEs were never checked; broadening reaches both.
    ("bootstrap_se",
     r"bootstrap (?:SE|standard error)[^\d\n]{0,12}(\d+\.\d+)", 5e-5,
     "NSCLC bootstrap SE", "PARA:NSCLC!!pan-cancer"),
    ("pancancer_bootstrap_se",
     r"bootstrap (?:SE|standard error)[^\d\n]{0,12}(\d+\.\d+)", 5e-5,
     "pan-cancer bootstrap SE", "PARA:pan-cancer!!NSCLC"),
    ("interval_understated_pct",
     r"(?:is|by)\s+\**(\d+\.\d+)\s*%\s*\**\s*too narrow", 5e-2,
     "how much the NSCLC reported interval understates, in percent",
     "PARA:NSCLC!!pan-cancer"),
    ("pancancer_interval_understated_pct",
     r"(?:is|by)\s+\**(\d+\.\d+)\s*%\s*\**\s*too narrow", 5e-2,
     "how much the pan-cancer reported interval understates, in percent",
     "PARA:pan-cancer!!NSCLC"),
    ("partition_var_pct",
     r"partition choice accounts for\s+\**(\d+\.\d+)\s*%", 5e-2,
     "share of total variance from partition choice, NSCLC",
     "PARA:NSCLC!!pan-cancer"),
    ("pancancer_partition_var_pct",
     r"partition choice accounts for\s+\**(\d+\.\d+)\s*%", 5e-2,
     "share of total variance from partition choice, pan-cancer",
     "PARA:pan-cancer!!NSCLC"),
    # A2: the rotation-null family p. Guarded to lines that are about the
    # rotation null, so the pattern can be loose enough to catch "p = 0.5".
    ("rotation_p", r"\bp\s*(?:=|&nbsp;=&nbsp;)\s*\**\s*(\d?\.\d+)", 2e-6,
     "rotation-null family p-value",
     r"rotation|B\s*=\s*1,?000|1,000 draws|B = 1000"),
    ("rotation_observed",
     r"[Oo]bserved family mean r[^\d\n]{0,12}\**(\d?\.\d+)", 5e-5,
     "observed family mean r", None),
    # NOT CHECKED: the rotation-null draw count B. `pipeline/README.md:110`
    # legitimately discusses "at B=100 the smallest attainable p (1/101)" as a
    # counterfactual, and no guard separates that from a claim without becoming
    # a regex for one line of one file. B is a round integer in a single
    # sentence; it is not where drift happens.
    # A3: the effective number of independent tests.
    ("m_eff",
     r"M(?:_|<sub>)?eff(?:</sub>)?\s*(?:=|&nbsp;=&nbsp;)\s*\**\s*(\d?\.\d+)",
     5e-3, "effective number of independent tests", None),
    # A4/A5 gap closures.
    # FOUND DEAD 2026-09-05, by giving it the self-test row it never had. The
    # pattern was `calibrat\w+...` guarded on /[Cc]ontrol C/, and the value it
    # is meant to check lives in a markdown table row of 14-SCIENCE-AUDIT.md:
    #
    #     | **Calibrated** | **0.0410** (p = 0.00995) |
    #
    # which does not contain the words "Control C" -- those are in the section
    # heading. The guard excluded the only line carrying the number, so this
    # pattern matched NOTHING, anywhere, and had done so since it was written.
    # It is the same line-guard-versus-wrapping defect handoff #24 recorded for
    # the small-panel bracket, sitting undetected in a DIFFERENT pattern,
    # invisible because the claim count went up either way.
    ("control_c_calibrated_r2",
     r"\*\*Calibrated\*\*\s*\|\s*\*\*(\d?\.\d{3,})\*\*", 5e-4,
     "Control C, calibrated by permutation (R^2)", None),
    # Guarded on `k=10` specifically, not on "panel size": 09-PAPER-DRAFT.md:337
    # quotes a DIFFERENT Spearman (rho = -0.50, the gap-versus-panel-size
    # correlation) two lines above the small-panel bracket's -1.000. A looser
    # guard flagged that correct number as drift. Note the character class
    # accepts U+2212 MINUS SIGN as well as ASCII hyphen -- the documents use the
    # typographic minus, so `-?` alone captured "0.50" from "-0.50".
    ("small_panel_spearman",
     r"Spearman[^\d\n−-]{0,14}([−-]?\d?\.\d+)", 5e-4,
     "small-panel monotonicity (Spearman)", r"k\s*=\s*10"),
    # Cohort sizes, each anchored on the phrase that names the cohort.
    ("pancancer_n", r"([\d,]{3,7})\s+TCGA patients", 0.5,
     "pan-TCGA cohort size", None),
    ("pancancer_n", r"Pan-TCGA[^\n]{0,12}?\(?\s*n\s*=\s*([\d,]{3,7})", 0.5,
     "pan-TCGA cohort size", None),
    ("pancancer_n", r"Pan-TCGA\s+([\d,]{5,7})\s+patients", 0.5,
     "pan-TCGA cohort size", None),
    ("pancancer_n_types", r"([\d,]{1,3})\s+(?:cancer|tumour|tumor)\s+types", 0.5,
     "number of cancer types", None),
    ("nsclc_n", r"NSCLC[^\n]{0,12}?\(?\s*n\s*=\s*([\d,]{3,7})", 0.5,
     "NSCLC cohort size", None),
    ("nsclc_n", r"([\d,]{3,7})\s+non.{0,2}small[- ]cell lung", 0.5,
     "NSCLC cohort size", None),
    ("pancancer_n_events", r"([\d,]{3,7})\s+PFI events", 0.5,
     "pan-TCGA PFI events", None),
    # "619 tissue source sites" -- stated in 09-PAPER-DRAFT.md and poster.html,
    # sourced by nothing until now. Recomputed live; see `_site_counts()`.
    ("pancancer_n_sites", r"([\d,]{2,6})\s+tissue source sites", 0.5,
     "pan-TCGA distinct tissue source sites", None),
    # Mechanically measured rather than stored: these are the two numbers that
    # have actually drifted in this project's history.
    # `\w+\s+` allows an adjective: 09-PAPER-DRAFT.md writes "78 automated
    # tests", which the adjacent-word-only pattern silently did not match. That
    # left the manuscript's own test count unchecked -- and stale.
    ("test_count", r"\**(\d{2,4})\**\s+(?:\w+\s+)?tests\b", 0.5,
     "tests pytest collects right now", r"test"),
    ("abstract_chars", r"\**([\d,]{3,6})\**\s*/\s*\**2,?600\**", 0.5,
     "abstract characters, counted by 08_count_abstract.py just now", None),
    ("abstract_headroom", r"(\d{1,4})\s+characters? of headroom", 0.5,
     "abstract characters remaining (limit minus live count)", None),
    # The public snapshot's size, from a live rebuild. Added 2026-09-07 because
    # this pair drifted four times in two sessions while three documents quoted
    # it and nothing checked it. Both patterns are anchored on the word
    # "stages"/"staged" and the unit, NOT on the digits -- anchoring on the
    # digits is the mistake that made three A7 patterns vacuous.
    #
    # The LIVE pushed figures (169 files / 1,447,396 bytes) must NOT match
    # these: they describe a different tree, one that no rebuild reproduces, and
    # they are correct where they appear. Both are excluded by requiring the
    # verb "stages"/"a rebuild stages" rather than matching any nearby number.
    # Both are WRAP-scoped: in BOTH documents that carry them the count and the
    # byte figure are split across a line break, so a line-scoped pattern would
    # match nothing and pass while checking nothing -- the exact failure this
    # file has shipped four times.
    ("snapshot_files", r"rebuild stages \**(\d{2,4})\**\s+files", 0.5,
     "files a snapshot rebuild stages right now", "WRAP:rebuild stages"),
    #
    # THE BYTE SIZE IS DELIBERATELY NOT CHECKED HERE, and this is a decision,
    # not an omission. It was written as a claim first and removed after one
    # measurement: `pipeline/scripts/` and `09-PAPER-DRAFT.md` are BOTH staged,
    # so the staged byte total moves when any source file or the manuscript
    # changes by a single character -- including this comment, and including the
    # line that would record the number. A gate on it would go red after almost
    # every edit in the project and demand a manuscript change to clear, which
    # is how a check earns being ignored. The FILE COUNT has the property the
    # byte total lacks: it moves only when a file is added or removed, which is
    # a real event and exactly the drift that went unnoticed (169 -> 170 -> 172
    # -> 173). The byte figure is instead written in the documents as a DATED
    # observation, the same fix applied to the two hand-maintained countdowns:
    # a stale dated observation reads as history, a stale bare claim reads as a
    # fact. `authority()` still measures it, so it is one line from being a gate
    # if that judgement ever changes.
    # A7. The split-half reliabilities are the load-bearing numbers of the
    # scorer result -- the whole "the estimand is undefined under ssGSEA"
    # argument rests on the null reliability collapsing from 0.799 to 0.258 --
    # and both were quoted in prose and checked by nothing.
    #
    # EVERY pattern here is anchored on the SURROUNDING WORDS and captures any
    # decimal. The first draft anchored on the digits instead (`0\.2[56]\d+`)
    # and `--self-test` proved three of them vacuous: injecting a wrong value
    # moved it outside the pattern that was supposed to catch it, so the row
    # reported nothing and passed by construction. That is the same defect the
    # module docstring records for `\b78 tests\b`, reintroduced. A pattern must
    # not depend on the value being right.
    ("a7_meanz_rel_null", r"from\s+(\d?\.\d{3,})\s+under mean-z", 5e-4,
     "A7 split-half reliability, mean-z null", None),
    ("a7_ssgsea_rel_null", r"under mean-z to\s+(\d?\.\d{3,})", 5e-4,
     "A7 split-half reliability, ssGSEA null", None),
    ("a7_ssgsea_rel_null", r"null reliability\s+(\d?\.\d{3,})", 5e-4,
     "A7 split-half reliability, ssGSEA null (mechanism sentence)", None),
    ("a7_ssgsea_rel_obs", r"split-half\s+(\d?\.\d{3,})", 5e-4,
     "A7 split-half reliability, ssGSEA observed", r"ssGSEA"),
    ("a7_ssgsea_rel_obs", r"against\s+(\d?\.\d{3,})\s+observed", 5e-4,
     "A7 split-half reliability, ssGSEA observed (mechanism sentence)", None),
    # The disattenuation multiplier. If the reliability moves this must move
    # with it, and nothing tied them together before.
    ("a7_ssgsea_rel_null", r"1\s*/\s*√\s*(\d?\.\d{3,})", 5e-4,
     "A7 disattenuation divisor (ssGSEA null reliability)", None),
    # The two cross-scorer agreement statistics that reached the manuscript.
    ("a7_agree_unc_pearson", r"Pearson r\s*=\s*(\d?\.\d{3,})", 5e-4,
     "A7 cross-scorer agreement, uncorrected excess (Pearson)",
     r"scorer|agreement|ssGSEA"),
    ("a7_agree_unc_spearman", r"Spearman ρ\s*=\s*(\d?\.\d{3,})", 5e-4,
     "A7 cross-scorer agreement, uncorrected excess (Spearman)",
     r"scorer|agreement|ssGSEA"),
    ("a7_agree_resid_pearson", r"agree closely \(r\s*=\s*(\d?\.\d{3,})", 5e-4,
     "A7 cross-scorer agreement, residualised observed r (Pearson)", None),
    ("a7_agree_resid_absdiff", r"mean \|difference\|\s*(\d?\.\d{3,})", 5e-4,
     "A7 cross-scorer agreement, residualised observed r (mean |diff|)", None),

    # ---- added 2026-09-04: the soft spots handoff #23 listed ----------------
    #
    # Median image-signature correlation. GUARDED on "image-signature", because
    # the Results section quotes an identically shaped sentence 25 lines later
    # for the image-AXIS correlation ("0.067 pan-cancer and -0.023 in NSCLC").
    # Without the guard the axis numbers are captured and reported as drift in
    # the signature claim. The paper line was reflowed so the label and both
    # numbers share one line; this scanner is line-based.
    ("pancancer_median_r",
     r"correlation was\s+(\d?\.\d{2,})\s+pan-cancer", 5e-4,
     "median image-signature correlation, pan-TCGA", r"image[–-]signature"),
    ("nsclc_median_r",
     r"pan-cancer and\s+(\d?\.\d{2,})\s+in NSCLC", 5e-4,
     "median image-signature correlation, NSCLC", r"image[–-]signature"),
    # Median split-robustness Delta r. The patterns are deliberately tight
    # rather than guarded: the SAME LINE also carries the Delta-MAE triples
    # (0.0088 [0.0078, 0.0098]), which a loose pattern would swallow.
    ("pancancer_median_delta_r",
     r"median Δr\s*=\s*([+−-]?\d?\.\d{2,})", 5e-4,
     "median split-robustness Delta r, pan-TCGA", None),
    ("nsclc_median_delta_r",
     r"median Δr\s*=\s*[+−-]?\d?\.\d{2,}\s+and\s+([+−-]?\d?\.\d{2,})", 5e-4,
     "median split-robustness Delta r, NSCLC", None),
    # The ancestry arm's cohort size and its two identifiability statistics.
    ("ancestry_n", r"ancestry arm \(([\d,]{3,7})\s+patients", 0.5,
     "ancestry arm cohort size", None),
    ("ancestry_cramers_v", r"Cramér's V\s*=\s*(\d?\.\d{2,})", 5e-4,
     "ancestry x site Cramer's V", None),
    # WRAP-scoped: the paper wraps this as "with 95th / percentile 0.088", so
    # a per-line scan never saw it. The `\n?` in the old pattern is the tell --
    # it was written for a scan that had no newline to give it, and the claim
    # went unchecked from the day it was added.
    ("ancestry_perm_p95", r"95th\s+percentile\s+(\d?\.\d{2,})", 5e-4,
     "ancestry x site Cramer's V permutation null, 95th percentile",
     "WRAP:Cram"),
    # Reliability quoted at three decimals in the manuscript prose. The
    # full-precision copies in 14-SCIENCE-AUDIT.md are already covered by the
    # LONG_DECIMAL rule; these rounded ones were not covered by anything.
    ("a7_ssgsea_rel_obs", r"split-half\s+(\d?\.\d{3,})", 5e-4,
     "ssGSEA observed split-half reliability", None),
    ("a7_meanz_rel_null", r"from\s+(\d?\.\d{3,})\s+under mean-z", 5e-4,
     "mean-z null reliability (prose, 3dp)", None),
    ("a7_ssgsea_rel_null", r"under mean-z to\s+(\d?\.\d{3,})", 5e-4,
     "ssGSEA null reliability (prose, 3dp)", None),

    # ---- added 2026-09-05: reliability, decomposition, label-side, immune p --
    #
    # EVERY anchor below was chosen by READING THE DOCUMENT'S LINE BREAKS, not
    # by picking the most natural-sounding phrase. Handoff #24's defect was a
    # guard whose anchor and number sat on different lines, and this scanner is
    # line-based, so an anchor is only usable if it CANNOT be separated from its
    # number by a wrap. Where the paper's own wrapping made that impossible the
    # paragraph was reflowed (see 09-PAPER-DRAFT.md "fell only from 0.984 to
    # 0.976", which used to break between "to" and "0.976").
    #
    # The alpha LEVELS are medians; the alpha GAPS are means. See authority().
    ("pancancer_alpha_null_raw", r"median α\s+(\d?\.\d{2,})", 5e-4,
     "median raw Cronbach alpha of random sets, pan-TCGA", None),
    ("pancancer_alpha_null_resid", r"null's α to\s+(\d?\.\d{2,})", 5e-4,
     "median residualised alpha of random sets, pan-TCGA", None),
    # The manuscript sentence was REWORDED to make these two checkable. It read
    # "...curated signatures fell only from 0.984 to 0.976." and wrapped
    # between "only" and "from", putting the anchor on one line and both
    # numbers on the next; `--self-test` duly reported both patterns MISSED.
    # The prose now reads "their α fell from 0.984 to 0.976" on a single line,
    # so the anchor cannot be separated from its numbers by a reflow.
    ("pancancer_alpha_obs_raw", r"α fell from\s+(\d?\.\d{2,})", 5e-4,
     "median raw alpha of curated signatures, pan-TCGA", None),
    ("pancancer_alpha_obs_resid",
     r"α fell from\s+\d?\.\d{2,}\s+to\s+(\d?\.\d{2,})", 5e-4,
     "median residualised alpha of curated signatures, pan-TCGA", None),
    # The raw-versus-residualised reliability gaps, and the fold difference the
    # sentence draws from them. These four numbers are the quantitative core of
    # the paper's methodological claim and nothing checked them.
    ("pancancer_rel_gap_raw", r"understated:\s*\+?(\d?\.\d{2,})\s+pan-cancer",
     5e-4, "mean raw-score reliability gap, pan-TCGA", None),
    ("nsclc_rel_gap_raw",
     r"understated:[^\n]*?pan-cancer and\s+\+?(\d?\.\d{2,})\s+NSCLC", 5e-4,
     "mean raw-score reliability gap, NSCLC", None),
    ("pancancer_rel_gap_resid",
     r"\*\*\+?(\d?\.\d{2,})\s+and\s+\+?\d?\.\d{2,}\*\*\s+when computed", 5e-4,
     "mean residualised reliability gap, pan-TCGA", None),
    ("nsclc_rel_gap_resid",
     r"\*\*\+?\d?\.\d{2,}\s+and\s+\+?(\d?\.\d{2,})\*\*\s+when computed", 5e-4,
     "mean residualised reliability gap, NSCLC", None),
    ("pancancer_rel_gap_fold", r"a\s+(\d{1,3}(?:\.\d+)?)-fold and", 5e-2,
     "residualised-over-raw reliability gap ratio, pan-TCGA", None),
    ("nsclc_rel_gap_fold", r"-fold and\s+(\d{1,3}(?:\.\d+)?)-fold", 5e-2,
     "residualised-over-raw reliability gap ratio, NSCLC", None),
    ("pancancer_rel_gap_angio", r"reaches\s+\+?(\d?\.\d{2,})\s+for the smallest",
     5e-4, "reliability gap of the 36-gene angiogenesis set, pan-TCGA", None),
    # Median residualised correlations and their null means.
    # SPELLING-AGNOSTIC as of 2026-09-06. This pattern was anchored on the
    # British "residualised"; converting the manuscript to US spelling silently
    # killed it -- the claim count fell 345 -> 344 and --self-test reported the
    # row MISSED, which is the only reason it was noticed. Any pattern whose
    # anchor is a word with a British/US variant must accept both, or the next
    # spelling decision breaks it again.
    ("pancancer_r_resid_refit",
     r"residuali[sz]ed correlations were\s+(\d?\.\d{2,})", 5e-4,
     "median residualised image-signature correlation, pan-TCGA", None),
    ("nsclc_r_resid_refit", r"and\s+(\d?\.\d{2,})\s+against null means", 5e-4,
     "median residualised image-signature correlation, NSCLC", None),
    ("pancancer_null_mean_r", r"null means of\s+(\d?\.\d{2,})", 5e-4,
     "median null-set correlation, pan-TCGA", None),
    ("nsclc_null_mean_r",
     r"null means of\s+\d?\.\d{2,}\s+and\s+(\d?\.\d{2,})", 5e-4,
     "median null-set correlation, NSCLC", None),
    # Median image-AXIS correlation. Anchored on "image-axis" so it cannot be
    # confused with the image-SIGNATURE sentence 25 lines earlier, which has an
    # identical shape and is checked by its own pair of patterns above.
    ("pancancer_r_axis",
     r"image[–-]axis correlation was\s+([−+-]?\d?\.\d{2,})", 5e-4,
     "median image-axis correlation, pan-TCGA", None),
    ("nsclc_r_axis",
     r"image[–-]axis correlation was[^\n]*?pan-cancer and\s+\*{0,2}([−+-]?\d?\.\d{2,})",
     5e-4, "median image-axis correlation, NSCLC", None),
    # Purity decomposition and the covariate-only baseline.
    ("pancancer_r_partial", r"was\s+(\d?\.\d{2,})\s+pan-cancer\s*\(against", 5e-4,
     "median purity-partial correlation, pan-TCGA", None),
    ("pancancer_r_unadjusted", r"\(against\s+(\d?\.\d{2,})\s+unadjusted", 5e-4,
     "median unadjusted image-signature correlation, pan-TCGA", None),
    ("pancancer_incr_r2", r"stage was\s+(\d?\.\d{2,})\s+and", 5e-4,
     "incremental adjusted R^2 of the image, pan-TCGA", None),
    ("nsclc_incr_r2", r"stage was\s+\d?\.\d{2,}\s+and\s+(\d?\.\d{2,})", 5e-4,
     "incremental adjusted R^2 of the image, NSCLC", None),
    ("pancancer_r_covariates", r"reached median r\s*=\s*(\d?\.\d{2,})", 5e-4,
     "median covariate-only baseline correlation, pan-TCGA", None),
    # "**-0.000** over it" -- anchored on the trailing "over it" because the
    # preceding "the embedding added" wraps onto the previous line.
    ("pancancer_emb_over_cov", r"\*{0,2}([−+-]?\d\.\d{3,})\*{0,2}\s+over it",
     5e-4, "embedding gain over covariates, pan-TCGA", None),
    ("nsclc_emb_over_cov", r"embedding added\s+\+?([−+-]?\d?\.\d{2,})", 5e-4,
     "embedding gain over covariates, NSCLC", None),
    # The label-side site-variance control, quoted as percentages.
    ("label_site_variance_pct",
     r"median\s+(\d{1,3}(?:\.\d+)?)\s*%\s+of the label", 5e-2,
     "site-only R^2 of the label, pan-TCGA (percent)", None),
    ("label_site_given_type_pct",
     r"reduces it to\s+(\d{1,3}(?:\.\d+)?)\s*%", 5e-2,
     "site-given-cancer-type R^2 of the label (percent)", None),
    ("control_c_calibrated_pct",
     r"gives\s+(\d{1,3}(?:\.\d+)?)\s*%", 5e-2,
     "Control C, permutation-calibrated (percent)", r"calibrat"),
    ("label_plate_within_site", r"explains\s+(\d\.\d{3})", 5e-4,
     "plate-within-site R^2 of the label", None),
    # The immune-versus-other outcome contrast. Recomputed, not stored; see
    # `_immune_contrast()` for why both tails of both tests are in the table.
    # REPOINTED 2026-09-06 from the one-sided authority to the two-sided one,
    # when the author chose two-sided for both tests. The patterns did not
    # change; only the value they are held to did, and every one of them went
    # RED the moment the documents were edited -- which is the check working.
    ("immune_mwu_p_2s", r"Mann[–—-]Whitney p\s*=\s*(\d?\.\d{2,})", 5e-4,
     "immune-vs-other outcome contrast, two-sided Mann-Whitney p", None),
    ("immune_fisher_p_2s", r"Fisher p\s*=\s*(\d?\.\d{2,})", 5e-4,
     "immune-vs-other outcome contrast, two-sided Fisher p", None),
    ("immune_mwu_p_2s", r"\*P\*\s*=\s*(\d?\.\d{2,})", 5e-4,
     "immune-vs-other outcome contrast, two-sided Mann-Whitney p (abstract)",
     r"MSigDB-immune|signatures;"),
    ("immune_mwu_p_2s", r"contrast \(p\s*=\s*(\d?\.\d{2,})\s*→", 5e-4,
     "immune-vs-other contrast p, IL6 excluded (Figure 3 caption)", None),
    ("immune_mwu_p_2s_il6", r"→\s*(\d?\.\d{2,})\)", 5e-4,
     "immune-vs-other contrast p, IL6 counted immune (Figure 3 caption)", None),

    # ---- 07-ABSTRACT-DRAFT.md's own numbers ---------------------------------
    #
    # The abstract is the SUBMITTED artefact, and `--coverage` reported it at
    # 6% checked -- the worst of the eleven documents. Its headline ISI is now
    # covered by PAREN_TRIPLE; these four cover the rest of its Results
    # sentence. The paragraph was REFLOWED (words and order untouched, line
    # breaks moved) so that "0.318 (0.261-0.376)" and "0.45 at 10" each sit on
    # one line; before that they straddled a break and no line-based pattern
    # could reach them. The character count is unchanged and is re-verified
    # live by `_abstract_chars()` on every run.
    ("immune_mwu_p_2s", r"\bp\s*=\s*(\d?\.\d{2,})", 5e-4,
     "immune-vs-other outcome contrast p (abstract, plain form)",
     r"MSigDB-immune"),
    ("pancancer_rel_gap_fold",
     r"understated it\s+(\d{1,3}(?:\.\d+)?)-fold", 5e-2,
     "residualised-over-raw reliability gap ratio (abstract)", None),
    # WIDENED 2026-09-06. These two were anchored on the literal word "alpha",
    # which the ABSTRACT writes but the PAPER does not -- 09-PAPER-DRAFT.md:61
    # and :103 both write "Cronbach's α from 0.97 to 0.80". So the manuscript's
    # reliability headline, the number the whole Reliability paragraph turns on,
    # was reported UNCHECKED by --coverage while the abstract's copy of the same
    # claim was checked. Accepting either spelling covers both. Verified against
    # every occurrence in the tree first: 07:62, 09:61, 09:103 and
    # 13-NEXTGEN:80 all quote the SAME pan-cancer pair, so one authority is
    # correct for all of them and no false mismatch is possible.
    ("pancancer_alpha_null_raw", r"(?:alpha|α) from\s+(\d?\.\d{2,})", 5e-3,
     "median raw alpha of random sets (abstract + paper, 2dp)", None),
    ("pancancer_alpha_null_resid",
     r"(?:alpha|α) from\s+\d?\.\d{2,}\s+to\s+(\d?\.\d{2,})", 5e-3,
     "median residualised alpha of random sets (abstract + paper, 2dp)", None),

    # ---- added 2026-09-06, worked from `--coverage`, not from reading ----
    # Every pattern below was written against the document's ACTUAL line
    # wrapping, checked by printing the target lines first. That is the defect
    # that shipped twice before: an anchor that sits on the far side of a line
    # break from its number matches nothing and reports nothing.
    #
    # The MDE. A pattern anchored on "minimum detectable effect" already
    # existed and matches 09-PAPER-DRAFT.md:457 -- but 07-ABSTRACT-DRAFT.md
    # wraps between "minimum" and "detectable" (":17-18"), so the abstract's
    # copy was invisible to it. Anchoring on the line-local "detectable effect
    # is ... in C-index" reaches both.
    ("nsclc_mde", r"detectable effect is\s+(\d?\.\d+)\s+in C-index", 5e-4,
     "NSCLC minimum detectable effect, C-index (line-local form)", None),
    # Cohort sizes and type count, in the Methods sentence of the abstract.
    ("pancancer_n", r"([\d,]+)\s+patients across\s+\d+\s+cancer types", 0.5,
     "pan-cancer cohort size (\"N patients across M cancer types\")", None),
    # NOT ADDED: a `pancancer_n_types` twin of the pattern above. The existing
    # `([\d,]{1,3}) (?:cancer|tumour|tumor) types` already reaches every
    # occurrence, including this one -- confirmed by the fact that --coverage
    # reports 07-ABSTRACT-DRAFT.md:47's "7,168" as unchecked but NOT its "31".
    # A second pattern for the same literal would add a row to the claim table
    # and check nothing new, which is precisely the inflation this file warns
    # about.
    ("nsclc_n", r"\(n\s*=\s*([\d,]+)\)\s+as a second cohort", 0.5,
     "NSCLC cohort size, as stated in the abstract's Methods", None),
    # The global-axis cost. 2dp in the abstract against a stored 0.0673, so the
    # decimals-aware tolerance (0.005) carries it; the 5e-4 floor never binds.
    ("pancancer_r_axis", r"image-axis r\s*=\s*(\d?\.\d+)", 5e-4,
     "pan-cancer image-to-global-axis correlation", None),
    # The small-panel scaling, quoted at 2dp in the abstract and 3dp in the
    # paper and cover letter. One pattern serves all three because _tol_for
    # derives the tolerance from the decimals each document actually wrote.
    ("small_panel_k160", r"(\d?\.\d{2,})\s+at 160 genes", 5e-4,
     "reliability gap at k=160", None),
    ("small_panel_k10", r"at 160 genes to\s+(\d?\.\d{2,})\s+at 10", 5e-4,
     "reliability gap at k=10", None),
    # Anchored to "at 10 (Spearman", NOT to a bare "Spearman". Two unscanned
    # documents (00-PROJECT-BRIEF.md:255, 08-READINESS-ASSESSMENT.md:102) quote
    # a DIFFERENT Spearman -- the set-size scaling, -0.50/-0.59 -- and a loose
    # anchor would report those as drift the day either file joins DOCS.
    # 09-PAPER-DRAFT.md:107 wraps directly after "(Spearman", so only the
    # abstract's copy is reachable by a line-based pattern; that is a limit of
    # the approach, not a miss.
    ("small_panel_spearman", r"at 10\s*\(Spearman\s+([−-]?\d?\.\d+)", 5e-4,
     "Spearman of reliability gap against panel size", None),
    # A9's split-determinism finding. These are claims about a RESULT -- the
    # tie structure is what makes the partition platform-dependent, and the
    # number of patients that move is the size of the defect. Both are read
    # from results/split_determinism_darwin_arm64.json, not from memory.
    # SPELLING-AGNOSTIC anchors are unnecessary here (no British/US variants),
    # but the anchors deliberately avoid the digits of the correct value.
    ("splitdet_macos_sites_tied", r"(\d{2,3}) of the 68 NSCLC sites", 0.5,
     "NSCLC sites sharing a size with another site", None),
    # WRAP-scoped for the same reason: Limitations 8 wraps this as "216 of /
    # 944 patients change fold".
    ("splitdet_macos_moved", r"(\d{2,3}) of\s+944 patients change fold", 0.5,
     "patients assigned a different fold under the stable sort",
     "WRAP:change fold"),
    # A10, the ssGSEA tie convention. Claims about a RESULT: the tie structure
    # is what exposes the rank table to the defect, and the score difference is
    # its size. Read from results/sort_audit_darwin_arm64.json.
    #
    # Every anchor below is prose, never the digits of the correct value, and
    # every one is checked against where the number actually SITS after line
    # wrapping -- a wrap-spanning anchor has silently matched nothing three
    # times in this file's history, so `--self-test` is the acceptance test for
    # these rows, not a reading of the regex.
    ("sortaudit_macos_tied_elements",
     r"([\d,]+) of 38,747,424\s+rank entries", 0.5,
     "tied rank entries in the NSCLC expression matrix", None),
    ("sortaudit_macos_cells_diff",
     r"different scores for ([\d,]+) of 15,104", 0.5,
     "patient-signature cells whose ssGSEA score depends on the tie order", None),
    ("sortaudit_macos_max_diff",
     r"difference between the two orderings is ([\d.]+)", 0.005,
     "max ssGSEA score difference between tie conventions", None),
    ("sortaudit_macos_score_sd",
     r"per-signature score standard deviation of ([\d.]+)", 0.005,
     "mean per-signature ssGSEA score sd", None),

    # A10's ssGSEA re-run under the pinned ordering. The point estimate appears
    # in three places with three different framings, so each is matched on its
    # own anchor rather than on a bare "+0.174", which would also match a stale
    # copy of the same number somewhere else. `WRAP:` matches the
    # whitespace-normalised paragraph: the paper wraps "+0.174 [0.133, 0.210]"
    # across a line break, which is precisely how two earlier patterns went
    # vacuous.
    ("a10ss_ssgsea_unc",
     r"moves from \+0\.180 to\s+\+(\d?\.\d{2,})", 5e-4,
     "A10 ssGSEA uncorrected contrast, pinned tie ordering (Limitations 8)",
     r"WRAP:recomputed under the pinned\s+ordering"),
    ("a10ss_ssgsea_unc",
     r"does not propagate to this contrast:\s*\+(\d?\.\d{2,})",
     5e-4, "A10 ssGSEA uncorrected contrast, pinned (Results)",
     "WRAP:does not propagate to this contrast"),
    ("a10ss_ssgsea_unc",
     r"\+0\.180; \+(\d?\.\d{2,}) pinned", 5e-4,
     "A10 ssGSEA uncorrected contrast, pinned (poster panel 6)", None),
    ("a10ss_shift",
     r"a shift of (\d?\.\d{2,}), or", 5e-4,
     "A10 ssGSEA shift under the pinned ordering", None),
    ("a10ss_ssgsea_unc",
     r"and at \+(\d?\.\d{2,}) under the pinned tie ordering", 5e-4,
     "A10 pinned ssGSEA contrast (Limitations 9)",
     "WRAP:and at"),
    ("a10ss_shift_pct_of_halfwidth",
     r"or (\d{1,2})% of the pre-fix interval", 1.0,
     "A10 ssGSEA shift as a percent of the pre-fix CI half-width", None),
    ("a10ss_frac_of_meanz_pct",
     r"mean-z magnitude \((\d{2})% against", 1.0,
     "A10 pinned ssGSEA as a percent of the mean-z uncorrected contrast", None),
    # The pre-fix fraction, in the same sentence and previously unchecked. Both
    # halves of a "54% against 56%" comparison have to be anchored, or the
    # sentence can drift on the half nobody is watching.
    ("a7_frac_of_meanz_pct",
     r"mean-z magnitude \(\d{2}% against (\d{2})%", 1.0,
     "A7 pre-fix ssGSEA as a percent of the mean-z uncorrected contrast", None),
    ("a7_frac_of_meanz_pct",
     r"survives the change\s+of scorer at roughly (\d{2})% of its magnitude", 1.0,
     "A7 pre-fix ssGSEA fraction of mean-z, as stated in the Results heading",
     "WRAP:survives the change of scorer at roughly"),
    ("a10ss_beats_ssgsea_unc",
     r"same (\d+) of 16 signatures beating their null", 0.5,
     "A10 pinned ssGSEA signatures beating null (Results)", None),
    ("a10ss_beats_ssgsea_unc",
     r"\+0\.174, (\d+) of 16 signatures beat their null", 0.5,
     "A10 pinned ssGSEA signatures beating null (Limitations 8)", None),
    # The A10 comparison table in 14-SCIENCE-AUDIT.md. Worded "above their null"
    # rather than "beating their null" ON PURPOSE: the FRACTIONS net below fires
    # on "beat|exceed|beating" and requires 16/16, so the ssGSEA arm's genuine
    # 15/16 read as a mismatch there. Rewording keeps the coarse net strict for
    # the primary instead of teaching it to accept 15/16 everywhere, and these
    # two patterns then check the values the rewording moved out of its reach.
    ("a7_beats_ssgsea_unc", r"(\d+) of 16 pre-fix", 0.5,
     "A7 pre-fix ssGSEA signatures above null (audit A10 table)", None),
    ("a10ss_beats_ssgsea_unc", r"(\d+) of 16 pinned", 0.5,
     "A10 pinned ssGSEA signatures above null (audit A10 table)", None),
    # A11's corrected NSCLC Control C. Anchored on "median R²" rather than on a
    # bare decimal, because the pan-cancer value sits a few lines away in the
    # same document and a loose pattern would happily check one against the
    # other -- which is the exact confusion A11 was.
    ("nsclc_label_site_variance_r2",
     r"median R² (\d?\.\d{3,})", 5e-4,
     "A11 corrected NSCLC label-side site variance, median R^2", None),
    ("nsclc_gaps_n_patients_used",
     r"([\d,]+) of 7,168 patients, 31 sites", 0.5,
     "A11 corrected NSCLC science-gaps patient count", None),
    # The pan-cancer outcome effect on the C-index scale. Anchored on the
    # phrase, not the digits, and deliberately NOT on the neighbouring 0.046
    # MDE -- the two sit in the same sentence and a loose anchor would match
    # the wrong one.
    # ANCHOR AND NUMBER ON ONE LINE. The first version of this pattern read
    # "real pan-cancer effect has a\s+median of (...)" and spanned a line wrap,
    # so it matched nothing and --self-test reported it VACUOUS. That is the
    # third time a wrap-spanning anchor has been shipped in this file.
    # The poster's ssGSEA arm. Added when the poster's "ssGSEA sensitivity not
    # run" line was found STALE -- A7 closed 2026-09-03 and the run is frozen.
    # A stale claim is exactly what this checker cannot see unless the number
    # it replaces is given a pattern, so the replacement got one.
    ("a7_ssgsea_unc", r"uncorrected contrast \(\+(\d?\.\d{2,})", 5e-4,
     "ssGSEA uncorrected excess, A7 scorer sensitivity", None),
    ("pancancer_outcome_delta_c_median",
     r"median C-index excess of\s+(\d?\.\d{2,})", 5e-4,
     "median C-index excess over null, signatures that beat null, pan-TCGA",
     None),
]

# Number-WORDS. Every pattern above matches digits, so "all thirteen
# re-verified" survived in the manuscript for a whole session after the
# reference list dropped to twelve. A word is a claim exactly as much as a
# digit is.
#
# (key, regex whose capture is a number-word, description, line-guard or None)
NUMBER_WORDS: dict[str, int] = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
    "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11, "twelve": 12,
    "thirteen": 13, "fourteen": 14, "fifteen": 15, "sixteen": 16,
    "seventeen": 17, "eighteen": 18, "nineteen": 19, "twenty": 20,
    "thirty": 30, "thirty-one": 31, "thirty-two": 32,
}

WORD_CLAIMS: list[tuple[str, str, str, str | None]] = [
    # The reference count. This is the exact claim that went stale: the list
    # dropped from thirteen to twelve and the prose still said "all thirteen".
    # DELIBERATELY UNGUARDED. The obvious guard is /[Rr]eference/, and it is
    # wrong for exactly the reason handoff #24 recorded: the sentence is
    # "**Bibliographic metadata for all twelve re-verified against PubMed on
    # 2026-09-03**", and the nearest word "References" is the SECTION HEADING
    # two lines above it. A line-guard cannot see it. The pattern is made
    # specific by requiring the "(re-)verified" verb instead, and an unknown
    # capture is skipped rather than failed.
    ("n_references", r"\b(?:all|the)\s+([a-z]+(?:-[a-z]+)?)\s+(?:re-)?verified",
     "reference count stated as a word", None),
    ("n_references", r"\b([a-z]+(?:-[a-z]+)?)\s+references,\s*all cited",
     "reference count stated as a word", None),
    # REMOVED 2026-09-05, one commit after it was added, because it produced a
    # FALSE POSITIVE on correct prose. The row read:
    #
    #   ("pancancer_n_sigs", r"\b([a-z]+)\s+...signatures", guard=Hallmark|MSigDB)
    #
    # and it fired on "Six signatures are MSigDB-immune", reporting 6 against
    # the expected 16. The sentence is right; the pattern cannot tell "the 16
    # Hallmark signatures" from "six immune signatures", and no guard fixes
    # that, because both sentences legitimately contain the same anchor words.
    #
    # It is deleted rather than tightened. A pattern narrowed until it matches
    # only the phrasing that exists today matches NOTHING the moment the prose
    # is reworded, and a row that matches nothing inflates the apparent size of
    # this table while checking as much as the vacuous patterns this file's
    # docstring warns about. The reference-count rows above stay: they target a
    # claim with exactly one phrasing and a self-test row proves each can fail.

]

# Values that must equal ONE OF a set, when a line quotes several cohorts at
# once and the order is not guaranteed -- "AUROC 0.998 pan-cancer and 0.992 in
# NSCLC". Checking each captured value against the SET is order-independent and
# still fails on a wrong number, which per-cohort anchoring could not do here.
#
# (keys, regex with one capture, tolerance, description, line-guard or None)
SET_CLAIMS: list[tuple[tuple[str, ...], str, float, str, str | None]] = [
    (("pancancer_site_auroc", "nsclc_site_auroc"),
     r"(?:AUROC|auroc)[^\d\n]{0,24}?(\d?\.\d{2,})", 5e-4,
     "median site AUROC (must match one cohort's median)", None),
    (("pancancer_mde", "nsclc_mde"),
     r"(?:MDE|minimum detectable effect)[^\d\n]{0,26}?(\d?\.\d{2,})", 5e-4,
     "minimum detectable effect in C-index", None),
    # The small-panel bracket: five median gaps quoted in one sentence, in
    # DESCENDING k, so per-key anchoring would need five near-identical
    # patterns. The set form is order-independent and still fails on a wrong
    # number.
    #
    # DELIBERATELY UNGUARDED, and this is the second time the reason has had to
    # be learned. The first version guarded on `k\s*=\s*160`, the sentence's
    # opening term. `--self-test` reported MISSED: the bracket sentence WRAPS,
    # so "0.084 at k=160" and "**0.453 at k=10**" land on different lines, and
    # this scanner is line-based -- the guard excluded the very line carrying
    # the injected value. A line-guard is only safe when it cannot be separated
    # from its number by a line break, which is a property of the DOCUMENT's
    # wrapping, not of the claim. The pattern below is instead made specific by
    # its shape: exactly three decimals (so 0.3182 cannot match) immediately
    # followed by " at <digit>" or " at k=<digit>".
    (("small_panel_k10", "small_panel_k20", "small_panel_k40",
      "small_panel_k80", "small_panel_k160"),
     # The trailing size is pinned to the five PANEL SIZES, not left as any
     # digit. Unpinned, `\d\.\d{3} at \d` also matched "p = 0.001 at 1,000
     # draws" -- the rotation null, in two different documents. Note what is
     # anchored: k is a fixed design constant of the bracket, whereas the gap
     # value is captured and compared. Anchoring on the GAP's digits is the
     # vacuous-pattern defect; anchoring on k is not.
     # Two OR three decimals: the manuscript quotes the bracket at three
     # ("0.084 at k=160") and the abstract, which is character-limited, at two
     # ("0.08 at 160 genes"). `_tol_for` widens the tolerance to the rounding
     # half-width of whatever was written, so a 2-decimal quote is held to
     # 2 decimals and no tighter.
     r"(?<![\d.])(\d\.\d{2,3})(?![\d])\s+at\s+(?:k\s*=\s*)?(?:10|20|40|80|160)\b",
     5e-4,
     "small-panel bracket median gap (must match one of k=10/20/40/80/160)",
     None),
    # The four ancestry group counts, quoted as "EUR 5,312, AFR 594, ASIAN 502,
    # AMR 172". Set form for the same reason.
    (("ancestry_eur", "ancestry_afr", "ancestry_asian", "ancestry_amr"),
     r"(?:EUR|AFR|ASIAN|AMR)\s+([\d,]{3,6})", 0.5,
     "ancestry group count (must match one of EUR/AFR/ASIAN/AMR)", None),
    # Dashed RANGES: "excess 0.165-0.390 and 0.148-0.440", and the rotation
    # null's "(95% range 0.098-0.122)". Set form because one line carries four
    # bounds in a fixed order that per-key anchoring cannot address without
    # four near-identical patterns.
    #
    # The rotation-null bounds are in the SET deliberately, not because they
    # are the same quantity, but because the only guard that reaches BOTH
    # sentences is /excess|range/ -- the manuscript's provenance note wraps
    # between "mean excess" and its numbers, so an /excess/-only guard cannot
    # see line 21. Widening the guard to /range/ pulls in the rotation
    # sentence, so its bounds are listed rather than excluded. Every value in
    # the set is still read from a frozen artefact, so a wrong number in either
    # sentence still fails.
    (("pancancer_excess_min", "pancancer_excess_max",
      "nsclc_excess_min", "nsclc_excess_max",
      "rotation_null_lo", "rotation_null_hi"),
     r"(?<![\d.])(\d\.\d{3})(?![\d])\s*[–—-]\s*\d\.\d{3}(?![\d])", 5e-4,
     "low bound of a dashed range (excess spread or rotation null)",
     r"excess|range"),
    (("pancancer_excess_min", "pancancer_excess_max",
      "nsclc_excess_min", "nsclc_excess_max",
      "rotation_null_lo", "rotation_null_hi"),
     r"(?<![\d.])\d\.\d{3}(?![\d])\s*[–—-]\s*(\d\.\d{3})(?![\d])", 5e-4,
     "high bound of a dashed range (excess spread or rotation null)",
     r"excess|range"),
]

# Fractions written as N/M -- "16/16 signatures beat their null", "10/32".
# Each fraction found on a guarded line must equal one of the listed
# (numerator_key, denominator_key) pairs. A fraction matching none of them is a
# MISMATCH, so 15/16 or 9/32 fails.
#
# (line-guard, [(num_key, den_key), ...], description)
FRACTIONS: list[tuple[str, list[tuple[str, str]], str]] = [
    (r"signatures? (?:beat|exceed|beating)|ISI\s*\||beats null",
     [("pancancer_sigs_beat", "pancancer_n_sigs"),
      ("nsclc_sigs_beat", "nsclc_n_sigs")],
     "signatures beating their null"),
    (r"[Oo]utcome (?:arm|\(PFI)|PFI, within",
     [("pancancer_outcome_beat", "pancancer_n_outcome_tests"),
      ("nsclc_outcome_beat", "nsclc_n_outcome_tests")],
     "outcome-arm signatures beating their null"),
]

# Bracketed point-and-interval triples, e.g. `0.2911 [0.2692, 0.3133]`. This is
# the highest-value pattern in the project: it is how every headline is written,
# and it catches a CI bound that drifted away from its point estimate.
#
# The point estimate identifies WHICH quantity is being quoted; the two bounds
# are then checked against it. A triple whose point estimate matches nothing
# known is REPORTED, not failed -- HPC4 diagnostics and ssGSEA arms are
# legitimately not in the authority table.
# NOTE the U+2212 MINUS SIGN in every sign class. `[+-]?` alone matched no
# negative triple the documents actually contain, because they are typeset with
# the typographic minus, not the ASCII hyphen -- so the FIRST negative interval
# to reach the manuscript (A7's ssGSEA reconstruction) was invisible to the
# highest-value pattern in this file. `--self-test` caught it. The identical
# defect was found and fixed in the small-panel Spearman pattern earlier and
# was never propagated here; `_num` has always normalised U+2212, so only the
# regex was wrong.
TRIPLE = re.compile(
    r"(?<![\d.])([+\-−]?\d\.\d{3,})\s*\[\s*([+\-−]?\d\.\d{3,})\s*,\s*"
    r"([+\-−]?\d\.\d{3,})\s*\]")

# THE SAME TRIPLE, WRITTEN THE WAY THE ABSTRACT AND THE RESULTS ACTUALLY WRITE
# IT: `0.291 (95% CI, 0.269-0.313)`. Found by `--coverage` on 2026-09-05.
#
# This is the exact hole the coverage report was built to expose. `TRIPLE`
# requires SQUARE BRACKETS and a COMMA, and the manuscript's two most important
# sentences -- the abstract's headline and the Results' headline -- use round
# brackets and an en-dash instead. So the ISI and both its confidence bounds,
# in the two places a reviewer reads first, matched no pattern at all. The
# claim count never showed it, because the same numbers ARE checked in their
# square-bracket form elsewhere in the file; coverage is per-literal, and that
# is why it saw what a count could not.
#
# The separator class carries EN DASH and EM DASH as well as ASCII hyphen: the
# documents are typeset with U+2013, and `-` alone matched none of them.
PAREN_TRIPLE = re.compile(
    r"(?<![\d.])([+\-−]?\d\.\d{3,})\s*\(\s*(?:95%\s*CI,?\s*)?"
    r"([+\-−]?\d\.\d{3,})\s*[–—-]\s*([+\-−]?\d\.\d{3,})\s*\)")

# (point_key, lo_key, hi_key, label)
TRIPLES: list[tuple[str, str, str, str]] = [
    ("nsclc_isi", "nsclc_ci_lo", "nsclc_ci_hi", "NSCLC ISI"),
    ("pancancer_isi", "pancancer_ci_lo", "pancancer_ci_hi", "pan-TCGA ISI"),
    # These three were REPORTED AS UNMATCHED before their sources were traced:
    # two secondary-endpoint triples in 09-PAPER-DRAFT.md:293 and the rotation
    # null's family-mean interval in 14-SCIENCE-AUDIT.md:344.
    ("nsclc_dmae", "nsclc_dmae_lo", "nsclc_dmae_hi", "NSCLC secondary Delta-MAE"),
    ("pancancer_dmae", "pancancer_dmae_lo", "pancancer_dmae_hi",
     "pan-TCGA secondary Delta-MAE"),
    ("rotation_null_mean", "rotation_null_lo", "rotation_null_hi",
     "rotation-null family mean"),
    # A7 arms. Before these were listed, every scorer-sensitivity interval in
    # the manuscript was an UNMATCHED triple -- reported, never checked. The
    # control arm's triple is numerically identical to `nsclc_isi`, so it is
    # already covered; the other four were not.
    ("a7_meanz_unc", "a7_meanz_unc_lo", "a7_meanz_unc_hi",
     "A7 arm: mean-z, uncorrected"),
    ("a7_ssgsea_unc", "a7_ssgsea_unc_lo", "a7_ssgsea_unc_hi",
     "A7 arm: ssGSEA, uncorrected"),
    ("a7_recon_meanz", "a7_recon_meanz_lo", "a7_recon_meanz_hi",
     "A7 reconstruction: mean-z"),
    ("a7_recon_ssgsea_signed", "a7_recon_ssgsea_signed_lo",
     "a7_recon_ssgsea_signed_hi",
     "A7 reconstruction: ssGSEA (signed -- this interval is negative)"),
]

# FULL-PRECISION QUOTES, and why this rule rather than a per-value regex.
#
# The first attempt anchored each long decimal on its own correct prefix
# (`0\.31824\d{6,}`). That catches the realistic drift -- somebody retyping a
# 15-digit number and fumbling the tail -- but a value replaced wholesale slips
# through, because it no longer matches the pattern that would have checked it.
# Verified: substituting 0.399999999999999 for the NSCLC ISI was NOT caught.
#
# So the rule is positional instead of value-anchored: ANY decimal carrying 10+
# digits after the point is claiming to be a frozen artefact value, because
# nothing else in these documents is written to that precision. It must equal
# one of the authoritative values. Anything else is reported by value and line.
LONG_DECIMAL = re.compile(r"(?<![\d.])(\d\.\d{10,})(?![\d])")

# Claims with no summary.json to read, each carrying its derivation so it is
# auditable rather than asserted.
DERIVED: list[tuple[float, str, float, str, str | None]] = [
    (0.0218, r"(?:gap|difference)[^\d\n]{0,24}(0\.0\d{3})(?![\d])", 5e-5,
     "the A9 platform gap: the frozen macOS ISI minus HPC4's 0.2964", None),
]

DOCS = [
    "09-PAPER-DRAFT.md", "07-ABSTRACT-DRAFT.md", "14-SCIENCE-AUDIT.md", "15-CHECKLIST.md",
    "16-YOUR-TASKS.md", "pipeline/README.md",
    "poster/poster.html", "submission/SUBMISSION-CHECKLIST.md",
    "submission/COVER-LETTER.md", "pipeline/results/README.md",
    # CITATION.cff quotes the cohort size and cancer-type count in its abstract
    # and was not scanned by any earlier version.
    "CITATION.cff",
    # Added 2026-09-06 with S3's derived NSCLC multiplicity values. It states
    # results, so it belongs in the scan even though it is an internal document.
    "18-SUPPLEMENTARY-INVENTORY.md",
    # Added 2026-09-07. It is a LIVE working document, not a historical record,
    # and it carried a "Test suite | 78 tests | done" row that had been wrong for
    # roughly thirty tests without anything noticing -- because this file had
    # never been scanned.
    "11-SUBMISSION-PACK.md",
    # Added 2026-09-07 so that EXEMPT below is LOAD-BEARING. Both are skipped by
    # the exemption, so scanning them changes no claim count; what changes is
    # that the exemption now actually fires and says so in `skipped:`. Until
    # today neither file was in this list, which made every EXEMPT entry dead
    # code -- a guard for a door nobody walked through. See EXEMPT's note.
    "05-PRE-REGISTRATION.md",
    "02-PROJECT-DECISION.md",
]


def _num(s: str) -> float:
    return float(s.replace(",", "").replace("&nbsp;", "").replace("−", "-"))


def _tol_for(literal: str, floor: float) -> float:
    """Tolerance implied by how many decimals the document actually wrote.

    A document writing `0.386` is claiming 0.3855 <= x < 0.3865 and nothing
    more, so holding it to 5e-5 against a stored 0.38571423... reports drift
    that does not exist -- which is how the poster's correctly-rounded family
    mean r came up as a MISMATCH. Rounding half-width is 0.5 * 10^-decimals;
    the configured tolerance is a floor beneath it, never a ceiling.
    """
    lit = literal.replace(",", "").replace("&nbsp;", "")
    # SCIENTIFIC NOTATION needs the half-width SCALED by the exponent, or the
    # check goes vacuous. "1.669e-41" under the decimal rule below implies a
    # tolerance of 5e-8 against a value of 1.669e-41, which accepts every small
    # number ever written -- a pattern that cannot fail. Scaling gives 5e-45,
    # which is the interval the document actually claims. Added 2026-09-06,
    # when the first e-notation claims (S3's NSCLC p and q) were introduced.
    mant, sep, exp = lit.partition("e") if "e" in lit else lit.partition("E")
    if sep:
        ndec = len(mant.partition(".")[2])
        half = 0.5 * 10 ** -ndec if ndec else 0.5
        try:
            return half * 10 ** int(exp)
        except ValueError:
            pass
    frac = lit.partition(".")[2]
    ndec = len(frac)
    return max(floor, 0.5 * 10 ** -ndec + 1e-12) if ndec else max(floor, 0.5)


# An explicit, visible opt-out for lines that quote a number on purpose and
# must NOT be corrected. Three exist today, and each is a real category:
#
#   * a direct QUOTATION from a cited paper (Schmauch's 28 cancer types),
#   * a deliberate RECORD OF PRIOR DRIFT ("the previous 57 tests had drifted"),
#   * a RESOLVED, struck-through checklist item.
#
# A marker carries its reason in the document, where a reader sees it. The
# alternative -- inferring intent from strikethrough or from phrases like
# "previously" -- would silently swallow real drift that happened to be phrased
# that way, which is the failure mode this whole script exists to prevent.
IGNORE_MARKER = re.compile(r"<!--\s*numcheck:\s*ignore\b")


# Numeric literals that are never claims about a result, so listing them as
# "unchecked" in --coverage would bury the real gaps in noise. Each is matched
# against the WHOLE literal-with-context, and each is here for a reason:
# dates and version strings are identifiers; markdown list markers and heading
# levels are structure; a bare year is a citation.
COVERAGE_NOISE = re.compile(
    r"""(?x)
    ^20\d\d$ | ^19\d\d$                 # years / citation dates
  | ^\d{1,2}$                           # list markers, small counts, headings
  | ^0$ | ^1$
""")


def _numeric_literals(line: str) -> list[tuple[int, int, str]]:
    """Every number-shaped token in a line, as (start, end, literal)."""
    return [(m.start(), m.end(), m.group(0))
            for m in re.finditer(r"(?<![\w.])\d[\d,]*(?:\.\d+)?(?![\w])", line)]


def _compile_guard(guard: str | None
                   ) -> tuple[re.Pattern | None, bool, re.Pattern | None]:
    """Compile a guard, returning (regex, paragraph_scoped, line_exclusion).

    "PARA:<regex>" opts into paragraph scope; anything else is line-scoped,
    which is the historical default and stays byte-for-byte unchanged.

    A trailing "!!<regex>" adds a LINE-level exclusion: the line must NOT match
    it, whatever the scope of the positive part. This is what makes paragraph
    scope safe in text that discusses both cohorts at once. A paragraph
    comparing NSCLC with pan-cancer matches /NSCLC/ and /pan-cancer/ both, so
    on its own a paragraph guard would let the NSCLC pattern capture the
    pan-cancer number and report it as drift -- measured 2026-09-05, four such
    false mismatches across `15-CHECKLIST.md` and `poster/poster.html`. Writing
    "PARA:NSCLC!!pan-cancer" says what is actually meant: the paragraph is
    about NSCLC, and this particular line is not the one about the other
    cohort.
    """
    if guard is None:
        return None, False, None
    neg = None
    if "!!" in guard:
        guard, _, tail = guard.partition("!!")
        neg = re.compile(tail, re.IGNORECASE)
    if guard.startswith("PARA:"):
        return re.compile(guard[len("PARA:"):], re.IGNORECASE), True, neg
    if guard.startswith("WRAP:"):
        # WRAP implies paragraph scope AND paragraph MATCHING; the caller
        # handles the second half. An empty guard body means "no guard, just
        # let the pattern see across the line break".
        body = guard[len("WRAP:"):]
        return (re.compile(body, re.IGNORECASE) if body else None), True, neg
    return re.compile(guard, re.IGNORECASE), False, neg


def _paragraph_of(lines: list[str]) -> dict[int, str]:
    """Map each 1-based line number to the text of its blank-line paragraph.

    Markdown table rows and list items are contiguous non-blank lines, so they
    land in one paragraph together, which is the intended behaviour: a guard
    naming a cohort should reach every row of that cohort's table.
    """
    out: dict[int, str] = {}
    start = 1
    block: list[str] = []

    def flush(end: int) -> None:
        text = "\n".join(block)
        for k in range(start, end):
            out[k] = text

    for i, line in enumerate(lines, 1):
        if line.strip():
            block.append(line)
        else:
            flush(i)
            out[i] = ""
            block, start = [], i + 1
    flush(len(lines) + 1)
    return out


def _stylesheet_lines(lines: list[str]) -> set[int]:
    """1-based line numbers falling inside a <style> or <script> block.

    WHY, and why this is an honest exclusion rather than a flattering one.
    `poster/poster.html` reported 4% coverage (15 of 353 literals), which read
    as the worst-covered submission-facing document in the project. Measured
    2026-09-05: 175 of those 353 literals are stylesheet content -- `700` in
    `font-weight: 700`, `1.06` in `line-height: 1.06`, and the arithmetic in
    CSS comments like `/* 330 + 20 + 445 + 20 + 330 = 1145mm */`. None of them
    is a claim about a result, and none of them can be, because a stylesheet
    cannot state one. Counting them made the poster's coverage figure an
    artefact of the denominator rather than a measure of what is unverified.

    The exclusion is reported, not silent: `--coverage` prints the count it
    removed per document, so the smaller denominator can be audited rather
    than taken on trust. It touches ONLY the coverage denominator -- claim
    scanning is unchanged, so a pattern that matched inside a <style> block
    before still matches now, and a real claim hidden there would still be
    compared. The one exclusion worth naming: the A0 sheet size (1189 x 841 mm)
    appears only in the layout comment inside <style>, so it is no longer
    counted as an unchecked literal. It is a paper-size standard rather than a
    result, and nothing in this repository derives from it.
    """
    inside: set[int] = set()
    depth = 0
    for i, line in enumerate(lines, 1):
        low = line.lower()
        opens = len(re.findall(r"<(?:style|script)\b", low))
        closes = len(re.findall(r"</(?:style|script)>", low))
        if depth or opens:
            inside.add(i)
        depth = max(0, depth + opens - closes)
    return inside


def scan(root: Path, table: dict[str, float], every_run: dict[str, float],
         *, verbose: bool = False, quiet: bool = False,
         coverage: dict[str, list] | None = None,
         match_counts: dict[str, int] | None = None,
         match_files: dict[str, set[str]] | None = None
         ) -> tuple[int, int, list[str], list[str], list[str]]:
    """Scan `root`'s documents.

    Returns (checked, failures, unmatched, skipped_files, skipped_lines).

    If `coverage` is given it is filled in with, per document, which numeric
    literals were actually compared against an authoritative value and which
    were merely seen. See `--coverage` and the note on why a claim COUNT is not
    a coverage measure.

    If `match_files` is given it records, per pattern description, the SET of
    documents that pattern actually compared something in -- its REACH. This is
    recorded inside `compare()` rather than at the pattern loops, so it cannot
    drift away from what was really checked: every comparison in every claim
    family passes through that one function.
    """
    checked = failures = 0
    unmatched: list[str] = []
    skipped_files: list[str] = []
    skipped_lines: list[str] = []

    def say(*a):
        if not quiet:
            print(*a)

    for rel in DOCS:
        path = root / rel
        if not path.exists():
            skipped_files.append(f"{rel} (not found)")
            continue
        if path.name in EXEMPT:
            skipped_files.append(f"{rel} (exempt: {EXEMPT[path.name]})")
            continue
        lines = path.read_text(errors="replace").splitlines()
        paras = _paragraph_of(lines)

        def guarded(grx, para, neg, i, line):
            """Does the guard admit line `i`? Line scope unless PARA-scoped."""
            if neg is not None and neg.search(line):
                return False
            if grx is None:
                return True
            return bool(grx.search(paras[i] if para else line))
        ignored = {i for i, line in enumerate(lines, 1)
                   if IGNORE_MARKER.search(line)}
        if ignored:
            skipped_lines.append(f"{rel}: {len(ignored)} line(s) marked "
                                 "numcheck:ignore")

        # Which (line, literal) pairs were actually compared. Recorded even when
        # `coverage` is None so the bookkeeping cannot drift between modes.
        covered: set[tuple[int, str]] = set()

        def compare(expected, literal, i, desc, floor):
            nonlocal checked, failures
            got = _num(literal)
            tol = _tol_for(literal, floor)
            checked += 1
            covered.add((i, literal.strip()))
            if match_files is not None:
                match_files.setdefault(desc, set()).add(rel)
            if abs(got - expected) > tol:
                failures += 1
                say(f"  MISMATCH {rel}:{i}")
                say(f"    {desc}: document says {got!r}, expected {expected!r}")
            elif verbose:
                say(f"  ok       {rel}:{i}  {desc} = {got}")

        for key, pattern, tol, desc, guard in CLAIMS:
            if key not in table:
                continue
            rx = re.compile(pattern, re.IGNORECASE)
            grx, para, neg = _compile_guard(guard)
            if match_counts is not None:
                match_counts.setdefault(desc, 0)
            if isinstance(guard, str) and guard.startswith("WRAP:"):
                # WRAP: the number may sit on a DIFFERENT LINE from its anchor.
                # Matching per line cannot see those claims at all, so they
                # pass while checking nothing -- `ancestry_perm_p95` even
                # carried a `\n?` in its pattern, written for a scan that
                # never had a newline to offer it. Here the paragraph is
                # whitespace-normalised and matched ONCE, so a claim that
                # wraps is still a claim. Anchor the pattern on enough words
                # to stay unambiguous: a paragraph is a wide net.
                seen_para: set[int] = set()
                for i, line in enumerate(lines, 1):
                    if i in ignored:
                        continue
                    text = paras.get(i)
                    if not text or id(text) in seen_para:
                        continue
                    seen_para.add(id(text))
                    flat = " ".join(text.split())
                    if neg is not None and neg.search(flat):
                        continue
                    if grx is not None and not grx.search(flat):
                        continue
                    for m in rx.finditer(flat):
                        compare(table[key], m.group(1), i, desc, tol)
                        if match_counts is not None:
                            match_counts[desc] += 1
                continue
            for i, line in enumerate(lines, 1):
                if i in ignored or not guarded(grx, para, neg, i, line):
                    continue
                for m in rx.finditer(line):
                    compare(table[key], m.group(1), i, desc, tol)
                    if match_counts is not None:
                        match_counts[desc] += 1

        for keys, pattern, tol, desc, guard in SET_CLAIMS:
            live = [k for k in keys if k in table]
            if not live:
                continue
            rx = re.compile(pattern, re.IGNORECASE)
            grx, para, neg = _compile_guard(guard)
            for i, line in enumerate(lines, 1):
                if i in ignored or not guarded(grx, para, neg, i, line):
                    continue
                for m in rx.finditer(line):
                    got = _num(m.group(1))
                    eff = _tol_for(m.group(1), tol)
                    checked += 1
                    covered.add((i, m.group(1).strip()))
                    if any(abs(table[k] - got) <= eff for k in live):
                        if verbose:
                            say(f"  ok       {rel}:{i}  {desc} = {got}")
                    else:
                        failures += 1
                        say(f"  MISMATCH {rel}:{i}")
                        say(f"    {desc}: document says {got!r}, which matches none of "
                            + ", ".join(f"{k}={table[k]!r}" for k in live))

        for guard, pairs, desc in FRACTIONS:
            live = [(n, d) for n, d in pairs if n in table and d in table]
            if not live:
                continue
            # FRACTIONS guards stay LINE-scoped and case-SENSITIVE: a fraction
            # like "10/32" is only interpretable from the text beside it, and
            # the character class is load-bearing (see the note above). Assert
            # rather than silently treating "PARA:" as literal text, so a
            # future paragraph guard added here fails loudly instead of
            # matching nothing.
            assert not guard.startswith("PARA:"), (
                f"FRACTIONS guard {guard!r} is line-scoped only")
            grx = re.compile(guard)
            for i, line in enumerate(lines, 1):
                if i in ignored or not grx.search(line):
                    continue
                for m in re.finditer(r"(?<![\d./])(\d{1,3})\s*/\s*(\d{1,3})(?![\d/])",
                                     line):
                    got = (int(m.group(1)), int(m.group(2)))
                    checked += 1
                    covered.add((i, m.group(1)))
                    covered.add((i, m.group(2)))
                    if any(got == (int(table[n]), int(table[d])) for n, d in live):
                        if verbose:
                            say(f"  ok       {rel}:{i}  {desc} = {got[0]}/{got[1]}")
                    else:
                        failures += 1
                        say(f"  MISMATCH {rel}:{i}")
                        say(f"    {desc}: document says {got[0]}/{got[1]}, which matches "
                            "none of " + ", ".join(
                                f"{int(table[n])}/{int(table[d])}" for n, d in live))

        for expected, pattern, tol, desc, guard in DERIVED:
            rx = re.compile(pattern, re.IGNORECASE)
            grx, para, neg = _compile_guard(guard)
            for i, line in enumerate(lines, 1):
                if i in ignored or not guarded(grx, para, neg, i, line):
                    continue
                for m in rx.finditer(line):
                    compare(expected, m.group(1), i, desc, tol)

        # Full-precision quotes: must equal SOME authoritative value exactly.
        for i, line in enumerate(lines, 1):
            if i in ignored:
                continue
            for m in LONG_DECIMAL.finditer(line):
                checked += 1
                covered.add((i, m.group(1)))
                got = float(m.group(1))
                hit = next((k for k, v in table.items()
                            if abs(v - got) <= 1e-15), None)
                if hit is None:
                    # A full-precision number may legitimately quote a
                    # SUPERSEDED run -- results/README.md does exactly that,
                    # to say what each stored run produced. That is provenance,
                    # not a headline claim, so it is allowed as long as some
                    # stored run actually produced it.
                    hit = next((k for k, v in every_run.items()
                                if abs(v - got) <= 1e-15), None)
                    if hit is not None:
                        hit += " (superseded run, provenance)"
                if hit is None:
                    failures += 1
                    say(f"  MISMATCH {rel}:{i}")
                    say(f"    {got!r} is written to full precision but "
                        "matches NO frozen value.")
                    near = min(table.items(), key=lambda kv: abs(kv[1] - got))
                    say(f"    nearest is {near[0]} = {near[1]!r} "
                        f"(off by {abs(near[1] - got):.3e})")
                elif verbose:
                    say(f"  ok       {rel}:{i}  full precision = {hit}")

        # Bracketed triples.
        for i, line in enumerate(lines, 1):
            if i in ignored:
                continue
            for m in list(TRIPLE.finditer(line)) + list(
                    PAREN_TRIPLE.finditer(line)):
                # SIGNED, not magnitudes. This used to `lstrip("+-")` and
                # compare `abs()`, which had two consequences: a triple written
                # with U+2212 crashed `float()` once the regex was widened to
                # match it, and -- worse -- a sign-flipped interval compared
                # equal to its own negation. For A7's ssGSEA reconstruction the
                # sign IS the result, so magnitude-only checking would have
                # accepted the one error most worth catching.
                pt_lit, lo_lit, hi_lit = m.groups()
                pt = _num(pt_lit)
                for pk, lk, hk, label in TRIPLES:
                    if pk not in table:
                        continue
                    if abs(pt - table[pk]) <= _tol_for(pt_lit, 5e-5):
                        # The point estimate is what IDENTIFIED the quantity,
                        # so it is verified even though `compare` is only
                        # called on the two bounds.
                        checked += 1
                        covered.add((i, pt_lit.strip()))
                        compare(table[lk], lo_lit, i,
                                f"{label} CI low, in the triple beside "
                                f"{m.group(1)}", 5e-5)
                        compare(table[hk], hi_lit, i,
                                f"{label} CI high, in the triple beside "
                                f"{m.group(1)}", 5e-5)
                        break
                else:
                    unmatched.append(f"{rel}:{i}  {m.group(0)}")

        # Number-WORDS. Same shape as CLAIMS, but the capture is a word and is
        # translated through NUMBER_WORDS before comparing. A word the table
        # does not know is ignored rather than failed -- the patterns are loose
        # by necessity (English), so an unknown capture means "not a number",
        # not "wrong number".
        for key, pattern, desc, guard in WORD_CLAIMS:
            if key not in table:
                continue
            rx = re.compile(pattern, re.IGNORECASE)
            grx, para, neg = _compile_guard(guard)
            for i, line in enumerate(lines, 1):
                if i in ignored or not guarded(grx, para, neg, i, line):
                    continue
                for m in rx.finditer(line):
                    word = m.group(1).lower()
                    if word not in NUMBER_WORDS:
                        continue
                    got = NUMBER_WORDS[word]
                    checked += 1
                    covered.add((i, m.group(1)))
                    if got != int(table[key]):
                        failures += 1
                        say(f"  MISMATCH {rel}:{i}")
                        say(f"    {desc}: document says {word!r} "
                            f"({got}), expected {int(table[key])}")
                    elif verbose:
                        say(f"  ok       {rel}:{i}  {desc} = {word}")

        if coverage is not None:
            style = _stylesheet_lines(lines)
            seen, in_style, unchecked = 0, 0, []
            for i, line in enumerate(lines, 1):
                if i in ignored:
                    continue
                lits = _numeric_literals(line)
                if i in style:
                    # Counted and reported separately, never silently dropped:
                    # a shrinking denominator that nobody can see is how a
                    # coverage figure stops meaning anything.
                    in_style += len(lits)
                    continue
                for _s, _e, lit in lits:
                    seen += 1
                    if (i, lit) in covered:
                        continue
                    if COVERAGE_NOISE.match(lit.replace(",", "")):
                        continue
                    unchecked.append((i, lit, line.strip()[:96]))
            # Count only claims compared OUTSIDE the excluded regions, or a
            # pattern that happened to match inside a <style> block would
            # inflate the numerator against a denominator that excludes it.
            n_cov = len([1 for (i, _lit) in covered if i not in style])
            coverage[rel] = [seen, n_cov, unchecked, in_style]

    return checked, failures, unmatched, skipped_files, skipped_lines


# Each entry is (description, substitution) where substitution is applied to
# every document. The self-test asserts each one produces at least one MISMATCH.
# If a pattern ever goes vacuous -- as v1's did, wholesale -- its row here goes
# red instead of the whole script passing by construction.
SELF_TEST: list[tuple[str, str, str]] = [
    ("NSCLC ISI, full precision, replaced wholesale",
     r"0\.318240180054566", "0.399999999999999"),
    # A9 split determinism. The injections target the PAPER's phrasing, which
    # wraps between "216 of" and "944 patients" -- the pattern uses \s+ for
    # exactly that reason, and this row is what proves the wrap is handled.
    ("A9 tied-site count drifted",
     r"51 of the 68 NSCLC sites", "37 of the 68 NSCLC sites"),
    ("A9 moved-patient count drifted",
     r"216 of\s+944 patients change fold", "409 of 944 patients change fold"),
    # A10, the ssGSEA tie convention. Unlike the A9 rows above, every anchor
    # here sits on ONE line with its number -- the paragraph was reworded so it
    # does rather than relying on `\n?` in the pattern. Three wrap-spanning
    # anchors have silently matched nothing in this file's history; these four
    # rows are what prove these do not.
    ("A10 tied rank-entry count drifted",
     r"35,394,448 of 38,747,424", "35,394,449 of 38,747,424"),
    ("A10 differing-cell count drifted",
     r"different scores for 15,060 of 15,104",
     "different scores for 15,061 of 15,104"),
    ("A10 max ssGSEA score difference drifted",
     r"difference between the two orderings is 0\.0128",
     "difference between the two orderings is 0.0193"),
    ("A10 ssGSEA score sd drifted",
     r"per-signature score standard deviation of 0\.0211",
     "per-signature score standard deviation of 0.0299"),
    ("poster ssGSEA uncorrected arm drifted",
     r"uncorrected contrast \(\+0\.180", "uncorrected contrast (+0.311"),
    ("pan-cancer outcome C-index effect drifted",
     r"median C-index excess of 0\.0257", "median C-index excess of 0.0413"),
    ("a CI bound drifted away from its point estimate",
     r"\[\s*\+?0\.2614\s*,", "[+0.2999,"),
    # Injected against the phrasing the documents ACTUALLY use. The first
    # version of this row targeted "partition-choice sd of 0.0089", which
    # appears nowhere -- the self-test reported NO SUCH TEXT rather than
    # passing, which is the whole point of counting substitutions.
    ("partition sd", r"[Pp]artition sd(\*{0,2})\s*(\|?\s*\*{0,2})0\.0089",
     r"Partition sd\1\g<2>0.0111"),
    ("bootstrap SE", r"[Bb]ootstrap SE(\*{0,2})(\s*\|?\s*(?:of )?\*{0,2})0\.0294",
     r"Bootstrap SE\1\g<2>0.0333"),
    ("interval-understated percent", r"4\.5%\s*\*{0,2}\s*too narrow",
     "9.9% too narrow"),
    # A5 pan-cancer, landed 2026-09-05. Each row injects into the phrasing the
    # documents actually use, and each is worth a row for a separate reason:
    # the sd and SE prove the paragraph-scoped guards reach the pan-cancer
    # paragraph at all, and the percent row proves the `(\d+\.\d+)` capture
    # widened this session actually sees a TWO-DIGIT percent -- under the old
    # `(\d\.\d+)` this row would have reported MISSED, because the pattern
    # could not match "16.0" in the first place.
    ("pan-cancer partition sd",
     r"[Pp]artition sd(\*{0,2})\s*(\|?\s*(?:is |of )?\*{0,2})0\.0065",
     r"Partition sd\1\g<2>0.0099"),
    ("pan-cancer bootstrap SE",
     r"bootstrap (SE|standard error)(\s*\|?\s*(?:of )?\*{0,2})0\.0111",
     r"bootstrap \1\g<2>0.0222"),
    ("pan-cancer interval-understated percent (two digits)",
     r"16\.0%\s*\*{0,2}\s*too narrow", "26.0% too narrow"),
    ("pan-cancer partition share of total variance",
     r"partition choice accounts for 25\.7%",
     "partition choice accounts for 35.7%"),
    ("NSCLC partition share of total variance",
     r"partition choice accounts for 8\.4%",
     "partition choice accounts for 9.4%"),
    # S3's derived NSCLC multiplicity values. Both are in SCIENTIFIC NOTATION,
    # which is why they are here: under the pre-2026-09-06 `_tol_for` a claim of
    # 1.669e-41 carried a 5e-8 ABSOLUTE tolerance, so every wrong value these
    # rows could inject would have been accepted and both rows would report
    # MISSED. They are the regression guard for the exponent-scaling rule.
    ("frozen NSCLC smallest per-signature p (scientific notation)",
     r"min p = 1\.669e-41", "min p = 2.669e-41"),
    ("frozen NSCLC largest BH q (scientific notation)",
     r"max q = 5\.971e-06", "max q = 6.971e-06"),
    # The A9/A10 sensitivity run. The shift row matters most: it is DERIVED from
    # the two summaries, so this proves the derived key is compared and not just
    # computed. Both are the numbers a reader would use to judge whether the
    # correction changed the finding, which is exactly the kind of number this
    # project has previously let go stale.
    ("stable-sort sensitivity ISI", r"ISI of 0\.2922", "ISI of 0.3922"),
    ("stable-sort sensitivity shift", r"moves by 0\.0260", "moves by 0.0460"),
    ("signatures beating null (16/16 -> 15/16)", r"\b16/16\b", "15/16"),
    ("outcome arm (10/32 -> 11/32)", r"\b10/32\b", "11/32"),
    ("pan-cancer cohort size", r"7,168 TCGA patients", "7,169 TCGA patients"),
    ("number of cancer types", r"31 cancer types", "32 cancer types"),
    ("site AUROC", r"median AUROC 0\.998", "median AUROC 0.888"),
    ("minimum detectable effect", r"minimum detectable effect is 0\.046",
     "minimum detectable effect is 0.055"),
    ("rotation-null p, replaced wholesale", r"p = 0\.000999", "p = 0.5"),
    ("M_eff", r"M_eff = 5\.00", "M_eff = 7.00"),
    ("secondary Delta-MAE triple", r"0\.0088 \[0\.0078, 0\.0098\]",
     "0.0088 [0.0071, 0.0098]"),
    ("rotation-null family-mean triple", r"0\.1097 \[0\.0977, 0\.1223\]",
     "0.1097 [0.0911, 0.1223]"),
    # A7, one row per family the scorer result introduced. Each is anchored on
    # a FROZEN artefact value, per the note below.
    ("A7 ssGSEA uncorrected arm triple",
     r"\+0\.180 \[0\.139, 0\.216\]", "+0.180 [0.139, 0.244]"),
    ("A7 mean-z uncorrected arm triple",
     r"\+0\.319 \[0\.266, 0\.373\]", "+0.319 [0.266, 0.399]"),
    ("A7 ssGSEA reconstruction triple (negative value)",
     r"−0\.1333 \[−0\.1831, −0\.0835\]", "−0.1333 [−0.1911, −0.0835]"),
    ("A7 mean-z reconstruction triple",
     r"\+0\.3177 \[0\.2822,\s*\n?0\.3533\]", "+0.3177 [0.2999, 0.3533]"),
    ("A7 ssGSEA null split-half reliability",
     r"0\.799 under mean-z to 0\.258", "0.799 under mean-z to 0.318"),
    ("A7 disattenuation divisor",
     r"1/√0\.258", "1/√0.358"),
    ("A7 full-precision reliability (audit table)",
     r"0\.9619780095308824", "0.9619780095309999"),
    ("A7 full-precision arm value (audit table)",
     r"0\.1799529379910222", "0.1799529379919999"),
    # A10's pinned re-run. Every one of these substitutions has to be caught by
    # a DIFFERENT pattern from the pre-fix rows above, because the whole point
    # of the pair is that the two numbers can move independently. A single
    # authority covering both would pass while checking one of them.
    ("A10 pinned ssGSEA value (Results)",
     r"contrast:\s*\n?\+0\.174", "contrast:\n+0.199"),
    ("A10 pinned ssGSEA value (Limitations 8)",
     r"moves from \+0\.180 to\s*\n?\+0\.174", "moves from +0.180 to\n+0.211"),
    ("A10 pinned ssGSEA value (poster panel 6)",
     r"\+0\.180; \+0\.174 pinned", "+0.180; +0.161 pinned"),
    ("A10 pinned ssGSEA value (Limitations 9)",
     r"and at \+0\.174 under the pinned", "and at +0.188 under the pinned"),
    # A11's corrected NSCLC values. The substitution for the median R^2 is the
    # PAN-CANCER value on purpose: the failure being pinned is not "a digit
    # changed", it is "this cell was filled from the other cohort", and a row
    # that injects a random number would not prove the checker can tell those
    # two apart.
    ("A11 NSCLC Control C median R^2 replaced by pan-cancer's",
     r"median R² 0\.14533984970354008", "median R² 0.44379065488261865"),
    ("A11 NSCLC science-gaps patient count",
     r"944 of 7,168 patients, 31 sites", "7,168 of 7,168 patients, 31 sites"),
    ("A10 shift between the pre-fix and pinned arms",
     r"a shift of 0\.006, or", "a shift of 0.031, or"),
    ("A10 shift as a percent of the pre-fix half-width",
     r"or 16% of the pre-fix interval", "or 42% of the pre-fix interval"),
    ("A10 pinned arm as a percent of the mean-z arm",
     r"magnitude \(54% against", "magnitude (61% against"),
    ("A7 pre-fix arm as a percent of the mean-z arm",
     r"against 56%\)", "against 71%)"),
    ("A10 pinned signatures above null (audit table)",
     r"15 of 16 pinned", "14 of 16 pinned"),
    ("A7 pre-fix signatures above null (audit table)",
     r"15 of 16 pre-fix", "13 of 16 pre-fix"),
    ("A10 full-precision pinned arm value (audit table)",
     r"0\.17377773241035724", "0.17377773241039999"),
    ("A10 full-precision shift (audit table)",
     r"0\.00617520558066495", "0.00617520558069999"),
    # NOTE: the two LIVE-MEASURED claims -- the pytest collection count and the
    # abstract character count -- are appended by `_self_test_rows()` below,
    # built from the authority table rather than written as literals here. A
    # literal goes stale the moment the value it names legitimately changes:
    # this row was `\b78 tests\b`, four tests were added, the documents were
    # updated to 82, and the injection then matched nothing and reported NO
    # SUCH TEXT. Everything above is anchored on a FROZEN artefact value, which
    # by definition does not move.
    #
    # ---- added 2026-09-04, one per newly covered claim family --------------
    #
    # Each injection moves the number to a value the NEW pattern must still
    # match, so the pattern fires and the COMPARISON fails. Anchoring on the
    # correct value's digits in a way that makes the wrong value unmatchable is
    # the vacuous-pattern defect this file has now shipped four times; these
    # rows exist to prove these ten patterns do not have it.
    ("median image-signature correlation, pan-TCGA",
     r"correlation was 0\.630 pan-cancer", "correlation was 0.930 pan-cancer"),
    ("median image-signature correlation, NSCLC",
     r"pan-cancer and 0\.400 in NSCLC", "pan-cancer and 0.700 in NSCLC"),
    ("median split-robustness Delta r, pan-TCGA",
     r"median Δr = \+0\.020", "median Δr = +0.070"),
    ("median split-robustness Delta r, NSCLC",
     r"(median Δr = \+0\.0\d\d and )\+0\.023", r"\g<1>+0.077"),
    ("ancestry arm cohort size",
     r"ancestry arm \(6,580 patients", "ancestry arm (6,999 patients"),
    ("ancestry x site Cramer's V",
     r"Cramér's V = 0\.537", "Cramér's V = 0.937"),
    ("ancestry group counts (EUR)", r"EUR 5,312", "EUR 5,999"),
    ("small-panel bracket (k=10 median gap)",
     r"\*\*0\.453 at k=10\*\*", "**0.953 at k=10**"),
    ("ssGSEA observed split-half reliability",
     r"split-half 0\.934", "split-half 0.834"),
    ("mean-z null reliability quoted at 3dp",
     r"from 0\.799 under mean-z", "from 0.699 under mean-z"),

    # ---- added 2026-09-05 ---------------------------------------------------
    #
    # Three of these cover patterns that ALREADY EXISTED and had no self-test
    # row, which is its own kind of vacuous: an untested pattern is a pattern
    # nobody has seen go red. Auditing the CLAIMS table against SELF_TEST found
    # `control_c_calibrated_r2`, `rotation_observed` and the three A7 agreement
    # statistics in that state.
    ("Control C calibrated R^2 (pre-existing pattern, was DEAD)",
     r"\*\*Calibrated\*\* \| \*\*0\.0410\*\*",
     "**Calibrated** | **0.0610**"),
    ("observed family mean r (pre-existing pattern, never self-tested)",
     r"family mean r(\*{0,2})[= ]{1,3}(\*{0,2})0\.386",
     r"family mean r\1=\g<2>0.986"),
    ("A7 cross-scorer agreement, Pearson (pre-existing, never self-tested)",
     r"Pearson r = 0\.633", "Pearson r = 0.733"),
    ("A7 cross-scorer agreement, Spearman (pre-existing, never self-tested)",
     r"Spearman ρ = 0\.582", "Spearman ρ = 0.682"),
    ("A7 residualised agreement, Pearson (pre-existing, never self-tested)",
     r"agree closely \(r = 0\.766", "agree closely (r = 0.866"),
    ("A7 residualised agreement, mean |diff| (pre-existing, never self-tested)",
     r"mean \|difference\| 0\.049", "mean |difference| 0.149"),
    # The reliability section.
    ("median raw alpha of random sets", r"median α 0\.974", "median α 0.874"),
    ("median residualised alpha of random sets",
     r"null's α to 0\.801", "null's α to 0.901"),
    ("median raw alpha of curated signatures",
     r"α fell from 0\.984", "α fell from 0.884"),
    ("median residualised alpha of curated signatures",
     r"α fell from 0\.984 to 0\.976", "α fell from 0.984 to 0.876"),
    ("mean raw reliability gap, pan-TCGA",
     r"understated: \+0\.014", "understated: +0.114"),
    ("mean raw reliability gap, NSCLC",
     r"pan-cancer and \+0\.032 NSCLC", "pan-cancer and +0.132 NSCLC"),
    ("mean residualised reliability gap, pan-TCGA",
     r"\*\*\+0\.184 and", "**+0.284 and"),
    ("mean residualised reliability gap, NSCLC",
     r"and \+0\.167\*\*", "and +0.267**"),
    ("reliability gap fold difference, pan-TCGA",
     r"a 13-fold and", "a 17-fold and"),
    ("reliability gap fold difference, NSCLC",
     r"-fold and 5\.3-fold", "-fold and 8.3-fold"),
    ("angiogenesis reliability gap",
     r"reaches \+0\.420 for the smallest", "reaches +0.520 for the smallest"),
    # Residualised correlations and null means.
    ("median residualised correlation, pan-TCGA",
     r"correlations were 0\.394", "correlations were 0.494"),
    ("median residualised correlation, NSCLC",
     r"and 0\.372 against null means", "and 0.472 against null means"),
    ("median null-set correlation, pan-TCGA",
     r"null means of 0\.106", "null means of 0.206"),
    ("median null-set correlation, NSCLC",
     r"null means of 0\.106 and 0\.051", "null means of 0.106 and 0.151"),
    ("median image-axis correlation, pan-TCGA",
     r"image–axis correlation was 0\.067", "image–axis correlation was 0.167"),
    ("median image-axis correlation, NSCLC",
     r"pan-cancer and \*\*−0\.023\*\* in NSCLC",
     "pan-cancer and **−0.123** in NSCLC"),
    # Purity decomposition and covariate baseline.
    ("purity-partial correlation, pan-TCGA",
     r"was 0\.547 pan-cancer \(against", "was 0.647 pan-cancer (against"),
    ("unadjusted correlation, pan-TCGA",
     r"\(against 0\.627 unadjusted", "(against 0.727 unadjusted"),
    ("incremental adjusted R^2, pan-TCGA",
     r"stage was 0\.036 and", "stage was 0.136 and"),
    ("incremental adjusted R^2, NSCLC",
     r"stage was 0\.036 and 0\.054", "stage was 0.036 and 0.154"),
    ("covariate-only baseline correlation, pan-TCGA",
     r"reached median r = 0\.655", "reached median r = 0.755"),
    ("embedding gain over covariates, pan-TCGA",
     r"\*\*−0\.000\*\* over it", "**−0.100** over it"),
    ("embedding gain over covariates, NSCLC",
     r"embedding added \+0\.031", "embedding added +0.131"),
    # Label-side site variance, all four quoted as percentages or 3dp.
    ("site-only R^2 of the label (percent)",
     r"median 44% of the label", "median 54% of the label"),
    ("site-given-type R^2 of the label (percent)",
     r"reduces it to 5\.5%", "reduces it to 8.5%"),
    ("Control C calibrated, quoted as a percent",
     r"gives 4\.1%", "gives 9.1%"),
    ("plate-within-site R^2 of the label",
     r"explains 0\.000", "explains 0.500"),
    # The immune-versus-other contrast, all three places it is quoted.
    ("immune-vs-other Mann-Whitney p (Results)",
     r"Mann–Whitney p = 0\.031", "Mann–Whitney p = 0.046"),
    ("immune-vs-other Fisher p (Results)",
     r"Fisher p = 0\.093", "Fisher p = 0.193"),
    ("immune-vs-other p (abstract)", r"\*P\* = 0\.031", "*P* = 0.046"),
    # The abstract's remaining headroom. Both operands were already checked and
    # the SUBTRACTION between them was not, so the document read
    # "2,515 / 2,600 with 83 characters of headroom" -- each count right, the
    # difference wrong, on the sentence that says how much room is left before
    # a hard deadline. Injecting 83 is therefore the exact historical error.
    ("abstract headroom (limit minus live count)",
     r"85 characters of headroom", "83 characters of headroom"),
    ("immune-vs-other p, IL6 counterfactual (Figure 3 caption)",
     r"\(p = 0\.031 → 0\.016\)", "(p = 0.031 → 0.048)"),
    # Number-WORDS. The exact defect: the reference list dropped to twelve and
    # the prose still said thirteen. Every other pattern in this file matches
    # digits and could not see it.
    ("reference count stated as a word",
     r"\btwelve references, all cited", "thirteen references, all cited"),
    ("reference count stated as a word (verification sentence)",
     r"for all twelve re-verified", "for all thirteen re-verified"),
    # The parenthetical CI form, which `--coverage` showed was unparsed. One
    # row per cohort per bound, injected into the ABSTRACT's phrasing, because
    # that is the sentence that had no checker at all.
    ("parenthetical CI, pan-TCGA low bound",
     r"\(95% CI,? 0\.269", "(95% CI, 0.229"),
    ("parenthetical CI, pan-TCGA high bound",
     r"(\(95% CI,? 0\.269[–-])0\.313", r"\g<1>0.353"),
    ("parenthetical CI, NSCLC low bound",
     r"\(0\.261([–-]0\.376\))", r"(0.201\1"),
    ("parenthetical CI, NSCLC high bound",
     r"(\(0\.261[–-])0\.376\)", r"\g<1>0.316)"),
    # Per-signature excess ranges, quoted in the Results and in the provenance
    # note at the top of the manuscript.
    ("per-signature excess range, pan-TCGA low",
     r"excess 0\.165", "excess 0.135"),
    ("per-signature excess range, NSCLC high",
     r"0\.148[–-]0\.440", "0.148–0.470"),
    # 07-ABSTRACT-DRAFT.md's own Results sentence.
    ("immune-vs-other p (abstract draft, plain form)",
     r"others, p=0\.031", "others, p=0.046"),
    ("reliability gap fold difference (abstract draft)",
     r"understated it 13-fold", "understated it 19-fold"),
    ("median raw alpha of random sets (abstract draft, 2dp)",
     r"alpha from 0\.97 to", "alpha from 0.87 to"),
    ("median residualised alpha of random sets (abstract draft, 2dp)",
     r"alpha from 0\.97 to 0\.80", "alpha from 0.97 to 0.90"),
    ("small-panel bracket quoted at 2dp (abstract draft)",
     r"0\.08 at 160 genes", "0.28 at 160 genes"),

    # ---- rows for the patterns added 2026-09-06 from `--coverage` ----
    # Each injection targets the phrasing of the document the NEW pattern was
    # written to reach, not a phrasing some other pattern already covered. A row
    # that a pre-existing pattern would catch anyway proves nothing about the
    # new one.
    #
    # The abstract's MDE. Deliberately anchored on the abstract-only tail
    # "; the real" -- 09-PAPER-DRAFT.md:457 writes "0.046 in C-index, while"
    # and is already reached by the older "minimum detectable effect" row, so
    # substituting there too would let this row pass on the paper's back.
    ("MDE in the abstract, whose line wraps before \"detectable\"",
     r"0\.046 in C-index; the real", "0.077 in C-index; the real"),
    ("pan-cancer cohort size, \"N patients across M cancer types\" form",
     r"7,168 patients across", "7,169 patients across"),
    ("NSCLC cohort size, abstract Methods form (wraps after \"NSCLC\")",
     r"\(n=944\) as a second cohort", "(n=955) as a second cohort"),
    ("pan-cancer image-to-global-axis correlation (abstract phrasing)",
     r"image-axis r = 0\.07", "image-axis r = 0.77"),
    ("small-panel k=10 gap quoted at 2dp (abstract draft)",
     r"to 0\.45 at 10", "to 0.95 at 10"),
    ("small-panel Spearman in the abstract, which writes \"at 10\" not \"k=10\"",
     r"\(Spearman -1\.00\)", "(Spearman -0.10)"),
    # Proves the alpha patterns really did widen to the paper's "α", rather
    # than continuing to pass on the abstract's "alpha" alone.
    ("Cronbach's alpha written as the GREEK letter (paper, not abstract)",
     r"Cronbach's α from 0\.97", "Cronbach's α from 0.57"),
]


def _self_test_rows(table) -> list[tuple[str, str, str]]:
    """SELF_TEST plus the rows that must track the live-measured values."""
    rows = list(SELF_TEST)
    if "test_count" in table:
        n = int(table["test_count"])
        rows.append(("test count (live pytest collection)",
                     rf"\b{n} tests\b", f"{n + 1} tests"))
    if "pancancer_n_sites" in table:
        n = int(table["pancancer_n_sites"])
        rows.append(("tissue source sites (live recount from meta.parquet)",
                     rf"\b{n} tissue source sites",
                     f"{n + 1} tissue source sites"))
    if "abstract_chars" in table and "abstract_limit" in table:
        n, lim = int(table["abstract_chars"]), int(table["abstract_limit"])
        rows.append(("abstract character count (live recount)",
                     rf"{n // 1000},?{n % 1000:03d}\s*/\s*{lim // 1000},?{lim % 1000:03d}",
                     f"{n - 100:,} / {lim:,}"))
    # The snapshot pair, added 2026-09-07. Both are live-measured, so like the
    # rows above they cannot be written as constants. Note the file count and
    # the byte figure are separated by a LINE BREAK in
    # `submission/SUBMISSION-CHECKLIST.md` and not in `09-PAPER-DRAFT.md`;
    # injection runs against whole-file text, so `\s+` covers both without the
    # row having to know which document it is rewriting.
    if "snapshot_files" in table:
        n = int(table["snapshot_files"])
        rows.append(("public snapshot file count (live rebuild)",
                     rf"rebuild stages \*\*{n} files",
                     f"rebuild stages **{n + 1} files"))
    return rows


def self_test(table, every_run) -> int:
    """Prove each family of claim can still FAIL. See the module docstring."""
    rows = _self_test_rows(table)
    print("SELF-TEST: injecting a wrong value for each family of claim.")
    print("Every row must report CAUGHT. A row reporting MISSED means that "
          "pattern has\ngone vacuous and is passing by construction.\n")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "repo"
        for rel in DOCS:
            src = REPO / rel
            if src.exists():
                dst = tmp / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)

        base_checked, base_fail, _, _, _ = scan(tmp, table, every_run, quiet=True)
        print(f"  baseline (unmodified copy): {base_checked} claims, "
              f"{base_fail} mismatch(es)")
        if base_fail:
            print("  !! the unmodified copy already fails; fix that first")

        originals = {rel: (tmp / rel).read_text(errors="replace")
                     for rel in DOCS if (tmp / rel).exists()}
        missed = []
        for desc, pattern, replacement in rows:
            rx = re.compile(pattern)
            hits = 0
            for rel, text in originals.items():
                new, n = rx.subn(replacement, text)
                hits += n
                (tmp / rel).write_text(new)
            _, fail, _, _, _ = scan(tmp, table, every_run, quiet=True)
            for rel, text in originals.items():
                (tmp / rel).write_text(text)

            if hits == 0:
                verdict = "NO SUCH TEXT (nothing to inject into)"
                missed.append(desc)
            elif fail > base_fail:
                verdict = f"CAUGHT ({fail - base_fail} mismatch(es))"
            else:
                verdict = "!! MISSED -- the pattern is vacuous"
                missed.append(desc)
            print(f"  {verdict:38s} {desc}  [{hits} substitution(s)]")

    if missed:
        print(f"\nSELF-TEST FAILED: {len(missed)} injected error(s) not caught:")
        for d in missed:
            print(f"  - {d}")
        return 1
    print(f"\nSELF-TEST PASSED: all {len(rows)} injected errors were caught.")
    return 0


def _reach_report(counts: dict[str, int], reach: dict[str, set[str]]) -> None:
    """Per pattern: WHICH documents it actually checked, and how many lines.

    WHY REPORTING RATHER THAN GATING, decided 2026-09-07 and recorded so it is
    not re-litigated. `--strict` is per-pattern and GLOBAL: a pattern that used
    to check three documents and now checks one still passes. The obvious fix --
    gate on an expected per-document match count -- was rejected. It would fire
    on every ordinary prose edit that adds or removes a mention, and a check that
    goes red on almost every edit gets ignored, which is how the snapshot byte
    total failed as a gate. 21 of 29 guarded patterns already match exactly one
    line in one file, so such a gate would mostly pin accidents of wording.

    That this is worth printing at all is not hypothetical. While rewriting a
    checklist item on 2026-09-07 I turned "rebuild stages **173 files**" into
    "rebuild now stages **173**" -- and `--strict` stayed GREEN, because the
    pattern still matched the manuscript. One document silently lost its guard.
    Nothing reported it; the count in the summary line simply went down by one,
    which reads as prose churn. This section makes that visible.
    """
    if not counts:
        return
    print("PATTERN REACH -- which documents each pattern actually checked.\n")
    print("A pattern that used to check three documents and now checks one "
          "still\npasses `--strict`, which only asks whether it matched "
          "ANYTHING. Shrinkage\nis reported here rather than gated; see this "
          "function's docstring for why.\n")
    single, none = 0, 0
    for desc in sorted(counts):
        files = sorted(reach.get(desc, ()))
        n = counts[desc]
        if not files:
            none += 1
            print(f"  (0 files)  {desc}   <-- VACUOUS: matched nothing")
            continue
        if len(files) == 1:
            single += 1
        print(f"  {len(files)} file(s), {n:3d} line(s)  {desc}")
        for f in files:
            print(f"        {f}")
    print(f"\n  {len(counts)} pattern(s): {none} matching nothing, {single} "
          f"reaching exactly one document,\n  "
          f"{len(counts) - none - single} reaching two or more. A pattern that "
          "drops a document\n  keeps passing -- compare this list against the "
          "last session's before trusting it.\n")


def coverage_report(table, every_run, strict: bool = False) -> int:
    """Per document: which numbers are checked, and which are merely present.

    WHY THIS EXISTS, and why a claim COUNT was never coverage. The check prints
    "checked N claims across M files", which sounds like a coverage figure and
    is not one: it counts the claims the patterns FOUND, and says nothing about
    the numbers no pattern looks at. That is exactly the hole "78 automated
    tests" hid in -- a stale number sitting in the manuscript, in a file the
    checker scanned, reported by nothing, because no pattern was pointed at it.
    221 checked claims looked like coverage while that was true.

    This lists the gap instead of implying it is empty. It is a REPORT, not a
    gate: it always exits 0. Most unchecked literals are legitimately not
    claims about a result -- tolerances, gene counts, dates in prose, figure
    dimensions -- and failing on them would make the tool ignorable. The point
    is that a human can now read the list and see what is unguarded.
    """
    cov: dict[str, list] = {}
    counts: dict[str, int] = {}
    reach: dict[str, set[str]] = {}
    checked, failures, _, skipped_files, _ = scan(
        REPO, table, every_run, quiet=True, coverage=cov,
        match_counts=counts, match_files=reach)

    print("COVERAGE REPORT -- which numeric literals each document has "
          "checked.\n")
    print("A claim count is not a coverage measure. This is the list of "
          "numbers that\nappear in the documents and that NO pattern in this "
          "file examines.\n")

    tot_seen = tot_cov = tot_style = 0
    for rel, (seen, ncov, unchecked, in_style) in cov.items():
        tot_seen += seen
        tot_cov += ncov
        tot_style += in_style
        pct = (100.0 * ncov / seen) if seen else 100.0
        print(f"{rel}")
        print(f"    {seen:4d} numeric literals, {ncov:3d} checked "
              f"({pct:.0f}%), {len(unchecked)} unchecked and not filtered "
              "as structural")
        if in_style:
            print(f"         (+{in_style} more inside <style>/<script>, "
                  "excluded from the denominator: a stylesheet cannot state "
                  "a result)")
        for i, lit, ctx in unchecked[:14]:
            print(f"      :{i:<5d} {lit:<12s} {ctx}")
        if len(unchecked) > 14:
            print(f"      ... and {len(unchecked) - 14} more")
        print()

    _reach_report(counts, reach)

    print(f"TOTAL: {tot_cov} of {tot_seen} numeric literals checked "
          f"({100.0 * tot_cov / max(tot_seen, 1):.0f}%), "
          f"{checked} claim comparisons, {failures} mismatch(es).")
    if tot_style:
        print(f"       {tot_style} further literal(s) excluded as stylesheet "
              "or script content.")
    for entry in skipped_files:
        print(f"  skipped: {entry}")
    if not strict:
        print("\nThis report never fails without --strict. Read it and decide "
              "which unchecked\nnumbers are claims about a result; those are "
              "the ones worth a pattern and a\nself-test row.")
        return 0
    return _coverage_floor_verdict(cov)


# How many literals each document must still have CHECKED. Recorded 2026-09-07
# at the exact live counts, which is safe for a reason that was measured rather
# than assumed -- see `_coverage_floor_verdict`.
COVERAGE_FLOOR = {
    "09-PAPER-DRAFT.md": 157,
    "07-ABSTRACT-DRAFT.md": 25,
    "14-SCIENCE-AUDIT.md": 103,
    "15-CHECKLIST.md": 17,
    "16-YOUR-TASKS.md": 9,
    "pipeline/README.md": 21,
    "poster/poster.html": 19,
    "submission/SUBMISSION-CHECKLIST.md": 8,
    "submission/COVER-LETTER.md": 7,
    "pipeline/results/README.md": 58,
    "CITATION.cff": 2,
    "18-SUPPLEMENTARY-INVENTORY.md": 8,
    "11-SUBMISSION-PACK.md": 3,
}


def _coverage_floor_verdict(cov: dict) -> int:
    """Fail when a document has FEWER literals checked than it used to.

    WHY A COUNT AND NOT A PERCENTAGE. The obvious gate is a floor on each
    document's checked FRACTION, and it is the wrong one: the denominator moves
    for reasons that have nothing to do with losing a guard. Measured across one
    session of ordinary documentation work on 2026-09-07 -- adding a measured
    table to `pipeline/README.md` and a log row to `results/README.md`, both
    dense with exit codes and counts:

        percentages   2 of 13 documents fell by a point
        checked COUNTS  all 13 unchanged, exactly

    Writing prose about numbers lowers the fraction while removing nothing. A
    fraction floor would therefore have fired twice in one session for no
    defect at all, and a gate that goes red on ordinary work gets ignored --
    the same reasoning that made the snapshot byte total a dated observation
    instead of a check. The count only falls when a pattern really has stopped
    matching something it used to match.

    This is the AGGREGATE guard. Its per-pattern counterpart is the reach
    report, which catches the sharper case: a pattern that still matches
    somewhere, so `--strict` stays green, while quietly dropping a document.
    That is a real defect, caused in this repository on 2026-09-07 by rewording
    a checklist item.
    """
    lost = []
    for rel, floor in sorted(COVERAGE_FLOOR.items()):
        entry = cov.get(rel)
        if entry is None:
            lost.append(f"  {rel}: not scanned at all (floor {floor})")
            continue
        ncov = entry[1]
        if ncov < floor:
            lost.append(f"  {rel}: {ncov} checked, was {floor} "
                        f"(-{floor - ncov})")
    if lost:
        print("\nCOVERAGE REGRESSED -- these documents have fewer literals "
              "checked than\nthey did when the floor was recorded:")
        print("\n".join(lost))
        print("\nEither a guard stopped matching (find it in the reach report "
              "above), or the\nclaim was deliberately removed -- in which case "
              "lower the floor in\nCOVERAGE_FLOOR, in the same commit, with "
              "the reason.")
        return 1
    print(f"\nOK: all {len(COVERAGE_FLOOR)} documents still have at least as "
          "many literals\nchecked as when the floor was recorded.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--list", action="store_true",
                    help="print the authority table and exit")
    ap.add_argument("--verbose", action="store_true",
                    help="print every matched claim, not only mismatches")
    ap.add_argument("--self-test", action="store_true",
                    help="inject wrong values and assert each is caught")
    ap.add_argument("--fast", action="store_true",
                    help="skip the pytest collection and abstract recount")
    ap.add_argument("--coverage", action="store_true",
                    help="list, per document, which numbers are NOT checked")
    ap.add_argument("--strict", action="store_true",
                    help="also fail if any CLAIMS pattern matched NOTHING "
                         "anywhere -- a vacuous pattern passes silently")
    args = ap.parse_args()

    table = authority(slow=not args.fast)
    every_run = every_stored_run()
    if args.list:
        print("AUTHORITATIVE VALUES, read from the frozen artefacts just now\n")
        for k, v in sorted(table.items()):
            print(f"  {k:28s} {v!r}")
        return 0
    if args.self_test:
        return self_test(table, every_run)
    if args.coverage:
        return coverage_report(table, every_run, strict=args.strict)

    print("Checking documents against the frozen results.")
    print(f"Authority: {len(table)} values from results/nsclc_v3, "
          f"results/pancancer_v3,\n           results/partition_variance_nsclc, "
          f"pancancer_v3/science_gaps.json,\n           plus a live pytest "
          f"collection and abstract recount; and {len(every_run)}\n"
          f"           values across all stored runs, for provenance quotes.\n")

    matches: dict[str, int] = {}
    checked, failures, unmatched, skipped_files, skipped_lines = scan(
        REPO, table, every_run, verbose=args.verbose, match_counts=matches)
    # A pattern that matches nothing PASSES. That is not hypothetical: three
    # consecutive sessions shipped one, and each time a green run said
    # "OK: every checked claim agrees" while the pattern checked no line in
    # the repository. The only way it was ever caught was a human grepping
    # --verbose for the description. `--strict` makes it mechanical.
    vacuous = sorted(d for d, n in matches.items() if n == 0)
    files_seen = len(DOCS) - len(skipped_files)

    print(f"\nchecked {checked} claims across {files_seen} files "
          f"({len(CLAIMS)} scalar + {len(SET_CLAIMS)} set + {len(FRACTIONS)} "
          f"fraction + {len(DERIVED)} derived\npatterns, {len(TRIPLES)} "
          "bracketed-triple sources, plus the full-precision rule)")
    if unmatched:
        print(f"\n{len(unmatched)} bracketed triple(s) matched no known "
              "quantity -- NOT checked, and not\nnecessarily wrong (HPC4 "
              "diagnostics and ssGSEA arms are legitimately absent\nfrom the "
              "authority table):")
        for u in unmatched[:12]:
            print(f"  {u}")
        if len(unmatched) > 12:
            print(f"  ... and {len(unmatched) - 12} more")
    if vacuous:
        print(f"\n{len(vacuous)} CLAIMS pattern(s) matched NOTHING anywhere. "
              "Each passes by\nconstruction and checks no document:")
        for d in vacuous:
            print(f"  VACUOUS  {d}")
        print("  Fix the pattern, or delete it. A claim nobody makes needs no "
              "check;\n  a claim made in prose the pattern cannot see is the "
              "dangerous case.")
    if skipped_files or skipped_lines:
        print("skipped:")
        for entry in skipped_files + skipped_lines:
            print(f"  {entry}")
    print("\nCOVERAGE: only the patterns in the tables above are checked. A "
          "number not\nlisted there is NOT verified by this script -- adding "
          "one is a two-line edit,\nand `--self-test` proves the existing ones "
          "can still fail.")

    if failures:
        print(f"\nFAILED: {failures} mismatch(es).")
        return 1
    if args.strict and vacuous:
        print(f"\nFAILED (--strict): {len(vacuous)} pattern(s) matched "
              "nothing. Every claim agrees,\nbut some of them agree about "
              "no text at all.")
        return 1
    print("\nOK: every checked claim agrees with the frozen results.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
