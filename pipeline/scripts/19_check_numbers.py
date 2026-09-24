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
import importlib.util
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
    # Added 2026-09-17 (F3.3). Both are dated records of an earlier state: the
    # brief was written on 2026-08-11 before any data were analysed, and the
    # code audit on 2026-08-18. Scanning them un-exempted flagged eleven and four
    # "mismatches" that are what was true then (78 tests, NSCLC's 68 sites, a
    # 0.0315 reliability gap later corrected) plus the brief quoting a prior
    # paper's 28 cancer types.
    "00-PROJECT-BRIEF.md":
        "historical record written 2026-08-11, before any data were analysed",
    "04-CODE-AUDIT.md":
        "historical record of the 2026-08-18 code audit",
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


def _count_beats(path: Path) -> int:
    """How many signatures in a run's `immune_excess.csv` beat their null."""
    return sum(r["beats_null"].strip().lower() in ("true", "1") for r in _rows(path))


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


def _a16_extent() -> dict[str, float] | None:
    """A16's extent: how many signature scores describe tissue other than the tumour.

    Session 57. The manuscript's tumor-only paragraph, `pipeline/README.md` and
    the audit all state these counts -- 51 of 944 NSCLC (48 adjacent normal,
    2 recurrence, 1 normal only), 578 of 7,168 pan-cancer profiles averaged and
    40 with no tumour sample, 618 in all -- and NO artefact stores them: they
    lived only as text, in those documents and in `00_build_interim.py`'s
    docstring. A mutation sweep of the README moved them and nothing failed,
    and the two-digit ones never even reached the coverage list, which drops
    `^\\d{1,2}$` as structural.

    Recomputed from the interim expression INDEX only (`columns=[]`, ~2 s),
    with the builder's own rules: NSCLC takes each patient's FIRST TOIL sample
    in TOIL column order; pan-cancer averages every sample of a patient, and
    `_sample_mix` there defines "more than one sample" and "no primary-tumour
    sample". Measured 2026-09-22: every documented count reproduced exactly,
    and 578 + 40 = 618 is both the sum and the union (no patient is in both).

    Session 58 added the breakdown the audit's A16 block states and nothing
    checked: 6,550 single tumour profiles; 520 / 27 / 20 / 11 inside the 578;
    28 / 8 / 4 inside the 40; the five most affected types; NSCLC's 851 and 90.
    Measured 2026-09-23 against the same index: every one reproduced exactly.

    Returns None if the interim data is absent, which is reported rather than
    silently skipped.
    """
    interim = ROOT / "data" / "interim"
    expr_p, meta_p = interim / "expression_hugo.parquet", interim / "meta.parquet"
    nsclc_p = interim / "cohort_nsclc.parquet"
    if not (expr_p.exists() and meta_p.exists() and nsclc_p.exists()):
        return None
    try:
        import numpy as np
        import pandas as pd

        sys.path.insert(0, str(ROOT / "src"))
        from aacr27 import data as _data

        ids = pd.Index(pd.read_parquet(expr_p, columns=[]).index).astype(str)
        codes = pd.Series(ids.str[13:15], index=ids.str[:12])
        out: dict[str, float] = {}

        ns = set(pd.read_parquet(nsclc_p)["patient_id"])
        c_ns = codes[codes.index.isin(ns)]
        first = c_ns[~c_ns.index.duplicated(keep="first")]
        has01 = c_ns.groupby(level=0).agg(lambda c: "01" in set(c))
        off = first[first != "01"]
        out["a16_ns_not_tumour"] = float(len(off))
        out["a16_ns_adjacent_normal"] = float(
            sum(off[p] == "11" and has01[p] for p in off.index))
        out["a16_ns_recurrence"] = float((off == "02").sum())
        out["a16_ns_normal_only"] = float(sum(not has01[p] for p in off.index))
        # Session 58: the audit's A16 narrative block breaks these down.
        ns_n = c_ns.groupby(level=0).size()
        ns_set = c_ns.groupby(level=0).agg(frozenset)
        out["a16_ns_single_tumour"] = float(
            ((ns_n == 1) & ns_set.map(lambda s: s == {"01"})).sum())
        out["a16_ns_tumour_normal"] = float(ns_set.map(lambda s: {"01", "11"} <= s).sum())

        meta = pd.read_parquet(meta_p)
        pat = _data.collapse_to_patient(meta, np.zeros((len(meta), 1)))[0]
        pan = set(pat["patient_id"])
        per = codes[codes.index.isin(pan)].groupby(level=0).agg(lambda c: tuple(sorted(c)))
        multi = per.map(len) > 1
        no01 = ~per.map(lambda t: "01" in t)
        out["a16_pan_averaged"] = float(multi.sum())
        out["a16_pan_no_tumour"] = float(no01.sum())
        out["a16_pan_total"] = float((multi | no01).sum())
        # Session 58: the audit's A16 narrative block -- the mix inside the
        # averaged and no-tumour groups, and the five most affected cancer
        # types, each named, so a swapped name fails as well as a moved count.
        # "Others" is the averaged remainder, so the four parts sum to the
        # whole by construction; the no-tumour parts are each counted, so a
        # fourth kind appearing in the data breaks the document's 28 + 8 + 4.
        mix = per.map(frozenset)
        out["a16_pan_single_tumour"] = float((~multi & ~no01).sum())
        for key, code in (("normal", "11"), ("recurrence", "02"), ("metastasis", "06")):
            out[f"a16_pan_tumour_{key}"] = float((multi & mix.map({"01", code}.__eq__)).sum())
            out[f"a16_pan_{key}_only"] = float((no01 & mix.map({code}.__eq__)).sum())
        out["a16_pan_multi_other"] = out["a16_pan_averaged"] - sum(
            out[f"a16_pan_tumour_{k}"] for k in ("normal", "recurrence", "metastasis"))
        typ = pat.set_index("patient_id")["cancer type abbreviation"]
        hit = typ.reindex(per.index[(multi | no01).to_numpy()]).value_counts()
        for t in A16_TOP_TYPES:
            out[f"a16_pan_type_{t.lower()}"] = float(hit.get(t, 0))
        return out
    except Exception:
        return None


def _site_counts() -> dict[str, float] | None:
    """Distinct tissue source sites per cohort, recomputed from the interim data.

    WHY THIS IS LIVE-MEASURED RATHER THAN READ FROM AN ARTEFACT. The manuscript
    and the poster both state "619 tissue source sites", and NO frozen artefact
    contains that number: `site_control.csv` holds only the sites that were
    *evaluable* in the negative control (98 pan-cancer, 15 NSCLC), which is a
    different quantity. Until now 619 was an unsourced count in a manuscript --
    exactly the failure mode this script exists to prevent.

    It is recomputed the same way `12_partition_variance.py`'s `main()` prints it --
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
            ns_tss = pd.read_parquet(nsclc)["tss"]
            out["nsclc_n_sites"] = float(ns_tss.nunique())
            # Session 60 (H89): E15's site attrition -- what Control C's
            # min_site_n=10 drops -- recomputed here rather than recalled.
            sz = ns_tss.value_counts()
            out["e15_ns_sites_dropped"] = float((sz < 10).sum())
            out["e15_ns_patients_dropped"] = float(sz[sz < 10].sum())
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
            sz = meta.loc[keep, "tss"].value_counts()
            out["e15_pan_sites_dropped"] = float((sz < 10).sum())
            out["e15_pan_patients_dropped"] = float(sz[sz < 10].sum())
            out["e15_pan_retained"] = float(keep.sum() - sz[sz < 10].sum())
        return out or None
    except Exception:
        return None


def _abstract_counts() -> dict[str, int] | None:
    """Every figure 08_count_abstract.py prints, run now rather than recalled.

    Keys: `total` and `limit` (always, or None is returned), and `title`,
    `body`, `naive` (title + body INCLUDING spaces) and `words` when the counter
    printed them. Session 51 (2026-09-17) found the dry run quoting a naive
    length of 2,920 that the counter had not printed since its closing-rule fix
    (it prints 2,915), and a title count and a word count nothing checked. All
    four come from the same run of the counter, so they cannot disagree with
    the total the way a second definition of "what is counted" could.
    """
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
    counts = {"total": int(m.group(1).replace(",", "")),
              "limit": int(m.group(2).replace(",", ""))}
    for key, pat in (("title", r"^\s*title\s+(\d+)\s*$"),
                     ("body", r"^\s*body\s+(\d+)\s*$"),
                     ("naive", r"INCLUDING spaces would read (\d+)"),
                     ("words", r"\(body is (\d+) words\)")):
        f = re.search(pat, out.stdout, re.MULTILINE)
        if f:
            counts[key] = int(f.group(1))
    return counts


def _abstract_chars() -> tuple[int, int] | None:
    """(counted, limit) from 08_count_abstract.py, run now rather than recalled."""
    c = _abstract_counts()
    return None if c is None else (c["total"], c["limit"])


def _cover_letter_words() -> dict[str, float]:
    """The cover letter's body length, measured as it would actually be pasted.

    The letter states a body length and a fallback ("delete paragraph 4 and N
    remain"), and the author uses both at submission to decide whether it sets
    at one page. Session 52 found the stated 614 was RIGHT WHEN WRITTEN -- 609
    by this measure at `713d0be`, 2026-09-02 -- and 794 by the time anyone
    re-derived it: every dated correction note since had been appended inside
    the measured span while the number stood still. Nothing checked it, so
    nothing said so. Recomputed here on every run for that reason.

    "As pasted" is the words between the `[Date]` line and the signature line
    that carries the corresponding author's ORCID, with HTML comments removed:
    comments are invisible in rendered markdown and are not part of what goes
    into the cover-letter field. The closing anchor is the FIRST `orcid.org/`
    line after `[Date]` -- the letter's own contact line -- and deliberately
    not the address on it: session 52's first draft anchored on the email
    literal and the deposit's leak audit refused the staged tree, correctly,
    because this script is published. If either anchor moves, this returns
    nothing rather than measure the wrong span: a silently mis-scoped count is
    worse than an unchecked one.
    """
    path = REPO / "submission" / "COVER-LETTER.md"
    try:
        lines = path.read_text(encoding="utf-8").split("\n")
    except OSError:
        return {}
    i0 = next((i for i, ln in enumerate(lines) if ln.strip() == "[Date]"), None)
    if i0 is None:
        return {}
    i1 = next((i for i, ln in enumerate(lines[i0 + 1:], start=i0 + 1)
               if "orcid.org/" in ln.lower()), None)
    if i1 is None or i1 <= i0:
        return {}
    body = re.sub(r"<!--.*?-->", " ", "\n".join(lines[i0:i1]), flags=re.S)
    total = len(body.split())
    out = {"cover_letter_body_words": float(total)}
    paras = re.split(r"\n\s*\n", body)
    for para in paras:
        if para.startswith("The correction that makes this measurable"):
            n4 = len(para.split())
            out["cover_letter_para4_words"] = float(n4)
            out["cover_letter_body_minus_para4"] = float(total - n4)
            break

    # Session 55: the letter ALSO enumerates its remaining paragraphs
    # ("48 / 111 / 112 / 122 words of argument and 215 of required
    # boilerplate"), and the author reads that list to decide which part of
    # the argument to lose (B-24). Every member was exact when measured here,
    # but the enumeration sums to 608 against the 647 the header says remain
    # after paragraph 4 goes -- a 39-word gap a reader cannot resolve, because
    # the sentence never said it counts only the PROSE paragraphs and not the
    # date, address, salutation, sign-off and byline. Measured, not asserted:
    # the paragraphs sum to the body exactly, so the gap is entirely that
    # scaffolding. Each member is bound here so the list cannot drift the way
    # the body length did, and the two sums are bound separately because the
    # project's own headroom defect was every part right and the arithmetic
    # between them wrong.
    #
    # Anchored on each paragraph's opening words -- structure, never personal
    # data, since this script is published. A paragraph whose opening moves
    # contributes NOTHING rather than a wrong span, exactly as the body span
    # above refuses to measure between anchors it cannot find.
    openers = (
        ("cover_letter_para_submit_words", "We submit "),
        ("cover_letter_para_finding_words",
         "**The finding is about tumor biology"),
        ("cover_letter_para_immune_words",
         "**The immune signal itself is real"),
        ("cover_letter_para_notwhat_words", "**What this work is not.**"),
        ("cover_letter_para_scope_words",
         "The manuscript falls within the journal's"),
    )
    enumerated = 0.0
    found_all = True
    for key, opener in openers:
        hit = next((p for p in paras if p.startswith(opener)), None)
        if hit is None:
            found_all = False
            continue
        n = float(len(hit.split()))
        out[key] = n
        enumerated += n
    if found_all and "cover_letter_para4_words" in out:
        out["cover_letter_enumerated_words"] = enumerated
        # Everything in the body that is NOT one of the six prose paragraphs:
        # the date line, the address block, the salutation, the sign-off and
        # the byline. Derived by subtraction from measured parts, so it cannot
        # disagree with them.
        out["cover_letter_scaffold_words"] = (
            float(total) - enumerated - out["cover_letter_para4_words"])
    return out


def _paper_body_words() -> dict[str, float]:
    """The manuscript's body word count and figure count, which AACR requires
    the COVER LETTER to state, and which nothing had ever computed.

    The cover letter has carried `[X] words` and `[N] figures` as placeholders
    since it was written. They are not optional -- the submission pack records
    that AACR asks the cover letter to state both -- and they are a count one
    document makes of ANOTHER, which is the shape session 53 found goes stale
    silently (H30). So they are derived here rather than typed.

    The span is the letter's own definition: everything except the cover page,
    the abstract, Materials and Methods, tables and references -- that is
    Introduction + Results + Discussion, with Data/Code Availability, the
    authorship sections, Acknowledgments and References all excluded as back
    matter. HTML comments and tags are stripped, so the `<details>` markup of
    the folded Supplementary Note contributes nothing, and the measure is the
    same word-character token count `_paper_abstract_words()` uses, so the two
    figures in the letter are commensurable.

    `paper_suppnote_words` is reported beside them because it is the size of
    the decision in B-21: the folded note is inside the counted span today, and
    lifting it into a real supplementary file removes that many words from the
    manuscript.
    """
    path = REPO / "09-PAPER-DRAFT.md"
    try:
        lines = path.read_text(encoding="utf-8").split("\n")
    except OSError:
        return {}

    def _span(name: str) -> tuple[int, int] | None:
        i = next((i for i, ln in enumerate(lines)
                  if ln.startswith("## " + name)), None)
        if i is None:
            return None
        j = next((j for j, ln in enumerate(lines[i + 1:], start=i + 1)
                  if ln.startswith("## ")), len(lines))
        return i + 1, j

    def _words(text: str) -> int:
        text = re.sub(r"<!--.*?-->", " ", text, flags=re.S)
        text = re.sub(r"<[^>]+>", " ", text)
        text = re.sub(r"[*_`]", "", text)
        return len([w for w in text.split() if re.search(r"\w", w)])

    spans = [_span(n) for n in ("Introduction", "Results", "Discussion")]
    if any(s is None for s in spans):
        return {}
    body = sum(_words("\n".join(lines[a:b])) for a, b in spans)  # type: ignore

    figs = {int(m) for m in re.findall(r"(?:Fig\.|Figure)\s+(\d+)",
                                       "\n".join(lines))}
    out = {"paper_body_words": float(body)}
    if figs:
        out["paper_figure_count"] = float(len(figs))

    # Supplementary Note 1. Until 2026-09-24 it was a folded <details> block
    # inside Limitation 8, measured here from the manuscript, and its 1,593
    # words sat inside the counted body. B-21 lifted it out word for word into
    # its own file; the body above no longer contains it, and its size is now
    # measured from that file -- the lifted text after the note's first
    # italic-marked subsection heading, which is where the old block's body
    # began (the old count also included the block's <summary>, which the
    # file's own introduction replaces).
    note = REPO / "SUPPLEMENTARY-NOTE-1.md"
    if note.exists():
        nl = note.read_text(encoding="utf-8").split("\n")
        start = next((i for i, ln in enumerate(nl)
                      if ln.startswith("*The mechanism,")), None)
        if start is not None:
            out["paper_suppnote_words"] = float(_words("\n".join(nl[start:])))
    return out


def _paper_abstract_words() -> dict[str, float]:
    """The manuscript abstract's own word counts, and the sum of them.

    09-PAPER-DRAFT.md states that its abstract "measures 234 + 28 = 262 words"
    and is therefore "inside the observed ceiling" of 131-264 words that the
    same paragraph reports across 40 published articles. That is the cover
    letter's shape again (session 52, S52-1/S52-2): a document counting ITSELF,
    with a conclusion resting on the count. Session 53 checked the history
    before assuming drift and found NONE -- the body measures the same today as
    at `509be54`, the commit that wrote the numbers, and 234/28 are exactly
    reproduced by the measure below. The defect is not staleness; it is that
    nothing re-derived any of the three, and one pending author decision (B-30)
    edits a line INSIDE the measured span.

    The measure: whitespace-delimited tokens containing a word character,
    between the body's opening line and the `**Significance:**` line, and from
    that line to the `<details>` that folds the superseded structured version.
    The word-character test drops four tokens a journal would not call words --
    three spaced em-dashes and the `=` of a reported statistic. A naive
    whitespace split counts those four and gives 238 + 28 = 266, which is why
    the paragraph now states its method: the ceiling claim is true under this
    measure and not under the naive one, and a reader is entitled to know which.

    Anchors are structural. If any moves this returns nothing rather than
    measure the wrong span.
    """
    path = REPO / "09-PAPER-DRAFT.md"
    try:
        lines = path.read_text(encoding="utf-8").split("\n")
    except OSError:
        return {}
    i0 = next((i for i, ln in enumerate(lines)
               if ln.startswith("Deep learning models predict")), None)
    i1 = next((i for i, ln in enumerate(lines)
               if ln.startswith("**Significance:**")), None)
    if i0 is None or i1 is None or i1 <= i0:
        return {}
    i2 = next((i for i, ln in enumerate(lines[i1 + 1:], start=i1 + 1)
               if ln.startswith("<details")), None)
    if i2 is None:
        return {}

    def _words(text: str) -> int:
        text = re.sub(r"[*_`]", "", text)
        return len([w for w in text.split() if re.search(r"\w", w)])

    body = _words("\n".join(lines[i0:i1]))
    sig = _words(re.sub(r"^\s*\*\*Significance:\*\*\s*", "",
                        "\n".join(lines[i1:i2])))
    out = {
        "paper_abstract_body_words": float(body),
        "paper_abstract_significance_words": float(sig),
        "paper_abstract_total_words": float(body + sig),
    }
    # Session 54: the folded structured abstract's length, which the same
    # paragraph states twice. It said 461 -- a naive whitespace split of the
    # block -- in the paragraph that defines "word" as the measure above,
    # under which the block is 456. Counted here by that measure, from the
    # line after its <summary> to its </details>.
    s0 = next((i for i, ln in enumerate(lines[i2:], start=i2)
               if ln.startswith("<summary>Superseded structured abstract")), None)
    s1 = None if s0 is None else next(
        (i for i, ln in enumerate(lines[s0:], start=s0) if ln.startswith("</details>")),
        None)
    if s0 is not None and s1 is not None:
        out["paper_structured_abstract_words"] = float(
            _words("\n".join(lines[s0 + 1:s1])))
    return out


# The dry run's pre-measured abstract variants (session 51). Each entry is
# (key, prefix, suffix): the claim pattern is prefix + "(number)" + suffix, and
# the self-test row rewrites the same cell. Prefixes and suffixes must not
# contain a capturing group. "(?m)" lets the self-test, which substitutes over
# whole-file text, anchor on a table row the way the line scan does.
_NUM = r"[−-]?[\d,]+"
_B = r"\**"  # an optional bold marker around a cell
ABSTRACT_VARIANT_SPEC: list[tuple[str, str, str]] = []
for _n in (1, 2, 3, 4):
    _row = rf"(?m)^\| {_n} \| [^|]+ \| "
    ABSTRACT_VARIANT_SPEC += [
        (f"ladder_cut{_n}_saves", _row, r" \|"),
        (f"ladder_cut{_n}_total", _row + r"\d+ \| ", r" \|"),
        (f"ladder_cut{_n}_headroom", _row + r"\d+ \| [\d,]+ \| " + _B, _B + r" \|"),
    ]
# Session 58: the AI-sentence, cut-4a and B-30 variant tables were retired on
# 2026-09-23, when the author's decisions went into the draft -- every row had
# been defined against the 2,512-character draft, which no longer exists. The
# ladder above is what remains live. (They are in git history at `dc461c4`.)

# Session 57: the dry run's PROSE that re-quotes its own tables. The tables
# above have been checked claims since session 51; the sentences that quote
# them back were not, and that is how the ladder's closing sentence kept a
# pre-fix "226" after the tables were re-derived on 2026-09-17 (session 56,
# H55). A mutation sweep of the dry run -- each of the 42 literals `--coverage`
# lists as unchecked moved one step in its last printed place, `--strict`
# without `--fast`, restored by sha256, behind a positive control on a ladder
# cell that WAS caught -- came back 41 UNCAUGHT. Nine are live figures in such
# sentences. Two more on the same line (the "/ 88" and "/ 24" headrooms) were
# never in that list at all: the coverage filter drops them as structural, and
# only mutating the line showed they were unguarded too. The rest are ORCID
# digit groups, HPC4 job ids, AACR's own fees and table cost, two tick sizes
# with no authority, and dated notes that must keep their old values.
# Same (key, prefix, suffix) shape as ABSTRACT_VARIANT_SPEC; keys repeat
# where one authority is quoted twice.
DRYRUN_PROSE_SPEC: list[tuple[str, str, str]] = [
    ("abstract_chars", r"cumulative from ", r" and were re-derived"),
    ("ladder_cut3_headroom", r"portal demands more than ", r" characters"),
]


def _abstract_variants() -> dict[str, float]:
    """The dry run's cut ladder, recomputed from the draft.

    Session 51. `submission/ABSTRACT-PORTAL-DRY-RUN.md` §4 carries a cut ladder
    and (until 2026-09-23) tables of totals for the 2027 AI sentences, all
    measured by hand with the counter, and every one goes stale the moment the
    abstract changes.
    The rows DEFINE their variants -- the backticked text a cut removes, the
    quoted sentences a variant adds -- and the NUMBERS in the rows are the
    claims. So the definitions are read from the rows, applied to the live
    draft, and counted with 08_count_abstract.py's own `extract` and `count`
    (imported, not re-implemented). A cut whose text no longer occurs exactly
    once in the body yields no authority, so `--strict` fails on the vanished
    keys: when the abstract changes, the ladder has to be re-derived.
    """
    dry_p = REPO / "submission" / "ABSTRACT-PORTAL-DRY-RUN.md"
    draft_p = REPO / "07-ABSTRACT-DRAFT.md"
    script = ROOT / "scripts" / "08_count_abstract.py"
    if not (dry_p.exists() and draft_p.exists() and script.exists()):
        return {}
    spec = importlib.util.spec_from_file_location("_count_abstract", script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    title, body = mod.extract(draft_p.read_text())
    base, limit = mod.count(title) + mod.count(body), mod.LIMIT
    dry = dry_p.read_text()

    def saves(old: str, new: str = "") -> int | None:
        # Whitespace-tolerant, because the draft wraps where the table does not.
        rx = r"\s+".join(re.escape(w) for w in old.split())
        if len(re.findall(rx, body)) != 1:
            return None
        return mod.count(old) - mod.count(new)

    out: dict[str, float] = {}
    cut: dict[str, int] = {}
    for m in re.finditer(r"^\| ([1-4]) \| `([^`]+)`(?: → `([^`]+)`)? \|", dry, re.M):
        s = saves(m.group(2), m.group(3) or "")
        if s is not None:
            cut[m.group(1)] = s

    def put(key: str, total: int) -> None:
        out[f"{key}_total"], out[f"{key}_headroom"] = float(total), float(limit - total)

    running = base
    for n in ("1", "2", "3", "4"):
        if n not in cut:
            break  # cumulative: a missing rung invalidates every rung below it
        running -= cut[n]
        out[f"ladder_cut{n}_saves"] = float(cut[n])
        put(f"ladder_cut{n}", running)
    return out


def _snapshot_size() -> tuple[int, int] | None:
    """(staged file count, staged bytes) from a LIVE rebuild, not from memory.

    Three documents quote this pair and nothing verified it, so it drifted
    every time anyone edited a staged file: 169 -> 170 -> 172 -> 173 files and
    1,447,396 -> 1,493,616 -> 1,756,296 -> 1,764,887 bytes inside two sessions,
    twice leaving a stale pair in the manuscript. It is derived by RUNNING
    `22_build_public_snapshot.py` into a temp directory and reading the numbers
    it prints, rather than by re-implementing its staging rules here -- a second
    definition of "what is staged" would be free to disagree with the first, and
    this project has already been bitten by exactly that (S1's two definitions
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
        # Added 2026-09-16 (ledger D4): the number of patients the pooled index
        # actually resamples once incomplete cases are dropped, which the run
        # measures and records only in these notes. (The same notes also carry
        # "alone reaches C=0.676" -- but that is a LITERAL the code prints for
        # every cohort, so it is read from results/type_only_concordance.json
        # below instead.)
        m = re.search(r"(\d+) of \d+ patients dropped as incomplete", notes)
        if m:
            table[f"{prefix}_isi_dropped"] = float(m.group(1))
        m = re.search(r"(\d+) complete-case patients resampled", notes)
        if m:
            table[f"{prefix}_isi_complete_n"] = float(m.group(1))

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
            # Session 58: "BH-FDR changes nothing: 16/16 and 10/32 survive" is a
            # claim about the stored BH column, which only pan-cancer carries.
            if "beats_null_bh" in rows[0]:
                table[f"{prefix}_outcome_bh_beat"] = float(
                    sum(x["beats_null_bh"] == "True" for x in rows))
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
        # Session 58: Limitation 8's "All 16 of 16 signatures still exceed
        # their nulls" after the corrected orderings, read from the same run.
        table["sortaudit_beats"] = float(_count_beats(p.parent / "immune_excess.csv"))
    else:
        missing.append(str(p.relative_to(REPO)))
    # Session 59: HPC4's corrected-ordering run of the same configuration (job
    # 129151), copied into the diagnostics directory by A6 and hashed since
    # 2026-09-24. Its ISI agrees with macOS to one ULP; its Delta-MAE does not
    # (A12, the platform-dependent random-patient split), and the audit and the
    # manuscript quote this one. NOT registered as a TRIPLE, measured: the
    # triple rule identifies by point estimate with a 5e-5 floor, so it claimed
    # the tumour-only 0.004844... [0.00151, 0.00819] in the results README and
    # failed its bounds. The BRACKETED copies were already guarded -- a triple
    # must match some stored run at printed precision or fail -- so only the
    # two UNBRACKETED copies in the audit needed binding, below.
    p = RESULTS / "session48_diagnostics" / "job129151_nsclc_stablesort_hpc4" / "summary.json"
    if p.exists():
        dm = json.loads(p.read_text())["secondary_delta_mae"]
        # (Assigned directly: `_triple` is defined further down this function.)
        table["hpc4_sortaudit_dmae"] = float(dm["value"])
        table["hpc4_sortaudit_dmae_lo"] = float(dm["ci_lo"])
        table["hpc4_sortaudit_dmae_hi"] = float(dm["ci_hi"])
    else:
        missing.append(str(p.relative_to(REPO)))
    # Session 60 (H86): the PINNED-split NSCLC run -- the author's A12 decision
    # (option b) made permanent by `splits.stable_group_kfold`. Its ISI is the
    # corrected-ordering run's bit for bit (the pin never reaches the
    # preserved-site path), so only the random-patient secondaries are new: the
    # Delta-MAE, which reproduces ledger D1's in-process patch exactly, and the
    # site-prediction control's median AUROC, which D1 never stored. Both are
    # reported BESIDE the unpinned values, never replacing them. Scalars, not a
    # registered TRIPLE (#59: a triple identifies by point estimate and would
    # claim the other Delta-MAE intervals within printed tolerance); bracketed
    # copies are guarded by the stored-run rule, which reads this summary.
    # The pan-cancer pinned run (macOS, 2026-09-24; the VPN was down) is read
    # the same way once it exists; NSCLC's is required.
    for coh, run, required in (("nsclc", "nsclc_v3_pinned", True),
                               ("pancancer", "pancancer_v3_pinned", False)):
        d = RESULTS / run
        if not ((d / "summary.json").exists() and (d / "site_control.csv").exists()):
            if required:
                missing.append(str((d / "summary.json").relative_to(REPO)))
            continue
        dm = json.loads((d / "summary.json").read_text())["secondary_delta_mae"]
        table[f"pinned_{coh}_dmae"] = float(dm["value"])
        table[f"pinned_{coh}_dmae_lo"] = float(dm["ci_lo"])
        table[f"pinned_{coh}_dmae_hi"] = float(dm["ci_hi"])
        table[f"pinned_{coh}_site_auroc"] = float(
            _median([float(x["auroc"]) for x in _rows(d / "site_control.csv")]))
        table[f"pinned_{coh}_n_sites"] = float(len(_rows(d / "site_control.csv")))
        ps = _rows(d / "per_signature.csv")
        table[f"pinned_{coh}_median_dr"] = float(_median([float(x["delta_r"]) for x in ps]))
        table[f"pinned_{coh}_emb_over_cov"] = float(
            _median([float(x["embedding_over_covariates"]) for x in ps]))
        table[f"pinned_{coh}_cov_r"] = float(_median([float(x["r_covariates"]) for x in ps]))
        table[f"pinned_{coh}_emb_over_cov_abs"] = abs(table[f"pinned_{coh}_emb_over_cov"])

    # The SAME sensitivity run for pan-cancer, produced on Linux/x86_64 because
    # the pan-cancer working set never fitted on the macOS machine. Reported
    # beside `pancancer_v3/`, never replacing it -- option B, as for NSCLC.
    # The shift is DERIVED from the pair and taken as a MAGNITUDE: NSCLC's
    # corrected value moved DOWN and pan-cancer's moved UP, and the prose in
    # both places says "moves by", so a signed authority would make one of the
    # two documents wrong for a reason that has nothing to do with the science.
    p = RESULTS / "pancancer_v3_stablesort" / "summary.json"
    if p.exists():
        j = json.loads(p.read_text())["primary_excess"]
        table["pcsortaudit_isi"] = float(j["value"])
        table["pcsortaudit_ci_lo"] = float(j["ci_lo"])
        table["pcsortaudit_ci_hi"] = float(j["ci_hi"])
        frozen = RESULTS / "pancancer_v3" / "summary.json"
        if frozen.exists():
            f = json.loads(frozen.read_text())["primary_excess"]
            table["pcsortaudit_shift"] = abs(float(j["value"]) - float(f["value"]))
        table["pcsortaudit_beats"] = float(_count_beats(p.parent / "immune_excess.csv"))
    else:
        missing.append(str(p.relative_to(REPO)))

    # The NSCLC outcome arm's power at given C-index advantages (ledger E10),
    # written by 30_outcome_power.py from the frozen outcome table, in percent.
    p = RESULTS / "outcome_power_nsclc.json"
    if p.exists():
        curve = json.loads(p.read_text())["power_by_c_index_advantage"]
        for key, name in (("0.02", "002"), ("0.03", "003"), ("0.04", "004")):
            table[f"nsclc_power_{name}_pct"] = 100 * float(curve[key]["median"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # Cancer type ALONE on pooled pan-TCGA PFI (ledger E13): computed by
    # 28_type_only_concordance.py from the cohort, because the 0.676 in the
    # run notes is a constant string, not a measurement.
    p = RESULTS / "type_only_concordance.json"
    if p.exists():
        table["pancancer_type_only_c"] = float(
            json.loads(p.read_text())["c_index"]["event_proportion"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # S1's multiplicity columns for the FROZEN NSCLC cohort. DERIVED here, from
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
        # Asserted, not stored: an authority key nothing reads is the defect
        # `models.site_prediction_control`'s dead `.attrs` were (removed 2026-09-06). The 16/16 count itself is checked by the
        # FRACTIONS mechanism wherever a document writes it.
        assert int(_np.nansum(_q < 0.05)) == len(_d), (
            "frozen NSCLC no longer gives 16/16 after BH")
        # Session 58: the comment above was only half true. The FRACTIONS net
        # fires on "signatures beat|exceed|beating", and the manuscript's own
        # BH sentence writes the count as "N of 16", which no net reads -- so
        # the count is now an authority, bound where it is written out.
        table["nsclc_bh_beat"] = float(
            int(((_q < 0.05) & (_d["excess_z"].to_numpy() > 0)).sum()))
    else:
        missing.append(str(ie.relative_to(REPO)))
    # Session 58: the same count pan-cancer, from the stored BH columns.
    ie = RESULTS / "pancancer_v3" / "immune_excess.csv"
    if ie.exists():
        table["pancancer_bh_beat"] = float(sum(
            r["beats_null_bh"].strip().lower() in ("true", "1") for r in _rows(ie)))
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
        # Session 54: the chi-square interval for the between-partition sd and
        # the variance shares it implies, which the Results quote for these
        # six- and five-partition sweeps. They predate `27_partition_sweep.py`
        # and never stored it, so it is derived here with that script's own
        # formula (df = partitions - 1) rather than typed.
        import math
        from scipy import stats as sps
        df = vals["n_partitions_used"] - 1
        sd_lo = sd * math.sqrt(df / sps.chi2.ppf(0.975, df))
        sd_hi = sd * math.sqrt(df / sps.chi2.ppf(0.025, df))
        table[f"{prefix}_partition_sd_chi2_lo"] = sd_lo
        table[f"{prefix}_partition_sd_chi2_hi"] = sd_hi
        table[f"{prefix}_partition_share_lo_pct"] = 100.0 * sd_lo**2 / (sd_lo**2 + se**2)
        table[f"{prefix}_partition_share_hi_pct"] = 100.0 * sd_hi**2 / (sd_hi**2 + se**2)
    if "nsclc_bootstrap_se" in table and "pancancer_bootstrap_se" in table:
        table["bootstrap_se_ratio"] = (table["nsclc_bootstrap_se"]
                                       / table["pancancer_bootstrap_se"])
        table["partition_sd_ratio"] = (table["nsclc_partition_sd"]
                                       / table["pancancer_partition_sd"])

    # B9 (session 49): the 24-partition sweeps under the pinned ordering and
    # their comparison, from `27_partition_sweep.py derive` / `compare`. Kept
    # under their own `pv24_` keys: the frozen six- and five-partition values
    # above stay reported beside them (Option B), so neither may stand in for
    # the other.
    for cohort in ("nsclc", "pancancer"):
        p = RESULTS / f"partition_variance_{cohort}_stable24" / "derived.json"
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
            continue
        j = json.loads(p.read_text())
        pre = f"pv24_{cohort}_"
        table[pre + "n"] = float(j["n_partitions_used"])
        table[pre + "sd"] = float(j["isi_sd_across_partitions"])
        table[pre + "sd_lo"], table[pre + "sd_hi"] = map(float, j["isi_sd_chi2_95ci"])
        table[pre + "understated_pct"] = 100.0 * float(j["interval_understated_by"])
        table[pre + "share_pct"] = 100.0 * float(j["partition_share_of_total_variance"])
        lo, hi = j["partition_share_from_sd_chi2_95ci"]
        table[pre + "share_lo_pct"], table[pre + "share_hi_pct"] = 100.0 * lo, 100.0 * hi
        # Session 60: the 24-partition runs' own bootstrap and combined SE,
        # from the same run's summary. Until then the results README's
        # "Bootstrap SE 0.011073" (24 partitions) was CLAIMED by the frozen
        # six-partition pattern below and passed only by its 5e-5 floor
        # (0.011100 against 0.011073) -- a different quantity, unguarded.
        sp = RESULTS / f"partition_variance_{cohort}_stable24" / "summary.json"
        if sp.exists():
            sj = json.loads(sp.read_text())
            table[pre + "bootstrap_se"] = float(sj["bootstrap_se"])
            table[pre + "combined_se"] = float(sj["combined_se"])
        else:
            missing.append(str(sp.relative_to(REPO)))
    p = RESULTS / "partition_variance_scaling_stable24.json"
    if p.exists():
        j = json.loads(p.read_text())
        table["pv24_bootstrap_ratio"] = float(j["bootstrap_component_ratio"])
        table["pv24_sqrt_patients"] = float(j["sqrt_patients_prediction"])
        table["pv24_sqrt_sites"] = float(j["sqrt_sites_prediction"])
        table["pv24_partition_ratio"] = float(j["partition_component_ratio"])
        table["pv24_partition_ratio_lo"], table["pv24_partition_ratio_hi"] = map(
            float, j["partition_component_ratio_f_95ci"])
        # The prose states these three verdicts in words; if an input moved
        # enough to flip one, the words would be wrong while every number
        # still matched, so they are asserted here.
        want = {"sqrt_patients_inside_ratio_ci": False, "sqrt_sites_inside_ratio_ci": False,
                "one_inside_ratio_ci": True, "sd_chi2_95ci_overlap": True}
        got = {k: j[k] for k in want}
        if got != want:
            raise SystemExit(f"B8 VERDICT CHANGED: {got} (the Results paragraph says {want})")
    else:
        missing.append(str(p.relative_to(REPO)))

    # E14 (session 49): Control C's within-type calibration for all 16
    # signatures. Script 11's A4b calibrates ONE probe; these records rerun it
    # for every signature, and each must reproduce that probe bit for bit.
    for cohort in ("pancancer", "nsclc"):
        p = RESULTS / "control_c_calibration" / f"{cohort}.json"
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
            continue
        j = json.loads(p.read_text())
        if not j["probe_reproduces_stored_a4b_bitwise"]:
            raise SystemExit(f"E14 BROKEN: {p.name} does not reproduce script 11's probe")
        table[f"e14_{cohort}_median_pct"] = 100.0 * float(j["median_calibrated_r2"])
        table[f"e14_{cohort}_min_pct"] = 100.0 * float(j["min_calibrated_r2"])
        table[f"e14_{cohort}_max_pct"] = 100.0 * float(j["max_calibrated_r2"])
        table[f"e14_{cohort}_n_sig"] = float(j["n_perm_p_at_or_below_0.05"])
        table[f"e14_{cohort}_median_r2"] = float(j["median_calibrated_r2"])

    # Session 49's sensitivity analyses (E1, E3, E4, E5, A16). Each reads the
    # record the run wrote; the medians are recomputed from its tables.
    def _triple(prefix, v, lo, hi):
        table[prefix], table[prefix + "_lo"], table[prefix + "_hi"] = float(v), float(lo), float(hi)

    def _med(path, col):
        return _median(float(r[col]) for r in _rows(path))

    p = RESULTS / "nsclc_sensitivity" / "summary.json"
    if p.exists():
        j = json.loads(p.read_text())
        if not j["reference_reproduced_bitwise"]:
            raise SystemExit("E1-E5 BROKEN: the NSCLC sensitivity reference is not bitwise")
        var = {v["variant"]: v for v in j["variants"]}
        for name in ("reference", "global_axis", "unmatched_null", "null_boot_400",
                     "null_boot_1000", "k3", "k10"):
            _triple(f"sens_ns_{name}", var[name]["isi"], var[name]["lo"], var[name]["hi"])
            table[f"sens_ns_{name}_beats"] = float(var[name]["signatures_beating_null"])
        table["sens_ns_unmatched_shift"] = float(var["unmatched_null"]["shift_from_reference"])
        ie = RESULTS / "nsclc_sensitivity" / "global_axis" / "immune_excess.csv"
        table["sens_ns_global_null_r"] = _med(ie, "null_mean_r")
        table["sens_ns_global_curated_r"] = _med(ie, "r_residual_refit")
        table["sens_ns_global_degenerate"] = float(sum(float(r["degenerate"]) for r in _rows(ie)))
    else:
        missing.append(str(p.relative_to(REPO)))
    p = RESULTS / "pancancer_sensitivity" / "summary.json"
    if p.exists():
        var = {v["variant"]: v for v in json.loads(p.read_text())["variants"]}
        _triple("sens_pan_reference", var["reference"]["isi"], var["reference"]["lo"],
                var["reference"]["hi"])
        g = var["global_axis"]
        _triple("sens_pan_global_axis", g["isi"], g["lo"], g["hi"])
        table["sens_pan_global_axis_beats"] = float(g["signatures_beating_null"])
        ie = RESULTS / "pancancer_sensitivity" / "global_axis" / "immune_excess.csv"
        table["sens_pan_global_curated_r"] = _med(ie, "r_residual_refit")
        table["sens_pan_global_null_r"] = _med(ie, "null_mean_r")
        ie = RESULTS / "pancancer_v3_stablesort" / "immune_excess.csv"
        table["sens_pan_ref_curated_r"] = _med(ie, "r_residual_refit")
        table["sens_pan_ref_null_r"] = _med(ie, "null_mean_r")
    else:
        missing.append(str(p.relative_to(REPO)))

    for name, path in (("ns", RESULTS / "nsclc_v3_tumouronly" / "summary.json"),
                       ("pan", RESULTS / "pancancer_tumouronly" / "summary.json")):
        if not path.exists():
            missing.append(str(path.relative_to(REPO)))
            continue
        j = json.loads(path.read_text())
        pe = j.get("primary_excess") or {"value": j["isi"], "ci_lo": j["lo"], "ci_hi": j["hi"]}
        _triple(f"a16_{name}", pe["value"], pe["ci_lo"], pe["ci_hi"])
        if name == "ns" and "cohort" in j:
            table["a16_ns_patients"] = float(j["cohort"]["n_patients"])
    p = RESULTS / "pancancer_tumouronly" / "summary.json"
    if p.exists():
        table["a16_pan_patients"] = float(json.loads(p.read_text())["n_patients"])
    # Session 54: the one outcome comparison the tumour-only NSCLC run gains,
    # quoted as a triple in the audit. Selected by signature AND adjustment,
    # and required to be unique, because each signature has a raw row too.
    p = RESULTS / "nsclc_v3_tumouronly" / "outcome_arm.csv"
    if p.exists():
        hx = [r for r in _rows(p) if r["signature"] == "SIG_HALLMARK_HYPOXIA"
              and r["adjustment"] == "axis_residualised"]
        if len(hx) != 1:
            raise SystemExit(f"{p}: expected one axis_residualised hypoxia "
                             f"row, found {len(hx)}")
        _triple("a16_ns_hypoxia_outcome", float(hx[0]["excess_z"]),
                float(hx[0]["excess_lo"]), float(hx[0]["excess_hi"]))
    # Session 59: the rest of the audit's A16 NSCLC paragraph -- the evidence
    # the author reads for B-18 -- was unbound: the shift, the per-signature
    # comparison, Delta-MAE and the outcome count. The comparison is against
    # the corrected-ordering run the paragraph names (`nsclc_v3_stablesort/`),
    # and was re-derived before binding: against the FROZEN `nsclc_v3/` it
    # gives 11 rises, TGF-beta +0.115 and angiogenesis -0.073, so a recount
    # against the wrong baseline would have "corrected" a right paragraph.
    tu = RESULTS / "nsclc_v3_tumouronly"
    ss = RESULTS / "nsclc_v3_stablesort"
    if (tu / "immune_excess.csv").exists() and (ss / "immune_excess.csv").exists():
        ref = {r["signature"]: float(r["excess_z"]) for r in _rows(ss / "immune_excess.csv")}
        d = {r["signature"]: float(r["excess_z"]) - ref[r["signature"]]
             for r in _rows(tu / "immune_excess.csv")}
        if set(d) != set(ref) or len(d) != 16:
            raise SystemExit("A16: the tumour-only and corrected-ordering runs no "
                             "longer carry the same 16 signatures")
        top = max(d, key=d.get)
        if "TGF_BETA" not in top:
            raise SystemExit(f"A16 CHANGED: the largest per-signature rise is {top}, "
                             "and the audit names TGF-beta")
        table["a16_ns_beats"] = float(_count_beats(tu / "immune_excess.csv"))
        table["a16_ns_sigs_rise"] = float(sum(v > 0 for v in d.values()))
        table["a16_ns_tgf_rise"] = d[top]
        table["a16_ns_angio_fall"] = -d["SIG_HALLMARK_ANGIOGENESIS"]
    if (tu / "summary.json").exists():
        dm = json.loads((tu / "summary.json").read_text())["secondary_delta_mae"]
        # A scalar, NOT a registered triple: registering it made the triple
        # rule identify HPC4's corrected-ordering 0.00481 [0.00199, 0.00748]
        # and the manuscript's 0.0048 [0.002, 0.0075] as this quantity, by
        # point estimate alone, and fail their bounds (measured, session 59).
        table["a16_ns_dmae"] = float(dm["value"])
    for name, run in (("a16_ns_outcome_beat", tu), ("sortaudit_outcome_beat", ss)):
        if (run / "outcome_arm.csv").exists():
            table[name] = float(sum(r["beats_null"].strip().lower() in ("true", "1")
                                    for r in _rows(run / "outcome_arm.csv")))
    if "a16_ns" in table and "sortaudit_isi" in table:
        table["a16_ns_shift"] = table["a16_ns"] - table["sortaudit_isi"]
    if "a16_pan" in table and "pcsortaudit_isi" in table:
        # A FALL, stored positive: the documents write the sign in the text.
        table["a16_pan_fall"] = table["pcsortaudit_isi"] - table["a16_pan"]

    # E6 (session 49): PLAGE, the third scorer, and E8: omega beside alpha.
    p = RESULTS / "scorer_sensitivity_plage" / "summary.json"
    if p.exists():
        j = json.loads(p.read_text())
        for key, arm in (("e6_meanz_unc", "mean_z__uncorrected"),
                         ("e6_plage_unc", "plage__uncorrected")):
            a_ = j["arms"][arm]
            table[key], table[key + "_lo"], table[key + "_hi"] = (
                float(a_["isi"]), float(a_["lo"]), float(a_["hi"]))
            table[key + "_beats"] = float(a_["beats_null"])
        table["e6_plage_recon"] = float(j["reconstructed"]["plage"]["isi"])
        table["e6_plage_recon_beats"] = float(j["reconstructed"]["plage"]["beats_null"])
        # Session 59: the results README quotes the mean-z reconstruction's count too.
        table["e6_meanz_recon_beats"] = float(j["reconstructed"]["mean_z"]["beats_null"])
        if j["arms"]["plage__disattenuated"]["isi"] is not None:
            raise SystemExit("E6 CHANGED: the disattenuated PLAGE arm now returns a value; "
                             "Results and Limitation 9 say it is undefined")
    else:
        missing.append(str(p.relative_to(REPO)))
    p = RESULTS / "omega_reliability" / "summary.json"
    if p.exists():
        j = json.loads(p.read_text())
        if j["alpha_gate_max_abs_diff"] > 1e-10:
            raise SystemExit("E8 BROKEN: its alpha no longer reproduces the registered one")
        for k in ("median_alpha_observed", "median_omega_observed", "median_alpha_null",
                  "median_omega_null", "median_gap_alpha", "median_gap_omega"):
            table["e8_" + k] = float(j[k])
        for key, block in (("e8_recon_alpha", "reconstructed_alpha"),
                           ("e8_recon_omega", "reconstructed_omega")):
            table[key] = float(j[block]["isi"])
            table[key + "_lo"], table[key + "_hi"] = float(j[block]["lo"]), float(j[block]["hi"])
        table["e8_beats_null_omega"] = float(j["beats_null_omega"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # E16 (session 49): the site-classification AUROC on patient subsamples.
    # A fraction below 1 was drawn three times; the aggregate quoted in Results
    # is the MEDIAN over those repeats, computed here rather than stored.
    for cohort, short in (("pancancer", "pan"), ("nsclc", "ns")):
        for arm in ("registered", "small"):
            p = RESULTS / "site_auroc_stress" / f"{cohort}_{arm}.json"
            if not p.exists():
                missing.append(str(p.relative_to(REPO)))
                continue
            for run in json.loads(p.read_text())["runs"]:
                if run["median_auroc"] is None:
                    continue
                tag = f"e16_{short}_{arm}_f{round(run['fraction'] * 100):03d}"
                table.setdefault(tag + "_medians", []).append(float(run["median_auroc"]))
                if run["fraction"] == 0.25:
                    # Session 54: the Results quote the quarter's patient count.
                    n_q = table.setdefault(f"e16_{short}_quarter_n", float(run["n_patients"]))
                    if n_q != float(run["n_patients"]):
                        raise SystemExit(f"{p}: quarter repeats differ in size")
                    # Session 58: "10 to 14 evaluable sites" at a quarter.
                    table.setdefault(f"e16_{short}_{arm}_quarter_sites", []).append(
                        float(run["n_sites_evaluated"]))
                if run["fraction"] == 1.0:
                    table[tag + "_sites"] = float(run["n_sites_evaluated"])
                    table[tag + "_frac_above_0.9"] = float(run["frac_above_0.9"])
    for key in [k for k in table if k.endswith("_quarter_sites")]:
        vals = table.pop(key)
        table[key + "_min"], table[key + "_max"] = min(vals), max(vals)
    for key in [k for k in table if k.endswith("_frac_above_0.9")]:
        table[key[: -len("_frac_above_0.9")] + "_frac_pct"] = 100.0 * table[key]
    for key in [k for k in table if k.endswith("_medians")]:
        vals = sorted(table.pop(key))
        base = key[: -len("_medians")]
        table[base + "_median"] = _median(vals)
        table[base + "_lowest"] = vals[0]
        table[base + "_highest"] = vals[-1]

    # E11 (session 49): the outcome arm's set-level composition test.
    p = RESULTS / "outcome_category_test.json"
    if p.exists():
        j = json.loads(p.read_text())
        table["e11_p_all_inside"] = float(j["p_all_winners_inside_group"])
        table["e11_relabellings"] = float(j["relabellings_enumerated"])
        table["e11_relabellings_at_least"] = (
            float(j["p_excess_difference_one_sided"]) * float(j["relabellings_enumerated"]))
    else:
        missing.append(str(p.relative_to(REPO)))

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
        # Session 59: the ssGSEA reconstruction's count ("-0.1333 0/16").
        table["a7_recon_ssgsea_beats"] = float(j["reconstructed"]["ssgsea"]["beats_null"])
        # Session 58: the mean-z arm's count beside it ("16/16 under mean-z").
        table["a7_beats_meanz_unc"] = float(j["beats_null"]["mean_z__uncorrected"])
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
            # Session 58: the poster table's "Sets" column -- how many curated
            # signatures are large enough to be cut to k.
            table[f"small_panel_k{int(row['k'])}_sets"] = float(row["n_sets"])
        # Session 54: the poster prints the bracket's end-to-end ratio in bold
        # ("5.4x larger at k = 10 than at 160"). Derived from the two stored
        # values above, never typed, so it moves if either does.
        if "small_panel_k10" in table and "small_panel_k160" in table:
            table["small_panel_ratio_k10_k160"] = (
                table["small_panel_k10"] / table["small_panel_k160"])
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

    # ---- added 2026-09-08: the two audit numbers that DECIDE A1's verdict ----
    #
    # `14-SCIENCE-AUDIT.md` argues that the missing BH correction is harmless for
    # the ISI and possibly decisive for the outcome arm, and it rests that
    # argument on exactly two numbers: the smallest pan-cancer `excess_lo` (far
    # from zero, so multiplicity cannot flip it) and hypoxia's outcome-arm
    # `excess_lo` (a hair above zero, so it might). Both were unchecked -- the
    # document that states the audit's conclusions was the third-least covered
    # file in the tree at 10%.
    #
    # They come from two DIFFERENT artefacts, which is why one key would not do:
    # the ISI figure is a minimum over `immune_excess.csv`, the outcome figure is
    # one cell of `outcome_arm.csv` selected by BOTH signature and adjustment.
    # Hypoxia appears twice there (`axis_residualised` 0.006885 and `raw`
    # 0.018480); the audit quotes the residualised arm, and taking the first
    # matching row rather than naming the adjustment would silently check the
    # wrong one whenever the file's row order changed.
    p = RESULTS / "pancancer_v3" / "immune_excess.csv"
    if p.exists():
        table["pancancer_min_excess_lo"] = min(
            float(r["excess_lo"]) for r in _rows(p))
    else:
        missing.append(str(p.relative_to(REPO)))

    p = RESULTS / "pancancer_v3" / "outcome_arm.csv"
    if p.exists():
        hyp = [r for r in _rows(p)
               if r["signature"] == "SIG_HALLMARK_HYPOXIA"
               and r["adjustment"] == "axis_residualised"]
        if len(hyp) == 1:
            table["pancancer_hypoxia_excess_lo"] = float(hyp[0]["excess_lo"])
        else:
            missing.append(f"{p.relative_to(REPO)}:HYPOXIA/axis_residualised "
                           f"matched {len(hyp)} row(s), expected exactly 1")

        # Session 53: the five per-signature excesses the POSTER draws. These
        # are the outcome arm's whole result as a reader at the board sees it,
        # and every one was unchecked -- the poster was the least-covered
        # document in the tree at 10%, and these six numbers are the ones a
        # reviewer would challenge.
        #
        # The adjustment is named for the reason the hypoxia row above names
        # it, and the reason is sharper here: under `raw` the SAME five
        # signatures beat their null but the ORDER changes -- EMT leads at
        # 0.0623 and G2M falls to 0.0520 -- so a row-order-dependent read would
        # not fail loudly, it would quietly check a different arm's numbers
        # against the residualised ones the poster prints.
        for _key, _sig in (
            ("g2m", "SIG_HALLMARK_G2M_CHECKPOINT"),
            ("e2f", "SIG_HALLMARK_E2F_TARGETS"),
            ("angio", "SIG_HALLMARK_ANGIOGENESIS"),
            ("emt", "SIG_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
            ("hypoxia", "SIG_HALLMARK_HYPOXIA"),
        ):
            _hit = [r for r in _rows(p)
                    if r["signature"] == _sig
                    and r["adjustment"] == "axis_residualised"]
            if len(_hit) == 1:
                table[f"pancancer_outcome_excess_{_key}"] = float(
                    _hit[0]["excess_z"])
            else:
                missing.append(
                    f"{p.relative_to(REPO)}:{_sig}/axis_residualised matched "
                    f"{len(_hit)} row(s), expected exactly 1")
    else:
        missing.append(str(p.relative_to(REPO)))

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
        # Session 54: the Results sentence that follows ("Spearman rho = -0.50,
        # p = 0.048 pan-cancer; rho = -0.59, p = 0.015 in NSCLC: it averages
        # +0.141 across the 200-gene sets"), computed the way figure 2 computes
        # it, from the deposited post-filtering panel sizes. "The 200-gene
        # sets" are the eleven Hallmark sets of nominal size 200, which keep
        # 198-200 genes after filtering; the six that keep exactly 200 average
        # 0.152 pan-cancer, so the method is recorded here, not guessed.
        gp = RESULTS / "gene_set_sizes.csv"
        if gp.exists():
            from scipy import stats as sps
            sizes = {r["signature"]: int(r["n_genes"]) for r in _rows(gp)}
            ks = [sizes[r["signature"].removeprefix("SIG_HALLMARK_")] for r in rows]
            rho, pv = sps.spearmanr(ks, gap_resid)
            table[f"{prefix}_relgap_size_rho"] = float(rho)
            table[f"{prefix}_relgap_size_p"] = float(pv)
            # Session 58: the magnitude, for text that prints the sign apart.
            table[f"{prefix}_relgap_size_rho_abs"] = abs(float(rho))
            full = [g for k, g in zip(ks, gap_resid) if k >= 198]
            if len(full) != 11:
                raise SystemExit(f"{gp}: expected 11 nominal-200 sets, found {len(full)}")
            table[f"{prefix}_rel_gap_fullsize_mean"] = sum(full) / len(full)
            # Session 59 (B-34, applied 2026-09-23): the sentence now names the
            # eleven and their post-filtering range, "(198-200 genes after
            # filtering)"; both ends are bound, and the eleven is asserted above.
            fs = [k for k in ks if k >= 198]
            table["fullsize_set_min_genes"] = float(min(fs))
            table["fullsize_set_max_genes"] = float(max(fs))
            # Session 58: "the smallest, 36-gene angiogenesis set". Asserted to
            # BE the smallest, so the word is checked along with the number.
            if min(sizes, key=sizes.get) != "ANGIOGENESIS":
                raise SystemExit(f"{gp}: angiogenesis is no longer the smallest set")
            table["gene_set_size_angiogenesis"] = float(sizes["ANGIOGENESIS"])
        else:
            missing.append(str(gp.relative_to(REPO)))

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
        # Session 54: the poster states the same median as a fraction
        # ("0.055 given cancer type"), which the percent key cannot match.
        table["label_site_given_type_r2"] = (
            table["label_site_given_type_pct"] / 100.0)
        # E14: what the NSCLC/pan-cancer gap is made of (Results).
        table["lsv_pan_adj_pct"] = 100.0 * _median(float(r["adj_r2_site_only"]) for r in lsv)
        table["lsv_pan_type_pct"] = 100.0 * _median(float(r["r2_type_only"]) for r in lsv)
        table["label_plate_within_site"] = _median(
            float(r["r2_plate_within_site"]) for r in lsv)
    else:
        missing.append(str(p.relative_to(REPO)))
    p = RESULTS / "nsclc_v3_stablesort" / "label_site_variance.csv"
    if p.exists():
        lsv = _rows(p)
        table["lsv_nsclc_adj_pct"] = 100.0 * _median(float(r["adj_r2_site_only"]) for r in lsv)
        table["lsv_nsclc_type_pct"] = 100.0 * _median(float(r["r2_type_only"]) for r in lsv)
        table["lsv_nsclc_site_given_type_pct"] = 100.0 * _median(
            float(r["r2_site_given_type"]) for r in lsv)
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
    # Session 54: "only 0.6-1.0% of sites" is the range, across the three
    # comparison groups, of the share of sites where that group appears that
    # carry >= 10 patients of it and of the reference. A recount over all 619
    # sites gives 0.5-1.0%, so the denominator is recorded, not assumed.
    p = RESULTS / "ancestry" / "site_identifiability.csv"
    if p.exists():
        fr = [100.0 * float(r["frac_sites_supporting"]) for r in _rows(p)]
        table["ancestry_site_frac_min_pct"] = min(fr)
        table["ancestry_site_frac_max_pct"] = max(fr)
    else:
        missing.append(str(p.relative_to(REPO)))

    # Session 50 (ledger F7.4): numbers the supplementary CAPTIONS state, read
    # from the same source files the tables are built from. Before this,
    # CAPTIONS.md had 1 of its 58 numeric literals checked.
    for tag, rel in (("pan", "pancancer_v3/label_site_variance.csv"),
                     ("nsclc", "nsclc_v3_stablesort/label_site_variance.csv")):
        p = RESULTS / rel
        if p.exists():
            rows = _rows(p)
            ns = {r["n"] for r in rows}
            sites = {r["n_sites"] for r in rows}
            if len(ns) != 1 or len(sites) != 1:
                raise SystemExit(f"{rel}: n or n_sites varies across signatures")
            table[f"cap_{tag}_lsv_n"] = float(ns.pop())
            table[f"cap_{tag}_lsv_sites"] = float(sites.pop())
        else:
            missing.append(str(p.relative_to(REPO)))
    p = RESULTS / "nsclc_v3" / "site_control_combat.csv"
    if p.exists():
        table["cap_combat_auroc_median"] = _median([float(r["auroc"]) for r in _rows(p)])
    # Session 57: the deposited README's demo table ("From `results/demo/`"),
    # which a reader who runs `02_run_audit.py --demo` compares against.
    p = RESULTS / "demo" / "immune_excess.csv"
    if p.exists():
        for r in _rows(p):
            name = r["signature"].removeprefix("SIG_immune_")
            for col, short in DEMO_TABLE_COLUMNS:
                table[f"demo_{name}_{short}"] = float(r[col])
    else:
        missing.append(str(p.relative_to(REPO)))
    p = RESULTS / "pancancer_v3" / "decomposition.csv"
    if p.exists():
        emt = [r for r in _rows(p)
               if r["signature"] == "SIG_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
        if len(emt) != 1:
            raise SystemExit("pancancer_v3/decomposition.csv: expected one EMT row")
        r = emt[0]
        base = r.get("M4_+stage") or r["M3_+site"]
        table["cap_emt_incr_unadj"] = float(r["Mfull_+image"]) - float(base)
        table["cap_emt_incr_adj"] = float(r["incremental_image_r2"])
    else:
        missing.append(str(p.relative_to(REPO)))

    # Session 50 (ledger F1.1): how much of each partition the A9 fix changed,
    # and every yardstick ratio Limitation 8 quotes for it. Written by
    # `34_sort_fix_fold_change.py` from the stored per-patient folds and
    # summaries. Before this, Limitation 8's 2.4/2.9/3.4/0.87/0.79 ratios had no
    # authority at all, and two of them were wrong at one decimal (2.4 for
    # 2.457; 3.4 for 3.347).
    p = RESULTS / "sort_fix_fold_change.json"
    if p.exists():
        j = json.loads(p.read_text())
        for c, tag in (("nsclc", "ns"), ("pancancer", "pc")):
            r = j[c]
            table[f"sffc_{tag}_sites_tied"] = float(r["n_sites_tied"])
            table[f"sffc_{tag}_pct_tied"] = float(r["pct_patients_in_tied_sites"])
            table[f"sffc_{tag}_moved"] = float(r["n_patients_fold_changed"])
            table[f"sffc_{tag}_pct_moved"] = float(r["pct_patients_fold_changed"])
            table[f"sffc_{tag}_ari"] = float(r["adjusted_rand_index"])
            table[f"sffc_{tag}_move_sds_small"] = float(r["move_in_small_sweep_sds"])
            table[f"sffc_{tag}_move_sds_24"] = float(r["move_in_24_partition_sds"])
            table[f"sffc_{tag}_n_pairs"] = float(r["n_partition_pairs_24"])
            table[f"sffc_{tag}_pairs_ge_move"] = float(
                r["n_pairs_differing_by_at_least_the_move"])
        table["sffc_gap_sds_small"] = float(j["nsclc_platform_gap"]["gap_in_small_sweep_sds"])
        table["sffc_gap_sds_24"] = float(j["nsclc_platform_gap"]["gap_in_24_partition_sds"])
        table["sffc_contrast_small"] = float(j["contrast"]["nsclc_over_pancancer_in_small_sweep_sds"])
        table["sffc_contrast_24"] = float(j["contrast"]["nsclc_over_pancancer_in_24_partition_sds"])
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
            # Session 58: the task list quotes it as a percentage.
            table[f"sortaudit_{tag}_frac_sd_pct"] = 100.0 * table[f"sortaudit_{tag}_frac_sd"]
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

    # A12's two-platform record (ledger C9), asserted the same two-sided way:
    # `26_sort_audit.py --library-order` on each machine. The STABLE GroupKFold
    # fold assignment must agree across machines and the LIBRARY one must not
    # -- if the library ones agreed, the A12 write-up would overstate what was
    # measured; if the stable ones disagreed, pinning would not fix it.
    gko = {}
    for tag, fname in (("macos", "groupkfold_order_darwin_arm64.json"),
                       ("hpc4", "groupkfold_order_linux_x86_64_login.json")):
        p = RESULTS / fname
        if p.exists():
            gko[tag] = json.loads(p.read_text())
        else:
            missing.append(str(p.relative_to(REPO)))
    if len(gko) == 2:
        if gko["macos"]["fold_hash_stable"] != gko["hpc4"]["fold_hash_stable"]:
            raise SystemExit(
                "A12 BROKEN: the STABLE GroupKFold fold assignment differs across "
                f"platforms ({gko['macos']['fold_hash_stable']} vs "
                f"{gko['hpc4']['fold_hash_stable']}), so pinning would not fix it.")
        if gko["macos"]["fold_hash_library"] == gko["hpc4"]["fold_hash_library"]:
            raise SystemExit(
                "A12 UNSUPPORTED: the LIBRARY GroupKFold fold assignments agree "
                f"across platforms (both {gko['macos']['fold_hash_library']}); "
                "check that the two files really came from different machines.")

    # Session 60 (H89): results the science audit states that no key held until
    # its unchecked literals were swept. Every source below is a hashed run
    # (the null-reliability records were hashed for this, PROVENANCE 233 ->
    # 236). Derived quantities -- widths, shifts, ratios -- are computed from
    # the same stored values the audit's tables print beside them.
    try:
        sens = RESULTS / "nsclc_sensitivity"
        ref_ie, glob_ie = sens / "reference" / "immune_excess.csv", sens / "global_axis" / "immune_excess.csv"
        unc = "excess_z_UNCORRECTED_do_not_report"
        table["sens_ns_ref_curated_r"] = _med(ref_ie, "r_residual_refit")
        table["sens_ns_ref_null_r"] = _med(ref_ie, "null_mean_r")
        table["sens_ns_ref_alpha_null"] = _med(ref_ie, "alpha_null_mean")
        table["sens_ns_ref_unc_excess"] = _med(ref_ie, unc)
        table["sens_ns_global_alpha_null"] = _med(glob_ie, "alpha_null_mean")
        table["sens_ns_global_unc_excess_abs"] = abs(_med(glob_ie, unc))
        degen = [r for r in _rows(glob_ie) if float(r["degenerate"])]
        table["sens_ns_degen_null_r_min"] = min(float(r["null_mean_r"]) for r in degen)
        table["sens_ns_degen_null_r_max"] = max(float(r["null_mean_r"]) for r in degen)
        table["sens_ns_degen_null_sd_min"] = min(float(r["null_sd_r"]) for r in degen)
        table["sens_ns_degen_null_sd_max"] = max(float(r["null_sd_r"]) for r in degen)
        ref = ("sens_ns_reference", "sens_ns_reference_lo", "sens_ns_reference_hi")
        table["sens_ns_k3_shift_abs"] = abs(table["sens_ns_k3"] - table[ref[0]])
        table["sens_ns_k10_shift"] = table["sens_ns_k10"] - table[ref[0]]
        for arm in ("null_boot_400", "null_boot_1000"):
            table[f"sens_ns_{arm}_move"] = max(
                abs(table[f"sens_ns_{arm}_lo"] - table[ref[1]]),
                abs(table[f"sens_ns_{arm}_hi"] - table[ref[2]]))
        pref_ie = RESULTS / "pancancer_v3_stablesort" / "immune_excess.csv"
        pglob_ie = RESULTS / "pancancer_sensitivity" / "global_axis" / "immune_excess.csv"
        for tag, ie in (("ref", pref_ie), ("global", pglob_ie)):
            table[f"sens_pan_{tag}_null_sd"] = _med(ie, "null_sd_r")
            table[f"sens_pan_{tag}_alpha_null"] = _med(ie, "alpha_null_mean")
            table[f"sens_pan_{tag}_unc_excess"] = _med(ie, unc)
        pv = {v["variant"]: v for v in json.loads(
            (RESULTS / "pancancer_sensitivity" / "summary.json").read_text())["variants"]}
        table["sens_pan_global_axis_var_pct"] = 100.0 * float(pv["global_axis"]["axis_variance_explained"])
        table["sens_pan_global_shift"] = float(pv["global_axis"]["shift_from_reference"])
        table["sens_pan_reference_width"] = table["sens_pan_reference_hi"] - table["sens_pan_reference_lo"]
        table["sens_pan_global_axis_width"] = table["sens_pan_global_axis_hi"] - table["sens_pan_global_axis_lo"]
        # The A12 memo's two NSCLC runs, both hashed: corrected ordering, and pinned.
        ss, pn = RESULTS / "nsclc_v3_stablesort", RESULTS / "nsclc_v3_pinned"
        table["sortaudit_dmae"] = float(json.loads((ss / "summary.json").read_text())
                                        ["secondary_delta_mae"]["value"])
        dr_ss = {r["signature"]: float(r["delta_r"]) for r in _rows(ss / "per_signature.csv")}
        dr_pn = {r["signature"]: float(r["delta_r"]) for r in _rows(pn / "per_signature.csv")}
        assert dr_ss.keys() == dr_pn.keys()
        table["sortaudit_mean_dr"] = sum(dr_ss.values()) / len(dr_ss)
        table["pinned_nsclc_mean_dr"] = sum(dr_pn.values()) / len(dr_pn)
        diffs = [abs(dr_pn[k] - dr_ss[k]) for k in dr_ss]
        table["pinned_vs_ss_dr_max"] = max(diffs)
        table["pinned_vs_ss_dr_mean"] = sum(diffs) / len(diffs)
        # A5's partition tables: each partition's ISI, the intervals' widths and
        # the observed/predicted ratios printed beside them.
        for coh, run in (("nsclc", "partition_variance_nsclc"), ("pancancer", "partition_variance_pancancer")):
            for r in _rows(RESULTS / run / "per_partition.csv"):
                table[f"pv_{coh}_seed{int(r['split_seed'])}"] = float(r["isi"])
            for kind in ("reported", "honest"):
                table[f"{coh}_{kind}_ci_width"] = (table[f"{coh}_{kind}_ci_hi"]
                                                   - table[f"{coh}_{kind}_ci_lo"])
        sqrt_n = (table["pancancer_n"] / table["nsclc_n"]) ** 0.5   # = cohort_sqrt_n_ratio, below
        table["bootstrap_ratio_over_sqrt_n"] = table["bootstrap_se_ratio"] / sqrt_n
        table["partition_ratio_over_sqrt_n"] = table["partition_sd_ratio"] / sqrt_n
        table["partition_ratio_over_sqrt_sites"] = table["partition_sd_ratio"] / table["pv24_sqrt_sites"]
        table["pv24_bootstrap_over_sqrt_n"] = table["pv24_bootstrap_ratio"] / table["pv24_sqrt_patients"]
        # Control C's A4b probe and E14's wider frame, from the calibration records.
        cc = {c: json.loads((RESULTS / "control_c_calibration" / f"{c}.json").read_text())
              for c in ("nsclc", "pancancer")}
        probe = cc["pancancer"]["per_signature"][cc["pancancer"]["probe_signature"]]
        table["control_c_probe_observed_r2"] = float(probe["observed_r2"])
        table["control_c_probe_perm_median_r2"] = float(probe["perm_median_r2"])
        table["control_c_probe_perm_p"] = float(probe["perm_p"])
        for c in ("nsclc", "pancancer"):
            table[f"e15_{c}_median_observed_r2"] = float(cc[c]["median_observed_r2"])
        # E14's table: the label-side decomposition's medians, both cohorts.
        for c, run in (("pancancer", "pancancer_v3"), ("nsclc", "nsclc_v3_stablesort")):
            lsv = _rows(RESULTS / run / "label_site_variance.csv")
            for col in ("adj_r2_site_only", "r2_type_only", "r2_site_given_type"):
                table[f"lsv_{c}_{col}"] = _median(float(r[col]) for r in lsv)
        # Third wave: four more derived or stored values the audit prints.
        table["sortaudit_macos_elements"] = (table["sortaudit_macos_samples"]
                                             * table["sortaudit_macos_genes"])
        table["lsv_crude_gap"] = table["label_site_variance_r2"] - table["nsclc_label_site_variance_r2"]
        oc = json.loads((RESULTS / "outcome_category_test.json").read_text())
        table["e11_excess_difference"] = float(oc["mean_excess_difference_group_minus_rest"])
        pl = json.loads((RESULTS / "scorer_sensitivity_plage" / "summary.json").read_text())
        table["e6_meanz_recon"] = float(pl["reconstructed"]["mean_z"]["isi"])
        # E9: the null-reliability records.
        for c in ("nsclc", "pancancer"):
            nr = json.loads((RESULTS / "null_reliability" / f"{c}.json").read_text())
            ps = nr["per_signature"].values()
            table[f"e9_{c}_alpha_min"] = float(nr["overall_alpha_null_min"])
            table[f"e9_{c}_p01_min"] = min(float(v["alpha_null_p01"]) for v in ps)
            table[f"e9_{c}_ratio_max"] = float(nr["overall_max_abs_r_over_sqrt_alpha"])
            table[f"e9_{c}_floor_multiple"] = float(nr["overall_alpha_null_min"]) / 0.01
    except (OSError, KeyError, ValueError) as exc:
        missing.append(f"session-60 audit authorities: {exc!r}")

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

    # Session 51: the dry run's cut ladder and AI-sentence variants. Pure
    # Python over two files, so it belongs in the fast table too.
    table.update(_abstract_variants())

    # Session 52: the cover letter's own length, for the same reason -- it is a
    # submission-day number the author acts on, and it had drifted ~30%.
    table.update(_cover_letter_words())

    # Session 53: the manuscript abstract's own counts and their sum, for the
    # same reason. B-30 edits a line inside the measured span, so these three
    # would go stale the moment the author chooses a wording.
    table.update(_paper_abstract_words())

    # Session 53: the manuscript's body and figure counts, which AACR asks the
    # cover letter to state and which had been `[X]` and `[N]` placeholders.
    table.update(_paper_body_words())

    if slow:
        n = _test_count()
        if n is not None:
            table["test_count"] = float(n)
        counts = _abstract_counts()
        ab = None if counts is None else (counts["total"], counts["limit"])
        if counts is not None:
            # Session 51: the other figures the counter prints, from the same
            # run. Live-measured, so they sit beside abstract_chars here and in
            # the fast table's expected absences.
            for key, name in (("title", "abstract_title_chars"),
                              ("naive", "abstract_naive_len"),
                              ("words", "abstract_body_words")):
                if key in counts:
                    table[name] = float(counts[key])
        if ab is not None:
            table["abstract_chars"], table["abstract_limit"] = float(ab[0]), float(ab[1])
            # DERIVED, not stored: the two counts were checked and the
            # subtraction between them was not, so 07-ABSTRACT-DRAFT.md read
            # "2,515 / 2,600 with 83 characters of headroom" -- 2,515 correct,
            # 2,600 correct, and the difference between them wrong by two, on
            # the line that says how much room is left before a hard deadline.
            # Found 2026-09-07 by running the counter, not by reading.
            table["abstract_headroom"] = float(ab[1] - ab[0])
            # Session 58: the dry run's "(97%)" of the limit and its naive
            # len()'s "16% higher" -- two-digit, so the coverage list dropped
            # them as structural and a mutation sweep moved both unnoticed.
            table["abstract_pct_of_limit"] = 100.0 * ab[0] / ab[1]
            if "abstract_naive_len" in table:
                table["abstract_naive_excess_pct"] = (
                    100.0 * (table["abstract_naive_len"] / ab[0] - 1.0))
        sc = _site_counts()
        if sc is not None:
            table.update(sc)
        a16x = _a16_extent()
        if a16x is not None:
            table.update(a16x)
        snap = _snapshot_size()
        if snap is not None:
            table["snapshot_files"], table["snapshot_bytes"] = (
                float(snap[0]), float(snap[1]))
    # Session 54: ratios and transforms the Results print, derived from values
    # read above so they move when those do.
    import math
    if "pancancer_n" in table and "nsclc_n" in table:
        table["cohort_n_ratio"] = table["pancancer_n"] / table["nsclc_n"]
        table["cohort_sqrt_n_ratio"] = math.sqrt(table["cohort_n_ratio"])
    for side in ("null", "obs"):
        if f"a7_ssgsea_rel_{side}" in table:
            table[f"a7_ssgsea_{side}_inflation"] = 1.0 / math.sqrt(
                table[f"a7_ssgsea_rel_{side}"])
    if "nsclc_label_site_variance_r2" in table:
        table["nsclc_label_site_variance_pct"] = (
            100.0 * table["nsclc_label_site_variance_r2"])
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
    # RE-SCOPED 2026-09-16. These were anchored on "corrected ordering", which
    # was unique to the NSCLC paragraph until the pan-cancer re-run landed and
    # its paragraph used the same phrase -- at which point both guards started
    # checking pan-cancer numbers against NSCLC authorities and went red. The
    # anchors are now the sentence that opens each paragraph, and `WRAP:`
    # rather than `PARA:` because both cross a line wrap.
    ("sortaudit_isi", r"ISI of\s+\**(\d\.\d+)", 5e-5,
     "NSCLC ISI under the corrected sorts",
     "WRAP:That separately reported re-run now exists"),
    ("sortaudit_shift", r"moves by\s+\**(\d\.\d+)", 5e-5,
     "how far the corrected sorts move the NSCLC point estimate",
     "WRAP:That separately reported re-run now exists"),
    ("pcsortaudit_isi", r"ISI of\s+\**(\d\.\d+)", 5e-5,
     "pan-cancer ISI under the corrected sorts",
     "WRAP:equivalent re-run for the pan-cancer cohort"),
    ("pcsortaudit_shift", r"moves by\s+\**(\d\.\d+)", 5e-5,
     "how far the corrected sorts move the pan-cancer point estimate",
     "WRAP:equivalent re-run for the pan-cancer cohort"),
    # S1's DERIVED NSCLC multiplicity values. These are the project's first
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
    # Floors 5e-7 since 2026-09-24, so the PRINTED precision governs (a floor
    # is a minimum tolerance, never a ceiling): at 5e-5 a six-decimal quote of
    # the 24-partition run's SE passed against this six-partition authority.
    # The 24-partition rows are excluded here and bound to their own `pv24_`
    # keys below.
    ("bootstrap_se",
     r"bootstrap (?:SE|standard error)[^\d\n]{0,12}(\d+\.\d+)", 5e-7,
     "NSCLC bootstrap SE", "PARA:NSCLC!!pan-cancer|_stable24"),
    ("pancancer_bootstrap_se",
     r"bootstrap (?:SE|standard error)[^\d\n]{0,12}(\d+\.\d+)", 5e-7,
     "pan-cancer bootstrap SE", "PARA:pan-cancer!!NSCLC|_stable24"),
    ("pv24_nsclc_bootstrap_se",
     r"nsclc_stable24/` \|.*bootstrap SE (\d\.\d{4}), combined \d\.\d{4}", 5e-7,
     "24-partition NSCLC bootstrap SE (results README)", None),
    ("pv24_nsclc_combined_se",
     r"nsclc_stable24/` \|.*bootstrap SE \d\.\d{4}, combined (\d\.\d{4})", 5e-7,
     "24-partition NSCLC combined SE (results README)", None),
    ("pv24_pancancer_bootstrap_se",
     r"stable24/` \|.*Bootstrap SE (\d\.\d{6}) and combined", 5e-7,
     "24-partition pan-cancer bootstrap SE (results README)", None),
    ("pv24_pancancer_combined_se",
     r"stable24/` \|.*Bootstrap SE \d\.\d{6} and combined (\d\.\d{6})", 5e-7,
     "24-partition pan-cancer combined SE (results README)", None),
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
    # Session 49 sensitivity analyses: what the triples do not carry.
    ("sens_ns_global_axis_beats", r"0\.0466\], with\s+(\d+) of 16 signatures above", 0.5,
     "E1 NSCLC signatures above null", "WRAP:"),
    ("sens_ns_global_null_r", r"curated ones\s+\(median r (\d\.\d{3}) against", 5e-4,
     "E1 NSCLC median random-set r", "WRAP:"),
    ("sens_ns_global_curated_r", r"curated ones\s+\(median r \d\.\d{3} against (\d\.\d{3})\)", 5e-4,
     "E1 NSCLC median curated r", "WRAP:"),
    ("sens_ns_global_degenerate", r"and (\d+) of the 16\s+nulls become degenerate", 0.5,
     "E1 NSCLC degenerate nulls", "WRAP:"),
    ("sens_pan_global_axis_beats", r"0\.3189\], with (\d+) of 16\s+signatures", 0.5,
     "E1 pan-TCGA signatures above null", "WRAP:"),
    ("sens_pan_global_curated_r", r"more predictable \(median r (\d\.\d{3}) and", 5e-4,
     "E1 pan-TCGA median curated r, global axis", "WRAP:"),
    ("sens_pan_global_null_r", r"\(median r \d\.\d{3} and (\d\.\d{3}), against", 5e-4,
     "E1 pan-TCGA median random-set r, global axis", "WRAP:"),
    ("sens_pan_ref_curated_r", r"\d\.\d{3}, against (\d\.\d{3}) and \d\.\d{3}\)", 5e-4,
     "pan-TCGA median curated r, registered", "WRAP:"),
    ("sens_pan_ref_null_r", r"\d\.\d{3}, against \d\.\d{3} and (\d\.\d{3})\)", 5e-4,
     "pan-TCGA median random-set r, registered", "WRAP:"),
    ("sens_ns_unmatched_shift", r"matching lowers the NSCLC index by (\d\.\d{3})", 5e-4,
     "E3 shift", "WRAP:"),
    ("sens_ns_null_boot_400_lo", r"interval nearly so: \[(\d\.\d{4}),", 5e-5,
     "E4 400 draws, CI low", "WRAP:"),
    ("sens_ns_null_boot_400_hi", r"interval nearly so: \[\d\.\d{4}, (\d\.\d{4})\]", 5e-5,
     "E4 400 draws, CI high", "WRAP:"),
    ("sens_ns_null_boot_1000_lo", r"\] and \[(\d\.\d{4}), \d\.\d{4}\]\.", 5e-5,
     "E4 1,000 draws, CI low", "WRAP:interval nearly so"),
    ("sens_ns_null_boot_1000_hi", r"\] and \[\d\.\d{4}, (\d\.\d{4})\]\.", 5e-5,
     "E4 1,000 draws, CI high", "WRAP:interval nearly so"),
    ("sens_ns_k3_beats", r"0\.2722\],\s+with (\d+) of 16 signatures", 0.5,
     "E5 3 folds, signatures above null", "WRAP:"),
    ("a16_pan_patients", r"0\.3149\] \(([\d,]+) patients\)", 0.5,
     "A16 pan-TCGA tumor-only patients", "WRAP:"),
    # E6 and E8: what their triples do not carry.
    ("e6_plage_unc_beats", r"0\.1834\] with (\d+)\s*\n?\s*of 16 signatures above their null", 0.5,
     "E6 PLAGE signatures above null", "WRAP:"),
    ("e6_meanz_unc_beats", r"\[0\.2442, 0\.3514\] and (\d+) of 16", 0.5,
     "E6 mean-z control signatures above null", "WRAP:"),
    ("e6_plage_recon", r"PLAGE stays positive\s*\n?\s*\(\+(\d\.\d{4})", 5e-5,
     "E6 PLAGE reconstruction", "WRAP:"),
    ("e6_plage_recon_beats", r"\(\+\d\.\d{4}, (\d+) of 16 above null\)", 0.5,
     "E6 PLAGE reconstruction, signatures above null", "WRAP:"),
    ("e8_median_omega_observed", r"curated sets almost unchanged\s+\(median ω (\d\.\d{3}) against α", 5e-4,
     "E8 median omega, curated", "WRAP:"),
    ("e8_median_alpha_observed", r"curated sets almost unchanged\s+\(median ω \d\.\d{3} against α (\d\.\d{3})\)", 5e-4,
     "E8 median alpha, curated", "WRAP:"),
    ("e8_median_omega_null", r"random sets \(median ω (\d\.\d{3}) against", 5e-4,
     "E8 median omega, random sets", "WRAP:"),
    ("e8_median_alpha_null", r"random sets \(median ω \d\.\d{3} against α (\d\.\d{3})\)", 5e-4,
     "E8 median alpha, random sets", "WRAP:"),
    ("e8_median_gap_alpha", r"gap widens from (\d\.\d{3}) to", 5e-4,
     "E8 reliability gap under alpha", "WRAP:"),
    ("e8_median_gap_omega", r"gap widens from \d\.\d{3} to (\d\.\d{3})", 5e-4,
     "E8 reliability gap under omega", "WRAP:"),
    ("e8_beats_null_omega", r"with α on the same estimator, and (\d+) of 16 signatures", 0.5,
     "E8 signatures above null under omega", "WRAP:"),
    # E16: the subsample AUROC figures quoted in Negative controls.
    ("e16_pan_registered_f100_median", r"median stays at (\d\.\d{3}) at the full 7,168", 5e-4,
     "E16 pan-cancer median AUROC, full cohort", "WRAP:"),
    ("e16_pan_registered_f050_median", r"at the full 7,168 patients,\s+(\d\.\d{3}) at half", 5e-4,
     "E16 pan-cancer median AUROC, half", "WRAP:"),
    ("e16_pan_registered_f025_median", r"at half and (\d\.\d{3}) at a quarter", 5e-4,
     "E16 pan-cancer median AUROC, quarter", "WRAP:"),
    ("e16_ns_registered_f100_median", r"in NSCLC, (\d\.\d{3}), \d\.\d{3} and", 5e-4,
     "E16 NSCLC median AUROC, full cohort", "WRAP:"),
    ("e16_ns_registered_f050_median", r"in NSCLC, \d\.\d{3}, (\d\.\d{3}) and", 5e-4,
     "E16 NSCLC median AUROC, half", "WRAP:"),
    ("e16_ns_registered_f025_lowest", r"and (\d\.\d{3}) to \d\.\d{3} over the same steps", 5e-4,
     "E16 NSCLC median AUROC, quarter, lowest repeat", "WRAP:"),
    ("e16_ns_registered_f025_highest", r"and \d\.\d{3} to (\d\.\d{3}) over the same steps", 5e-4,
     "E16 NSCLC median AUROC, quarter, highest repeat", "WRAP:"),
    ("e16_ns_small_f100_sites", r"evaluates (\d+) NSCLC sites", 0.5,
     "E16 NSCLC sites at the 5-patient threshold", "WRAP:"),
    ("e16_ns_small_f100_median", r"NSCLC sites \(median (\d\.\d{3}), \d+%", 5e-4,
     "E16 NSCLC median AUROC at the 5-patient threshold", "WRAP:"),
    ("e16_ns_small_f100_frac_pct", r"NSCLC sites \(median \d\.\d{3}, (\d+)% above", 0.5,
     "E16 NSCLC percent above 0.9 at the 5-patient threshold", "WRAP:"),
    ("e16_pan_small_f100_median", r"pan-cancer sites \(median (\d\.\d{3}), \d+%", 5e-4,
     "E16 pan-cancer median AUROC at the 5-patient threshold", "WRAP:"),
    ("e16_pan_small_f100_frac_pct", r"pan-cancer sites \(median \d\.\d{3}, (\d+)% above", 0.5,
     "E16 pan-cancer percent above 0.9 at the 5-patient threshold", "WRAP:"),
    ("e16_pan_small_f100_sites", r"and (\d+) pan-cancer sites \(median", 0.5,
     "E16 pan-cancer sites at the 5-patient threshold", "WRAP:"),
    # E11: the set-level composition of the pan-cancer outcome winners.
    ("e11_p_all_inside", r"all land there with\s+probability (\d\.\d+)", 5e-5,
     "P(all five outcome winners inside the six-signature group)", "WRAP:"),
    ("e11_relabellings_at_least", r"these six have \((\d+) of [\d,]+\)", 0.5,
     "relabellings at least as extreme", "WRAP:"),
    ("e11_relabellings", r"these six have \(\d+ of ([\d,]+)\)", 0.5,
     "relabellings enumerated, C(16, 6)", "WRAP:"),
    # B9: the 24-partition sweeps. Each number has an anchor of its own in the
    # Results paragraph, none of which the frozen patterns above can reach.
    ("pv24_nsclc_sd",
     r"In NSCLC the standard deviation across the 24\s+partitions is (\d\.\d+)",
     5e-5, "sd across 24 NSCLC partitions", "WRAP:"),
    ("pv24_nsclc_sd_lo",
     r"24\s+partitions is \d\.\d+, with a\s+chi-square 95% interval of \[(\d\.\d+),",
     5e-5, "chi-square lower bound, NSCLC 24-partition sd", "WRAP:"),
    ("pv24_nsclc_sd_hi",
     r"24\s+partitions is \d\.\d+, with a\s+chi-square 95% interval of \[\d\.\d+,\s+(\d\.\d+)\]",
     5e-5, "chi-square upper bound, NSCLC 24-partition sd", "WRAP:"),
    ("pv24_nsclc_understated_pct",
     r"the NSCLC interval understates total uncertainty by\s+(\d+\.\d+)%",
     5e-2, "NSCLC 24-partition understatement, percent", "WRAP:"),
    ("pv24_nsclc_share_pct",
     r"makes\s+up\s+(\d+\.\d+)% \[[\d.]+%,\s+[\d.]+%\]\s+of\s+the\s+NSCLC total",
     5e-2, "NSCLC 24-partition share of variance", "WRAP:"),
    ("pv24_nsclc_share_lo_pct",
     r"makes\s+up\s+[\d.]+% \[(\d+\.\d+)%,\s+[\d.]+%\]\s+of\s+the\s+NSCLC total",
     5e-2, "NSCLC 24-partition share, lower bound", "WRAP:"),
    ("pv24_nsclc_share_hi_pct",
     r"makes\s+up\s+[\d.]+% \[[\d.]+%,\s+(\d+\.\d+)%\]\s+of\s+the\s+NSCLC total",
     5e-2, "NSCLC 24-partition share, upper bound", "WRAP:"),
    ("pv24_pancancer_sd",
     r"Pan-cancer, the\s+standard deviation across the 24 partitions is (\d\.\d+)",
     5e-5, "sd across 24 pan-cancer partitions", "WRAP:"),
    ("pv24_pancancer_sd_lo",
     r"across the 24 partitions is \d\.\d+ \[(\d\.\d+),",
     5e-5, "chi-square lower bound, pan-cancer 24-partition sd", "WRAP:"),
    ("pv24_pancancer_sd_hi",
     r"across the 24 partitions is \d\.\d+ \[\d\.\d+,\s+(\d\.\d+)\]",
     5e-5, "chi-square upper bound, pan-cancer 24-partition sd", "WRAP:"),
    ("pv24_pancancer_understated_pct",
     r"pan-cancer interval understates total uncertainty by\s+(\d+\.\d+)%",
     5e-2, "pan-cancer 24-partition understatement, percent", "WRAP:"),
    ("pv24_pancancer_share_pct",
     r"makes\s+up\s+(\d+\.\d+)% \[[\d.]+%,\s+[\d.]+%\]\s+of\s+the\s+pan-cancer total",
     5e-2, "pan-cancer 24-partition share of variance", "WRAP:"),
    ("pv24_pancancer_share_lo_pct",
     r"makes\s+up\s+[\d.]+% \[(\d+\.\d+)%,\s+[\d.]+%\]\s+of\s+the\s+pan-cancer total",
     5e-2, "pan-cancer 24-partition share, lower bound", "WRAP:"),
    ("pv24_pancancer_share_hi_pct",
     r"makes\s+up\s+[\d.]+% \[[\d.]+%,\s+(\d+\.\d+)%\]\s+of\s+the\s+pan-cancer total",
     5e-2, "pan-cancer 24-partition share, upper bound", "WRAP:"),
    # ...and the poster's one-line version of the same (panel 4).
    ("pv24_pancancer_sd", r"the sds agree \((\d\.\d+), \d\.\d+\)", 5e-5,
     "sd across 24 pan-cancer partitions (poster)", "WRAP:"),
    ("pv24_nsclc_sd", r"the sds agree \(\d\.\d+, (\d\.\d+)\)", 5e-5,
     "sd across 24 NSCLC partitions (poster)", "WRAP:"),
    ("pv24_pancancer_share_pct", r"shares\s+separate, (\d+\.\d+)% vs", 5e-2,
     "pan-cancer 24-partition share (poster)", "WRAP:"),
    ("pv24_nsclc_share_pct", r"shares\s+separate, \d+\.\d+% vs (\d+\.\d+)%", 5e-2,
     "NSCLC 24-partition share (poster)", "WRAP:"),
    ("pv24_bootstrap_ratio",
     r"falls with cohort size, by (\d\.\d+) against the", 5e-4,
     "bootstrap SE ratio between cohorts, 24 partitions", "WRAP:"),
    ("pv24_sqrt_patients",
     r"falls with cohort size, by \d\.\d+ against the (\d\.\d+) predicted", 5e-4,
     "sqrt(patients) prediction", "WRAP:"),
    ("pv24_partition_ratio",
     r"its ratio\s+between the cohorts is (\d\.\d+), with an F-based", 5e-4,
     "partition sd ratio between cohorts, 24 partitions", "WRAP:"),
    ("pv24_partition_ratio_lo",
     r"F-based 95% interval of \[(\d\.\d+),", 5e-4,
     "F interval lower bound, partition sd ratio", "WRAP:"),
    ("pv24_partition_ratio_hi",
     r"F-based 95% interval of \[\d\.\d+,\s+(\d\.\d+)\]", 5e-4,
     "F interval upper bound, partition sd ratio", "WRAP:"),
    ("pv24_sqrt_sites",
     r"and the site-count\s+prediction \((\d\.\d+)\)", 5e-4,
     "sqrt(sites) prediction", "WRAP:"),
    # A2: the rotation-null family p. Guarded to lines that are about the
    # rotation null, so the pattern can be loose enough to catch "p = 0.5".
    # Session 59: the floor was 2e-6 until 2026-09-23, believed necessary so
    # that the same pattern also accepts "0.001". It was not: `_tol_for` gives
    # "0.001" its own printed half-width (5e-4), whatever the floor. What 2e-6
    # did was let "0.000999" drift to 0.000997-0.001001 unseen (measured: a
    # one-unit move, 0.000999 -> 0.000998, passed `--strict`). The floor is now
    # the printed half-width of "0.000999", and "0.001" still passes.
    ("rotation_p", r"\bp\s*(?:=|&nbsp;=&nbsp;)\s*\**\s*(\d?\.\d+)", 5e-7,
     "rotation-null family p-value",
     r"rotation|B\s*=\s*1,?000|1,000 draws|B = 1000"),
    ("rotation_observed",
     r"[Oo]bserved family mean r[^\d\n]{0,12}\**(\d?\.\d+)", 5e-5,
     "observed family mean r", None),
    # NOT CHECKED: the rotation-null draw count B. `pipeline/README.md` (Multiplicity)
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
    # Guarded on `k=10` specifically, not on "panel size": the manuscript's Reliability
    # subsection
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
    # E10, 2026-09-17: the NSCLC outcome arm's power curve. One sentence states
    # all three and wraps, so each is anchored on its own words.
    ("nsclc_power_003_pct", r"the NSCLC arm's power is (\d+)%", 0.5,
     "NSCLC outcome-arm power at a C-index advantage of 0.03", None),
    ("nsclc_power_002_pct", r"power is \d+% \((\d+)% at 0\.02", 0.5,
     "NSCLC outcome-arm power at 0.02", "WRAP:the NSCLC arm's power is"),
    ("nsclc_power_004_pct", r"(\d+)% at 0\.04\)", 0.5,
     "NSCLC outcome-arm power at 0.04", "WRAP:the NSCLC arm's power is"),
    # E13 and D4, 2026-09-16: read from the frozen pan-cancer summary's notes.
    ("pancancer_type_only_c", r"cancer type alone reaches C = (\d\.\d{3})", 5e-4,
     "concordance of cancer type alone on pooled pan-TCGA PFI", None),
    ("pancancer_isi_dropped", r"Of these, (\d+) patients \(16 from", 0.5,
     "pan-cancer patients without a cancer type, dropped from the index", None),
    ("pancancer_isi_complete_n", r"use the ([\d,]+) complete cases", 0.5,
     "pan-cancer patients the index and its bootstrap resample", "WRAP:Of these, "),
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
    # Session 51: the dry run's other counter figures. Each is anchored on its
    # own wording, not on digits, and each has a live self-test row. The naive
    # length read 2,920 for days after the counter stopped printing it.
    ("abstract_headroom",
     r"2,?600\**\s*\(\d{1,3}%\),\s+headroom\s+(\d{1,4})\b", 0.5,
     "abstract headroom, dry-run form (limit minus live count)", "WRAP:"),
    ("abstract_title_chars",
     r"\*\*Title\*\*\s+\((\d{2,3}) characters under the AACR rule\)", 0.5,
     "abstract title characters, counted by 08_count_abstract.py just now", "WRAP:"),
    ("abstract_naive_len",
     r"`len\(\)` including spaces reads ([\d,]{4,6})\b", 0.5,
     "abstract title + body INCLUDING spaces, from 08_count_abstract.py", "WRAP:"),
    ("abstract_body_words",
     r"headroom\s+\d{1,4},\s+body\s+(\d{2,4})\s+words", 0.5,
     "abstract body word count, from 08_count_abstract.py", "WRAP:"),
    # Session 52: the cover letter's length and its one-paragraph fallback.
    # Anchored on the wording, not the digits. The stated 614 was accurate on
    # 2026-09-02 and 180 words stale by 2026-09-18 because every correction
    # note landed inside the measured span; these three rows are why that
    # cannot happen quietly again.
    ("cover_letter_body_words",
     r"\*\*Body length:\s+(\d{3,4})\s+words\*\*", 0.5,
     "cover letter body words AS PASTED, recounted just now", None),
    ("cover_letter_para4_words",
     r"\*\*(\d{2,4})\s+words measured\*\*", 0.5,
     "cover letter paragraph 4 word count, recounted just now", None),
    ("cover_letter_body_minus_para4",
     r"words measured\*\*\s+—\s+and\s+(\d{3,4})\s+remain", 0.5,
     "cover letter body minus paragraph 4, recounted just now", "WRAP:"),
    # Session 55: the paragraph enumeration the author reads to decide WHICH
    # part of the argument to cut. Every member was exact when first measured;
    # what was missing is that the list sums to 608 while the header says 647
    # remain, and nothing said the other 39 words are the letterhead
    # scaffolding. All seven are WRAP-scoped: the sentence wraps between the
    # boilerplate figure and its noun, which is how the project's last
    # line-scoped claim (a Mann-Whitney P in this same file) hid a defect.
    ("cover_letter_para_submit_words",
     r"remaining paragraphs are (\d{2,3}) / \d{2,3} / \d{2,3} / \d{2,3}"
     r" words of argument", 0.5,
     "cover letter submission paragraph words, recounted just now",
     "WRAP:remaining paragraphs are"),
    ("cover_letter_para_finding_words",
     r"remaining paragraphs are \d{2,3} / (\d{2,3}) / \d{2,3} / \d{2,3}"
     r" words of argument", 0.5,
     "cover letter finding paragraph words, recounted just now",
     "WRAP:remaining paragraphs are"),
    ("cover_letter_para_immune_words",
     r"remaining paragraphs are \d{2,3} / \d{2,3} / (\d{2,3}) / \d{2,3}"
     r" words of argument", 0.5,
     "cover letter immune-signal paragraph words, recounted just now",
     "WRAP:remaining paragraphs are"),
    ("cover_letter_para_notwhat_words",
     r"remaining paragraphs are \d{2,3} / \d{2,3} / \d{2,3} / (\d{2,3})"
     r" words of argument", 0.5,
     "cover letter scope-limits paragraph words, recounted just now",
     "WRAP:remaining paragraphs are"),
    ("cover_letter_para_scope_words",
     r"words of argument and (\d{2,3}) of required boilerplate", 0.5,
     "cover letter boilerplate paragraph words, recounted just now",
     "WRAP:remaining paragraphs are"),
    ("cover_letter_enumerated_words",
     r"required boilerplate — (\d{3}) in all", 0.5,
     "cover letter enumerated paragraphs, summed just now",
     "WRAP:remaining paragraphs are"),
    ("cover_letter_scaffold_words",
     r"in all, the other (\d{2}) being the date", 0.5,
     "cover letter letterhead scaffolding words, derived just now",
     "WRAP:remaining paragraphs are"),
    # Session 55: the cover letter's SCIENCE, found by mutation sweep rather
    # than by reading the coverage list. Before this, the only numbers the
    # letter had bound were its own word counts, one Mann-Whitney P and the
    # two ISI triples; a sweep of all 47 unchecked literals -- each moved one
    # step in its last printed place against the real file, `--strict` without
    # `--fast`, restored by sha256 -- came back 38 UNCAUGHT, and nine of those
    # are results an editor is invited to check. This is the outward-facing
    # document, so it is the first place to sweep, not the last.
    ("nsclc_n",
     r"in a separate ([\d,]{3,5})-patient NSCLC", 0.5,
     "cover letter NSCLC cohort size", "WRAP:immune-specificity index"),
    ("rotation_p",
     r"between-signature correlation gives \*P\* = (\d+\.\d+)", 0.0005,
     "cover letter family-level rotation p", "WRAP:immune-specificity index"),
    ("rotation_B",
     r"\*P\* = [\d.]+ at ([\d,]+) draws", 0.5,
     "cover letter rotation null draws", "WRAP:immune-specificity index"),
    # The reliability sentence: "their Cronbach's alpha falls from 0.97 to
    # 0.80 while curated signatures hold above 0.97". The first two ARE
    # measured quantities and are bound here. The third 0.97 is NOT, and the
    # attempt to bind it is worth recording: it was written as a point
    # estimate against `pancancer_alpha_obs_resid`, and `--strict` refused it
    # at a tolerance of half the last printed place, because that median is
    # 0.9764 and rounds to 0.98. The number is not a rounded median -- it is
    # an INEQUALITY BOUND the author chose, and 0.9764 > 0.97 is true. A
    # threshold is not a measurement, and a rule that checks it as one reports
    # a defect that does not exist. It stays unbound, deliberately. B-30 was
    # answered on 2026-09-23 and the clause now says "the median curated set"
    # (or "signature") in the four live summaries -- which silently took them
    # OUT of these two patterns: `--strict` still passed, because the folded
    # superseded abstract kept the old wording and the patterns matched there,
    # and only the cover letter's self-test rows noticed. Both wordings are
    # admitted now.
    ("pancancer_alpha_null_raw",
     r"from (\d+\.\d+) to \d+\.\d+ while (?:the median )?curated", 0.005,
     "random-set alpha before residualisation, all four summaries", None),
    ("pancancer_alpha_null_resid",
     r"from \d+\.\d+ to (\d+\.\d+) while (?:the median )?curated", 0.005,
     "random-set alpha after residualisation, all four summaries", None),
    ("nsclc_rel_gap_fold",
     r"gap is then (\d+\.\d+)-fold \(NSCLC\)", 0.05,
     "cover letter NSCLC reliability-gap fold", None),
    ("pancancer_rel_gap_fold",
     r"-fold \(NSCLC\) to (\d+)-fold \(pan-cancer\)", 0.5,
     "cover letter pan-cancer reliability-gap fold", None),
    ("small_panel_k160",
     r"scales inversely with panel size — (\d+\.\d+) at 160 genes", 0.0005,
     "cover letter reliability gap at 160 genes", "WRAP:Random gene sets"),
    ("small_panel_k10",
     r"at 160 genes rising to (\d+\.\d+) at 10", 0.0005,
     "cover letter reliability gap at 10 genes", "WRAP:Random gene sets"),
    # The reproducibility promise. 0.3182 is bound; its partner 0.2964 is
    # deliberately NOT, and that is a recorded decision rather than an
    # oversight: it is the pre-fix index the HPC4 run returned, it lives in
    # `results/scorer_sensitivity_hpc4/run.log`, and this file already states
    # that directory's numbers are diagnostic and not quotable. Binding it
    # would mean deriving an authority from a .log, which F6.8 decided against.
    ("nsclc_isi",
     r"the (\d+\.\d+) frozen on macOS / arm64", 0.00005,
     "cover letter frozen NSCLC index, the reproducibility promise", None),
    # Session 53: the manuscript abstract's two counts and the SUM of them, the
    # same triple the cover letter carries and for the same reason. The sum has
    # its own row deliberately: a conclusion rests on it ("inside the observed
    # ceiling"), and September's abstract-headroom defect was exactly the shape
    # of both parts right and the arithmetic between them wrong. B-30 edits a
    # line inside the measured span, so all three move when the author chooses.
    # Session 53: the two figures AACR requires the COVER LETTER to state, both
    # of them counts of ANOTHER document, which is the shape that goes stale
    # unwatched. The word count is WRAP-scoped because the sentence wraps
    # between the number and its qualifying clause.
    ("paper_body_words",
     r"\*\*([\d,]{4,6})\s+words\*\* excluding the cover page", 0.5,
     "cover letter: manuscript body word count", "WRAP:"),
    ("paper_figure_count",
     r"with \*\*(\d{1,2})\s+figures and no tables\*\*", 0.5,
     "cover letter: manuscript figure count", "WRAP:"),
    ("paper_abstract_body_words",
     r"measures \*\*(\d{3})\s+words of body", 0.5,
     "manuscript abstract body words, recounted just now", None),
    ("paper_abstract_significance_words",
     r"words of body and (\d{2})\s+of Significance", 0.5,
     "manuscript abstract Significance words, recounted just now", None),
    ("paper_abstract_total_words",
     r"of Significance,\s+(\d{3})\s+in total\*\*", 0.5,
     "manuscript abstract body + Significance, recounted just now", None),
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
    # Session 60 (H86): the pinned-split NSCLC re-run's secondaries, which
    # Results reports beside the frozen ones in a paragraph of its own. Worded
    # "of" rather than "=" so the frozen pan-TCGA pattern above cannot claim
    # them; WRAP-scoped on the paragraph's opening words.
    ("pinned_nsclc_median_dr", r"NSCLC gives a median Δr of \+(\d?\.\d{3})", 5e-4,
     "pinned split: NSCLC median Delta r", "WRAP:Pinned random-patient split"),
    ("pinned_nsclc_dmae", r"NSCLC gives a median Δr of \+[\d.]+, Δ-MAE (\d?\.\d{4}) \[", 5e-5,
     "pinned split: NSCLC Delta-MAE", "WRAP:Pinned random-patient split"),
    ("pinned_nsclc_emb_over_cov",
     r"embedding advantage over covariates of \+(\d?\.\d{3})", 5e-4,
     "pinned split: NSCLC embedding over covariates", "WRAP:Pinned random-patient split"),
    ("pinned_nsclc_site_auroc",
     r"covariates of \+[\d.]+ and a median site AUROC of (\d?\.\d{3})", 5e-4,
     "pinned split: NSCLC median site AUROC", "WRAP:Pinned random-patient split"),
    # ...and the pan-cancer half, run on macOS (the VPN was down) and reported
    # in the same paragraph; each anchored on the cohort's own words.
    ("pinned_pancancer_median_dr", r"pan-cancer gives a median Δr of \+(\d?\.\d{3})", 5e-4,
     "pinned split: pan-cancer median Delta r", "WRAP:Pinned random-patient split"),
    ("pinned_pancancer_dmae",
     r"pan-cancer gives a median Δr of \+[\d.]+, Δ-MAE (\d?\.\d{4}) \[", 5e-5,
     "pinned split: pan-cancer Delta-MAE", "WRAP:Pinned random-patient split"),
    ("pinned_pancancer_emb_over_cov_abs",
     r"an embedding advantage of [−-](\d?\.\d{3}) and", 5e-4,
     "pinned split: pan-cancer embedding over covariates (magnitude)",
     "WRAP:Pinned random-patient split"),
    ("pinned_pancancer_site_auroc",
     r"advantage of [−-][\d.]+ and a\s+median site AUROC of (\d?\.\d{3})", 5e-4,
     "pinned split: pan-cancer median site AUROC", "WRAP:Pinned random-patient split"),
    # The two numbers A1's verdict rests on. Both patterns are anchored on the
    # surrounding prose rather than on the digits, so the injection rows in
    # `--self-test` can move the value and still be matched -- the vacuity
    # defect this file has shipped four times.
    ("pancancer_min_excess_lo",
     r"`excess_lo` is \+(\d?\.\d+), nowhere near the boundary", 5e-3,
     "smallest pan-cancer excess_lo, the ISI's multiplicity margin", None),
    ("pancancer_hypoxia_excess_lo",
     r"Hypoxia beats its null with `excess_lo = \+(\d?\.\d+)`", 5e-5,
     "hypoxia outcome-arm excess_lo, axis-residualised", None),
    # Session 53: the same six numbers where the POSTER prints them. All five
    # excesses sit on one line of poster.html, so a line-scoped pattern reaches
    # them; each is anchored on the signature name that PRECEDES it, never on
    # the digits. Tolerance is 5e-4 because the poster rounds to three decimals
    # (0.0686961... printed as 0.069) -- tight enough that a transcription slip
    # of one in the last printed place still fails.
    ("pancancer_outcome_excess_g2m",
     r"G2M (\d?\.\d{2,}), E2F", 5e-4,
     "poster panel: G2M outcome excess, axis-residualised", None),
    ("pancancer_outcome_excess_e2f",
     r"E2F (\d?\.\d{2,}), angiogenesis", 5e-4,
     "poster panel: E2F outcome excess, axis-residualised", None),
    ("pancancer_outcome_excess_angio",
     r"angiogenesis (\d?\.\d{2,}), EMT", 5e-4,
     "poster panel: angiogenesis outcome excess, axis-residualised", None),
    ("pancancer_outcome_excess_emt",
     r"EMT (\d?\.\d{2,}), hypoxia", 5e-4,
     "poster panel: EMT outcome excess, axis-residualised", None),
    ("pancancer_outcome_excess_hypoxia",
     r"hypoxia (\d?\.\d{2,}) \(Fisher exact", 5e-4,
     "poster panel: hypoxia outcome excess, axis-residualised", None),
    ("pancancer_hypoxia_excess_lo",
     r"narrowest at excess_lo = \+(\d?\.\d+)", 5e-5,
     "poster panel: hypoxia excess_lo, the multiplicity margin", None),
    # Session 53: pipeline/README.md's headline results table. It restates the
    # frozen numbers for BOTH cohorts and it is DEPOSITED -- it ships in the
    # public snapshot -- so if it drifts, the published artefact contradicts the
    # paper in the one file a reader of the repository opens first. 21 of its
    # 459 literals were checked and none of these were among them. Each row is
    # one line, so line scope reaches them; every pattern is anchored on the row
    # label or the cohort name, never on the digits.
    ("pancancer_n", r"Pan-TCGA \(n=([\d,]+), 31 types\)", 0.5,
     "deposited README results table: pan-TCGA cohort size", None),
    ("nsclc_n", r"\| NSCLC \(n=([\d,]+)\) \|", 0.5,
     "deposited README results table: NSCLC cohort size", None),
    ("pancancer_isi", r"\| ISI \| \*\*(\d?\.\d+) \[", 5e-5,
     "deposited README results table: pan-TCGA ISI", None),
    ("pancancer_ci_lo", r"\| ISI \| \*\*\d?\.\d+ \[(\d?\.\d+),", 5e-5,
     "deposited README results table: pan-TCGA CI low", None),
    ("pancancer_ci_hi", r"\| ISI \| \*\*\d?\.\d+ \[\d?\.\d+, (\d?\.\d+)\]", 5e-5,
     "deposited README results table: pan-TCGA CI high", None),
    ("nsclc_isi", r"\]\*\* \| (\d?\.\d+) \[", 5e-5,
     "deposited README results table: NSCLC ISI", None),
    ("nsclc_ci_lo", r"\]\*\* \| \d?\.\d+ \[(\d?\.\d+),", 5e-5,
     "deposited README results table: NSCLC CI low", None),
    ("nsclc_ci_hi", r"\]\*\* \| \d?\.\d+ \[\d?\.\d+, (\d?\.\d+)\]", 5e-5,
     "deposited README results table: NSCLC CI high", None),
    ("nsclc_mde", r"0/32 \(MDE (\d?\.\d+)\)", 5e-4,
     "deposited README results table: NSCLC minimum detectable effect", None),
    ("pancancer_site_auroc", r"\| Site AUROC \| (\d?\.\d+) \|", 5e-4,
     "deposited README results table: pan-TCGA site AUROC", None),
    ("nsclc_site_auroc", r"\| Site AUROC \| \d?\.\d+ \| (\d?\.\d+) \|", 5e-4,
     "deposited README results table: NSCLC site AUROC", None),
    # ...and the manuscript's copy of the same five, which the poster must agree
    # with. Verified 2026-09-20 that all five already agree with each other and
    # with `outcome_arm.csv`; these rows are what keeps that true. Every one is
    # WRAP-scoped and uses `\s+` between tokens, because this Results paragraph
    # ALREADY wraps between a signature's name and its value -- "**E2F targets**"
    # ends one line and "(0.062)" opens the next -- which is precisely the shape
    # a line-scoped pattern cannot see (the defect of session 50's F3.8).
    ("pancancer_outcome_excess_g2m",
     r"G2M checkpoint\*\*\s+\(excess\s+(\d?\.\d{2,})\)", 5e-4,
     "manuscript Results: G2M outcome excess", "WRAP:"),
    ("pancancer_outcome_excess_e2f",
     r"E2F targets\*\*\s+\((\d?\.\d{2,})\)", 5e-4,
     "manuscript Results: E2F outcome excess (wraps)", "WRAP:"),
    ("pancancer_outcome_excess_angio",
     r"\*\*angiogenesis\*\*\s+\((\d?\.\d{2,})\)", 5e-4,
     "manuscript Results: angiogenesis outcome excess", "WRAP:"),
    ("pancancer_outcome_excess_emt",
     r"mesenchymal transition\*\*\s+\((\d?\.\d{2,})\)", 5e-4,
     "manuscript Results: EMT outcome excess", "WRAP:"),
    ("pancancer_outcome_excess_hypoxia",
     r"\*\*hypoxia\*\*\s+\((\d?\.\d{2,})\)", 5e-4,
     "manuscript Results: hypoxia outcome excess", "WRAP:"),
    # Session 54: THE REST OF THE POSTER. Session 53 bound the outcome panel
    # and called it the least-checked thing in the project; it was not the
    # only one. A mutation sweep (each literal moved by one digit in its last
    # printed place, against the real file, `--strict` run each time) found
    # the poster's LEAD SENTENCE -- both cohorts' ISI and intervals, printed
    # largest on the board -- guarded by nothing, along with every number in
    # the reliability panel, the small-panel table, the per-signature ranges,
    # the rotation null, both widened intervals and most of panel 5. They were
    # not wrong. They were unchecked, which on an A0 board that cannot be
    # corrected once it is up is the same risk. The poster writes HTML
    # entities (`&ndash;`, `&minus;`, `&rsquo;`) where the manuscript writes
    # the characters, which is why the manuscript's patterns never reached it;
    # these are anchored on the poster's own prose, never on the digits.
    ("pancancer_isi", r"ISI = (\d?\.\d+) \(95% CI [\d.]+&ndash;[\d.]+\) pan-TCGA", 5e-4,
     "poster lead: pan-TCGA ISI", None),
    ("pancancer_ci_lo", r"ISI = [\d.]+ \(95% CI (\d*\.\d+)&ndash;", 5e-4,
     "poster lead: pan-TCGA CI low", None),
    ("pancancer_ci_hi", r"\(95% CI [\d.]+&ndash;(\d*\.\d+)\) pan-TCGA", 5e-4,
     "poster lead: pan-TCGA CI high", None),
    ("nsclc_isi", r"pan-TCGA and (\d?\.\d+) \([\d.]+&ndash;[\d.]+\) in NSCLC", 5e-4,
     "poster lead: NSCLC ISI", None),
    ("nsclc_ci_lo", r"pan-TCGA and [\d.]+ \((\d*\.\d+)&ndash;", 5e-4,
     "poster lead: NSCLC CI low", None),
    ("nsclc_ci_hi", r"&ndash;(\d*\.\d+)\) in NSCLC\.", 5e-4,
     "poster lead: NSCLC CI high", None),
    ("pancancer_excess_min", r"Per-signature excess (\d*\.\d+)&ndash;[\d.]+ pan-cancer", 5e-4,
     "poster panel 4: smallest pan-TCGA per-signature excess", None),
    ("pancancer_excess_max", r"Per-signature excess [\d.]+&ndash;(\d*\.\d+) pan-cancer", 5e-4,
     "poster panel 4: largest pan-TCGA per-signature excess", None),
    ("nsclc_excess_min", r"pan-cancer, (\d*\.\d+)&ndash;[\d.]+ NSCLC\. Family-level", 5e-4,
     "poster panel 4: smallest NSCLC per-signature excess", None),
    ("nsclc_excess_max", r"pan-cancer, [\d.]+&ndash;(\d*\.\d+) NSCLC\. Family-level", 5e-4,
     "poster panel 4: largest NSCLC per-signature excess", None),
    ("rotation_null_mean", r"against a null mean of (\d*\.\d+) \(95% range", 5e-4,
     "poster panel 4: rotation-null family mean r", None),
    ("rotation_null_lo", r"\(95% range (\d*\.\d+)&ndash;[\d.]+\)", 5e-4,
     "poster panel 4: rotation-null 95% range, low", None),
    ("rotation_null_hi", r"\(95% range [\d.]+&ndash;(\d*\.\d+)\)", 5e-4,
     "poster panel 4: rotation-null 95% range, high", None),
    ("pancancer_r_axis", r"Median image&ndash;axis r = (\d*\.\d+) and", 5e-4,
     "poster panel 4: pan-TCGA median image-axis r", None),
    ("nsclc_r_axis", r"Median image&ndash;axis r = [\d.]+ and (&minus;\d*\.\d+),", 5e-4,
     "poster panel 4: NSCLC median image-axis r (negative)", None),
    ("nsclc_honest_ci_lo", r"widens the interval to \[(\d*\.\d+), [\d.]+\]", 5e-5,
     "poster panel 4: NSCLC partition-widened interval, low", None),
    ("nsclc_honest_ci_hi", r"widens the interval to \[[\d.]+, (\d*\.\d+)\]", 5e-5,
     "poster panel 4: NSCLC partition-widened interval, high", None),
    ("nsclc_interval_understated_pct",
     r"widens the interval to \[[\d.]+, [\d.]+\], (\d*\.\d+)% wider than reported", 0.05,
     "poster panel 4: NSCLC interval widening, percent", None),
    ("pancancer_honest_ci_lo", r"widens it to \[(\d*\.\d+), [\d.]+\]", 5e-5,
     "poster panel 4: pan-TCGA partition-widened interval, low", None),
    ("pancancer_honest_ci_hi", r"widens it to \[[\d.]+, (\d*\.\d+)\]", 5e-5,
     "poster panel 4: pan-TCGA partition-widened interval, high", None),
    ("pancancer_interval_understated_pct",
     r"widens it to \[[\d.]+, [\d.]+\], <b>(\d*\.\d+)% wider</b>", 0.05,
     "poster panel 4: pan-TCGA interval widening, percent",
     "WRAP:Pan-cancer partition sd"),
    ("pancancer_partition_var_pct", r"larger share pan-cancer \((\d*\.\d+)% vs [\d.]+%\)", 0.05,
     "poster panel 4: partition share of variance, pan-TCGA", None),
    ("nsclc_partition_var_pct", r"larger share pan-cancer \([\d.]+% vs (\d*\.\d+)%\)", 0.05,
     "poster panel 4: partition share of variance, NSCLC", None),
    # Panel 2, reliability. The paragraph wraps twice inside its own claims.
    ("pancancer_alpha_null_raw", r"median Cronbach&rsquo;s &alpha; (\d*\.\d+) pan-cancer", 5e-4,
     "poster panel 2: median random-set alpha, raw, pan-TCGA", None),
    ("pancancer_alpha_null_resid", r"drops the null&rsquo;s &alpha; to <b>(\d*\.\d+)</b>", 5e-4,
     "poster panel 2: median random-set alpha, residualized, pan-TCGA",
     "WRAP:Random gene sets are internally consistent"),
    ("pancancer_alpha_obs_resid", r"curated signatures hold at <b>(\d*\.\d+)</b>", 5e-4,
     "poster panel 2: median curated alpha, residualized, pan-TCGA", None),
    ("pancancer_rel_gap_raw", r"gap is \+(\d*\.\d+) pan-cancer and \+[\d.]+ NSCLC on raw", 5e-4,
     "poster panel 2: curated-random alpha gap, raw, pan-TCGA",
     "WRAP:Random gene sets are internally consistent"),
    ("nsclc_rel_gap_raw", r"\+[\d.]+ pan-cancer and \+(\d*\.\d+) NSCLC on raw scores", 5e-4,
     "poster panel 2: curated-random alpha gap, raw, NSCLC", None),
    ("pancancer_rel_gap_resid", r"against <b>\+(\d*\.\d+)</b> and <b>\+[\d.]+</b> on the", 5e-4,
     "poster panel 2: curated-random alpha gap, residualized, pan-TCGA", None),
    ("nsclc_rel_gap_resid", r"against <b>\+[\d.]+</b> and <b>\+(\d*\.\d+)</b> on the", 5e-4,
     "poster panel 2: curated-random alpha gap, residualized, NSCLC", None),
    ("small_panel_k10", r'<td class="hi">10</td><td class="hi">(\d*\.\d+)</td>', 5e-5,
     "poster panel 2 table: alpha gap at k = 10", None),
    ("small_panel_k20", r"<td>20</td><td>(\d*\.\d+)</td>", 5e-5,
     "poster panel 2 table: alpha gap at k = 20", None),
    ("small_panel_k40", r"<td>40</td><td>(\d*\.\d+)</td>", 5e-5,
     "poster panel 2 table: alpha gap at k = 40", None),
    ("small_panel_k80", r"<td>80</td><td>(\d*\.\d+)</td>", 5e-5,
     "poster panel 2 table: alpha gap at k = 80", None),
    ("small_panel_k160", r"<td>160</td><td>(\d*\.\d+)</td>", 5e-5,
     "poster panel 2 table: alpha gap at k = 160", None),
    ("small_panel_spearman", r"Spearman &rho; = (&minus;\d*\.\d+), p &lt;", 5e-4,
     "poster panel 2: Spearman(k, gap)", None),
    ("small_panel_ratio_k10_k160", r"<b>(\d*\.\d+)&times;</b> larger at", 0.05,
     "poster panel 2: alpha gap at k = 10 over k = 160", None),
    # Panel 3, methods, and panel 5, limitations.
    ("nsclc_n", r"NSCLC ([\d,]+) \(\d+ events\)", 0.5,
     "poster panel 3: NSCLC cohort size", None),
    ("nsclc_n_events", r"NSCLC [\d,]+ \((\d+) events\)", 0.5,
     "poster panel 3: NSCLC PFI events", None),
    ("sortaudit_macos_genes", r"Xena TOIL expression, ([\d,]+) genes", 0.5,
     "poster panel 3: genes in the TOIL expression matrix (sort audit's count)", None),
    ("nsclc_site_auroc", r"<b>AUROC [\d.]+</b> pan-cancer, (\d*\.\d+) NSCLC\.", 5e-4,
     "poster panel 5: NSCLC site AUROC, median",
     "WRAP:Site is almost perfectly recoverable"),
    ("pancancer_emb_over_cov", r"&Delta;r = (&minus;\d*\.\d+) over site", 5e-4,
     "poster panel 5: pan-TCGA embedding over covariates, Delta r", None),
    ("nsclc_emb_over_cov", r"stage; NSCLC \+(\d*\.\d+)\) &mdash; hence", 5e-4,
     "poster panel 5: NSCLC embedding over covariates, Delta r", None),
    ("label_site_variance_r2", r"site R&sup2; (\d*\.\d+) crude", 5e-4,
     "poster panel 5: label-side site R^2, crude, median", None),
    ("label_site_given_type_r2", r"crude, (\d*\.\d+) given cancer type", 5e-4,
     "poster panel 5: label-side site R^2 given cancer type, median", None),
    ("nsclc_isi", r"\(NSCLC ISI (\d*\.\d+) vs [\d.]+\)", 5e-5,
     "poster panel 5: the frozen NSCLC ISI beside HPC4's", None),
    # The outcome panel's two p-values and the pan-cancer MDE.
    ("immune_mwu_p_2s",
     r'<div class="n">p = (\d*\.\d+)</div>\s*<div class="l">Mann&ndash;Whitney across process', 5e-4,
     "poster chip: Mann-Whitney p across process categories, two-sided", "WRAP:"),
    ("immune_fisher_p_2s", r"\(Fisher exact p = (\d*\.\d+)\)", 5e-4,
     "poster panel: Fisher exact p, immune vs other winners, two-sided", None),
    ("pancancer_mde", r"pan-cancer MDE (\d*\.\d+)\)", 5e-4,
     "poster panel: pan-TCGA minimum detectable effect",
     "WRAP:Pan-TCGA, progression-free interval"),
    # The deposited landing page, which restates two of the same results.
    ("nsclc_site_auroc", r"\| Site AUROC \(median\) \| [\d.]+ \| (\d*\.\d+) \|", 5e-4,
     "snapshot README results table: NSCLC site AUROC", None),
    ("a7_meanz_rel_null", r"random-set reliability from (\d*\.\d+) to [\d.]+", 5e-4,
     "snapshot README: NSCLC random-set reliability, mean-z", None),
    ("a7_ssgsea_rel_null", r"random-set reliability from [\d.]+ to (\d*\.\d+)", 5e-4,
     "snapshot README: NSCLC random-set reliability, ssGSEA", None),
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
    # E14: the composition of the NSCLC/pan-cancer gap, each on its own anchor.
    ("lsv_pan_adj_pct", r"to an adjusted (\d+\.\d)% pan-cancer", 5e-2,
     "adjusted site R^2 of the label, pan-TCGA median (percent)", "WRAP:"),
    ("lsv_nsclc_adj_pct", r"adjusted \d+\.\d% pan-cancer and (\d+\.\d)%\s+in NSCLC", 5e-2,
     "adjusted site R^2 of the label, NSCLC median (percent)", "WRAP:"),
    ("lsv_pan_type_pct", r"Cancer type alone accounts for a median (\d+\.\d)% of label", 5e-2,
     "type-only R^2 of the label, pan-TCGA median (percent)", "WRAP:"),
    ("lsv_nsclc_type_pct", r"of label variance\s+pan-cancer and (\d+\.\d)% across the two NSCLC", 5e-2,
     "type-only R^2 of the label, NSCLC median (percent)", "WRAP:"),
    ("lsv_nsclc_site_given_type_pct", r"Conditioned on type, site explains a median\s+(\d+\.\d)% in NSCLC",
     5e-2, "site-given-type R^2 of the label, NSCLC median (percent)", "WRAP:"),
    ("label_site_given_type_pct", r"median\s+\d+\.\d% in NSCLC and (\d+\.\d)% pan-cancer", 5e-2,
     "site-given-type R^2 of the label, pan-TCGA (percent), beside NSCLC", "WRAP:"),
    # The single-probe calibration (script 11's A4b), now named as such.
    ("control_c_calibrated_pct",
     r"allograft rejection alone,\s+gave\s+(\d{1,3}(?:\.\d+)?)\s*%", 5e-2,
     "Control C, permutation-calibrated, one probe (percent)", "WRAP:"),
    # E14: the same calibration over all 16 signatures, both cohorts.
    ("e14_pancancer_median_pct", r"leaves a median of (\d\.\d)% across the 16", 5e-2,
     "Control C calibrated, pan-cancer median (percent)", "WRAP:"),
    ("e14_pancancer_min_pct", r"across the 16\s+signatures \(range (\d\.\d)% to", 5e-2,
     "Control C calibrated, pan-cancer minimum (percent)", "WRAP:"),
    ("e14_pancancer_max_pct", r"\(range \d\.\d% to (\d\.\d)%; \d+ of 16", 5e-2,
     "Control C calibrated, pan-cancer maximum (percent)", "WRAP:"),
    ("e14_pancancer_n_sig", r"to \d\.\d%; (\d+) of 16 at p", 0.5,
     "Control C calibrated, pan-cancer signatures at p <= 0.05", "WRAP:"),
    ("e14_nsclc_median_pct", r"with a median of (\d\.\d)% across the 16 signatures in NSCLC", 5e-2,
     "Control C calibrated, NSCLC median (percent)", "WRAP:"),
    ("e14_nsclc_n_sig", r"signatures in NSCLC \((\d+) of 16 at", 0.5,
     "Control C calibrated, NSCLC signatures at p <= 0.05", "WRAP:"),
    ("e14_pancancer_median_r2", r"a median (\d\.\d{3}) after\s+within-type permutation", 5e-4,
     "Control C calibrated, pan-cancer median (poster)", "WRAP:"),
    # WITHDRAWN 2026-09-17 (E17): "plate within site explains 0.000". The plate
    # field is empty in every input row, so that zero was true by construction
    # and said nothing. The authority `label_plate_within_site` is still read so
    # a future non-zero value would be seen, but no document may quote it as
    # evidence.
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
    # Session 50 (F3.8): the cover letter wraps between "Mann–Whitney" and
    # "*P*", so the line-scoped row above (guard on the SAME line) never saw its
    # value -- which was 0.016, the IL6-counted-immune p, beside the 0-of-6 vs
    # 5-of-10 counts that give 0.031. Anchored on the counts, across the wrap.
    ("immune_mwu_p_2s",
     r"MSigDB-immune versus \d+ of \d+ other signatures, Mann.Whitney \*P\*\s*=\s*(\d?\.\d{2,})",
     5e-4, "immune-vs-other Mann-Whitney p beside its counts, across a line wrap",
     "WRAP:"),
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
    # which the ABSTRACT writes but the PAPER does not -- the manuscript's abstract
    # and its folded superseded abstract both write "Cronbach's α from 0.97 to 0.80". So the manuscript's
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
    # existed and matches the manuscript's outcome-arm Results sentence -- but
    # 07-ABSTRACT-DRAFT.md wraps between "minimum" and "detectable" (in "What
    # changed from draft 1"), so the abstract's
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
    # reports the "7,168" in 07-ABSTRACT-DRAFT.md's BODY Methods sentence as
    # unchecked but NOT its "31".
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
    # documents (00-PROJECT-BRIEF.md section 9, 08-READINESS-ASSESSMENT.md's "What
    # is genuinely strong") quote
    # a DIFFERENT Spearman -- the set-size scaling, -0.50/-0.59 -- and a loose
    # anchor would report those as drift the day either file joins DOCS.
    # The manuscript's folded superseded abstract wraps directly after
    # "(Spearman", so only the
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
    # Session 50 (F3.8): five correct values the line-scoped rows could not see
    # because each wraps across a line break, and one whose guard sits on
    # another line. Each is anchored on its own sentence.
    ("pancancer_n", r"Across ([\d,]{3,7}) TCGA patients in 31 cancer types", 0.5,
     "pan-TCGA cohort size, cover letter, across a wrap", "WRAP:"),
    ("nsclc_n", r"TCGA NSCLC \(n=([\d,]{3,7})\) as a second cohort", 0.5,
     "NSCLC cohort size, abstract draft, across a wrap", "WRAP:"),
    ("nsclc_n", r"with ([\d,]{3,7}) non-small cell lung cancer patients as a second", 0.5,
     "NSCLC cohort size, NextGen abstract, across a wrap", "WRAP:"),
    ("splitdet_macos_sites_tied", r"(\d{2,3}) of the 68 NSCLC sites sit in a tied size group", 0.5,
     "NSCLC tied sites, audit A9, across a wrap", "WRAP:"),
    ("rotation_p", r"rotation null p=(\d?\.\d+) · M_eff", 5e-7,
     "rotation-null family p-value, checklist summary, across a wrap", "WRAP:"),
    ("small_panel_spearman", r"Spearman\(k, gap\) = ([−-]?\d\.\d+)", 5e-4,
     "small-panel Spearman, audit A6 heading", None),
    # Session 50 (F7.4): the supplementary captions. CAPTIONS.md writes each
    # caption as one line, so these are line-scoped on anchors that occur only
    # there.
    ("nsclc_site_auroc", r"median AUROC \d\.\d{3} pan-TCGA, (\d\.\d{3}) NSCLC", 5e-4,
     "NSCLC site AUROC median, S7 caption", None),
    ("pancancer_n_sites_evaluable", r"rows are (\d+) pan-TCGA sites plus", 0.5,
     "pan-TCGA evaluable sites, S7 caption", None),
    ("nsclc_n_sites_evaluable", r"pan-TCGA sites plus (\d+) NSCLC sites twice", 0.5,
     "NSCLC evaluable sites, S7 caption", None),
    ("cap_combat_auroc_median", r"below chance \(median (\d\.\d{3})\)", 5e-4,
     "post-ComBat NSCLC site AUROC median, S7 caption", None),
    ("cap_pan_lsv_sites", r"enter \((\d+) of 619 pan-TCGA", 0.5,
     "pan-TCGA Control C sites with 10+ patients, S4/S5 captions", None),
    ("cap_nsclc_lsv_sites", r"pan-TCGA, (\d+) of 68 NSCLC\)", 0.5,
     "NSCLC Control C sites with 10+ patients, S4/S5 captions", None),
    ("nsclc_label_site_variance_r2", r"median R-squared of (\d\.\d{3}), against", 5e-4,
     "NSCLC Control C median R^2, S4 caption", None),
    ("label_site_variance_r2", r"median R-squared of \d\.\d{3}, against (\d\.\d{3}) here", 5e-4,
     "pan-TCGA Control C median R^2, S4 caption", None),
    ("nsclc_label_site_variance_r2", r"R-squared is (\d\.\d{3}) here against", 5e-4,
     "NSCLC Control C median R^2, S5 caption", None),
    ("label_site_variance_r2", r"R-squared is \d\.\d{3} here against (\d\.\d{3}) there", 5e-4,
     "pan-TCGA Control C median R^2, S5 caption", None),
    ("cap_nsclc_lsv_n", r"NSCLC contributes (\d+) patients across", 0.5,
     "NSCLC Control C patients, S5 caption", None),
    ("cap_nsclc_lsv_sites", r"patients across (\d+) usable sites against", 0.5,
     "NSCLC Control C sites, S5 caption", None),
    ("cap_pan_lsv_n", r"against pan-TCGA's ([\d,]+) across", 0.5,
     "pan-TCGA Control C patients, S5 caption", None),
    ("cap_pan_lsv_sites", r"pan-TCGA's [\d,]+ across (\d+),", 0.5,
     "pan-TCGA Control C sites, S5 caption", None),
    ("cap_emt_incr_unadj", r"transition: (\d\.\d{3}) unadjusted", 5e-4,
     "EMT image gain, ordinary R^2, S9 caption", None),
    ("cap_emt_incr_adj", r"unadjusted, (\d\.\d{3}) adjusted", 5e-4,
     "EMT image gain, adjusted R^2, S9 caption", None),
    ("small_panel_spearman", r"\(Spearman rho = (-?\d\.\d+)\)", 5e-3,
     "small-panel Spearman, S3 caption", None),
    # Session 50 (F1.1). The 216 is macOS's shipped-against-stable count; the
    # sentence used to read as a macOS-against-Linux count, which nothing
    # measured. Linux's own shipped-against-stable count sits beside it now.
    ("splitdet_hpc4_moved", r"change fold on macOS and (\d{2,3}) on Linux", 0.5,
     "patients assigned a different fold under the stable sort, Linux", "WRAP:"),
    # Limitation 8's yardstick ratios, read from `sort_fix_fold_change.json`.
    # None had an authority before session 50, and two were wrong at one
    # decimal. Every anchor is prose plus the NEIGHBOURING number, never the
    # digits being checked.
    ("sffc_gap_sds_small", r"0\.0218 is (\d\.\d) times the six-partition", 0.0,
     "A9 platform gap in six-partition sds", "WRAP:"),
    ("sffc_gap_sds_24", r"deviation of 0\.0089 \((\d\.\d) times the 24-partition value\)",
     0.0, "A9 platform gap in 24-partition sds", "WRAP:"),
    ("sffc_ns_move_sds_small", r"moves by 0\.0260, which is (\d\.\d) times", 0.0,
     "NSCLC sort-fix move in six-partition sds", "WRAP:"),
    ("sffc_ns_move_sds_24",
     r"\((\d\.\d) times the 24-partition value measured under the pinned ordering\)",
     0.0, "NSCLC sort-fix move in 24-partition sds", "WRAP:"),
    ("sffc_ns_n_pairs", r"of the (\d+) pairs among the 24 partitions, \d+ differ by as much",
     0.0, "pairs of NSCLC 24 partitions", "WRAP:"),
    ("sffc_ns_pairs_ge_move", r"pairs among the 24 partitions, (\d+) differ by as much",
     0.0, "NSCLC partition pairs differing by at least the sort-fix move", "WRAP:"),
    ("sffc_pc_move_sds_small", r"\*\*(\d\.\d+) times the pan-cancer between-partition", 0.0,
     "pan-cancer sort-fix move in five-partition sds", "WRAP:"),
    ("sffc_pc_move_sds_24", r"\((\d\.\d+) times the 24-partition value\) — that is", 0.0,
     "pan-cancer sort-fix move in 24-partition sds", "WRAP:"),
    ("sffc_pc_pairs_ge_move", r"since (\d+) of the \d+ pairs among the 24 partitions differ by more",
     0.0, "pan-cancer partition pairs differing by more than the sort-fix move", "WRAP:"),
    ("sffc_pc_n_pairs", r"since \d+ of the (\d+) pairs among the 24 partitions differ by more",
     0.0, "pairs of pan-cancer 24 partitions", "WRAP:"),
    ("sffc_ns_move_sds_small", r"the same fix moved the estimate by (\d\.\d) partition standard", 0.0,
     "NSCLC sort-fix move in six-partition sds, contrast sentence", "WRAP:"),
    ("sffc_pc_move_sds_small", r"partition standard deviations and here by (\d\.\d), a factor", 0.0,
     "pan-cancer sort-fix move in five-partition sds, contrast sentence", "WRAP:"),
    ("sffc_contrast_small", r"a factor of (\d\.\d) \(\d\.\d in 24-partition units\)", 0.0,
     "NSCLC-over-pan-cancer sort-fix move, six/five-partition units", "WRAP:"),
    ("sffc_contrast_24", r"a factor of \d\.\d \((\d\.\d) in 24-partition units\)", 0.0,
     "NSCLC-over-pan-cancer sort-fix move, 24-partition units", "WRAP:"),
    ("sffc_pc_sites_tied", r"(\d{3}) of 619 holding", 0.0,
     "pan-cancer sites sharing a size with another site", "WRAP:"),
    ("sffc_pc_pct_tied", r"of 619 holding (\d+\.\d)% of its patients", 0.0,
     "pan-cancer patients in tied sites, percent", "WRAP:"),
    ("sffc_ns_sites_tied", r"against (\d{2}) of 68 holding", 0.0,
     "NSCLC sites sharing a size with another site", "WRAP:"),
    ("sffc_ns_pct_tied", r"of 68 holding (\d+\.\d)% in NSCLC", 0.0,
     "NSCLC patients in tied sites, percent", "WRAP:"),
    ("sffc_pc_moved", r"the fix moved (\d,\d{3}) of 7,168 pan-cancer", 0.0,
     "pan-cancer patients the sort fix moved to another fold", "WRAP:"),
    ("sffc_pc_pct_moved", r"of 7,168 pan-cancer patients \((\d+\.\d)%\)", 0.0,
     "pan-cancer patients the sort fix moved, percent", "WRAP:"),
    ("sffc_ns_moved", r"against (\d{3}) of 944 \(\d+\.\d%\) in NSCLC", 0.0,
     "NSCLC patients the sort fix moved to another fold", "WRAP:"),
    ("sffc_ns_pct_moved", r"of 944 \((\d+\.\d)%\) in NSCLC", 0.0,
     "NSCLC patients the sort fix moved, percent", "WRAP:"),
    ("sffc_pc_ari", r"two partitions is (0\.\d+) pan-cancer and", 0.0,
     "adjusted Rand index, pan-cancer frozen vs corrected partition", "WRAP:"),
    ("sffc_ns_ari", r"pan-cancer and (0\.\d+) in NSCLC\. In pan-cancer the fix", 0.0,
     "adjusted Rand index, NSCLC frozen vs corrected partition", "WRAP:"),
    # The Linux site-control median next to the macOS one it differs from
    # (Limitation 8 and Negative controls both say why).
    ("e16_ns_registered_f100_median", r"full-cohort NSCLC median reads (\d\.\d{3})", 5e-4,
     "E16 NSCLC full-cohort site AUROC, Linux", "WRAP:"),
    ("nsclc_site_auroc", r"where the macOS run gives (\d\.\d{3})", 5e-4,
     "frozen NSCLC site AUROC median, macOS", "WRAP:"),
    ("nsclc_site_auroc", r"reads \d\.\d{3} here and (\d\.\d{3}) above", 5e-4,
     "frozen NSCLC site AUROC median, macOS, Negative controls", "WRAP:"),
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
    # Floors 5e-7 since 2026-09-24 (they were 0.005, a hundred times the
    # printed half-width): at 0.005 the note's macOS score sd, printed 0.0211,
    # passed against 0.021154 -- which rounds to 0.0212.
    ("sortaudit_macos_max_diff",
     r"difference between the two orderings is ([\d.]+)", 5e-7,
     "max ssGSEA score difference between tie conventions", None),
    ("sortaudit_macos_score_sd",
     r"per-signature score standard deviation of ([\d.]+)", 5e-7,
     "mean per-signature ssGSEA score sd", None),
    # E20: the Linux values beside them, so the prose names both platforms.
    ("sortaudit_hpc4_max_diff",
     r"\(macOS; ([\d.]+) and\s+[\d.]+ on Linux", 5e-5,
     "max ssGSEA score difference between tie conventions, Linux", r"WRAP:\(macOS; "),
    ("sortaudit_hpc4_score_sd",
     r"\(macOS; [\d.]+ and\s+([\d.]+) on Linux", 5e-5,
     "mean per-signature ssGSEA score sd, Linux", r"WRAP:\(macOS; "),
    ("sortaudit_hpc4_frac_sd",
     r"movement in the per-cell scores \(([\d.]+) on Linux\)", 5e-3,
     "max ssGSEA difference as a fraction of the score sd, Linux", "WRAP:movement in the per-cell"),

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

# Session 51: every number in the dry run's cut ladder and AI-variant tables.
# Line-scoped; one claim per table cell. See `_abstract_variants()`.
CLAIMS += [(key, prefix + "(" + _NUM + ")" + suffix, 0.5,
            f"dry-run abstract variant {key}, recomputed from the draft", None)
           for key, prefix, suffix in ABSTRACT_VARIANT_SPEC]

# Session 57: the dry run's prose re-quotations of those same tables. One
# claim per quotation, numbered so each has its own PATTERN REACH line even
# where two quote the same authority. See DRYRUN_PROSE_SPEC.
CLAIMS += [(key, prefix + "(" + _NUM + ")" + suffix, 0.5,
            f"dry-run prose re-quotes {key} (quotation {i + 1})", None)
           for i, (key, prefix, suffix) in enumerate(DRYRUN_PROSE_SPEC)]

# Session 57: A16's extent -- the patients whose signature scores describe
# tissue other than the slide's tumour -- against `_a16_extent()`, a live
# recount from the interim expression index. Each (key, what, shape) gives the
# claim regex (shape with {v} as the capture) and the self-test row (shape with
# {v} as the live value, pushed off it), so neither can go stale.
A16_EXTENT_SPEC: list[tuple[str, str, str]] = [
    ("a16_ns_not_tumour", "manuscript: NSCLC scored on non-tumour tissue",
     r"(In NSCLC, ){v}( of [\d,]+ patients were paired with a profile)"),
    ("a16_ns_adjacent_normal", "manuscript: NSCLC adjacent-normal profiles",
     r"(tissue: ){v}( adjacent normal, \d+ recurrence)"),
    ("a16_ns_recurrence", "manuscript: NSCLC recurrence profiles",
     r"(adjacent normal, ){v}( recurrence, and \d+ patient)"),
    ("a16_ns_normal_only", "manuscript: NSCLC patients with only a normal sample",
     r"(recurrence, and ){v}( patient with only a normal)"),
    ("a16_pan_averaged", "manuscript: pan-cancer profiles averaged",
     r"(Pan-cancer, ){v}( of [\d,]+ profiles averaged several samples)"),
    ("a16_pan_no_tumour", "manuscript: pan-cancer patients with no tumour sample",
     r"(profiles averaged several samples, and ){v}(\b)"),
    # Session 59: the results README's row for `nsclc_v3_tumouronly/`.
    ("a16_ns_not_tumour", "results README: NSCLC patients re-profiled tumour-only",
     r"(only, so the ){v}( patients whose frozen profile was an adjacent normal)"),
    ("a16_ns_adjacent_normal", "results README: NSCLC adjacent-normal profiles",
     r"(was an adjacent normal \(){v}(\), a recurrence)"),
    ("a16_ns_recurrence", "results README: NSCLC recurrence profiles",
     r"(a recurrence \(){v}(\) or normal-only)"),
    ("a16_ns_normal_only", "results README: NSCLC normal-only patients",
     r"(or normal-only \(){v}(, dropped\))"),
    ("a16_pan_no_tumour", "results README: pan-cancer patients without a tumour sample",
     r"(patients \(the ){v}( without a tumour sample dropped\))"),
    ("a16_ns_not_tumour", "README: NSCLC scored on non-tumour tissue",
     r"(`14-SCIENCE-AUDIT\.md`: ){v}( of [\d,]+ NSCLC and)"),
    ("nsclc_n", "README: NSCLC cohort in the A16 sentence",
     r"(`14-SCIENCE-AUDIT\.md`: \d+ of ){v}( NSCLC and)"),
    ("a16_pan_total", "README: pan-cancer profiles not a single tumour sample",
     r"(NSCLC and ){v}( of [\d,]+ pan-TCGA patients)"),
    ("pancancer_n", "README: pan-cancer cohort in the A16 sentence",
     r"(NSCLC and [\d,]+ of ){v}( pan-TCGA patients)"),
    ("a16_ns_not_tumour", "audit: NSCLC scored on non-tumour tissue",
     r"(in extent \(){v}( of [\d,]+ NSCLC, )"),
    ("nsclc_n", "audit: NSCLC cohort in the A16 row",
     r"(in extent \(\d+ of ){v}( NSCLC, )"),
    ("a16_pan_total", "audit: pan-cancer profiles not a single tumour sample",
     r"(NSCLC, ){v}( of [\d,]+ pan-TCGA\))"),
    ("pancancer_n", "audit: pan-cancer cohort in the A16 row",
     r"(NSCLC, [\d,]+ of ){v}( pan-TCGA\))"),
    # Session 58: the audit's A16 narrative block, which states the breakdown
    # of the counts above. Its two-digit figures never reached the coverage
    # list (`^\d{1,2}$` is dropped as structural), so they were found by
    # reading the block, not by the sweep. "Two" and "one" in its NSCLC half
    # are words and stay unbound.
    ("pancancer_n", "audit block: pan-cancer cohort",
     r"(takes the mean\. Of the ){v}( patients,)"),
    ("a16_pan_single_tumour", "audit block: pan-cancer single tumour profiles",
     r"(\b){v}( have exactly one TOIL sample and it is the tumour)"),
    ("a16_pan_averaged", "audit block: pan-cancer profiles averaged",
     r"(it is the tumour\. ){v}( are averages of)"),
    ("a16_pan_tumour_normal", "audit block: pan-cancer tumour plus adjacent normal",
     r"(several samples \(){v}( tumour plus adjacent normal)"),
    ("a16_pan_tumour_recurrence", "audit block: pan-cancer tumour plus recurrence",
     r"(adjacent normal, ){v}( tumour plus recurrence)"),
    ("a16_pan_tumour_metastasis", "audit block: pan-cancer tumour plus metastasis",
     r"(\b){v}( tumour plus metastasis, \d+ others\))"),
    ("a16_pan_multi_other", "audit block: pan-cancer other averaged mixes",
     r"(tumour plus metastasis, ){v}( others\))"),
    ("a16_pan_no_tumour", "audit block: pan-cancer patients with no tumour sample",
     r"(others\), and ){v}( have no tumour sample at all)"),
    ("a16_pan_normal_only", "audit block: pan-cancer adjacent normal only",
     r"(\(){v}( adjacent normal only, \d+ recurrence only)"),
    ("a16_pan_recurrence_only", "audit block: pan-cancer recurrence only",
     r"(adjacent normal only, ){v}( recurrence only, )"),
    ("a16_pan_metastasis_only", "audit block: pan-cancer metastasis only",
     r"(recurrence only, ){v}( metastasis only\))"),
    ("a16_pan_type_thca", "audit block: most affected type THCA",
     r"(affected types are THCA \(){v}(\), KIRC)"),
    ("a16_pan_type_kirc", "audit block: most affected type KIRC",
     r"(\), KIRC \(){v}(\), PRAD)"),
    ("a16_pan_type_prad", "audit block: most affected type PRAD",
     r"(\), PRAD \(){v}(\), LUAD and LIHC)"),
    ("a16_pan_type_luad", "audit block: most affected type LUAD",
     r"(LUAD and LIHC \(){v}( each\))"),
    ("a16_pan_type_lihc", "audit block: most affected type LIHC",
     r"(LUAD and LIHC \(){v}( each\))"),
    ("a16_ns_single_tumour", "audit block: NSCLC single tumour profiles",
     r"(column order\. For ){v}( of [\d,]+ patients that is the only sample)"),
    ("nsclc_n", "audit block: NSCLC cohort, single-sample sentence",
     r"(column order\. For [\d,]+ of ){v}( patients that is the only sample)"),
    ("a16_ns_tumour_normal", "audit block: NSCLC tumour plus normal",
     r"(Of the ){v}( with a tumour and a normal sample)"),
    ("a16_ns_adjacent_normal", "audit block: NSCLC received the normal profile",
     r"(with a tumour and a normal sample, ){v}( received the)"),
    ("a16_ns_not_tumour", "audit block: NSCLC scored on non-tumour tissue",
     r"(So ){v}( of [\d,]+ NSCLC signature)"),
    ("nsclc_n", "audit block: NSCLC cohort, closing sentence",
     r"(So [\d,]+ of ){v}( NSCLC signature)"),
]

# Session 58: the five types the audit's A16 block names as most affected,
# each a live count in `_a16_extent()`.
A16_TOP_TYPES = ("THCA", "KIRC", "PRAD", "LUAD", "LIHC")


def _a16_claim(shape: str) -> str:
    """The claim regex for a shape "(prefix){v}(suffix)": the number is group 1.

    The shape's two outer groups exist for the self-test's \\g<1>/\\g<2>; left
    in the claim they would make the prefix group 1 and compare it instead.
    """
    cut = shape.index("){v}(")
    assert shape[0] == "(" and shape[-1] == ")", shape
    return shape[1:cut] + r"([\d,]+)" + shape[cut + 5:-1]


CLAIMS += [(key, _a16_claim(shape), 0.5, f"A16 extent, {what}", None)
           for key, what, shape in A16_EXTENT_SPEC]

# Session 58: small counts the coverage list never showed. `--coverage` drops
# `^\d{1,2}$` as structural, which is right for list markers and wrong for
# "16 of 16 signatures exceeded their null" in the conference abstract. A
# mutation sweep of the two-digit literals sharing a line with a checked claim
# in the five outward-facing documents moved 111 of them behind a caught
# positive control; these are the ones that state a result or a count the
# results support. Same (key, what, shape) form as A16_EXTENT_SPEC, and the same
# generated self-test row, pushed 2 off the rounded live value.
SMALL_COUNT_SPEC: list[tuple[str, str, str]] = [
    # -- the conference abstract draft (07-ABSTRACT-DRAFT.md) --
    ("pancancer_sigs_beat", "abstract: signatures beating null, both cohorts (pan-cancer)",
     r"(; ){v}( of \d+ signatures exceeded their null in both)"),
    ("nsclc_sigs_beat", "abstract: signatures beating null, both cohorts (NSCLC)",
     r"(; ){v}( of \d+ signatures exceeded their null in both)"),
    ("pancancer_n_sigs", "abstract: signature family size",
     r"(; \d+ of ){v}( signatures exceeded their null in both)"),
    ("other_beat", "abstract: non-immune outcome winners",
     r"(MSigDB-immune versus ){v}( of \d+ others, )"),
    ("other_n", "abstract: non-immune signatures",
     r"(MSigDB-immune versus \d+ of ){v}( others, )"),
    ("immune_beat", "immune outcome winners, in the 'did (0 of 6' form",
     r"(signature did \(){v}( of \d+)"),
    ("immune_n", "immune signatures, in the 'did (0 of 6' form",
     r"(signature did \(\d+ of ){v}(\b)"),
    ("pancancer_n_types", "pan-TCGA cancer types, in the '(n=7,168, 31 types)' form",
     r"(\(n=[\d,]+, ){v}( types\))"),
    ("immune_beat", "abstract draft note: immune winners, hyphenated",
     r"(it is ){v}(-of-\d+ vs \d+-of-\d+, Mann)"),
    ("immune_n", "abstract draft note: immune signatures, hyphenated",
     r"(it is \d+-of-){v}( vs \d+-of-\d+, Mann)"),
    ("other_beat", "abstract draft note: non-immune winners, hyphenated",
     r"(it is \d+-of-\d+ vs ){v}(-of-\d+, Mann)"),
    ("other_n", "abstract draft note: non-immune signatures, hyphenated",
     r"(it is \d+-of-\d+ vs \d+-of-){v}(, Mann)"),
    # -- the manuscript (09-PAPER-DRAFT.md) --
    ("pancancer_sigs_beat", "manuscript abstract: all signatures exceed (pan-cancer)",
     r"(in NSCLC, all ){v}( signatures exceeding their)"),
    ("nsclc_sigs_beat", "manuscript abstract: all signatures exceed (NSCLC)",
     r"(in NSCLC, all ){v}( signatures exceeding their)"),
    ("other_beat", "non-immune outcome winners, 'versus 5 of 10 other signatures'",
     r"(MSigDB-immune versus ){v}( of \d+ other signatures)"),
    ("other_n", "non-immune signatures, 'versus 5 of 10 other signatures'",
     r"(MSigDB-immune versus \d+ of ){v}( other signatures)"),
    ("immune_beat", "manuscript Results: immune outcome winners",
     r"(\(\d+\), ){v}( of \d+ immune signatures beat their null)"),
    ("immune_n", "manuscript Results: immune signatures",
     r"(\(\d+\), \d+ of ){v}( immune signatures beat their null)"),
    ("other_beat", "manuscript Results: non-immune outcome winners",
     r"(against ){v}( of \d+ others \(Mann)"),
    ("other_n", "manuscript Results: non-immune signatures",
     r"(against \d+ of ){v}( others \(Mann)"),
    ("nsclc_n_sites", "manuscript Methods: NSCLC sites, tied-size sentence",
     r"(are common — \d+ of the ){v}( NSCLC sites share)"),
    ("nsclc_n_sites", "manuscript Limitation 8: NSCLC sites, tied-size sentence",
     r"(\d+ of the ){v}( NSCLC sites share a size with at least one other site)"),
    ("nsclc_n_sites", "manuscript: NSCLC sites in the site-count ratio",
     r"(for [\d,]+ sites against ){v}(\)\.)"),
    ("gene_set_size_angiogenesis", "manuscript: the smallest set's size",
     r"(the smallest, ){v}(-gene angiogenesis set)"),
    ("nsclc_sigs_beat", "manuscript: acceptance check, signatures beating null",
     r"(0\.3763\], ){v}(/\d+ beating null\), confirming)"),
    ("nsclc_n_sigs", "manuscript: acceptance check, family size",
     r"(0\.3763\], \d+/){v}( beating null\), confirming)"),
    ("a7_beats_ssgsea_unc", "manuscript: ssGSEA uncorrected, signatures beating null",
     r"(under mean-z, with ){v}( of \d+ signatures still)"),
    ("a7_beats_meanz_unc", "manuscript: mean-z uncorrected, signatures beating null",
     r"(beating their null \(){v}(/\d+ under mean-z\))"),
    ("nsclc_bh_beat", "manuscript: signatures beating null after BH (NSCLC)",
     r"(neither headline: ){v}( of \d+ signatures still exceeded their null on the index)"),
    ("pancancer_bh_beat", "manuscript: signatures beating null after BH (pan-cancer)",
     r"(neither headline: ){v}( of \d+ signatures still exceeded their null on the index)"),
    ("cap_nsclc_lsv_sites", "manuscript: Control C usable NSCLC sites",
     r"(label across ){v}( usable sites and)"),
    ("e16_pan_registered_quarter_sites_min", "manuscript: E16 quarter, fewest evaluable sites",
     r"(patients, ){v}( to \d+ evaluable sites)"),
    ("e16_pan_registered_quarter_sites_max", "manuscript: E16 quarter, most evaluable sites",
     r"(patients, \d+ to ){v}( evaluable sites)"),
    ("sortaudit_beats", "manuscript Limitation 8: NSCLC corrected ordering, signatures beating null",
     r"(differ by as much\. All ){v}( of \d+)"),
    ("pcsortaudit_beats", "manuscript Limitation 8: pan-cancer corrected ordering, signatures beating null",
     r"(analysis derived\. All ){v}( of \d+ signatures still exceed)"),
    ("a7_beats_ssgsea_unc", "manuscript Limitation 9: ssGSEA uncorrected, signatures beating null",
     r"(vs \+[\d.]+, ){v}(/\d+ vs \d+/\d+ beating)"),
    ("a7_beats_meanz_unc", "manuscript Limitation 9: mean-z uncorrected, signatures beating null",
     r"(vs \+[\d.]+, \d+/\d+ vs ){v}(/\d+ beating)"),
    # -- the poster --
    ("pancancer_rel_gap_fold", "poster: pan-cancer reliability-gap fold",
     r"(<b>){v}(-fold and [\d.]+-fold</b>)"),
    ("small_panel_k10_sets", "poster: sets evaluable at k=10",
     r'(<td class="hi">10</td><td class="hi">[\d.]+</td><td>){v}(</td>)'),
    ("small_panel_k20_sets", "poster: sets evaluable at k=20",
     r"(<td>20</td><td>[\d.]+</td><td>){v}(</td>)"),
    ("small_panel_k40_sets", "poster: sets evaluable at k=40",
     r"(<td>40</td><td>[\d.]+</td><td>){v}(</td>)"),
    ("small_panel_k80_sets", "poster: sets evaluable at k=80",
     r"(<td>80</td><td>[\d.]+</td><td>){v}(</td>)"),
    ("small_panel_k160_sets", "poster: sets evaluable at k=160",
     r"(<td>160</td><td>[\d.]+</td><td>){v}(</td>)"),
    ("pancancer_sigs_beat", "poster lead: signatures beating null, both cohorts (pan-cancer)",
     r"(\b){v}( of \d+ signatures exceeded their own null in both cohorts)"),
    ("nsclc_sigs_beat", "poster lead: signatures beating null, both cohorts (NSCLC)",
     r"(\b){v}( of \d+ signatures exceeded their own null in both cohorts)"),
    ("pancancer_n_sigs", "poster lead: signature family size",
     r"(\b\d+ of ){v}( signatures exceeded their own null in both cohorts)"),
    ("pancancer_bh_beat", "poster: signatures beating null after BH",
     r"(q&nbsp;&lt;&nbsp;0\.05: ){v}(/\d+;)"),
    ("pancancer_n_sigs", "poster: BH family size",
     r"(q&nbsp;&lt;&nbsp;0\.05: \d+/){v}(;)"),
    ("a10ss_beats_ssgsea_unc", "poster: pinned ssGSEA, signatures beating null",
     r"(\+[\d.]+ pinned, ){v}(/\d+ signatures)"),
    # -- the cover letter --
    ("pancancer_sigs_beat", "cover letter: signatures beating null (pan-cancer)",
     r"(that ){v}( of \d+ signatures exceed their null\.)"),
    ("nsclc_sigs_beat", "cover letter: signatures beating null (NSCLC)",
     r"(that ){v}( of \d+ signatures exceed their null\.)"),
    ("nsclc_sigs_beat", "cover letter: A7 leaves the NSCLC count untouched",
     r"(significance and ){v}(/\d+ beating null are untouched)"),
    # -- session 59: the "N of M signatures" word form, swept on every line that
    # writes it (12, line-scoped, across the 28 scanned documents). Three were
    # unbound: the S1 caption's BH count, and E8's omega count in the audit and
    # in Limitation 9. The "16" of "15 or 16 of 16" is the alpha
    # reconstruction's count, which `omega_reliability/summary.json` does not
    # store, so it stays unbound; so do the "of 16" design denominators. --
    ("pancancer_bh_beat", "S1 caption: signatures beating null after BH (pan-cancer)",
     r"(cohorts give ){v}( of \d+ signatures beating their null after)"),
    ("nsclc_bh_beat", "S1 caption: signatures beating null after BH (NSCLC)",
     r"(cohorts give ){v}( of \d+ signatures beating their null after)"),
    ("e8_beats_null_omega", "audit E8: signatures above their null under omega",
     r"(\b){v}( of \d+ signatures still exceed their null under omega)"),
    ("e8_beats_null_omega", "audit E8: the omega half of '15 or 16 of 16'",
     r"(agree on direction and on ){v}( or \d+ of \d+)"),
    ("e8_beats_null_omega", "manuscript Limitation 9: the omega half of '15 or 16 of 16'",
     r"(same ){v}( or \d+ of \d+ signatures above their null)"),
    # -- session 59: the two-digit literals sharing a line with a checked claim
    # in the deposited documents and the audit, swept behind caught controls
    # (66 moved in the audit, 62 UNCAUGHT; 23 in nine other documents, 22
    # UNCAUGHT -- most of them dates, design constants and run records). These
    # are the ones that state a result. --
    ("pancancer_n_types", "deposited README headline table: cancer types",
     r"(\(n = [\d,]+, ){v}( types\))"),
    ("nsclc_n_sites", "S-captions: NSCLC sites in the Control C site count",
     r"(pan-TCGA, \d+ of ){v}( NSCLC\))"),
    ("pancancer_n_sigs", "MANIFEST: S1 rows per cohort (pan-cancer)",
     r"(all ){v}( rows per cohort belong in the supplement)"),
    ("nsclc_n_sigs", "MANIFEST: S1 rows per cohort (NSCLC)",
     r"(all ){v}( rows per cohort belong in the supplement)"),
    ("pancancer_n_types", "NextGen abstract: cancer types",
     r"(pooled across ){v}( diseases)"),
    ("nsclc_n_sites", "audit A5: NSCLC sites",
     r"(NSCLC, n=[\d,]+, ){v}( sites:)"),
    ("sens_ns_reference_beats", "audit E1 table: NSCLC within-type row, signatures above null",
     r"(\| within type \(registered\) \| 0\.2922 \[[^\]]+\] \| ){v}(/16 \|)"),
    ("sens_ns_global_axis_beats", "audit E1 table: NSCLC global-axis row, signatures above null",
     r"(\| global, both types together \| [^|]+ \| ){v}(/16 \|)"),
    ("pcsortaudit_beats", "audit E1 table: pan-cancer within-type row, signatures above null",
     r"(\| within type \(registered\) \| 0\.2968 \[[^\]]+\] \| ){v}(/16 \|)"),
    ("sens_pan_global_axis_beats", "audit E1 table: pan-cancer global-axis row, signatures above null",
     r"(\| global, all \d+ types together \| [^|]+ \| ){v}(/16 \|)"),
    ("pancancer_n_types", "audit E1 table: pan-cancer global-axis row, cancer types",
     r"(\| global, all ){v}( types together \|)"),
    ("sens_ns_unmatched_null_beats", "audit E3 row: signatures above null, size-only null",
     r"(\| E3, null matched on size only \| [^|]+ \| ){v}(/16 \|)"),
    ("sens_ns_k3_beats", "audit E5 row: signatures above null, 3 folds",
     r"(\| E5, 3 folds \| [^|]+ \| ){v}(/16 \|)"),
    ("sens_ns_k10_beats", "audit E5 row: signatures above null, 10 folds",
     r"(\| E5, 10 folds \| [^|]+ \| ){v}(/16 \|)"),
    ("e6_meanz_unc_beats", "audit E6 table: mean-z uncorrected, signatures above null",
     r"(\| mean-z, uncorrected \(control, same run\) \| [^|]+ \| ){v}(/16 \|)"),
    ("e6_plage_unc_beats", "audit E6 table: PLAGE uncorrected, signatures above null",
     r"(\| PLAGE, uncorrected \| [^|]+ \| ){v}(/16 \|)"),
    ("a10ss_beats_ssgsea_unc", "audit E6 table: pinned ssGSEA uncorrected, signatures above null",
     r"(\| ssGSEA, uncorrected \(stored, pinned ordering\) \| [^|]+ \| ){v}(/16 \|)"),
    ("a10ss_beats_ssgsea_unc", "audit A-item table: A10, signatures above null either way",
     r"(half-width; ){v}(/16 either way;)"),
    # -- session 59: the results README's two-digit literals, swept behind a
    # caught control; these state results (the rest are dates, run records,
    # partition and fold counts by design, mode tallies and list numbers). --
    ("pancancer_n_types", "results README: pancancer_v3 row, cancer types",
     r"(The pan-TCGA headline, n=[\d,]+, ){v}( types\.)"),
    ("pancancer_n_partitions_used", "results README: A5 pan-cancer, usable partitions",
     r"(partitions requested, \*\*){v}( usable, \d+ refused\*\*)"),
    ("pancancer_n_partitions_refused", "results README: A5 pan-cancer, refused partitions",
     r"(partitions requested, \*\*\d+ usable, ){v}( refused\*\*)"),
    ("pcsortaudit_beats", "results README: pan-cancer corrected ordering, signatures beating null",
     r"(0\.31872055896523865\]`\. \*\*){v}(/16 signatures still beat their null)"),
    ("sortaudit_beats", "results README: NSCLC corrected ordering, signatures beating null",
     r"(anticipated a move of this size\. \*\*){v}(/16 signatures still beat their null)"),
    ("a10ss_beats_ssgsea_unc", "results README: A10 pinned ssGSEA, signatures above null",
     r"(0\.21018161815588407\]\*\*, ){v}(/16 signatures above their null, against the pre-fix)"),
    ("a7_beats_ssgsea_unc", "results README: A7 pre-fix ssGSEA, signatures above null",
     r"(0\.2159114981000094\]\*\*, also ){v}(/16\.)"),
    ("e6_meanz_unc_beats", "results README: E6 mean-z uncorrected, signatures above null",
     r"(Uncorrected: mean_z \+[\d.]+ \[[^\]]+\] ){v}(/16, plage)"),
    ("e6_plage_unc_beats", "results README: E6 PLAGE uncorrected, signatures above null",
     r"(plage \+[\d.]+ \[[^\]]+\] ){v}(/16 \(stored pinned ssGSEA)"),
    ("a10ss_beats_ssgsea_unc", "results README: E6 row, stored pinned ssGSEA, signatures above null",
     r"(stored pinned ssGSEA \+[\d.]+, ){v}(/16\))"),
    ("e6_meanz_recon_beats", "results README: E6 mean-z reconstruction, signatures above null",
     r"(Split-half reconstruction: mean_z \+[\d.]+ ){v}(/16, plage)"),
    ("e6_plage_recon_beats", "results README: E6 PLAGE reconstruction, signatures above null",
     r"(plage \*\*\+[\d.]+ ){v}(/16\*\*, against)"),
    ("a7_recon_ssgsea_beats", "results README: A7 ssGSEA reconstruction, signatures above null",
     r"(against ssGSEA's [−-][\d.]+ ){v}(/16, so)"),
    ("e8_beats_null_omega", "results README: E8 omega reconstruction, signatures above null",
     r"(0\.2554\], ){v}(/16 still above null)"),
    ("sens_ns_global_axis_beats", "results README: E1 NSCLC global axis, signatures above null",
     r"(global axis [−-][\d.]+ \[[^\]]+\], ){v}(/16, with \d+ of 16 nulls degenerate)"),
    ("sens_ns_global_degenerate", "results README: E1 NSCLC global axis, degenerate nulls",
     r"(/16, with ){v}( of 16 nulls degenerate)"),
    ("sens_ns_unmatched_null_beats", "results README: E3 size-only null, signatures above null",
     r"(size-only null [\d.]+ \[[^\]]+\], ){v}(/16 \()"),
    ("sens_ns_k3_beats", "results README: E5 3 folds, signatures above null",
     r"(3 folds [\d.]+ \[[^\]]+\], ){v}(/16 \()"),
    ("sens_ns_k10_beats", "results README: E5 10 folds, signatures above null",
     r"(10 folds [\d.]+ \[[^\]]+\], ){v}(/16 \()"),
    ("sens_pan_global_axis_beats", "results README: E1 pan-cancer global axis, signatures above null",
     r"(0\.4353\], ){v}(/16 above null, \d+ degenerate nulls)"),
    ("splitdet_macos_sites_tied", "results README: fold census, NSCLC tied sites",
     r"(against NSCLC's ){v}( of \d+ \([\d.]+%\), and the fix moved)"),
    ("nsclc_n_sites", "results README: fold census, NSCLC sites",
     r"(against NSCLC's \d+ of ){v}( \([\d.]+%\), and the fix moved)"),
    ("nsclc_n_sites", "results README: fold census reference, NSCLC sites",
     r"(`split_determinism_darwin_arm64\.json` \(){v}( sites, \d+ tied, \d+ moved\))"),
    ("splitdet_macos_sites_tied", "results README: fold census reference, NSCLC tied sites",
     r"(`split_determinism_darwin_arm64\.json` \(\d+ sites, ){v}( tied, \d+ moved\))"),
    ("splitdet_macos_sites_tied", "results README: split-determinism row, NSCLC tied sites",
     r"(\*\*){v}( of \d+ NSCLC sites are tied on size\*\*)"),
    ("nsclc_n_sites", "results README: split-determinism row, NSCLC sites",
     r"(\*\*\d+ of ){v}( NSCLC sites are tied on size\*\*)"),
    ("fullsize_set_min_genes", "manuscript Results (B-34): fewest genes in a nominal-200 set",
     r"(Hallmark sets \(){v}(–\d+ genes after filtering\))"),
    ("fullsize_set_max_genes", "manuscript Results (B-34): most genes in a nominal-200 set",
     r"(Hallmark sets \(\d+–){v}( genes after filtering\))"),
    # -- the portal dry run --
    ("abstract_pct_of_limit", "dry run: abstract length as a percentage of the limit",
     r"(\*\*[\d,]+ / [\d,]+ \(){v}(%\), headroom)"),
    ("abstract_naive_excess_pct", "dry run: how much higher a naive len() reads",
     r"(reads [\d,]+ — ){v}(% higher)"),
]

CLAIMS += [(key, _a16_claim(shape), 0.5, f"small count, {what}", None)
           for key, what, shape in SMALL_COUNT_SPEC]

# Session 58: the operational checklists, swept by mutation behind a caught
# positive control. `15-CHECKLIST.md` and `16-YOUR-TASKS.md` are what the
# author works down, and their completed items quote the current science --
# partition sds, honest intervals, the tumour-only indices, the AI-variant
# arithmetic -- with nothing checking it. Records of past runs (PIDs, job ids,
# wall clocks of one run, file counts at a dated commit) stay unbound. Each
# entry is (key, what, shape, printed format); the claim holds the document to
# its printed precision, and the self-test row moves the value four units in
# its last printed place.
CHECKLIST_SPEC: list[tuple[str, str, str, str]] = [
    # -- 15-CHECKLIST.md --
    ("nsclc_honest_ci_lo", "checklist: NSCLC honest CI, low",
     r"(Honest CI \[){v}(, [\d.]+\])", "{:.4f}"),
    ("nsclc_honest_ci_hi", "checklist: NSCLC honest CI, high",
     r"(Honest CI \[[\d.]+, ){v}(\])", "{:.4f}"),
    ("pancancer_bootstrap_se", "checklist: pan-cancer bootstrap SE",
     r"(\b){v}(; the reported interval is \*\*[\d.]+% too narrow)", "{:.4f}"),
    ("pancancer_honest_ci_lo", "checklist: pan-cancer honest CI, low",
     r"(\[){v}(, [\d.]+\]; partition choice accounts for)", "{:.4f}"),
    ("pancancer_honest_ci_hi", "checklist: pan-cancer honest CI, high",
     r"(\[[\d.]+, ){v}(\]; partition choice accounts for)", "{:.4f}"),
    ("pv24_nsclc_sd", "checklist: NSCLC 24-partition sd",
     r"(was computed\. NSCLC: sd ){v}(,)", "{:.4f}"),
    ("pv24_nsclc_understated_pct", "checklist: NSCLC 24-partition understatement",
     r"(interval understated by ){v}(%, share [\d.]+%\. Pan-cancer)", "{:.1f}"),
    ("pv24_nsclc_share_pct", "checklist: NSCLC 24-partition share",
     r"(%, share ){v}(%\. Pan-cancer: sd)", "{:.1f}"),
    ("pv24_pancancer_sd", "checklist: pan-cancer 24-partition sd",
     r"(%\. Pan-cancer: sd ){v}(,)", "{:.4f}"),
    ("pv24_pancancer_understated_pct", "checklist: pan-cancer 24-partition understatement",
     r"(understated by ){v}(%, share [\d.]+%\. The sds do not differ)", "{:.1f}"),
    ("pv24_pancancer_share_pct", "checklist: pan-cancer 24-partition share",
     r"(%, share ){v}(%\. The sds do not differ)", "{:.1f}"),
    ("cohort_n_ratio", "checklist: pan-cancer cohort as a multiple of NSCLC",
     r"(despite it being ){v}(× the size)", "{:.1f}"),
    ("nsclc_partition_sd_chi2_lo", "checklist: six-partition NSCLC sd interval, low",
     r"(the chi-square intervals are \[){v}(, [\d.]+\] and)", "{:.4f}"),
    ("nsclc_partition_sd_chi2_hi", "checklist: six-partition NSCLC sd interval, high",
     r"(the chi-square intervals are \[[\d.]+, ){v}(\] and)", "{:.4f}"),
    ("pancancer_partition_sd_chi2_lo", "checklist: five-partition pan-cancer sd interval, low",
     r"(\[){v}(, [\d.]+\], and the variance shares)", "{:.4f}"),
    ("pancancer_partition_sd_chi2_hi", "checklist: five-partition pan-cancer sd interval, high",
     r"(\[[\d.]+, ){v}(\], and the variance shares)", "{:.4f}"),
    ("nsclc_partition_share_lo_pct", "checklist: six-partition NSCLC share interval, low",
     r"(and the variance shares \[){v}(%, [\d.]+%\] and)", "{:.1f}"),
    ("nsclc_partition_share_hi_pct", "checklist: six-partition NSCLC share interval, high",
     r"(and the variance shares \[[\d.]+%, ){v}(%\] and)", "{:.1f}"),
    ("pancancer_partition_share_lo_pct", "checklist: five-partition pan-cancer share interval, low",
     r"(\[){v}(%, [\d.]+%\]\. Both pairs overlap)", "{:.1f}"),
    ("pancancer_partition_share_hi_pct", "checklist: five-partition pan-cancer share interval, high",
     r"(\[[\d.]+%, ){v}(%\]\. Both pairs overlap)", "{:.1f}"),
    ("cap_pan_lsv_sites", "checklist: Control C usable pan-cancer sites",
     r"(\(){v}( and \d+, which are Control C's \*usable\* sites\))", "{:,.0f}"),
    ("cap_nsclc_lsv_sites", "checklist: Control C usable NSCLC sites",
     r"(\(\d+ and ){v}(, which are Control C's \*usable\* sites\))", "{:,.0f}"),
    ("pancancer_n_sites", "checklist: pan-cancer sites that govern fold assignment",
     r"(govern fold assignment are ){v}( and \d+\.)", "{:,.0f}"),
    ("nsclc_n_sites", "checklist: NSCLC sites that govern fold assignment",
     r"(govern fold assignment are \d+ and ){v}(\.)", "{:,.0f}"),
    ("a7_ssgsea_unc", "checklist: A7 ssGSEA uncorrected contrast",
     r"(magnitude \(\+){v}( vs \+[\d.]+, \d+/\d+ vs \d+/\d+\))", "{:.3f}"),
    ("a7_meanz_unc", "checklist: A7 mean-z uncorrected contrast",
     r"(magnitude \(\+[\d.]+ vs \+){v}(, \d+/\d+ vs \d+/\d+\))", "{:.3f}"),
    ("a7_beats_ssgsea_unc", "checklist: A7 ssGSEA signatures beating null",
     r"(magnitude \(\+[\d.]+ vs \+[\d.]+, ){v}(/\d+ vs \d+/\d+\))", "{:,.0f}"),
    ("a7_beats_meanz_unc", "checklist: A7 mean-z signatures beating null",
     r"(magnitude \(\+[\d.]+ vs \+[\d.]+, \d+/\d+ vs ){v}(/\d+\))", "{:,.0f}"),
    ("a7_recon_ssgsea", "checklist: A7 reconstruction, magnitude of the reversed sign",
     r"(reconstruction reverses sign \(−){v}(\))", "{:.4f}"),
    ("a7_ssgsea_obs_inflation", "checklist: A7 disattenuation of the observed",
     r"(≈ [\d.]+× against ){v}(× for)", "{:.2f}"),
    ("splitdet_macos_moved", "checklist: patients the sort fix moved, macOS",
     r"(0 of \d+ patients differ, was ){v}(\))", "{:,.0f}"),
    ("nsclc_site_auroc", "checklist: NSCLC site AUROC in the poster-panel record",
     r"(site AUROC [\d.]+ / ){v}(; Δr)", "{:.3f}"),
    ("paper_abstract_body_words", "checklist: manuscript abstract body words",
     r"(unstructured paragraph of \*\*){v}( words plus a)", "{:,.0f}"),
    ("paper_abstract_significance_words", "checklist: manuscript Significance words",
     r"(\b){v}(-word Significance = \d+\*\*)", "{:,.0f}"),
    ("paper_abstract_total_words", "checklist: manuscript abstract total words",
     r"(-word Significance = ){v}(\*\*)", "{:,.0f}"),
    ("paper_structured_abstract_words", "checklist: folded structured abstract words",
     r"(version is ){v}(, which the manuscript now says)", "{:,.0f}"),
    ("sortaudit_macos_tied_rows", "checklist: A10 samples carrying rank ties",
     r"(itself contradicted\. All ){v}( samples carry rank ties)", "{:,.0f}"),
    ("sortaudit_macos_frac_sd", "checklist: A10 score move in score sds",
     r"(ssGSEA scores move by ){v}( of a score sd)", "{:.2f}"),
    ("nsclc_interval_understated_pct", "checklist: NSCLC six-partition understatement",
     r"(\b){v}(% too narrow, not the 27% once quoted)", "{:.1f}"),
    ("pancancer_bh_beat", "checklist summary: signatures surviving BH",
     r"(changes nothing: ){v}(/\d+ and \d+/\d+ survive)", "{:,.0f}"),
    ("pancancer_outcome_bh_beat", "checklist summary: outcome comparisons surviving BH",
     r"(changes nothing: \d+/\d+ and ){v}(/\d+ survive)", "{:,.0f}"),
    # -- submission/SUBMISSION-CHECKLIST.md --
    ("pancancer_relgap_size_rho_abs", "submission checklist: deposited-snapshot size-gap rho",
     r"(reproducing ρ = −){v}(, p = [\d.]+;)", "{:.3f}"),
    ("pancancer_relgap_size_p", "submission checklist: deposited-snapshot size-gap p",
     r"(reproducing ρ = −[\d.]+, p = ){v}(;)", "{:.4f}"),
    ("nsclc_honest_ci_lo", "submission checklist: NSCLC honest CI, low",
     r"(\[){v}(, [\d.]+\]\. The manuscript reports the measurement)", "{:.4f}"),
    ("nsclc_honest_ci_hi", "submission checklist: NSCLC honest CI, high",
     r"(\[[\d.]+, ){v}(\]\. The manuscript reports the measurement)", "{:.4f}"),
    ("pv24_nsclc_understated_pct", "submission checklist: NSCLC 24-partition understatement",
     r"(total uncertainty by ){v}(% in NSCLC and by)", "{:.1f}"),
    ("pv24_pancancer_understated_pct", "submission checklist: pan-cancer 24-partition understatement",
     r"(% in NSCLC and by ){v}(% in pan-TCGA)", "{:.1f}"),
    ("pancancer_isi", "submission checklist: the v2/v3 mean excess",
     r"(mean excess ){v}(, range)", "{:.6f}"),
    ("paper_abstract_body_words", "submission checklist: manuscript abstract body words",
     r"(measuring \*\*){v}( \+ \d+ = \d+ words\*\*)", "{:,.0f}"),
    ("paper_abstract_significance_words", "submission checklist: manuscript Significance words",
     r"(measuring \*\*\d+ \+ ){v}( = \d+ words\*\*)", "{:,.0f}"),
    ("paper_abstract_total_words", "submission checklist: manuscript abstract total words",
     r"(measuring \*\*\d+ \+ \d+ = ){v}( words\*\*)", "{:,.0f}"),
    ("paper_structured_abstract_words", "submission checklist: folded structured abstract words",
     r"(\b){v}( by the word measure the manuscript states)", "{:,.0f}"),
    ("paper_abstract_total_words", "submission checklist: abstract total against the ceiling",
     r"(\b){v}( sits exactly on the observed ceiling)", "{:,.0f}"),
    # -- 16-YOUR-TASKS.md --
    ("abstract_chars", "task list: the abstract's count with the author's decisions in",
     r"(The abstract counts ){v}( characters, \d+ under the limit)", "{:,.0f}"),
    ("abstract_headroom", "task list: headroom with the author's decisions in",
     r"(The abstract counts [\d,]+ characters, ){v}( under the limit)", "{:,.0f}"),
    ("a16_ns_not_tumour", "task list: A16 NSCLC not scored on tumour",
     r"(not always the tumour\.\*\* ){v}( of [\d,]+ NSCLC)", "{:,.0f}"),
    ("nsclc_n", "task list: NSCLC cohort in the A16 item",
     r"(not always the tumour\.\*\* \d+ of ){v}( NSCLC)", "{:,.0f}"),
    ("a16_pan_total", "task list: A16 pan-cancer not a single tumour sample",
     r"(and ){v}( of [\d,]+ pan-cancer patients\. With tumour-only)", "{:,.0f}"),
    ("pancancer_n", "task list: pan-cancer cohort in the A16 item",
     r"(and [\d,]+ of ){v}( pan-cancer patients\. With tumour-only)", "{:,.0f}"),
    ("sortaudit_isi", "task list: NSCLC corrected-ordering ISI before tumour-only",
     r"(index RISES from ){v}( to [\d.]+ \(still)", "{:.4f}"),
    ("a16_ns", "task list: NSCLC tumour-only ISI",
     r"(index RISES from [\d.]+ to ){v}( \(still)", "{:.4f}"),
    ("a16_pan", "task list: pan-cancer tumour-only ISI",
     r"(since finished: ){v}(, against [\d.]+ under the corrected ordering)", "{:.4f}"),
    ("pcsortaudit_isi", "task list: pan-cancer corrected-ordering ISI",
     r"(since finished: [\d.]+, against ){v}( under the corrected ordering)", "{:.4f}"),
    ("sortaudit_macos_tied_rows", "task list: A10 samples carrying rank ties",
     r"(It was implicated: ){v}( of \d+ samples)", "{:,.0f}"),
    ("sortaudit_macos_frac_sd_pct", "task list: A10 score move as a percentage of an sd",
     r"(scores move by ){v}(% of a score sd)", "{:.1f}"),
    ("pancancer_bh_beat", "task list summary: signatures surviving BH",
     r"(BH-FDR \(){v}(/\d+ and \d+/\d+ both survive)", "{:,.0f}"),
    ("pancancer_outcome_bh_beat", "task list summary: outcome comparisons surviving BH",
     r"(BH-FDR \(\d+/\d+ and ){v}(/\d+ both survive)", "{:,.0f}"),
    ("nsclc_interval_understated_pct", "task list summary: NSCLC six-partition understatement",
     r"(partition variance \(interval ){v}(% too narrow, not 27%\))", "{:.1f}"),
]

CLAIMS += [(key, _a16_claim(shape).replace(r"([\d,]+)", r"(\d[\d,]*(?:\.\d+)?)", 1), 0.0,
            f"checklist figure, {what}", None)
           for key, what, shape, _fmt in CHECKLIST_SPEC]

# Session 59: the audit's A16 NSCLC paragraph and the two other places that
# quote its figures (the audit's A-item table, the results README's row for
# `nsclc_v3_tumouronly/`). Same form and same generated self-test row as
# CHECKLIST_SPEC: held to printed precision, moved four units in the last place.
A16_NARRATIVE_SPEC: list[tuple[str, str, str, str]] = [
    ("a16_ns_shift", "audit A16: NSCLC tumour-only shift",
     r"(\+){v}( \(about seven 24-partition)", "{:.4f}"),
    ("a16_ns_shift", "audit A-item table: NSCLC tumour-only shift",
     r"(0\.4079\] \(\+){v}(, job 130002\))", "{:.4f}"),
    ("a16_pan_fall", "audit A-item table: pan-cancer tumour-only fall",
     r"(0\.3149\] \([−-]){v}(, job 130155\))", "{:.4f}"),
    ("a16_ns_shift", "results README: NSCLC tumour-only shift",
     r"(a shift of \+){v}(, with \d+ of 16 per-signature)", "{:.4f}"),
    ("a16_ns_beats", "audit A16: signatures still beating null, tumour-only",
     r"(overlap\)\. ){v}(/16 signatures still beat)", "{:,.0f}"),
    ("a16_ns_beats", "results README: signatures above null, tumour-only",
     r"(\]\*\*, all ){v}( signatures above their null, against the corrected)", "{:,.0f}"),
    ("a16_ns_sigs_rise", "audit A16: per-signature excesses that rise",
     r"(their null, and ){v}( of 16 per-signature excesses rise)", "{:,.0f}"),
    ("a16_ns_sigs_rise", "results README: per-signature excesses moving up",
     r"(, with ){v}( of 16 per-signature excesses moving up)", "{:,.0f}"),
    ("a16_ns_tgf_rise", "audit A16: the largest rise, TGF-beta",
     r"(TGF-β, \+){v}(; angiogenesis falls by)", "{:.3f}"),
    ("a16_ns_angio_fall", "audit A16: angiogenesis's fall",
     r"(angiogenesis falls by ){v}(\)\. Adjacent-normal)", "{:.3f}"),
    ("a16_ns_dmae", "audit A16: tumour-only Delta-MAE",
     r"(Δ-MAE barely moves \(){v}()", "{:.5f}"),
    ("a16_ns_outcome_beat", "audit A16: outcome comparisons beating null, tumour-only",
     r"(\]\): ){v}(/\d+ against \d+/\d+\.)", "{:,.0f}"),
    ("hpc4_sortaudit_dmae", "audit A16: the HPC4 corrected-ordering Delta-MAE it is compared with",
     r"(against ){v}(\)\. The outcome arm gains one signature)", "{:.5f}"),
    ("hpc4_sortaudit_dmae", "audit A12 table: HPC4 unpinned Delta-MAE",
     r"(\| HPC4 unpinned vs macOS unpinned \(A12 above\) \| [^|]+ \| [^|]+ \| ){v}( vs [\d.]+ \|)", "{:.5f}"),
    ("sortaudit_outcome_beat", "audit A16: outcome comparisons beating null, corrected ordering",
     r"(\]\): \d+/\d+ against ){v}(/\d+\.)", "{:,.0f}"),
]

CLAIMS += [(key, _a16_claim(shape).replace(r"([\d,]+)", r"(\d[\d,]*(?:\.\d+)?)", 1), 0.0,
            f"audit A16 figure, {what}", None)
           for key, what, shape, _fmt in A16_NARRATIVE_SPEC]

# B-21, 2026-09-24 (after session 60): the manuscript's Results summaries.
# Five passages and the site control's stress test were lifted word for word
# into SUPPLEMENTARY-NOTE-2.md to bring the counted body under 5,000 words, and
# a summary was left in each place. Every figure a summary restates was checked
# in the original text and is still checked there, in the note; these bind the
# restatements to the same authorities, at printed precision, each with a
# generated self-test row. Shapes are each literal's own line context.
B21_SUMMARY_SPEC: list[tuple[str, str, str, str]] = [
    ("e6_plage_unc", "PLAGE uncorrected contrast",
     r"(as ssGSEA does \(\+){v}( \[)", "{:.4f}"),
    ("e6_plage_recon", "PLAGE scorer-agnostic reconstruction",
     r"(reconstruction \(\+){v}(\), so the reversal)", "{:.4f}"),
    ("nsclc_label_site_variance_pct", "NSCLC label-side site R^2",
     r"(same control attributes ){v}(% of the label to site)", "{:.1f}"),
    ("lsv_nsclc_site_given_type_pct", "NSCLC site-given-type R^2",
     r"(the median, ){v}(% once cancer type is conditioned on)", "{:.1f}"),
    ("e14_nsclc_median_pct", "NSCLC Control C calibrated median",
     r"(conditioned on, and ){v}(% under the)", "{:.1f}"),
    ("pv24_nsclc_sd", "NSCLC 24-partition sd",
     r"(0\.0089 over six partitions and ){v}( over)", "{:.4f}"),
    ("nsclc_honest_ci_lo", "NSCLC partition-widened interval, low",
     r"(\[){v}(, [\d.]+\]: the reported interval is 4\.5%)", "{:.4f}"),
    ("nsclc_honest_ci_hi", "NSCLC partition-widened interval, high",
     r"(\[[\d.]+, ){v}(\]: the reported interval is 4\.5%)", "{:.4f}"),
    ("pv24_nsclc_understated_pct", "NSCLC 24-partition understatement",
     r"(4\.5% too narrow \(){v}(% at 24\))", "{:.1f}"),
    ("pv24_nsclc_share_pct", "NSCLC 24-partition share",
     r"(8\.4% of total variance \(){v}(% at 24\))", "{:.1f}"),
    ("pv24_pancancer_sd", "pan-cancer 24-partition sd",
     r"(could be evaluated and ){v}( over 24,)", "{:.4f}"),
    ("pancancer_bootstrap_se", "pan-cancer bootstrap SE",
     r"(\b){v}(\. The interval widens from)", "{:.4f}"),
    ("pancancer_reported_ci_lo", "pan-cancer reported interval, low",
     r"(widens from \[){v}(, [\d.]+\] to \[)", "{:.4f}"),
    ("pancancer_reported_ci_hi", "pan-cancer reported interval, high",
     r"(widens from \[[\d.]+, ){v}(\] to \[)", "{:.4f}"),
    ("pancancer_honest_ci_lo", "pan-cancer partition-widened interval, low",
     r"(0\.3152\] to \[){v}(, [\d.]+\]: the)", "{:.4f}"),
    ("pancancer_honest_ci_hi", "pan-cancer partition-widened interval, high",
     r"(0\.3152\] to \[[\d.]+, ){v}(\]: the)", "{:.4f}"),
    ("pv24_pancancer_understated_pct", "pan-cancer 24-partition understatement",
     r"(16\.0% too narrow \(){v}(% at 24\))", "{:.1f}"),
    ("pv24_pancancer_share_pct", "pan-cancer 24-partition share",
     r"(25\.7% of total variance \(){v}(% at 24\))", "{:.1f}"),
    ("pv24_bootstrap_ratio", "bootstrap SE ratio between cohorts",
     r"(sampling theory says it should, by ){v}( against the)", "{:.3f}"),
    ("pv24_sqrt_patients", "square-root patient-count prediction",
     r"(should, by [\d.]+ against the ){v}( predicted by the)", "{:.3f}"),
    ("sens_ns_k3", "NSCLC index at 3 folds",
     r"(3 folds the NSCLC index is ){v}( and with 10)", "{:.4f}"),
    ("sens_ns_k10", "NSCLC index at 10 folds",
     r"(and with 10 it is ){v}(, so its value)", "{:.4f}"),
    ("sens_ns_unmatched_null", "NSCLC index, size-only null",
     r"(size and mean expression, gives ){v}(, so expression)", "{:.4f}"),
    ("e16_pan_registered_f050_median", "pan-cancer site AUROC at half the cohort",
     r"(the median stays at ){v}( at half the cohort)", "{:.3f}"),
    ("e16_pan_registered_f025_median", "pan-cancer site AUROC at a quarter",
     r"(half the cohort and ){v}( at a quarter)", "{:.3f}"),
]

CLAIMS += [(key, _a16_claim(shape).replace(r"([\d,]+)", r"(\d[\d,]*(?:\.\d+)?)", 1), 0.0,
            f"B-21 summary figure, {what}", None)
           for key, what, shape, _fmt in B21_SUMMARY_SPEC]

# Session 60 (H89): the science audit's unchecked literals, swept. Every one of
# its 1,072 unchecked, non-structural literals was moved one step in its last
# printed digit (in-process, the audit rescanned alone, behind a caught
# positive control): 7 CAUGHT, 1,065 not. These are the ones that state a
# CURRENT result an existing authority already holds -- the sort audit's
# counts, the sensitivity arms, Control C, E8, E11, E16, the A5 partition
# tables, the small-panel bracket -- each bound at printed precision to that
# authority, with a generated self-test row. Shapes are the literal's own
# line context, generated and checked to match exactly one line in every
# scanned document. What is NOT here, and why, is classified in the audit's
# "Session 60: the 1,072 unchecked numbers" block.
AUDIT_SWEEP_SPEC: list[tuple[str, str, str, str]] = [
    ('sortaudit_macos_samples', 'audit sweep :690 sortaudit_macos_samples',
     '(py`\\ on\\ the\\ real\\ NSCLC\\ matrix\\ \\(){v}(\\ ×\\ 41,046\\):)', '{:.0f}'),
    ('sortaudit_macos_genes', 'audit sweep :690 sortaudit_macos_genes',
     '(\\ the\\ real\\ NSCLC\\ matrix\\ \\(944\\ ×\\ ){v}(\\):)', '{:,.0f}'),
    ('sortaudit_macos_tied_rows', 'audit sweep :694 sortaudit_macos_tied_rows',
     '(\\ least\\ one\\ exact\\ rank\\ tie\\ \\|\\ \\*\\*){v}(\\ of\\ 944\\*\\*\\ \\|)', '{:.0f}'),
    ('sortaudit_macos_samples', 'audit sweep :694 sortaudit_macos_samples',
     '(one\\ exact\\ rank\\ tie\\ \\|\\ \\*\\*944\\ of\\ ){v}(\\*\\*\\ \\|)', '{:.0f}'),
    ('sortaudit_macos_tied_elements', 'audit sweep :695 sortaudit_macos_tied_elements',
     '(\\|\\ tied\\ rank\\ entries\\ \\|\\ \\*\\*){v}(\\ of\\ 38,747,424\\ \\(91)', '{:,.0f}'),
    ('sortaudit_macos_cells_diff', 'audit sweep :698 sortaudit_macos_cells_diff',
     '(differing\\ between\\ the\\ two\\ \\|\\ \\*\\*){v}(\\ of\\ 15,104\\*\\*\\ \\|)', '{:,.0f}'),
    ('sortaudit_macos_cells_total', 'audit sweep :698 sortaudit_macos_cells_total',
     '(between\\ the\\ two\\ \\|\\ \\*\\*15,060\\ of\\ ){v}(\\*\\*\\ \\|)', '{:,.0f}'),
    ('sortaudit_macos_frac_sd', 'audit sweep :702 sortaudit_macos_frac_sd',
     '(\\ as\\ a\\ fraction\\ of\\ that\\ sd\\ \\|\\ \\*\\*){v}(\\*\\*\\ \\|)', '{:.3f}'),
    ('sortaudit_macos_frac_sd', 'audit sweep :708 sortaudit_macos_frac_sd',
     '(happens\\ constantly\\.\\ \\*\\*){v}(\\ of\\ a\\ standard\\ dev)', '{:.3f}'),
    ('a10ss_shift_pct_of_halfwidth', 'audit sweep :733 a10ss_shift_pct_of_halfwidth',
     '(ves\\ \\*\\*0\\.00617520558066495\\*\\*\\ —\\ ){v}(%\\ of\\ the\\ pre\\-fix)', '{:.1f}'),
    ('sortaudit_macos_frac_sd', 'audit sweep :735 sortaudit_macos_frac_sd',
     '(ct\\ moved\\ the\\ NSCLC\\ primary\\.\\ A\\ ){v}(\\-of\\-a\\-score\\-sd\\ mov)', '{:.3f}'),
    ('sortaudit_macos_cells_diff', 'audit sweep :735 sortaudit_macos_cells_diff',
     '(606\\-of\\-a\\-score\\-sd\\ movement\\ in\\ ){v}(\\ of)', '{:,.0f}'),
    ('sortaudit_macos_cells_total', 'audit sweep :736 sortaudit_macos_cells_total',
     '(){v}(\\ per\\-cell\\ scores\\ t)', '{:,.0f}'),
    ('sortaudit_macos_tied_rows', 'audit sweep :782 sortaudit_macos_tied_rows',
     '(\\|\\ samples\\ with\\ ties\\ \\|\\ ){v}(\\ of\\ 944\\ \\|\\ 944\\ of\\ 9)', '{:.0f}'),
    ('sortaudit_macos_samples', 'audit sweep :782 sortaudit_macos_samples',
     '(\\|\\ samples\\ with\\ ties\\ \\|\\ 944\\ of\\ ){v}(\\ \\|\\ 944\\ of\\ 944\\ \\|\\ id)', '{:.0f}'),
    ('sortaudit_hpc4_tied_rows', 'audit sweep :782 sortaudit_hpc4_tied_rows',
     '(ples\\ with\\ ties\\ \\|\\ 944\\ of\\ 944\\ \\|\\ ){v}(\\ of\\ 944\\ \\|\\ identica)', '{:.0f}'),
    ('sortaudit_hpc4_samples', 'audit sweep :782 sortaudit_hpc4_samples',
     '(th\\ ties\\ \\|\\ 944\\ of\\ 944\\ \\|\\ 944\\ of\\ ){v}(\\ \\|\\ identical\\ \\|)', '{:.0f}'),
    ('sortaudit_macos_tied_elements', 'audit sweep :783 sortaudit_macos_tied_elements',
     '(\\|\\ tied\\ rank\\ entries\\ \\|\\ ){v}(\\ \\|\\ 35,394,448\\ \\|\\ id)', '{:,.0f}'),
    ('sortaudit_hpc4_tied_elements', 'audit sweep :783 sortaudit_hpc4_tied_elements',
     '(d\\ rank\\ entries\\ \\|\\ 35,394,448\\ \\|\\ ){v}(\\ \\|\\ identical\\ \\|)', '{:,.0f}'),
    ('sortaudit_macos_frac_sd', 'audit sweep :787 sortaudit_macos_frac_sd',
     '(\\|\\ as\\ a\\ fraction\\ of\\ score\\ sd\\ \\|\\ ){v}(\\ \\|\\ 0\\.618\\ \\|\\ differs)', '{:.3f}'),
    ('sortaudit_hpc4_frac_sd', 'audit sweep :787 sortaudit_hpc4_frac_sd',
     '(raction\\ of\\ score\\ sd\\ \\|\\ 0\\.606\\ \\|\\ ){v}(\\ \\|\\ differs\\ \\(quicks)', '{:.3f}'),
    ('sortaudit_macos_tied_elements', 'audit sweep :789 sortaudit_macos_tied_elements',
     '(\\ both\\ machines\\ count\\ the\\ same\\ ){v}(\\ tied)', '{:,.0f}'),
    ('a16_ns_patients', 'audit sweep :1196 a16_ns_patients',
     '(y`\\ on\\ `\\-\\-tumour\\-only`\\ inputs;\\ ){v}(\\ patients,\\ the\\ nor)', '{:.0f}'),
    ('a16_pan_patients', 'audit sweep :1211 a16_pan_patients',
     '(\\(){v}(\\ patients\\)\\ gives\\ a)', '{:,.0f}'),
    ('a16_pan_patients', 'audit sweep :1217 a16_pan_patients',
     '(pancancer_tumouronly/`\\)\\.\\ Over\\ ){v}(\\ patients\\ and\\ 619\\ )', '{:,.0f}'),
    ('pancancer_n_sites', 'audit sweep :1217 pancancer_n_sites',
     '(y/`\\)\\.\\ Over\\ 7,128\\ patients\\ and\\ ){v}(\\ sites,)', '{:.0f}'),
    ('a16_pan_fall', 'audit sweep :1219 a16_pan_fall',
     '(\\.2749,\\ 0\\.3189\\]\\.\\ The\\ shift\\ is\\ −){v}(,\\ less\\ than\\ half\\ t)', '{:.4f}'),
    ('pv24_pancancer_sd', 'audit sweep :1220 pv24_pancancer_sd',
     '(standard\\ deviation\\ \\(){v}(\\)\\.\\ In\\ NSCLC\\ the\\ co)', '{:.4f}'),
    ('sens_ns_global_curated_r', 'audit sweep :1260 sens_ns_global_curated_r',
     '(96\\ \\[−0\\.0630,\\ 0\\.0466\\]\\ \\|\\ 1/16\\ \\|\\ ){v}(\\ \\|\\ 0\\.361\\ \\|\\ 0\\.944\\ \\|)', '{:.3f}'),
    ('sens_ns_global_null_r', 'audit sweep :1260 sens_ns_global_null_r',
     '(630,\\ 0\\.0466\\]\\ \\|\\ 1/16\\ \\|\\ 0\\.351\\ \\|\\ ){v}(\\ \\|\\ 0\\.944\\ \\|\\ −0\\.004\\ )', '{:.3f}'),
    ('sens_pan_ref_curated_r', 'audit sweep :1284 sens_pan_ref_curated_r',
     '(68\\ \\[0\\.2749,\\ 0\\.3189\\]\\ \\|\\ 16/16\\ \\|\\ ){v}(\\ \\|\\ 0\\.100\\ \\|\\ 0\\.025\\ \\|)', '{:.3f}'),
    ('sens_pan_ref_null_r', 'audit sweep :1284 sens_pan_ref_null_r',
     '(49,\\ 0\\.3189\\]\\ \\|\\ 16/16\\ \\|\\ 0\\.391\\ \\|\\ ){v}(\\ \\|\\ 0\\.025\\ \\|\\ 0\\.801\\ \\|)', '{:.3f}'),
    ('sens_pan_global_curated_r', 'audit sweep :1285 sens_pan_global_curated_r',
     '(12\\ \\[0\\.2692,\\ 0\\.4353\\]\\ \\|\\ 15/16\\ \\|\\ ){v}(\\ \\|\\ 0\\.452\\ \\|\\ 0\\.069\\ \\|)', '{:.3f}'),
    ('sens_pan_global_null_r', 'audit sweep :1285 sens_pan_global_null_r',
     '(92,\\ 0\\.4353\\]\\ \\|\\ 15/16\\ \\|\\ 0\\.725\\ \\|\\ ){v}(\\ \\|\\ 0\\.069\\ \\|\\ 0\\.805\\ \\|)', '{:.3f}'),
    ('sens_pan_ref_null_r', 'audit sweep :1289 sens_pan_ref_null_r',
     '(et\\.\\ Random\\ sets\\ rise\\ from\\ r\\ =\\ ){v}()', '{:.3f}'),
    ('sens_pan_global_null_r', 'audit sweep :1290 sens_pan_global_null_r',
     '(to\\ ){v}(,\\ but\\ the\\ curated\\ )', '{:.3f}'),
    ('sens_pan_ref_curated_r', 'audit sweep :1290 sens_pan_ref_curated_r',
     '(signatures\\ rise\\ further,\\ from\\ ){v}(\\ to\\ 0\\.725\\.\\ The)', '{:.3f}'),
    ('sens_pan_global_curated_r', 'audit sweep :1290 sens_pan_global_curated_r',
     '(s\\ rise\\ further,\\ from\\ 0\\.391\\ to\\ ){v}(\\.\\ The)', '{:.3f}'),
    ('sens_ns_unmatched_shift', 'audit sweep :1366 sens_ns_unmatched_shift',
     '(\\ matching\\ lowers\\ the\\ index\\ by\\ ){v}(\\ \\(what\\ matching\\ bu)', '{:.3f}'),
    ('sens_ns_null_boot_400_lo', 'audit sweep :1367 sens_ns_null_boot_400_lo',
     '(he\\ bootstrap\\ \\|\\ unchanged;\\ CI\\ \\[){v}(,\\ 0\\.3512\\]\\ \\|\\ 16/16\\ )', '{:.4f}'),
    ('sens_ns_null_boot_400_hi', 'audit sweep :1367 sens_ns_null_boot_400_hi',
     '(trap\\ \\|\\ unchanged;\\ CI\\ \\[0\\.2336,\\ ){v}(\\]\\ \\|\\ 16/16\\ \\|\\ interv)', '{:.4f}'),
    ('sens_ns_null_boot_1000_lo', 'audit sweep :1368 sens_ns_null_boot_1000_lo',
     '(\\ 1,000\\ draws\\ \\|\\ unchanged;\\ CI\\ \\[){v}(,\\ 0\\.3511\\]\\ \\|\\ 16/16\\ )', '{:.4f}'),
    ('sens_ns_null_boot_1000_hi', 'audit sweep :1368 sens_ns_null_boot_1000_hi',
     '(raws\\ \\|\\ unchanged;\\ CI\\ \\[0\\.2333,\\ ){v}(\\]\\ \\|\\ 16/16\\ \\|\\ at\\ mos)', '{:.4f}'),
    ('sens_ns_k3', 'audit sweep :1372 sens_ns_k3',
     '(The\\ index\\ rises\\ with\\ K:\\ ){v}(,\\ 0\\.292\\ and\\ 0\\.315\\ )', '{:.3f}'),
    ('sens_ns_reference', 'audit sweep :1372 sens_ns_reference',
     '(he\\ index\\ rises\\ with\\ K:\\ 0\\.219,\\ ){v}(\\ and\\ 0\\.315\\ at\\ 3,\\ 5)', '{:.3f}'),
    ('sens_ns_k10', 'audit sweep :1372 sens_ns_k10',
     '(ises\\ with\\ K:\\ 0\\.219,\\ 0\\.292\\ and\\ ){v}(\\ at\\ 3,\\ 5\\ and\\ 10\\ fo)', '{:.3f}'),
    ('e11_p_all_inside', 'audit sweep :1406 e11_p_all_inside',
     '(\\(all\\ five\\ inside\\)\\ =\\ 6/4,368\\ =\\ ){v}(\\.)', '{:.4f}'),
    ('e11_relabellings', 'audit sweep :1408 e11_relabellings',
     '(\\ \\ \\ No\\ other\\ of\\ the\\ C\\(16,\\ 6\\)\\ =\\ ){v}(\\ relabellings\\ does)', '{:,.0f}'),
    ('e11_relabellings', 'audit sweep :1409 e11_relabellings',
     '(\\ \\ \\ \\ P\\ =\\ 1/){v}(\\)\\.)', '{:,.0f}'),
    ('immune_mwu_p_2s', 'audit sweep :1413 immune_mwu_p_2s',
     '(immune\\ vs\\ other:\\ Mann–Whitney\\ ){v}(,\\ Fisher)', '{:.3f}'),
    ('immune_fisher_p_2s', 'audit sweep :1414 immune_fisher_p_2s',
     '(\\ \\ \\ \\ ){v}(\\)\\ is\\ unchanged\\ and)', '{:.3f}'),
    ('sortaudit_macos_max_diff', 'audit sweep :1424 sortaudit_macos_max_diff',
     '(\\ \\ Limitations\\ 8\\ now\\ says\\ its\\ ){v}(\\ /\\ 0\\.0211\\ /\\ 0\\.61\\ a)', '{:.4f}'),
    ('sortaudit_macos_frac_sd', 'audit sweep :1424 sortaudit_macos_frac_sd',
     '(ow\\ says\\ its\\ 0\\.0128\\ /\\ 0\\.0211\\ /\\ ){v}(\\ are\\ macOS\\ values)', '{:.2f}'),
    ('sortaudit_hpc4_max_diff', 'audit sweep :1425 sortaudit_hpc4_max_diff',
     "(nd\\ gives\\ Linux's\\ beside\\ them\\ \\(){v}(\\ /\\ 0\\.0211\\ /\\ 0\\.62\\)\\.)", '{:.4f}'),
    ('sortaudit_hpc4_score_sd', 'audit sweep :1425 sortaudit_hpc4_score_sd',
     "(Linux's\\ beside\\ them\\ \\(0\\.0130\\ /\\ ){v}(\\ /\\ 0\\.62\\)\\.\\ All\\ thre)", '{:.4f}'),
    ('sortaudit_hpc4_frac_sd', 'audit sweep :1425 sortaudit_hpc4_frac_sd',
     '(eside\\ them\\ \\(0\\.0130\\ /\\ 0\\.0211\\ /\\ ){v}(\\)\\.\\ All\\ three\\ Linux)', '{:.2f}'),
    ('nsclc_label_site_variance_pct', 'audit sweep :1460 nsclc_label_site_variance_pct',
     '(14,\\ the\\ NSCLC/pan\\-cancer\\ gap\\ \\(){v}(%\\ against\\ 44%\\)\\.\\*\\*\\ )', '{:.1f}'),
    ('label_site_variance_r2', 'audit sweep :1466 label_site_variance_r2',
     '(\\|\\ site\\ R²,\\ crude\\ \\|\\ ){v}(\\ \\|\\ 0\\.145\\ \\|)', '{:.3f}'),
    ('nsclc_label_site_variance_r2', 'audit sweep :1466 nsclc_label_site_variance_r2',
     '(\\|\\ site\\ R²,\\ crude\\ \\|\\ 0\\.444\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('label_site_given_type_r2', 'audit sweep :1469 label_site_given_type_r2',
     '(\\|\\ site\\ given\\ cancer\\ type\\ \\|\\ ){v}(\\ \\|\\ 0\\.099\\ \\|)', '{:.3f}'),
    ('e14_pancancer_median_r2', 'audit sweep :1470 e14_pancancer_median_r2',
     '(ll\\ sites,\\ 200\\ permutations\\)\\ \\|\\ ){v}(\\ \\(6/16\\ at\\ p\\ ≤\\ 0\\.05)', '{:.3f}'),
    ('e14_nsclc_median_r2', 'audit sweep :1470 e14_nsclc_median_r2',
     '(\\ \\|\\ 0\\.047\\ \\(6/16\\ at\\ p\\ ≤\\ 0\\.05\\)\\ \\|\\ ){v}(\\ \\(15/16\\)\\ \\|)', '{:.3f}'),
    ('cap_pan_lsv_sites', 'audit sweep :1490 cap_pan_lsv_sites',
     '(\\|\\ pan\\-cancer\\ \\|\\ ){v}(\\ of\\ 619\\ \\|\\ 417\\ \\(lar)', '{:.0f}'),
    ('pancancer_n_sites', 'audit sweep :1490 pancancer_n_sites',
     '(\\|\\ pan\\-cancer\\ \\|\\ 202\\ of\\ ){v}(\\ \\|\\ 417\\ \\(largest\\ dr)', '{:.0f}'),
    ('cap_pan_lsv_n', 'audit sweep :1490 cap_pan_lsv_n',
     '(he\\ complete\\-case\\ mask,\\ giving\\ ){v}(\\ \\|)', '{:,.0f}'),
    ('cap_nsclc_lsv_n', 'audit sweep :1491 cap_nsclc_lsv_n',
     '(NSCLC\\ \\|\\ 31\\ of\\ 68\\ \\|\\ 37\\ \\|\\ 148\\ \\|\\ ){v}(\\ \\|)', '{:.0f}'),
    ('e6_plage_recon', 'audit sweep :1575 e6_plage_recon',
     '(\\ split\\-half\\ reconstruction\\ \\|\\ \\+){v}(\\ \\|\\ 11/16\\ \\|)', '{:.4f}'),
    ('a7_recon_ssgsea', 'audit sweep :1576 a7_recon_ssgsea',
     '(\\ split\\-half\\ reconstruction\\ \\|\\ −){v}(\\ \\|\\ 0/16\\ \\|)', '{:.4f}'),
    ('e8_median_alpha_observed', 'audit sweep :1597 e8_median_alpha_observed',
     '(\\|\\ curated\\ sets\\ \\|\\ ){v}(\\ \\|\\ 0\\.9773\\ \\|)', '{:.4f}'),
    ('e8_median_omega_observed', 'audit sweep :1597 e8_median_omega_observed',
     '(\\|\\ curated\\ sets\\ \\|\\ 0\\.9773\\ \\|\\ ){v}(\\ \\|)', '{:.4f}'),
    ('e8_median_alpha_null', 'audit sweep :1598 e8_median_alpha_null',
     '(random\\ sets\\ \\(25\\ draws\\ each\\)\\ \\|\\ ){v}(\\ \\|\\ 0\\.2270\\ \\|)', '{:.4f}'),
    ('e8_median_omega_null', 'audit sweep :1598 e8_median_omega_null',
     '(ts\\ \\(25\\ draws\\ each\\)\\ \\|\\ 0\\.8358\\ \\|\\ ){v}(\\ \\|)', '{:.4f}'),
    ('e8_median_gap_alpha', 'audit sweep :1599 e8_median_gap_alpha',
     '(\\|\\ curated\\ −\\ random\\ gap\\ \\|\\ ){v}(\\ \\|\\ 0\\.7333\\ \\|)', '{:.4f}'),
    ('e8_median_gap_omega', 'audit sweep :1599 e8_median_gap_omega',
     '(rated\\ −\\ random\\ gap\\ \\|\\ 0\\.1264\\ \\|\\ ){v}(\\ \\|)', '{:.4f}'),
    ('e16_pan_registered_f100_median', 'audit sweep :1626 e16_pan_registered_f100_median',
     '(cancer,\\ sites\\ ≥\\ 20\\ patients\\ \\|\\ ){v}(\\ \\(98\\ sites\\)\\ \\|\\ 0\\.99)', '{:.4f}'),
    ('e16_pan_registered_f050_median', 'audit sweep :1626 e16_pan_registered_f050_median',
     '(atients\\ \\|\\ 0\\.9985\\ \\(98\\ sites\\)\\ \\|\\ ){v}(\\ \\(40–42\\)\\ \\|\\ 0\\.9940\\ )', '{:.4f}'),
    ('e16_pan_registered_f025_median', 'audit sweep :1626 e16_pan_registered_f025_median',
     '(\\(98\\ sites\\)\\ \\|\\ 0\\.9976\\ \\(40–42\\)\\ \\|\\ ){v}(\\ \\(10–14\\)\\ \\|)', '{:.4f}'),
    ('e16_ns_registered_f100_median', 'audit sweep :1627 e16_ns_registered_f100_median',
     '(\\ NSCLC,\\ sites\\ ≥\\ 20\\ patients\\ \\|\\ ){v}(\\ \\(15\\ sites\\)\\ \\|\\ 0\\.98)', '{:.4f}'),
    ('e16_ns_registered_f050_median', 'audit sweep :1627 e16_ns_registered_f050_median',
     '(atients\\ \\|\\ 0\\.9894\\ \\(15\\ sites\\)\\ \\|\\ ){v}(\\ \\(4–5\\)\\ \\|\\ 0\\.894–0\\.9)', '{:.4f}'),
    ('e16_ns_registered_f025_lowest', 'audit sweep :1627 e16_ns_registered_f025_lowest',
     '(4\\ \\(15\\ sites\\)\\ \\|\\ 0\\.9805\\ \\(4–5\\)\\ \\|\\ ){v}(–0\\.972\\ \\(1–2\\)\\ \\|)', '{:.3f}'),
    ('e16_ns_registered_f025_highest', 'audit sweep :1627 e16_ns_registered_f025_highest',
     '(sites\\)\\ \\|\\ 0\\.9805\\ \\(4–5\\)\\ \\|\\ 0\\.894–){v}(\\ \\(1–2\\)\\ \\|)', '{:.3f}'),
    ('e16_pan_small_f100_median', 'audit sweep :1628 e16_pan_small_f100_median',
     '(\\-cancer,\\ sites\\ ≥\\ 5\\ patients\\ \\|\\ ){v}(\\ \\(320\\ sites,\\ 97%\\ a)', '{:.4f}'),
    ('e16_pan_small_f100_sites', 'audit sweep :1628 e16_pan_small_f100_sites',
     '(\\ sites\\ ≥\\ 5\\ patients\\ \\|\\ 0\\.9981\\ \\(){v}(\\ sites,\\ 97%\\ above\\ )', '{:.0f}'),
    ('e16_pan_small_f025_lowest', 'audit sweep :1628 e16_pan_small_f025_lowest',
     '(0\\ sites,\\ 97%\\ above\\ 0\\.9\\)\\ \\|\\ —\\ \\|\\ ){v}(–0\\.9969\\ \\(114–123\\)\\ )', '{:.4f}'),
    ('e16_pan_small_f025_highest', 'audit sweep :1628 e16_pan_small_f025_highest',
     '(,\\ 97%\\ above\\ 0\\.9\\)\\ \\|\\ —\\ \\|\\ 0\\.9963–){v}(\\ \\(114–123\\)\\ \\|)', '{:.4f}'),
    ('e16_pan_small_quarter_sites_min', 'audit sweep :1628 e16_pan_small_quarter_sites_min',
     '(ove\\ 0\\.9\\)\\ \\|\\ —\\ \\|\\ 0\\.9963–0\\.9969\\ \\(){v}(–123\\)\\ \\|)', '{:.0f}'),
    ('e16_pan_small_quarter_sites_max', 'audit sweep :1628 e16_pan_small_quarter_sites_max',
     '(0\\.9\\)\\ \\|\\ —\\ \\|\\ 0\\.9963–0\\.9969\\ \\(114–){v}(\\)\\ \\|)', '{:.0f}'),
    ('e16_ns_small_f100_median', 'audit sweep :1629 e16_ns_small_f100_median',
     '(\\|\\ NSCLC,\\ sites\\ ≥\\ 5\\ patients\\ \\|\\ ){v}(\\ \\(47\\ sites,\\ 94%\\ ab)', '{:.4f}'),
    ('e16_ns_small_f025_lowest', 'audit sweep :1629 e16_ns_small_f025_lowest',
     '(7\\ sites,\\ 94%\\ above\\ 0\\.9\\)\\ \\|\\ —\\ \\|\\ ){v}(–0\\.9783\\ \\(17–20\\)\\ \\|)', '{:.4f}'),
    ('e16_ns_small_f025_highest', 'audit sweep :1629 e16_ns_small_f025_highest',
     '(,\\ 94%\\ above\\ 0\\.9\\)\\ \\|\\ —\\ \\|\\ 0\\.9623–){v}(\\ \\(17–20\\)\\ \\|)', '{:.4f}'),
    ('e16_ns_small_f100_median', 'audit sweep :1634 e16_ns_small_f100_median',
     '(threshold\\ the\\ median\\ is\\ still\\ ){v}(\\ and\\ 0\\.998,\\ so\\ sit)', '{:.3f}'),
    ('e16_pan_small_f100_median', 'audit sweep :1634 e16_pan_small_f100_median',
     '(the\\ median\\ is\\ still\\ 0\\.987\\ and\\ ){v}(,\\ so\\ sites)', '{:.3f}'),
    ('e16_ns_registered_f100_median', 'audit sweep :1644 e16_ns_registered_f100_median',
     '(l\\-cohort\\ NSCLC\\ median\\ here\\ is\\ ){v}(\\ on\\ Linux\\ against\\ )', '{:.4f}'),
    ('nsclc_site_auroc', 'audit sweep :1645 nsclc_site_auroc',
     '(\\ \\ ){v}(\\ on\\ macOS\\.\\ That\\ is)', '{:.4f}'),
    ('nsclc_isi', 'audit sweep :2174 nsclc_isi',
     '(stored\\ \\+){v}(\\ —\\ \\*\\*failed\\*\\*,\\ ret)', '{:.4f}'),
    ('bisect_hpc4_pooled', 'audit sweep :2243 bisect_hpc4_pooled',
     '(eir\\ mean\\-z\\ arms\\ return\\ 0\\.2964/){v}(\\ —\\ that\\ is\\ A9,\\ not)', '{:.4f}'),
    ('pancancer_hypoxia_excess_lo', 'audit sweep :2267 pancancer_hypoxia_excess_lo',
     '(eat\\ its\\ null\\ at\\ `excess_lo\\ =\\ \\+){v}(`\\.\\ It)', '{:.4f}'),
    ('label_site_variance_r2', 'audit sweep :2286 label_site_variance_r2',
     '(\\ Site\\ on\\ the\\ label,\\ crude\\ \\|\\ \\*\\*){v}(\\*\\*\\ \\|)', '{:.4f}'),
    ('label_site_given_type_r2', 'audit sweep :2287 label_site_given_type_r2',
     '(ite\\ \\*\\*given\\ cancer\\ type\\*\\*\\ \\|\\ \\*\\*){v}(\\*\\*\\ \\|)', '{:.4f}'),
    ('label_site_given_type_pct', 'audit sweep :2293 label_site_given_type_pct',
     '(on\\ cancer\\ type\\ it\\ falls\\ to\\ ){v}(%,\\ and\\ \\*\\*plate\\ —\\ b)', '{:.1f}'),
    ('label_site_given_type_pct', 'audit sweep :2319 label_site_given_type_pct',
     '(ite\\ given\\ cancer\\ type"\\ route\\ \\(){v}(%\\)\\.\\ Two\\ different\\ )', '{:.1f}'),
    ('label_site_given_type_pct', 'audit sweep :2330 label_site_given_type_pct',
     '("agreement\\ with\\ ){v}(%"\\ argument\\ still\\ )', '{:.1f}'),
    ('control_c_calibrated_pct', 'audit sweep :2324 control_c_calibrated_pct',
     '(16\\.\\ The\\ manuscript\\ quoted\\ its\\ ){v}(%\\ beside\\ two\\ 16\\-si)', '{:.1f}'),
    ('e14_pancancer_median_pct', 'audit sweep :2326 e14_pancancer_median_pct',
     '(ation/`\\):\\ pan\\-cancer\\ median\\ \\*\\*){v}(%\\*\\*\\ \\(range)', '{:.1f}'),
    ('e14_pancancer_min_pct', 'audit sweep :2327 e14_pancancer_min_pct',
     '(){v}(–6\\.0%,\\ 6\\ of\\ 16\\ at\\ )', '{:.1f}'),
    ('e14_pancancer_max_pct', 'audit sweep :2327 e14_pancancer_max_pct',
     '(2\\.8–){v}(%,\\ 6\\ of\\ 16\\ at\\ p\\ ≤\\ )', '{:.1f}'),
    ('e14_nsclc_median_pct', 'audit sweep :2327 e14_nsclc_median_pct',
     '(\\ at\\ p\\ ≤\\ 0\\.05\\),\\ NSCLC\\ median\\ \\*\\*){v}(%\\*\\*\\ \\(range\\ 1\\.8–10\\.)', '{:.1f}'),
    ('e14_nsclc_min_pct', 'audit sweep :2327 e14_nsclc_min_pct',
     '(\\ NSCLC\\ median\\ \\*\\*7\\.0%\\*\\*\\ \\(range\\ ){v}(–10\\.5%,\\ 15\\ of)', '{:.1f}'),
    ('e14_nsclc_max_pct', 'audit sweep :2327 e14_nsclc_max_pct',
     '(LC\\ median\\ \\*\\*7\\.0%\\*\\*\\ \\(range\\ 1\\.8–){v}(%,\\ 15\\ of)', '{:.1f}'),
    ('small_panel_k10', 'audit sweep :2339 small_panel_k10',
     '(\\|\\ \\*\\*10\\*\\*\\ \\|\\ \\*\\*){v}(\\*\\*\\ \\|\\ 16\\ \\|)', '{:.4f}'),
    ('small_panel_k20', 'audit sweep :2340 small_panel_k20',
     '(\\|\\ 20\\ \\|\\ ){v}(\\ \\|\\ 16\\ \\|)', '{:.4f}'),
    ('small_panel_k40', 'audit sweep :2341 small_panel_k40',
     '(\\|\\ 40\\ \\|\\ ){v}(\\ \\|\\ 15\\ \\|)', '{:.4f}'),
    ('small_panel_k80', 'audit sweep :2342 small_panel_k80',
     '(\\|\\ 80\\ \\|\\ ){v}(\\ \\|\\ 14\\ \\|)', '{:.4f}'),
    ('small_panel_k160', 'audit sweep :2343 small_panel_k160',
     '(\\|\\ 160\\ \\|\\ ){v}(\\ \\|\\ 11\\ \\|)', '{:.4f}'),
    ('small_panel_ratio_k10_k160', 'audit sweep :2346 small_panel_ratio_k10_k160',
     '(range\\.\\ At\\ k=10\\ the\\ gap\\ is\\ \\*\\*){v}(×\\ larger\\*\\*\\ than\\ at)', '{:.1f}'),
    ('rotation_p', 'audit sweep :2367 rotation_p',
     '(\\|\\ \\*\\*p\\*\\*\\ \\|\\ \\*\\*){v}(\\*\\*\\ \\(B\\ =\\ 1000\\)\\ \\|)', '{:.6f}'),
    ('nsclc_reported_ci_lo', 'audit sweep :2405 nsclc_reported_ci_lo',
     '(\\|\\ Reported\\ interval\\ \\|\\ \\[){v}(,\\ 0\\.3695\\],\\ width\\ 0)', '{:.4f}'),
    ('nsclc_reported_ci_hi', 'audit sweep :2405 nsclc_reported_ci_hi',
     '(\\ Reported\\ interval\\ \\|\\ \\[0\\.2542,\\ ){v}(\\],\\ width\\ 0\\.1153\\ \\|)', '{:.4f}'),
    ('nsclc_honest_ci_lo', 'audit sweep :2406 nsclc_honest_ci_lo',
     '(\\|\\ \\*\\*Honest\\ interval\\*\\*\\ \\|\\ \\*\\*\\[){v}(,\\ 0\\.3720\\],\\ width\\ 0)', '{:.4f}'),
    ('nsclc_honest_ci_hi', 'audit sweep :2406 nsclc_honest_ci_hi',
     '(onest\\ interval\\*\\*\\ \\|\\ \\*\\*\\[0\\.2516,\\ ){v}(\\],\\ width\\ 0\\.1204\\*\\*\\ )', '{:.4f}'),
    ('nsclc_interval_understated_pct', 'audit sweep :2408 nsclc_interval_understated_pct',
     '(\\*\\*The\\ reported\\ interval\\ is\\ ){v}(%\\ too\\ narrow\\ —\\ not)', '{:.1f}'),
    ('pancancer_isi', 'audit sweep :2419 pancancer_isi',
     '(\\|\\ 0\\ \\|\\ ){v}(\\ \\|\\ 6,266\\ s\\ \\|)', '{:.4f}'),
    ('pancancer_reported_ci_lo', 'audit sweep :2436 pancancer_reported_ci_lo',
     '(\\|\\ Reported\\ interval\\ \\|\\ \\[){v}(,\\ 0\\.3152\\],\\ width\\ 0)', '{:.4f}'),
    ('pancancer_reported_ci_hi', 'audit sweep :2436 pancancer_reported_ci_hi',
     '(\\ Reported\\ interval\\ \\|\\ \\[0\\.2717,\\ ){v}(\\],\\ width\\ 0\\.0435\\ \\|)', '{:.4f}'),
    ('pancancer_honest_ci_lo', 'audit sweep :2437 pancancer_honest_ci_lo',
     '(\\|\\ \\*\\*Honest\\ interval\\*\\*\\ \\|\\ \\*\\*\\[){v}(,\\ 0\\.3187\\],\\ width\\ 0)', '{:.4f}'),
    ('pancancer_honest_ci_hi', 'audit sweep :2437 pancancer_honest_ci_hi',
     '(onest\\ interval\\*\\*\\ \\|\\ \\*\\*\\[0\\.2683,\\ ){v}(\\],\\ width\\ 0\\.0505\\*\\*\\ )', '{:.4f}'),
    ('pancancer_interval_understated_pct', 'audit sweep :2439 pancancer_interval_understated_pct',
     '(\\*\\*The\\ reported\\ interval\\ is\\ ){v}(%\\ too\\ narrow\\.\\*\\*\\ Pa)', '{:.1f}'),
    ('pancancer_partition_var_pct', 'audit sweep :2461 pancancer_partition_var_pct',
     '(ty\\ in\\ the\\ \\*\\*larger\\*\\*\\ cohort\\ —\\ ){v}(%\\ pan\\-cancer\\ again)', '{:.1f}'),
    ('nsclc_partition_var_pct', 'audit sweep :2461 nsclc_partition_var_pct',
     '(rt\\ —\\ 25\\.7%\\ pan\\-cancer\\ against\\ ){v}(%\\ in)', '{:.1f}'),
    ('cohort_n_ratio', 'audit sweep :2462 cohort_n_ratio',
     '(NSCLC,\\ on\\ a\\ cohort\\ ){v}(×\\ the\\ size\\.\\ The\\ tw)', '{:.1f}'),
    ('nsclc_bootstrap_se', 'audit sweep :2463 nsclc_bootstrap_se',
     '(\\ standard\\ error\\ falls\\ with\\ n,\\ ){v}(\\ →\\ 0\\.0111,\\ while\\ t)', '{:.4f}'),
    ('pancancer_bootstrap_se', 'audit sweep :2463 pancancer_bootstrap_se',
     '(\\ error\\ falls\\ with\\ n,\\ 0\\.0294\\ →\\ ){v}(,\\ while\\ the)', '{:.4f}'),
    ('nsclc_partition_sd', 'audit sweep :2464 nsclc_partition_sd',
     '(en\\-partition\\ sd\\ barely\\ moves,\\ ){v}(\\ →\\ 0\\.0065\\.\\ Partiti)', '{:.4f}'),
    ('pancancer_partition_sd', 'audit sweep :2464 pancancer_partition_sd',
     '(ion\\ sd\\ barely\\ moves,\\ 0\\.0089\\ →\\ ){v}(\\.\\ Partition\\ varian)', '{:.4f}'),
    ('nsclc_interval_understated_pct', 'audit sweep :2468 nsclc_interval_understated_pct',
     '(is\\ audit\\ item\\ posed\\ —\\ "is\\ the\\ ){v}(%\\ widening)', '{:.2f}'),
    ('nsclc_interval_understated_pct', 'audit sweep :2478 nsclc_interval_understated_pct',
     '(r,\\ and\\ the\\ measured\\ answer\\ is\\ ){v}(%\\.)', '{:.1f}'),
    ('pancancer_partition_var_pct', 'audit sweep :2487 pancancer_partition_var_pct',
     '(bution\\*\\*:\\ partition\\ choice\\ is\\ ){v}(%\\ of\\ total\\ ISI\\ var)', '{:.2f}'),
    ('nsclc_partition_var_pct', 'audit sweep :2488 nsclc_partition_var_pct',
     '(pan\\-cancer\\ against\\ ){v}(%\\ in\\ NSCLC,\\ on\\ 7\\.6)', '{:.2f}'),
    ('nsclc_partition_sd', 'audit sweep :2499 nsclc_partition_sd',
     '(\\|\\ NSCLC\\ \\|\\ 6\\ \\|\\ ){v}(\\ \\|\\ \\[0\\.005550,\\ 0\\.02)', '{:.6f}'),
    ('nsclc_partition_sd_chi2_lo', 'audit sweep :2499 nsclc_partition_sd_chi2_lo',
     '(\\|\\ NSCLC\\ \\|\\ 6\\ \\|\\ 0\\.008892\\ \\|\\ \\[){v}(,\\ 0\\.021807\\]\\ \\|\\ 3\\.93)', '{:.6f}'),
    ('nsclc_partition_sd_chi2_hi', 'audit sweep :2499 nsclc_partition_sd_chi2_hi',
     '(C\\ \\|\\ 6\\ \\|\\ 0\\.008892\\ \\|\\ \\[0\\.005550,\\ ){v}(\\]\\ \\|\\ 3\\.93x\\ \\|)', '{:.6f}'),
    ('pancancer_partition_sd', 'audit sweep :2500 pancancer_partition_sd',
     '(\\|\\ pan\\-cancer\\ \\|\\ 5\\ \\|\\ ){v}(\\ \\|\\ \\[0\\.003908,\\ 0\\.01)', '{:.6f}'),
    ('pancancer_partition_sd_chi2_lo', 'audit sweep :2500 pancancer_partition_sd_chi2_lo',
     '(\\ pan\\-cancer\\ \\|\\ 5\\ \\|\\ 0\\.006523\\ \\|\\ \\[){v}(,\\ 0\\.018745\\]\\ \\|\\ 4\\.80)', '{:.6f}'),
    ('pancancer_partition_sd_chi2_hi', 'audit sweep :2500 pancancer_partition_sd_chi2_hi',
     '(r\\ \\|\\ 5\\ \\|\\ 0\\.006523\\ \\|\\ \\[0\\.003908,\\ ){v}(\\]\\ \\|\\ 4\\.80x\\ \\|)', '{:.6f}'),
    ('nsclc_partition_var_pct', 'audit sweep :2510 nsclc_partition_var_pct',
     '(\\|\\ NSCLC\\ \\|\\ ){v}(%\\ \\|\\ \\[3\\.44%,\\ 35\\.48%)', '{:.2f}'),
    ('nsclc_partition_share_lo_pct', 'audit sweep :2510 nsclc_partition_share_lo_pct',
     '(\\|\\ NSCLC\\ \\|\\ 8\\.38%\\ \\|\\ \\[){v}(%,\\ 35\\.48%\\]\\ \\|)', '{:.2f}'),
    ('nsclc_partition_share_hi_pct', 'audit sweep :2510 nsclc_partition_share_hi_pct',
     '(\\|\\ NSCLC\\ \\|\\ 8\\.38%\\ \\|\\ \\[3\\.44%,\\ ){v}(%\\]\\ \\|)', '{:.2f}'),
    ('pancancer_partition_var_pct', 'audit sweep :2511 pancancer_partition_var_pct',
     '(\\|\\ pan\\-cancer\\ \\|\\ ){v}(%\\ \\|\\ \\[11\\.03%,\\ 74\\.04)', '{:.2f}'),
    ('pancancer_partition_share_lo_pct', 'audit sweep :2511 pancancer_partition_share_lo_pct',
     '(\\|\\ pan\\-cancer\\ \\|\\ 25\\.67%\\ \\|\\ \\[){v}(%,\\ 74\\.04%\\]\\ \\|)', '{:.2f}'),
    ('pancancer_partition_share_hi_pct', 'audit sweep :2511 pancancer_partition_share_hi_pct',
     '(an\\-cancer\\ \\|\\ 25\\.67%\\ \\|\\ \\[11\\.03%,\\ ){v}(%\\]\\ \\|)', '{:.2f}'),
    ('pancancer_partition_var_pct', 'audit sweep :2513 pancancer_partition_var_pct',
     '(Those\\ overlap\\ too\\.\\ The\\ ){v}(%\\-vs\\-8\\.38%\\ contras)', '{:.2f}'),
    ('nsclc_partition_var_pct', 'audit sweep :2513 nsclc_partition_var_pct',
     '(se\\ overlap\\ too\\.\\ The\\ 25\\.67%\\-vs\\-){v}(%\\ contrast\\ is\\ a\\ co)', '{:.2f}'),
    ('bootstrap_se_ratio', 'audit sweep :2527 bootstrap_se_ratio',
     '(\\|\\ bootstrap\\ component\\ \\|\\ ){v}(\\ \\|\\ 2\\.7556\\ \\|\\ 0\\.961\\ )', '{:.4f}'),
    ('cohort_sqrt_n_ratio', 'audit sweep :2527 cohort_sqrt_n_ratio',
     '(ootstrap\\ component\\ \\|\\ 2\\.6490\\ \\|\\ ){v}(\\ \\|\\ 0\\.961\\ \\|)', '{:.4f}'),
    ('partition_sd_ratio', 'audit sweep :2528 partition_sd_ratio',
     '(\\|\\ partition\\ component\\ \\|\\ ){v}(\\ \\|\\ 2\\.7556\\ \\|\\ 0\\.495\\ )', '{:.4f}'),
    ('cohort_sqrt_n_ratio', 'audit sweep :2528 cohort_sqrt_n_ratio',
     '(artition\\ component\\ \\|\\ 1\\.3630\\ \\|\\ ){v}(\\ \\|\\ 0\\.495\\ \\|)', '{:.4f}'),
    ('partition_sd_ratio', 'audit sweep :2529 partition_sd_ratio',
     '(\\|\\ partition\\ component\\ \\|\\ ){v}(\\ \\|\\ 3\\.0171\\ \\(sqrt\\ of)', '{:.4f}'),
    ('pv24_sqrt_sites', 'audit sweep :2529 pv24_sqrt_sites',
     '(artition\\ component\\ \\|\\ 1\\.3630\\ \\|\\ ){v}(\\ \\(sqrt\\ of\\ the\\ site)', '{:.4f}'),
    ('pv24_sqrt_sites', 'audit sweep :2556 pv24_sqrt_sites',
     '(\\ prediction\\ is\\ sqrt\\(619/68\\)\\ =\\ ){v}(\\ and\\ observed/pred)', '{:.4f}'),
    ('pv24_nsclc_sd', 'audit sweep :2595 pv24_nsclc_sd',
     '(\\|\\ NSCLC\\ \\|\\ 24\\ \\|\\ ){v}(\\ \\|\\ \\[0\\.005952,\\ 0\\.01)', '{:.6f}'),
    ('pv24_nsclc_sd_lo', 'audit sweep :2595 pv24_nsclc_sd_lo',
     '(\\|\\ NSCLC\\ \\|\\ 24\\ \\|\\ 0\\.007658\\ \\|\\ \\[){v}(,\\ 0\\.010743\\]\\ \\|\\ 1\\.80)', '{:.6f}'),
    ('pv24_nsclc_sd_hi', 'audit sweep :2595 pv24_nsclc_sd_hi',
     '(\\ \\|\\ 24\\ \\|\\ 0\\.007658\\ \\|\\ \\[0\\.005952,\\ ){v}(\\]\\ \\|\\ 1\\.80x\\ \\|\\ 6\\.37%\\ )', '{:.6f}'),
    ('pv24_nsclc_share_pct', 'audit sweep :2595 pv24_nsclc_share_pct',
     '(0\\.005952,\\ 0\\.010743\\]\\ \\|\\ 1\\.80x\\ \\|\\ ){v}(%\\ \\[3\\.95%,\\ 11\\.80%\\]\\ )', '{:.2f}'),
    ('pv24_nsclc_share_lo_pct', 'audit sweep :2595 pv24_nsclc_share_lo_pct',
     '(2,\\ 0\\.010743\\]\\ \\|\\ 1\\.80x\\ \\|\\ 6\\.37%\\ \\[){v}(%,\\ 11\\.80%\\]\\ \\|\\ 3\\.34%)', '{:.2f}'),
    ('pv24_nsclc_share_hi_pct', 'audit sweep :2595 pv24_nsclc_share_hi_pct',
     '(0743\\]\\ \\|\\ 1\\.80x\\ \\|\\ 6\\.37%\\ \\[3\\.95%,\\ ){v}(%\\]\\ \\|\\ 3\\.34%\\ \\|)', '{:.2f}'),
    ('pv24_nsclc_understated_pct', 'audit sweep :2595 pv24_nsclc_understated_pct',
     '(80x\\ \\|\\ 6\\.37%\\ \\[3\\.95%,\\ 11\\.80%\\]\\ \\|\\ ){v}(%\\ \\|)', '{:.2f}'),
    ('pv24_pancancer_sd', 'audit sweep :2596 pv24_pancancer_sd',
     '(\\|\\ pan\\-cancer\\ \\|\\ 24\\ \\|\\ ){v}(\\ \\|\\ \\[0\\.005616,\\ 0\\.01)', '{:.6f}'),
    ('pv24_pancancer_sd_lo', 'audit sweep :2596 pv24_pancancer_sd_lo',
     '(pan\\-cancer\\ \\|\\ 24\\ \\|\\ 0\\.007226\\ \\|\\ \\[){v}(,\\ 0\\.010136\\]\\ \\|\\ 1\\.80)', '{:.6f}'),
    ('pv24_pancancer_sd_hi', 'audit sweep :2596 pv24_pancancer_sd_hi',
     '(\\ \\|\\ 24\\ \\|\\ 0\\.007226\\ \\|\\ \\[0\\.005616,\\ ){v}(\\]\\ \\|\\ 1\\.80x\\ \\|\\ 29\\.86%)', '{:.6f}'),
    ('pv24_pancancer_share_pct', 'audit sweep :2596 pv24_pancancer_share_pct',
     '(0\\.005616,\\ 0\\.010136\\]\\ \\|\\ 1\\.80x\\ \\|\\ ){v}(%\\ \\[20\\.46%,\\ 45\\.59%\\])', '{:.2f}'),
    ('pv24_pancancer_share_lo_pct', 'audit sweep :2596 pv24_pancancer_share_lo_pct',
     '(,\\ 0\\.010136\\]\\ \\|\\ 1\\.80x\\ \\|\\ 29\\.86%\\ \\[){v}(%,\\ 45\\.59%\\]\\ \\|\\ 19\\.41)', '{:.2f}'),
    ('pv24_pancancer_share_hi_pct', 'audit sweep :2596 pv24_pancancer_share_hi_pct',
     '(36\\]\\ \\|\\ 1\\.80x\\ \\|\\ 29\\.86%\\ \\[20\\.46%,\\ ){v}(%\\]\\ \\|\\ 19\\.41%\\ \\|)', '{:.2f}'),
    ('pv24_pancancer_understated_pct', 'audit sweep :2596 pv24_pancancer_understated_pct',
     '(x\\ \\|\\ 29\\.86%\\ \\[20\\.46%,\\ 45\\.59%\\]\\ \\|\\ ){v}(%\\ \\|)', '{:.2f}'),
    ('pv24_bootstrap_ratio', 'audit sweep :2600 pv24_bootstrap_ratio',
     '(\\|\\ bootstrap\\ component\\ \\|\\ ){v}(\\ \\|\\ —\\ \\|\\ 2\\.7556\\ \\(obs)', '{:.4f}'),
    ('pv24_sqrt_patients', 'audit sweep :2600 pv24_sqrt_patients',
     '(trap\\ component\\ \\|\\ 2\\.6520\\ \\|\\ —\\ \\|\\ ){v}(\\ \\(observed/predict)', '{:.4f}'),
    ('pv24_partition_ratio', 'audit sweep :2601 pv24_partition_ratio',
     '(\\|\\ partition\\ component\\ \\|\\ ){v}(\\ \\|\\ \\[0\\.6971,\\ 1\\.6115)', '{:.4f}'),
    ('pv24_partition_ratio_lo', 'audit sweep :2601 pv24_partition_ratio_lo',
     '(rtition\\ component\\ \\|\\ 1\\.0599\\ \\|\\ \\[){v}(,\\ 1\\.6115\\]\\ \\(F,\\ 23\\ a)', '{:.4f}'),
    ('pv24_partition_ratio_hi', 'audit sweep :2601 pv24_partition_ratio_hi',
     '(component\\ \\|\\ 1\\.0599\\ \\|\\ \\[0\\.6971,\\ ){v}(\\]\\ \\(F,\\ 23\\ and\\ 23\\ df)', '{:.4f}'),
    ('pv24_sqrt_patients', 'audit sweep :2601 pv24_sqrt_patients',
     '(,\\ 1\\.6115\\]\\ \\(F,\\ 23\\ and\\ 23\\ df\\)\\ \\|\\ ){v}(\\ —\\ outside\\ \\|\\ 3\\.017)', '{:.4f}'),
    ('pv24_sqrt_sites', 'audit sweep :2601 pv24_sqrt_sites',
     '(d\\ 23\\ df\\)\\ \\|\\ 2\\.7556\\ —\\ outside\\ \\|\\ ){v}(\\ —\\ outside\\ \\|)', '{:.4f}'),
    ('partition_sd_ratio', 'audit sweep :2603 partition_sd_ratio',
     '(\\ point\\ estimates\\ above\\ \\(ratio\\ ){v}(,\\ "about\\ half\\ the)', '{:.4f}'),
    ('pv24_nsclc_understated_pct', 'audit sweep :2622 pv24_nsclc_understated_pct',
     '(erstates\\ total\\ uncertainty\\ by\\ ){v}(%\\ in\\ NSCLC\\ and\\ 19\\.)', '{:.1f}'),
    ('pv24_pancancer_understated_pct', 'audit sweep :2622 pv24_pancancer_understated_pct',
     '(ertainty\\ by\\ 3\\.3%\\ in\\ NSCLC\\ and\\ ){v}(%)', '{:.1f}'),
    # Second wave: results bound to authorities added for them (session 60).
    ('sortaudit_dmae', 'audit sweep :988 sortaudit_dmae',
     '(\\ \\|\\ 0\\.034\\ /\\ 0\\.014\\ \\|\\ 0\\.00481\\ vs\\ ){v}(\\ \\|)', '{:.5f}'),
    ('sortaudit_mean_dr', 'audit sweep :993 sortaudit_mean_dr',
     '(split\\ would;\\ mean\\ Δr\\ is\\ ){v}(\\ unpinned,\\ 0\\.0168\\ )', '{:.4f}'),
    ('pinned_nsclc_mean_dr', 'audit sweep :993 pinned_nsclc_mean_dr',
     '(;\\ mean\\ Δr\\ is\\ 0\\.0205\\ unpinned,\\ ){v}(\\ pinned\\.\\ \\(2\\)\\ Pinne)', '{:.4f}'),
    ('pinned_vs_ss_dr_max', 'audit sweep :987 pinned_vs_ss_dr_max',
     '(\\ macOS\\ unpinned\\ \\|\\ identical\\ \\|\\ ){v}(\\ /\\ 0\\.016\\ \\|\\ 0\\.00428)', '{:.3f}'),
    ('pinned_vs_ss_dr_mean', 'audit sweep :987 pinned_vs_ss_dr_mean',
     '(npinned\\ \\|\\ identical\\ \\|\\ 0\\.043\\ /\\ ){v}(\\ \\|\\ 0\\.00428\\ \\[0\\.0013)', '{:.3f}'),
    ('sens_ns_ref_curated_r', 'audit sweep :1259 sens_ns_ref_curated_r',
     '(22\\ \\[0\\.2355,\\ 0\\.3527\\]\\ \\|\\ 16/16\\ \\|\\ ){v}(\\ \\|\\ 0\\.079\\ \\|\\ 0\\.835\\ \\|)', '{:.3f}'),
    ('sens_ns_ref_null_r', 'audit sweep :1259 sens_ns_ref_null_r',
     '(55,\\ 0\\.3527\\]\\ \\|\\ 16/16\\ \\|\\ 0\\.373\\ \\|\\ ){v}(\\ \\|\\ 0\\.835\\ \\|\\ 0\\.300\\ \\|)', '{:.3f}'),
    ('sens_ns_ref_alpha_null', 'audit sweep :1259 sens_ns_ref_alpha_null',
     '(27\\]\\ \\|\\ 16/16\\ \\|\\ 0\\.373\\ \\|\\ 0\\.079\\ \\|\\ ){v}(\\ \\|\\ 0\\.300\\ \\|)', '{:.3f}'),
    ('sens_ns_ref_unc_excess', 'audit sweep :1259 sens_ns_ref_unc_excess',
     '(/16\\ \\|\\ 0\\.373\\ \\|\\ 0\\.079\\ \\|\\ 0\\.835\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('sens_ns_global_alpha_null', 'audit sweep :1260 sens_ns_global_alpha_null',
     '(466\\]\\ \\|\\ 1/16\\ \\|\\ 0\\.351\\ \\|\\ 0\\.361\\ \\|\\ ){v}(\\ \\|\\ −0\\.004\\ \\|)', '{:.3f}'),
    ('sens_ns_global_unc_excess_abs', 'audit sweep :1260 sens_ns_global_unc_excess_abs',
     '(16\\ \\|\\ 0\\.351\\ \\|\\ 0\\.361\\ \\|\\ 0\\.944\\ \\|\\ −){v}(\\ \\|)', '{:.3f}'),
    ('sens_ns_degen_null_r_min', 'audit sweep :1269 sens_ns_degen_null_r_min',
     '(ce\\)\\.\\ Their\\ null\\ means\\ are\\ r\\ =\\ ){v}(–0\\.37\\ with\\ sd\\ 0\\.01)', '{:.2f}'),
    ('sens_ns_degen_null_r_max', 'audit sweep :1269 sens_ns_degen_null_r_max',
     '(Their\\ null\\ means\\ are\\ r\\ =\\ 0\\.34–){v}(\\ with\\ sd\\ 0\\.018–0\\.0)', '{:.2f}'),
    ('sens_ns_degen_null_sd_min', 'audit sweep :1269 sens_ns_degen_null_sd_min',
     '(ans\\ are\\ r\\ =\\ 0\\.34–0\\.37\\ with\\ sd\\ ){v}(–0\\.020,\\ the)', '{:.3f}'),
    ('sens_ns_degen_null_sd_max', 'audit sweep :1269 sens_ns_degen_null_sd_max',
     '(e\\ r\\ =\\ 0\\.34–0\\.37\\ with\\ sd\\ 0\\.018–){v}(,\\ the)', '{:.3f}'),
    ('sens_pan_ref_null_sd', 'audit sweep :1284 sens_pan_ref_null_sd',
     '(89\\]\\ \\|\\ 16/16\\ \\|\\ 0\\.391\\ \\|\\ 0\\.100\\ \\|\\ ){v}(\\ \\|\\ 0\\.801\\ \\|\\ 0\\.306\\ \\|)', '{:.3f}'),
    ('sens_pan_ref_alpha_null', 'audit sweep :1284 sens_pan_ref_alpha_null',
     '(/16\\ \\|\\ 0\\.391\\ \\|\\ 0\\.100\\ \\|\\ 0\\.025\\ \\|\\ ){v}(\\ \\|\\ 0\\.306\\ \\|)', '{:.3f}'),
    ('sens_pan_ref_unc_excess', 'audit sweep :1284 sens_pan_ref_unc_excess',
     '(391\\ \\|\\ 0\\.100\\ \\|\\ 0\\.025\\ \\|\\ 0\\.801\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('sens_pan_global_null_sd', 'audit sweep :1285 sens_pan_global_null_sd',
     '(53\\]\\ \\|\\ 15/16\\ \\|\\ 0\\.725\\ \\|\\ 0\\.452\\ \\|\\ ){v}(\\ \\|\\ 0\\.805\\ \\|\\ 0\\.424\\ \\|)', '{:.3f}'),
    ('sens_pan_global_alpha_null', 'audit sweep :1285 sens_pan_global_alpha_null',
     '(/16\\ \\|\\ 0\\.725\\ \\|\\ 0\\.452\\ \\|\\ 0\\.069\\ \\|\\ ){v}(\\ \\|\\ 0\\.424\\ \\|)', '{:.3f}'),
    ('sens_pan_global_unc_excess', 'audit sweep :1285 sens_pan_global_unc_excess',
     '(725\\ \\|\\ 0\\.452\\ \\|\\ 0\\.069\\ \\|\\ 0\\.805\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('sens_pan_global_axis_var_pct', 'audit sweep :1288 sens_pan_global_axis_var_pct',
     '(global\\ first\\ component\\ \\(){v}(%\\ of\\ expression\\ va)', '{:.1f}'),
    ('sens_pan_global_shift', 'audit sweep :1291 sens_pan_global_shift',
     '(excess\\ grows\\ by\\ ){v}(,\\ its\\ interval\\ wid)', '{:.3f}'),
    ('sens_pan_reference_width', 'audit sweep :1291 sens_pan_reference_width',
     '(054,\\ its\\ interval\\ widens\\ from\\ ){v}(\\ to\\ 0\\.166,\\ one\\ sig)', '{:.3f}'),
    ('sens_pan_global_axis_width', 'audit sweep :1291 sens_pan_global_axis_width',
     '(interval\\ widens\\ from\\ 0\\.044\\ to\\ ){v}(,\\ one\\ signature)', '{:.3f}'),
    ('sens_ns_null_boot_400_move', 'audit sweep :1367 sens_ns_null_boot_400_move',
     '(6\\ \\|\\ interval\\ moves\\ by\\ at\\ most\\ ){v}(\\ \\|)', '{:.4f}'),
    ('sens_ns_null_boot_1000_move', 'audit sweep :1368 sens_ns_null_boot_1000_move',
     '(33,\\ 0\\.3511\\]\\ \\|\\ 16/16\\ \\|\\ at\\ most\\ ){v}(;\\ the\\ registered\\ 2)', '{:.4f}'),
    ('sens_ns_k3_shift_abs', 'audit sweep :1369 sens_ns_k3_shift_abs',
     '(8\\ \\[0\\.1673,\\ 0\\.2722\\]\\ \\|\\ 15/16\\ \\|\\ −){v}(\\ \\|)', '{:.3f}'),
    ('sens_ns_k10_shift', 'audit sweep :1370 sens_ns_k10_shift',
     '(4\\ \\[0\\.2607,\\ 0\\.3711\\]\\ \\|\\ 16/16\\ \\|\\ \\+){v}(\\ \\|)', '{:.3f}'),
    ('e9_nsclc_alpha_min', 'audit sweep :1387 e9_nsclc_alpha_min',
     '(LC\\ \\(task\\ 0,\\ 258\\ s\\)\\ \\|\\ 16,000\\ \\|\\ ){v}(\\ \\|\\ 0\\.226\\ \\|\\ 0\\.553\\ \\|)', '{:.4f}'),
    ('e9_nsclc_p01_min', 'audit sweep :1387 e9_nsclc_p01_min',
     '(0,\\ 258\\ s\\)\\ \\|\\ 16,000\\ \\|\\ 0\\.0586\\ \\|\\ ){v}(\\ \\|\\ 0\\.553\\ \\|\\ 0\\ \\|\\ 0\\ \\|)', '{:.3f}'),
    ('e9_nsclc_ratio_max', 'audit sweep :1387 e9_nsclc_ratio_max',
     '(\\)\\ \\|\\ 16,000\\ \\|\\ 0\\.0586\\ \\|\\ 0\\.226\\ \\|\\ ){v}(\\ \\|\\ 0\\ \\|\\ 0\\ \\|)', '{:.3f}'),
    ('e9_pancancer_alpha_min', 'audit sweep :1388 e9_pancancer_alpha_min',
     '(\\ \\(task\\ 1,\\ 4,148\\ s\\)\\ \\|\\ 16,000\\ \\|\\ ){v}(\\ \\|\\ 0\\.243\\ \\|\\ 0\\.390\\ \\|)', '{:.3f}'),
    ('e9_pancancer_p01_min', 'audit sweep :1388 e9_pancancer_p01_min',
     '(,\\ 4,148\\ s\\)\\ \\|\\ 16,000\\ \\|\\ 0\\.172\\ \\|\\ ){v}(\\ \\|\\ 0\\.390\\ \\|\\ 0\\ \\|\\ 0\\ \\|)', '{:.3f}'),
    ('e9_pancancer_ratio_max', 'audit sweep :1388 e9_pancancer_ratio_max',
     '(s\\)\\ \\|\\ 16,000\\ \\|\\ 0\\.172\\ \\|\\ 0\\.243\\ \\|\\ ){v}(\\ \\|\\ 0\\ \\|\\ 0\\ \\|)', '{:.3f}'),
    ('e9_nsclc_floor_multiple', 'audit sweep :1390 e9_nsclc_floor_multiple',
     '(The\\ nearest\\ draw\\ is\\ ){v}(\\ times\\ the\\ reliabi)', '{:.1f}'),
    ('lsv_pancancer_adj_r2_site_only', 'audit sweep :1467 lsv_pancancer_adj_r2_site_only',
     '(ted\\ for\\ its\\ number\\ of\\ terms\\ \\|\\ ){v}(\\ \\|\\ 0\\.112\\ \\|)', '{:.3f}'),
    ('lsv_nsclc_adj_r2_site_only', 'audit sweep :1467 lsv_nsclc_adj_r2_site_only',
     '(its\\ number\\ of\\ terms\\ \\|\\ 0\\.424\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('lsv_pancancer_r2_type_only', 'audit sweep :1468 lsv_pancancer_r2_type_only',
     '(\\|\\ cancer\\ type\\ alone\\ \\|\\ ){v}(\\ \\|\\ 0\\.036\\ \\|)', '{:.3f}'),
    ('lsv_nsclc_r2_type_only', 'audit sweep :1468 lsv_nsclc_r2_type_only',
     '(\\|\\ cancer\\ type\\ alone\\ \\|\\ 0\\.385\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('lsv_nsclc_r2_site_given_type', 'audit sweep :1469 lsv_nsclc_r2_site_given_type',
     '(e\\ given\\ cancer\\ type\\ \\|\\ 0\\.055\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('e15_pancancer_median_observed_r2', 'audit sweep :1483 e15_pancancer_median_observed_r2',
     '(\\ 10\\ patients:\\ median\\ crude\\ R²\\ ){v}(\\ and\\ 0\\.177\\ on\\ that)', '{:.3f}'),
    ('e15_nsclc_median_observed_r2', 'audit sweep :1483 e15_nsclc_median_observed_r2',
     '(ts:\\ median\\ crude\\ R²\\ 0\\.485\\ and\\ ){v}(\\ on\\ that\\ wider\\ fra)', '{:.3f}'),
    ('control_c_probe_observed_r2', 'audit sweep :2311 control_c_probe_observed_r2',
     '(\\|\\ Observed\\ \\|\\ ){v}(\\ \\|)', '{:.4f}'),
    ('control_c_probe_perm_median_r2', 'audit sweep :2312 control_c_probe_perm_median_r2',
     '(\\|\\ Permuted\\ median\\ \\|\\ ){v}(\\ \\|)', '{:.4f}'),
    ('control_c_probe_perm_p', 'audit sweep :2313 control_c_probe_perm_p',
     '(alibrated\\*\\*\\ \\|\\ \\*\\*0\\.0410\\*\\*\\ \\(p\\ =\\ ){v}(\\)\\ \\|)', '{:.5f}'),
    ('nsclc_reported_ci_width', 'audit sweep :2405 nsclc_reported_ci_width',
     '(val\\ \\|\\ \\[0\\.2542,\\ 0\\.3695\\],\\ width\\ ){v}(\\ \\|)', '{:.4f}'),
    ('nsclc_honest_ci_width', 'audit sweep :2406 nsclc_honest_ci_width',
     '(\\*\\ \\|\\ \\*\\*\\[0\\.2516,\\ 0\\.3720\\],\\ width\\ ){v}(\\*\\*\\ \\|)', '{:.4f}'),
    ('pv_pancancer_seed1', 'audit sweep :2420 pv_pancancer_seed1',
     '(\\|\\ 1\\ \\|\\ ){v}(\\ \\|\\ 3,541\\ s\\ \\|)', '{:.4f}'),
    ('pv_pancancer_seed2', 'audit sweep :2421 pv_pancancer_seed2',
     '(\\|\\ 2\\ \\|\\ ){v}(\\ \\|\\ 11,404\\ s\\ \\|)', '{:.4f}'),
    ('pv_pancancer_seed3', 'audit sweep :2422 pv_pancancer_seed3',
     '(\\|\\ 3\\ \\|\\ ){v}(\\ \\|\\ 3,374\\ s\\ \\|)', '{:.4f}'),
    ('pv_pancancer_seed4', 'audit sweep :2423 pv_pancancer_seed4',
     '(\\|\\ 4\\ \\|\\ ){v}(\\ \\|\\ 3,364\\ s\\ \\|)', '{:.4f}'),
    ('pancancer_reported_ci_width', 'audit sweep :2436 pancancer_reported_ci_width',
     '(val\\ \\|\\ \\[0\\.2717,\\ 0\\.3152\\],\\ width\\ ){v}(\\ \\|)', '{:.4f}'),
    ('pancancer_honest_ci_width', 'audit sweep :2437 pancancer_honest_ci_width',
     '(\\*\\ \\|\\ \\*\\*\\[0\\.2683,\\ 0\\.3187\\],\\ width\\ ){v}(\\*\\*\\ \\|)', '{:.4f}'),
    ('bootstrap_ratio_over_sqrt_n', 'audit sweep :2527 bootstrap_ratio_over_sqrt_n',
     '(component\\ \\|\\ 2\\.6490\\ \\|\\ 2\\.7556\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('partition_ratio_over_sqrt_n', 'audit sweep :2528 partition_ratio_over_sqrt_n',
     '(component\\ \\|\\ 1\\.3630\\ \\|\\ 2\\.7556\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('partition_ratio_over_sqrt_sites', 'audit sweep :2529 partition_ratio_over_sqrt_sites',
     '(the\\ site\\ counts,\\ 619\\ vs\\ 68\\)\\ \\|\\ ){v}(\\ \\|)', '{:.3f}'),
    ('pv24_bootstrap_over_sqrt_n', 'audit sweep :2600 pv24_bootstrap_over_sqrt_n',
     '(\\ \\|\\ 2\\.7556\\ \\(observed/predicted\\ ){v}(\\)\\ \\|\\ —\\ \\|)', '{:.3f}'),
    # Third wave (session 60).
    ('sortaudit_shift', 'audit sweep :734 sortaudit_shift',
     '(h,\\ and\\ about\\ a\\ quarter\\ of\\ the\\ ){v}(\\ that\\ the\\ same\\ cla)', '{:.3f}'),
    ('a10ss_shift', 'audit sweep :736 a10ss_shift',
     '(res\\ therefore\\ aggregates\\ to\\ a\\ ){v}(\\ movement\\ in\\ the\\ i)', '{:.3f}'),
    ('sortaudit_macos_elements', 'audit sweep :695 sortaudit_macos_elements',
     '(ank\\ entries\\ \\|\\ \\*\\*35,394,448\\ of\\ ){v}(\\ \\(91%\\)\\*\\*\\ \\|)', '{:,.0f}'),
    ('e11_excess_difference', 'audit sweep :1407 e11_excess_difference',
     '(mean\\ axis\\-residualised\\ excess\\ ){v}(\\ above\\ the\\ other\\ t)', '{:.4f}'),
    ('lsv_crude_gap', 'audit sweep :1473 lsv_crude_gap',
     '(nd\\ it\\ cannot\\ explain\\ a\\ gap\\ of\\ ){v}(\\.)', '{:.2f}'),
    ('control_c_probe_observed_r2', 'audit sweep :2315 control_c_probe_observed_r2',
     '(The\\ crude\\ ){v}(\\ is\\ almost\\ entirel)', '{:.2f}'),
    ('pv24_pancancer_sd', 'audit sweep :2607 pv24_pancancer_sd',
     '(between\\-partition\\ sd\\ is\\ about\\ ){v}(\\ to\\ 0\\.008\\ in\\ both)', '{:.3f}'),
    ('pv24_nsclc_sd', 'audit sweep :2607 pv24_nsclc_sd',
     '(artition\\ sd\\ is\\ about\\ 0\\.007\\ to\\ ){v}(\\ in\\ both)', '{:.3f}'),
    ('e6_meanz_recon', 'audit sweep :1574 e6_meanz_recon',
     '(\\ split\\-half\\ reconstruction\\ \\|\\ \\+){v}(\\ \\|\\ 16/16\\ \\|)', '{:.4f}'),
    # Fourth wave: current results inside sections otherwise classified as records.
    ('sffc_pc_pct_moved', 'audit sweep :1785 sffc_pc_pct_moved',
     '(\\ \\ \\ the\\ fix\\ moved\\ \\*\\*){v}(%\\*\\*\\ of\\ pan\\-cancer\\ )', '{:.1f}'),
    ('sffc_pc_moved', 'audit sweep :1785 sffc_pc_moved',
     '(\\.2%\\*\\*\\ of\\ pan\\-cancer\\ patients\\ \\(){v}(\\ of\\ 7,168;\\ adjuste)', '{:,.0f}'),
    ('pancancer_n', 'audit sweep :1785 pancancer_n',
     '(pan\\-cancer\\ patients\\ \\(3,812\\ of\\ ){v}(;\\ adjusted)', '{:,.0f}'),
    ('sffc_pc_ari', 'audit sweep :1786 sffc_pc_ari',
     '(\\ \\ \\ Rand\\ index\\ ){v}(,\\ a\\ nearly\\ unrelat)', '{:.3f}'),
    ('sffc_ns_pct_moved', 'audit sweep :1786 sffc_ns_pct_moved',
     '(nrelated\\ partition\\)\\ against\\ \\*\\*){v}(%\\*\\*\\ of\\ NSCLC)', '{:.1f}'),
    ('sffc_ns_moved', 'audit sweep :1787 sffc_ns_moved',
     '(\\ \\ \\ patients\\ \\(){v}(\\ of\\ 944;\\ ARI\\ 0\\.537)', '{:.0f}'),
    ('nsclc_n', 'audit sweep :1787 nsclc_n',
     '(\\ \\ \\ patients\\ \\(216\\ of\\ ){v}(;\\ ARI\\ 0\\.537\\)\\.\\ The\\ )', '{:.0f}'),
    ('sffc_ns_ari', 'audit sweep :1787 sffc_ns_ari',
     '(\\ \\ \\ patients\\ \\(216\\ of\\ 944;\\ ARI\\ ){v}(\\)\\.\\ The\\ mechanism\\ i)', '{:.3f}'),
    ('sffc_pc_sites_tied', 'audit sweep :1789 sffc_pc_sites_tied',
     '(th\\ cohorts\\),\\ and\\ pan\\-TCGA\\ has\\ ){v}(\\ of\\ 619\\ sites)', '{:.0f}'),
    ('pancancer_n_sites', 'audit sweep :1789 pancancer_n_sites',
     '(rts\\),\\ and\\ pan\\-TCGA\\ has\\ 600\\ of\\ ){v}(\\ sites)', '{:.0f}'),
    ('sffc_pc_pct_tied', 'audit sweep :1790 sffc_pc_pct_tied',
     '(\\ \\ \\ tied,\\ holding\\ ){v}(%\\ of\\ patients,\\ aga)', '{:.1f}'),
    ('sffc_ns_pct_tied', 'audit sweep :1790 sffc_ns_pct_tied',
     '(atients,\\ against\\ 51\\ of\\ 68\\ and\\ ){v}(%\\.\\ So\\ the\\ larger)', '{:.1f}'),
    ('sortaudit_shift', 'audit sweep :1794 sortaudit_shift',
     '(\\ \\ \\ move\\ of\\ ){v}(\\)\\ overstated\\ the\\ c)', '{:.4f}'),
    ('sffc_ns_n_pairs', 'audit sweep :1794 sffc_ns_n_pairs',
     '(stated\\ the\\ comparison\\.\\ Of\\ the\\ ){v}(\\ pairs\\ among\\ the\\ 2)', '{:.0f}'),
    ('sffc_pc_pairs_ge_move', 'audit sweep :1796 sffc_pc_pairs_ge_move',
     '(cer\\ move\\ \\(0\\.0057\\)\\ is\\ typical:\\ ){v}(\\ of\\ 276\\ pairs\\ diff)', '{:.0f}'),
    ('sffc_pc_n_pairs', 'audit sweep :1796 sffc_pc_n_pairs',
     '(e\\ \\(0\\.0057\\)\\ is\\ typical:\\ 153\\ of\\ ){v}(\\ pairs\\ differ\\ by)', '{:.0f}'),
    ('sffc_gap_sds_small', 'audit sweep :1800 sffc_gap_sds_small',
     '(uthority:\\ the\\ platform\\ gap\\ is\\ ){v}(\\ standard\\ deviatio)', '{:.3f}'),
    ('sffc_contrast_small', 'audit sweep :1802 sffc_contrast_small',
     '(C\\-over\\-pan\\-cancer\\ contrast\\ is\\ ){v}(,\\ where\\ it\\ wrote\\ 3)', '{:.3f}'),
    ('splitdet_hpc4_moved', 'audit sweep :1815 splitdet_hpc4_moved',
     '(\\ \\ \\ own\\ is\\ ){v}(,\\ and\\ no\\ machine\\-a)', '{:.0f}'),
    ('e16_ns_registered_f100_median', 'audit sweep :1819 e16_ns_registered_f100_median',
     '(its\\ full\\-cohort\\ NSCLC\\ median\\ \\(){v}(\\)\\ differs)', '{:.3f}'),
    ('nsclc_site_auroc', 'audit sweep :1820 nsclc_site_auroc',
     '(\\ \\ \\ from\\ the\\ frozen\\ macOS\\ ){v}(\\ because\\ the\\ site\\ )', '{:.3f}'),
    ('pancancer_alpha_obs_resid', 'audit sweep :2086 pancancer_alpha_obs_resid',
     '(the\\ median\\ curated\\ signature:\\ ){v}(\\ pan\\-cancer,\\ 0\\.977)', '{:.3f}'),
    ('nsclc_alpha_obs_resid', 'audit sweep :2086 nsclc_alpha_obs_resid',
     '(\\ signature:\\ 0\\.976\\ pan\\-cancer,\\ ){v}(\\ NSCLC\\.\\ It\\ is\\ not)', '{:.3f}'),
    ('pancancer_outcome_excess_g2m', 'audit sweep :2911 pancancer_outcome_excess_g2m',
     '(G2M\\ falls\\ to\\ 0\\.0520,\\ against\\ ){v}(\\ and\\ 0\\.0382\\ residu)', '{:.4f}'),
    ('pancancer_outcome_excess_emt', 'audit sweep :2911 pancancer_outcome_excess_emt',
     '(to\\ 0\\.0520,\\ against\\ 0\\.0687\\ and\\ ){v}(\\ residualised\\.\\ A\\ r)', '{:.4f}'),
    ('splitdet_macos_moved', 'audit sweep :654 splitdet_macos_moved',
     '(le\\ is\\ \\*\\*0\\ of\\ 944\\*\\*,\\ down\\ from\\ ){v}(\\.\\ On\\ the\\ milliseco)', '{:.0f}'),
    ('splitdet_macos_moved', 'audit sweep :674 splitdet_macos_moved',
     '(hash\\ `e55ee4e8887b5167`,\\ and\\ "){v}(\\ of\\ 944\\ patients\\ a)', '{:.0f}'),
    ('nsclc_n', 'audit sweep :1117 nsclc_n',
     '(axis\\ \\(0\\ of\\ ){v}(\\)\\.\\ The\\ same\\ run\\ re)', '{:.0f}'),
    # E15's attrition row, bound to the live recount in `_site_counts`.
    ("e15_pan_sites_dropped", "audit E15: pan-cancer sites under 10 patients",
     r"(\| pan-cancer \| 202 of 619 \| ){v}( \(largest dropped site)", "{:.0f}"),
    ("e15_pan_patients_dropped", "audit E15: pan-cancer patients in those sites",
     r"(\(largest dropped site: 9 patients\) \| ){v}( \| 17 of the)", "{:,.0f}"),
    ("e15_pan_retained", "audit E15: pan-cancer patients retained",
     r"(\| 17 of the ){v}( retained patients have no cancer type)", "{:,.0f}"),
    ("e15_ns_patients_dropped", "audit E15: NSCLC patients in sites under 10",
     r"(\| NSCLC \| 31 of 68 \| 37 \| ){v}( \| 796 \|)", "{:.0f}"),
]

CLAIMS += [(key, _a16_claim(shape).replace(r"([\d,]+)", r"(\d[\d,]*(?:\.\d+)?)", 1), 0.0,
            f"audit figure, {what}", None)
           for key, what, shape, _fmt in AUDIT_SWEEP_SPEC]

# Session 60 (H89): A5's NSCLC table of six partitions, whose rows ("| 1 |
# 0.3255 |") no line can tell apart from another two-column table. Scoped to
# the paragraph whose header is exactly `| split_seed | ISI |` -- the pan-cancer
# table carries a runtime column, and 17-CPTAC-PLAN.md's slide table another
# header -- and held to printed precision.
CLAIMS += [(f"pv_nsclc_seed{_s}", rf"^\| {_s} \| (\d\.\d{{4}}) \|$", 5e-7,
            f"audit A5: NSCLC partition {_s} ISI",
            r"PARA:split_seed \| ISI \|\s+\|---\|---\|\s")
           for _s in range(6)]

# Session 57: the deposited README printed ComBat figures -- raw 1.000,
# in-sample 0.049, out-of-fold 0.020 -- as "Measured:", with no scope. They are
# the SYNTHETIC measurement in `models.site_prediction_control_oof_combat`'s
# docstring; no stored artefact carries any of them, and the stored out-of-fold
# runs give 0.005 (NSCLC 0.00508, demo 0.00494). The passage now scopes them
# and states the stored figure, which is bound here to the same authority the
# captions use.
# Session 57: the deposited README's demo table, bound to the artefact it
# names. A mutation sweep of the README moved each of its 15 numbers and none
# failed. Each claim anchors on its own row and column; the self-test pushes
# each cell 0.010 off its live value.
DEMO_TABLE_COLUMNS = [("r_residual", "r"), ("null_mean_r", "null"),
                      ("excess_z", "excess"), ("excess_lo", "lo"), ("excess_hi", "hi")]
_DCELL = r"[+−-]?\d\.\d{3}"
DEMO_TABLE_SPEC: list[tuple[str, str, str]] = []
for _name in ("clean", "siteconf", "purity"):
    _row = rf"^\| `SIG_immune_{_name}` \| "
    DEMO_TABLE_SPEC += [
        (f"demo_{_name}_r", _row, r" \| "),
        (f"demo_{_name}_null", _row + _DCELL + r" \| ", r" \| "),
        (f"demo_{_name}_excess", _row + _DCELL + r" \| " + _DCELL + r" \| \**", r"\** \| \["),
        (f"demo_{_name}_lo", _row + _DCELL + r" \| " + _DCELL + r" \| \**" + _DCELL
         + r"\** \| \[", r", "),
        (f"demo_{_name}_hi", _row + _DCELL + r" \| " + _DCELL + r" \| \**" + _DCELL
         + r"\** \| \[" + _DCELL + r", ", r"\] \|"),
    ]
CLAIMS += [(key, prefix + "(" + _DCELL + ")" + suffix, 5e-4,
            f"README demo table {key}", None)
           for key, prefix, suffix in DEMO_TABLE_SPEC]

# Session 57: the deposited README's sort-fix paragraph. "The gap is 0.0218"
# was checked; the ratio and the partition sd beside it, and the ratio again
# at full precision, were not -- each moved in the sweep without failing.
README_SORTFIX_SPEC: list[tuple[str, str, str, str]] = [
    ("sffc_gap_sds_small", r"It is ", r" times the standard deviation of partition", "{:.1f}"),
    ("partition_sd", r"choice over six partitions \(", r"\)", "{:.4f}"),
    ("sffc_gap_sds_small", r"interval for; ", r" to full precision", "{:.3f}"),
]
CLAIMS += [(key, prefix + r"(\d\.\d+)" + suffix, 0.0,
            f"README sort-fix paragraph, {key} as {fmt}", None)
           for key, prefix, suffix, fmt in README_SORTFIX_SPEC]

CLAIMS += [("cap_combat_auroc_median",
            r"On the frozen NSCLC run the out-of-fold median is (\d\.\d{3})", 5e-4,
            "README: stored out-of-fold ComBat site AUROC", None)]

# Session 54: THE MANUSCRIPT'S RESULTS, SWEPT THE WAY THE POSTER WAS. A mutation
# sweep of every literal `--coverage` listed as unchecked in the abstract,
# Results and Discussion found some sixty that moved without failing
# `--strict`: the partition-widened intervals and every step of their
# arithmetic, the Control C sample sizes, the size-gap Spearman statistics,
# the rotation-null family means, the multiplicity margin and M_eff, the
# tumor-only and site-AUROC counts, and the live abstract's cohort size and
# fold range, which wrap. (The bracketed triples among them -- the ssGSEA,
# PLAGE, reconstruction and global-axis arms -- are guarded by the triple rule,
# which now fails a point estimate that identifies nothing.) Each row is
# (key, claim regex, tolerance, description, scope, and the self-test's
# prefix, suffix, printed format, injected step, and whether the page prints
# the magnitude beside a separate sign). The self-test prefix and suffix are
# written with single spaces and matched with `\s+`, because most of these
# claims wrap in the file.
_MANUSCRIPT_S54: list[tuple] = [
    ("pancancer_n", r"In ([\d,]+) TCGA patients across", 0.5,
     "manuscript abstract: pan-TCGA cohort size", "WRAP:",
     r"In ", r" TCGA patients", "{:,.0f}", 100, False),
    ("nsclc_rel_gap_fold", r"gap (\d+\.\d)- to \d+-fold larger", 0.05,
     "manuscript abstract: curated-random gap fold, NSCLC", "WRAP:",
     r"random gap ", r"- to", "{:.1f}", 0.6, False),
    ("pancancer_rel_gap_fold", r"gap \d+\.\d- to (\d+)-fold larger", 0.5,
     "manuscript abstract: curated-random gap fold, pan-TCGA", "WRAP:",
     r"- to ", r"-fold larger", "{:.0f}", 2, False),
    ("nsclc_rel_gap_fold", r"gap was (\d+\.\d)-fold \(NSCLC\)", 0.05,
     "superseded structured abstract: gap fold, NSCLC", "WRAP:",
     r"gap was ", r"-fold \(NSCLC\)", "{:.1f}", 0.6, False),
    ("pancancer_rel_gap_fold", r"\(NSCLC\) to (\d+)-fold \(pan-cancer\)", 0.5,
     "superseded structured abstract: gap fold, pan-TCGA", "WRAP:",
     r"\(NSCLC\) to ", r"-fold \(pan-cancer\)", "{:.0f}", 2, False),
    ("small_panel_spearman", r"at 10 \(Spearman ρ = ([−-]\d\.\d+), p", 5e-4,
     "superseded structured abstract: small-panel Spearman", "WRAP:",
     r"at 10 \(Spearman ρ = −", r", p <", "{:.3f}", -0.006, True),
    ("pancancer_relgap_size_rho", r"\(Fig\. 2A; Spearman ρ = ([−-]\d\.\d+),", 5e-3,
     "manuscript Results: Spearman(size, gap), pan-TCGA", "WRAP:",
     r"\(Fig\. 2A; Spearman ρ = −", r", p =", "{:.2f}", -0.06, True),
    ("pancancer_relgap_size_p",
     r"\(Fig\. 2A; Spearman ρ = [−-]\d\.\d+, p = (\d\.\d+) pan-cancer", 5e-4,
     "manuscript Results: Spearman(size, gap) p, pan-TCGA", "WRAP:",
     r"p = ", r" pan-cancer; ρ", "{:.3f}", 0.006, False),
    ("nsclc_relgap_size_rho", r"pan-cancer; ρ = ([−-]\d\.\d+), p = \d\.\d+ in NSCLC", 5e-3,
     "manuscript Results: Spearman(size, gap), NSCLC", "WRAP:",
     r"pan-cancer; ρ = −", r", p = [\d.]+ in NSCLC", "{:.2f}", -0.06, True),
    ("nsclc_relgap_size_p", r"pan-cancer; ρ = [−-]\d\.\d+, p = (\d\.\d+) in NSCLC", 5e-4,
     "manuscript Results: Spearman(size, gap) p, NSCLC", "WRAP:",
     r"pan-cancer; ρ = −[\d.]+, p = ", r" in NSCLC\)", "{:.3f}", 0.006, False),
    ("pancancer_rel_gap_fullsize_mean", r"it averages \+(\d\.\d{3}) across the", 5e-4,
     "manuscript Results: mean gap over the nominal-200 sets, pan-TCGA", "WRAP:",
     r"it averages \+", r" across the", "{:.3f}", 0.006, False),
    ("small_panel_ratio_k10_k160", r"p < 0\.0001\) — (\d+\.\d)-fold", 0.05,
     "manuscript Results: gap at k = 10 over k = 160", "WRAP:",
     r"p < 0\.0001\) — ", r"-fold", "{:.1f}", 0.6, False),
    ("sortaudit_macos_frac_sd", r"up to (\d\.\d{2}) of a score standard deviation", 5e-3,
     "manuscript Results: largest tie-order score move, in score sds (macOS)", "WRAP:",
     r"up to ", r" of a score standard", "{:.2f}", 0.06, False),
    ("a7_ssgsea_rel_null", r"multiplies the null by 1/√(\d\.\d{3,})", 5e-4,
     "manuscript Results: ssGSEA random-set reliability in the inflation", "WRAP:",
     r"null by 1/√", r" ≈", "{:.3f}", 0.006, False),
    ("a7_ssgsea_null_inflation", r"1/√\d\.\d{3,} ≈ (\d\.\d{2})", 5e-3,
     "manuscript Results: disattenuation inflation of the ssGSEA null", "WRAP:",
     r" ≈ ", r" and the observed", "{:.2f}", 0.06, False),
    ("a7_ssgsea_obs_inflation", r"the observed by only (\d\.\d{2})\.", 5e-3,
     "manuscript Results: disattenuation inflation of the ssGSEA observed", "WRAP:",
     r"observed by only ", r"\. Dividing", "{:.2f}", 0.06, False),
    ("rotation_observed", r"family mean correlation of (\d\.\d{3}) against", 5e-4,
     "manuscript Results: rotation null, observed family mean r", "WRAP:",
     r"family mean correlation of ", r" against a null", "{:.3f}", 0.006, False),
    ("rotation_null_mean", r"against a null family mean of (\d\.\d{3})", 5e-4,
     "manuscript Results: rotation null, null family mean r", "WRAP:",
     r"null family mean of ", r" \(95% range", "{:.3f}", 0.006, False),
    ("pancancer_hypoxia_excess_lo", r"narrowest margin \(excess_lo = \+(\d\.\d+)\)", 5e-5,
     "manuscript Results: hypoxia excess_lo, the multiplicity margin", "WRAP:",
     r"narrowest margin \(excess_lo = \+", r"\)", "{:.4f}", 0.0006, False),
    ("m_eff", r"\(11\) is (\d\.\d{2}) for the 16 signatures", 5e-3,
     "manuscript Results: effective number of tests (Li and Ji)", "WRAP:",
     r"\(11\) is ", r" for the 16 signatures", "{:.2f}", 0.4, False),
    ("cap_pan_lsv_sites", r"regression uses the (\d+) sites that contribute", 0.5,
     "manuscript Results: Control C sites, pan-TCGA", "WRAP:",
     r"regression uses the ", r" sites that contribute", "{:.0f}", 7, False),
    ("cap_pan_lsv_n", r"at least 10 patients \(([\d,]+) patients, after", 0.5,
     "manuscript Results: Control C patients, pan-TCGA", "WRAP:",
     r"at least 10 patients \(", r" patients, after", "{:,.0f}", 30, False),
    ("cap_pan_lsv_sites", r"dummy-coding (\d+) sites", 0.5,
     "manuscript Results: Control C sites dummy-coded, pan-TCGA", "WRAP:",
     r"dummy-coding ", r" sites", "{:.0f}", 7, False),
    ("nsclc_label_site_variance_pct", r"a median (\d+\.\d)% of the label across", 0.05,
     "manuscript Results: label-side site R^2, NSCLC, percent", "WRAP:",
     r"third the size: a median ", r"% of the", "{:.1f}", 0.6, False),
    ("cap_nsclc_lsv_n", r"usable sites and (\d+) patients, against", 0.5,
     "manuscript Results: Control C patients, NSCLC", "WRAP:",
     r"usable sites and ", r" patients, against", "{:.0f}", 9, False),
    ("cap_pan_lsv_sites", r"against 44% across (\d+) sites", 0.5,
     "manuscript Results: Control C sites, pan-TCGA, beside NSCLC", "WRAP:",
     r"44% across ", r" sites pan-cancer", "{:.0f}", 7, False),
    ("nsclc_reported_ci_lo", r"widens the interval from \[(\d\.\d{4}), \d\.\d{4}\] to", 5e-5,
     "manuscript Results: NSCLC reported interval, low", "WRAP:",
     r"widens the interval from \[", r", [\d.]+\] to", "{:.4f}", 0.0006, False),
    ("nsclc_reported_ci_hi", r"widens the interval from \[\d\.\d{4}, (\d\.\d{4})\] to", 5e-5,
     "manuscript Results: NSCLC reported interval, high", "WRAP:",
     r"widens the interval from \[[\d.]+, ", r"\] to", "{:.4f}", 0.0006, False),
    ("nsclc_honest_ci_lo", r"\] to \[(\d\.\d{4}), \d\.\d{4}\]: in NSCLC", 5e-5,
     "manuscript Results: NSCLC partition-widened interval, low", "WRAP:",
     r"\] to \[", r", [\d.]+\]: in NSCLC", "{:.4f}", 0.0006, False),
    ("nsclc_honest_ci_hi", r"\] to \[\d\.\d{4}, (\d\.\d{4})\]: in NSCLC", 5e-5,
     "manuscript Results: NSCLC partition-widened interval, high", "WRAP:",
     r"\] to \[[\d.]+, ", r"\]: in NSCLC", "{:.4f}", 0.0006, False),
    ("pancancer_reported_ci_lo", r"combination widens \[(\d\.\d{4}), \d\.\d{4}\] to", 5e-5,
     "manuscript Results: pan-TCGA reported interval, low", "WRAP:",
     r"combination widens \[", r", [\d.]+\] to", "{:.4f}", 0.0006, False),
    ("pancancer_reported_ci_hi", r"combination widens \[\d\.\d{4}, (\d\.\d{4})\] to", 5e-5,
     "manuscript Results: pan-TCGA reported interval, high", "WRAP:",
     r"combination widens \[[\d.]+, ", r"\] to", "{:.4f}", 0.0006, False),
    ("pancancer_honest_ci_lo",
     r"combination widens \[\d\.\d{4}, \d\.\d{4}\] to \[(\d\.\d{4}), ", 5e-5,
     "manuscript Results: pan-TCGA partition-widened interval, low", "WRAP:",
     r"combination widens \[[\d.]+, [\d.]+\] to \[", r", [\d.]+\]: the reported", "{:.4f}",
     0.0006, False),
    ("pancancer_honest_ci_hi",
     r"combination widens \[\d\.\d{4}, \d\.\d{4}\] to \[\d\.\d{4}, (\d\.\d{4})\]", 5e-5,
     "manuscript Results: pan-TCGA partition-widened interval, high", "WRAP:",
     r"to \[[\d.]+, ", r"\]: the reported interval is", "{:.4f}", 0.0006, False),
    ("pancancer_partition_var_pct",
     r"choice — (\d+\.\d)% pan-cancer against \d+\.\d% in NSCLC", 0.05,
     "manuscript Results: partition share of variance, pan-TCGA", "WRAP:",
     r"choice — ", r"% pan-cancer against", "{:.1f}", 0.6, False),
    ("nsclc_partition_var_pct",
     r"choice — \d+\.\d% pan-cancer against (\d+\.\d)% in NSCLC", 0.05,
     "manuscript Results: partition share of variance, NSCLC", "WRAP:",
     r"% pan-cancer against ", r"% in NSCLC, in the cohort", "{:.1f}", 0.6, False),
    ("cohort_n_ratio", r"in the cohort (\d+\.\d) times larger", 0.05,
     "manuscript Results: pan-TCGA over NSCLC cohort size", "WRAP:",
     r"in the cohort ", r" times larger", "{:.1f}", 0.6, False),
    ("nsclc_partition_sd_chi2_lo", r"spans \[(\d\.\d{4}), \d\.\d{4}\] in NSCLC", 5e-5,
     "manuscript Results: six-partition sd chi-square interval, NSCLC, low", "WRAP:",
     r"spans \[", r", [\d.]+\] in NSCLC", "{:.4f}", 0.0006, False),
    ("nsclc_partition_sd_chi2_hi", r"spans \[\d\.\d{4}, (\d\.\d{4})\] in NSCLC", 5e-5,
     "manuscript Results: six-partition sd chi-square interval, NSCLC, high", "WRAP:",
     r"spans \[[\d.]+, ", r"\] in NSCLC", "{:.4f}", 0.0006, False),
    ("pancancer_partition_sd_chi2_lo", r"in NSCLC and \[(\d\.\d{4}), \d\.\d{4}\] pan-cancer", 5e-5,
     "manuscript Results: five-partition sd chi-square interval, pan-TCGA, low", "WRAP:",
     r"\] in NSCLC and \[", r", [\d.]+\] pan-cancer", "{:.4f}", 0.0006, False),
    ("pancancer_partition_sd_chi2_hi", r"in NSCLC and \[\d\.\d{4}, (\d\.\d{4})\] pan-cancer", 5e-5,
     "manuscript Results: five-partition sd chi-square interval, pan-TCGA, high", "WRAP:",
     r"\] in NSCLC and \[[\d.]+, ", r"\] pan-cancer", "{:.4f}", 0.0006, False),
    ("nsclc_partition_share_lo_pct", r"become \[(\d+\.\d)%, \d+\.\d%\] and", 0.05,
     "manuscript Results: partition share interval, NSCLC, low", "WRAP:",
     r"those become \[", r"%, [\d.]+%\] and", "{:.1f}", 0.6, False),
    ("nsclc_partition_share_hi_pct", r"become \[\d+\.\d%, (\d+\.\d)%\] and", 0.05,
     "manuscript Results: partition share interval, NSCLC, high", "WRAP:",
     r"those become \[[\d.]+%, ", r"%\] and", "{:.1f}", 0.6, False),
    ("pancancer_partition_share_lo_pct", r"become \[\d+\.\d%, \d+\.\d%\] and \[(\d+\.\d)%, ", 0.05,
     "manuscript Results: partition share interval, pan-TCGA, low", "WRAP:",
     r"%\] and \[", r"%, [\d.]+%\]\. Each", "{:.1f}", 0.6, False),
    ("pancancer_partition_share_hi_pct", r"and \[\d+\.\d%, (\d+\.\d)%\]\. Each standard", 0.05,
     "manuscript Results: partition share interval, pan-TCGA, high", "WRAP:",
     r"%\] and \[[\d.]+%, ", r"%\]\. Each standard", "{:.1f}", 0.6, False),
    ("bootstrap_se_ratio", r"falls by a factor of (\d\.\d{4}), against", 5e-5,
     "manuscript Results: bootstrap SE, NSCLC over pan-TCGA", "WRAP:",
     r"falls by a factor of ", r", against the", "{:.4f}", 0.0006, False),
    ("cohort_sqrt_n_ratio", r"against the (\d\.\d{4}) predicted by the ratio of the square roots", 5e-5,
     "manuscript Results: square root of the cohort-size ratio", "WRAP:",
     r"against the ", r" predicted by the ratio", "{:.4f}", 0.0006, False),
    ("partition_sd_ratio", r"partition component falls by only (\d\.\d{4}),", 5e-5,
     "manuscript Results: partition sd, NSCLC over pan-TCGA", "WRAP:",
     r"falls by only ", r", roughly half", "{:.4f}", 0.0006, False),
    ("cohort_sqrt_n_ratio", r"by patient count \((\d\.\d{4})\)", 5e-5,
     "manuscript Results: square root of the cohort-size ratio, restated", "WRAP:",
     r"by patient count \(", r"\) or by", "{:.4f}", 0.0006, False),
    ("pv24_sqrt_sites", r"tissue-source-site count \((\d\.\d{4}), for", 5e-5,
     "manuscript Results: square root of the site-count ratio", "WRAP:",
     r"tissue-source-site count \(", r", for", "{:.4f}", 0.0006, False),
    ("pancancer_n_sites", r"\(\d\.\d{4}, for (\d+) sites against", 0.5,
     "manuscript Results: pan-TCGA tissue source sites", "WRAP:",
     r", for ", r" sites against", "{:.0f}", 7, False),
    ("nsclc_interval_understated_pct", r"total uncertainty, by (\d+\.\d)% in NSCLC and", 0.05,
     "manuscript Results: interval understatement, NSCLC", "WRAP:",
     r"uncertainty, by ", r"% in NSCLC and", "{:.1f}", 0.6, False),
    ("pancancer_interval_understated_pct",
     r"uncertainty, by \d+\.\d% in NSCLC and (\d+\.\d)% pan-cancer", 0.05,
     "manuscript Results: interval understatement, pan-TCGA", "WRAP:",
     r"% in NSCLC and ", r"% pan-cancer —", "{:.1f}", 0.6, False),
    ("nsclc_partition_var_pct", r"six-partition NSCLC share of (\d+\.\d)%", 0.05,
     "manuscript Results: six-partition NSCLC share, restated", "WRAP:",
     r"six-partition NSCLC share of ", r"% nor", "{:.1f}", 0.6, False),
    ("a16_ns_patients", r"\[\d\.\d{4}, \d\.\d{4}\] \((\d+) patients\), so", 0.5,
     "manuscript Results: NSCLC tumor-only cohort", "WRAP:",
     r"0\.4079\] \(", r" patients\), so", "{:.0f}", 5, False),
    ("nsclc_site_auroc", r"pan-cancer and (\d\.\d{3}) in NSCLC, with \d+% of evaluable", 5e-4,
     "manuscript Results: NSCLC site AUROC, median", "WRAP:",
     r"\*\* pan-cancer and ", r" in NSCLC, with", "{:.3f}", -0.006, False),
    ("e16_ns_registered_f100_frac_pct", r"in NSCLC, with (\d+)% of evaluable sites above", 0.5,
     "manuscript Results: NSCLC sites above 0.9, percent", "WRAP:",
     None, None, None, None, False),
    ("e16_pan_registered_f100_frac_pct", r"in NSCLC, with (\d+)% of evaluable sites above", 0.5,
     "manuscript Results: pan-TCGA sites above 0.9, percent", "WRAP:",
     r"in NSCLC, with ", r"% of evaluable sites", "{:.0f}", -3, False),
    ("cap_combat_auroc_median", r"anti-predictive \(measured (\d\.\d{3})\)", 5e-4,
     "manuscript Results: post-ComBat site AUROC, median", "WRAP:",
     r"anti-predictive \(measured ", r"\)", "{:.3f}", 0.006, False),
    ("e16_pan_quarter_n", r"at a quarter \(([\d,]+) patients, \d+ to \d+ evaluable", 0.5,
     "manuscript Results: pan-TCGA quarter-cohort patients", "WRAP:",
     r"at a quarter \(", r" patients, 10", "{:,.0f}", 40, False),
    ("ancestry_amr", r"AMR (\d+)\) is reported", 0.5,
     "manuscript Results: ancestry arm, AMR patients", "WRAP:",
     r"AMR ", r"\) is reported", "{:.0f}", 5, False),
    ("paper_structured_abstract_words", r"The structured (\d+)-word version", 0.5,
     "manuscript: the folded structured abstract's own length", "WRAP:",
     r"The structured ", r"-word version", "{:.0f}", 5, False),
    ("paper_structured_abstract_words", r"Superseded structured abstract \((\d+) words\)", 0.5,
     "manuscript: the folded structured abstract's length, in its summary line", None,
     r"Superseded structured abstract \(", r" words\)", "{:.0f}", 5, False),
    ("nsclc_n", r"In NSCLC, \d+ of (\d+) patients were paired", 0.5,
     "manuscript Results: NSCLC cohort, tumor-only denominator", "WRAP:",
     r"In NSCLC, 51 of ", r" patients were", "{:.0f}", 10, False),
    ("pancancer_n", r"Pan-cancer, \d+ of ([\d,]+) profiles averaged", 0.5,
     "manuscript Results: pan-TCGA cohort, tumor-only denominator", "WRAP:",
     r"Pan-cancer, 578 of ", r" profiles averaged", "{:,.0f}", 100, False),
    ("ancestry_site_frac_min_pct", r"Only (\d\.\d)–\d\.\d% of sites carry", 0.05,
     "manuscript Results: ancestry, smallest share of sites supporting a contrast",
     "WRAP:", r"Only ", r"–[\d.]+% of sites", "{:.1f}", 0.6, False),
    ("ancestry_site_frac_max_pct", r"Only \d\.\d–(\d\.\d)% of sites carry", 0.05,
     "manuscript Results: ancestry, largest share of sites supporting a contrast",
     "WRAP:", r"Only [\d.]+–", r"% of sites", "{:.1f}", 0.6, False),
    ("a7_meanz_rel_null", r"rank-walk statistic does, from (\d\.\d{3}) to", 5e-4,
     "manuscript Discussion: random-set reliability, mean-z", "WRAP:",
     r"statistic does, from ", r" to [\d.]+ here", "{:.3f}", 0.006, False),
    ("a7_ssgsea_rel_null", r"rank-walk statistic does, from \d\.\d{3} to (\d\.\d{3})", 5e-4,
     "manuscript Discussion: random-set reliability, ssGSEA", "WRAP:",
     r"statistic does, from [\d.]+ to ", r" here", "{:.3f}", 0.006, False),
]
CLAIMS += [row[:5] for row in _MANUSCRIPT_S54]

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
    # The pinned NSCLC run's median joined the set on 2026-09-24 (H86): it is a
    # reported cohort median too, bound exactly by `pinned_nsclc_site_auroc`.
    (("pancancer_site_auroc", "nsclc_site_auroc", "pinned_nsclc_site_auroc",
      "pinned_pancancer_site_auroc"),
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

def _stored_triples(every_run: dict[str, float]) -> list[tuple[str, float, float, float]]:
    """(run.block, value, ci_lo, ci_hi) for every stored run's triple."""
    out = []
    for k, v in every_run.items():
        if k.endswith(".value"):
            base = k[: -len(".value")]
            lo, hi = every_run.get(base + ".ci_lo"), every_run.get(base + ".ci_hi")
            if lo is not None and hi is not None:
                out.append((base, v, lo, hi))
    return out


# Triples quoted as a record of a measurement this repository never stored,
# each EXACTLY as written, so a changed digit is still a failure. Every entry
# says where the number lives. Keep this list short: the right fix for a new
# triple is a TRIPLES row with a stored authority.
UNSTORED_TRIPLES: dict[tuple[str, str, str], str] = {
    ("0.00481", "0.00199", "0.00748"):
        "A12: HPC4's secondary Delta-MAE; its summary.json is stored and hashed "
        "since 2026-09-24 under session48_diagnostics/, which every_stored_run() "
        "does not scan, so the literal stays listed here",
    ("0.0048", "0.0020", "0.0075"):
        "A12: the same HPC4 secondary Delta-MAE, rounded for the manuscript",
    # ("0.00428", "0.00138", "0.00703") -- the macOS pinned-ordering Delta-MAE --
    # was listed here until 2026-09-24, as "a diagnostic run whose summary was
    # not stored". The shipped pinned run stores it (results/nsclc_v3_pinned/),
    # so it now passes only by matching that stored run at printed precision.
}


# (point_key, lo_key, hi_key, label)
TRIPLES: list[tuple[str, str, str, str]] = [
    ("nsclc_isi", "nsclc_ci_lo", "nsclc_ci_hi", "NSCLC ISI"),
    ("pancancer_isi", "pancancer_ci_lo", "pancancer_ci_hi", "pan-TCGA ISI"),
    # These three were REPORTED AS UNMATCHED before their sources were traced:
    # two secondary-endpoint triples in the manuscript's Results paragraph on
    # preserved-site splitting, and the rotation null's family-mean interval in
    # 14-SCIENCE-AUDIT.md's A2 section.
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
    # Session 54, when an unidentified point became a failure: the pinned
    # ssGSEA arm (A10), the 24-partition pan-cancer sd with its chi-square
    # interval, and the tumour-only NSCLC run's one outcome gain.
    ("a10ss_ssgsea_unc", "a10ss_ssgsea_unc_lo", "a10ss_ssgsea_unc_hi",
     "A10 arm: ssGSEA, uncorrected, pinned ordering"),
    ("pv24_pancancer_sd", "pv24_pancancer_sd_lo", "pv24_pancancer_sd_hi",
     "pan-TCGA partition sd over 24 partitions"),
    ("pv24_nsclc_sd", "pv24_nsclc_sd_lo", "pv24_nsclc_sd_hi",
     "NSCLC partition sd over 24 partitions"),
    ("a16_ns_hypoxia_outcome", "a16_ns_hypoxia_outcome_lo",
     "a16_ns_hypoxia_outcome_hi",
     "A16 NSCLC tumor-only: hypoxia outcome excess, axis-residualised"),
    # Session 49 sensitivity analyses (Results, "Sensitivity to other analysis
    # choices", and Limitation 8's corrected-ordering values).
    ("sens_ns_reference", "sens_ns_reference_lo", "sens_ns_reference_hi",
     "NSCLC corrected-ordering ISI"),
    ("sens_pan_reference", "sens_pan_reference_lo", "sens_pan_reference_hi",
     "pan-TCGA corrected-ordering ISI"),
    ("sens_ns_global_axis", "sens_ns_global_axis_lo", "sens_ns_global_axis_hi",
     "E1 NSCLC, global axis"),
    ("sens_pan_global_axis", "sens_pan_global_axis_lo", "sens_pan_global_axis_hi",
     "E1 pan-TCGA, global axis"),
    ("sens_ns_unmatched_null", "sens_ns_unmatched_null_lo", "sens_ns_unmatched_null_hi",
     "E3 NSCLC, size-only null matching"),
    ("sens_ns_k3", "sens_ns_k3_lo", "sens_ns_k3_hi", "E5 NSCLC, 3 folds"),
    ("sens_ns_k10", "sens_ns_k10_lo", "sens_ns_k10_hi", "E5 NSCLC, 10 folds"),
    ("e6_plage_unc", "e6_plage_unc_lo", "e6_plage_unc_hi", "E6 PLAGE, uncorrected"),
    ("e6_meanz_unc", "e6_meanz_unc_lo", "e6_meanz_unc_hi", "E6 mean-z control, uncorrected"),
    ("e8_recon_omega", "e8_recon_omega_lo", "e8_recon_omega_hi", "E8 reconstruction under omega"),
    ("e8_recon_alpha", "e8_recon_alpha_lo", "e8_recon_alpha_hi", "E8 reconstruction under alpha"),
    ("a16_ns", "a16_ns_lo", "a16_ns_hi", "A16 NSCLC, tumor-only expression"),
    ("a16_pan", "a16_pan_lo", "a16_pan_hi", "A16 pan-TCGA, tumor-only expression"),
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
    # Added 2026-09-24 (session 59, B-21): Limitation 8's forensic note, lifted
    # out of the manuscript word for word. Every claim it carried was checked
    # inside 09-PAPER-DRAFT.md and is checked here now.
    "SUPPLEMENTARY-NOTE-1.md",
    # Added 2026-09-24 (B-21, after session 60): five Results passages and the
    # site control's stress test, lifted out of the manuscript word for word to
    # bring its counted body under 5,000 words. Every claim they carried was
    # checked inside 09-PAPER-DRAFT.md and is checked here now.
    "SUPPLEMENTARY-NOTE-2.md",
    "16-YOUR-TASKS.md", "pipeline/README.md",
    "poster/poster.html", "submission/SUBMISSION-CHECKLIST.md",
    "submission/COVER-LETTER.md", "pipeline/results/README.md",
    # CITATION.cff quotes the cohort size and cancer-type count in its abstract
    # and was not scanned by any earlier version.
    "CITATION.cff",
    # Added 2026-09-06 with S1's derived NSCLC multiplicity values. It states
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
    # Added 2026-09-08. It is the operational script for pasting the abstract
    # into the AACR portal and it QUOTES THE CHARACTER COUNT, which is the one
    # number on which submission is hard-blocked. 11-SUBMISSION-PACK.md carried a
    # test count that was wrong for ~thirty tests precisely because nothing
    # scanned it; a submission-day document is the last place to repeat that.
    "submission/ABSTRACT-PORTAL-DRY-RUN.md",
    # Added 2026-09-09. It is the PUBLIC deposit's root README -- the first page
    # a reader of `immune-specificity-index` sees -- and it states all NINE
    # headline values (both ISIs and both CI bounds, both site AUROCs, both
    # cohort sizes, the MDE). Nothing checked one of them. It was found by
    # asking the ninth-axis question of `DOCS`: this list records a reason for
    # every document ADDED, which makes it look complete, while eighteen prose
    # documents sat outside it with no decision recorded either way.
    "snapshot/README.md",
    # Added 2026-09-17 (F3.3): every prose document that states a result now has
    # a recorded decision. Scanned with no mismatches on first addition:
    "17-CPTAC-PLAN.md", "08-READINESS-ASSESSMENT.md", "10-PROJECT-REFRAME.md",
    "pipeline/ENVIRONMENT.md", "pipeline/results/supplementary/CAPTIONS.md",
    "pipeline/results/supplementary/MANIFEST.md",
    # Scanned after a live stale number was fixed in each (two quoted a 78-test
    # suite and now state no count; the third quoted a prior paper's 28 cancer
    # types in the grammar the cancer-type pattern reads, now reworded):
    "12-MENTOR-OUTREACH.md", "submission/PREPRINT-PLAN.md",
    "13-NEXTGEN-EXTENDED-ABSTRACT.md",
    # Scanned so that EXEMPT above is load-bearing for them (see its note).
    "00-PROJECT-BRIEF.md", "04-CODE-AUDIT.md",
]


def _num(s: str) -> float:
    # `&minus;` is how poster.html writes a negative number; session 54 added it
    # so a poster pattern can capture the sign instead of dropping it.
    return float(s.replace(",", "").replace("&nbsp;", "").replace("−", "-")
                 .replace("&minus;", "-"))


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
    # when the first e-notation claims (S1's NSCLC p and q) were introduced.
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


def _wrap_line(lines: list[str], first: int, flat: str, pos: int) -> int:
    """The 1-based line on which a WRAP match at `flat[pos]` really sits.

    A WRAP pattern matches the paragraph flattened to one line, so the only
    line number to hand is the paragraph's FIRST. Until session 54 every WRAP
    claim was recorded there, which put MISMATCH reports on the wrong line and
    made `--coverage` list a guarded literal as unchecked on its real line --
    five poster numbers read as unguarded that a mutation sweep proved were
    caught. `flat` is the paragraph's lines joined by single spaces, so the
    literal's line follows from counting whitespace-separated tokens.
    """
    head = flat[:pos]
    token = len(head.split()) - (1 if head and not head[-1].isspace() else 0)
    for j in range(first, len(lines) + 1):
        n = len(lines[j - 1].split())
        if token < n:
            return j
        token -= n
    return first


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
         match_files: dict[str, set[str]] | None = None,
         only: set[str] | None = None,
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
        if only is not None and rel not in only:
            continue
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
            # Recorded UNSIGNED, because `_numeric_literals` never includes a
            # sign: a claim captured as "−0.023" used to be booked under a
            # string no literal had, and `--coverage` listed that guarded
            # number as unchecked. Session 54.
            covered.add((i, re.sub(r"^(?:[−+-]|&minus;)+", "", literal.strip())))
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
                        compare(table[key], m.group(1),
                                _wrap_line(lines, i, flat, m.start(1)), desc, tol)
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

        # Bracketed triples, matched over each PARAGRAPH flattened to one line
        # and attributed back to the line the point estimate sits on. Until
        # session 54 this loop ran per line, so a triple the typesetting had
        # wrapped -- "+0.3177 [0.2822," ending one line and "0.3533]" opening
        # the next -- was neither checked nor reported: invisible.
        seen_para: set[int] = set()
        found: list[tuple[int, re.Match]] = []
        for i0 in range(1, len(lines) + 1):
            text = paras.get(i0)
            if not text or id(text) in seen_para:
                continue
            seen_para.add(id(text))
            flat = " ".join(text.split())
            for m in list(TRIPLE.finditer(flat)) + list(
                    PAREN_TRIPLE.finditer(flat)):
                found.append((_wrap_line(lines, i0, flat, m.start(1)), m))
        for i, m in found:
            if i in ignored:
                continue
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
                # A superseded or diagnostic run's own stored triple is
                # provenance, not drift -- the same whitelist the
                # full-precision rule uses -- but only when all three numbers
                # match it at the precision written.
                lo_v, hi_v = _num(lo_lit), _num(hi_lit)
                prov = next((name for name, sv, slo, shi in _stored_triples(every_run)
                             if abs(pt - sv) <= _tol_for(pt_lit, 5e-5)
                             and abs(lo_v - slo) <= _tol_for(lo_lit, 5e-5)
                             and abs(hi_v - shi) <= _tol_for(hi_lit, 5e-5)),
                            None)
                literal = tuple(x.lstrip("+") for x in (pt_lit, lo_lit, hi_lit))
                if prov is not None or literal in UNSTORED_TRIPLES:
                    checked += 1
                    covered.update({(i, _x.lstrip("+−-")) for _x in literal})
                    if verbose:
                        say(f"  ok       {rel}:{i}  triple {m.group(0)} = "
                            f"{prov or UNSTORED_TRIPLES[literal]}")
                    continue
                # A POINT ESTIMATE THAT IDENTIFIES NOTHING IS A FAILURE.
                # It used to be reported and passed, on the grounds that
                # HPC4 diagnostics and ssGSEA arms were absent from the
                # table. They have since been added, and session 54 found
                # what the exemption was costing: the point is what
                # identifies a triple, so a DRIFTED point matched nothing,
                # was merely reported, and exited 0 -- the bounds were
                # guarded and the headline number beside them was not. A
                # mutation sweep moved seven Results points (the ssGSEA,
                # mean-z, PLAGE, reconstruction and global-axis arms) and
                # --strict stayed green for every one. A new triple now
                # needs a TRIPLES row before it can be quoted.
                unmatched.append(f"{rel}:{i}  {m.group(0)}")
                failures += 1
                say(f"  MISMATCH {rel}:{i}")
                say(f"    bracketed triple {m.group(0)!r}: its point "
                    "estimate matches no quantity in TRIPLES -- either it "
                    "drifted, or the quantity needs a TRIPLES row")

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
            # Session 57. A list comma after a number is part of the token
            # `_numeric_literals` returns ("2,611," in "2,611, 2,674 and"), and
            # a pattern's capture may or may not swallow it -- `[\d,]+` does,
            # a capture followed by a literal ", " does not. Compared raw, a
            # GUARDED number was listed as unguarded whenever the two
            # disagreed. Measured before fixing: changing the tokenizer instead
            # would have traded those false entries for new ones wherever a
            # greedy capture had kept the comma, so both sides drop a trailing
            # comma here and neither the tokenizer nor any count moves.
            covered_bare = {(j, w.rstrip(",")) for (j, w) in covered}
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
                    if (i, lit) in covered or (i, lit.rstrip(",")) in covered_bare:
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
    ('E16 pan-cancer full-cohort median drifted',
     r"median stays at 0\.998 at the full", 'median stays at 0.898 at the full'),
    ('E16 pan-cancer half-cohort median drifted',
     r"patients,\s+0\.998 at half", 'patients, 0.898 at half'),
    ('E16 pan-cancer quarter median drifted',
     r"at half and 0\.994 at a quarter", 'at half and 0.894 at a quarter'),
    ('E16 NSCLC full-cohort median drifted',
     r"in NSCLC, 0\.989, 0\.980", 'in NSCLC, 0.889, 0.980'),
    ('E16 NSCLC half-cohort median drifted',
     r"in NSCLC, 0\.989, 0\.980", 'in NSCLC, 0.989, 0.880'),
    ('E16 NSCLC quarter range low drifted',
     r"and 0\.894 to 0\.972 over", 'and 0.794 to 0.972 over'),
    ('E16 NSCLC quarter range high drifted',
     r"and 0\.894 to 0\.972 over", 'and 0.894 to 0.992 over'),
    ('E16 NSCLC small-threshold site count drifted',
     r"evaluates 47 NSCLC sites", 'evaluates 41 NSCLC sites'),
    ('E16 NSCLC small-threshold median drifted',
     r"NSCLC sites \(median 0\.987", 'NSCLC sites (median 0.887'),
    ('E16 NSCLC small-threshold percent drifted',
     r"\(median 0\.987, 94%", '(median 0.987, 84%'),
    ('E16 pan-cancer small-threshold site count drifted',
     r"and 320 pan-cancer sites", 'and 300 pan-cancer sites'),
    ('E16 pan-cancer small-threshold median drifted',
     r"pan-cancer sites \(median 0\.998", 'pan-cancer sites (median 0.898'),
    ('E16 pan-cancer small-threshold percent drifted',
     r"\(median 0\.998, 97% above", '(median 0.998, 87% above'),
    ("E6 PLAGE uncorrected CI drifted",
     r"\+0\.1353 \[0\.0889, 0\.1834\]", "+0.1353 [0.0989, 0.1834]"),
    ("E6 mean-z control CI drifted",
     r"\+0\.2963 \[0\.2442, 0\.3514\]", "+0.2963 [0.2542, 0.3514]"),
    ("E6 PLAGE reconstruction drifted", r"\(\+0\.0772, 11 of 16", "(+0.0972, 11 of 16"),
    ("E6 PLAGE reconstruction count drifted", r"\(\+0\.0772, 11 of 16", "(+0.0772, 14 of 16"),
    ("E8 omega curated median drifted", r"median ω 0\.977 against α", "median ω 0.877 against α"),
    ("E8 omega null median drifted", r"\(median ω 0\.227 against", "(median ω 0.327 against"),
    ("E8 alpha gap drifted", r"from 0\.126 to 0\.733", "from 0.226 to 0.733"),
    ("E8 omega gap drifted", r"from 0\.126 to 0\.733", "from 0.126 to 0.633"),
    ("E8 omega reconstruction CI drifted",
     r"0\.2147 \[0\.1739, 0\.2554\]", "0.2147 [0.1839, 0.2554]"),
    ("E8 alpha reconstruction CI drifted",
     r"0\.2917 \[0\.2584, 0\.3251\]", "0.2917 [0.2684, 0.3251]"),
    ("E14 pan-cancer adjusted median drifted", r"to an adjusted 42\.4%", "to an adjusted 44.4%"),
    ("E14 NSCLC adjusted median drifted", r"pan-cancer and 11\.2%", "pan-cancer and 14.2%"),
    ("E14 pan-cancer type-only median drifted", r"a median 38\.5% of label", "a median 48.5% of label"),
    ("E14 NSCLC type-only median drifted", r"and 3\.6% across the two", "and 13.6% across the two"),
    ("E14 NSCLC site-given-type median drifted", r"median\s+9\.9% in NSCLC", "median 5.9% in NSCLC"),
    ("E14 pan-cancer site-given-type beside NSCLC drifted", r"in NSCLC and 5\.5% pan-cancer",
     "in NSCLC and 9.5% pan-cancer"),
    ('E1 NSCLC global-axis CI drifted',
     r'−0\.0096 \[−0\.0630, 0\.0466\]', '−0.0096 [−0.0630, 0.0766]'),
    ('E1 pan-cancer global-axis CI drifted',
     r'0\.3512 \[0\.2692, 0\.4353\]', '0.3512 [0.2692, 0.4053]'),
    ('E3 size-only null CI drifted',
     r'0\.3041 \[0\.2434, 0\.3536\]', '0.3041 [0.2534, 0.3536]'),
    ('E5 3-fold CI drifted',
     r'0\.2188 \[0\.1673, 0\.2722\]', '0.2188 [0.1973, 0.2722]'),
    ('E5 10-fold CI drifted',
     r'0\.3154 \[0\.2607, 0\.3711\]', '0.3154 [0.2607, 0.3911]'),
    ('A16 NSCLC tumor-only CI drifted',
     r'0\.3448 \[0\.2731, 0\.4079\]', '0.3448 [0.2731, 0.4379]'),
    ('A16 pan-cancer tumor-only CI drifted',
     r'0\.2937 \[0\.2701, 0\.3149\]', '0.2937 [0.2801, 0.3149]'),
    ('E1 NSCLC signatures above null drifted',
     r'0\.0466\], with\s+1 of 16', '0.0466], with 4 of 16'),
    ('E1 NSCLC random-set r drifted',
     r'\(median r 0\.361 against', '(median r 0.261 against'),
    ('E1 NSCLC curated r drifted',
     r'against 0\.351\)', 'against 0.451)'),
    ('E1 NSCLC degenerate count drifted',
     r'and 10 of the 16', 'and 6 of the 16'),
    ('E1 pan-cancer signatures above null drifted',
     r'0\.3189\], with 15 of 16', '0.3189], with 12 of 16'),
    ('E1 pan-cancer curated r drifted',
     r'\(median r 0\.725 and', '(median r 0.625 and'),
    ('E1 pan-cancer random-set r drifted',
     r'and 0\.452, against', 'and 0.352, against'),
    ('pan-cancer registered curated r drifted',
     r'against 0\.391 and', 'against 0.491 and'),
    ('pan-cancer registered random-set r drifted',
     r'and 0\.100\)', 'and 0.200)'),
    ('E3 shift drifted',
     r'NSCLC index by 0\.012', 'NSCLC index by 0.021'),
    ('E4 400-draw CI low drifted',
     r'nearly so: \[0\.2336,', 'nearly so: [0.2436,'),
    ('E4 400-draw CI high drifted',
     r'0\.2336, 0\.3512\]', '0.2336, 0.3612]'),
    ('E4 1000-draw CI low drifted',
     r'and \[0\.2333, 0\.3511\]\.', 'and [0.2433, 0.3511].'),
    ('E4 1000-draw CI high drifted',
     r'and \[0\.2333, 0\.3511\]\.', 'and [0.2333, 0.3611].'),
    ('E5 3-fold signatures above null drifted',
     r'0\.2722\],\s+with 15 of 16', '0.2722], with 11 of 16'),
    ('A16 pan-cancer patients drifted',
     r'\(7,128 patients\)', '(7,168 patients)'),
    ("E11 composition probability drifted",
     r"probability 0\.0014, and", "probability 0.0140, and"),
    ("E11 extreme-relabelling count drifted",
     r"these six have \(1 of 8,008\)", "these six have (3 of 8,008)"),
    ("E11 relabelling total drifted",
     r"these six have \(1 of 8,008\)", "these six have (1 of 4,368)"),
    # B9: the 24-partition sweeps, one row per pattern.
    ("B9 NSCLC 24-partition sd drifted",
     r"24 partitions is 0\.0077, with a", "24 partitions is 0.0097, with a"),
    ("B9 NSCLC 24-partition sd lower bound drifted",
     r"interval of \[0\.0060, 0\.0107\]", "interval of [0.0040, 0.0107]"),
    ("B9 NSCLC 24-partition sd upper bound drifted",
     r"interval of \[0\.0060, 0\.0107\]", "interval of [0.0060, 0.0197]"),
    ("B9 NSCLC understatement drifted",
     r"uncertainty by 3\.3%, and partition", "uncertainty by 4.3%, and partition"),
    ("B9 NSCLC share drifted",
     r"makes up 6\.4% \[3\.9%", "makes up 8.4% [3.9%"),
    ("B9 NSCLC share lower bound drifted",
     r"6\.4% \[3\.9%, 11\.8%\]", "6.4% [2.9%, 11.8%]"),
    ("B9 NSCLC share upper bound drifted",
     r"6\.4% \[3\.9%, 11\.8%\]", "6.4% [3.9%, 35.5%]"),
    ("B9 pan-cancer 24-partition sd drifted",
     r"24 partitions is 0\.0072 \[", "24 partitions is 0.0065 ["),
    ("B9 pan-cancer sd lower bound drifted",
     r"0\.0072 \[0\.0056, 0\.0101\]", "0.0072 [0.0039, 0.0101]"),
    ("B9 pan-cancer sd upper bound drifted",
     r"0\.0072 \[0\.0056, 0\.0101\]", "0.0072 [0.0056, 0.0187]"),
    ("B9 pan-cancer understatement drifted",
     r"uncertainty by 19\.4%", "uncertainty by 16.0%"),
    ("B9 pan-cancer share drifted",
     r"makes up 29\.9% \[", "makes up 25.7% ["),
    ("B9 pan-cancer share lower bound drifted",
     r"29\.9% \[20\.5%, 45\.6%\]", "29.9% [11.0%, 45.6%]"),
    ("B9 pan-cancer share upper bound drifted",
     r"29\.9% \[20\.5%, 45\.6%\]", "29.9% [20.5%, 74.0%]"),
    ("B9 bootstrap ratio drifted",
     r"by 2\.652 against the", "by 2.649 against the"),
    ("B9 sqrt(patients) prediction drifted",
     r"against the 2\.756 predicted", "against the 2.556 predicted"),
    ("B9 partition ratio drifted",
     r"cohorts is 1\.060, with", "cohorts is 1.363, with"),
    ("B9 partition ratio lower bound drifted",
     r"interval of \[0\.697, 1\.611\]", "interval of [0.797, 1.611]"),
    ("B9 partition ratio upper bound drifted",
     r"interval of \[0\.697, 1\.611\]", "interval of [0.697, 2.611]"),
    ("B9 poster pan-cancer sd drifted", r"sds agree \(0\.0072,", "sds agree (0.0065,"),
    ("B9 poster NSCLC sd drifted", r", 0\.0077\) while", ", 0.0089) while"),
    ("B9 poster pan-cancer share drifted", r"separate, 29\.9% vs", "separate, 25.7% vs"),
    ("B9 poster NSCLC share drifted", r"% vs 6\.4%, as a floor", "% vs 8.4%, as a floor"),
    ("B9 sqrt(sites) prediction drifted",
     r"prediction \(3\.017\)", "prediction (3.217)"),
    # E20: the Linux ssGSEA tie-order figures quoted beside macOS's.
    ("A10 Linux max score difference drifted",
     r"\(macOS; 0\.0130 and", "(macOS; 0.0190 and"),
    ("A10 Linux score sd drifted",
     r"0\.0211 on Linux", "0.0291 on Linux"),
    ("A10 Linux fraction of an sd drifted",
     r"scores \(0\.62 on Linux\)", "scores (0.72 on Linux)"),
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
    # 0.0211 until 2026-09-24, when the note's misrounding was corrected (H92).
    ("A10 ssGSEA score sd drifted",
     r"per-signature score standard deviation of 0\.0212",
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
    # S1's derived NSCLC multiplicity values. Both are in SCIENTIFIC NOTATION,
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
    # The pan-cancer counterparts, added 2026-09-16 when HPC4 job 129159 closed
    # the project's longest-blocked item. The two pairs share their REGEXES --
    # "ISI of" and "moves by" -- and are told apart only by the paragraph each
    # pattern is scoped to, so these rows are also what proves that scoping
    # works: if the WRAP: anchors ever collapse back onto one paragraph, one of
    # these four injections stops being caught.
    # E10 (2026-09-17).
    ("NSCLC power at 0.03", r"the NSCLC arm's power is 45%", "the NSCLC arm's power is 65%"),
    ("NSCLC power at 0.04", r"68% at 0\.04", "78% at 0.04"),
    # E13/D4 (2026-09-16).
    ("cancer type alone on pooled PFI", r"alone reaches C = 0\.676", "alone reaches C = 0.776"),
    ("pan-cancer patients without a type", r"Of these, 19 patients", "Of these, 18 patients"),
    ("pan-cancer complete cases for the index", r"use the 7,149 complete", "use the 7,150 complete"),
    ("pan-cancer stable-sort sensitivity ISI",
     r"ISI of 0\.2968", "ISI of 0.3968"),
    ("pan-cancer stable-sort sensitivity shift",
     r"moves by 0\.0057", "moves by 0.0257"),
    ("signatures beating null (16/16 -> 15/16)", r"\b16/16\b", "15/16"),
    ("outcome arm (10/32 -> 11/32)", r"\b10/32\b", "11/32"),
    ("pan-cancer cohort size", r"7,168 TCGA patients", "7,169 TCGA patients"),
    ("number of cancer types", r"31 cancer types", "32 cancer types"),
    ("site AUROC", r"median AUROC 0\.998", "median AUROC 0.888"),
    ("minimum detectable effect", r"minimum detectable effect is 0\.046",
     "minimum detectable effect is 0.055"),
    ("rotation-null p, replaced wholesale", r"p = 0\.000999", "p = 0.5"),
    ("rotation-null p, one unit in its last printed place (the 2e-6 floor let this pass)",
     r"p = 0\.000999", "p = 0.000998"),
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
    # ---- added 2026-09-08 ---------------------------------------------------
    # The injected value must stay MATCHABLE by the pattern, or the row proves
    # only that the text exists. Both patterns anchor on prose, so any digits
    # survive the substitution and the COMPARISON is what fails.
    ("smallest pan-cancer excess_lo, the ISI's multiplicity margin",
     r"`excess_lo` is \+0\.14, nowhere", "`excess_lo` is +0.94, nowhere"),
    ("hypoxia outcome-arm excess_lo, axis-residualised",
     r"`excess_lo = \+0\.0069`", "`excess_lo = +0.0099`"),
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
     r"alone, gave 4\.1%", "alone, gave 9.1%"),
    ("E14 pan-cancer calibrated median drifted", r"leaves a median of 4\.7%", "leaves a median of 8.7%"),
    ("E14 pan-cancer calibrated minimum drifted", r"\(range 2\.8% to", "(range 1.8% to"),
    ("E14 pan-cancer calibrated maximum drifted", r"to 6\.0%; 6 of 16", "to 9.0%; 6 of 16"),
    ("E14 pan-cancer significant count drifted", r"6\.0%; 6 of 16", "6.0%; 9 of 16"),
    ("E14 NSCLC calibrated median drifted", r"median of 7\.0% across", "median of 4.0% across"),
    ("E14 NSCLC significant count drifted", r"in NSCLC \(15 of 16", "in NSCLC (12 of 16"),
    ("E14 poster calibrated median drifted", r"a median 0\.047 after", "a median 0.041 after"),
    # The immune-versus-other contrast, all three places it is quoted.
    ("immune-vs-other Mann-Whitney p (Results)",
     r"Mann–Whitney p = 0\.031", "Mann–Whitney p = 0.046"),
    ("immune-vs-other Fisher p (Results)",
     r"Fisher p = 0\.093", "Fisher p = 0.193"),
    ("immune-vs-other p (abstract)", r"\*P\* = 0\.031", "*P* = 0.046"),
    ("immune-vs-other p, IL6 counterfactual (Figure 3 caption)",
     r"\(p = 0\.031 → 0\.016\)", "(p = 0.031 → 0.048)"),
    # Number-WORDS. The exact defect: the reference list dropped to twelve and
    # the prose still said thirteen. Every other pattern in this file matches
    # digits and could not see it.
    ("reference count stated as a word",
     r"\bthirteen references,\s+all cited", "twelve references, all cited"),
    ("reference count stated as a word (verification sentence)",
     r"all thirteen verified at source", "all twelve verified at source"),
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
    # "; the real" -- the manuscript's outcome-arm Results sentence writes
    # "0.046 in C-index, while"
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

    # ---- session 50 (F1.1): Limitation 8's yardsticks and the fold census ----
    # Two of these rows inject the exact values the manuscript carried before
    # session 50 ("2.4" and "3.4"), which is the direction the real defect took.
    ("A9 Linux moved-patient count drifted",
     r"change fold on macOS and 180 on Linux", "change fold on macOS and 280 on Linux"),
    ("A9 platform gap in six-partition sds reverted to the old 2.4",
     r"0\.0218 is 2\.5 times", "0.0218 is 2.4 times"),
    ("A9 platform gap in 24-partition sds drifted",
     r"\(2\.9 times the 24-partition", "(2.4 times the 24-partition"),
    ("NSCLC sort-fix move in six-partition sds drifted",
     r"0\.0260, which is 2\.9", "0.0260, which is 2.4"),
    ("NSCLC sort-fix move in 24-partition sds drifted",
     r"\(3\.4 times the\s+24-partition value measured",
     "(3.9 times the 24-partition value measured"),
    ("NSCLC partition pairs at least as far apart as the sort-fix move drifted",
     r"24 partitions, 3 differ by as much", "24 partitions, 13 differ by as much"),
    ("NSCLC number of partition pairs drifted",
     r"of the 276 pairs among the 24 partitions, 3",
     "of the 267 pairs among the 24 partitions, 3"),
    ("pan-cancer sort-fix move in five-partition sds drifted",
     r"\*\*0\.87 times the pan-cancer", "**0.97 times the pan-cancer"),
    ("pan-cancer sort-fix move in 24-partition sds drifted",
     r"\(0\.79 times the 24-partition value\) — that is",
     "(0.89 times the 24-partition value) — that is"),
    ("pan-cancer partition pairs farther apart than the sort-fix move drifted",
     r"since 153 of the 276", "since 135 of the 276"),
    ("pan-cancer number of partition pairs drifted",
     r"since 153 of the 276 pairs", "since 153 of the 277 pairs"),
    ("NSCLC sort-fix move in the contrast sentence drifted",
     r"moved the estimate by\s+2\.9 partition", "moved the estimate by 2.4 partition"),
    ("pan-cancer sort-fix move in the contrast sentence drifted",
     r"here by 0\.9, a factor", "here by 0.7, a factor"),
    ("sort-fix contrast reverted to the old 3.4",
     r"a factor of 3\.3 \(", "a factor of 3.4 ("),
    ("sort-fix contrast in 24-partition units drifted",
     r"\(4\.3 in\s+24-partition units\)", "(4.9 in 24-partition units)"),
    ("pan-cancer tied sites drifted", r"600 of 619 holding", "601 of 619 holding"),
    ("pan-cancer tied-site patient share drifted",
     r"holding 79\.5% of its", "holding 75.9% of its"),
    ("NSCLC tied sites in the census sentence drifted",
     r"against 51 of 68\s+holding", "against 52 of 68 holding"),
    ("NSCLC tied-site patient share drifted",
     r"holding 38\.1% in NSCLC", "holding 31.8% in NSCLC"),
    ("pan-cancer patients moved by the sort fix drifted",
     r"moved 3,812 of 7,168", "moved 3,821 of 7,168"),
    ("pan-cancer share moved by the sort fix drifted",
     r"patients \(53\.2%\)", "patients (35.2%)"),
    ("NSCLC patients moved by the sort fix drifted",
     r"against 216 of 944 \(22\.9%\)", "against 261 of 944 (22.9%)"),
    ("NSCLC share moved by the sort fix drifted",
     r"of 944 \(22\.9%\) in NSCLC", "of 944 (29.2%) in NSCLC"),
    ("pan-cancer adjusted Rand index drifted",
     r"is 0\.13\s+pan-cancer and", "is 0.31 pan-cancer and"),
    ("NSCLC adjusted Rand index drifted",
     r"pan-cancer and 0\.54 in NSCLC", "pan-cancer and 0.45 in NSCLC"),
    ("Linux full-cohort NSCLC site AUROC drifted",
     r"median reads\s+0\.989", "median reads 0.899"),
    ("macOS NSCLC site AUROC beside the Linux one drifted (Limitation 8)",
     r"where the macOS run gives 0\.992", "where the macOS run gives 0.929"),
    ("macOS NSCLC site AUROC beside the Linux one drifted (Negative controls)",
     r"here and 0\.992 above", "here and 0.929 above"),
    # The exact regression session 50 found in the cover letter (F3.8): the
    # IL6-counted p beside the IL6-excluded counts, across a line wrap.
    # F3.8: the wrapped occurrences, injected in their raw (wrapped) form.
    ("cover letter pan-TCGA cohort size across a wrap drifted",
     r"Across 7,168 TCGA\npatients", "Across 7,169 TCGA\npatients"),
    # Session 55: one row per pattern added for the cover letter's science.
    # Anchors are letter-specific wherever the sentence is the letter's own;
    # the two reliability rows deliberately are NOT, because that sentence is
    # shared with the conference abstract, the NextGen abstract and the
    # manuscript, and injecting into all of them is the honest proof that the
    # pattern fires everywhere it reaches rather than only where it was aimed.
    ("cover letter NSCLC cohort size drifted",
     r"in a separate 944-patient NSCLC", "in a separate 954-patient NSCLC"),
    ("cover letter family-level rotation p drifted",
     r"correlation gives \*P\* = 0\.001",
     "correlation gives *P* = 0.010"),
    ("cover letter rotation null draws drifted",
     r"1,000 draws\. The image reads", "1,100 draws. The image reads"),
    ("cover letter random-set alpha before residualisation drifted",
     r"falls from 0\.97 to 0\.80 while the median curated",
     "falls from 0.87 to 0.80 while the median curated"),
    ("cover letter random-set alpha after residualisation drifted",
     r"falls from 0\.97 to 0\.80 while the median curated",
     "falls from 0.97 to 0.90 while the median curated"),
    ("cover letter NSCLC reliability-gap fold drifted",
     r"gap is then 5\.3-fold", "gap is then 5.9-fold"),
    ("cover letter pan-cancer reliability-gap fold drifted",
     r"5\.3-fold \(NSCLC\) to 13-fold", "5.3-fold (NSCLC) to 15-fold"),
    ("cover letter small-panel gap at 160 genes drifted",
     r"scales inversely with panel size — 0\.084",
     "scales inversely with panel size — 0.094"),
    ("cover letter small-panel gap at 10 genes drifted",
     r"rising to 0\.453 at 10", "rising to 0.493 at 10"),
    ("cover letter frozen NSCLC index, the reproducibility promise, drifted",
     r"the 0\.3182 frozen on macOS", "the 0.3172 frozen on macOS"),
    ("abstract draft NSCLC cohort size across a wrap drifted",
     r"TCGA NSCLC\n\(n=944\) as a second", "TCGA NSCLC\n(n=945) as a second"),
    ("NextGen NSCLC cohort size across a wrap drifted",
     r"with 944\nnon-small cell lung", "with 945\nnon-small cell lung"),
    ("audit A9 tied-site count across a wrap drifted",
     r"\*\*51 of the 68 NSCLC sites sit in a\ntied", "**52 of the 68 NSCLC sites sit in a\ntied"),
    ("checklist rotation p across a wrap drifted",
     r"rotation null\np=0\.000999", "rotation null\np=0.000899"),
    ("checklist rotation p across a wrap, one unit in its last printed place",
     r"rotation null\np=0\.000999", "rotation null\np=0.000998"),
    ("audit A6 small-panel Spearman heading drifted",
     r"Spearman\(k, gap\) = −1\.000", "Spearman(k, gap) = −0.100"),
    # F7.4: one row per caption authority, each on CAPTIONS.md's own phrasing.
    ("S7 caption NSCLC site AUROC drifted",
     r"median AUROC 0\.998 pan-TCGA, 0\.992 NSCLC", "median AUROC 0.998 pan-TCGA, 0.929 NSCLC"),
    ("S7 caption pan-TCGA evaluable sites drifted", r"rows are 98 pan-TCGA", "rows are 89 pan-TCGA"),
    ("S7 caption NSCLC evaluable sites drifted", r"plus 15 NSCLC sites twice", "plus 51 NSCLC sites twice"),
    ("S7 caption ComBat median drifted", r"below chance \(median 0\.005\)", "below chance (median 0.050)"),
    ("S4/S5 caption pan-TCGA site count drifted", r"enter \(202 of 619", "enter (220 of 619"),
    ("S4/S5 caption NSCLC site count drifted", r"31 of 68 NSCLC\)", "13 of 68 NSCLC)"),
    ("S4 caption NSCLC median R^2 drifted", r"R-squared of 0\.145, against", "R-squared of 0.154, against"),
    ("S4 caption pan-TCGA median R^2 drifted", r"against 0\.444 here", "against 0.404 here"),
    ("S5 caption NSCLC median R^2 drifted", r"R-squared is 0\.145 here", "R-squared is 0.415 here"),
    ("S5 caption pan-TCGA median R^2 drifted", r"against 0\.444 there", "against 0.404 there"),
    ("S5 caption NSCLC patients drifted", r"contributes 796 patients", "contributes 769 patients"),
    ("S5 caption NSCLC sites drifted", r"across 31 usable sites", "across 13 usable sites"),
    ("S5 caption pan-TCGA patients drifted", r"pan-TCGA's 5,775 across", "pan-TCGA's 5,757 across"),
    ("S5 caption pan-TCGA sites drifted", r"5,775 across 202,", "5,775 across 220,"),
    ("S9 caption EMT unadjusted gain drifted", r"transition: 0\.059 unadjusted", "transition: 0.095 unadjusted"),
    ("S9 caption EMT adjusted gain drifted", r"unadjusted, 0\.065 adjusted", "unadjusted, 0.056 adjusted"),
    ("S3 caption Spearman drifted", r"\(Spearman rho = -1\.00\)", "(Spearman rho = -0.10)"),
    ("cover letter's Mann-Whitney p reverted to the IL6-counted 0.016",
     r"Mann–Whitney\n\*P\* = 0\.031", "Mann–Whitney\n*P* = 0.016"),
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
    # The abstract's remaining headroom. Both operands were already checked and
    # the SUBTRACTION between them was not, so the document once read "2,515 /
    # 2,600 with 83 characters of headroom" -- each count right, the difference
    # wrong by two. The row injects that same two-character error. It is LIVE
    # since 2026-09-17: it hard-coded "85 characters", the counter fix of that
    # day made every document say 88, and the row matched nothing -- the
    # self-test failed on its own stale literal (full_suite_20260917_run50b.log).
    if "abstract_headroom" in table:
        h = int(table["abstract_headroom"])
        rows.append(("abstract headroom (limit minus live count)",
                     rf"\b{h} characters of headroom", f"{h - 2} characters of headroom"))
        # Session 51: the dry run's own form of the same number.
        rows.append(("abstract headroom, dry-run form",
                     rf"\),\s+headroom\s+{h},", f"), headroom {h - 2},"))
    # Session 51: each dry-run variant cell, pushed 7 off its live value.
    for key, prefix, suffix in ABSTRACT_VARIANT_SPEC:
        if key in table:
            rows.append((f"dry-run abstract variant {key}",
                         "(?m)(" + prefix.removeprefix("(?m)") + ")" + _NUM
                         + "(" + suffix + ")",
                         r"\g<1>" + str(int(table[key]) + 7) + r"\g<2>"))
    # Session 57: each prose re-quotation of those tables, pushed 7 off its
    # live value -- its own row, so a quotation cannot pass on the back of the
    # table cell it copies (session 53's both-copies rule).
    for i, (key, prefix, suffix) in enumerate(DRYRUN_PROSE_SPEC):
        if key in table:
            rows.append((f"dry-run prose re-quotes {key} (quotation {i + 1})",
                         "(?m)(" + prefix.removeprefix("(?m)") + ")" + _NUM
                         + "(" + suffix + ")",
                         r"\g<1>" + str(int(table[key]) + 7) + r"\g<2>"))
    # Session 57: the README's demo table, each cell pushed 0.010 off.
    for key, prefix, suffix in DEMO_TABLE_SPEC:
        if key in table:
            rows.append((f"README demo table drifted, {key}",
                         "(?m)(" + prefix + ")" + _DCELL + "(" + suffix + ")",
                         r"\g<1>" + f"{table[key] + 0.010:.3f}" + r"\g<2>"))
    # Session 57: the README's sort-fix ratio and partition sd, each pushed
    # four printed places off its live value.
    for key, prefix, suffix, fmt in README_SORTFIX_SPEC:
        if key in table:
            v = table[key]
            step = 4 * 10 ** -int(fmt[3])
            rows.append((f"README sort-fix paragraph drifted, {key} as {fmt}",
                         "(" + prefix + ")" + fmt.format(v).replace(".", r"\.") + "(" + suffix + ")",
                         r"\g<1>" + fmt.format(v + step) + r"\g<2>"))
    # Session 57: the README's stored out-of-fold ComBat median.
    if "cap_combat_auroc_median" in table:
        c = table["cap_combat_auroc_median"]
        rows.append(("README stored ComBat median drifted",
                     rf"(out-of-fold median is ){c:.3f}",
                     r"\g<1>" + f"{c + 0.004:.3f}"))
    # Session 57: A16's extent, each claim pushed 2 off its live value.
    for key, what, shape in A16_EXTENT_SPEC:
        if key in table:
            v = int(table[key])
            rows.append((f"A16 extent drifted, {what}",
                         shape.replace("{v}", f"{v:,}"),
                         r"\g<1>" + f"{v + 2:,}" + r"\g<2>"))
    # Session 58: the checklists' figures, each moved four units in its last
    # printed place.
    for key, what, shape, fmt in CHECKLIST_SPEC:
        if key in table:
            v = table[key]
            places = int(fmt.split(".")[1].rstrip("f}")) if "." in fmt else 0
            step = 4 * 10 ** -places
            rows.append((f"checklist figure drifted, {what}",
                         shape.replace("{v}", re.escape(fmt.format(v))),
                         r"\g<1>" + fmt.format(v + step) + r"\g<2>"))
    # Session 60: A5's NSCLC partition rows, each moved four units.
    for _s in range(6):
        if f"pv_nsclc_seed{_s}" in table:
            _v = f"{table[f'pv_nsclc_seed{_s}']:.4f}"
            rows.append((f"audit A5 partition row drifted, NSCLC seed {_s}",
                         rf"(?m)(^\| {_s} \| ){re.escape(_v)}( \|$)",
                         r"\g<1>" + f"{float(_v) + 0.0004:.4f}" + r"\g<2>"))
    # Session 60: the audit sweep's bindings, the same way.
    for key, what, shape, fmt in AUDIT_SWEEP_SPEC:
        if key in table:
            v = table[key]
            places = int(fmt.split(".")[1].rstrip("f}")) if "." in fmt else 0
            step = 4 * 10 ** -places
            rows.append((f"audit sweep figure drifted, {what}",
                         shape.replace("{v}", re.escape(fmt.format(v))),
                         r"\g<1>" + fmt.format(v + step) + r"\g<2>"))
    # Session 59: the audit's A16 paragraph, the same way.
    for key, what, shape, fmt in A16_NARRATIVE_SPEC:
        if key in table:
            v = table[key]
            places = int(fmt.split(".")[1].rstrip("f}")) if "." in fmt else 0
            step = 4 * 10 ** -places
            rows.append((f"audit A16 figure drifted, {what}",
                         shape.replace("{v}", re.escape(fmt.format(v))),
                         r"\g<1>" + fmt.format(v + step) + r"\g<2>"))
    # B-21: the manuscript's Results summaries, the same way.
    for key, what, shape, fmt in B21_SUMMARY_SPEC:
        if key in table:
            v = table[key]
            places = int(fmt.split(".")[1].rstrip("f}")) if "." in fmt else 0
            step = 4 * 10 ** -places
            rows.append((f"B-21 summary figure drifted, {what}",
                         shape.replace("{v}", re.escape(fmt.format(v))),
                         r"\g<1>" + fmt.format(v + step) + r"\g<2>"))
    # Session 58: the small counts, each pushed 2 off its ROUNDED live value
    # (three are derived ratios or percentages printed as integers).
    for key, what, shape in SMALL_COUNT_SPEC:
        if key in table:
            v = round(table[key])
            rows.append((f"small count drifted, {what}",
                         shape.replace("{v}", f"{v:,}"),
                         r"\g<1>" + f"{v + 2:,}" + r"\g<2>"))
    # Session 52: the cover letter's length, its paragraph-4 fallback and the
    # subtraction between them. The subtraction gets its own row because that
    # is the shape the abstract's headroom defect took -- both counts right,
    # the difference wrong -- and here the difference is what the author acts
    # on at submission.
    if "cover_letter_body_words" in table:
        w = int(table["cover_letter_body_words"])
        rows.append(("cover letter body length (live recount)",
                     rf"\*\*Body length:\s+{w} words\*\*",
                     f"**Body length: {w - 11} words**"))
    if "cover_letter_para4_words" in table:
        p = int(table["cover_letter_para4_words"])
        rows.append(("cover letter paragraph 4 length (live recount)",
                     rf"\*\*{p} words measured\*\*",
                     f"**{p + 6} words measured**"))
    if "cover_letter_body_minus_para4" in table:
        d = int(table["cover_letter_body_minus_para4"])
        rows.append(("cover letter body minus paragraph 4 (the subtraction)",
                     rf"(words measured\*\*\s+—\s+and\s+){d}(\s+remain)",
                     r"\g<1>" + str(d - 2) + r"\g<2>"))
    # Session 55: one row per member of the paragraph enumeration, plus one
    # for its sum and one for the scaffolding the sum leaves over. Each
    # anchors on its own POSITION in the list, so a row proves its own
    # pattern fires rather than borrowing a neighbour's -- the "both copies
    # need their own row" rule from session 53, applied within one sentence.
    # WHITESPACE IN THESE SHAPES IS `\s+`, since 2026-09-24. They used literal
    # spaces, and H83's re-wrap of the letter (session 59, `eeb95f7`) put line
    # breaks after "paragraphs are" and inside "in all": six of the seven rows
    # then found NO SUCH TEXT while the WRAP-scoped patterns they exist to test
    # kept passing `--strict`. Session 60's baseline suite caught it; nothing
    # between `eeb95f7` and that suite ran the self-test.
    _enum_rows = (
        ("cover_letter_para_submit_words", "submission paragraph",
         r"(remaining\s+paragraphs\s+are\s+){v}(\s+/\s+)"),
        ("cover_letter_para_finding_words", "finding paragraph",
         r"(remaining\s+paragraphs\s+are\s+\d+\s+/\s+){v}(\s+/\s+)"),
        ("cover_letter_para_immune_words", "immune-signal paragraph",
         r"(remaining\s+paragraphs\s+are\s+\d+\s+/\s+\d+\s+/\s+){v}(\s+/\s+)"),
        ("cover_letter_para_notwhat_words", "scope-limits paragraph",
         r"(remaining\s+paragraphs\s+are\s+\d+\s+/\s+\d+\s+/\s+\d+\s+/\s+){v}(\s+words\s+of)"),
        ("cover_letter_para_scope_words", "boilerplate paragraph",
         r"(words\s+of\s+argument\s+and\s+){v}(\s+of)"),
        ("cover_letter_enumerated_words", "enumerated-paragraph sum",
         r"(required\s+boilerplate\s+—\s+){v}(\s+in\s+all)"),
        ("cover_letter_scaffold_words", "letterhead scaffolding",
         r"(in\s+all,\s+the\s+other\s+){v}(\s+being\s+the\s+date)"),
    )
    for key, what, shape in _enum_rows:
        if key not in table:
            continue
        v = int(table[key])
        rows.append((f"cover letter {what} words (live recount)",
                     shape.replace("{v}", str(v)),
                     r"\g<1>" + str(v + 3) + r"\g<2>"))
    # Session 53: the manuscript abstract's two counts and their SUM. The sum
    # is injected in the direction the real defect takes -- upward, past the
    # 131-264 ceiling the sentence claims to sit inside -- so a MISSED here
    # would mean the ceiling claim had gone unguarded.
    # Session 53: the cover letter's two AACR-required counts of the manuscript.
    if "paper_body_words" in table:
        _bw = int(table["paper_body_words"])
        rows.append(("cover letter: manuscript body word count",
                     rf"(\*\*){_bw:,}(\s+words\*\* excluding the cover page)",
                     r"\g<1>" + f"{_bw + 250:,}" + r"\g<2>"))
    if "paper_figure_count" in table:
        _fc = int(table["paper_figure_count"])
        rows.append(("cover letter: manuscript figure count",
                     rf"(with \*\*){_fc}(\s+figures and no tables\*\*)",
                     r"\g<1>" + str(_fc + 1) + r"\g<2>"))
    if "paper_abstract_body_words" in table:
        b = int(table["paper_abstract_body_words"])
        rows.append(("manuscript abstract body words (live recount)",
                     rf"(measures \*\*){b}(\s+words of body)",
                     r"\g<1>" + str(b + 5) + r"\g<2>"))
    if "paper_abstract_significance_words" in table:
        s = int(table["paper_abstract_significance_words"])
        rows.append(("manuscript abstract Significance words (live recount)",
                     rf"(words of body and ){s}(\s+of Significance)",
                     r"\g<1>" + str(s - 3) + r"\g<2>"))
    if "paper_abstract_total_words" in table:
        t2 = int(table["paper_abstract_total_words"])
        rows.append(("manuscript abstract total (the sum, pushed past the ceiling)",
                     rf"(of Significance,\s+){t2}(\s+in total\*\*)",
                     r"\g<1>" + str(t2 + 4) + r"\g<2>"))
    # Session 53: the poster's five drawn outcome excesses and its multiplicity
    # margin. Each injection is +0.006 on the PRINTED (3-decimal) form, which is
    # an order of magnitude above the 5e-4 tolerance and still small enough that
    # a row passing here would mean the pattern had gone vacuous rather than
    # that the tolerance was generous.
    for _k, _pre, _post in (
        ("g2m", "G2M ", ", E2F"),
        ("e2f", "E2F ", ", angiogenesis"),
        ("angio", "angiogenesis ", ", EMT"),
        ("emt", "EMT ", ", hypoxia"),
        ("hypoxia", "hypoxia ", r" \(Fisher exact"),
    ):
        _key = f"pancancer_outcome_excess_{_k}"
        if _key in table:
            _shown = f"{table[_key]:.3f}"
            rows.append((f"poster outcome excess, {_k} (drawn on the board)",
                         rf"({_pre}){re.escape(_shown)}({_post})",
                         r"\g<1>" + f"{table[_key] + 0.006:.3f}" + r"\g<2>"))
    if "pancancer_hypoxia_excess_lo" in table:
        _lo = table["pancancer_hypoxia_excess_lo"]
        rows.append(("poster multiplicity margin, hypoxia excess_lo",
                     rf"(narrowest at excess_lo = \+){_lo:.4f}",
                     r"\g<1>" + f"{_lo + 0.0006:.4f}"))
    # Session 53: the deposited README's results table. Each row is injected in
    # the shape a transcription slip actually takes -- one digit of the reported
    # precision -- against the table as it is written, not against a copy.
    for _key, _pre, _post, _fmt, _d in (
        ("pancancer_n", r"Pan-TCGA \(n=", r", 31 types\)", "{:,.0f}", 100),
        ("nsclc_n", r"\| NSCLC \(n=", r"\) \|", "{:,.0f}", 10),
        ("pancancer_isi", r"\| ISI \| \*\*", r" \[", "{:.4f}", 4e-4),
        ("pancancer_ci_lo", r"\| ISI \| \*\*[\d.]+ \[", r",", "{:.4f}", 4e-4),
        ("nsclc_mde", r"0/32 \(MDE ", r"\)", "{:.3f}", 6e-3),
        ("pancancer_site_auroc", r"\| Site AUROC \| ", r" \|", "{:.3f}", -3e-3),
    ):
        if _key in table:
            _now = _fmt.format(table[_key])
            _new = _fmt.format(table[_key] + _d)
            rows.append((f"deposited README results table, {_key}",
                         rf"({_pre}){re.escape(_now)}({_post})",
                         r"\g<1>" + _new + r"\g<2>"))
    # The manuscript's copy of the same five needs its OWN rows. A row that
    # only rewrites the poster would still report CAUGHT -- both documents read
    # the same authority, so the poster's mismatch alone fails the run -- and
    # the manuscript pattern would never once be proven to fire. That is the
    # vacuity this file has shipped before, wearing a green row.
    for _k, _pre, _post in (
        ("g2m", r"G2M checkpoint\*\* \(excess ", r"\)"),
        ("e2f", r"E2F targets\*\*\n\(", r"\)"),
        ("angio", r"\*\*angiogenesis\*\* \(", r"\)"),
        ("emt", r"mesenchymal transition\*\* \(", r"\)"),
        ("hypoxia", r"\*\*hypoxia\*\* \(", r"\)"),
    ):
        _key = f"pancancer_outcome_excess_{_k}"
        if _key in table:
            _shown = f"{table[_key]:.3f}"
            rows.append((f"manuscript Results outcome excess, {_k}",
                         rf"({_pre}){re.escape(_shown)}({_post})",
                         r"\g<1>" + f"{table[_key] + 0.006:.3f}" + r"\g<2>"))
    # Session 54, the manuscript's Results: one row per claim in
    # `_MANUSCRIPT_S54`, built from that table so the claim and the row that
    # proves it can fail cannot drift apart. Spaces in the prefix and suffix
    # match any whitespace, because most of these claims wrap in the file.
    for (_key, _c, _t, _d, _s, _pre, _post, _fmt, _step, _absval) in _MANUSCRIPT_S54:
        if _pre is None or _key not in table:
            continue
        _v = abs(table[_key]) if _absval else table[_key]
        _pre_rx, _post_rx = _pre.replace(" ", r"\s+"), _post.replace(" ", r"\s+")
        rows.append((f"session 54 manuscript claim, {_key}: {_d}",
                     rf"({_pre_rx}){re.escape(_fmt.format(_v))}({_post_rx})",
                     r"\g<1>" + _fmt.format(_v + _step) + r"\g<2>"))
    # Session 54: one row per poster and landing-page claim added this session,
    # each injected against the file as written, in the shape a transcription
    # slip takes -- one step in the last printed place, well above the implied
    # tolerance. `absval` rows are the three printed with `&minus;`: the row
    # rewrites the magnitude and leaves the sign where the page put it.
    for _key, _pre, _post, _fmt, _d, _absval in (
        # The lead's ISI and the pan-cancer MDE are anchored on poster-only
        # markup: the abstract and the manuscript print the same words, and a
        # row rewriting all three would report CAUGHT on their patterns alone.
        ("pancancer_isi", r'<p class="lead">ISI = ', r" \(95% CI", "{:.3f}", 0.006, False),
        ("pancancer_ci_lo", r"\(95% CI ", r"&ndash;", "{:.3f}", 0.006, False),
        ("pancancer_ci_hi", r"\(95% CI [\d.]+&ndash;", r"\) pan-TCGA", "{:.3f}", 0.006, False),
        ("nsclc_isi", r"pan-TCGA and ", r" \(", "{:.3f}", 0.006, False),
        ("nsclc_ci_lo", r"pan-TCGA and [\d.]+ \(", r"&ndash;", "{:.3f}", 0.006, False),
        ("nsclc_ci_hi", r"&ndash;", r"\) in NSCLC\.", "{:.3f}", 0.006, False),
        ("pancancer_excess_min", r"Per-signature excess ", r"&ndash;", "{:.3f}", 0.006, False),
        ("pancancer_excess_max", r"Per-signature excess [\d.]+&ndash;", r" pan-cancer",
         "{:.3f}", 0.006, False),
        ("nsclc_excess_min", r"pan-cancer, ", r"&ndash;[\d.]+ NSCLC\. Family", "{:.3f}", 0.006, False),
        ("nsclc_excess_max", r"pan-cancer, [\d.]+&ndash;", r" NSCLC\. Family", "{:.3f}", 0.006, False),
        ("rotation_null_mean", r"null mean of ", r" \(95% range", "{:.3f}", 0.006, False),
        ("rotation_null_lo", r"\(95% range ", r"&ndash;", "{:.3f}", 0.006, False),
        ("rotation_null_hi", r"\(95% range [\d.]+&ndash;", r"\)", "{:.3f}", 0.006, False),
        ("pancancer_r_axis", r"axis r = ", r" and", "{:.3f}", 0.006, False),
        ("nsclc_r_axis", r"axis r = [\d.]+ and &minus;", r",", "{:.3f}", 0.006, True),
        ("nsclc_honest_ci_lo", r"widens the interval to \[", r",", "{:.4f}", 0.0006, False),
        ("nsclc_honest_ci_hi", r"widens the interval to \[[\d.]+, ", r"\]", "{:.4f}", 0.0006, False),
        ("nsclc_interval_understated_pct", r"\], ", r"% wider than reported", "{:.1f}", 0.6, False),
        ("pancancer_honest_ci_lo", r"widens it to \[", r",", "{:.4f}", 0.0006, False),
        ("pancancer_honest_ci_hi", r"widens it to \[[\d.]+, ", r"\]", "{:.4f}", 0.0006, False),
        ("pancancer_interval_understated_pct", r"<b>", r"% wider</b>", "{:.1f}", 0.6, False),
        ("pancancer_partition_var_pct", r"larger share pan-cancer \(", r"% vs", "{:.1f}", 0.6, False),
        ("nsclc_partition_var_pct", r"% vs ", r"%\)\. At 24", "{:.1f}", 0.6, False),
        ("pancancer_alpha_null_raw", r"Cronbach&rsquo;s &alpha; ", r" pan-cancer", "{:.3f}", 0.006, False),
        ("pancancer_alpha_null_resid", r"null&rsquo;s\s+&alpha; to <b>", r"</b>", "{:.3f}", 0.006, False),
        ("pancancer_alpha_obs_resid", r"hold at <b>", r"</b>", "{:.3f}", 0.006, False),
        ("pancancer_rel_gap_raw", r"gap is\s+\+", r" pan-cancer and", "{:.3f}", 0.006, False),
        ("nsclc_rel_gap_raw", r"pan-cancer and \+", r" NSCLC on raw", "{:.3f}", 0.006, False),
        ("pancancer_rel_gap_resid", r"against <b>\+", r"</b> and", "{:.3f}", 0.006, False),
        ("nsclc_rel_gap_resid", r"</b> and <b>\+", r"</b> on the", "{:.3f}", 0.006, False),
        ("small_panel_k10", r'<td class="hi">10</td><td class="hi">', r"</td>", "{:.4f}", 0.0006, False),
        ("small_panel_k20", r"<td>20</td><td>", r"</td>", "{:.4f}", 0.0006, False),
        ("small_panel_k40", r"<td>40</td><td>", r"</td>", "{:.4f}", 0.0006, False),
        ("small_panel_k80", r"<td>80</td><td>", r"</td>", "{:.4f}", 0.0006, False),
        ("small_panel_k160", r"<td>160</td><td>", r"</td>", "{:.4f}", 0.0006, False),
        ("small_panel_spearman", r"Spearman &rho; = &minus;", r", p &lt;", "{:.3f}", -0.006, True),
        ("small_panel_ratio_k10_k160", r"<b>", r"&times;</b> larger", "{:.1f}", 0.6, False),
        ("nsclc_n", r"NSCLC ", r" \(\d+ events\)", "{:,.0f}", 10, False),
        ("nsclc_n_events", r"NSCLC [\d,]+ \(", r" events\)", "{:,.0f}", 5, False),
        ("sortaudit_macos_genes", r"Xena TOIL expression, ", r" genes", "{:,.0f}", 100, False),
        ("nsclc_site_auroc", r"</b> pan-cancer,\s+", r" NSCLC\. Preserved", "{:.3f}", -0.006, False),
        ("pancancer_emb_over_cov", r"&Delta;r = &minus;", r" over site", "{:.3f}", 0.006, True),
        ("nsclc_emb_over_cov", r"stage; NSCLC \+", r"\) &mdash;", "{:.3f}", 0.006, False),
        ("label_site_variance_r2", r"site R&sup2; ", r" crude", "{:.3f}", 0.006, False),
        ("label_site_given_type_r2", r"crude, ", r" given cancer type", "{:.3f}", 0.006, False),
        ("nsclc_isi", r"\(NSCLC ISI ", r" vs", "{:.4f}", 0.0006, False),
        ("immune_mwu_p_2s", r'<div class="n">p = ', r'</div>\s+<div class="l">Mann',
         "{:.3f}", 0.006, False),
        ("immune_fisher_p_2s", r"\(Fisher exact p = ", r"\)", "{:.3f}", 0.006, False),
        ("pancancer_mde", r"about [\d.]+; pan-cancer MDE\s+", r"\)", "{:.3f}", 0.006, False),
        ("nsclc_site_auroc", r"\| Site AUROC \(median\) \| [\d.]+ \| ", r" \|", "{:.3f}", -0.006, False),
        ("a7_meanz_rel_null", r"random-set reliability from ", r" to", "{:.3f}", 0.006, False),
        ("a7_ssgsea_rel_null", r"random-set reliability from [\d.]+ to ", r"\.", "{:.3f}", 0.006, False),
    ):
        if _key in table:
            _v = abs(table[_key]) if _absval else table[_key]
            rows.append((f"session 54 poster/landing claim, {_key} after {_pre[:24]!r}",
                         rf"({_pre}){re.escape(_fmt.format(_v))}({_post})",
                         r"\g<1>" + _fmt.format(_v + _d) + r"\g<2>"))
    # Session 60 (H92): the 24-partition SEs the six-partition pattern used to
    # claim, and the macOS score sd whose 0.005 floor hid a misrounding -- each
    # moved one unit in its last printed place, which the old floors accepted.
    for _key, _pre, _post, _fmt, _d in (
        ("pv24_nsclc_bootstrap_se", r"bootstrap SE ", r", combined 0\.0303", "{:.4f}", 0.0001),
        ("pv24_nsclc_combined_se", r"bootstrap SE 0\.0294, combined ", r";", "{:.4f}", 0.0001),
        ("pv24_pancancer_bootstrap_se", r"Bootstrap SE ", r" and combined", "{:.6f}", 0.000001),
        ("pv24_pancancer_combined_se", r"Bootstrap SE 0\.011073 and combined ", r",", "{:.6f}", 0.000001),
        ("sortaudit_macos_score_sd", r"standard deviation of ", r" \(macOS;", "{:.4f}", 0.0001),
        ("sortaudit_macos_max_diff", r"orderings is ", r", against a mean", "{:.4f}", 0.0001),
    ):
        if _key in table:
            rows.append((f"session 60 printed-precision claim, {_key}",
                         rf"({_pre}){re.escape(_fmt.format(table[_key]))}({_post})",
                         r"\g<1>" + _fmt.format(table[_key] + _d) + r"\g<2>"))
    # Session 60 (H86): the pinned-split paragraph, one row per figure, each moved
    # in its last printed place by more than its tolerance.
    for _key, _pre, _post, _fmt, _d in (
        ("pinned_nsclc_median_dr", r"NSCLC gives a median Δr\s+of\s+\+", r",", "{:.3f}", 0.004),
        ("pinned_nsclc_dmae", r"Δ-MAE\s+", r"\s+\[0\.0014", "{:.4f}", 0.0004),
        ("pinned_nsclc_emb_over_cov", r"over\s+covariates\s+of\s+\+", r"\s+and", "{:.3f}", 0.004),
        ("pinned_nsclc_site_auroc", r"median\s+site\s+AUROC\s+of\s+", r";\s+pan-cancer gives",
         "{:.3f}", -0.004),
        ("pinned_pancancer_median_dr", r"pan-cancer gives a median Δr\s+of\s+\+", r",", "{:.3f}", 0.004),
        ("pinned_pancancer_dmae", r"Δ-MAE\s+", r"\s+\[0\.0083", "{:.4f}", 0.0004),
        ("pinned_pancancer_emb_over_cov_abs", r"embedding\s+advantage\s+of\s+[−-]", r"\s+and\s+a",
         "{:.3f}", 0.004),
        ("pinned_pancancer_site_auroc", r"median\s+site\s+AUROC\s+of\s+", r"\.\s+No conclusion",
         "{:.3f}", -0.004),
    ):
        if _key in table:
            rows.append((f"session 60 pinned-split claim, {_key}",
                         rf"({_pre}){re.escape(_fmt.format(table[_key]))}({_post})",
                         r"\g<1>" + _fmt.format(table[_key] + _d) + r"\g<2>"))
    # Session 51: the counter's other figures, each injected as a small drift.
    if "abstract_title_chars" in table:
        t = int(table["abstract_title_chars"])
        rows.append(("abstract title characters (live recount)",
                     rf"\(\s*{t} characters under the AACR rule",
                     f"({t + 3} characters under the AACR rule"))
    if "abstract_naive_len" in table:
        n = int(table["abstract_naive_len"])
        rows.append(("abstract length including spaces (live recount)",
                     rf"including spaces reads {n // 1000},?{n % 1000:03d}",
                     f"including spaces reads {n + 5:,}"))
    if "abstract_body_words" in table:
        w = int(table["abstract_body_words"])
        rows.append(("abstract body word count (live recount)",
                     rf"body\s+{w}\s+words", f"body {w + 1} words"))
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
        # ONLY THE DOCUMENTS A ROW CHANGED ARE RESCANNED. Until session 54 every
        # row rewrote and rescanned all of DOCS, so this step cost rows x
        # documents and was 80-88% of the whole suite, growing with every claim
        # this project added. `scan` keeps no state across documents -- each is
        # read from `root / rel` and compared alone -- so a document the row did
        # not touch contributes exactly its baseline, and the verdict below is
        # the same arithmetic on fewer terms. Measured, not argued: the first
        # run of this form reproduced the previous form's per-row verdicts and
        # mismatch counts for all 395 rows (14-SCIENCE-AUDIT.md, session 54).
        base_by_doc = {rel: scan(tmp, table, every_run, quiet=True,
                                 only={rel})[1] for rel in originals}
        assert sum(base_by_doc.values()) == base_fail, (
            "per-document baselines do not sum to the whole-tree baseline; "
            "scan has acquired cross-document state and this loop is no "
            "longer exact")
        missed = []
        spot_checks = 3
        for desc, pattern, replacement in rows:
            rx = re.compile(pattern)
            hits = 0
            changed: dict[str, str] = {}
            for rel, text in originals.items():
                new, n = rx.subn(replacement, text)
                if n:
                    hits += n
                    changed[rel] = new
            fail = 0
            if changed:
                for rel, new in changed.items():
                    (tmp / rel).write_text(new)
                _, fail, _, _, _ = scan(tmp, table, every_run, quiet=True,
                                        only=set(changed))
                if spot_checks:
                    # The exactness claim above, re-proved on every run for the
                    # first few rows that change anything: the whole-tree
                    # rescan the old loop did must give the same count.
                    spot_checks -= 1
                    _, whole, _, _, _ = scan(tmp, table, every_run, quiet=True)
                    assert whole - base_fail == fail - sum(
                        base_by_doc[rel] for rel in changed), (
                        f"row {desc!r}: rescanning only the changed documents "
                        "disagrees with rescanning the whole tree")
                for rel in changed:
                    (tmp / rel).write_text(originals[rel])
            delta = fail - sum(base_by_doc[rel] for rel in changed)

            if hits == 0:
                verdict = "NO SUCH TEXT (nothing to inject into)"
                missed.append(desc)
            elif delta > 0:
                verdict = f"CAUGHT ({delta} mismatch(es))"
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
    # Raised 157 -> 402 in session 54, which bound the Results' partition
    # arithmetic, Control C, the size-gap statistics and the rest, and made a
    # bracketed triple's unidentified point a failure.
    # Lowered 402 -> 386 in session 59 with the reason, as this table asks:
    # B-21 lifted Limitation 8's note (59 checked literals) into its own file,
    # floored below at 59; the manuscript and the note together hold 445.
    # Lowered 386 -> 290 on 2026-09-24 (B-21, after session 60), with the
    # reason: five Results passages and the site control's stress test were
    # lifted word for word into Supplementary Note 2, floored below at its
    # measured 170, and the summaries left in their place are bound by
    # B21_SUMMARY_SPEC. The manuscript held 398 checked before the lift; the
    # manuscript and the note together hold 460 after it, and no pattern lost
    # its last document (PATTERN REACH diffed before and after: 114 moved from
    # the manuscript to the note, none vanished).
    "09-PAPER-DRAFT.md": 290,
    # Session 59 (B-21): Limitation 8's note, lifted out of the manuscript on
    # 2026-09-24. Set at its measured count, 59 checked literals.
    "SUPPLEMENTARY-NOTE-1.md": 59,
    # B-21 (after session 60): the sensitivity analyses in full, lifted out of
    # the manuscript on 2026-09-24. Set at its measured count, 170.
    "SUPPLEMENTARY-NOTE-2.md": 170,
    "07-ABSTRACT-DRAFT.md": 25,
    # Raised 103 -> 546 in session 60 (H89): the sweep of all 1,072 unchecked
    # literals bound 288 to authorities; the rest are classified in
    # results/session60_diagnostics/audit_sweep_20260924.tsv.
    "14-SCIENCE-AUDIT.md": 546,
    "15-CHECKLIST.md": 17,
    "16-YOUR-TASKS.md": 9,
    "pipeline/README.md": 21,
    # Raised 19 -> 79 in session 54, which bound the rest of the board's
    # results (the lead sentence among them) after session 53 had bound its
    # outcome panel without moving this line. A floor left at 19 would let
    # sixty of those guards disappear without a word.
    "poster/poster.html": 79,
    "submission/SUBMISSION-CHECKLIST.md": 8,
    # Raised 7 -> 18 in session 54, by patterns written for the manuscript
    # that also reach the letter's copies of the same numbers.
    "submission/COVER-LETTER.md": 18,
    # Lowered 58 -> 56 on 2026-09-08, with the guard found rather than
    # silenced. Session 41 wrote a new run row for `_run4` that stated
    # "131 tests" in FOUR live-claim phrasings, and the sentence "an idle
    # figure at 131 tests still does not exist" in five rows. Those are
    # HISTORICAL statements -- a run row records what a named run reported
    # on a date -- but they were phrased so the checker read them as live,
    # so every future growth of the suite would go red and demand an edit
    # to a historical record. This file's own convention is the hyphenated
    # form (`131-test`), which the pattern cannot see. Restoring the
    # convention deliberately removes those literals from the checked set.
    "pipeline/results/README.md": 56,
    "CITATION.cff": 2,
    "18-SUPPLEMENTARY-INVENTORY.md": 8,
    "11-SUBMISSION-PACK.md": 3,
    # Added 2026-09-08 with the document itself. Only 3 of its 160 literals
    # are checked -- most are dates, fees and portal field numbers with no
    # frozen artefact behind them -- but one of the three is the abstract
    # character count, which is the number submission is hard-blocked on.
    "submission/ABSTRACT-PORTAL-DRY-RUN.md": 3,
    # Added 2026-09-09 with the document. 13 of its 35 literals are checked --
    # 37%, the highest fraction of any document here, because a deposit README
    # is almost entirely results. Note WHY this entry has to exist at all: the
    # loop below iterates over COVERAGE_FLOOR, not over DOCS, so a document
    # added to DOCS and forgotten here is never checked for LOST coverage and
    # `--coverage --strict` stays green while its guards disappear. That gap is
    # now closed by the reconciliation in `_check_coverage_floor_is_complete`.
    # Raised 13 -> 16 in session 54 with the three landing-page claims bound.
    "snapshot/README.md": 16,
    # Added 2026-09-17 (F3.3) at the counts measured when each joined DOCS.
    # PREPRINT-PLAN's only checked literal was a stale test count, reworded to
    # carry no number, so its floor is 0: it is scanned for mismatches, not for
    # coverage.
    "08-READINESS-ASSESSMENT.md": 2, "10-PROJECT-REFRAME.md": 1,
    "12-MENTOR-OUTREACH.md": 2, "13-NEXTGEN-EXTENDED-ABSTRACT.md": 12,
    "17-CPTAC-PLAN.md": 10, "pipeline/ENVIRONMENT.md": 2,
    "pipeline/results/supplementary/CAPTIONS.md": 1,
    "pipeline/results/supplementary/MANIFEST.md": 4,
    "submission/PREPRINT-PLAN.md": 0,
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
    # RECONCILE FIRST. The loop below ranges over COVERAGE_FLOOR, so a document
    # that is in DOCS but absent here is silently exempt from the floor: its
    # checked count can fall to zero and this function still prints OK. That is
    # not hypothetical -- `snapshot/README.md` was added to DOCS on 2026-09-09
    # and had no floor entry until this reconciliation was written, in the same
    # hour, for the same reason. The set a guard ranges over is a thing that has
    # to be checked, not assumed.
    unfloored = sorted(set(DOCS) - set(COVERAGE_FLOOR) - set(EXEMPT))
    if unfloored:
        print("\nCOVERAGE FLOOR INCOMPLETE -- these documents are scanned but "
              "have no floor,\nso losing every guard in them would not be "
              "noticed:")
        print("\n".join(f"  {r}" for r in unfloored))
        print("\nAdd each to COVERAGE_FLOOR at its currently measured count, "
              "or to EXEMPT\nwith the reason it is not checked.")
        return 1
    stale = sorted(set(COVERAGE_FLOOR) - set(DOCS))
    if stale:
        print("\nCOVERAGE FLOOR IS STALE -- these have a floor but are no "
              "longer scanned:")
        print("\n".join(f"  {r}" for r in stale))
        return 1

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


def _absent_claim_keys(table: dict[str, float]) -> list[str]:
    """Claim keys with no authority value -- claims `scan` skips WITHOUT a trace.

    Every claim family in `scan` begins with `if key not in table: continue`,
    and for CLAIMS that `continue` runs BEFORE the match count is recorded, so
    `--strict`'s vacuity check cannot see the row at all: a pattern whose
    authority vanished is neither checked nor reported. Measured 2026-09-17
    (session 50, ledger F3.4): zero such keys on the working tree without
    `--fast`; with `--fast`, exactly the five live-measured keys, by design.
    """
    keys = {k for k, *_ in CLAIMS} | {k for k, *_ in WORD_CLAIMS}
    keys |= {k for ks, *_ in SET_CLAIMS for k in ks}
    keys |= {k for row in TRIPLES for k in row[:3]}
    return sorted(keys - set(table))


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
    # `--fast` skips the live-measured authorities on purpose, so their claims
    # are expected to be absent there and only there.
    absent = _absent_claim_keys(table) if args.strict and not args.fast else []
    if absent:
        print(f"\nFAILED (--strict): {len(absent)} claim key(s) have no "
              "authority value, so their\npatterns were skipped without being "
              "counted as matching or not:")
        for k in absent:
            print(f"  NO AUTHORITY  {k}")
        return 1
    print("\nOK: every checked claim agrees with the frozen results.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
