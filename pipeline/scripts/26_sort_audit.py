#!/usr/bin/env python3
"""Audit every non-stable ordering in the pipeline for the A9 defect class.

WHY THIS SCRIPT EXISTS
----------------------
A9 was `splits.py:177`. `np.argsort` defaults to `kind="quicksort"` (introsort),
which is NOT STABLE, so exactly-tied keys are ordered by the sort's internal
tie-breaking -- implementation-defined, and different between numpy's arm64 and
x86_64 kernels. See `24_split_determinism.py` for the two-machine measurement.

Handoff #27 named two further non-stable sorts as "latent and NOT implicated on
this data" and rested that on an INFERENCE: downstream residualised scores agree
across platforms to 1.0e-14, therefore the sorts cannot be breaking ties
differently. That is an argument, not a measurement, and this project has been
burned twice by exactly that substitution. This script measures instead.

THE TWO QUESTIONS THAT DECIDE WHETHER THE DEFECT IS PRESENT
-----------------------------------------------------------
  1. Are there EXACT ties in the sort key on the real data?
     No ties  =>  the order is unique  =>  no sort algorithm can differ, on any
     platform.  This is a property of the data, checkable on one machine.
  2. If there are ties, does the OUTPUT change when the tie order changes?
     Some consumers are provably order-invariant -- a median, a mean, a set
     membership test.  Those are safe even with ties.

A case is DANGEROUS only when the answer is "ties exist" AND "output changes".
Both are measured here; neither is assumed.

MEASURED 2026-09-05, macOS arm64, numpy 2.1.3 / pandas 2.2.3, NSCLC inputs
--------------------------------------------------------------------------
  case                     ties?           output changes?   verdict
  A signatures.py:212      944/944 rows    YES, materially   WAS A LIVE BUG
  B stats.py:180 BH-FDR    forced, 32/32   no, even forced   SAFE, by algebra
  C barcodes.py:130        51 tied sites   yes (cum_frac)    DEAD CODE
  D models.py:339          0 of 15         no (median only)  SAFE, by consumer
  E splits.py:177          51 of 68        yes               A9, NOW FIXED

  [2026-09-22: the line numbers above are those of 2026-09-05. A and E keep
  theirs as the NAMES of defects A10 and A9. The cases below label B, C and D
  by function -- `stats.bh_fdr`, `barcodes.site_summary`,
  `models.site_prediction_control` -- because by 2026-09-22 D's line had
  drifted 85 lines, into `impute_median`. C's sort was pinned to
  kind="stable" on 2026-09-06; the case registers that through
  `shipped_pins_stable`, so its verdict above predates the pin.]

CASE A OVERTURNED WHAT HANDOFF #27 CARRIED FORWARD. It recorded this sort as
"empirically harmless here, since continuous expression evidently produces no
exact rank ties", inferred from downstream scores agreeing to 1.0e-14. The
inference was wrong and the data says so loudly: EVERY one of the 944 samples
carries ties and 35,394,448 of 38,747,424 rank entries (91%) are tied, because
log-expression has a floor and every gene sitting on it shares an average rank.
The tie comment in `signatures._ssgsea_sample_tables` said exactly this; the
handoff contradicted the source it was describing. [That comment was cited as
`signatures.py:206` until 2026-09-22, by which point it had moved to 239.]

The difference is not confined to the intermediate. `--deep` runs the real O(k)
reduction over the 16 real signatures under both tie conventions and measures
(macOS, 2026-09-05): 15,060 of 15,104 cells differ, all 16 signatures affected,
max |difference| 1.2824e-02 against a mean per-signature score sd of 2.1154e-02
-- 60.6% of a standard deviation. Tied genes carry identical weights, so a swap
only moves the score where exactly ONE of the pair is a set member; with 41,046
genes and ~200-gene sets that happens constantly.

Scope: this never touched the frozen headline results, which run
`scorer="mean_z"` (results/nsclc_v3/config.json) and never reach this function.
It did reach the ssGSEA arms of the A7 scorer-sensitivity analysis.

Case B is safe for a reason worth stating, because it is NOT "no ties": exactly
tied p-values receive adjacent ranks i and i+1, and the reverse running minimum
in the step-up assigns them the SAME q regardless of which came first. Swapping
them permutes `order` and `q` together, so `out_q` is bit-identical. This script
FORCES ties and confirms it rather than trusting the argument.

USAGE
-----
    python3 scripts/26_sort_audit.py              # audit, exit 1 if DANGEROUS
    python3 scripts/26_sort_audit.py --self-test  # prove each detector can fire

`--self-test` is the part that matters. A green audit you have never seen go red
is not evidence, so the self-test drives each detector on data whose answer is
known by construction and asserts it fires. The hash comparison is checked BOTH
ways -- it must differ on a tie-heavy matrix and agree on a tie-free one, since
a comparison that always differs is as useless as one that never does.

Verdicts are about the SHIPPED code, read from source by `pins_stable`, not
about the two algorithms. `output_changes` stays True for A and E after the fix:
quicksort and stable still disagree on the data, which is precisely the evidence
that pinning a convention mattered. What changed is that the shipped call now
pins one. Labels are derived from the source for a reason -- script 24's output
went on calling the quicksort order "shipped" after the source had been fixed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from aacr27 import stats as stats_mod  # noqa: E402

EXPR_NSCLC = REPO / "data" / "interim" / "expr_nsclc.parquet"
COHORT_NSCLC = REPO / "data" / "interim" / "cohort_nsclc.parquet"
SITE_CONTROL = REPO / "results" / "nsclc_v3" / "site_control.csv"
OUTCOME_ARM = REPO / "results" / "nsclc_v3" / "outcome_arm.csv"


def h(a: np.ndarray) -> str:
    """Stable 16-hex digest of an array's bytes, shape and dtype."""
    m = hashlib.sha256()
    m.update(str(a.shape).encode())
    m.update(str(a.dtype).encode())
    m.update(np.ascontiguousarray(a).tobytes())
    return m.hexdigest()[:16]


def count_row_ties(row: np.ndarray) -> int:
    """Number of elements that share their value with at least one other."""
    _, counts = np.unique(row, return_counts=True)
    return int(counts[counts > 1].sum())


def shipped_line(rel: str, needle: str) -> str:
    """The shipped source line containing `needle`, read from disk.

    Every verdict below that says "fixed" is derived from this rather than
    remembered. A label that is not read from the source is a label that can
    drift away from it -- which is exactly what happened to script 24's output
    labels, which went on calling the quicksort order "shipped" after the source
    had already been changed.
    """
    for line in (REPO / rel).read_text().splitlines():
        if needle in line:
            return line
    return f"<{needle!r} not found in {rel}>"


def pins_stable(rel: str, needle: str) -> bool:
    return 'kind="stable"' in shipped_line(rel, needle)


# --------------------------------------------------------------------------
# Case A -- signatures.py:212, the ssGSEA depth table.
# --------------------------------------------------------------------------
def case_a(expr: pd.DataFrame) -> dict:
    """`depth[i, np.argsort(-ranks[i])] = countdown`.

    Reproduces the call exactly as `signatures._rank_tables` makes it, on the
    real NSCLC expression matrix, and asks whether any row of the average-rank
    matrix contains an exact tie. Row by row: the full float64 rank matrix is
    944 x 41046 = 310 MB and there is a pan-cancer job resident.
    """
    ranks = expr.rank(axis=1, method="average").to_numpy(dtype=float)
    n_samples, n_genes = ranks.shape
    assert n_samples > 0 and n_genes > 0, f"degenerate expression matrix {ranks.shape}"

    countdown = np.arange(n_genes, 0, -1, dtype=np.int32)
    tied_rows, tied_elems = 0, 0
    mq = hashlib.sha256()
    ms = hashlib.sha256()
    for i in range(n_samples):
        t = count_row_ties(ranks[i])
        if t:
            tied_rows += 1
            tied_elems += t
        dq = np.empty(n_genes, dtype=np.int32)
        ds = np.empty(n_genes, dtype=np.int32)
        dq[np.argsort(-ranks[i])] = countdown                  # the shipped call
        ds[np.argsort(-ranks[i], kind="stable")] = countdown   # tie-order fixed
        mq.update(dq.tobytes())
        ms.update(ds.tobytes())

    hq, hs = mq.hexdigest()[:16], ms.hexdigest()[:16]
    return {
        "case": "A signatures.py:212  ssGSEA depth table",
        "n_samples": n_samples,
        "n_genes": n_genes,
        "tied_rows": tied_rows,
        "tied_elements": tied_elems,
        "h_quicksort": hq,
        "h_stable": hs,
        "output_changes": hq != hs,
        "ties_exist": tied_rows > 0,
        "shipped_pins_stable": pins_stable(
            "src/aacr27/signatures.py", "depth[i, np.argsort(-ranks[i]"),
    }


def case_a_deep(expr: pd.DataFrame) -> dict:
    """Does the depth-table difference SURVIVE into the 16 signature scores?

    Case A shows the intermediate differs. That is not yet a finding: tied genes
    carry identical weights, so a swap inside a tie group cancels in the O(k)
    reduction whenever BOTH or NEITHER of the pair is a set member. It moves the
    score only where exactly one is. This runs the real reduction both ways on
    the real gene sets and reports the difference against the score's own scale,
    because "1e-2 differs" means nothing without knowing the sd.

    Costs ~30 s on macOS. Not run by default; `--deep`.
    """
    from aacr27 import signatures as sig_mod

    gmt = REPO / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"
    if not gmt.exists():
        return {"case": "A-deep ssGSEA scores", "skipped": f"no {gmt.name}"}
    raw = sig_mod.SignatureSet.from_gmt(gmt)
    sigs = sig_mod.SignatureSet(name=raw.name,
                                sets={f"SIG_{k}": v for k, v in raw.sets.items()},
                                source=str(gmt))
    members = sig_mod._member_positions(expr.columns, sigs)
    assert members, "no signature genes matched the expression columns"

    n_samples, n_genes = expr.shape
    all_depth = n_genes * (n_genes + 1) / 2.0
    alpha = 0.25
    countdown = np.arange(n_genes, 0, -1, dtype=np.int32)
    out = {k: {n: np.empty(n_samples) for n in members} for k in ("quicksort", "stable")}
    block = max(1, sig_mod._SSGSEA_BLOCK_CELLS // max(n_genes, 1))

    for start in range(0, n_samples, block):
        stop = min(start + block, n_samples)
        ranks = expr.iloc[start:stop].rank(axis=1, method="average").to_numpy(dtype=float)
        w = np.abs(ranks) ** alpha
        for kind in ("quicksort", "stable"):
            depth = np.empty(ranks.shape, dtype=np.int32)
            for i in range(ranks.shape[0]):
                depth[i, np.argsort(-ranks[i], kind=kind)] = countdown
            for name, idx in members.items():
                n_miss = n_genes - len(idx)
                ww, dd = w[:, idx], depth[:, idx]
                total_hit = ww.sum(axis=1)
                with np.errstate(divide="ignore", invalid="ignore"):
                    sc = ((ww * dd).sum(axis=1) / total_hit
                          - (all_depth - dd.sum(axis=1, dtype=np.float64)) / n_miss) / n_genes
                out[kind][name][start:stop] = np.where(
                    (total_hit == 0) | (n_miss == 0), np.nan, sc)

    q = pd.DataFrame(out["quicksort"], index=expr.index)
    s = pd.DataFrame(out["stable"], index=expr.index)
    d = (q - s).abs().to_numpy()
    sd = float(np.nanmean(q.std().to_numpy()))
    return {
        "case": "A-deep  ssGSEA SCORES, both tie conventions",
        "n_signatures": q.shape[1],
        "cells_differing": int(np.nansum(d > 0)),
        "cells_total": int(d.size),
        "max_abs_diff": float(np.nanmax(d)),
        "mean_score_sd": sd,
        "max_diff_as_frac_of_sd": float(np.nanmax(d)) / sd if sd else float("nan"),
        "signatures_affected": int((np.nan_to_num(d) > 0).any(axis=0).sum()),
        "output_changes": bool(np.nanmax(d) > 0),
        "ties_exist": True,
        "shipped_pins_stable": pins_stable(
            "src/aacr27/signatures.py", "depth[i, np.argsort(-ranks[i]"),
    }


# --------------------------------------------------------------------------
# Case B -- stats.bh_fdr, Benjamini-Hochberg.
# --------------------------------------------------------------------------
def _bh_with_kind(p: np.ndarray, kind: str) -> np.ndarray:
    """`stats.bh_fdr`'s q-vector, with the sort algorithm swapped out.

    This is a COPY of the shipped body, so it can silently drift away from it and
    leave case B auditing a function nobody calls. `_assert_bh_copy_is_faithful`
    pins the copy to `stats.bh_fdr` and is run on every invocation.
    """
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p, kind=kind)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out_q = np.empty(n)
    out_q[order] = q
    return out_q


def _assert_bh_copy_is_faithful() -> None:
    """`_bh_with_kind(p, "quicksort")` must BE `stats.bh_fdr(p)`, not resemble it.

    Auditing a private copy of a function proves nothing about the function. The
    shipped `bh_fdr` uses the default (quicksort) sort, so that arm of the copy
    has to reproduce it bit for bit on tied and untied input alike.
    """
    rng = np.random.default_rng(7)
    for p in (rng.uniform(0, 1, 32),
              np.repeat(rng.uniform(0, 1, 8), 4),
              np.full(12, 0.03),
              np.array([0.001, 0.9])):
        _, q_ship = stats_mod.bh_fdr(p)
        q_copy = _bh_with_kind(p, "quicksort")
        assert np.array_equal(q_ship, q_copy), (
            "26_sort_audit's BH copy has drifted from stats.bh_fdr; case B is "
            f"auditing the wrong function.\n  shipped {q_ship}\n  copy    {q_copy}")


def case_b(*, inject_tie: bool = True) -> dict:
    """The real p-value families, plus a FORCED-tie adversarial family.

    `inject_tie=False` audits only the p-values the pipeline actually produced.
    The default is True because the interesting question is not "does this data
    happen to tie" but "would a tie break it", and the answer is no.
    """
    _assert_bh_copy_is_faithful()
    families: list[tuple[str, np.ndarray]] = []
    for name, path, col in (
        ("outcome_arm", OUTCOME_ARM, "p"),
        ("per_signature", REPO / "results" / "nsclc_v3" / "per_signature.csv", "p"),
    ):
        if path.exists():
            df = pd.read_csv(path)
            hit = [c for c in df.columns if c == col or c.startswith("p_")]
            for c in hit:
                v = pd.to_numeric(df[c], errors="coerce").dropna().to_numpy()
                if len(v) >= 2:
                    families.append((f"{name}:{c}", v))

    if inject_tie:
        # 32 p-values, every one of them duplicated: the worst case BH can see.
        rng = np.random.default_rng(0)
        base = rng.uniform(0, 1, 16)
        families.append(("FORCED all-tied (16 values x2)", np.repeat(base, 2)))
        families.append(("FORCED all-identical", np.full(32, 0.04)))

    worst_delta = 0.0
    ties_exist = False
    detail = []
    for name, v in families:
        t = count_row_ties(v)
        ties_exist |= t > 0
        qq = _bh_with_kind(v, "quicksort")
        qs = _bh_with_kind(v, "stable")
        d = float(np.max(np.abs(qq - qs))) if len(v) else 0.0
        worst_delta = max(worst_delta, d)
        detail.append((name, len(v), t, d))

    return {
        "case": "B stats.bh_fdr  Benjamini-Hochberg",
        "families": detail,
        "max_abs_q_delta": worst_delta,
        "output_changes": worst_delta > 0.0,
        "ties_exist": ties_exist,
    }


# --------------------------------------------------------------------------
# Case C -- barcodes.site_summary. DEAD CODE.
# --------------------------------------------------------------------------
def case_c() -> dict:
    """`.nunique().sort_values(ascending=False)` then a cumulative sum.

    pandas `Series.sort_values` also defaults to `kind="quicksort"`, and per-site
    patient COUNTS are massively tied -- far more tied than the sizes in A9. The
    `cum_frac` column is order-dependent by construction, so this is a live
    instance of the same defect.

    It is not a bug in the delivered analysis for one reason only: `site_summary`
    is never called. Measured by grep over every .py/.md/.sh/.ipynb in the repo:
    the sole occurrence is its own `def`. Recorded here so that "dead" is a
    measurement someone can re-run, not a memory.
    """
    called_from = []
    for pat in ("*.py", "*.md", "*.sh", "*.ipynb"):
        for f in REPO.rglob(pat):
            if ".git" in f.parts or f.name == "26_sort_audit.py":
                continue
            try:
                txt = f.read_text(errors="ignore")
            except OSError:
                continue
            for ln, line in enumerate(txt.splitlines(), 1):
                # A CALL, not a MENTION. Until 2026-09-16 this matched the bare
                # name, and the check defeated itself: `ENVIRONMENT.md` contains
                # the sentence "`site_summary` is never called", so writing down
                # that the function is dead was what made the detector report it
                # alive. The stored artefact said `is_dead: true` on 2026-09-05
                # and the next `--write` would have flipped it to false, silently
                # contradicting the manuscript's "is not called anywhere in the
                # pipeline" -- which was, and remains, TRUE.
                #
                # Requiring the open paren keeps the .md/.sh/.ipynb sweep, which
                # exists to catch a notebook or script that really does call it,
                # while no longer counting prose that merely names it. Measured
                # at the time of the change: one hit under the old rule (that
                # sentence), zero under this one.
                if "site_summary(" in line and not line.lstrip().startswith("def "):
                    called_from.append(f"{f.relative_to(REPO)}:{ln}")

    ties = -1
    changes = False
    if COHORT_NSCLC.exists():
        c = pd.read_csv(COHORT_NSCLC) if COHORT_NSCLC.suffix == ".csv" \
            else pd.read_parquet(COHORT_NSCLC)
        col = next((x for x in ("tss", "site") if x in c.columns), None)
        pid = next((x for x in ("patient_id", "patient") if x in c.columns), None)
        if col and pid:
            counts = c.groupby(col)[pid].nunique()
            ties = count_row_ties(counts.to_numpy())
            a = counts.sort_values(ascending=False, kind="quicksort")
            b = counts.sort_values(ascending=False, kind="stable")
            changes = not a.index.equals(b.index)

    return {
        "case": "C barcodes.site_summary  [DEAD CODE]",
        "call_sites": called_from,
        "tied_site_counts": ties,
        "output_changes": changes,
        "ties_exist": ties > 0,
        "dead": len(called_from) == 0,
        # Added 2026-09-06. This case reported no `shipped_pins_stable` at all,
        # so `verdict()` could only ever call it DANGEROUS -- pinning the sort
        # in barcodes.py left the audit's answer unchanged, and the audit went
        # on exiting 1 describing a defect that had been fixed. A case that
        # cannot register its own repair is not measuring the shipped code.
        "shipped_pins_stable": pins_stable(
            "src/aacr27/barcodes.py", ".sort_values(ascending=False"),
    }


# --------------------------------------------------------------------------
# Case D -- models.site_prediction_control, site AUROC ordering.
# --------------------------------------------------------------------------
def case_d() -> dict:
    """`pd.DataFrame(rows).sort_values("auroc", ascending=False)`.

    Row ORDER is non-stable under ties. What is derived from it is a median --
    `AuditResult.headline` recomputes it straight off the `auroc` column, and the
    reported figure is read from `site_control.csv` by `19_check_numbers.py` --
    and a median is permutation-invariant, so ties can exist without touching a
    reported number.

    [Corrected 2026-09-22. This docstring said, in the present tense, that "the
    `.attrs` written at models.py:341-342 are never read by anything". Nothing is
    written there any more: `median_auroc` and `frac_above_0.9` were removed from
    `site_prediction_control` on 2026-09-06 -- its own comment says why: dead on
    arrival, and `.attrs` does not survive `.to_csv()`, which is precisely how
    the rotation-null p-value was silently lost once. The docstring's line
    numbers had drifted too: by 2026-09-22 `models.py:339` landed in
    `impute_median` and `experiment.py:249` in a docstring.]
    """
    if not SITE_CONTROL.exists():
        return {"case": "D models.site_prediction_control  site AUROC order", "skipped": "no site_control.csv"}
    df = pd.read_csv(SITE_CONTROL)
    a = pd.to_numeric(df["auroc"], errors="coerce").dropna().to_numpy()
    ties = count_row_ties(a)
    q = np.sort(a, kind="quicksort")[::-1]
    s = np.sort(a, kind="stable")[::-1]
    return {
        "case": "D models.site_prediction_control  site AUROC order",
        "n_sites": len(a),
        "tied_aurocs": ties,
        "median_quicksort": float(np.median(q)),
        "median_stable": float(np.median(s)),
        "output_changes": float(np.median(q)) != float(np.median(s)),
        "ties_exist": ties > 0,
    }


# --------------------------------------------------------------------------
# Case E -- splits.py:177. The known live defect; standing positive control.
# --------------------------------------------------------------------------
def case_e() -> dict:
    """A9 itself. Must keep reporting DANGEROUS until splits.py:177 is fixed.

    If this case ever goes quiet while `splits.py:177` still reads
    `np.argsort(-(sizes + noise))`, the probe has broken, not the bug.
    """
    line = shipped_line("src/aacr27/splits.py", "np.argsort(-(sizes + noise)")
    fixed = pins_stable("src/aacr27/splits.py", "np.argsort(-(sizes + noise)")

    ties, changes = -1, False
    if COHORT_NSCLC.exists():
        c = pd.read_parquet(COHORT_NSCLC)
        col = next((x for x in ("tss", "site") if x in c.columns), None)
        pid = next((x for x in ("patient_id", "patient") if x in c.columns), None)
        if col and pid:
            sizes = c.groupby(col)[pid].nunique().to_numpy().astype(float)
            ties = count_row_ties(sizes)
            oq = np.argsort(-sizes)
            os_ = np.argsort(-sizes, kind="stable")
            changes = not np.array_equal(oq, os_)

    return {
        "case": "E splits.py:177  site-to-fold order  [A9]",
        "source_line": line.strip(),
        "tied_site_sizes": ties,
        "output_changes": changes,
        "ties_exist": ties > 0,
        "shipped_pins_stable": fixed,
    }


# --------------------------------------------------------------------------
# Case F -- scikit-learn's GroupKFold (A12). A library ordering, not ours.
# --------------------------------------------------------------------------
def _groupkfold_folds(groups: np.ndarray, *, stable: bool) -> np.ndarray:
    """GroupKFold(5) fold ids, with the library's own argsort or a stable one.

    The stable branch is sklearn 1.6.1's non-shuffle code with kind="stable"
    and nothing else changed (the same patch ledger item D1 measured).
    """
    from sklearn.model_selection import GroupKFold

    X = np.zeros((len(groups), 1))
    if not stable:
        fold = np.full(len(groups), -1)
        for k, (_, te) in enumerate(GroupKFold(n_splits=5).split(X, groups=groups)):
            fold[te] = k
        return fold
    _, group_idx = np.unique(groups, return_inverse=True)
    sizes = np.bincount(group_idx)
    order = np.argsort(sizes, kind="stable")[::-1]
    per_fold = np.zeros(5)
    to_fold = np.zeros(len(sizes))
    for i, w in enumerate(sizes[order]):
        j = np.argmin(per_fold)
        per_fold[j] += w
        to_fold[order[i]] = j
    return to_fold[group_idx].astype(int)


def case_f(splits_src: str | None = None, groups: np.ndarray | None = None) -> dict:
    """A12: `splits.random_patient_split` hands patients to sklearn's GroupKFold.

    GroupKFold orders groups with `np.argsort(n_samples_per_group)[::-1]` -- the
    default, unstable sort -- and after `collapse_to_patient` every group has
    size 1, so the whole order is the kernel's tie rule. The partition feeds the
    random-patient scheme (Delta r, Delta-MAE, the covariate baselines) and the
    site-prediction control, never the index. The library hashes are per
    platform by design.

    PINNED 2026-09-24 (the author's decision on A12, the memo's option b):
    `splits.stable_group_kfold` is GroupKFold's code with `kind="stable"`, and
    all three call sites use it, so the shipped verdict is FIXED. Until then it
    was DOCUMENTED, NOT PINNED. The case sees the pin by the helper's NAME in
    splits.py; the self-test checks the shipped (pinned) source, an injected
    unpinned source that must read as not pinned, and tie-free groups that must
    not change the folds.
    """
    import hashlib
    import inspect

    from sklearn.model_selection import GroupKFold

    lib_src = inspect.getsource(GroupKFold._iter_test_indices)
    lib_line = next((ln.strip() for ln in lib_src.splitlines() if "argsort(" in ln), "")
    if splits_src is None:
        splits_src = (REPO / "src" / "aacr27" / "splits.py").read_text()
    pinned = "GroupKFold" not in splits_src or "stable_group_kfold" in splits_src
    if groups is None:          # one row per patient, as after collapse_to_patient
        groups = np.random.default_rng(0).permutation(944)
    lib_folds = _groupkfold_folds(groups, stable=False)
    stable_folds = _groupkfold_folds(groups, stable=True)
    changes = not np.array_equal(lib_folds, stable_folds)
    hashes = {n: hashlib.sha256(np.argsort(np.ones(n, np.int64)).tobytes()).hexdigest()[:16]
              for n in (944, 7168)}
    return {
        "case": "F sklearn GroupKFold  random-patient group order  [A12]",
        "library_line": lib_line,
        "library_pins_stable": 'kind="stable"' in lib_line,
        "all_ties_argsort_hash_944": hashes[944],
        "all_ties_argsort_hash_7168": hashes[7168],
        # The fold assignment itself, both ways: the stable one must agree
        # across machines, the library one is free to differ.
        "fold_hash_library": hashlib.sha256(lib_folds.astype(np.int64).tobytes()).hexdigest()[:16],
        "fold_hash_stable": hashlib.sha256(stable_folds.astype(np.int64).tobytes()).hexdigest()[:16],
        "ties_exist": True,
        "output_changes": changes,
        "shipped_pins_stable": pinned,
        "documented": "A12 (pinned 2026-09-24 as splits.stable_group_kfold; the index "
                      "does not use this partition)",
    }


def verdict(r: dict) -> str:
    """The verdict is about the SHIPPED code, not about the two algorithms.

    `output_changes` stays True after a fix -- quicksort and stable still
    disagree on the data, which is the evidence that pinning one mattered. What
    changes is whether the shipped call pins a convention, so that is what
    decides the verdict.
    """
    if r.get("skipped"):
        return "SKIPPED"
    if r.get("dead"):
        return "DEAD CODE"
    if not r.get("ties_exist"):
        return "SAFE (no ties in the key)"
    if not r.get("output_changes"):
        return "SAFE (output is order-invariant)"
    if r.get("shipped_pins_stable"):
        return "FIXED (ties exist and matter; shipped pins kind=\"stable\")"
    if r.get("documented"):
        return f"DOCUMENTED, NOT PINNED -- {r['documented']}"
    return "DANGEROUS"


def report(results: list[dict]) -> int:
    bad = 0
    for r in results:
        v = verdict(r)
        print(f"\n{r['case']}")
        print(f"  verdict: {v}")
        for k, val in r.items():
            if k == "case":
                continue
            if k == "families":
                for name, n, t, d in val:
                    print(f"    family {name}: n={n} tied={t} max|dq|={d:.3e}")
            elif k == "call_sites":
                print(f"    call_sites: {val if val else 'NONE -- never called'}")
            else:
                print(f"    {k}: {val}")
        if v == "DANGEROUS":
            bad += 1
    print("\n" + "=" * 72)
    if bad:
        print(f"{bad} case(s) DANGEROUS: ties exist, they change the output, and")
        print("the shipped call does not pin a tie convention.")
    else:
        print("No dangerous non-stable ordering found.")
    return bad


def self_test() -> int:
    """Prove each detector can fire. Exits 1 if any detector is vacuous."""
    print("SELF-TEST -- inject a known tie and assert the detector reports it.\n")
    fails = 0

    # A: test the DETECTOR on data whose answer is known by construction, not by
    # perturbing the real matrix. An earlier version of this block injected one
    # tie into the real matrix and asserted the count went 0 -> 2. That premise
    # was false -- the real matrix already carries 35.4M ties -- so the assertion
    # could never hold and the self-test correctly reported itself VACUOUS. The
    # bug was in the test, and the finding it was meant to guard is the tie count
    # itself. Both halves are now checked.
    tie_free = np.arange(10, dtype=float)
    one_pair = np.array([0.0, 0.0, 1.0, 2.0, 3.0])
    all_same = np.zeros(7)
    counts_ok = (count_row_ties(tie_free) == 0
                 and count_row_ties(one_pair) == 2
                 and count_row_ties(all_same) == 7)
    print(f"  A tie counter       : tie-free={count_row_ties(tie_free)} "
          f"one-pair={count_row_ties(one_pair)} all-same={count_row_ties(all_same)}"
          f"  -> {'CAUGHT' if counts_ok else 'VACUOUS'}")
    fails += not counts_ok

    # And the hash comparison must be able to go BOTH ways: differ on a tie-heavy
    # matrix, agree on a tie-free one. A comparison that always differs is as
    # useless as one that never does.
    rng = np.random.default_rng(0)
    cols = [f"G{i}" for i in range(40)]
    heavy = pd.DataFrame(rng.normal(size=(6, 40)).round(0),
                         index=[f"S{i}" for i in range(6)], columns=cols)
    clean = pd.DataFrame(rng.normal(size=(6, 40)),
                         index=[f"S{i}" for i in range(6)], columns=cols)
    rh, rc = case_a(heavy), case_a(clean)
    hash_ok = rh["output_changes"] and not rc["output_changes"]
    print(f"  A hash comparison   : tie-heavy differs={rh['output_changes']} "
          f"tie-free differs={rc['output_changes']}"
          f"  -> {'CAUGHT' if hash_ok else 'VACUOUS'}")
    fails += not hash_ok

    # The finding itself, asserted so it cannot silently evaporate.
    expr = pd.read_parquet(EXPR_NSCLC)
    real = case_a(expr)
    real_ok = real["tied_rows"] == real["n_samples"] and real["output_changes"]
    print(f"  A on real NSCLC     : {real['tied_rows']}/{real['n_samples']} rows "
          f"tied, {real['tied_elements']} tied entries, depth differs="
          f"{real['output_changes']}  -> {'CAUGHT' if real_ok else 'VACUOUS'}")
    fails += not real_ok

    # B: the detector is `max_abs_q_delta`. Prove it can be non-zero at all by
    # feeding the same machinery an order-DEPENDENT reduction.
    v = np.array([0.01, 0.01, 0.5, 0.5])
    qq = _bh_with_kind(v, "quicksort")
    qs = _bh_with_kind(v, "stable")
    forced = case_b(inject_tie=True)
    tie_families = [f for f in forced["families"] if f[2] > 0]
    ok = len(tie_families) > 0 and np.allclose(qq, qs)
    print(f"  B tie detector      : {len(tie_families)} tied famil(ies) reached the "
          f"comparison -> {'CAUGHT' if ok else 'VACUOUS'}")
    fails += not ok

    # F: the library ordering must be SEEN to matter (ties exist and the fold
    # assignment changes under a stable sort), and the shipped splits.py --
    # pinned since 2026-09-24 (A12) -- must read as pinned, verdict FIXED. Until
    # then this check asserted the opposite (not pinned), which was the truth.
    f = case_f()
    ok = (f["ties_exist"] and f["output_changes"] and f["shipped_pins_stable"]
          and verdict(f).startswith("FIXED"))
    print(f"  F library ordering  : order-changes={f['output_changes']} "
          f"pinned={f['shipped_pins_stable']}  -> {'CAUGHT' if ok else 'VACUOUS'}")
    fails += not ok
    # ...and in the other two directions: a splits.py that still hands patients
    # to GroupKFold is reported as NOT pinned (so the pin detector is not stuck
    # at True), and groups of all-different sizes (no ties) change nothing.
    fu = case_f(splits_src="from sklearn.model_selection import GroupKFold\n")
    sizes = np.repeat(np.arange(20), np.arange(1, 21))   # group g has g+1 rows
    fn = case_f(groups=sizes)
    ok = (not fu["shipped_pins_stable"] and verdict(fu).startswith("DOCUMENTED")
          and not fn["output_changes"])
    print(f"  F other directions  : unpinned-source reads pinned={fu['shipped_pins_stable']} "
          f"tie-free changes={fn['output_changes']}  -> {'CAUGHT' if ok else 'VACUOUS'}")
    fails += not ok

    e = case_e()
    ok = e["ties_exist"] and e["output_changes"]
    print(f"  E positive control  : ties={e['tied_site_sizes']} order-changes="
          f"{e['output_changes']}  -> {'CAUGHT' if ok else 'VACUOUS'}")
    fails += not ok

    print()
    if fails:
        print(f"SELF-TEST FAILED: {fails} detector(s) could not be made to fire.")
    else:
        print("SELF-TEST PASSED: every detector fired on an injected tie.")
    return fails


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true",
                    help="inject a tie into each case and assert the detector fires")
    ap.add_argument("--deep", action="store_true",
                    help="also score the 16 real signatures under both tie "
                         "conventions (~30 s) -- the measurement that shows the "
                         "depth-table difference survives into the scores")
    ap.add_argument("--write", action="store_true",
                    help="record THIS machine's measurements to "
                         "results/sort_audit_<system>_<machine>.json, so the "
                         "A10 numbers quoted in the documents are checked "
                         "against an artefact rather than a transcript. "
                         "Implies --deep.")
    ap.add_argument("--library-order", type=Path, metavar="OUT",
                    help="run case F only (no expression needed) and record it, "
                         "with this machine's versions, to OUT -- a NEW file per "
                         "platform; never one of the sort_audit_*.json artefacts")
    args = ap.parse_args()

    if args.self_test:
        return 1 if self_test() else 0

    if args.library_order:
        import sklearn

        if args.library_order.name.startswith("sort_audit_"):
            print("refusing: --library-order does not write the A10 artefacts", file=sys.stderr)
            return 2
        f = case_f()
        f["verdict"] = verdict(f)
        f["platform"] = {"system": platform.system(), "machine": platform.machine(),
                         "python": platform.python_version(), "numpy": np.__version__,
                         "sklearn": sklearn.__version__}
        args.library_order.parent.mkdir(parents=True, exist_ok=True)
        args.library_order.write_text(json.dumps(f, indent=2) + "\n")
        report([f])
        print(f"wrote {args.library_order}")
        return 0

    if not EXPR_NSCLC.exists():
        print(f"missing {EXPR_NSCLC}", file=sys.stderr)
        return 2

    expr = pd.read_parquet(EXPR_NSCLC)
    assert expr.shape[0] > 0 and expr.shape[1] > 0, f"empty expression {expr.shape}"
    deep = args.deep or args.write
    results = [case_a(expr)]
    if deep:
        results.append(case_a_deep(expr))
    results += [case_b(), case_c(), case_d(), case_e(), case_f()]
    bad = report(results)

    if args.write:
        by = {r["case"].split()[0]: r for r in results}
        a, ad, b, c, e = by["A"], by["A-deep"], by["B"], by["C"], by["E"]
        out = (REPO / "results" /
               f"sort_audit_{platform.system().lower()}_"
               f"{platform.machine().lower()}.json")
        payload = {
            "platform": {
                "system": platform.system(),
                "machine": platform.machine(),
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
            # Case A -- the ssGSEA rank table. `hash_depth_stable` is the value
            # that must agree across machines; `hash_depth_quicksort` is the one
            # that must NOT, and both are recorded so the claim is falsifiable
            # in either direction.
            "n_samples": a["n_samples"],
            "n_genes": a["n_genes"],
            "tied_rows": a["tied_rows"],
            "tied_elements": a["tied_elements"],
            "hash_depth_quicksort": a["h_quicksort"],
            "hash_depth_stable": a["h_stable"],
            "signatures_pins_stable": a["shipped_pins_stable"],
            # Case A-deep -- does it survive into the scores.
            "ssgsea_cells_differing": ad["cells_differing"],
            "ssgsea_cells_total": ad["cells_total"],
            "ssgsea_signatures_affected": ad["signatures_affected"],
            "ssgsea_max_abs_diff": ad["max_abs_diff"],
            "ssgsea_mean_score_sd": ad["mean_score_sd"],
            "ssgsea_max_diff_as_frac_of_sd": ad["max_diff_as_frac_of_sd"],
            # Cases B, C, E -- the orderings that are safe, dead, or A9 itself.
            "bh_max_abs_q_delta": b["max_abs_q_delta"],
            "site_summary_is_dead": c["dead"],
            "site_summary_tied_counts": c["tied_site_counts"],
            # Added 2026-09-16. Cases A and E each persist their pin verdict;
            # case C COMPUTED one and threw it away, so a regression in
            # `barcodes.py`'s `kind="stable"` would have been invisible in this
            # artefact -- the exact failure these files exist to prevent. Three
            # cases carry a pin flag and only two were written down, with
            # nothing recording why.
            #
            # NOTE THE ASYMMETRY THIS CREATES, which is deliberate and is
            # recorded rather than tidied away: `sort_audit_darwin_arm64.json`
            # gains this key on its next `--write`, and
            # `sort_audit_linux_x86_64.json` does NOT. That file describes the
            # PRE-FIX cluster tree, where `signatures_pins_stable` and
            # `splits_pins_stable` are False -- and those two False values are
            # the evidence that the shipped cluster source does not pin, which
            # is what lets this script fire a true positive there. Re-writing it
            # would mean either overwriting a file in the HPC4 project directory
            # (forbidden) or running in the stable-sort checkout, whose output
            # lands on the SAME filename (it is derived from platform.system()
            # and platform.machine()) and would flip those flags to True. That
            # trades the evidence for a field, so it is declined.
            "site_summary_pins_stable": c["shipped_pins_stable"],
            "splits_tied_site_sizes": e["tied_site_sizes"],
            "splits_pins_stable": e["shipped_pins_stable"],
        }
        out.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"\nwrote {out.relative_to(REPO)}")

    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
