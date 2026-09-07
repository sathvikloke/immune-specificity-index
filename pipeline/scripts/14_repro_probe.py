#!/usr/bin/env python3
"""Which stage of the NSCLC ISI is not platform-stable?

    python scripts/14_repro_probe.py

WHY THIS EXISTS
===============
The scorer-sensitivity run reproduced the stored mean-z ISI on macOS/Python
3.13.9 (+0.3182) but returned +0.2964 on Einstein HPC4 (Linux/Python 3.12.14)
with byte-identical inputs and every package pinned to the same version --
numpy 2.1.3, pandas 2.2.3, scikit-learn 1.6.1, scipy 1.15.3. A 0.022 shift in a
headline estimate across platforms is a reproducibility defect, and the
repository claims one-command reproduction, so it has to be located rather than
absorbed.

This prints a fingerprint of each stage in order, so the FIRST line that differs
between two machines names the culprit. It computes no model and takes seconds.

Stages, in dependency order:
  1. inputs        -- are the parquet/npy files actually identical?
  2. column order  -- does the gene axis arrive in the same order?
  3. mean express  -- the vector that drives expression-matched binning
  4. qcut bins     -- the binning itself, which is rank-based and tie-sensitive
  5. null draws    -- the sampled random gene sets
  6. scores        -- mean-z scores for one real signature
Run it on both machines and diff the output.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aacr27 import signatures  # noqa: E402

INTERIM = ROOT / "data" / "interim"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"


def h(obj) -> str:
    """Stable short digest of an array or sequence of strings."""
    if isinstance(obj, np.ndarray):
        b = np.ascontiguousarray(obj).tobytes()
    else:
        b = "\n".join(map(str, obj)).encode()
    return hashlib.sha256(b).hexdigest()[:16]


def file_digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def main() -> int:
    print(f"python   {sys.version.split()[0]}")
    print(f"numpy    {np.__version__}")
    print(f"pandas   {pd.__version__}")
    try:
        import sklearn
        print(f"sklearn  {sklearn.__version__}")
    except ImportError:
        pass
    print(f"blas     {np.__config__.get_info('blas_opt_info') if hasattr(np.__config__, 'get_info') else 'n/a'}"[:120])
    print()

    # -- 1. inputs -------------------------------------------------------
    for name in ("cohort_nsclc.parquet", "X_nsclc.npy", "expr_nsclc.parquet"):
        print(f"1 input      {name:24s} {file_digest(INTERIM / name)}")
    print(f"1 input      {'h.all...tme.gmt':24s} {file_digest(GMT)}")

    cohort = pd.read_parquet(INTERIM / "cohort_nsclc.parquet")
    expr = pd.read_parquet(INTERIM / "expr_nsclc.parquet")
    if not expr.index.equals(pd.Index(cohort["patient_id"])):
        expr = expr.reindex(cohort["patient_id"])

    # -- 2. axis order ---------------------------------------------------
    print(f"2 columns    {expr.shape}                  {h(list(expr.columns))}")
    print(f"2 index      {'':24s} {h(list(expr.index))}")

    # -- 3. mean expression ----------------------------------------------
    mean_expr = expr.mean(axis=0)
    print(f"3 mean_expr  sum={mean_expr.sum():.12f}   {h(mean_expr.to_numpy())}")
    # Exact-tie count matters: qcut ranks with method='first', so ties are broken
    # by position and are stable, but a 1-ulp difference in a mean is not.
    print(f"3 mean_expr  n_exact_zero={int((mean_expr == 0).sum())}  "
          f"n_unique={mean_expr.nunique()}")

    # -- 4. binning ------------------------------------------------------
    pool = np.asarray(list(expr.columns))
    shared = mean_expr.reindex(pool).dropna()
    bins = pd.qcut(shared.rank(method="first"), 10, labels=False).astype(int)
    print(f"4 qcut bins  counts={np.bincount(bins.to_numpy()).tolist()}")
    print(f"4 qcut bins  {'':24s} {h(bins.to_numpy())}")

    # -- 5. null draws ---------------------------------------------------
    sigs = signatures.SignatureSet.from_gmt(GMT)
    sigs = signatures.SignatureSet(
        name=sigs.name, sets={f"SIG_{k}": v for k, v in sigs.sets.items()})
    filtered = sigs.filter_to(set(expr.columns))
    print(f"5 filtered   {len(filtered)} signatures        "
          f"{h([f'{k}:{len(v)}' for k, v in filtered.sets.items()])}")

    one = signatures.SignatureSet(
        name="probe", sets={"p": next(iter(filtered.sets.values()))})
    fam = signatures.random_gene_sets(
        list(expr.columns), one, n_per_signature=20, seed=0,
        match_expression=mean_expr)["p"]
    flat = [g for k in sorted(fam.sets) for g in fam.sets[k]]
    print(f"5 null draws 20 sets, {len(flat)} genes      {h(flat)}")
    print(f"5 null draw0 first 5: {fam.sets[sorted(fam.sets)[0]][:5]}")

    # -- 6. scores -------------------------------------------------------
    scored = signatures.score_mean_z(expr, filtered)
    col = sorted(scored.columns)[0]
    print(f"6 mean_z     {col}")
    print(f"6 mean_z     sum={scored[col].sum():.12f}  {h(scored[col].to_numpy())}")

    # -- 7. ridge alpha selection ----------------------------------------
    # THE SUSPECT. models.ALPHAS is np.logspace(-2, 5, 24) and RidgeCV picks per
    # target by leave-one-out generalised CV. Adjacent grid points differ by a
    # factor of ~1.96, so when two alphas score near-identically, a 1-ulp
    # difference in the target vector flips the selection to a penalty twice as
    # heavy -- and requirements-lock.txt already warns that "the ISI depends on
    # the selected penalty". If this line differs across machines while stages
    # 1-5 match, the ~0.02 ISI shift is explained: it is not drift, it is a
    # discrete jump.
    #
    # BUT READ A MATCH HERE NARROWLY. This is ONE fit on the FULL cohort over the
    # 16 observed signatures. The audit fits RidgeCV inside each of 5 folds, under
    # two split schemes, for the observed signatures AND for 16 x n_null random
    # sets -- of order 1e5 selections, against 16 here. Identical alphas at this
    # single fit therefore do NOT clear alpha selection; they only say the easiest
    # case is stable. Ruling it out needs per-fold alphas from the real path.
    # (2026-09-02: this line hashed identically on macOS and HPC4. That was
    # briefly over-read as exonerating the ridge. It does not.)
    from sklearn.linear_model import RidgeCV  # noqa: PLC0415

    from aacr27 import models  # noqa: PLC0415

    X = np.load(INTERIM / "X_nsclc.npy")
    cols = sorted(scored.columns)
    Y = scored[cols].to_numpy()
    Xs = (X - X.mean(0)) / X.std(0).clip(min=1e-12)
    est = RidgeCV(alphas=models.ALPHAS, alpha_per_target=True).fit(Xs, Y)
    alphas = np.atleast_1d(est.alpha_)
    print(f"7 X          {X.shape}                   {h(X)}")
    print(f"7 alphas     {h(alphas)}")
    for name, a in zip(cols, alphas):
        print(f"7   {name.replace('SIG_HALLMARK_', ''):42s} alpha={a:12.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
