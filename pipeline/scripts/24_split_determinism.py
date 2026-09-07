#!/usr/bin/env python3
"""A9, CLOSED: the cross-validation partition is not platform-invariant.

WHAT THIS SCRIPT ESTABLISHES, BY MEASUREMENT ON TWO MACHINES
------------------------------------------------------------
`splits._assign_sites_greedy` orders sites largest-first with

    order = per_site.index.to_numpy()[np.argsort(-(sizes + noise))]   # :177

`np.argsort` defaults to `kind="quicksort"` (introsort), which is NOT STABLE.
With the default `jitter=0.0` the `noise` term is exactly zero, so sites of
equal size are exactly tied and their relative order is decided entirely by the
sort's internal tie-breaking. That is implementation-defined: numpy ships
different SIMD kernels and pivot choices for different architectures, so arm64
and x86_64 order the tied block differently.

In NSCLC, 51 of 68 sites sit in a tied size group. A different order feeds the
greedy assignment a different sequence, producing a DIFFERENT site->fold
partition -- while preserving fold SIZES exactly, which is why every existing
guard (`_assert_site_disjoint`, the degenerate-partition check, the equal-byte
size of the two bisection dumps) passed and six sessions did not see it.

MEASURED, 2026-09-05, macOS arm64 py3.13.9 vs Einstein HPC4 x86_64 py3.12.14,
numpy 2.1.3 / pandas 2.2.3 / sklearn 1.6.1 / scipy 1.15.3 on BOTH:

    H_X_bytes            b5d70e2e64ef4ad6   b5d70e2e64ef4ad6   identical
    H_tss                5292af3ef56ccbef   5292af3ef56ccbef   identical
    H_sizes              1085bc3b226f462f   1085bc3b226f462f   identical
    H_order_quicksort    9fd91de2c6f61b19   0e7c4c7488097d58   DIFFERS
    H_order_stable       654c68f3e3aa46a2   654c68f3e3aa46a2   identical
    H_FOLD (shipped)     e55ee4e8887b5167   a40c862b12d9a391   DIFFERS
    H_FOLD (stable)      9f31e33c2caaf1e3   9f31e33c2caaf1e3   identical

And on the millisecond reproducer (`tests/_synthetic_isi`, whose 20 synthetic
sites hold exactly 10 patients each, so EVERY site is tied):

                        macOS                HPC4                 gap
    shipped             1.036173699004058    0.9752094760755513   5.89%
    stable sort         0.9752094760755506   0.9752094760755513   7e-16

Applying the stable sort ON MACOS reproduces the HPC4 value to 15 significant
figures. The entire 5.89% "platform gap" was the sort order. Genuine
cross-platform floating-point noise is ~7e-16 -- six orders of magnitude below
what was attributed to it.

WHAT THIS RETIRES
-----------------
A9 was framed for six sessions as "what mechanism turns a 1e-14 input
difference into an O(0.1) prediction difference". That premise is FALSE and was
never measured. Measured here:

  * cond(K + alpha*I) for the per-fold ridge system is 18-77, not ~1e10. The
    raw kernel is exactly singular (column-standardising X puts the all-ones
    vector in null(X^T)), but the REGULARISED system that is actually solved is
    well conditioned.
  * Perturbing the inputs by 1e-14 relative moves the predictions by 6e-15 --
    an amplification factor of ~0.25. The ridge solve is a CONTRACTION.
  * The observed cross-platform output difference of 4.245e-01 against a
    1.010e-14 input difference implies an amplification of 4.2e+13, which this
    system cannot produce. The machines were not solving the same problem.

So `models.cross_val_predict_multi` is exonerated, and so is BLAS. The alpha
flips (11 of 80, each one grid step) are a CONSEQUENCE of the different
partition -- different training sets select different penalties -- not a cause.

THE FIX, VERIFIED BUT NOT APPLIED
---------------------------------
`kind="stable"` at splits.py:177 makes the partition byte-identical across both
architectures (H_FOLD 9f31e33c2caaf1e3 on both). It is NOT applied here,
because changing the partition changes every frozen number in
`results/nsclc_v3/` and `results/pancancer_v3/`. That is the maintainer's call.

Two further non-stable sorts are latent, not implicated on this data:
  * `signatures.py:212` `np.argsort(-ranks[i])` -- its own comment already notes
    tied genes are ordered by argsort's tie-breaking. Empirically harmless here:
    the residualised scores downstream (S2) agree across platforms to 1.0e-14,
    so continuous expression evidently produces no exact rank ties.
  * `stats.py:180` `np.argsort(p)` in the BH-FDR path -- tied p-values.
`outcome.py:173,181` already pass `kind="stable"` deliberately.

USAGE
-----
    python3 scripts/24_split_determinism.py            # report every hash
    python3 scripts/24_split_determinism.py --check    # assert invariance
    python3 scripts/24_split_determinism.py --write    # durable JSON artefact

`--check` asserts the STABLE-sort fold hash equals the value measured on both
machines. It is designed to be able to fail: change splits.py's assignment logic
and it goes red. It deliberately does NOT assert the shipped hash, because that
value is platform-dependent -- which is the whole finding.

`--write` records what THIS machine measures to
`results/split_determinism_<system>_<machine>.json`. One file per platform, each
produced by an actual run on that platform -- so the cross-platform table above
is quoted from artefacts rather than from memory, and `19_check_numbers.py` can
check it. Run it on both machines and keep both files.
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

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
INTERIM = REPO / "data" / "interim"

from aacr27 import splits  # noqa: E402

# Measured on macOS arm64 AND Einstein HPC4 x86_64. Platform-invariant.
STABLE_FOLD_HASH = "9f31e33c2caaf1e3"
STABLE_ORDER_HASH = "654c68f3e3aa46a2"


def h(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]


def _src_line() -> str:
    """The shipped `splits.py:177`, read from disk.

    Read rather than remembered: the printed labels in this script once claimed
    the quicksort order was "shipped" after the source had already been fixed.
    A label derived from the source cannot drift away from it.
    """
    src = (REPO / "src" / "aacr27" / "splits.py").read_text().splitlines()
    for line in src:
        if "np.argsort(-(sizes + noise)" in line:
            return line
    return "<the argsort at splits.py:177 was not found>"


def _src_is_fixed() -> bool:
    return 'kind="stable"' in _src_line()


def _stable_greedy(per_site, n_folds, *, seed=0, jitter=0.0):
    """`_assign_sites_greedy` with the ONE change: kind="stable".

    Kept as a faithful copy rather than a monkeypatch of numpy, so what is
    being compared is exactly the shipped algorithm with a different tie-break.
    """
    rng = np.random.default_rng(seed)
    sizes = per_site["n"].to_numpy(dtype=float)
    noise = rng.uniform(0, jitter * max(sizes.max(), 1.0), size=len(sizes))
    order = per_site.index.to_numpy()[np.argsort(-(sizes + noise), kind="stable")]
    has_outcome = "outcome" in per_site.columns
    grand_mean = float(per_site["outcome"].mean()) if has_outcome else 0.0
    acc, sums, counts = (np.zeros(n_folds) for _ in range(3))
    assignment: dict[str, int] = {}
    for site in order:
        n = float(per_site.at[site, "n"])
        cost = acc + n
        cost = (cost - cost.mean()) / (cost.std() + 1e-9)
        if has_outcome:
            y = float(per_site.at[site, "outcome"])
            new_means = (sums + y * n) / np.maximum(counts + n, 1e-9)
            oc = np.abs(new_means - grand_mean)
            cost = cost + (oc - oc.mean()) / (oc.std() + 1e-9)
        best = np.flatnonzero(cost == cost.min())
        choice = int(best[0]) if len(best) == 1 else int(rng.choice(best))
        assignment[site] = choice
        acc[choice] += n
        counts[choice] += n
        if has_outcome:
            sums[choice] += float(per_site.at[site, "outcome"]) * n
    return assignment


def synthetic_isi(stable: bool) -> float:
    """The millisecond A9 reproducer, optionally with the stable-sort fix.

    A faithful copy of `tests/test_pipeline.py::_synthetic_isi`, which is what
    the frozen constants there pin. Its 20 sites hold exactly 10 patients each,
    so EVERY site is tied on size -- which is why this probe reproduces A9
    despite n=200 > p=40, nowhere near underdetermined.
    """
    from aacr27 import experiment, models, stats

    rng = np.random.default_rng(0)
    n, p, n_null = 200, 40, 8
    X = rng.standard_normal((n, p))
    beta = rng.standard_normal(p)
    obs = X @ beta + 2.0 * rng.standard_normal(n)
    nulls = rng.standard_normal((n, n_null)) + 0.15 * (X @ beta)[:, None]
    patients = pd.Series([f"P{i:04d}" for i in range(n)])
    sites = pd.Series([f"S{i % 20:02d}" for i in range(n)], dtype="string")

    orig = splits._assign_sites_greedy
    if stable:
        splits._assign_sites_greedy = _stable_greedy
    try:
        sp = splits.preserved_site_split(patients, sites, n_folds=5, seed=0)
    finally:
        splits._assign_sites_greedy = orig

    alpha = experiment._select_alpha_per_fold(X, obs, sp)
    pred = models.cross_val_predict_multi(
        X, np.column_stack([obs, nulls]), sp, patients.to_numpy(),
        fixed_alpha=alpha)
    r_obs = stats.corr_ci(obs, pred[:, 0]).value
    r_null = np.asarray([stats.corr_ci(nulls[:, j], pred[:, j + 1]).value
                         for j in range(n_null)], dtype=float)
    return float(stats.fisher_z(r_obs) - np.nanmean(stats.fisher_z(r_null)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="assert the stable-sort partition matches the "
                         "value measured on both platforms; exit 1 if not")
    ap.add_argument("--write", action="store_true",
                    help="write this platform's measurements to "
                         "results/split_determinism_<system>_<machine>.json")
    args = ap.parse_args()

    print(f"platform  {platform.system()} {platform.machine()} "
          f"py{platform.python_version()}")
    print(f"numpy {np.__version__}  pandas {pd.__version__}\n")

    cohort = INTERIM / "cohort_nsclc.parquet"
    if not cohort.exists():
        print(f"MISSING {cohort} -- this script needs the staged NSCLC inputs.")
        return 2

    frame = pd.read_parquet(cohort)
    sites = frame["tss"].astype("string")
    patients = frame["patient_id"]
    assert len(frame) == 944, f"expected 944 NSCLC rows, got {len(frame)}"

    per_site = pd.DataFrame({"site": sites, "patient": patients}).groupby(
        "site").agg(n=("patient", "nunique"))
    sizes = per_site["n"].to_numpy(dtype=float)
    vals, counts = np.unique(sizes, return_counts=True)
    tied = int(counts[counts > 1].sum())
    print(f"sites {len(per_site)}   distinct sizes {len(vals)}   "
          f"in a TIED size group: {tied}")
    print(f"H_site_index      {h(per_site.index.to_numpy().astype('U16'))}")
    print(f"H_sizes           {h(sizes)}")

    # Both orders are still computed and printed, even though splits.py:177 now
    # pins stable: they are the evidence that the two tie conventions really do
    # disagree on this data, which is what makes the fix load-bearing rather
    # than cosmetic. `_src_is_fixed` reads the shipped line so the labels can
    # never drift out of step with the source the way they did once already.
    fixed = _src_is_fixed()
    q = np.argsort(-sizes)
    s = np.argsort(-sizes, kind="stable")
    print(f"\nH_order_quicksort {h(q.astype(np.int64))}   <- platform-DEPENDENT"
          f"{'' if fixed else '  (SHIPPED)'}")
    print(f"H_order_stable    {h(s.astype(np.int64))}   <- invariant"
          f"{'  (SHIPPED)' if fixed else ''}")
    print(f"quicksort == stable here? {bool(np.array_equal(q, s))}")

    shipped = splits.preserved_site_split(
        patients, sites, n_folds=5, seed=0).fold.to_numpy().astype(np.int64)
    orig = splits._assign_sites_greedy
    splits._assign_sites_greedy = _stable_greedy
    try:
        stable = splits.preserved_site_split(
            patients, sites, n_folds=5, seed=0).fold.to_numpy().astype(np.int64)
    finally:
        splits._assign_sites_greedy = orig

    h_shipped, h_stable = h(shipped), h(stable)
    moved = int((shipped != stable).sum())
    print(f"\nH_FOLD shipped    {h_shipped}   sizes "
          f"{[int((shipped == k).sum()) for k in range(5)]}")
    print(f"H_FOLD stable     {h_stable}   sizes "
          f"{[int((stable == k).sum()) for k in range(5)]}")
    print(f"patients whose fold id differs between the two: {moved} of "
          f"{len(shipped)}")

    isi_shipped = synthetic_isi(stable=False)
    isi_stable = synthetic_isi(stable=True)
    print(f"\nsynthetic ISI shipped  {isi_shipped!r}")
    print(f"synthetic ISI stable   {isi_stable!r}")

    if args.write:
        out = (REPO / "results" /
               f"split_determinism_{platform.system().lower()}_"
               f"{platform.machine().lower()}.json")
        payload = {
            "platform": {
                "system": platform.system(),
                "machine": platform.machine(),
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
            "n_sites": int(len(per_site)),
            "n_sites_tied": tied,
            "n_patients": int(len(shipped)),
            "n_patients_fold_changed": moved,
            "hash_site_index": h(per_site.index.to_numpy().astype("U16")),
            "hash_sizes": h(sizes),
            "hash_order_quicksort": h(q.astype(np.int64)),
            "hash_order_stable": h(s.astype(np.int64)),
            "hash_fold_shipped": h_shipped,
            "hash_fold_stable": h_stable,
            "synthetic_isi_shipped": isi_shipped,
            "synthetic_isi_stable": isi_stable,
        }
        out.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"\nwrote {out.relative_to(REPO)}")

    if not args.check:
        if not args.write:
            print("\n(run with --check to assert platform invariance)")
        return 0

    ok = True
    if h(s.astype(np.int64)) != STABLE_ORDER_HASH:
        print(f"\nFAILED: stable site ORDER hash {h(s.astype(np.int64))} != "
              f"{STABLE_ORDER_HASH}")
        ok = False
    if h_stable != STABLE_FOLD_HASH:
        print(f"\nFAILED: stable FOLD hash {h_stable} != {STABLE_FOLD_HASH}\n"
              "The site-to-fold assignment changed. If that was intended, "
              "re-measure on BOTH machines and update the constant.")
        ok = False
    # The regression guard proper. Until 2026-09-05 this was a NOTE, because
    # splits.py:177 was deliberately unfixed and a hard failure would have made
    # the script unusable in reproduce.sh. The fix is applied now, so a shipped
    # partition that stops matching the stable one means the fix was reverted --
    # exactly the regression this script exists to catch. Fail on it.
    if h_shipped != h_stable:
        print(f"\nFAILED: the SHIPPED partition ({h_shipped}) differs from the "
              f"stable-sort partition ({h_stable}).\n"
              f"{int((shipped != stable).sum())} of {len(shipped)} patients are "
              "in a different fold. splits.py:177 must pass kind=\"stable\"; "
              f"the shipped line currently reads:\n  {_src_line().strip()}")
        ok = False
    if ok:
        print(f"\nOK: the stable-sort partition matches {STABLE_FOLD_HASH}, "
              "the value measured on macOS arm64 and HPC4 x86_64, and the "
              "SHIPPED partition is identical to it.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
