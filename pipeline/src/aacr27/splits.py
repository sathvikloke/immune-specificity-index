"""Cross-validation splits — random vs preserved-site.

The single most important module here. The study's primary endpoint is the
*difference* between performance under a naive random split and performance
under a site-disjoint split, so both must be constructed correctly and
comparably (same fold count, same patient-level grouping, same seed handling).

Preserved-site cross-validation follows Howard et al., Nat Commun 2021;12:4423
(https://www.nature.com/articles/s41467-021-24698-1), which showed that TCGA
tissue source site is trivially learnable from H&E and biases downstream
prediction. Their remedy assigns whole sites to folds via quadratic
programming so that outcome distribution stays balanced across folds.

We implement a greedy balanced assignment that is deterministic and dependency
free, and optionally an exact QP when `cvxpy` is installed. For the cohort
sizes here the greedy solution is typically within a few percent of optimal on
outcome balance; `fold_balance_report()` lets you verify that rather than
assume it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold


@dataclass(frozen=True)
class SplitResult:
    """Fold assignment plus the metadata needed to audit it."""

    fold: pd.Series  # index aligned to input, values 0..n_folds-1
    scheme: str
    n_folds: int
    balance: pd.DataFrame

    def indices(self, k: int) -> tuple[np.ndarray, np.ndarray]:
        """(train_idx, test_idx) positional arrays for fold k."""
        test = np.where(self.fold.to_numpy() == k)[0]
        train = np.where(self.fold.to_numpy() != k)[0]
        return train, test


def stable_group_kfold(groups, n_splits: int):
    """Yield (train_idx, test_idx) as scikit-learn 1.6.1's GroupKFold(n_splits).split does,
    with the group order pinned by a stable sort (A12, pinned 2026-09-24).

    GroupKFold's non-shuffle branch orders groups by
    `np.argsort(n_samples_per_group)[::-1]` -- the default, unstable sort. After
    `collapse_to_patient` every group has one row, so every key ties and the
    fold assignment is the platform's tie rule: macOS arm64 and Linux x86_64 give
    different random-patient folds for the same seed (14-SCIENCE-AUDIT.md, A12).
    This is that code with `kind="stable"` and nothing else changed -- the same
    patch ledger item D1 measured in-process
    (`results/session49_diagnostics/d1_pinned_groupkfold/`), under which the two
    platforms agree on every secondary quantity to round-off. Ties now fall in
    the reverse of `np.unique`'s order, on every platform.

    The name matters: `26_sort_audit.py` case F reads this module and reports the
    split as pinned only if `stable_group_kfold` appears in it.
    """
    if groups is None:
        raise ValueError("The 'groups' parameter should not be None.")
    groups = np.asarray(groups)
    if groups.ndim != 1:
        raise ValueError(f"groups must be one-dimensional, got shape {groups.shape}")
    n_samples = len(groups)
    if n_splits > n_samples:
        raise ValueError(
            f"Cannot have number of splits n_splits={n_splits} greater than the "
            f"number of samples: n_samples={n_samples}."
        )
    unique_groups, group_idx = np.unique(groups, return_inverse=True)
    if n_splits > len(unique_groups):
        raise ValueError(
            f"Cannot have number of splits n_splits={n_splits} greater than the "
            f"number of groups: {len(unique_groups)}."
        )
    # Distribute the most frequent groups first, each to the lightest fold.
    n_samples_per_group = np.bincount(group_idx)
    order = np.argsort(n_samples_per_group, kind="stable")[::-1]
    n_samples_per_fold = np.zeros(n_splits)
    group_to_fold = np.zeros(len(unique_groups))
    for group_index, weight in enumerate(n_samples_per_group[order]):
        lightest_fold = np.argmin(n_samples_per_fold)
        n_samples_per_fold[lightest_fold] += weight
        group_to_fold[order[group_index]] = lightest_fold
    sample_fold = group_to_fold[group_idx]
    indices = np.arange(n_samples)
    for f in range(n_splits):
        test = sample_fold == f
        yield indices[~test], indices[test]


def random_patient_split(
    patients: pd.Series,
    *,
    n_folds: int = 5,
    seed: int = 0,
    stratify: pd.Series | None = None,
) -> SplitResult:
    """Patient-disjoint but site-agnostic — i.e. what the literature usually does.

    Note this is already stricter than the worst practice (slide-level random
    splits, which leak a patient across train and test). We hold patients out
    properly and *still* expect to see inflation relative to site-disjoint;
    that gap is the finding.
    """
    frame = pd.DataFrame({"patient": patients.to_numpy()}, index=patients.index)
    fold = pd.Series(-1, index=patients.index, dtype=int)

    if stratify is not None:
        splitter = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=seed)
        iterator = splitter.split(frame, stratify.to_numpy(), groups=frame["patient"])
    else:
        # GroupKFold is deterministic and ignores random_state, so shuffle the
        # group order ourselves to make `seed` meaningful. The group order it
        # then applies is pinned (`stable_group_kfold`, A12): with one row per
        # patient every group ties, and the library's unstable sort made these
        # folds platform-dependent. Runs before 2026-09-24 used the library's.
        rng = np.random.default_rng(seed)
        uniq = frame["patient"].unique()
        perm = {p: i for i, p in enumerate(rng.permutation(uniq))}
        shuffled = frame["patient"].map(perm)
        iterator = stable_group_kfold(shuffled.to_numpy(), n_splits=n_folds)

    for k, (_, test_idx) in enumerate(iterator):
        fold.iloc[test_idx] = k

    if (fold < 0).any():
        raise RuntimeError("some rows were never assigned a fold")

    return SplitResult(
        fold=fold,
        scheme="random_patient",
        n_folds=n_folds,
        balance=_balance(fold, stratify),
    )


def preserved_site_split(
    patients: pd.Series,
    sites: pd.Series,
    *,
    n_folds: int = 5,
    outcome: pd.Series | None = None,
    seed: int = 0,
    exact: bool = False,
    jitter: float = 0.0,
) -> SplitResult:
    """Site-disjoint folds: no tissue source site appears in more than one fold.

    Sites are assigned whole. If `outcome` is given, assignment additionally
    balances the mean outcome across folds, which prevents the degenerate case
    where one fold happens to collect all the high-signal sites.

    Set `exact=True` to solve the assignment as a QP via cvxpy when available;
    otherwise a deterministic greedy heuristic is used (largest sites first,
    each placed in the fold that minimises a joint size+outcome imbalance).
    """
    if not patients.index.equals(sites.index):
        raise ValueError("patients and sites must share an index")

    per_site = pd.DataFrame({"site": sites, "patient": patients}).groupby("site").agg(
        n=("patient", "nunique")
    )
    if outcome is not None:
        per_site["outcome"] = (
            pd.DataFrame({"site": sites, "y": outcome}).groupby("site")["y"].mean()
        )

    n_sites = len(per_site)
    if n_sites < n_folds:
        raise ValueError(
            f"cannot build {n_folds} site-disjoint folds from {n_sites} sites. "
            "Reduce n_folds, or pool cohorts to increase site diversity."
        )

    if exact:
        assignment = _assign_sites_qp(per_site, n_folds)
    else:
        assignment = _assign_sites_greedy(per_site, n_folds, seed=seed, jitter=jitter)

    fold = sites.map(assignment).astype(int)

    # HARD GUARD. The outcome-balancing branch z-scores the size cost and the
    # outcome cost to unit variance independently, so a meaningless outcome
    # spread can outweigh a real size imbalance. On adverse draws that has
    # produced folds like [3, 213, 336, 0, 401] — including an EMPTY fold, which
    # makes a CV estimate silently undefined. Refuse rather than return it.
    counts = fold.value_counts().reindex(range(n_folds), fill_value=0)
    target = len(fold) / n_folds
    if (counts == 0).any() or (counts < 0.5 * target).any():
        raise RuntimeError(
            f"degenerate site-disjoint partition: fold sizes {counts.tolist()} "
            f"(target ~{target:.0f} each). This usually means a few sites hold "
            "most patients, or that outcome balancing overwhelmed size balancing. "
            "Retry with a different seed, drop `outcome=`, or reduce n_folds."
        )

    result = SplitResult(
        fold=fold,
        scheme="preserved_site",
        n_folds=n_folds,
        balance=_balance(fold, outcome),
    )
    _assert_site_disjoint(fold, sites)
    return result


def _assign_sites_greedy(
    per_site: pd.DataFrame, n_folds: int, *, seed: int = 0, jitter: float = 0.0
) -> dict[str, int]:
    """Largest-first greedy with joint size/outcome balancing.

    Classic LPT scheduling, extended so that when two folds are near-equal in
    size we prefer the one whose outcome mean moves least. Deterministic given
    `seed`, which breaks exact ties in the per-fold COST (see `rng.choice`
    below). It does NOT reorder equal-sized sites unless `jitter > 0` -- see the
    correction on that point at the sort.
    """
    rng = np.random.default_rng(seed)

    # Largest-first. `kind="stable"` is LOAD-BEARING, not tidiness -- this line
    # was audit item A9.
    #
    # `np.argsort` defaults to `kind="quicksort"` (introsort), which is not
    # stable, so exactly-tied keys come out in whatever order the sort's internal
    # tie-breaking produces. That is implementation-defined: numpy ships
    # different SIMD kernels for arm64 and x86_64 and they disagree. In NSCLC 51
    # of 68 sites sit in a tied size group, so macOS and HPC4 built DIFFERENT
    # cross-validation partitions from byte-identical inputs -- while preserving
    # fold sizes exactly, which is why `_assert_site_disjoint` and every other
    # guard passed and six sessions did not see it. Measured on two machines by
    # `scripts/24_split_determinism.py`; `scripts/26_sort_audit.py` sweeps for
    # the same defect elsewhere.
    #
    # An earlier comment here claimed ties were "ordered RANDOMLY per seed". They
    # were not, and could not be: with the default `jitter=0.0` the `noise` term
    # is IDENTICALLY ZERO, so `sizes + noise` is `sizes` exactly and the seed
    # never reaches the ordering at all. That is also the real explanation for
    # the "between-partition sd = 0.000" the old comment recorded as a puzzle.
    # With `jitter > 0` the noise does separate tied sites, no exact ties remain,
    # and stable and quicksort necessarily agree -- so this fix changes behaviour
    # only on the `jitter=0.0` path, which is the default and the frozen one.
    sizes = per_site["n"].to_numpy(dtype=float)
    noise = rng.uniform(0, jitter * max(sizes.max(), 1.0), size=len(sizes))
    order = per_site.index.to_numpy()[np.argsort(-(sizes + noise), kind="stable")]

    has_outcome = "outcome" in per_site.columns
    grand_mean = float(per_site["outcome"].mean()) if has_outcome else 0.0

    sizes = np.zeros(n_folds)
    sums = np.zeros(n_folds)
    counts = np.zeros(n_folds)
    assignment: dict[str, int] = {}

    for site in order:
        n = float(per_site.at[site, "n"])
        size_cost = sizes + n
        size_cost = (size_cost - size_cost.mean()) / (size_cost.std() + 1e-9)

        if has_outcome:
            y = float(per_site.at[site, "outcome"])
            new_means = (sums + y * n) / np.maximum(counts + n, 1e-9)
            out_cost = np.abs(new_means - grand_mean)
            out_cost = (out_cost - out_cost.mean()) / (out_cost.std() + 1e-9)
            cost = size_cost + out_cost
        else:
            cost = size_cost

        best = np.flatnonzero(cost == cost.min())
        choice = int(best[0]) if len(best) == 1 else int(rng.choice(best))

        assignment[site] = choice
        sizes[choice] += n
        counts[choice] += n
        if has_outcome:
            sums[choice] += float(per_site.at[site, "outcome"]) * n

    return assignment


def _assign_sites_qp(per_site: pd.DataFrame, n_folds: int) -> dict[str, int]:
    """Exact assignment via mixed-integer QP. Requires cvxpy; falls back if absent."""
    try:
        import cvxpy as cp
    except ImportError:  # pragma: no cover - optional dependency
        import warnings

        warnings.warn("cvxpy not installed; falling back to greedy assignment", stacklevel=2)
        return _assign_sites_greedy(per_site, n_folds)

    sites = per_site.index.to_numpy()
    n = per_site["n"].to_numpy(dtype=float)
    target = n.sum() / n_folds

    x = cp.Variable((len(sites), n_folds), boolean=True)
    fold_sizes = n @ x
    objective = cp.Minimize(cp.sum_squares(fold_sizes - target))
    constraints = [cp.sum(x, axis=1) == 1]

    problem = cp.Problem(objective, constraints)
    problem.solve()
    if x.value is None:
        import warnings

        warnings.warn("QP did not solve; falling back to greedy", stacklevel=2)
        return _assign_sites_greedy(per_site, n_folds)

    return {s: int(np.argmax(x.value[i])) for i, s in enumerate(sites)}


def _assert_site_disjoint(fold: pd.Series, sites: pd.Series) -> None:
    """Hard guarantee. A silent violation here would invalidate the primary result."""
    per_site_folds = pd.DataFrame({"fold": fold, "site": sites}).groupby("site")["fold"].nunique()
    leaky = per_site_folds[per_site_folds > 1]
    if len(leaky):
        raise RuntimeError(f"site(s) spanning multiple folds: {list(leaky.index)[:10]}")


def _balance(fold: pd.Series, outcome: pd.Series | None) -> pd.DataFrame:
    rows = fold.value_counts().sort_index().rename("n").to_frame()
    if outcome is not None:
        rows["outcome_mean"] = outcome.groupby(fold).mean()
        rows["outcome_std"] = outcome.groupby(fold).std()
    return rows


def fold_balance_report(result: SplitResult) -> str:
    """Human-readable balance check. Put this in the supplement."""
    b = result.balance
    lines = [
        f"scheme={result.scheme}  n_folds={result.n_folds}",
        f"fold sizes: min={b['n'].min()} max={b['n'].max()} "
        f"imbalance={b['n'].max() / max(b['n'].min(), 1):.2f}x",
    ]
    if "outcome_mean" in b.columns:
        spread = b["outcome_mean"].max() - b["outcome_mean"].min()
        lines.append(f"outcome mean across folds: spread={spread:.4f}")
    return "\n".join(lines)


def repeated_preserved_site_splits(
    patients: pd.Series,
    sites: pd.Series,
    *,
    n_folds: int = 5,
    repeats: int = 10,
    outcome: pd.Series | None = None,
    seed: int = 0,
) -> list[SplitResult]:
    """Several site-disjoint partitions, each with a different site ordering.

    Why this matters for the reported interval
    ------------------------------------------
    The greedy assignment is deterministic given a seed, so a single call hides
    how much the result depends on WHICH partition you happened to draw.
    Measured across 15 valid partitions of the NSCLC cohort: Delta_r mean 0.210
    with sd 0.019, against a patient-bootstrap half-width of 0.048. So roughly
    27% of the honest interval is missing when partition variance is ignored.

    Report the between-partition sd alongside the bootstrap CI, or combine them.
    Degenerate partitions are skipped (with the seed recorded) rather than
    silently retried, so the caller can see how often the geometry failed.
    """
    out: list[SplitResult] = []
    skipped: list[int] = []
    for r in range(repeats):
        try:
            out.append(
                preserved_site_split(
                    patients, sites, n_folds=n_folds, outcome=outcome,
                    seed=seed + r,
                    # r=0 reproduces the canonical deterministic partition; later
                    # repeats perturb the site ordering so the spread is real.
                    jitter=0.0 if r == 0 else 0.25,
                )
            )
        except (RuntimeError, ValueError):
            skipped.append(seed + r)
    if not out:
        raise RuntimeError(
            f"no valid site-disjoint partition in {repeats} attempts; "
            "the site-size distribution cannot support this n_folds"
        )
    if skipped:
        out[0].balance.attrs["skipped_seeds"] = skipped
    return out


def partition_variance(values: list[float]) -> dict[str, float]:
    """Between-partition spread of a statistic. Report next to the bootstrap CI."""
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if len(v) < 2:
        return {"n_partitions": len(v)}
    return {
        "n_partitions": int(len(v)),
        "mean": float(v.mean()),
        "sd": float(v.std(ddof=1)),
        "min": float(v.min()),
        "max": float(v.max()),
    }
