"""Statistics: correlation CIs, paired deltas, FDR, and power.

Every number that goes on the abstract should come out of this module, so that
the confidence interval is computed the same way everywhere. Two rules encoded
here, both from the landscape audit's failure-mode catalogue:

1. Report the *difference* with a CI, not two separate point estimates. A
   reviewer comparing 0.61 and 0.58 without a CI on the gap has learned nothing
   (audit section 8.5, item 3).
2. Bootstrap at the PATIENT level. Resampling slides inflates precision because
   slides from one patient are not independent.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats as sps


@dataclass(frozen=True)
class Estimate:
    """A point estimate with a confidence interval."""

    value: float
    lo: float
    hi: float
    n: int
    method: str = ""

    def __str__(self) -> str:
        return f"{self.value:.4f} [{self.lo:.4f}, {self.hi:.4f}] (n={self.n})"

    @property
    def width(self) -> float:
        return self.hi - self.lo


def fisher_z(r: np.ndarray | float) -> np.ndarray | float:
    """Fisher z-transform. Clipped to avoid infinities at |r|=1."""
    r = np.clip(np.asarray(r, dtype=float), -0.999999, 0.999999)
    return np.arctanh(r)


def inv_fisher_z(z: np.ndarray | float) -> np.ndarray | float:
    return np.tanh(np.asarray(z, dtype=float))


def corr_ci(
    x: np.ndarray, y: np.ndarray, *, method: str = "pearson", alpha: float = 0.05
) -> Estimate:
    """Correlation with a Fisher-z confidence interval.

    Analytic, so it assumes independent observations. When observations are
    patients this is fine; when they are slides use `bootstrap_corr` instead.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 4:
        return Estimate(np.nan, np.nan, np.nan, n, f"{method}-fisher")

    r = sps.pearsonr(x, y)[0] if method == "pearson" else sps.spearmanr(x, y)[0]
    se = 1.0 / np.sqrt(n - 3)
    crit = sps.norm.ppf(1 - alpha / 2)
    z = fisher_z(r)
    return Estimate(float(r), float(inv_fisher_z(z - crit * se)),
                    float(inv_fisher_z(z + crit * se)), n, f"{method}-fisher")


def bootstrap_corr(
    x: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    *,
    method: str = "pearson",
    n_boot: int = 2000,
    alpha: float = 0.05,
    seed: int = 0,
) -> Estimate:
    """Cluster (patient-level) bootstrap of a correlation.

    `groups` is the patient id per row. We resample patients with replacement
    and take all of each sampled patient's rows, which propagates within-patient
    correlation into the interval instead of pretending it away.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    groups = np.asarray(groups)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, groups = x[ok], y[ok], groups[ok]

    uniq = np.unique(groups)
    index_by_group = {g: np.flatnonzero(groups == g) for g in uniq}
    rng = np.random.default_rng(seed)
    fn = (lambda a, b: sps.pearsonr(a, b)[0]) if method == "pearson" else (
        lambda a, b: sps.spearmanr(a, b)[0]
    )

    point = fn(x, y)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        picked = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([index_by_group[g] for g in picked])
        draws[i] = fn(x[idx], y[idx]) if len(idx) > 3 else np.nan

    lo, hi = np.nanpercentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return Estimate(float(point), float(lo), float(hi), len(uniq),
                    f"{method}-cluster-bootstrap")


def paired_delta(
    a: np.ndarray,
    b: np.ndarray,
    groups: np.ndarray,
    *,
    n_boot: int = 2000,
    alpha: float = 0.05,
    seed: int = 0,
) -> Estimate:
    """Bootstrap CI on a paired difference of means (a - b), clustered by group.

    Use for "random-split score minus site-disjoint score" where both scores are
    computed on the same units. This is the primary endpoint's estimator.
    """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    groups = np.asarray(groups)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b, groups = a[ok], b[ok], groups[ok]

    uniq = np.unique(groups)
    index_by_group = {g: np.flatnonzero(groups == g) for g in uniq}
    rng = np.random.default_rng(seed)

    point = float(np.mean(a - b))
    draws = np.empty(n_boot)
    for i in range(n_boot):
        picked = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([index_by_group[g] for g in picked])
        draws[i] = np.mean(a[idx] - b[idx])

    lo, hi = np.percentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return Estimate(point, float(lo), float(hi), len(uniq), "paired-cluster-bootstrap")


def delta_r_ci(
    r1: float, n1: int, r2: float, n2: int, *, alpha: float = 0.05
) -> Estimate:
    """CI on the difference between two INDEPENDENT correlations, in Fisher z.

    For the ancestry contrast (e.g. r in EUR vs r in AFR), where the two groups
    are disjoint sets of patients.
    """
    z1, z2 = fisher_z(r1), fisher_z(r2)
    se = np.sqrt(1.0 / max(n1 - 3, 1) + 1.0 / max(n2 - 3, 1))
    crit = sps.norm.ppf(1 - alpha / 2)
    d = float(z1 - z2)
    return Estimate(d, float(d - crit * se), float(d + crit * se),
                    min(n1, n2), "fisher-z-difference")


def min_detectable_delta_z(n1: int, n2: int, *, alpha: float = 0.05,
                           power: float = 0.80) -> float:
    """Minimum detectable difference in Fisher z between two groups.

    Report this next to every subgroup result. If the observed gap is smaller
    than the MDE, the honest statement is "we could not detect a difference of
    at least X", not "there was no difference".
    """
    se = np.sqrt(1.0 / max(n1 - 3, 1) + 1.0 / max(n2 - 3, 1))
    return float((sps.norm.ppf(1 - alpha / 2) + sps.norm.ppf(power)) * se)


def bh_fdr(pvals: np.ndarray, alpha: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Benjamini-Hochberg. Returns (rejected, qvalues)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    out_q = np.empty(n)
    out_q[order] = q
    return out_q <= alpha, out_q


def effective_tests(corr_matrix: np.ndarray) -> float:
    """Effective number of independent tests (Li & Ji, Heredity 2005;95:221-227).

    Meff = sum_i [ I(|lambda_i| >= 1) + (|lambda_i| - floor(|lambda_i|)) ]

    DESCRIPTIVE ONLY — do NOT use as a multiplicity correction.
    Meff is a family-wise-error device: you divide alpha by it. This project
    controls FDR with Benjamini-Hochberg, which is valid under positive
    dependence and already accounts for correlation. Dividing alpha by Meff *and*
    applying BH is a category error; the two do not compose. Report Meff purely
    as a description of how redundant the signature panel is.

    Three corrections against the previous implementation:
      - I(|lambda| >= 1), not (lambda > 1). On an exact identity the old form
        returned 0.0 instead of k, because every eigenvalue is exactly 1 and
        `> 1` is false while the fractional part is 0.
      - NaN input raises instead of being silently coerced to r = 0, which
        turned a missing correlation into an assertion of independence.
      - Non-PSD input is projected to the nearest PSD matrix. Pairwise-complete
        correlations across signatures with different missingness are routinely
        non-PSD, which produces negative eigenvalues and a nonsense Meff.
    """
    C = np.asarray(corr_matrix, dtype=float)
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise ValueError(f"corr_matrix must be square, got shape {C.shape}")
    if not np.isfinite(C).all():
        raise ValueError(
            "corr_matrix contains NaN/inf. Compute it with complete cases, or "
            "impute deliberately — silently treating a missing correlation as 0 "
            "asserts independence that was never observed."
        )

    C = (C + C.T) / 2.0                       # enforce exact symmetry
    eig = np.linalg.eigvalsh(C)
    if (eig < -1e-8).any():                   # project to nearest PSD
        eig = np.clip(eig, 0.0, None)

    a = np.abs(eig)
    return float(np.sum((a >= 1.0).astype(float) + (a - np.floor(a))))


def rotation_null(
    null_r: np.ndarray, observed_family_mean: float
) -> dict[str, float]:
    """Family-level p-value that respects between-signature correlation.

    The per-signature nulls are NOT independent — every signature is scored
    against the same images and the same global axis, so a naive sigma/sqrt(m)
    interval around the family mean is far too narrow. In simulation the naive
    formula claimed 50 +/- 10.5 where the true 95% band was 13.9-89.0.

    Construction: `null_r` is (n_signatures x B). Each column j is treated as a
    synthetic family drawn under the null with the SAME cross-signature
    correlation structure as the observed family, because it reuses the same
    fits. The family mean of column j is one draw of the null family statistic.

        p = (1 + #{draws >= observed}) / (B + 1)

    The +1 numerator and denominator make it a valid (non-zero) permutation
    p-value; never report p = 0 from a finite null.
    """
    R = np.asarray(null_r, dtype=float)
    if R.ndim != 2:
        raise ValueError(f"null_r must be 2-D (n_signatures x B), got {R.shape}")

    family_means = np.nanmean(R, axis=0)
    family_means = family_means[np.isfinite(family_means)]
    B = len(family_means)
    if B < 2:
        return {"p": np.nan, "B": B}

    n_ge = int(np.sum(family_means >= observed_family_mean))
    return {
        "observed": float(observed_family_mean),
        "null_family_mean": float(np.mean(family_means)),
        "null_family_sd": float(np.std(family_means, ddof=1)),
        "ci_lo": float(np.percentile(family_means, 2.5)),
        "ci_hi": float(np.percentile(family_means, 97.5)),
        "p": float((1 + n_ge) / (B + 1)),
        "B": B,
    }


def summarise(estimates: dict[str, Estimate]) -> pd.DataFrame:
    """Tidy table of named estimates, ready for a supplement."""
    return pd.DataFrame(
        [
            {"name": k, "estimate": e.value, "ci_lo": e.lo, "ci_hi": e.hi,
             "ci_width": e.width, "n": e.n, "method": e.method}
            for k, e in estimates.items()
        ]
    )
