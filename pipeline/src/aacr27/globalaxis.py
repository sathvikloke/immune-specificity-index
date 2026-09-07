"""The global expression axis — and the corrected estimand built on it.

WHY THIS MODULE EXISTS
======================
The first version of this project's primary endpoint was "where does the curated
signature fall in its own random-gene-set distribution". A red-team lens showed
that endpoint cannot support the conclusion it was written to support, for two
reasons that are worth stating in full because they shape everything here.

**1. The inference did not discriminate.** Two hypotheses predict the identical
observation:

  H1  the image model learned nothing immune-specific;
  H2  bulk expression has a dominant global axis — purity, stromal fraction,
      proliferation, tissue of origin — that essentially every gene set of
      appreciable size loads onto, and the image reads that axis faithfully.

H2 is a *positive* result about H&E. Comparing a curated signature to a random
set does not separate them, because both load on the same axis. The separation
requires conditioning the SIGNATURE on the global axis (not, as the old Control B
did, conditioning the image PREDICTION on purity) and asking whether the residual
is still image-predictable above a residualised null.

**2. The null was near-degenerate.** Under mean-z scoring, two *independent*
random gene sets of size k with mean pairwise gene-gene correlation rho have
score correlation

    k * rho / (1 + (k - 1) * rho)

which at k=29, rho=0.10 is **0.763**, and at k=100, rho=0.10 is **0.917**. At the
sizes real TME signatures use, every random draw measures the same latent axis,
so the null distribution of r collapses toward a point. A percentile computed on
a near-degenerate null is a hair trigger: an excess of +0.01 lands at the 99th
percentile and -0.01 at the 1st. Percentiles are therefore replaced throughout by
**excess over the null mean in Fisher z, with a confidence interval**, and
`null_degeneracy()` is reported alongside as a mandatory diagnostic.

THE CORRECTED PRIMARY ENDPOINT
==============================
    immune-specific excess
      = z[ r(image_pred, signature | global axis) ]
      - mean_j z[ r(image_pred, random_set_j | global axis) ]

with a patient-clustered bootstrap CI. Pre-register the axis definition before
running: `method="pc1"` (PC1 of within-cancer-type log expression) or
`method="meanz"` (mean z-score over all expressed genes). They are not
interchangeable and the choice must be fixed in advance.

Interpretation, both branches:
  excess CI excludes 0   ->  immune-specific signal exists beyond composition.
                             Honest headline: "H&E reads composition, AND carries
                             an immune-specific residual of size X."
  excess CI covers 0     ->  what looked immune-specific is the global axis.
                             Honest headline names the axis rather than claiming
                             the model learned nothing.

Both are publishable. Neither requires a clinical claim the data cannot support.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import statsmodels.api as sm

from .stats import Estimate, fisher_z


@dataclass(frozen=True)
class GlobalAxis:
    """The dominant expression axis, plus what it is made of."""

    values: pd.Series
    method: str
    variance_explained: float | None
    within: str | None
    top_loadings: pd.Series | None = None

    def __str__(self) -> str:
        ve = "n/a" if self.variance_explained is None else f"{self.variance_explained:.1%}"
        return (
            f"GlobalAxis(method={self.method}, within={self.within}, "
            f"variance_explained={ve}, n={len(self.values)})"
        )


def compute_global_axis(
    expr: pd.DataFrame,
    *,
    method: str = "pc1",
    cancer_type: pd.Series | None = None,
    n_top_loadings: int = 30,
) -> GlobalAxis:
    """The dominant axis of variation in bulk expression.

    Parameters
    ----------
    expr
        samples x genes, log-transformed.
    method
        "pc1"   first principal component (captures the dominant axis directly)
        "meanz" mean z-score across all expressed genes (simpler, assumption-free)

        THESE ARE NOT INTERCHANGEABLE, AND "meanz" IS NOT A FREE SENSITIVITY
        ANALYSIS. Measured on NSCLC (n=944, 41,046 genes), residualising on each
        and then re-measuring the reliability of random 200-gene sets:

            axis                     r(axis, gene mean)   residual rho   alpha
            global PC1                      0.441            +0.031      0.866
            within-type PC1 (registered)    0.808            +0.017      0.749
            meanz                           1.000            +0.005      0.025

        "meanz" IS the gene-wise mean by construction, so residualising it forces
        the average residual gene-gene covariance to nearly zero. Random sets
        then have essentially no internal consistency, Spearman's correction
        divides by sqrt(0.025), and the disattenuated estimand becomes undefined
        for most draws. Switching axis method after seeing results would
        therefore not be a robustness check — it would change whether the primary
        endpoint exists at all. Fix it in advance and report the choice.
    cancer_type
        If given, the axis is computed WITHIN each cancer type and concatenated.
        Strongly recommended pan-cancer: otherwise PC1 is tissue of origin and
        every downstream residual is really a within-vs-between-tissue contrast.

    Sign convention: PC1's sign is arbitrary, so it is fixed so that the axis
    correlates positively with mean expression. Without this the residuals are
    unchanged but the loadings flip between runs, which makes reporting confusing.
    """
    if method not in {"pc1", "meanz"}:
        raise ValueError("method must be 'pc1' or 'meanz'")

    if cancer_type is not None:
        parts, ve = [], []
        for group, idx in cancer_type.groupby(cancer_type).groups.items():
            sub = expr.loc[expr.index.intersection(idx)]
            if len(sub) < 10:
                continue
            axis = compute_global_axis(sub, method=method, n_top_loadings=0)
            parts.append(axis.values)
            if axis.variance_explained is not None:
                ve.append(axis.variance_explained)
        if not parts:
            raise ValueError("no cancer type had >=10 samples")
        return GlobalAxis(
            values=pd.concat(parts).reindex(expr.index),
            method=method,
            variance_explained=float(np.mean(ve)) if ve else None,
            within="cancer_type",
        )

    z = (expr - expr.mean(axis=0)) / expr.std(axis=0).replace(0, np.nan)
    z = z.dropna(axis=1, how="all").fillna(0.0)

    if method == "meanz":
        return GlobalAxis(
            values=z.mean(axis=1), method="meanz", variance_explained=None, within=None
        )

    matrix = z.to_numpy(dtype=float)
    matrix -= matrix.mean(axis=0)
    # Economy SVD: we only need the leading component.
    U, S, Vt = np.linalg.svd(matrix, full_matrices=False)
    pc1 = U[:, 0] * S[0]
    ve = float(S[0] ** 2 / np.sum(S**2))

    loadings = pd.Series(Vt[0], index=z.columns)
    if np.corrcoef(pc1, matrix.mean(axis=1))[0, 1] < 0:
        pc1, loadings = -pc1, -loadings

    top = (
        loadings.reindex(loadings.abs().sort_values(ascending=False).index)
        .head(n_top_loadings)
        if n_top_loadings
        else None
    )
    return GlobalAxis(
        values=pd.Series(pc1, index=z.index),
        method="pc1",
        variance_explained=ve,
        within=None,
        top_loadings=top,
    )


def residualise(
    scores: pd.DataFrame | pd.Series,
    axis: GlobalAxis | pd.Series,
    *,
    within: pd.Series | None = None,
) -> pd.DataFrame | pd.Series:
    """Remove the global axis from signature scores by OLS.

    `within` optionally fits the removal separately per group (e.g. cancer type),
    which matters pan-cancer because the axis-to-signature slope differs by tissue
    and a single pooled slope would leave tissue structure in the residual.
    """
    axis_values = axis.values if isinstance(axis, GlobalAxis) else axis
    was_series = isinstance(scores, pd.Series)
    frame = scores.to_frame() if was_series else scores

    if within is not None:
        parts = []
        for _, idx in within.groupby(within).groups.items():
            common = frame.index.intersection(idx)
            if len(common) < 5:
                continue
            parts.append(_residualise_block(frame.loc[common], axis_values.reindex(common)))
        out = pd.concat(parts).reindex(frame.index)
    else:
        out = _residualise_block(frame, axis_values.reindex(frame.index))

    return out.iloc[:, 0] if was_series else out


def _residualise_block(frame: pd.DataFrame, axis: pd.Series) -> pd.DataFrame:
    a = pd.to_numeric(axis, errors="coerce").to_numpy(dtype=float)
    design = sm.add_constant(a.reshape(-1, 1), has_constant="add")
    out = {}
    for col in frame.columns:
        y = pd.to_numeric(frame[col], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(y) & np.isfinite(a)
        resid = np.full(len(y), np.nan)
        if ok.sum() >= 5:
            fit = sm.OLS(y[ok], design[ok]).fit()
            resid[ok] = y[ok] - fit.predict(design[ok])
        out[col] = resid
    return pd.DataFrame(out, index=frame.index)


def residualise_matrix(
    Y: np.ndarray,
    axis: np.ndarray,
    *,
    groups: np.ndarray | None = None,
) -> np.ndarray:
    """Group-wise OLS residuals of EVERY column of `Y` on [1, axis], vectorised.

    Numerically equivalent to calling `residualise()` column by column (verified
    to ~1e-12 in the test suite), but solves one least-squares system per group
    for all columns at once instead of one statsmodels fit per column.

    This is not only a speed convenience. Two places need it for correctness:

    - The ISI null residualises a 1,000-column random-set frame per signature.
      Column-by-column that is 16,000 OLS fits per run.
    - Reliability of a residualised composite needs the residual variance of
      every GENE (41,046 columns here), which is unaffordable any other way.
      See `signatures.alpha_from_score_and_item_variances`.

    Rows whose axis value is non-finite, and columns with fewer than 3 usable
    rows in a group, are returned as NaN rather than silently zero-filled.
    """
    Y = np.asarray(Y, dtype=float)
    a = np.asarray(axis, dtype=float)
    if Y.ndim != 2:
        raise ValueError(f"Y must be 2-D (n_samples x n_targets), got {Y.shape}")
    if len(a) != len(Y):
        raise ValueError(f"axis has {len(a)} rows but Y has {len(Y)}")

    out = np.full(Y.shape, np.nan)
    g = np.zeros(len(Y), dtype=int) if groups is None else np.asarray(groups)

    for key in pd.unique(pd.Series(g)):
        m = (g == key) & np.isfinite(a)
        if m.sum() < 3:
            continue
        A = np.column_stack([np.ones(m.sum()), a[m]])
        block = Y[m]
        finite = np.isfinite(block)
        full = finite.all(axis=0)

        if full.any():
            coef, *_ = np.linalg.lstsq(A, block[:, full], rcond=None)
            out[np.ix_(m, full)] = block[:, full] - A @ coef

        # Columns with missing values are fitted on their own complete cases; a
        # single NaN must not drop the column, and must not contaminate the
        # shared solve above.
        for j in np.where(~full)[0]:
            ok = finite[:, j]
            if ok.sum() < 3:
                continue
            Aj = A[ok]
            cj, *_ = np.linalg.lstsq(Aj, block[ok, j], rcond=None)
            idx = np.where(m)[0][ok]
            out[idx, j] = block[ok, j] - Aj @ cj

    return out


def residual_item_variances(
    items: pd.DataFrame,
    axis: pd.Series,
    *,
    within: pd.Series | None = None,
) -> pd.Series:
    """Variance of every column of `items` AFTER residualising on the axis.

    `items` is the z-scored gene matrix. The returned per-gene variances are the
    `sum_j var(item_j)` term of the general-form Cronbach alpha, which is what
    makes reliability of a RESIDUALISED composite computable in O(k) per gene
    set instead of O(k^2).
    """
    R = residualise_matrix(
        items.to_numpy(dtype=float),
        pd.to_numeric(axis, errors="coerce").to_numpy(dtype=float),
        groups=None if within is None else within.to_numpy(),
    )
    # A column that residualised to all-NaN has no variance to report. Return
    # NaN for it explicitly rather than letting nanvar warn and emit NaN anyway.
    usable = np.isfinite(R).sum(axis=0) >= 2
    out = np.full(R.shape[1], np.nan)
    if usable.any():
        out[usable] = np.nanvar(R[:, usable], axis=0, ddof=1)
    return pd.Series(out, index=items.columns)


def disattenuate(r: float | np.ndarray, reliability: float | np.ndarray) -> float | np.ndarray:
    """Spearman's correction for attenuation: r_true = r_obs / sqrt(reliability).

    Only the TARGET side is corrected here; the predictor (the image model) is
    common to every comparison, so its reliability cancels in the contrast.

    Two guards, both of which return NaN rather than a misleading number:

    - reliability <= 0.01: a gene set with no internal consistency is not a scale
      and its correlation cannot be disattenuated.
    - |r / sqrt(rel)| >= 1: sampling error in the reliability estimate has pushed
      the corrected value outside the valid correlation range. Clipping to 1
      would be worse than dropping it, because Fisher z of 1 is unbounded and a
      handful of such draws will dominate a null mean. Callers MUST report how
      many draws this removes — see `excess_over_null`, which does.
    """
    rel = np.asarray(reliability, dtype=float)
    rel = np.where(rel > 0.01, rel, np.nan)
    out = np.asarray(r, dtype=float) / np.sqrt(rel)
    return np.where(np.abs(out) < 1.0, out, np.nan)


def excess_over_null(
    r_observed: float,
    r_null: np.ndarray,
    *,
    n: int,
    alpha: float = 0.05,
    reliability_observed: float | None = None,
    reliability_null: np.ndarray | float | None = None,
) -> Estimate:
    """THE CORRECTED ESTIMAND: excess over the null mean, in Fisher z.

    Replaces the percentile, which is uninterpretable on a near-degenerate null.
    The returned value is a difference of Fisher-z correlations, so it has a real
    effect-size meaning and a CI that does not collapse when the null tightens.

    The interval combines two sources of uncertainty: sampling error in the
    observed correlation (1/sqrt(n-3)) and the spread of the null distribution
    itself. Report `null_degeneracy()` alongside — if the null SD is tiny, a
    large excess in z may still be a trivial excess in r.
    """
    r_null = np.asarray(r_null, dtype=float)

    # MANDATORY when curated and random sets differ in internal consistency.
    # Without this the curated set wins on measurement precision alone; in
    # simulation the uncorrected contrast had a 100% false-positive rate.
    method = "excess-over-null-z"
    if reliability_observed is not None and reliability_null is not None:
        n_before = int(np.isfinite(r_null).sum())
        r_observed = float(disattenuate(r_observed, reliability_observed))
        r_null = np.asarray(disattenuate(r_null, reliability_null), dtype=float)
        n_dropped = n_before - int(np.isfinite(r_null).sum())
        method = "excess-over-null-z-disattenuated"
        if n_dropped:
            # NO SILENT CAPS: a dropped draw is a draw whose reliability estimate
            # was too noisy to correct. Surfacing the count lets a reader judge
            # whether the null is trustworthy. A high drop rate means the gene
            # sets are too small or the sample too few to estimate alpha stably.
            method += f" [{n_dropped}/{n_before} null draws dropped: |r/sqrt(alpha)|>=1]"

    r_null = r_null[np.isfinite(r_null)]
    if len(r_null) < 5 or not np.isfinite(r_observed):
        return Estimate(np.nan, np.nan, np.nan, n, method)

    z_obs = float(fisher_z(r_observed))
    z_null = fisher_z(r_null)
    excess = z_obs - float(np.mean(z_null))

    se_obs = 1.0 / np.sqrt(max(n - 3, 1))
    se_null = float(np.std(z_null, ddof=1)) / np.sqrt(len(z_null))
    se = np.sqrt(se_obs**2 + se_null**2)

    from scipy import stats as sps

    crit = sps.norm.ppf(1 - alpha / 2)
    return Estimate(excess, excess - crit * se, excess + crit * se, n, method)


def null_degeneracy(r_null: np.ndarray) -> dict[str, float]:
    """Diagnostic: how concentrated is the null? Report this every time.

    Reports the EMPIRICAL spread of the null only. The theoretical concentration
    k*rho/(1+(k-1)*rho) is not recoverable from scores alone — use
    `expected_random_set_correlation()` for that. A null SD below ~0.02 in r means
    a percentile would have been a hair trigger and the excess-in-z estimand is
    doing real work.
    """
    r = np.asarray(r_null, dtype=float)
    r = r[np.isfinite(r)]
    if len(r) < 3:
        return {"n": len(r)}
    z = fisher_z(r)
    return {
        "n": float(len(r)),
        "r_mean": float(np.mean(r)),
        "r_sd": float(np.std(r, ddof=1)),
        "r_min": float(np.min(r)),
        "r_max": float(np.max(r)),
        "r_range": float(np.max(r) - np.min(r)),
        "z_sd": float(np.std(z, ddof=1)),
        "degenerate": float(np.std(r, ddof=1) < 0.02),
    }


def expected_random_set_correlation(k: int, rho: float) -> float:
    """k*rho / (1 + (k-1)*rho) — correlation between two independent size-k sets.

    Use this to show a reviewer, before any data, that the null MUST be
    concentrated at realistic signature sizes. At k=29, rho=0.10 -> 0.763.
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    return float(k * rho / (1 + (k - 1) * rho))


def three_way_report(
    image_pred: pd.Series,
    signature: pd.Series,
    axis: GlobalAxis | pd.Series,
    *,
    within: pd.Series | None = None,
) -> dict[str, float]:
    """The arm that separates H1 from H2. Run this for every signature.

    Returns:
      r_raw        image vs the raw signature       (what the literature reports)
      r_axis       image vs the global axis alone   (how much is composition)
      r_residual   image vs the axis-residualised signature  (what is left)

    If r_residual is materially above the residualised random-set null, immune-
    specific signal exists. If it collapses while r_axis stays high, the model was
    reading composition — which is a finding about what H&E can and cannot do,
    not a claim that the model is worthless.
    """
    axis_values = axis.values if isinstance(axis, GlobalAxis) else axis
    idx = image_pred.index.intersection(signature.index).intersection(axis_values.index)

    p = pd.to_numeric(image_pred.reindex(idx), errors="coerce")
    s = pd.to_numeric(signature.reindex(idx), errors="coerce")
    a = pd.to_numeric(axis_values.reindex(idx), errors="coerce")

    resid = residualise(s, a, within=None if within is None else within.reindex(idx))

    def safe_corr(x: pd.Series, y: pd.Series) -> float:
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 10:
            return float("nan")
        return float(np.corrcoef(x[ok], y[ok])[0, 1])

    return {
        "r_raw": safe_corr(p, s),
        "r_axis": safe_corr(p, a),
        "r_residual": safe_corr(p, resid),
        "n": int(len(idx)),
    }
