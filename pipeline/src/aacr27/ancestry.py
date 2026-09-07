"""The ancestry arm — does the image-to-signature relationship transfer equally?

WHAT THIS ARM CAN AND CANNOT SAY
================================
This is a **non-inferiority / equity** analysis, not a discovery analysis, and it
is scoped that way deliberately. 02-PROJECT-DECISION.md §9.4 cut it from the
abstract for five reasons that all still hold:

  1. it is fatal in NSCLC by the project's own numbers (ASIAN n≈20 against ~137
     needed), so it exists pan-cancer or not at all;
  2. ancestry is heavily confounded with tissue source site;
  3. it yields a non-inferiority statement, not a positive finding;
  4. it adds a fourth multiplicity family;
  5. it collides with the AACR CFA's terminology caution on race/ancestry.

It survives as **a supplementary table**. Nothing in this module should be
promoted to a headline, and every contrast is reported next to its minimum
detectable effect so that "we could not detect a difference of at least X" is
never silently rewritten as "there was no difference".

THE ESTIMAND
============
One model is trained on the whole cohort under preserved-site cross-validation.
Out-of-fold predictions are then correlated with the signature score **within**
each ancestry group, and each group is contrasted against EUR in Fisher z:

    delta_z(g) = z[r(pred, signature | ancestry = g)]
               - z[r(pred, signature | ancestry = EUR)]

That is the deployment question — a model fitted on a cohort that is 81% EUR,
does it work as well for everyone else? — rather than "is ancestry biologically
associated with the signature", which this design cannot answer.

TWO CONFOUNDS THAT MUST BE REPORTED, NOT ADJUSTED AWAY QUIETLY
==============================================================
**Cancer type.** The groups have very different disease mixes. Measured at
patient level in this cohort: ASIAN is 29.2% LIHC against EUR's 2.6%, an 11-fold
enrichment; AMR is 16.6% THCA against EUR's 5.9%. Signature predictability
varies by tumour type, so a crude between-group contrast is substantially a
between-disease contrast. `stratified_contrast()` therefore computes the
contrast WITHIN each cancer type and pools by inverse variance, and that
stratified number — not the crude one — is the reportable estimate. Both are
returned so the gap between them is visible.

**Tissue source site.** Ancestry and TSS are strongly associated, so a residual
difference may be a site effect. `cramers_v()` quantifies it. Note the bias
correction: with 619 TSS levels and ~7,000 patients the uncorrected statistic is
inflated by construction, so Bergsma's correction is applied and both values are
returned.

COUNTS MUST BE DERIVED AT PATIENT LEVEL
=======================================
An earlier revision quoted AFR 717 / EAS 535 / AMR 249. Those were **slide**
counts. Collapsed to patients this cohort gives AFR 665 / ASIAN 539 / **AMR 181**.
Quoting the slide numbers would overstate the smallest group by 38%.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

from .stats import Estimate, fisher_z, min_detectable_delta_z

REFERENCE_GROUP = "EUR"
ANALYSIS_GROUPS = ("EUR", "AFR", "ASIAN", "AMR")

# A Fisher-z correlation needs a usable n. Below this a cell contributes noise
# and an unstable 1/(n-3) weight, so it is dropped and the drop is reported.
MIN_CELL_N = 20


def cramers_v(x: pd.Series, y: pd.Series) -> dict[str, float]:
    """Association between two categorical variables, with and without correction.

    The uncorrected statistic is badly biased upward when either variable has
    many levels relative to n — which is exactly the ancestry-by-TSS case (619
    sites, ~7,000 patients). Bergsma's correction (J Korean Stat Soc 2013;42:323)
    removes the leading bias term. Report the corrected value; the raw one is
    returned only so a reader can see how large the correction was.
    """
    table = pd.crosstab(x, y)
    if table.size == 0 or table.shape[0] < 2 or table.shape[1] < 2:
        return {"v": np.nan, "v_uncorrected": np.nan, "n": int(table.values.sum())}

    chi2 = sps.chi2_contingency(table, correction=False)[0]
    n = int(table.values.sum())
    r, k = table.shape
    phi2 = chi2 / n

    v_raw = float(np.sqrt(phi2 / max(min(r - 1, k - 1), 1)))

    phi2_corr = max(0.0, phi2 - (r - 1) * (k - 1) / max(n - 1, 1))
    r_corr = r - (r - 1) ** 2 / max(n - 1, 1)
    k_corr = k - (k - 1) ** 2 / max(n - 1, 1)
    denom = max(min(r_corr - 1, k_corr - 1), 1e-12)

    return {
        "v": float(np.sqrt(phi2_corr / denom)),
        "v_uncorrected": v_raw,
        "n": n,
        "n_rows": r,
        "n_cols": k,
        "chi2": float(chi2),
    }


def cramers_v_calibrated(
    x: pd.Series, y: pd.Series, *, n_perm: int = 200, seed: int = 0
) -> dict[str, float]:
    """Cramér's V against a PERMUTATION null. Report `v_calibrated`.

    Bergsma's analytic correction removes the leading bias term but is badly
    behaved when the table is very sparse, and the ancestry-by-site table is
    exactly that: 4 groups x 608 sites over ~6,600 patients is 2.7 patients per
    cell. Measured on independent variables with 120 levels at n=400, the
    *corrected* statistic still returned 0.25.

    So the analytic value cannot be trusted on its own. Permuting one label
    destroys the association while preserving both margins and the sparsity, and
    the difference between observed and permuted V is the part that is real. This
    is the same device `decomposition.permutation_calibration` applies to
    site-R^2, and for the same reason.
    """
    obs = cramers_v(x, y)
    rng = np.random.default_rng(seed)
    xv = pd.Series(x).reset_index(drop=True)
    yv = pd.Series(y).reset_index(drop=True)

    draws = []
    for _ in range(n_perm):
        permuted = pd.Series(rng.permutation(xv.to_numpy()))
        v = cramers_v(permuted, yv)["v"]
        if np.isfinite(v):
            draws.append(v)

    if not draws:
        return {**obs, "v_calibrated": np.nan, "n_perm": 0}

    d = np.asarray(draws)
    return {
        **obs,
        "v_perm_median": float(np.median(d)),
        "v_perm_p95": float(np.percentile(d, 95)),
        "v_calibrated": float(obs["v"] - np.median(d)),
        "perm_p": float((1 + np.sum(d >= obs["v"])) / (len(d) + 1)),
        "n_perm": int(len(d)),
    }


def _safe_r(a: np.ndarray, b: np.ndarray) -> tuple[float, int]:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 4:
        return float("nan"), int(ok.sum())
    x, y = a[ok], b[ok]
    if x.std() == 0 or y.std() == 0:
        return float("nan"), int(ok.sum())
    return float(np.corrcoef(x, y)[0, 1]), int(ok.sum())


def group_correlations(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    groups: pd.Series,
    *,
    which: tuple[str, ...] = ANALYSIS_GROUPS,
    min_n: int = MIN_CELL_N,
) -> pd.DataFrame:
    """Within-group correlation between the signature and its OOF image prediction."""
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    g = pd.Series(groups).to_numpy()

    rows = []
    for name in which:
        mask = g == name
        if mask.sum() < min_n:
            rows.append({"group": name, "n": int(mask.sum()), "r": np.nan,
                         "z": np.nan, "se_z": np.nan, "below_min_n": True})
            continue
        r, n = _safe_r(y_true[mask], y_pred[mask])
        rows.append({
            "group": name, "n": n, "r": r,
            "z": float(fisher_z(r)) if np.isfinite(r) else np.nan,
            "se_z": 1.0 / np.sqrt(n - 3) if n > 3 else np.nan,
            "below_min_n": False,
        })
    return pd.DataFrame(rows)


def contrast_vs_reference(
    per_group: pd.DataFrame,
    *,
    reference: str = REFERENCE_GROUP,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """CRUDE contrast of each group against the reference, in Fisher z.

    Crude means unadjusted for cancer type, and in this cohort that is a real
    problem rather than a technicality — see the module docstring. Use
    `stratified_contrast` for the reportable number and read the two together.
    """
    ref = per_group[per_group["group"] == reference]
    if not len(ref) or not np.isfinite(ref["z"].iloc[0]):
        return pd.DataFrame()
    z_ref, n_ref = float(ref["z"].iloc[0]), int(ref["n"].iloc[0])

    crit = sps.norm.ppf(1 - alpha / 2)
    rows = []
    for _, row in per_group.iterrows():
        if row["group"] == reference or not np.isfinite(row["z"]):
            continue
        n_g = int(row["n"])
        se = np.sqrt(1.0 / max(n_g - 3, 1) + 1.0 / max(n_ref - 3, 1))
        d = float(row["z"]) - z_ref
        rows.append({
            "group": row["group"], "n": n_g, "n_reference": n_ref,
            "r": row["r"], "r_reference": float(ref["r"].iloc[0]),
            "delta_z": d, "delta_lo": d - crit * se, "delta_hi": d + crit * se,
            "se": se,
            # Report this next to every contrast. A null result means "we could
            # not detect a difference of at least MDE", never "no difference".
            "mde_80pct": min_detectable_delta_z(n_g, n_ref, alpha=alpha),
            "detectable": bool(abs(d) >= min_detectable_delta_z(n_g, n_ref, alpha=alpha)),
        })
    return pd.DataFrame(rows)


def stratified_contrast(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    groups: pd.Series,
    strata: pd.Series,
    *,
    reference: str = REFERENCE_GROUP,
    which: tuple[str, ...] = ANALYSIS_GROUPS,
    min_cell_n: int = MIN_CELL_N,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """THE REPORTABLE CONTRAST: computed within stratum, pooled by inverse variance.

    Fixed-effect meta-analysis over cancer types. Within stratum t,

        delta_z(g,t) = z[r(g,t)] - z[r(ref,t)],  var = 1/(n_gt-3) + 1/(n_reft-3)

    pooled with weights 1/var. Because every comparison happens inside one cancer
    type, the groups' different disease mixes cannot drive the result — which for
    this cohort is the difference between a real contrast and an artefact.

    Heterogeneity (Cochran's Q, I^2) is returned because a pooled estimate over
    strata that disagree is not meaningful, and `n_strata_used` / `n_strata_dropped`
    are returned because strata too small for a stable correlation are excluded
    and that exclusion must be visible.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    g = pd.Series(groups).reset_index(drop=True)
    s = pd.Series(strata).reset_index(drop=True).astype("string")

    crit = sps.norm.ppf(1 - alpha / 2)
    rows = []
    for name in which:
        if name == reference:
            continue
        deltas, weights, used, dropped, n_g_total, n_r_total = [], [], [], 0, 0, 0

        for stratum in s.dropna().unique():
            # `.fillna(False)` is load-bearing: `s` is a nullable "string" dtype,
            # so a missing stratum compares to pd.NA rather than False and the
            # resulting mask raises "boolean value of NA is ambiguous" on .sum().
            in_s = (s == stratum).fillna(False).to_numpy(dtype=bool)
            m_g = in_s & (g == name).to_numpy(dtype=bool)
            m_r = in_s & (g == reference).to_numpy(dtype=bool)
            if m_g.sum() < min_cell_n or m_r.sum() < min_cell_n:
                dropped += 1
                continue
            r_g, n_g = _safe_r(y_true[m_g], y_pred[m_g])
            r_r, n_r = _safe_r(y_true[m_r], y_pred[m_r])
            if not (np.isfinite(r_g) and np.isfinite(r_r)) or n_g <= 3 or n_r <= 3:
                dropped += 1
                continue
            var = 1.0 / (n_g - 3) + 1.0 / (n_r - 3)
            deltas.append(float(fisher_z(r_g) - fisher_z(r_r)))
            weights.append(1.0 / var)
            used.append(stratum)
            n_g_total += n_g
            n_r_total += n_r

        if not deltas:
            rows.append({"group": name, "delta_z": np.nan, "n_strata_used": 0,
                         "n_strata_dropped": dropped})
            continue

        d = np.asarray(deltas)
        w = np.asarray(weights)
        pooled = float((w * d).sum() / w.sum())
        se = float(np.sqrt(1.0 / w.sum()))
        Q = float((w * (d - pooled) ** 2).sum())
        df = len(d) - 1
        i2 = float(max(0.0, (Q - df) / Q) * 100) if Q > 0 and df > 0 else 0.0

        rows.append({
            "group": name,
            "delta_z": pooled, "delta_lo": pooled - crit * se, "delta_hi": pooled + crit * se,
            "se": se,
            "n_group": n_g_total, "n_reference": n_r_total,
            "n_strata_used": len(used), "n_strata_dropped": dropped,
            "strata": ", ".join(map(str, sorted(used))),
            "Q": Q, "df": df, "I2_pct": i2,
            "heterogeneous": bool(df > 0 and sps.chi2.sf(Q, df) < 0.10),
            "mde_80pct": float((sps.norm.ppf(1 - alpha / 2) + sps.norm.ppf(0.80)) * se),
        })
    return pd.DataFrame(rows)


def pooled_over_signatures(
    per_signature: pd.DataFrame,
    *,
    group: str,
    value_col: str = "delta_z",
    n_boot: int = 2000,
    seed: int = 0,
    alpha: float = 0.05,
) -> Estimate:
    """Mean contrast across the signature family, bootstrapped over signatures.

    The signature panel is correlated, so this interval is NOT a
    sqrt(m)-narrowing independent-sample mean; resampling the family is the
    honest cheap option here because — unlike the ISI — each signature's
    contrast already carries its own patient-level standard error and the
    quantity being pooled is a between-signature average.
    """
    vals = pd.to_numeric(
        per_signature.loc[per_signature["group"] == group, value_col], errors="coerce"
    ).dropna().to_numpy()
    if len(vals) < 2:
        return Estimate(float(vals[0]) if len(vals) else np.nan, np.nan, np.nan,
                        len(vals), f"mean-delta-z-{group}")

    rng = np.random.default_rng(seed)
    draws = np.array([
        np.mean(rng.choice(vals, size=len(vals), replace=True)) for _ in range(n_boot)
    ])
    lo, hi = np.percentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return Estimate(float(np.mean(vals)), float(lo), float(hi), len(vals),
                    f"mean-delta-z-{group}-signature-bootstrap")


def site_identifiability(
    ancestry: pd.Series,
    site: pd.Series,
    cancer_type: pd.Series | None = None,
    *,
    reference: str = REFERENCE_GROUP,
    which: tuple[str, ...] = ANALYSIS_GROUPS,
    min_each: int = 10,
) -> pd.DataFrame:
    """CAN the ancestry contrast be separated from site at all? Usually: no.

    A between-ancestry difference is only attributable to ancestry if the two
    groups are observed at the SAME sites. In TCGA they largely are not — sites
    are near-perfectly segregated by ancestry, because a source site is a
    hospital and hospitals serve local populations.

    Measured in this cohort at patient level:

        ASIAN    4 of 542 sites have >= 10 of both ASIAN and EUR
        AFR      6 of 576
        AMR      3 of 534

    That is the arm's real finding, and it is a stronger and more defensible
    statement than any group difference the crude or type-stratified contrast
    produces: **in TCGA, "ancestry" and "site" are not separable for this
    question.** Report this table before reporting any contrast, so nobody reads
    a between-group difference as a biological one.

    Returns one row per group with the number of supporting sites, the number of
    supporting (cancer type, site) cells, and how many patients those cover.
    """
    frame = pd.DataFrame({"anc": ancestry, "site": site}).reset_index(drop=True)
    if cancer_type is not None:
        frame["type"] = pd.Series(cancer_type).reset_index(drop=True)

    rows = []
    for name in which:
        if name == reference:
            continue
        sub = frame[frame["anc"].isin([name, reference])]
        by_site = sub.groupby("site")["anc"].agg(
            g=lambda s: int((s == name).sum()), e=lambda s: int((s == reference).sum())
        )
        ok_site = by_site[(by_site["g"] >= min_each) & (by_site["e"] >= min_each)]

        row = {
            "group": name,
            "n_group_total": int((frame["anc"] == name).sum()),
            "n_sites_total": int(by_site.shape[0]),
            "n_sites_supporting": int(len(ok_site)),
            "frac_sites_supporting": float(len(ok_site) / max(len(by_site), 1)),
            "n_group_in_supporting_sites": int(ok_site["g"].sum()),
            "n_reference_in_supporting_sites": int(ok_site["e"].sum()),
        }
        if cancer_type is not None:
            by_cell = sub.groupby(["type", "site"])["anc"].agg(
                g=lambda s: int((s == name).sum()), e=lambda s: int((s == reference).sum())
            )
            ok_cell = by_cell[(by_cell["g"] >= min_each) & (by_cell["e"] >= min_each)]
            row["n_type_site_cells_supporting"] = int(len(ok_cell))
            row["n_patients_in_supporting_cells"] = int(ok_cell["g"].sum() + ok_cell["e"].sum())
        rows.append(row)
    return pd.DataFrame(rows)


def composition_report(
    ancestry: pd.Series, cancer_type: pd.Series, site: pd.Series
) -> dict:
    """Everything a reader needs to judge whether a contrast is confounded."""
    frame = pd.DataFrame({"anc": ancestry, "type": cancer_type, "site": site})
    frame = frame[frame["anc"].isin(ANALYSIS_GROUPS)]
    counts = frame["anc"].value_counts()

    share = pd.crosstab(frame["anc"], frame["type"], normalize="index")
    ref = share.loc[REFERENCE_GROUP] if REFERENCE_GROUP in share.index else None
    enrich = None
    if ref is not None:
        enrich = (share / ref.replace(0, np.nan)).replace([np.inf, -np.inf], np.nan)

    return {
        "counts": counts.to_dict(),
        "n_total": int(len(frame)),
        "site_association": cramers_v(frame["anc"], frame["site"]),
        "type_association": cramers_v(frame["anc"], frame["type"]),
        "type_share": share,
        "type_enrichment_vs_reference": enrich,
    }
