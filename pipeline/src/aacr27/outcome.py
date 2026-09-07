"""The outcome arm — Venet's original design, run correctly.

WHY THIS ARM EXISTS
===================
Running the random-gene-set null on *expression* invited a fatal objection: the
image and the RNA are two measurements of one physical specimen, so a random gene
set being predictable from a picture of that specimen is the EXPECTED result, not
a discovery. Venet's finding was surprising precisely because survival is
**distal** to expression.

This module restores the distance. It asks Venet's actual question, of TME
signatures, in two layers nobody has run:

  LAYER 1 (RNA side, Venet's exact original)
      Does a curated TME signature predict outcome better than a size-matched
      random gene set? Venet found 28 of 47 published breast signatures did not.
      Never asked of the TME/immune signature family.

  LAYER 2 (image side, new)
      Does the IMAGE-PREDICTED signature predict outcome better than
      IMAGE-PREDICTED random sets? This is the quantity that matters if H&E-
      inferred TME scores are to guide anything.

Both layers are additionally run on global-axis-residualised scores, which is the
direct analogue of Venet's decisive control: adjusting for a proliferation
metagene (meta-PCNA) abrogated nearly all signature-outcome associations, and
unadjusted HR correlated with meta-PCNA at R^2 = 0.9.

SURVIVAL-ANALYSIS RULES ENCODED HERE
====================================
Two documented failure modes from the landscape audit are prevented structurally
rather than by convention:

1. **No dichotomisation.** Scores enter Cox models as continuous covariates.
   "Optimal" cutpoint searching drives type I error to ~40% rather than 5%
   (Altman, Lausen, Sauerbrei & Schumacher, JNCI 1994;86:829-835). There is no
   code path here that dichotomises a score.

2. **No C-hacking.** The reduction from model to scalar is fixed in advance:
   Harrell's C on the continuous linear predictor. Sonabend, Bender & Vollmer
   (Bioinformatics 2022;38:4178-4184) named the practice of picking whichever
   reduction or time point maximises C-index.

3. **Endpoint choice is explicit.** Liu et al., Cell 2018;173:400-416 (TCGA-CDR)
   recommends PFI for all but 4 of 33 TCGA types and OS for only 23. `Endpoint`
   defaults to PFI and REQUIRES you to acknowledge the choice, because using OS
   in a type where TCGA-CDR advises against it is a free win for a hostile
   reviewer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .globalaxis import excess_over_null, null_degeneracy

# TCGA-CDR default. See Liu et al., Cell 2018;173:400-416.
DEFAULT_ENDPOINT = "PFI"
VALID_ENDPOINTS = ("PFI", "OS", "DFI", "DSS")


@dataclass(frozen=True)
class Survival:
    """A right-censored endpoint, kept immutable so it cannot drift mid-analysis."""

    time: np.ndarray
    event: np.ndarray
    name: str = DEFAULT_ENDPOINT

    def __post_init__(self) -> None:
        if len(self.time) != len(self.event):
            raise ValueError("time and event must be the same length")
        if self.name not in VALID_ENDPOINTS:
            raise ValueError(f"endpoint must be one of {VALID_ENDPOINTS}, got {self.name!r}")

    @property
    def n_events(self) -> int:
        return int(np.nansum(self.event))

    @property
    def valid(self) -> np.ndarray:
        return np.isfinite(self.time) & np.isfinite(self.event) & (self.time > 0)

    def subset(self, mask: np.ndarray) -> "Survival":
        return Survival(self.time[mask], self.event[mask], self.name)


def from_frame(
    frame: pd.DataFrame,
    *,
    endpoint: str = DEFAULT_ENDPOINT,
    time_col: str | None = None,
    event_col: str | None = None,
) -> Survival:
    """Build a Survival from a TCGA-CDR style frame.

    TCGA-CDR column convention is `PFI` (event, 0/1) and `PFI.time` (days). The
    embeddings parquet ships TCGA-CDR merged, so these columns are usually
    already present.
    """
    event_col = event_col or endpoint
    time_col = time_col or f"{endpoint}.time"
    for col in (event_col, time_col):
        if col not in frame.columns:
            raise KeyError(
                f"{col!r} not found. Available survival-ish columns: "
                f"{[c for c in frame.columns if any(e in str(c).upper() for e in VALID_ENDPOINTS)][:12]}"
            )
    return Survival(
        time=pd.to_numeric(frame[time_col], errors="coerce").to_numpy(dtype=float),
        event=pd.to_numeric(frame[event_col], errors="coerce").to_numpy(dtype=float),
        name=endpoint,
    )


# ------------------------------------------------------------- concordance ---

def concordance_index(
    score: np.ndarray, surv: Survival, strata: np.ndarray | None = None
) -> float:
    """Harrell's C for a continuous risk score. Higher score = higher risk.

    A pair is comparable when the one with the shorter time had an event. Ties in
    score count a half. Returns NaN when there are no comparable pairs.

    `strata` RESTRICTS COMPARABLE PAIRS TO WITHIN-STRATUM, and pan-cancer that is
    not optional. Progression-free interval differs enormously between tumour
    types, so an unstratified C-index pooled across 31 diseases mostly measures
    "does this score know which cancer this is". Any score correlated with tissue
    of origin — which most expression signatures are — then looks prognostic.
    Passing cancer type here asks the question that was actually intended:
    does the score rank patients WITHIN a disease?

    Measured on pan-TCGA (n=7,168): unstratified, 94% of signatures beat their
    random-set null; stratified by cancer type, that collapses. The unstratified
    number was tissue of origin, not prognosis.
    """
    score = np.asarray(score, dtype=float)
    ok = surv.valid & np.isfinite(score)
    if strata is not None:
        strata = np.asarray(strata)
        ok &= pd.notna(strata)
    if ok.sum() < 10:
        return float("nan")

    t, e, s = surv.time[ok], surv.event[ok], score[ok]
    g = None if strata is None else np.asarray(strata)[ok]

    def block(tb: np.ndarray, eb: np.ndarray, sb: np.ndarray) -> tuple[float, float, float]:
        """Concordant / permissible / tied within one already-sorted block."""
        c = p = ti = 0.0
        for i in range(len(tb)):
            if eb[i] != 1:
                continue
            # Sorted by time, so everything strictly later is a contiguous tail.
            # Slicing that tail is what makes the stratified path affordable:
            # masking the whole cohort per event instead took the test suite from
            # 90 s to 17.5 min, and would have dominated the pan-cancer run.
            j = np.searchsorted(tb, tb[i], side="right")
            if j >= len(tb):
                continue
            others = sb[j:]
            p += len(others)
            c += float(np.sum(others < sb[i]))
            ti += float(np.sum(others == sb[i]))
        return c, p, ti

    concordant = permissible = tied = 0.0
    if g is None:
        order = np.argsort(t, kind="stable")
        concordant, permissible, tied = block(t[order], e[order], s[order])
    else:
        for key in pd.unique(pd.Series(g)):
            m = g == key
            if m.sum() < 2:
                continue
            tb, eb, sb = t[m], e[m], s[m]
            order = np.argsort(tb, kind="stable")
            c, p, ti = block(tb[order], eb[order], sb[order])
            concordant += c
            permissible += p
            tied += ti

    if permissible == 0:
        return float("nan")
    return float((concordant + 0.5 * tied) / permissible)


def cox_score(
    score: np.ndarray, surv: Survival, covariates: np.ndarray | None = None
) -> dict[str, float]:
    """Univariate (or covariate-adjusted) Cox on a CONTINUOUS score.

    The score is standardised so the hazard ratio is per standard deviation,
    which is comparable across signatures of different scales. There is no
    dichotomisation path.
    """
    from statsmodels.duration.hazard_regression import PHReg

    score = np.asarray(score, dtype=float)
    ok = surv.valid & np.isfinite(score)
    if covariates is not None:
        covariates = np.asarray(covariates, dtype=float)
        ok &= np.isfinite(covariates).all(axis=1)
    if ok.sum() < 20 or surv.subset(ok).n_events < 5:
        return {"hr": np.nan, "coef": np.nan, "se": np.nan, "p": np.nan,
                "c_index": np.nan, "n": int(ok.sum()), "n_events": 0}

    s = score[ok]
    sd = s.std()
    s = (s - s.mean()) / (sd if sd > 0 else 1.0)
    X = s.reshape(-1, 1) if covariates is None else np.hstack([s.reshape(-1, 1), covariates[ok]])

    try:
        fit = PHReg(surv.time[ok], X, status=surv.event[ok]).fit(disp=False)
        coef, se, p = float(fit.params[0]), float(fit.bse[0]), float(fit.pvalues[0])
    except Exception:  # noqa: BLE001 — a non-converging fit is a result, not a crash
        coef = se = p = np.nan

    return {
        "hr": float(np.exp(coef)) if np.isfinite(coef) else np.nan,
        "coef": coef,
        "se": se,
        "p": p,
        "c_index": concordance_index(s, surv.subset(ok)),
        "n": int(ok.sum()),
        "n_events": surv.subset(ok).n_events,
    }


# -------------------------------------------------------- the outcome null ---

def outcome_null(
    observed_scores: pd.Series,
    null_scores: pd.DataFrame,
    surv: Survival,
    *,
    metric: str = "c_index",
    alpha: float = 0.05,
    strata: np.ndarray | None = None,
) -> dict[str, object]:
    """Venet's question, asked correctly: does the curated signature beat random?

    Parameters
    ----------
    observed_scores
        The curated signature's score per patient (RNA-derived for layer 1, or
        the out-of-fold image prediction for layer 2).
    null_scores
        One column per random gene set, same index, scored identically.
    metric
        "c_index" (default, pre-specified) or "coef" (Cox log-hazard per SD).

    The reported excess is |C - 0.5| based, because a signature that predicts
    outcome in the WRONG direction is not evidence of specificity — random sets
    scatter symmetrically around 0.5 and an unsigned comparison would reward
    anti-predictive signatures.
    """
    if metric not in {"c_index", "coef"}:
        raise ValueError("metric must be 'c_index' or 'coef'")

    def evaluate(values: np.ndarray) -> float:
        if metric == "c_index":
            c = concordance_index(values, surv, strata=strata)
            # Fold to a directionless discrimination measure in [0.5, 1].
            return np.nan if not np.isfinite(c) else 0.5 + abs(c - 0.5)
        return abs(cox_score(values, surv)["coef"])

    obs = evaluate(observed_scores.to_numpy(dtype=float))
    draws = np.array(
        [evaluate(null_scores[c].to_numpy(dtype=float)) for c in null_scores.columns],
        dtype=float,
    )
    draws = draws[np.isfinite(draws)]

    if metric == "c_index":
        # Map C in [0.5, 1] onto a correlation-like scale so the Fisher-z excess
        # machinery applies: r = 2C - 1 is Somers' D, the standard rescaling.
        obs_r = 2 * obs - 1
        draws_r = 2 * draws - 1
    else:
        obs_r, draws_r = np.tanh(obs), np.tanh(draws)

    est = excess_over_null(obs_r, draws_r, n=int(surv.valid.sum()), alpha=alpha)

    # MINIMUM DETECTABLE EFFECT. A null outcome result is uninterpretable without
    # it: "0 of 32 signatures beat their null" can mean the signatures carry no
    # prognostic information, or simply that the cohort was too small to see the
    # information they do carry. Reporting the MDE forces the honest phrasing,
    # "we could not detect an excess of at least X".
    #
    # This matters concretely here. NSCLC (n=944, 328 events) could not detect an
    # excess below ~0.09 in Fisher z, i.e. a C-index advantage under ~0.045 —
    # while the pan-TCGA advantage is around 0.03. The NSCLC null and the
    # pan-cancer signal are therefore NOT in conflict; the smaller cohort simply
    # cannot see an effect that size.
    from scipy import stats as _sps
    se = ((est.hi - est.lo) / (2 * _sps.norm.ppf(1 - alpha / 2))
          if np.isfinite(est.hi) and np.isfinite(est.lo) else np.nan)
    mde = float((_sps.norm.ppf(1 - alpha / 2) + _sps.norm.ppf(0.80)) * se) if np.isfinite(se) else np.nan

    return {
        "mde_80pct_z": mde,
        # Somers' D = 2C - 1, so a difference of `mde` in Fisher z is roughly
        # mde/2 in C-index near C = 0.5. Quoted for readability only.
        "mde_80pct_c_index": float(np.tanh(mde) / 2) if np.isfinite(mde) else np.nan,
        "observed": float(obs),
        "null_mean": float(np.mean(draws)) if len(draws) else np.nan,
        "null_sd": float(np.std(draws, ddof=1)) if len(draws) > 1 else np.nan,
        "n_null": int(len(draws)),
        "excess_z": est.value,
        "excess_lo": est.lo,
        "excess_hi": est.hi,
        "beats_null": bool(est.lo > 0) if np.isfinite(est.lo) else False,
        "metric": metric,
        "stratified": strata is not None,
        "endpoint": surv.name,
        "n_events": surv.n_events,
        "degeneracy": null_degeneracy(draws_r),
    }


def run_outcome_arm(
    signature_scores: pd.DataFrame,
    null_families: dict[str, pd.DataFrame],
    surv: Survival,
    *,
    global_axis: pd.Series | None = None,
    within: pd.Series | None = None,
    alpha: float = 0.05,
    strata: np.ndarray | None = None,
) -> pd.DataFrame:
    """Both layers, raw and axis-residualised, for every signature.

    `null_families` maps signature name -> DataFrame of that signature's random
    sets, scored the same way. Pass image-predicted scores to get layer 2.

    The residualised rows are the decisive ones: they are the analogue of Venet's
    meta-PCNA adjustment, which abrogated nearly all published signature-outcome
    associations.
    """
    from .globalaxis import residualise

    rows = []
    for sig in signature_scores.columns:
        if sig not in null_families:
            continue
        obs = signature_scores[sig]
        nulls = null_families[sig]

        rows.append({"signature": sig, "adjustment": "raw",
                     **outcome_null(obs, nulls, surv, alpha=alpha, strata=strata)})

        if global_axis is not None:
            obs_r = residualise(obs, global_axis, within=within)
            nulls_r = residualise(nulls, global_axis, within=within)
            rows.append({"signature": sig, "adjustment": "axis_residualised",
                         **outcome_null(obs_r, nulls_r, surv, alpha=alpha,
                                        strata=strata)})

    out = pd.DataFrame(rows)
    if len(out):
        out = out.drop(columns=["degeneracy"], errors="ignore")
        out = out.sort_values(["signature", "adjustment"])
    return out


def endpoint_warning(cancer_types: pd.Series, endpoint: str) -> str | None:
    """Flag an endpoint choice a reviewer would challenge.

    TCGA-CDR (Liu et al., Cell 2018) recommends PFI for all but 4 of 33 types and
    OS for only 23. This does not hard-code the per-type table — look it up in
    TCGA-CDR Table 1 for the specific types in your cohort — but it does refuse to
    let OS pass silently.
    """
    if endpoint == "PFI":
        return None
    n_types = cancer_types.nunique()
    return (
        f"Endpoint is {endpoint}, not PFI, across {n_types} cancer type(s). "
        "TCGA-CDR recommends PFI for all but 4 of 33 TCGA types and OS for only 23. "
        "Justify this per type against TCGA-CDR Table 1 (Liu et al., Cell 2018;173:400-416) "
        "or a reviewer will do it for you."
    )
