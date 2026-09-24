"""Variance decomposition — how much of an image-predicted signature is purity?

The second scientific core of the project, alongside the Venet null.

The claim under test: H&E models that appear to predict immune-microenvironment
signatures may largely be predicting tumour purity and stromal fraction, both of
which are directly visible in histology and both of which drive bulk-derived
immune scores mechanically (a sample with less tumour has proportionally more
immune cells, so the immune signature rises without any change in immune
biology).

Precedent for the purity leg: H&E-predicted ABSOLUTE purity reaches Spearman
0.418-0.655 across TCGA (Oner et al., Patterns 2021). If the image model
predicts purity that well, and purity explains most of a signature's variance,
then the signature's apparent predictability is largely purity in disguise.

Method
------
For each signature we fit nested linear models and report incremental R^2:

    M0: intercept only
    M1: + tumour purity
    M2: + cancer type
    M3: + tissue source site
    M4: + stage
    Mfull: + image prediction

The quantity of interest is the *incremental* R^2 of the image prediction over
M4 — i.e. what the image adds beyond purity, type, site and stage. If that is
near zero, the headline is written.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import statsmodels.api as sm


@dataclass(frozen=True)
class DecompositionResult:
    signature: str
    table: pd.DataFrame  # one row per nested model
    incremental_image_r2: float
    n: int

    @property
    def summary_line(self) -> str:
        return (
            f"{self.signature}: image adds {self.incremental_image_r2:.4f} R^2 "
            f"beyond purity+type+site+stage (n={self.n})"
        )


def _design(
    frame: pd.DataFrame, terms: list[str], categorical: set[str]
) -> np.ndarray | None:
    """Build a design matrix, dummy-coding categoricals and dropping constants."""
    if not terms:
        return None
    blocks = []
    for t in terms:
        if t not in frame.columns:
            continue
        col = frame[t]
        if t in categorical:
            # dummy_na=True: a missing site/stage becomes its own indicator
            # rather than being silently absorbed into the reference level,
            # which would attribute its effect to whichever category sorted first.
            dummies = pd.get_dummies(
                col.astype("string"), drop_first=True, dummy_na=True, dtype=float
            )
            if dummies.shape[1]:
                blocks.append(dummies.to_numpy())
        else:
            blocks.append(pd.to_numeric(col, errors="coerce").to_numpy().reshape(-1, 1))
    if not blocks:
        return None
    X = np.hstack(blocks)
    keep = np.nanstd(X, axis=0) > 0
    return X[:, keep] if keep.any() else None


def _r2(
    y: np.ndarray, X: np.ndarray | None, mask: np.ndarray | None = None
) -> tuple[float, float, int, int]:
    """OLS on a FIXED sample. Returns (r2, adj_r2, n_params, n_used).

    `mask` must be the SAME complete-case mask for every rung of the ladder.
    Previously each rung computed its own mask, so with ~300 of 956 purity values
    missing, M0 was fitted on n=956 and M1..M_full on n=656 — and the reported
    "incremental" R^2 was a difference between models fitted on different
    samples, which is not an increment of anything.

    Adjusted R^2 is returned alongside because site is dummy-coded and TCGA has
    hundreds of source sites: on the real parquet, 641 TSS over 8,535 patients
    inflates unadjusted R^2 by roughly p/(n-1) by construction. That inflation
    lands on the covariate rungs and therefore DEFLATES incremental_image_r2 —
    a bias pointing toward this project's hoped-for conclusion, which is exactly
    the kind of thing a reviewer looks for.
    """
    if mask is None:
        mask = np.isfinite(y)
    if X is None:
        n = int(mask.sum())
        return 0.0, 0.0, 1, n

    design = sm.add_constant(X, has_constant="add")
    ok = mask & np.isfinite(design).all(axis=1)
    p = design.shape[1]
    n = int(ok.sum())
    if n <= p + 1:
        return np.nan, np.nan, p, n
    try:
        model = sm.OLS(y[ok], design[ok]).fit()
    except np.linalg.LinAlgError as exc:
        # Session 59 (B-19c / D3, decided 2026-09-24). statsmodels' default
        # `pinv` fit goes through LAPACK gesdd, which fails to converge on some
        # rank-deficient site designs (pan-cancer split_seed 15, MYC targets,
        # rung `Mfull_+image`, HPC4 arrays 130202/130211). Retry the SAME
        # least-squares problem with a pivoted-QR solver that handles rank
        # deficiency. Only the failing calls take this path: every stored run
        # completed without it, so no reported number moves. Loud, not silent.
        import warnings

        import scipy.linalg as sla

        warnings.warn(f"OLS default SVD failed ({exc}); refit with pivoted QR "
                      "(scipy gelsy)", RuntimeWarning, stacklevel=2)
        Xd, yd = design[ok], y[ok]
        beta, _, rank, _ = sla.lstsq(Xd, yd, lapack_driver="gelsy")
        resid = yd - Xd @ beta
        tss = float(((yd - yd.mean()) ** 2).sum())
        r2 = 1.0 - float(resid @ resid) / tss
        df_resid = n - int(rank)
        adj = 1.0 - (n - 1) / df_resid * (1.0 - r2) if df_resid > 0 else np.nan
        return float(r2), float(adj), p, n
    return float(model.rsquared), float(model.rsquared_adj), p, n


def _ladder_mask(
    frame: pd.DataFrame, y: np.ndarray, terms: list[str], categorical: set[str]
) -> np.ndarray:
    """Complete cases over the UNION of every term used anywhere in the ladder."""
    ok = np.isfinite(y)
    for t in terms:
        if t and t in frame.columns:
            col = frame[t]
            ok &= (col.notna().to_numpy() if t in categorical
                   else np.isfinite(pd.to_numeric(col, errors="coerce").to_numpy()))
    return ok


def decompose_signature(
    frame: pd.DataFrame,
    signature_col: str,
    *,
    image_pred_col: str,
    purity_col: str = "purity",
    type_col: str = "cancer_type",
    site_col: str = "tss",
    stage_col: str | None = "stage",
) -> DecompositionResult:
    """Nested variance decomposition for one signature.

    `frame` must be one row per patient (not per slide) with the signature
    score, the cross-validated out-of-fold image prediction, and the covariates.
    """
    categorical = {type_col, site_col, stage_col} - {None}
    y = pd.to_numeric(frame[signature_col], errors="coerce").to_numpy()

    ladder: list[tuple[str, list[str]]] = [
        ("M0_intercept", []),
        ("M1_purity", [purity_col]),
        ("M2_+type", [purity_col, type_col]),
        ("M3_+site", [purity_col, type_col, site_col]),
    ]
    if stage_col:
        ladder.append(("M4_+stage", [purity_col, type_col, site_col, stage_col]))

    covariate_terms = ladder[-1][1]
    ladder.append(("Mfull_+image", covariate_terms + [image_pred_col]))
    ladder.append(("Image_only", [image_pred_col]))

    # ONE mask, over every term the ladder touches, so all rungs share a sample.
    all_terms = sorted({t for _, terms in ladder for t in terms if t})
    mask = _ladder_mask(frame, y, all_terms, categorical)

    rows = []
    r2_by_name: dict[str, float] = {}
    adj_by_name: dict[str, float] = {}
    for name, terms in ladder:
        r2, adj, k, n_used = _r2(y, _design(frame, terms, categorical), mask)
        r2_by_name[name] = r2
        adj_by_name[name] = adj
        rows.append({"model": name, "terms": ", ".join(terms) or "(intercept)",
                     "r2": r2, "adj_r2": adj, "n_params": k, "n_used": n_used})

    table = pd.DataFrame(rows)
    table["incremental_r2"] = table["r2"].diff().fillna(table["r2"])
    table["incremental_adj_r2"] = table["adj_r2"].diff().fillna(table["adj_r2"])

    # All rungs share `mask`, so n_used must be constant. If it is not, the
    # increment is meaningless and we refuse to report it.
    if table["n_used"].nunique() > 1:
        raise RuntimeError(
            f"ladder rungs used different sample sizes {sorted(table['n_used'].unique())}; "
            "the shared complete-case mask failed"
        )

    base = adj_by_name.get("M4_+stage", adj_by_name.get("M3_+site", 0.0))
    incremental = float(adj_by_name["Mfull_+image"] - base)

    return DecompositionResult(
        signature=signature_col,
        table=table,
        incremental_image_r2=incremental,
        n=int(mask.sum()),
    )


def decompose_all(
    frame: pd.DataFrame,
    signature_cols: list[str],
    pred_cols: dict[str, str],
    **kwargs,
) -> pd.DataFrame:
    """Run the decomposition across every signature and tidy the result."""
    rows = []
    for sig in signature_cols:
        pred = pred_cols.get(sig)
        if pred is None or pred not in frame.columns:
            continue
        res = decompose_signature(frame, sig, image_pred_col=pred, **kwargs)
        wide = res.table.set_index("model")["r2"].to_dict()
        wide.update(
            {
                "signature": sig,
                "incremental_image_r2": res.incremental_image_r2,
                "n": res.n,
            }
        )
        rows.append(wide)
    out = pd.DataFrame(rows)
    if "incremental_image_r2" in out.columns:
        out = out.sort_values("incremental_image_r2", ascending=False)
    return out


def label_site_variance(
    frame: pd.DataFrame,
    signature_cols: list[str],
    *,
    site_col: str = "tss",
    type_col: str | None = "cancer_type",
    purity_col: str | None = "purity",
    plate_col: str | None = "plate",
    min_site_n: int = 10,
) -> pd.DataFrame:
    """CONTROL C — how much of the signature LABEL is explained by site, image-free.

    This is the novel leg of the study and it is deliberately cheap: it touches
    no embeddings at all, only the transcriptome-derived signature scores and
    the tissue source site parsed from the barcode.

    The argument it enables
    -----------------------
    The site-confounding literature (Howard 2021, PathoROB, de Jong) is entirely
    about site being learnable from the IMAGE. If site also explains variance in
    the LABEL — because collection, fixation, ischaemia time and RNA handling
    vary by institution — then an image model's site-driven performance is
    partly a *label* artefact, not only an image artefact.

    SCOPE, honestly: site effects on TCGA expression are documented (Howard
    et al., Nat Commun 2021; Dawood et al. 2023), so the existence of label-side
    site variance is NOT a new claim. What is not established is the DIFFERENTIAL
    test — whether a curated immune signature carries MORE site variance than a
    size-matched random gene set does, which is why `random_null_r2` below runs
    the random-set null against the site-R^2 itself. Report this as a secondary
    analysis, never as a headline discovery.

    Interpretation
    --------------
    r2_site_only high  ->  the ground truth itself is institution-dependent, and
                           "predicting the signature" partly means "predicting
                           the institution's RNA handling"
    r2_site_only ~0    ->  labels are clean; site confounding lives purely in the
                           image, and Control C is dropped (week-3 gate)

    Site is reported both raw and after adjusting for cancer type, because
    site correlates with tumour type in TCGA and the raw number would otherwise
    absorb a type effect.
    """
    sites = frame[site_col].astype("string")
    counts = sites.value_counts()
    keep = sites.isin(counts[counts >= min_site_n].index)
    sub = frame.loc[keep].copy()
    sub_sites = sites[keep]

    if sub_sites.nunique() < 2:
        raise ValueError(
            f"only {sub_sites.nunique()} site(s) with >={min_site_n} samples; "
            "lower min_site_n or pool cohorts"
        )

    categorical = {site_col, type_col, plate_col} - {None}
    rows = []
    plate_note = None

    for sig in signature_cols:
        y = pd.to_numeric(sub[sig], errors="coerce").to_numpy()
        if np.isfinite(y).sum() < 30:
            continue

        # One shared mask across every rung of THIS analysis too.
        lmask = _ladder_mask(sub, y, [c for c in (site_col, type_col, purity_col) if c], categorical)
        r2_site, adj_site, _, n_site = _r2(y, _design(sub, [site_col], categorical), lmask)
        row = {"signature": sig, "n": n_site,
               "n_sites": int(sub_sites.nunique()),
               "r2_site_only": r2_site, "adj_r2_site_only": adj_site}

        if type_col and type_col in sub.columns:
            r2_type, _, _, _ = _r2(y, _design(sub, [type_col], categorical), lmask)
            r2_both, _, _, _ = _r2(y, _design(sub, [type_col, site_col], categorical), lmask)
            row["r2_type_only"] = r2_type
            row["r2_site_given_type"] = r2_both - r2_type

        # PLATE within SITE: the technical-only lower bound. Plate is batch
        # with no plausible biological interpretation, so whatever variance it
        # explains is a floor for how much "site" variance is pure processing.
        # Session 59 (B-19a, decided 2026-09-24): when the plate column carries
        # no information -- entirely missing, or a single level -- the rung is
        # NOT ESTIMABLE, and it is reported as NaN with a note rather than as
        # the structural zero `_design` produced (dummy-coding one NaN level
        # and dropping the constant). That zero was once read as "the label
        # artefact is not technical" and had to be withdrawn (E17). Frozen
        # `label_site_variance.csv` files are not rewritten; this governs any
        # future run of `11_close_science_gaps.py`.
        if plate_col and plate_col in sub.columns:
            if sub[plate_col].nunique(dropna=True) < 2:
                row["r2_site_plus_plate"] = np.nan
                row["r2_plate_within_site"] = np.nan
                plate_note = (f"'{plate_col}' has {sub[plate_col].nunique(dropna=True)} "
                              "non-missing level(s); plate-within-site is not "
                              "estimable and is NaN, not 0")
            else:
                r2_plate, _, _, _ = _r2(
                    y, _design(sub, [site_col, plate_col], categorical), lmask
                )
                row["r2_site_plus_plate"] = r2_plate
                row["r2_plate_within_site"] = r2_plate - r2_site

        if purity_col and purity_col in sub.columns:
            r2_pur, _, _, _ = _r2(y, _design(sub, [purity_col], categorical), lmask)
            row["r2_purity_only"] = r2_pur

        rows.append(row)

    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values("r2_site_only", ascending=False)
        out.attrs["median_r2_site"] = float(out["r2_site_only"].median())
    if plate_note is not None:
        # `.attrs` does not survive `.to_csv()`; the NaN in the plate columns
        # does, and is the durable record. The note is for in-process callers.
        out.attrs["plate_note"] = plate_note
    return out


def purity_share(
    frame: pd.DataFrame, signature_col: str, image_pred_col: str,
    purity_col: str = "purity"
) -> dict[str, float]:
    """How much of the image prediction's explanatory power is purity-shaped?

    Returns three numbers:
      r_image      correlation of image prediction with the true signature
      r_purity     correlation of purity with the true signature
      r_partial    partial correlation of image prediction with the signature,
                   controlling for purity

    If r_partial collapses toward zero while r_image was healthy, the image
    model was reading purity.
    """
    sub = frame[[signature_col, image_pred_col, purity_col]].apply(
        pd.to_numeric, errors="coerce"
    ).dropna()
    if len(sub) < 10:
        return {"r_image": np.nan, "r_purity": np.nan, "r_partial": np.nan, "n": len(sub)}

    y = sub[signature_col].to_numpy()
    p = sub[image_pred_col].to_numpy()
    u = sub[purity_col].to_numpy()

    def resid(a: np.ndarray, on: np.ndarray) -> np.ndarray:
        X = sm.add_constant(on.reshape(-1, 1), has_constant="add")
        return a - sm.OLS(a, X).fit().predict(X)

    r_partial = float(np.corrcoef(resid(y, u), resid(p, u))[0, 1])
    return {
        "r_image": float(np.corrcoef(y, p)[0, 1]),
        "r_purity": float(np.corrcoef(y, u)[0, 1]),
        "r_partial": r_partial,
        "n": len(sub),
    }


def permutation_calibration(
    frame: pd.DataFrame,
    signature_col: str,
    *,
    site_col: str = "tss",
    within_col: str | None = "cancer_type",
    n_perm: int = 1000,
    seed: int = 0,
) -> dict[str, float]:
    """Calibrate site-R^2 against permuted site labels. Report observed MINUS median.

    Why this is mandatory rather than optional
    ------------------------------------------
    Dummy-coded site inflates R^2 by roughly p/(n-1) purely by construction. On
    the real parquet that is 641 distinct TSS codes over 8,535 patients — about
    7.5% of R^2 handed over for free. Under a pure-noise label at n=956 with 68
    sites the raw statistic returned 0.0568 against an expected 67/955 = 0.0702.

    Crucially, that inflation lands on the COVARIATE rungs and therefore
    DEFLATES incremental_image_r2 — a bias pointing toward this project's hoped-
    for conclusion. Reviewers look for exactly this.

    Labels are permuted WITHIN `within_col` (cancer type by default), because
    site correlates with tumour type in TCGA and a global permutation would
    destroy that structure and understate the null.
    """
    rng = np.random.default_rng(seed)
    y = pd.to_numeric(frame[signature_col], errors="coerce").to_numpy()
    categorical = {site_col}
    mask = _ladder_mask(frame, y, [site_col], categorical)

    observed, obs_adj, _, n_used = _r2(y, _design(frame, [site_col], categorical), mask)

    sites = frame[site_col].astype("string")
    groups = (frame[within_col].astype("string") if within_col and within_col in frame.columns
              else pd.Series(["_all"] * len(frame), index=frame.index))

    draws = []
    work = frame.copy()
    for _ in range(n_perm):
        permuted = sites.copy()
        for _, idx in groups.groupby(groups).groups.items():
            vals = sites.loc[idx].to_numpy()
            permuted.loc[idx] = rng.permutation(vals)
        work[site_col] = permuted
        r2p, _, _, _ = _r2(y, _design(work, [site_col], categorical), mask)
        if np.isfinite(r2p):
            draws.append(r2p)

    if not draws:
        return {"observed_r2": observed, "n_perm": 0}

    d = np.asarray(draws)
    return {
        "observed_r2": float(observed),
        "observed_adj_r2": float(obs_adj),
        "perm_median_r2": float(np.median(d)),
        "calibrated_r2": float(observed - np.median(d)),
        "perm_p": float((1 + np.sum(d >= observed)) / (len(d) + 1)),
        "n_perm": int(len(d)),
        "n": int(n_used),
    }


def out_of_fold_r2(
    frame: pd.DataFrame,
    signature_col: str,
    terms: list[str],
    split,
    *,
    categorical: set[str] | None = None,
) -> float:
    """R^2 of out-of-fold predictions, using the SAME preserved-site folds.

    In-sample OLS R^2 always rises when you add parameters; out-of-fold R^2 does
    not. Report both — a large gap is itself the overfitting diagnostic.
    """
    from sklearn.linear_model import LinearRegression

    categorical = categorical or set()
    y = pd.to_numeric(frame[signature_col], errors="coerce").to_numpy()
    X = _design(frame, terms, categorical)
    if X is None:
        return float("nan")
    X = np.nan_to_num(X, nan=0.0)

    preds = np.full(len(y), np.nan)
    for k in range(split.n_folds):
        tr, te = split.indices(k)
        tr = tr[np.isfinite(y[tr])]
        if len(tr) < 10 or len(te) == 0:
            continue
        preds[te] = LinearRegression().fit(X[tr], y[tr]).predict(X[te])

    ok = np.isfinite(preds) & np.isfinite(y)
    if ok.sum() < 10:
        return float("nan")
    ss_res = float(np.sum((y[ok] - preds[ok]) ** 2))
    ss_tot = float(np.sum((y[ok] - y[ok].mean()) ** 2))
    return float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
