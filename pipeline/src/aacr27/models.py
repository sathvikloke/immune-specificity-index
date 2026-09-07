"""Models and the cross-validated evaluation harness.

The substrate is SLIDE-LEVEL embeddings (one vector per whole-slide image), so
the model is a regularised linear head, not attention-MIL. That is a feature,
not a compromise: it runs on CPU in minutes, it has no tuning surface to overfit
through, and — critically for an audit paper — it is the same class of head the
literature uses for linear probing, so a weak result cannot be dismissed as
"you used a bad model".

Baselines are first-class here. Audit section 8.5 item 3: a model that beats
nothing has established nothing. Every signature is evaluated against:

    site_only     one-hot tissue source site  (the confounder alone)
    type_only     one-hot cancer type
    purity_only   ABSOLUTE tumour purity
    covariates    site + type + purity + stage
    embedding     the image embedding
    emb+cov       image embedding on top of the covariates

If `site_only` matches `embedding`, the image is a site detector.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.preprocessing import StandardScaler

from .splits import SplitResult


@dataclass
class FoldPredictions:
    """Out-of-fold predictions for one target under one split scheme."""

    target: str
    scheme: str
    model: str
    y_true: np.ndarray
    y_pred: np.ndarray
    patient: np.ndarray
    fold: np.ndarray
    meta: dict = field(default_factory=dict)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "patient": self.patient,
                "fold": self.fold,
                "y_true": self.y_true,
                "y_pred": self.y_pred,
                "target": self.target,
                "scheme": self.scheme,
                "model": self.model,
            }
        )


ALPHAS = np.logspace(-2, 5, 24)


def cross_val_predict_multi(
    X: np.ndarray,
    Y: np.ndarray,
    split: SplitResult,
    patients: np.ndarray,
    *,
    target_names: list[str] | None = None,
    alphas: np.ndarray = ALPHAS,
    standardise: bool = True,
    fixed_alpha: float | None = None,
) -> np.ndarray:
    """Out-of-fold predictions for MANY targets at once. Returns (n_samples, n_targets).

    `RidgeCV(alpha_per_target=True)` selects a separate alpha per column while
    sharing one SVD of the design matrix per fold. Verified numerically identical
    to looping `cross_val_predict` over columns (max abs diff ~1e-15, identical
    selected alphas) and roughly 200x faster on a 29,000-target block.

    That speedup is what makes B=1000 null draws affordable: the loop form costs
    about 12 CPU-hours for the full null, this costs minutes.

    Imputation and standardisation are fitted INSIDE each training fold.

    `fixed_alpha` skips the per-target alpha search, holding regularisation
    constant across the observed-vs-null contrast instead of letting each side
    tune itself. Pass a scalar, or a per-fold sequence of length `n_folds` so the
    alpha for fold k can be selected on fold k's TRAINING partition only — the
    latter is what `experiment` uses, because a single alpha chosen on all
    patients is a (small) selection leak into the null.

    MEASURED DIRECTION, because the original rationale here was wrong. It claimed
    per-target alpha selection "gives the null a best-case advantage". It does
    not: on NSCLC the null's mean r was 0.0446 under a fixed alpha and 0.0398
    under per-target RidgeCV. Per-target selection picks HEAVIER regularisation
    for noise-like targets, shrinking predictions toward the mean and LOWERING
    the null. So a fixed alpha is the conservative choice — it hands the null
    more, not less — and that is the reason to prefer it.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    if Y.ndim == 1:
        Y = Y.reshape(-1, 1)
    preds = np.full(Y.shape, np.nan)

    per_fold_alpha = None
    if fixed_alpha is not None and np.ndim(fixed_alpha) > 0:
        per_fold_alpha = np.asarray(fixed_alpha, dtype=float)
        if len(per_fold_alpha) != split.n_folds:
            raise ValueError(
                f"fixed_alpha has {len(per_fold_alpha)} entries but the split has "
                f"{split.n_folds} folds"
            )

    for k in range(split.n_folds):
        train_idx, test_idx = split.indices(k)
        if len(train_idx) < 10 or len(test_idx) == 0:
            continue
        if per_fold_alpha is not None:
            fold_alpha = per_fold_alpha[k]
            fixed_alpha = None if not np.isfinite(fold_alpha) else float(fold_alpha)

        med = fit_impute_median(X[train_idx])
        Xtr = impute_median(X[train_idx], med)
        Xte = impute_median(X[test_idx], med)
        if standardise:
            scaler = StandardScaler().fit(Xtr)
            Xtr, Xte = scaler.transform(Xtr), scaler.transform(Xte)

        # Targets with any non-finite value in the training rows are fitted on
        # their own complete cases; a single NaN must not drop the whole block.
        Ytr = Y[train_idx]
        usable = np.isfinite(Ytr).all(axis=0)
        if usable.any():
            if fixed_alpha is not None:
                est = Ridge(alpha=fixed_alpha).fit(Xtr, Ytr[:, usable])
            else:
                est = RidgeCV(alphas=alphas, alpha_per_target=True).fit(Xtr, Ytr[:, usable])
            # sklearn RAVELS a single-column target, so with exactly one usable
            # column predict() returns (n,) and the 2-D assignment below raises
            # "shape mismatch". The multi-target paths in this project normally
            # pass hundreds of columns, which is why this never fired — but the
            # legacy Venet null can filter down to one, and a crash there would
            # take out the whole run.
            block = est.predict(Xte)
            preds[np.ix_(test_idx, np.where(usable)[0])] = (
                block.reshape(-1, 1) if block.ndim == 1 else block
            )

        for j in np.where(~usable)[0]:
            ok = np.isfinite(Ytr[:, j])
            if ok.sum() < 10:
                continue
            est = (Ridge(alpha=fixed_alpha) if fixed_alpha is not None
                   else RidgeCV(alphas=alphas)).fit(Xtr[ok], Ytr[ok, j])
            preds[test_idx, j] = est.predict(Xte)

    return preds


def cross_val_predict(
    X: np.ndarray,
    y: np.ndarray,
    split: SplitResult,
    patients: np.ndarray,
    *,
    target: str = "",
    model_name: str = "ridge",
    alphas: np.ndarray = ALPHAS,
    standardise: bool = True,
) -> FoldPredictions:
    """Out-of-fold predictions with the regularisation tuned INSIDE each fold.

    Tuning alpha on the full dataset and then cross-validating is a textbook
    leakage mode (Kapoor & Narayanan, Patterns 2023, taxonomy of eight leakage
    types). RidgeCV is fitted on the training partition only, every fold.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    preds = np.full(len(y), np.nan)

    for k in range(split.n_folds):
        train_idx, test_idx = split.indices(k)
        train_idx = train_idx[np.isfinite(y[train_idx])]
        if len(train_idx) < 10 or len(test_idx) == 0:
            continue

        # Impute and standardise using TRAINING-FOLD statistics only. Fitting
        # either on the full matrix leaks test information into training
        # (Kapoor & Narayanan, Patterns 2023, preprocessing leakage).
        med = fit_impute_median(X[train_idx])
        Xtr = impute_median(X[train_idx], med)
        Xte = impute_median(X[test_idx], med)
        if standardise:
            scaler = StandardScaler().fit(Xtr)
            Xtr, Xte = scaler.transform(Xtr), scaler.transform(Xte)

        estimator = RidgeCV(alphas=alphas).fit(Xtr, y[train_idx])
        preds[test_idx] = estimator.predict(Xte)

    return FoldPredictions(
        target=target,
        scheme=split.scheme,
        model=model_name,
        y_true=y,
        y_pred=preds,
        patient=np.asarray(patients),
        fold=split.fold.to_numpy(),
    )


def build_baseline_matrix(
    frame: pd.DataFrame,
    kind: str,
    *,
    site_col: str = "tss",
    type_col: str = "cancer_type",
    purity_col: str = "purity",
    stage_col: str | None = "stage",
) -> np.ndarray:
    """Design matrix for a named baseline. See module docstring for the menu."""
    def onehot(col: str) -> np.ndarray:
        return pd.get_dummies(
            frame[col].astype("string"), drop_first=True, dtype=float
        ).to_numpy()

    def numeric(col: str) -> np.ndarray:
        return pd.to_numeric(frame[col], errors="coerce").to_numpy().reshape(-1, 1)

    if kind == "cohort_mean":
        # The sanity floor: predict the cohort mean for everyone. A model that
        # does not beat this has established nothing at all.
        return np.zeros((len(frame), 1))
    if kind == "type_mean":
        # Predict the per-cancer-type mean. Pan-cancer, this is the baseline that
        # tissue-of-origin alone achieves, and it is the one most easily mistaken
        # for biology.
        return onehot(type_col)
    if kind == "site_only":
        return onehot(site_col)
    if kind == "type_only":
        return onehot(type_col)
    if kind == "purity_only":
        return numeric(purity_col)
    if kind == "covariates":
        blocks = [onehot(site_col), onehot(type_col), numeric(purity_col)]
        if stage_col and stage_col in frame.columns:
            blocks.append(onehot(stage_col))
        return np.hstack(blocks)
    raise ValueError(
        f"unknown baseline {kind!r}; expected one of cohort_mean, type_mean, "
        "site_only, type_only, purity_only, covariates"
    )


def impute_median(X: np.ndarray, medians: np.ndarray | None = None) -> np.ndarray:
    """Column-median imputation. Returns a COPY; never mutates the input.

    `np.asarray(X, dtype=float)` does not copy an array that is already float64,
    so the previous in-place assignment silently rewrote the caller's embedding
    matrix. Every downstream call then saw imputed values it never asked for.

    Pass `medians` to apply a training-fold statistic to a test fold. Computing
    medians over the full matrix and then cross-validating is preprocessing
    leakage (Kapoor & Narayanan, Patterns 2023) — the exact failure this module's
    docstring cites. Use `fit_impute_median` inside the fold instead.
    """
    X = np.array(X, dtype=float, copy=True)
    if medians is None:
        medians = fit_impute_median(X)
    idx = np.where(~np.isfinite(X))
    X[idx] = np.take(medians, idx[1])
    return X


def fit_impute_median(X: np.ndarray) -> np.ndarray:
    """Column medians from a TRAINING partition only."""
    medians = np.nanmedian(np.asarray(X, dtype=float), axis=0)
    return np.where(np.isfinite(medians), medians, 0.0)


def site_prediction_control(
    X: np.ndarray,
    sites: pd.Series,
    patients: pd.Series,
    *,
    n_folds: int = 5,
    seed: int = 0,
    min_site_n: int = 20,
) -> pd.DataFrame:
    """THE negative control: how well do the embeddings predict tissue source site?

    Report this next to every performance number. Per audit section 8.5 item 1,
    "AUROC 0.89, and site-prediction AUROC after harmonisation = 0.51" is far
    more credible than 0.94 with no site analysis at all.

    One-vs-rest AUROC per site with at least `min_site_n` patients, using
    patient-disjoint folds.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold

    X = np.asarray(X, dtype=float)
    counts = sites.value_counts()
    eligible = counts[counts >= min_site_n].index.tolist()

    rows = []
    for site in eligible:
        y = (sites == site).astype(int).to_numpy()
        if y.sum() < 5 or (1 - y).sum() < 5:
            continue

        oof = np.full(len(y), np.nan)
        splitter = GroupKFold(n_splits=min(n_folds, len(np.unique(patients))))
        for train_idx, test_idx in splitter.split(X, y, groups=patients.to_numpy()):
            if y[train_idx].sum() < 2 or (1 - y[train_idx]).sum() < 2:
                continue
            scaler = StandardScaler().fit(X[train_idx])
            clf = LogisticRegression(max_iter=2000, C=1.0).fit(
                scaler.transform(X[train_idx]), y[train_idx]
            )
            oof[test_idx] = clf.predict_proba(scaler.transform(X[test_idx]))[:, 1]

        ok = np.isfinite(oof)
        if ok.sum() < 20 or len(np.unique(y[ok])) < 2:
            continue
        rows.append(
            {
                "site": site,
                "n_patients": int(counts[site]),
                "auroc": float(roc_auc_score(y[ok], oof[ok])),
            }
        )

    # NOTHING is stashed in `.attrs` here, deliberately. Two summaries used to
    # be -- `median_auroc` and `frac_above_0.9` -- and both were dead on
    # arrival: measured 2026-09-06, nothing in the repository read either, and
    # `.attrs` does not survive `.to_csv()` (verified: an empty dict comes back
    # from a round-trip), so they could not have reached a caller even if one
    # had wanted them. That is exactly how A2 silently lost the rotation-null
    # p-value. The reported median site AUROC is computed where it is used, by
    # `19_check_numbers.py` reading `site_control.csv`, which is a channel that
    # survives being written to disk.
    return pd.DataFrame(rows).sort_values("auroc", ascending=False)


def fit_combat(
    X: np.ndarray, batch: pd.Series, *, covariates: np.ndarray | None = None
) -> dict:
    """Estimate location/scale batch parameters on a TRAINING partition."""
    X = np.asarray(X, dtype=float)
    batch = pd.Series(batch).astype("string").reset_index(drop=True)

    beta = None
    if covariates is not None:
        design = np.hstack([np.ones((len(X), 1)), np.asarray(covariates, dtype=float)])
        beta, *_ = np.linalg.lstsq(design, X, rcond=None)
        X = X - design @ beta + design[:, :1] @ beta[:1]

    grand_mean = X.mean(axis=0)
    grand_sd = X.std(axis=0)
    grand_sd[grand_sd == 0] = 1.0

    per_batch = {}
    for b in batch.unique():
        mask = (batch == b).to_numpy()
        if mask.sum() < 3:
            continue
        sd = X[mask].std(axis=0)
        sd[sd == 0] = 1.0
        per_batch[b] = (X[mask].mean(axis=0), sd)

    return {"grand_mean": grand_mean, "grand_sd": grand_sd,
            "per_batch": per_batch, "beta": beta}


def apply_combat(X: np.ndarray, batch: pd.Series, params: dict) -> tuple[np.ndarray, np.ndarray]:
    """Apply fitted batch parameters. Returns (corrected_X, was_estimable_mask).

    Batches absent from the fitted parameters CANNOT be corrected — there is no
    training data for them. Those rows are returned uncorrected and flagged.

    This is not an edge case, it is the deployment scenario: under preserved-site
    cross-validation a test site has zero training samples BY CONSTRUCTION, so
    ComBat is unestimable for every held-out site. Reporting that honestly is a
    better finding than a corrected number obtained by fitting on the full matrix.
    """
    X = np.array(X, dtype=float, copy=True)
    batch = pd.Series(batch).astype("string").reset_index(drop=True)
    estimable = np.zeros(len(X), dtype=bool)

    if params.get("beta") is not None:
        raise NotImplementedError(
            "covariate-adjusted apply is not implemented; fit without covariates "
            "or extend this function deliberately"
        )

    for b, (mu, sd) in params["per_batch"].items():
        mask = (batch == b).to_numpy()
        if not mask.any():
            continue
        X[mask] = (X[mask] - mu) / sd * params["grand_sd"] + params["grand_mean"]
        estimable |= mask

    return X, estimable


def combat_correct(
    X: np.ndarray, batch: pd.Series, *, covariates: np.ndarray | None = None
) -> np.ndarray:
    """Fit-and-apply on the SAME data. In-sample only — see the warning.

    WARNING: fitting batch parameters on the full matrix and then evaluating by
    cross-validation is leakage. Per-site centring forces each site's deviations
    to sum to zero, so a held-out member lands opposite its training-fold peers
    and site-prediction AUROC collapses far BELOW 0.5 — an artefact, not a
    correction. An honest out-of-fold correction lands near 0.50.

    Use `fit_combat` on the training partition and `apply_combat` on the test
    partition. This function is retained only for in-sample description, and for
    the demo's illustration of the artefact.

    This is a deliberately simple location/scale correction with no empirical
    Bayes shrinkage. Confirm any published conclusion against a reference
    implementation (e.g. `neurocombat`, `pycombat`).
    """
    params = fit_combat(X, batch, covariates=covariates)
    if covariates is not None:
        # Preserve the original single-shot behaviour for the covariate path.
        X = np.array(X, dtype=float, copy=True)
        design = np.hstack([np.ones((len(X), 1)), np.asarray(covariates, dtype=float)])
        beta = params["beta"]
        X = X - design @ beta + design[:, :1] @ beta[:1]
        params = fit_combat(X, batch)
    corrected, _ = apply_combat(X, batch, params)
    return corrected


def site_prediction_control_oof_combat(
    X: np.ndarray,
    sites: pd.Series,
    patients: pd.Series,
    *,
    n_folds: int = 5,
    seed: int = 0,
    min_site_n: int = 20,
) -> tuple[pd.DataFrame, dict]:
    """Site-prediction AUROC after an HONEST out-of-fold ComBat.

    MEASURED BEHAVIOUR — do not report post-ComBat AUROC as a success metric
    ------------------------------------------------------------------------
    On synthetic data with a planted multi-dimensional site fingerprint:

        raw embeddings          AUROC 1.000
        in-sample ComBat        AUROC 0.049
        out-of-fold ComBat      AUROC 0.020

    Out-of-fold fitting does NOT rescue the number, and an earlier expectation
    that an honest correction "lands near 0.50" is WRONG. ComBat is in fact
    working on its own terms — the spread of per-site means falls from 0.968 to
    0.033 — but the AUROC stays confidently backwards.

    The reason: ComBat imposes a within-site zero-sum constraint on the very
    samples the site classifier trains on. The classifier calibrates a boundary
    to that artificial structure and then reverses on samples that do not carry
    it. So post-ComBat site-AUROC is uninterpretable in EITHER direction: below
    0.5 does not mean over-correction and 0.5 would not mean success.

    Report instead: (a) the raw site AUROC, which is the real confound measure,
    and (b) the estimable fraction below, which answers the deployment question.

    The deeper point, which is a finding rather than a nuisance
    ----------------------------------------------------------
    Under PRESERVED-SITE cross-validation a test site has zero training samples
    by construction, so its batch parameters are unestimable — ComBat cannot be
    applied to a site it has never seen. That is exactly the deployment
    scenario for a model meeting a new hospital. Reporting "ComBat is
    unestimable for held-out sites" is more informative than any corrected
    number, and this function returns the estimable fraction so it can be said
    quantitatively.

    Returns (per_site_auroc_table, diagnostics).
    """
    from sklearn.model_selection import GroupKFold

    X = np.asarray(X, dtype=float)
    sites = pd.Series(sites).reset_index(drop=True)
    patients = pd.Series(patients).reset_index(drop=True)

    corrected = np.array(X, dtype=float, copy=True)
    estimable_all = np.zeros(len(X), dtype=bool)

    splitter = GroupKFold(n_splits=min(n_folds, patients.nunique()))
    for train_idx, test_idx in splitter.split(X, groups=patients.to_numpy()):
        params = fit_combat(X[train_idx], sites.iloc[train_idx])
        corrected[test_idx], est = apply_combat(
            X[test_idx], sites.iloc[test_idx], params
        )
        estimable_all[test_idx] = est

    table = site_prediction_control(
        corrected, sites, patients, n_folds=n_folds, seed=seed, min_site_n=min_site_n
    )
    diagnostics = {
        "estimable_fraction": float(estimable_all.mean()),
        "n_unestimable": int((~estimable_all).sum()),
        "auroc_median": float(table["auroc"].median()) if len(table) else float("nan"),
        "auroc_interpretable": False,
        "note": (
            "post-ComBat site AUROC is NOT interpretable (see docstring): the "
            "correction's within-site zero-sum constraint makes the classifier "
            "systematically anti-predictive. Use the raw AUROC as the confound "
            "measure and the estimable fraction as the deployment answer."
        ),
    }
    return table, diagnostics


def combat_estimable_fraction(sites: pd.Series, split) -> dict:
    """What fraction of held-out samples could ComBat even be applied to?

    This is the deployment question, and under preserved-site cross-validation
    the answer is 0% BY CONSTRUCTION: a held-out site contributes no training
    samples, so its batch parameters do not exist. A model meeting a new
    hospital is in exactly that position.

    That is a cleaner and more honest finding than any corrected AUROC.
    """
    sites = pd.Series(sites).reset_index(drop=True)
    estimable = np.zeros(len(sites), dtype=bool)
    for k in range(split.n_folds):
        train_idx, test_idx = split.indices(k)
        seen = set(sites.iloc[train_idx].unique())
        estimable[test_idx] = sites.iloc[test_idx].isin(seen).to_numpy()
    return {
        "estimable_fraction": float(estimable.mean()),
        "n_estimable": int(estimable.sum()),
        "n_total": int(len(estimable)),
        "scheme": split.scheme,
    }
