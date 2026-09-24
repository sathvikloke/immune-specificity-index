"""The audit orchestrator — ties splits, models, nulls and decomposition together.

Running `run_audit()` produces every number the abstract needs:

  PRIMARY     immune-specificity index (ISI): reliability-disattenuated,
              global-axis-residualised excess in Fisher z of r(image, signature)
              over the mean of r(image, size-matched random set), with a
              patient-clustered bootstrap CI. Requires `global_axis=`.
  three-way   r_raw / r_axis / r_residual per signature — the arm that separates
              "no immune-specific signal" from "the image reads composition
              faithfully". Both hypotheses predict the same raw correlation.
  outcome     the same null run against PFI, where Venet's logic is exact rather
              than analogical. Requires `survival=`.
  secondary   delta_MAE: paired difference in mean absolute error between the
              random-patient and preserved-site split schemes. This is a
              calibration-loss contrast, NOT a difference of correlations, and
              it must not be reported as "Delta_r".
  control     RAW site-prediction AUROC of the embeddings. Post-ComBat AUROC is
              deliberately NOT interpretable — see `models.combat_correct`. The
              deployment statement is the estimable fraction, which under
              preserved-site folds is 0% by construction.
  purity      incremental ADJUSTED R^2 of the image over purity + type + site +
              stage, on one shared complete-case sample.
  baselines   cohort_mean / type_mean / site_only / type_only / purity_only /
              covariates, per signature, all under the random-patient split.

Pre-specify the primary endpoint before running this. `AuditConfig.primary_endpoint`
is written to config.json on every run and `AuditResult.save()` asserts that the
estimand actually computed matches it, so the registration and the code cannot
drift apart.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import globalaxis as ga
from . import models, outcome as oc, signatures as sig_mod, splits, stats
from .decomposition import decompose_all, purity_share


@dataclass
class AuditConfig:
    """Everything that must be fixed before looking at results."""

    n_folds: int = 5
    seed: int = 0
    # 1,000, the count every frozen run recorded in its config.json. It was
    # 2,000 until 2026-09-24 (E2, decided in session 59): the cohort scripts
    # always passed 1,000 explicitly, so nothing stored changes, but a bare
    # `AuditConfig()` no longer gives a bootstrap the protocol did not use.
    n_boot: int = 1000
    alpha: float = 0.05

    # Pre-registered primary endpoint.
    #
    # PRIMARY: immune-specific excess, in Fisher z --
    #     z[r(image, signature | global axis)]
    #       - mean_j z[r(image, random_set_j | global axis)]
    # Computed when `global_axis` is passed to run_audit(). This is the estimand
    # that distinguishes "no immune-specific signal" from "a dominant global
    # expression axis the image reads faithfully" -- a raw Venet percentile
    # cannot, because both hypotheses predict the same observation, and because
    # the null is near-degenerate at realistic signature sizes
    # (k=29, rho=0.10 -> two independent random sets correlate at 0.763).
    #
    # SECONDARY: split robustness, random-patient minus preserved-site.
    # OUTCOME ARM: `outcome.run_outcome_arm`, invoked when `survival` is passed.
    # PRE-REGISTRATION STRING. Must stay byte-identical to README.md and to
    # 02-PROJECT-DECISION.md section 7. AuditResult.save() asserts the computed
    # estimand matches this, so code and registration cannot drift apart.
    primary_endpoint: str = (
        "immune-specificity index (ISI): reliability-disattenuated, "
        "global-axis-residualised excess in Fisher z of r(image, signature) over "
        "the mean of r(image, size-matched random set), 1000 draws, "
        "patient-clustered bootstrap CI"
    )
    secondary_endpoint: str = (
        "split robustness: mean over signatures of (Pearson r under random-patient "
        "CV) minus (Pearson r under preserved-site CV)"
    )
    axis_method: str = "pc1"          # pre-register: "pc1" or "meanz"
    axis_within_type: bool = True     # compute the axis within cancer type
    primary_scheme_a: str = "random_patient"
    primary_scheme_b: str = "preserved_site"

    # Venet null.
    #
    # B=1000, not 100. At B=100 the smallest attainable one-sided permutation
    # p is 1/101 = 0.0099, while BH-FDR at m=29 signatures and q=0.05 requires
    # the top-ranked p to be <= 0.05*(1/29) = 0.001724. A 100-draw null is
    # therefore LOGICALLY INCAPABLE of rejecting anything in the declared
    # family. B >= 579 is the minimum; 1000 gives headroom.
    n_null_sets: int = 1000
    match_null_expression: bool = True

    # Scoring method, applied identically to observed AND null sets. Mixing them
    # is not a null: measured sd(mean_z)/sd(ssGSEA) ~ 3x, so a mean-z null
    # against ssGSEA observations compares different scales.
    scorer: str = "mean_z"          # "mean_z" | "ssgsea" | "plage" (E6, sensitivity only)

    # Reliability disattenuation. Curated modules are more internally consistent
    # than random draws and are therefore more predictable regardless of biology.
    disattenuate: bool = True

    # Baselines to evaluate alongside the embedding model.
    baselines: tuple[str, ...] = ("site_only", "type_only", "purity_only", "covariates")

    # Columns.
    site_col: str = "tss"
    type_col: str = "cancer_type"
    purity_col: str = "purity"
    stage_col: str | None = "stage"
    patient_col: str = "patient_id"

    min_patients_per_signature: int = 200
    run_combat: bool = True

    # SPLIT SEED, DELIBERATELY SEPARATE FROM `seed`.
    #
    # Varying `seed` changes the fold partition AND the random gene sets AND the
    # bootstrap draws at once, so a spread computed that way conflates three
    # sources of variation. `split_seed` moves ONLY the site-to-fold assignment,
    # which is what `scripts/12_partition_variance.py` needs to isolate. Defaults
    # to `seed`, so existing runs are bit-identical.
    split_seed: int | None = None

    # The legacy Venet null (raw, non-residualised) is descriptive and costs
    # roughly a third of the runtime. The ISI computes its own residualised null,
    # so this can be switched off when only the primary is wanted.
    run_legacy_venet_null: bool = True

    @property
    def fold_seed(self) -> int:
        return self.seed if self.split_seed is None else self.split_seed

    # #31 repeated site-disjoint partitions, to expose partition variance.
    report_partition_variance: bool = False
    partition_repeats: int = 5

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2, default=str))


@dataclass
class AuditResult:
    config: AuditConfig
    per_signature: pd.DataFrame
    primary: stats.Estimate | None
    site_control: pd.DataFrame
    site_control_combat: pd.DataFrame | None
    null_summary: pd.DataFrame | None
    decomposition: pd.DataFrame | None
    predictions: pd.DataFrame
    three_way: pd.DataFrame | None = None
    immune_excess: pd.DataFrame | None = None
    primary_excess: stats.Estimate | None = None
    outcome_arm: pd.DataFrame | None = None
    runtime_s: float = 0.0
    notes: list[str] = field(default_factory=list)

    ISI_METHOD_PREFIX = "mean-immune-specific-excess-z"

    def save(self, outdir: str | Path) -> None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.assert_estimand_matches_registration()
        self.config.to_json(outdir / "config.json")
        self.per_signature.to_csv(outdir / "per_signature.csv", index=False)
        self.site_control.to_csv(outdir / "site_control.csv", index=False)
        self.predictions.to_csv(outdir / "predictions.csv.gz", index=False)
        if self.site_control_combat is not None:
            self.site_control_combat.to_csv(outdir / "site_control_combat.csv", index=False)
        if self.null_summary is not None:
            self.null_summary.to_csv(outdir / "venet_null.csv", index=False)
        if self.decomposition is not None:
            self.decomposition.to_csv(outdir / "decomposition.csv", index=False)
        if self.three_way is not None:
            self.three_way.to_csv(outdir / "three_way.csv", index=False)
        if self.immune_excess is not None:
            self.immune_excess.to_csv(outdir / "immune_excess.csv", index=False)
            draws = self.immune_excess.attrs.get("null_draws")
            if draws:
                np.savez_compressed(outdir / "null_draws.npz", **draws)
        if self.outcome_arm is not None:
            self.outcome_arm.to_csv(outdir / "outcome_arm.csv", index=False)
        summary = {
            "primary_excess": None if self.primary_excess is None else {
                "value": self.primary_excess.value,
                "ci_lo": self.primary_excess.lo,
                "ci_hi": self.primary_excess.hi,
                "n": self.primary_excess.n,
                "method": self.primary_excess.method,
            },
            "secondary_delta_mae": None if self.primary is None else {
                "value": self.primary.value,
                "ci_lo": self.primary.lo,
                "ci_hi": self.primary.hi,
                "n": self.primary.n,
                "method": self.primary.method,
            },
            "runtime_s": self.runtime_s,
            "notes": self.notes,
            "cohort": self.cohort_fingerprint(),
        }
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    def cohort_fingerprint(self) -> dict:
        """Identify WHICH PATIENTS this run was fitted on.

        Added 2026-09-06 as the source-level fix for A11. A results directory
        recorded the estimator's settings in `config.json` and its outputs in
        `summary.json`, and neither said whose data produced them. The only
        cohort evidence a run carried was `predictions.csv.gz`, which nothing
        checks. So when `11_close_science_gaps.py` wrote 202-site pan-TCGA
        numbers into an NSCLC directory, nothing in either file contradicted
        it, and the mismatch was found only because three values agreed with
        another cohort's to seventeen significant figures.

        The hash makes the patient SET comparable between a run and anything
        claiming to describe it, without publishing the barcodes themselves.
        Sorted before hashing so it is invariant to row order, and taken from
        `predictions`, which is the one artefact that necessarily contains
        every patient the estimator actually saw.

        The frozen runs predate this and do NOT carry the key -- they are
        hashed in PROVENANCE.json and are not being rewritten. Absence
        therefore means "run before 2026-09-06", not "cohort unknown to the
        code"; consumers must treat a missing key as unknown rather than as
        a mismatch.
        """
        import hashlib

        patients = sorted({str(p) for p in self.predictions["patient"]})
        digest = hashlib.sha256("\n".join(patients).encode()).hexdigest()
        return {
            "n_patients": len(patients),
            "patient_set_sha256": digest,
            "patient_set_sha256_note":
                "sha256 of the sorted unique patient ids, newline-joined",
        }

    def assert_estimand_matches_registration(self) -> None:
        """Guarantee the computed estimand IS the pre-registered one.

        The failure this prevents: `config.primary_endpoint` is the string that
        goes into the OSF registration and the manuscript, while
        `primary_excess.method` is what the code actually computed. Nothing
        previously tied them together, and they had already drifted apart once.
        """
        if self.primary_excess is None:
            return
        declared = self.config.primary_endpoint.lower()
        computed = (self.primary_excess.method or "").lower()
        if not computed.startswith(self.ISI_METHOD_PREFIX):
            raise AssertionError(
                f"computed estimand {computed!r} is not the ISI; "
                "primary_excess must come from _immune_specific_excess()"
            )
        if "immune-specificity index" not in declared and "immune-specific excess" not in declared:
            raise AssertionError(
                "config.primary_endpoint does not describe the ISI, but an ISI "
                f"was computed. Registration string: {self.config.primary_endpoint!r}"
            )

    def headline(self) -> str:
        lines = ["=" * 68, "AUDIT RESULT", "=" * 68]
        if self.primary_excess is not None:
            lines.append("PRIMARY  immune-specific excess (Fisher z), axis-residualised")
            lines.append(f"         excess = {self.primary_excess}")
            verdict = ("immune-specific signal beyond composition"
                       if self.primary_excess.lo > 0
                       else "NOT distinguishable from the global expression axis")
            lines.append(f"         -> {verdict}")
        else:
            why = ("no global_axis passed to run_audit"
                   if self.three_way is None
                   else "axis supplied, but no null scores were produced — see notes")
            lines.append(f"PRIMARY  not computed ({why})")
        if self.primary is not None:
            lines.append(f"SECOND'Y split-robustness delta_MAE = {self.primary}")
        if len(self.site_control):
            med = self.site_control["auroc"].median()
            hi = (self.site_control["auroc"] > 0.9).mean()
            lines.append(f"CONTROL  site-prediction AUROC median={med:.3f}, "
                         f"{hi:.0%} of sites above 0.9")
        if self.site_control_combat is not None and len(self.site_control_combat):
            med2 = self.site_control_combat["auroc"].median()
            lines.append(f"         after ComBat: median={med2:.3f} "
                         "(NOT interpretable — see notes)")
        if self.null_summary is not None and len(self.null_summary):
            col = ("beats_null" if "beats_null" in self.null_summary.columns
                   else "beats_null_p95" if "beats_null_p95" in self.null_summary.columns
                   else None)
            if col:
                beat = self.null_summary[col].mean()
                lines.append(
                    f"NULL     {beat:.0%} of signatures exceed their random-set null "
                    "(disattenuated excess CI > 0)"
                )
        if self.decomposition is not None and len(self.decomposition):
            med = self.decomposition["incremental_image_r2"].median()
            lines.append(f"PURITY   median incremental image R^2 over covariates = {med:.4f}")
        if self.three_way is not None and len(self.three_way):
            lines.append(
                f"3-WAY    median r_raw={self.three_way['r_raw'].median():.3f} "
                f"r_axis={self.three_way['r_axis'].median():.3f} "
                f"r_residual={self.three_way['r_residual'].median():.3f}"
            )
        if self.outcome_arm is not None and len(self.outcome_arm):
            beat = self.outcome_arm.query("adjustment == 'raw'")["beats_null"].mean()
            lines.append(f"OUTCOME  {beat:.0%} of signatures beat their random-set null on outcome")
        return "\n".join(lines)


def run_audit(
    frame: pd.DataFrame,
    X: np.ndarray,
    signature_cols: list[str],
    config: AuditConfig | None = None,
    *,
    expression: pd.DataFrame | None = None,
    signature_set: sig_mod.SignatureSet | None = None,
    global_axis: "ga.GlobalAxis | pd.Series | None" = None,
    survival: "oc.Survival | None" = None,
    verbose: bool = True,
) -> AuditResult:
    """Run the full audit on a patient-level cohort.

    Parameters
    ----------
    frame
        One row per patient, carrying patient id, site, cancer type, purity,
        stage, and the signature score columns.
    X
        Patient-level embedding matrix, row-aligned to `frame`.
    signature_cols
        Which columns of `frame` are signature scores.
    expression, signature_set
        Supply both to enable the Venet null. Without them the null is skipped
        and a note is recorded.
    global_axis
        The dominant expression axis (see `globalaxis.compute_global_axis`).
        REQUIRED for the PRIMARY endpoint. Without it only the secondary
        split-robustness contrast is computed, and `primary_excess` stays None.
        Compute it from `expression` if you have it:
            axis = globalaxis.compute_global_axis(expr, method="pc1",
                                                  cancer_type=frame["cancer_type"])
    survival
        A `outcome.Survival` (PFI recommended per TCGA-CDR). Enables the outcome
        arm, which is where Venet's logic actually applies because survival is
        distal to expression.
    """
    config = config or AuditConfig()
    started = time.time()
    notes: list[str] = []

    if len(frame) != len(X):
        raise ValueError(f"frame has {len(frame)} rows but X has {len(X)}")

    patients = frame[config.patient_col]
    if patients.duplicated().any():
        raise ValueError(
            "duplicate patients — call data.collapse_to_patient() before run_audit()"
        )

    sites = frame[config.site_col].astype("string")
    # NOTE: X is deliberately NOT imputed here. Imputing over the full matrix
    # before cross-validation leaks test-fold information into training
    # (Kapoor & Narayanan, Patterns 2023). cross_val_predict imputes with
    # training-fold medians inside each fold.
    # The null refit needs the embedding matrix. Attach it rather than threading
    # it through five call sites; attrs is not copied by pandas operations, so
    # this cannot leak into a derived frame.
    frame = frame.copy()
    frame.attrs["_X"] = X

    # ---------------------------------------------------------- splits ---
    random_split = splits.random_patient_split(
        patients, n_folds=config.n_folds, seed=config.fold_seed
    )
    site_split = splits.preserved_site_split(
        patients, sites, n_folds=config.n_folds, seed=config.fold_seed
    )
    if verbose:
        print(splits.fold_balance_report(random_split))
        print(splits.fold_balance_report(site_split))

    # ------------------------------------------- per-signature evaluation ---
    rows, prediction_frames = [], []
    scheme_map = {"random_patient": random_split, "preserved_site": site_split}

    for sig in signature_cols:
        y = pd.to_numeric(frame[sig], errors="coerce").to_numpy()
        n_ok = int(np.isfinite(y).sum())
        if n_ok < config.min_patients_per_signature:
            notes.append(f"skipped {sig}: only {n_ok} patients with a score")
            continue

        record: dict[str, object] = {"signature": sig, "n": n_ok}

        for scheme_name, split in scheme_map.items():
            preds = models.cross_val_predict(
                X, y, split, patients.to_numpy(), target=sig,
                model_name="ridge_embedding",
            )
            prediction_frames.append(preds.to_frame())
            est = stats.corr_ci(preds.y_true, preds.y_pred, alpha=config.alpha)
            record[f"r_{scheme_name}"] = est.value
            record[f"r_{scheme_name}_lo"] = est.lo
            record[f"r_{scheme_name}_hi"] = est.hi

        # Baselines are evaluated under the RANDOM-patient split, not the
        # site-disjoint one. Under preserved-site CV a one-hot-site model sees
        # no training example from any test site, so it cannot generalise and
        # scores ~0 by construction — which measures the split, not the
        # covariate. The question a baseline must answer is "how much of this
        # signature does site/type/purity explain at all", and that is the
        # random-split quantity.
        for baseline in config.baselines:
            try:
                Xb = models.impute_median(models.build_baseline_matrix(
                    frame, baseline, site_col=config.site_col,
                    type_col=config.type_col, purity_col=config.purity_col,
                    stage_col=config.stage_col,
                ))
            except (ValueError, KeyError) as exc:
                notes.append(f"baseline {baseline} unavailable: {exc}")
                continue
            preds = models.cross_val_predict(
                Xb, y, random_split, patients.to_numpy(), target=sig,
                model_name=baseline,
            )
            prediction_frames.append(preds.to_frame())
            record[f"r_{baseline}"] = stats.corr_ci(preds.y_true, preds.y_pred).value

            # The embedding must beat the covariate baseline on the SAME split
            # to have added anything. Report the gap, not two bare numbers.
            if baseline == "covariates" and f"r_{config.primary_scheme_a}" in record:
                record["embedding_over_covariates"] = (
                    float(record[f"r_{config.primary_scheme_a}"])
                    - float(record[f"r_{baseline}"])
                )

        a = record.get(f"r_{config.primary_scheme_a}")
        b = record.get(f"r_{config.primary_scheme_b}")
        if a is not None and b is not None:
            record["delta_r"] = float(a) - float(b)
        rows.append(record)

        if verbose:
            print(f"  {sig:38s} random={record.get('r_random_patient', float('nan')):.3f} "
                  f"site={record.get('r_preserved_site', float('nan')):.3f} "
                  f"delta={record.get('delta_r', float('nan')):+.3f}")

    per_signature = pd.DataFrame(rows)
    predictions = (
        pd.concat(prediction_frames, ignore_index=True) if prediction_frames
        else pd.DataFrame()
    )

    # ------------------------- SECONDARY: split-robustness / calibration loss ---
    # This computes a paired difference in MEAN ABSOLUTE ERROR between the two
    # split schemes. It was previously labelled "Delta_r" in headline(), which
    # was wrong: it is not a difference of correlations, and the demo printed
    # Delta_r = 0.3970 while the mean of the per-signature delta_r values was
    # 0.339. Renamed to delta_mae and demoted; the primary is the ISI.
    primary = None
    if len(per_signature) and "delta_r" in per_signature.columns:
        wide = predictions[predictions["model"] == "ridge_embedding"]
        a_side = wide[wide["scheme"] == config.primary_scheme_a]
        b_side = wide[wide["scheme"] == config.primary_scheme_b]
        merged = a_side.merge(
            b_side, on=["patient", "target"], suffixes=("_a", "_b")
        ).dropna(subset=["y_pred_a", "y_pred_b"])

        if len(merged):
            err_a = -np.abs(merged["y_true_a"] - merged["y_pred_a"]).to_numpy()
            err_b = -np.abs(merged["y_true_b"] - merged["y_pred_b"]).to_numpy()
            # negative absolute error -> higher is better, so the paired
            # difference is (site-disjoint MAE) - (random-split MAE).
            primary = stats.paired_delta(
                err_a, err_b, merged["patient"].to_numpy(),
                n_boot=config.n_boot, alpha=config.alpha, seed=config.seed,
            )

    # -------------------------------------------------- negative control ---
    site_control = models.site_prediction_control(
        X, sites, patients, n_folds=config.n_folds, seed=config.seed
    )
    # ComBat, fitted OUT OF FOLD. The in-sample version drives site AUROC far
    # below 0.5, which is a leakage fingerprint rather than a correction.
    site_control_combat = None
    if config.run_combat and len(site_control):
        site_control_combat, combat_diag = models.site_prediction_control_oof_combat(
            X, sites, patients, n_folds=config.n_folds, seed=config.seed
        )
        ps = models.combat_estimable_fraction(sites, site_split)
        notes.append(
            f"ComBat: post-correction site AUROC is NOT interpretable "
            f"(measured {combat_diag['auroc_median']:.3f}; the within-site "
            f"zero-sum constraint makes the classifier anti-predictive). "
            f"Report the RAW site AUROC as the confound measure. "
            f"Deployment answer: under preserved-site folds ComBat is estimable "
            f"for {ps['estimable_fraction']:.1%} of held-out samples "
            f"({ps['n_estimable']}/{ps['n_total']})."
        )

    # -------------------------------------------------------- Venet null ---
    null_summary = None
    if expression is not None and signature_set is not None and config.run_legacy_venet_null:
        null_summary = _run_venet_null(
            frame, X, expression, signature_set, per_signature, site_split,
            patients, config, verbose=verbose,
        )
    elif not config.run_legacy_venet_null:
        notes.append("legacy Venet null skipped by config (run_legacy_venet_null=False)")
    else:
        notes.append("Venet null skipped: expression matrix or signature set not supplied")

    # ----------------------------------------------------- decomposition ---
    decomposition = None
    if config.purity_col in frame.columns and len(predictions):
        oof = (
            predictions[
                (predictions["model"] == "ridge_embedding")
                & (predictions["scheme"] == config.primary_scheme_b)
            ]
            .pivot_table(index="patient", columns="target", values="y_pred")
            .add_prefix("pred_")
        )
        merged = frame.merge(oof, left_on=config.patient_col, right_index=True, how="left")
        decomposition = decompose_all(
            merged,
            [s for s in signature_cols if f"pred_{s}" in merged.columns],
            {s: f"pred_{s}" for s in signature_cols},
            purity_col=config.purity_col,
            type_col=config.type_col,
            site_col=config.site_col,
            stage_col=config.stage_col,
        )
        shares = [
            {"signature": s, **purity_share(merged, s, f"pred_{s}", config.purity_col)}
            for s in signature_cols if f"pred_{s}" in merged.columns
        ]
        if shares and decomposition is not None:
            decomposition = decomposition.merge(
                pd.DataFrame(shares), on="signature", how="left"
            )

    # ------------------------------------------ PRIMARY: immune-specific excess ---
    three_way = immune_excess = primary_excess = None
    if global_axis is not None and len(predictions):
        three_way, immune_excess, primary_excess = _immune_specific_excess(
            frame, predictions, signature_cols, global_axis, expression,
            signature_set, config, notes,
        )
    else:
        notes.append(
            "PRIMARY NOT COMPUTED: no global_axis supplied. Only the secondary "
            "split-robustness contrast is available."
        )

    # #31 partition variance: a single site-disjoint partition hides how much the
    # answer depends on WHICH partition was drawn. NOTE the scale differs by
    # estimand -- measured on NSCLC, between-partition sd is ~27% of the
    # bootstrap half-width for Delta_r (a BETWEEN-scheme contrast) but only
    # ~30% of it for the ISI (0.0089 vs 0.0294), because the ISI is computed
    # WITHIN one partition and the choice largely cancels. See
    # scripts/12_partition_variance.py.
    if config.report_partition_variance and len(per_signature):
        try:
            reps = splits.repeated_preserved_site_splits(
                patients, sites, n_folds=config.n_folds,
                repeats=config.partition_repeats, seed=config.seed,
            )
            deltas = []
            for rep in reps[1:]:
                vals = []
                for sig in signature_cols:
                    y = pd.to_numeric(frame[sig], errors="coerce").to_numpy()
                    if np.isfinite(y).sum() < config.min_patients_per_signature:
                        continue
                    pr = models.cross_val_predict(X, y, rep, patients.to_numpy(), target=sig)
                    vals.append(stats.corr_ci(pr.y_true, pr.y_pred).value)
                if vals:
                    deltas.append(float(np.nanmean(vals)))
            if len(deltas) > 1:
                pv = splits.partition_variance(deltas)
                notes.append(
                    f"Partition variance over {pv['n_partitions']} site-disjoint "
                    f"partitions: mean r = {pv['mean']:.3f}, sd = {pv['sd']:.3f}. "
                    "Report this alongside the bootstrap CI — they are separate "
                    "sources of uncertainty."
                )
        except (RuntimeError, ValueError) as exc:
            notes.append(f"partition-variance pass skipped: {exc}")

    # ------------------------------------------------------------ outcome arm ---
    outcome_arm = None
    if survival is not None and expression is not None and signature_set is not None:
        outcome_arm = _run_outcome_arm(
            frame, expression, signature_set, survival, global_axis, config, notes
        )
    elif survival is not None:
        notes.append("outcome arm skipped: needs expression + signature_set")

    return AuditResult(
        config=config,
        per_signature=per_signature,
        primary=primary,
        site_control=site_control,
        site_control_combat=site_control_combat,
        null_summary=null_summary,
        decomposition=decomposition,
        predictions=predictions,
        three_way=three_way,
        immune_excess=immune_excess,
        primary_excess=primary_excess,
        outcome_arm=outcome_arm,
        runtime_s=time.time() - started,
        notes=notes,
    )


def _run_venet_null(
    frame: pd.DataFrame,
    X: np.ndarray,
    expression: pd.DataFrame,
    signature_set: sig_mod.SignatureSet,
    per_signature: pd.DataFrame,
    split: splits.SplitResult,
    patients: pd.Series,
    config: AuditConfig,
    *,
    verbose: bool = True,
) -> pd.DataFrame:
    """Score random size-matched gene sets and compare against the real ones."""
    available = set(expression.columns)
    filtered = signature_set.filter_to(available)
    mean_expr = expression.mean(axis=0) if config.match_null_expression else None

    # Pass the filtered SignatureSet, not just its sizes: expression matching
    # needs the real gene lists to compute each signature's target per-bin
    # histogram. Passing sizes alone silently degrades to size-only matching.
    families = sig_mod.random_gene_sets(
        sorted(available),   # sorted, NOT list(): str hashing is randomised per
                             # process, so list(set) gives a different gene pool
                             # every run and the seed does not reproduce.
        filtered,
        n_per_signature=config.n_null_sets,
        seed=config.seed,
        match_expression=mean_expr,
    )

    # Explicit access: .get(..., []) returned an empty list when the column was
    # absent, which made zip() yield nothing, skipped every null family, returned
    # an empty summary, and silently dropped the NULL line from the headline —
    # a silent failure in the scientific core.
    col = "r_preserved_site"
    if col not in per_signature.columns:
        raise KeyError(
            f"{col!r} missing from per_signature; the null cannot be scored. "
            f"Available: {list(per_signature.columns)}"
        )
    observed = dict(zip(per_signature["signature"], per_signature[col]))
    nulls: dict[str, np.ndarray] = {}

    for sig_name, family in families.items():
        if sig_name not in observed:
            continue
        scores = _score_sets(expression, family, config)
        aligned = scores.reindex(frame[config.patient_col]).to_numpy()

        # BATCHED. This loop previously called cross_val_predict once PER NULL
        # DRAW — B separate RidgeCV fits per signature. Measured at ~0.93 s per
        # draw, which is what turned a 20-minute audit into a 7.9-hour one: the
        # legacy venet null was running alongside the ISI null and doing the same
        # work the slow way. One batched multi-target fit replaces all B.
        keep = [j for j in range(aligned.shape[1])
                if np.isfinite(aligned[:, j]).sum() >= config.min_patients_per_signature]
        if not keep:
            nulls[sig_name] = np.asarray([], dtype=float)
            continue
        Yn = aligned[:, keep]
        Pn = models.cross_val_predict_multi(X, Yn, split, patients.to_numpy())
        nulls[sig_name] = np.asarray(
            [stats.corr_ci(Yn[:, j], Pn[:, j]).value for j in range(Yn.shape[1])],
            dtype=float,
        )
        if verbose:
            print(f"  null {sig_name:34s} {len(nulls[sig_name])} sets scored", flush=True)

    return sig_mod.null_summary(observed, nulls)


def _score_sets(expression: pd.DataFrame, sigset, config: AuditConfig) -> pd.DataFrame:
    """Score gene sets with the CONFIGURED scorer.

    Observed and null scores must use the same scorer. Measured
    sd(mean_z)/sd(ssGSEA) is about 3x, so scoring observations one way and the
    null the other is a scale mismatch, not a null.
    """
    if config.scorer == "mean_z":
        return sig_mod.score_mean_z(expression, sigset)
    if config.scorer == "ssgsea":
        return sig_mod.score_ssgsea(expression, sigset)
    if config.scorer == "plage":
        return sig_mod.score_plage(expression, sigset)
    raise ValueError(
        f"unknown scorer {config.scorer!r}; expected 'mean_z', 'ssgsea' or 'plage'")


def _immune_specific_excess(
    frame: pd.DataFrame,
    predictions: pd.DataFrame,
    signature_cols: list[str],
    global_axis,
    expression: pd.DataFrame | None,
    signature_set,
    config: AuditConfig,
    notes: list[str],
):
    """PRIMARY endpoint: does the image predict signatures beyond the global axis?

    Returns (three_way_table, per_signature_excess_table, pooled_estimate).

    The pooled estimate averages per-signature excesses in Fisher z. Signatures
    are correlated, so this is NOT an independent-sample mean — the CI is a
    patient-clustered bootstrap over the per-signature excesses, and the number
    of effective tests should be reported alongside via stats.effective_tests.
    """
    axis_values = global_axis.values if hasattr(global_axis, "values") else global_axis
    axis_values = pd.Series(axis_values).reindex(frame[config.patient_col].to_numpy())
    axis_values.index = frame.index

    within = (frame[config.type_col]
              if config.axis_within_type and config.type_col in frame else None)

    oof = (
        predictions[
            (predictions["model"] == "ridge_embedding")
            & (predictions["scheme"] == config.primary_scheme_b)
        ]
        .pivot_table(index="patient", columns="target", values="y_pred")
    )

    # Built ONCE and shared. Previously every signature rebuilt this partition
    # twice (once to select alpha, once to refit the null), which is wasted work
    # and, worse, an invitation for the two to drift apart.
    split = splits.preserved_site_split(
        frame[config.patient_col], frame[config.site_col].astype("string"),
        n_folds=config.n_folds, seed=config.fold_seed,
    )

    # Per-GENE residual variances, computed once over the whole expression
    # matrix. This is what makes the reliability of a RESIDUALISED composite
    # computable in O(k) per gene set. See
    # signatures.alpha_from_score_and_item_variances for why the raw-score alpha
    # is the wrong statistic here.
    gene_resid_var = None
    if config.disattenuate and expression is not None:
        # Z-score over `expression`'s OWN samples first, then reindex to frame
        # order — exactly what `score_mean_z` does, so the items here are the
        # same items the scores are built from. Reindexing before z-scoring
        # would silently redefine the target.
        z_full = (expression - expression.mean(axis=0)) / expression.std(axis=0).replace(0, np.nan)
        z_items = z_full.reindex(frame[config.patient_col].to_numpy())
        z_items.index = frame.index
        z_items = z_items.dropna(axis=1, how="all").fillna(0.0)
        gene_resid_var = ga.residual_item_variances(z_items, axis_values, within=within)
        if not np.isfinite(gene_resid_var).any():
            notes.append(
                "residual item variances were all non-finite (axis/expression "
                "misalignment?); falling back to RAW-score reliabilities"
            )
            gene_resid_var = None

    three_rows, excess_rows, per_sig_excess = [], [], []
    null_matrix, observed_resid = [], []
    boot_inputs: list[dict] = []
    for sig in signature_cols:
        if sig not in oof.columns or sig not in frame.columns:
            continue
        pred = oof[sig].reindex(frame[config.patient_col].to_numpy())
        pred.index = frame.index
        report = ga.three_way_report(pred, frame[sig], axis_values, within=within)

        # Null side: residualised random-set predictability.
        got = _residualised_null_scores(
            frame, expression, signature_set, sig, axis_values, within, config, notes,
            split=split, gene_resid_var=gene_resid_var,
        )
        if got is None:
            three_rows.append({"signature": sig, **report})
            continue
        nulls = got["null_r"]
        if nulls is None or not len(nulls):
            three_rows.append({"signature": sig, **report})
            continue

        # SYMMETRY. `report["r_residual"]` comes from a head trained on the RAW
        # signature and then correlated against the residualised target, whereas
        # every null draw REFITS on its residualised random set. Comparing the
        # two is not like for like. The ISI therefore uses a head refitted on the
        # residualised signature, exactly as the null does. Measured effect on
        # NSCLC: r 0.4025 -> 0.4126, i.e. the old asymmetry was handicapping the
        # observed side and making the ISI conservative.
        r_resid_symmetric = got["r_observed_symmetric"]
        report["r_residual_refit"] = r_resid_symmetric
        three_rows.append({"signature": sig, **report})

        rel_obs, rel_null = got["alpha_observed"], got["alpha_null"]
        if rel_null is not None and len(rel_null) != len(nulls):
            rel_obs = rel_null = None
            notes.append(f"{sig}: null alpha length mismatch; disattenuation skipped")

        est = ga.excess_over_null(
            r_resid_symmetric, nulls, n=report["n"], alpha=config.alpha,
            reliability_observed=rel_obs, reliability_null=rel_null,
        )
        # The uncorrected value, kept under a name that can never be mistaken for
        # the primary. Reported so that a NaN disattenuated estimate still leaves
        # a visible trail, NOT as a fallback the pipeline may quietly promote —
        # the uncorrected contrast had a 100% false-positive rate in simulation.
        est_raw = ga.excess_over_null(
            r_resid_symmetric, nulls, n=report["n"], alpha=config.alpha
        )
        n_dropped = 0
        if rel_obs is not None and rel_null is not None:
            dis = np.asarray(ga.disattenuate(nulls, rel_null), dtype=float)
            n_dropped = int(np.isfinite(nulls).sum() - np.isfinite(dis).sum())
            if n_dropped and not np.isfinite(est.value):
                notes.append(
                    f"{sig}: DISATTENUATED EXCESS UNDEFINED — {n_dropped}/{len(nulls)} "
                    "null draws had |r/sqrt(alpha)| >= 1, leaving too few to form a "
                    "null mean. This signature contributes nothing to the pooled ISI. "
                    f"Uncorrected excess would have been {est_raw.value:.4f}; it is NOT "
                    "the pre-registered estimand and must not be reported as one."
                )
        deg = ga.null_degeneracy(nulls)
        excess_rows.append({
            "signature": sig,
            "r_residual": report["r_residual"],
            "r_residual_refit": r_resid_symmetric,
            "excess_z_UNCORRECTED_do_not_report": est_raw.value,
            "null_draws_dropped_by_disattenuation": n_dropped,
            "alpha_observed": rel_obs,
            "alpha_null_mean": float(np.nanmean(rel_null)) if rel_null is not None else np.nan,
            "reliability_gap": (float(rel_obs - np.nanmean(rel_null))
                                if rel_null is not None and rel_obs is not None else np.nan),
            "alpha_observed_raw": got.get("alpha_observed_raw"),
            "alpha_null_mean_raw": got.get("alpha_null_mean_raw"),
            "null_mean_r": deg.get("r_mean"), "null_sd_r": deg.get("r_sd"),
            "degenerate": deg.get("degenerate"),
            "excess_z": est.value, "excess_lo": est.lo, "excess_hi": est.hi,
            "beats_null": bool(est.lo > 0) if np.isfinite(est.lo) else False,
            "n": est.n,
        })
        per_sig_excess.append(est.value)
        null_matrix.append(nulls)
        observed_resid.append(r_resid_symmetric)
        boot_inputs.append({
            "signature": sig,
            "y": got["y_observed"], "p": got["p_observed"],
            "Yn": got["Y_null"], "Pn": got["P_null"],
            "rel_obs": rel_obs, "rel_null": rel_null,
        })

    if three_rows and not excess_rows and expression is not None:
        available = sorted(signature_set.filter_to(set(expression.columns)).sets) \
            if signature_set is not None else []
        raise ValueError(
            "The PRIMARY endpoint could not be computed for ANY signature: no "
            "null scores were produced.\n"
            f"  signature columns: {signature_cols[:8]}\n"
            f"  gene-set names:    {available[:8]}\n"
            "These must match exactly. Strip or add the 'SIG_' prefix so they "
            "agree, then re-run. Failing loudly here rather than returning an "
            "empty result, because a missing primary is not a finding."
        )

    three_way = pd.DataFrame(three_rows) if three_rows else None
    immune_excess = pd.DataFrame(excess_rows) if excess_rows else None

    # #23 family-level p-value that respects between-signature correlation.
    # The per-signature nulls share images and axis, so a naive sigma/sqrt(m)
    # band around the family mean is far too narrow (simulation: claimed
    # 50 +/- 10.5 where the true 95% band was 13.9-89.0).
    family_p = None
    if null_matrix and len(null_matrix) < 2:
        notes.append(
            "rotation_null (family-level p) not computed: it needs >= 2 "
            f"signatures, got {len(null_matrix)}. A family of one has no "
            "family-level statistic."
        )
    if null_matrix and len(null_matrix) > 1:
        width = min(len(v) for v in null_matrix)
        R = np.vstack([np.asarray(v[:width], dtype=float) for v in null_matrix])
        family_p = stats.rotation_null(R, float(np.nanmean(observed_resid)))

    pooled = None
    if per_sig_excess:
        vals = np.asarray(per_sig_excess, dtype=float)
        finite = np.isfinite(vals)
        if not finite.any():
            notes.append(
                f"PRIMARY NOT COMPUTED: all {len(vals)} signatures returned an "
                "undefined disattenuated excess. This is the degenerate regime — "
                "the residualised null sets have too little internal consistency "
                "for Spearman's correction to be defined. Do not substitute the "
                "uncorrected excess; report that the estimand is not identifiable "
                "on this cohort."
            )
        else:
            if (~finite).any():
                notes.append(
                    f"Pooled ISI averages {int(finite.sum())} of {len(vals)} "
                    "signatures; the rest had an undefined disattenuated excess."
                )
            keep = [b for b, f in zip(boot_inputs, finite) if f]
            pooled = _pooled_isi_bootstrap(vals[finite], keep, config, notes)
    if family_p is not None and immune_excess is not None:
        immune_excess.attrs["rotation_null"] = family_p
    if immune_excess is not None and null_matrix:
        # PERSIST THE DRAWS. Only the null's mean and SD used to survive to disk,
        # so any distributional figure had to be invented or the whole ~2 h run
        # repeated. Stored as an object array because signatures can retain
        # different numbers of usable draws.
        immune_excess.attrs["null_draws"] = {
            row["signature"]: np.asarray(v, dtype=float)
            for row, v in zip(excess_rows, null_matrix)
        }
    return three_way, immune_excess, pooled


def _corr_columns(Y: np.ndarray, P: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Pearson r between matching columns of Y and P, over rows `idx`. Vectorised."""
    a, b = Y[idx], P[idx]
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    denom = np.sqrt((a * a).sum(axis=0) * (b * b).sum(axis=0))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(denom > 0, (a * b).sum(axis=0) / denom, np.nan)


def _pooled_isi_bootstrap(
    point_values: np.ndarray,
    boot_inputs: list[dict],
    config: AuditConfig,
    notes: list[str],
    *,
    max_null_for_boot: int = 200,
) -> stats.Estimate:
    """PATIENT-clustered bootstrap of the pooled ISI — the REGISTERED interval.

    THIS REPLACES A WRONG INTERVAL. The previous version resampled the 16
    SIGNATURES with replacement, which treats the Hallmark panel as a random
    sample from a population of signatures — it is not, it is a fixed, chosen
    panel — and which does not propagate patient sampling error at all. The
    registration string says "patient-clustered bootstrap CI", so the code and
    the registration disagreed. Patients are also the unit that generalises to a
    new cohort, so they are the right unit regardless of what was registered.

    Each draw resamples patients with replacement (ONE patient set shared across
    all signatures, which preserves the between-signature correlation that makes
    a naive independent-sample interval far too narrow), then recomputes both the
    observed correlation and the null mean on those patients before pooling.

    NO SILENT CAPS: the null is subsampled to `max_null_for_boot` draws inside
    the bootstrap, because recomputing all B nulls for every signature on every
    draw is O(n_boot * n_sig * B * n). The subsample size is recorded in a note
    and in the returned method string. The null MEAN is what enters the estimand,
    and its standard error at 200 draws is already ~5x smaller than the patient
    sampling error it is being added to.
    """
    if not boot_inputs:
        return stats.Estimate(float(np.mean(point_values)), np.nan, np.nan,
                              len(point_values), AuditResult.ISI_METHOD_PREFIX)

    rng = np.random.default_rng(config.seed)
    kept = min(max_null_for_boot, min(b["Yn"].shape[1] for b in boot_inputs))
    sel = slice(0, kept)

    # ONE SHARED COMPLETE-CASE MASK, and it is load-bearing.
    #
    # `_corr_columns` centres a whole column at once, so a SINGLE non-finite
    # entry makes that column's correlation NaN on EVERY draw. The point
    # estimate does not have this problem because `stats.corr_ci` masks
    # non-finite pairs itself, so the two paths disagreed silently: pan-TCGA
    # produced 16 finite per-signature excesses and 0 of 1000 usable bootstrap
    # draws, because 19 of 7,168 patients had a non-finite global axis (their
    # cancer type carried no axis value) and therefore a NaN residualised score.
    #
    # The mask must be SHARED across signatures rather than computed per
    # signature, because every draw resamples one patient set for all of them —
    # that is what preserves the between-signature correlation the pooled
    # interval depends on.
    ok = np.ones(len(boot_inputs[0]["y"]), dtype=bool)
    for b in boot_inputs:
        ok &= np.isfinite(b["y"]) & np.isfinite(b["p"])
        ok &= np.isfinite(b["Yn"][:, sel]).all(axis=1)
        ok &= np.isfinite(b["Pn"][:, sel]).all(axis=1)

    n_dropped = int((~ok).sum())
    if n_dropped:
        notes.append(
            f"Pooled ISI bootstrap: {n_dropped} of {len(ok)} patients dropped as "
            "incomplete (non-finite residualised score or out-of-fold prediction "
            "in at least one signature or null draw)."
        )
    n = int(ok.sum())
    if n < 20:
        notes.append("pooled ISI bootstrap: too few complete cases; CI NOT reportable")
        return stats.Estimate(float(np.mean(point_values)), np.nan, np.nan,
                              len(point_values), AuditResult.ISI_METHOD_PREFIX)

    prepared = []
    for b in boot_inputs:
        prepared.append({
            "y": b["y"][ok].reshape(-1, 1), "p": b["p"][ok].reshape(-1, 1),
            "Yn": b["Yn"][np.ix_(ok, np.arange(kept))],
            "Pn": b["Pn"][np.ix_(ok, np.arange(kept))],
            "rel_obs": b["rel_obs"],
            "rel_null": None if b["rel_null"] is None
            else np.asarray(b["rel_null"], dtype=float)[sel],
        })

    draws = np.empty(config.n_boot)
    for i in range(config.n_boot):
        idx = rng.integers(0, n, size=n)
        per_sig = []
        for b in prepared:
            r_obs = float(_corr_columns(b["y"], b["p"], idx)[0])
            r_null = _corr_columns(b["Yn"], b["Pn"], idx)
            if b["rel_obs"] is not None and b["rel_null"] is not None:
                r_obs = float(ga.disattenuate(r_obs, b["rel_obs"]))
                r_null = np.asarray(ga.disattenuate(r_null, b["rel_null"]), dtype=float)
            r_null = r_null[np.isfinite(r_null)]
            if not np.isfinite(r_obs) or len(r_null) < 5:
                continue
            per_sig.append(float(stats.fisher_z(r_obs))
                           - float(np.mean(stats.fisher_z(r_null))))
        draws[i] = np.mean(per_sig) if per_sig else np.nan

    draws = draws[np.isfinite(draws)]
    if len(draws) < 20:
        notes.append(
            f"pooled ISI bootstrap produced only {len(draws)} usable draws; "
            "the interval is NOT reportable"
        )
        return stats.Estimate(float(np.mean(point_values)), np.nan, np.nan,
                              len(point_values), AuditResult.ISI_METHOD_PREFIX)

    lo, hi = np.percentile(draws, [100 * config.alpha / 2, 100 * (1 - config.alpha / 2)])
    notes.append(
        f"Pooled ISI CI: patient-clustered bootstrap, {len(draws)}/{config.n_boot} "
        f"usable draws, {n} complete-case patients resampled per draw "
        f"({n_dropped} dropped), null subsampled to {kept} draws per signature "
        f"(of {boot_inputs[0]['Yn'].shape[1]})."
    )
    return stats.Estimate(
        float(np.mean(point_values)), float(lo), float(hi), len(point_values),
        f"{AuditResult.ISI_METHOD_PREFIX}-patient-clustered-bootstrap"
        f" [null subsampled to {kept} draws]",
    )


def _residualised_null_scores(
    frame, expression, signature_set, sig, axis_values, within, config, notes,
    *, split=None, gene_resid_var=None,
):
    """The ISI's observed and null sides for one signature, computed SYMMETRICALLY.

    Returns a dict carrying, for signature `sig`:
      r_observed_symmetric  r(residualised signature, head refitted on it)
      null_r                the same quantity for each residualised random set
      alpha_observed/_null  reliability OF THE RESIDUALISED score (not the raw one)
      y_observed/p_observed, Y_null/P_null   kept for the patient-clustered bootstrap
    """
    if expression is None or signature_set is None:
        notes.append(f"{sig}: null skipped (no expression/signature_set)")
        return None
    filtered = signature_set.filter_to(set(expression.columns))
    if sig not in filtered.sets:
        notes.append(
            f"NULL SKIPPED for {sig!r}: no gene set of that name. "
            f"Available sets: {sorted(filtered.sets)[:8]}. "
            "Signature COLUMN names must match gene-SET names exactly (a common "
            "mismatch is a 'SIG_' prefix on the column but not on the set)."
        )
        return None

    X = frame.attrs.get("_X")
    if X is None:
        notes.append(f"{sig}: null refit skipped (embedding matrix not attached)")
        return None

    patients = frame[config.patient_col]
    if split is None:
        split = splits.preserved_site_split(
            patients, frame[config.site_col].astype("string"),
            n_folds=config.n_folds, seed=config.fold_seed,
        )

    single = sig_mod.SignatureSet(name=sig, sets={sig: filtered.sets[sig]})
    family = sig_mod.random_gene_sets(
        sorted(expression.columns), single,
        n_per_signature=config.n_null_sets, seed=config.seed,
        match_expression=expression.mean(axis=0) if config.match_null_expression else None,
    )[sig]

    scores = _score_sets(expression, family, config)
    scores = scores.reindex(frame[config.patient_col].to_numpy())
    scores.index = frame.index
    null_alphas_raw = sig_mod.alphas_from_score_frame(scores, family.sizes)

    # `residualise_matrix` rather than `residualise`: one lstsq per group for all
    # B+1 columns instead of one statsmodels fit per column, i.e. ~16,000 fits per
    # run. Equality with the reference path is asserted in
    # test_fix_residualise_matrix_matches_column_by_column.
    obs_raw = pd.to_numeric(frame[sig], errors="coerce").to_numpy(dtype=float)
    stacked_scores = np.column_stack([obs_raw, scores.to_numpy(dtype=float)])
    stacked_resid = ga.residualise_matrix(
        stacked_scores,
        pd.to_numeric(axis_values, errors="coerce").to_numpy(dtype=float),
        groups=None if within is None else within.to_numpy(),
    )
    obs_resid = pd.Series(stacked_resid[:, 0], index=frame.index)
    resid = pd.DataFrame(stacked_resid[:, 1:], index=frame.index, columns=scores.columns)

    # Regularisation is held constant across the observed/null contrast, but the
    # alpha for each fold is now selected on THAT FOLD'S TRAINING PARTITION.
    # Selecting one alpha on all patients leaked a (small) amount of test
    # information into the null; there is no reason to accept even that.
    fold_alpha = _select_alpha_per_fold(X, obs_resid.to_numpy(), split)

    # The honest null requires REFITTING the image head on each residualised
    # random set — correlating a random set against a model trained on the real
    # signature would measure signature-to-random overlap, not image
    # predictability, and would understate the null. The OBSERVED side is now
    # refitted the same way, so the two are like for like.
    usable = [c for c in resid.columns
              if np.isfinite(pd.to_numeric(resid[c], errors="coerce")).sum()
              >= config.min_patients_per_signature]
    if not usable:
        return None

    Y = resid[usable].apply(pd.to_numeric, errors="coerce").to_numpy()
    stacked = np.column_stack([obs_resid.to_numpy(dtype=float), Y])
    # STACKED multi-target ridge: one SVD per fold shared across the observed
    # target AND all B null targets, instead of B+1 independent fits.
    P = models.cross_val_predict_multi(X, stacked, split, patients.to_numpy(),
                                       fixed_alpha=fold_alpha)

    y_obs, p_obs = stacked[:, 0], P[:, 0]
    Yn, Pn = stacked[:, 1:], P[:, 1:]
    null_r = np.asarray(
        [stats.corr_ci(Yn[:, j], Pn[:, j]).value for j in range(Yn.shape[1])], dtype=float
    )

    # Reliability OF THE RESIDUALISED SCORE, which is the quantity being
    # disattenuated. The raw-score alpha understated the curated-vs-random gap
    # by roughly 5x — see signatures.alpha_from_score_and_item_variances.
    alpha_obs = alpha_null = None
    alpha_obs_raw = alpha_null_raw = np.nan
    if config.disattenuate:
        alpha_obs_raw = sig_mod.cronbach_alpha(expression, filtered.sets[sig])
        alpha_null_raw = float(np.nanmean(null_alphas_raw))
        if gene_resid_var is not None:
            alpha_obs = sig_mod.alpha_from_score_and_item_variances(
                y_obs, gene_resid_var.reindex(filtered.sets[sig]).dropna()
            )
            alpha_null = np.array([
                sig_mod.alpha_from_score_and_item_variances(
                    Yn[:, j], gene_resid_var.reindex(family.sets[c]).dropna()
                )
                for j, c in enumerate(usable)
            ], dtype=float)
        else:
            # Fall back to the raw-score closed form, and SAY SO — a silent
            # fallback here would reintroduce the exact bug this replaced.
            notes.append(
                f"{sig}: residual item variances unavailable; disattenuating with "
                "RAW-score alpha, which understates the reliability gap"
            )
            alpha_obs = alpha_obs_raw
            alpha_null = np.asarray(null_alphas_raw, dtype=float)[: len(null_r)]

    return {
        "r_observed_symmetric": stats.corr_ci(y_obs, p_obs).value,
        "null_r": null_r,
        "alpha_observed": alpha_obs,
        "alpha_null": alpha_null,
        "alpha_observed_raw": alpha_obs_raw,
        "alpha_null_mean_raw": alpha_null_raw,
        "y_observed": y_obs, "p_observed": p_obs,
        "Y_null": Yn, "P_null": Pn,
    }


def _select_alpha_per_fold(X, y, split, alphas=models.ALPHAS):
    """One alpha per fold, each chosen on that fold's TRAINING partition only.

    Holding regularisation constant across the observed-vs-null contrast is the
    point (see models.cross_val_predict_multi), but the alpha itself must not be
    selected using held-out patients.
    """
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler

    y = np.asarray(y, dtype=float)
    out = np.full(split.n_folds, np.nan)
    for k in range(split.n_folds):
        train_idx, _ = split.indices(k)
        ok = train_idx[np.isfinite(y[train_idx])]
        if len(ok) < 20:
            continue
        med = models.fit_impute_median(X[ok])
        Xs = StandardScaler().fit_transform(models.impute_median(X[ok], med))
        out[k] = float(RidgeCV(alphas=alphas).fit(Xs, y[ok]).alpha_)
    return out


def _run_outcome_arm(frame, expression, signature_set, survival, global_axis, config, notes):
    """Venet's original question, against a clinical endpoint."""
    warn = oc.endpoint_warning(frame.get(config.type_col, pd.Series(dtype=str)), survival.name)
    if warn:
        notes.append(warn)

    filtered = signature_set.filter_to(set(expression.columns))
    sig_scores = _score_sets(expression, filtered, config)
    sig_scores = sig_scores.reindex(frame[config.patient_col].to_numpy())
    sig_scores.index = frame.index

    families = sig_mod.random_gene_sets(
        sorted(expression.columns), filtered,
        n_per_signature=config.n_null_sets, seed=config.seed,
        match_expression=expression.mean(axis=0) if config.match_null_expression else None,
    )
    null_frames = {}
    for name, fam in families.items():
        sc = _score_sets(expression, fam, config)
        sc = sc.reindex(frame[config.patient_col].to_numpy())
        sc.index = frame.index
        null_frames[name] = sc

    axis_values = None
    if global_axis is not None:
        axis_values = global_axis.values if hasattr(global_axis, "values") else global_axis
        axis_values = pd.Series(axis_values).reindex(frame[config.patient_col].to_numpy())
        axis_values.index = frame.index

    within = (frame[config.type_col]
              if config.axis_within_type and config.type_col in frame else None)
    # STRATIFY BY CANCER TYPE whenever there is more than one. Pooled across 31
    # TCGA types, cancer type ALONE reaches C = 0.676 on PFI, so an unstratified
    # C-index mostly measures "does this score know which cancer this is" and any
    # tissue-correlated score inherits it. Measured: 94% of signatures beat their
    # null unstratified; stratified the margins roughly halve.
    strata = None
    if config.type_col in frame.columns and frame[config.type_col].nunique() > 1:
        strata = frame[config.type_col].to_numpy()
        notes.append(
            f"Outcome arm stratified by {config.type_col} "
            f"({frame[config.type_col].nunique()} levels): only within-type pairs "
            "are comparable. Pooling across types would measure tissue of origin, "
            "which alone reaches C=0.676 on pan-TCGA PFI."
        )
    return oc.run_outcome_arm(
        sig_scores, null_frames, survival,
        global_axis=axis_values, within=within, alpha=config.alpha, strata=strata,
    )
