"""End-to-end tests on synthetic data with a KNOWN planted confound.

The point of these tests is not coverage. It is that the pipeline must
*recover a site effect that we deliberately planted*, and must *not* report one
when none exists. If it cannot do that on synthetic data, no result it produces
on TCGA can be believed.

Run:  python -m pytest pipeline/tests -q
"""

from __future__ import annotations

import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aacr27 import barcodes, experiment, models, signatures, splits, stats  # noqa: E402


# ------------------------------------------------------------- synthetic ---

def make_cohort(
    *,
    n_patients: int = 600,
    n_sites: int = 20,
    dim: int = 64,
    site_effect: float = 0.0,
    biology_effect: float = 1.0,
    seed: int = 0,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Build a synthetic cohort where we control exactly how much signal is site.

    `site_effect=0` means the signature is pure biology; the embedding predicts
    it equally well under both split schemes. `site_effect` large means the
    signature is largely a site artefact, so preserved-site CV should collapse
    while random CV stays high. That gap is the primary endpoint.
    """
    rng = np.random.default_rng(seed)

    site_ids = [f"{i:02d}" for i in range(n_sites)]
    site_of = rng.choice(site_ids, size=n_patients)
    patients = [f"TCGA-{site_of[i]}-{i:04d}" for i in range(n_patients)]

    # Latent biology, shared by the signature and the image.
    biology = rng.normal(size=n_patients)

    # Each site has a MULTI-DIMENSIONAL fingerprint in embedding space — scanner,
    # stain protocol and fixation each shift many dimensions at once. Modelling
    # this as a single scalar would make sites with similar offsets artificially
    # inseparable and would understate a confound that is, in real TCGA data,
    # detectable at 93-99% accuracy (de Jong et al., arXiv:2501.18055).
    site_fingerprint = {s: rng.normal(scale=1.0, size=dim) for s in site_ids}
    site_signal = np.vstack([site_fingerprint[s] for s in site_of])

    # A scalar per-site shift that also leaks into the transcriptomic signature,
    # standing in for site-correlated collection and fixation practice.
    site_scalar_map = dict(zip(site_ids, rng.normal(scale=1.0, size=n_sites)))
    site_vec = np.array([site_scalar_map[s] for s in site_of])

    purity = np.clip(0.5 + 0.15 * rng.normal(size=n_patients), 0.05, 0.99)

    # Embedding carries biology, site fingerprint, purity and noise.
    load_bio = rng.normal(size=dim)
    load_pur = rng.normal(size=dim)
    X = (
        np.outer(biology, load_bio)
        + site_signal
        + np.outer(purity, load_pur)
        + rng.normal(scale=1.0, size=(n_patients, dim))
    )

    signature = (
        biology_effect * biology
        + site_effect * site_vec
        + 0.3 * purity
        + rng.normal(scale=0.5, size=n_patients)
    )

    frame = pd.DataFrame(
        {
            "barcode": [f"{p}-01A-01D-0001-01" for p in patients],
            "patient_id": patients,
            "tss": site_of,
            "cancer_type": rng.choice(["LUAD", "LUSC", "BRCA"], size=n_patients),
            "purity": purity,
            "stage": rng.choice(["I", "II", "III", "IV"], size=n_patients),
            "SIG_test": signature,
        }
    )
    return frame, X


# ------------------------------------------------------------- barcodes ---

def test_barcode_parse_extracts_site():
    b = barcodes.parse("TCGA-02-0001-01C-01D-0182-01")
    assert b.tss == "02"
    assert b.participant == "0001"
    assert b.patient_id == "TCGA-02-0001"
    assert b.is_primary_tumour


def test_barcode_rejects_garbage():
    with pytest.raises(ValueError):
        barcodes.parse("not-a-barcode")
    assert barcodes.try_parse("nope") is None


def test_barcode_handles_short_form():
    b = barcodes.parse("TCGA-AB-1234")
    assert b.patient_id == "TCGA-AB-1234"
    assert b.sample is None


# --------------------------------------------------------------- splits ---

def test_preserved_site_split_is_actually_site_disjoint():
    frame, _ = make_cohort(n_patients=400, n_sites=15, seed=1)
    result = splits.preserved_site_split(
        frame["patient_id"], frame["tss"], n_folds=5, seed=0
    )
    per_site = pd.DataFrame({"f": result.fold, "s": frame["tss"]}).groupby("s")["f"].nunique()
    assert (per_site == 1).all(), "a site leaked across folds"


def test_preserved_site_split_refuses_impossible_request():
    frame, _ = make_cohort(n_patients=100, n_sites=3, seed=2)
    with pytest.raises(ValueError, match="site-disjoint"):
        splits.preserved_site_split(frame["patient_id"], frame["tss"], n_folds=5)


def test_random_split_keeps_patients_whole():
    frame, _ = make_cohort(n_patients=300, seed=3)
    dup = pd.concat([frame, frame], ignore_index=True)  # two slides per patient
    result = splits.random_patient_split(dup["patient_id"], n_folds=5, seed=0)
    per_patient = pd.DataFrame(
        {"f": result.fold, "p": dup["patient_id"]}
    ).groupby("p")["f"].nunique()
    assert (per_patient == 1).all(), "a patient spanned folds"


# ---------------------------------------------------------------- stats ---

def test_paired_delta_ci_covers_zero_for_null_difference():
    rng = np.random.default_rng(0)
    a = rng.normal(size=500)
    est = stats.paired_delta(a, a.copy(), np.arange(500), n_boot=300, seed=0)
    assert est.lo <= 0 <= est.hi
    assert abs(est.value) < 1e-9


def test_paired_delta_detects_real_shift():
    rng = np.random.default_rng(1)
    a = rng.normal(size=800)
    b = a - 0.5
    est = stats.paired_delta(a, b, np.arange(800), n_boot=300, seed=0)
    assert est.lo > 0, "failed to detect a planted 0.5 shift"


def test_bh_fdr_matches_known_case():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205])
    rejected, q = stats.bh_fdr(p, alpha=0.05)
    assert rejected[0] and rejected[1]
    assert np.all(np.diff(q) >= -1e-12), "q-values must be monotone"


def test_min_detectable_effect_shrinks_with_n():
    assert stats.min_detectable_delta_z(50, 50) > stats.min_detectable_delta_z(700, 700)


# ------------------------------------------------------------- controls ---

def test_site_control_detects_planted_site_signal():
    frame, X = make_cohort(n_patients=600, n_sites=12, site_effect=0.0, seed=4)
    out = models.site_prediction_control(
        X, frame["tss"], frame["patient_id"], min_site_n=20
    )
    assert len(out), "no sites were evaluable"
    # The embedding was built with a site component, so this must be detectable.
    assert out["auroc"].median() > 0.7, f"median AUROC {out['auroc'].median():.3f}"


def test_combat_reduces_site_detectability():
    frame, X = make_cohort(n_patients=600, n_sites=12, seed=5)
    before = models.site_prediction_control(
        X, frame["tss"], frame["patient_id"], min_site_n=20
    )
    Xc = models.combat_correct(X, frame["tss"])
    after = models.site_prediction_control(
        Xc, frame["tss"], frame["patient_id"], min_site_n=20
    )
    assert after["auroc"].median() < before["auroc"].median(), "ComBat did not reduce site signal"


# ------------------------------------------------------- the real check ---

def test_audit_recovers_planted_site_confound():
    """The load-bearing test.

    A signature that is mostly site artefact must show a LARGER random-vs-site
    performance gap than one that is pure biology. If this fails, the primary
    endpoint is not measuring what it claims to.
    """
    cfg = experiment.AuditConfig(
        n_folds=5, n_boot=200, min_patients_per_signature=100, run_combat=False
    )

    clean, Xc = make_cohort(n_patients=700, n_sites=20, site_effect=0.0,
                            biology_effect=1.5, seed=10)
    confounded, Xf = make_cohort(n_patients=700, n_sites=20, site_effect=2.5,
                                 biology_effect=0.3, seed=10)

    res_clean = experiment.run_audit(clean, Xc, ["SIG_test"], cfg, verbose=False)
    res_conf = experiment.run_audit(confounded, Xf, ["SIG_test"], cfg, verbose=False)

    d_clean = float(res_clean.per_signature["delta_r"].iloc[0])
    d_conf = float(res_conf.per_signature["delta_r"].iloc[0])

    assert d_conf > d_clean, (
        f"site-confounded signature should degrade more under preserved-site CV "
        f"(confounded delta={d_conf:.3f} vs clean delta={d_clean:.3f})"
    )


def test_audit_runs_and_reports_baselines():
    frame, X = make_cohort(n_patients=500, n_sites=15, site_effect=1.0, seed=11)
    cfg = experiment.AuditConfig(n_boot=100, min_patients_per_signature=100)
    res = experiment.run_audit(frame, X, ["SIG_test"], cfg, verbose=False)

    row = res.per_signature.iloc[0]
    for baseline in ("r_site_only", "r_purity_only", "r_covariates"):
        assert baseline in row.index, f"missing baseline {baseline}"
    assert res.primary is not None
    assert isinstance(res.headline(), str)


def test_audit_rejects_duplicate_patients():
    frame, X = make_cohort(n_patients=200, seed=12)
    dup_frame = pd.concat([frame, frame], ignore_index=True)
    dup_X = np.vstack([X, X])
    with pytest.raises(ValueError, match="duplicate patients"):
        experiment.run_audit(dup_frame, dup_X, ["SIG_test"], verbose=False)


# -------------------------------------------------------------- scoring ---

def _ssgsea_reference(expr, sig, *, alpha=0.25):
    """The pre-optimisation `score_ssgsea` as the identity oracle.

    Kept HERE rather than in `signatures` so the shipped module carries no dead
    code, and so the thing being asserted against is literally the definition the
    frozen mean-z results were compared to, not a paraphrase of it.

    One deliberate departure from verbatim: `np.argsort(-r, kind="stable")`,
    matching the same change in `signatures._ssgsea_sample_tables`. Both sides
    previously used the default quicksort, which is NOT stable, so on the
    tie-heavy cases below this test pinned the optimised form to a
    PLATFORM-DEPENDENT target — it compared two implementations that shared the
    same undefined tie convention, and therefore passed on arm64 and x86_64
    while asserting different numbers on each. Pinning both to stable is what
    makes the assertion mean a fixed value. If you revert one side, revert both,
    or this test starts measuring the tie convention instead of the algebra.
    """
    ranks = expr.rank(axis=1, method="average")
    n_genes = expr.shape[1]
    results = {}

    for name, genes in sig.sets.items():
        cols = [g for g in genes if g in expr.columns]
        if not cols:
            continue
        member = expr.columns.isin(cols)
        scores = np.empty(len(expr))

        for i in range(len(expr)):
            r = ranks.iloc[i].to_numpy()
            order = np.argsort(-r, kind="stable")
            in_set = member[order]
            weights = np.abs(r[order]) ** alpha

            hits = np.cumsum(np.where(in_set, weights, 0.0))
            total_hit = hits[-1]
            misses = np.cumsum(np.where(in_set, 0.0, 1.0))
            total_miss = misses[-1]

            if total_hit == 0 or total_miss == 0:
                scores[i] = np.nan
                continue
            running = hits / total_hit - misses / total_miss
            scores[i] = float(np.sum(running) / n_genes)

        results[name] = scores

    return pd.DataFrame(results, index=expr.index)


def test_ssgsea_optimised_matches_reference_loop():
    """The fast ssGSEA must be the SAME statistic, not merely a similar one.

    `score_ssgsea` was rewritten from an O(n_samples * n_genes) walk per gene set
    into an O(n_samples * k) reduction, which is what made the scorer sensitivity
    analysis affordable (6.3 ms -> 0.004 ms per patient-set). A sensitivity
    analysis run with a subtly different ssGSEA answers a different question than
    the one it claims to, so the two forms are pinned to 1e-10 here.

    The tie-heavy case is the one that actually bites. Tied ranks are ordered by
    `np.argsort`'s tie-breaking, the cumulative sums depend on that order, and
    real log-expression has enormous tie groups at the floor. A rewrite that
    reproduces the intent but not the exact `argsort(-r)` call passes on
    continuous data and silently diverges on the real matrix.
    """
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(300)]
    samples = [f"S{i}" for i in range(40)]
    sig = signatures.SignatureSet(name="t", sets={
        "small": list(rng.choice(genes, size=8, replace=False)),
        "mid": list(rng.choice(genes, size=60, replace=False)),
        "big": list(rng.choice(genes, size=250, replace=False)),
    })

    continuous = pd.DataFrame(
        rng.normal(size=(40, 300)), index=samples, columns=genes)
    # Rounded to one decimal AND floored at zero: ~40% of entries become exact
    # zeros, so every sample has one huge tie group plus many small ones.
    tied = continuous.round(1).clip(lower=0.0)

    for label, expr in (("continuous", continuous), ("tied", tied)):
        fast = signatures.score_ssgsea(expr, sig)
        slow = _ssgsea_reference(expr, sig)
        assert list(fast.columns) == list(slow.columns), label
        assert fast.index.equals(slow.index), label
        np.testing.assert_allclose(
            fast.to_numpy(), slow.to_numpy(), rtol=0, atol=1e-10,
            err_msg=f"optimised ssGSEA diverges from the reference loop ({label})",
        )

    # A tie group spanning the WHOLE row is the degenerate corner: every gene
    # shares one rank, so the order is argsort's alone and nothing else pins it.
    flat = pd.DataFrame(np.zeros((5, 300)), index=samples[:5], columns=genes)
    np.testing.assert_allclose(
        signatures.score_ssgsea(flat, sig).to_numpy(),
        _ssgsea_reference(flat, sig).to_numpy(), rtol=0, atol=1e-10,
    )


def test_ssgsea_blocks_over_samples_without_changing_the_answer():
    """Sample blocking is a memory bound, not a change of estimand.

    Peak memory is held down by materialising the per-sample rank tables a block
    at a time. Ranks are computed within a row, so blocking is exact — but only
    as long as nothing in the reduction reaches across samples, which is what
    this pins.
    """
    rng = np.random.default_rng(3)
    genes = [f"G{i}" for i in range(120)]
    expr = pd.DataFrame(rng.normal(size=(37, 120)),
                        index=[f"S{i}" for i in range(37)], columns=genes)
    sig = signatures.SignatureSet(name="t", sets={
        "a": list(rng.choice(genes, size=15, replace=False))})

    one_block = signatures.score_ssgsea(expr, sig)
    original = signatures._SSGSEA_BLOCK_CELLS
    try:
        # 120 genes / 600 cells = 5 samples per block, so 37 rows spans 8 blocks
        # with a ragged last one.
        signatures._SSGSEA_BLOCK_CELLS = 600
        many_blocks = signatures.score_ssgsea(expr, sig)
    finally:
        signatures._SSGSEA_BLOCK_CELLS = original

    np.testing.assert_allclose(
        many_blocks.to_numpy(), one_block.to_numpy(), rtol=0, atol=1e-12)


def test_ssgsea_skips_sets_with_no_genes_in_the_matrix():
    """A set naming nothing in the matrix is dropped, not scored as NaN."""
    rng = np.random.default_rng(1)
    genes = [f"G{i}" for i in range(50)]
    expr = pd.DataFrame(rng.normal(size=(10, 50)),
                        index=[f"S{i}" for i in range(10)], columns=genes)
    sig = signatures.SignatureSet(name="t", sets={
        "real": genes[:10], "absent": ["ZZZ1", "ZZZ2"]})
    out = signatures.score_ssgsea(expr, sig)
    assert list(out.columns) == ["real"]
    assert np.isfinite(out["real"]).all()


# ----------------------------------------------------------- venet null ---

def test_random_gene_sets_match_requested_sizes():
    genes = [f"G{i}" for i in range(2000)]
    fams = signatures.random_gene_sets(genes, {"SIG_A": 40}, n_per_signature=10, seed=0)
    assert len(fams["SIG_A"]) == 10
    assert all(len(v) == 40 for v in fams["SIG_A"].sets.values())


def test_expression_matched_null_reproduces_signature_bin_histogram():
    """The null must match the REAL signature's expression distribution.

    Regression test for a bug where the draw spread genes uniformly across
    expression bins, which matches a uniform distribution rather than the
    signature's own — systematically wrong for immune signatures, whose genes
    cluster in the moderate-expression range.
    """
    genes = [f"G{i}" for i in range(1000)]
    # Mean expression increases with index, so bin membership is index-ordered.
    mean_expr = pd.Series(np.arange(1000, dtype=float), index=genes)

    # A signature drawn ONLY from the top decile of expression.
    real = signatures.SignatureSet(name="t", sets={"SIG_hi": genes[900:940]})

    fams = signatures.random_gene_sets(
        genes, real, n_per_signature=20, seed=0,
        match_expression=mean_expr, n_bins=10,
    )
    drawn = fams["SIG_hi"].sets
    assert len(drawn) == 20

    for members in drawn.values():
        assert len(members) == 40, f"expected 40 genes, got {len(members)}"
        idx = np.array([int(g[1:]) for g in members])
        # Every gene must come from the top decile, mirroring the real signature.
        assert idx.min() >= 900, (
            f"expression matching broken: drew gene index {idx.min()} "
            "from outside the signature's own expression bin"
        )


def test_size_only_null_ignores_expression_bins():
    """Without a gene list, matching must degrade to size-only and warn."""
    genes = [f"G{i}" for i in range(1000)]
    mean_expr = pd.Series(np.arange(1000, dtype=float), index=genes)
    with pytest.warns(UserWarning, match="SIZE-ONLY"):
        fams = signatures.random_gene_sets(
            genes, {"SIG_hi": 40}, n_per_signature=5, seed=0,
            match_expression=mean_expr,
        )
    spread = [int(g[1:]) for v in fams["SIG_hi"].sets.values() for g in v]
    assert min(spread) < 900, "size-only draw should span the whole pool"


def test_null_percentile_is_sane():
    draws = np.linspace(0, 1, 101)
    assert signatures.null_percentile(0.5, draws) == pytest.approx(50, abs=2)
    assert signatures.null_percentile(1.5, draws) == pytest.approx(100, abs=1)


def test_signature_filter_drops_uncovered_sets():
    sig = signatures.SignatureSet(name="t", sets={"big": [f"G{i}" for i in range(20)],
                                                  "tiny": ["G0", "ZZZ"]})
    filtered = sig.filter_to({f"G{i}" for i in range(20)}, min_genes=5)
    assert "big" in filtered.sets
    assert "tiny" not in filtered.sets


# -------------------------------------------------- Control C: label side ---

def test_label_site_variance_detects_site_driven_labels():
    """A label built from site must show high r2_site_only; a clean one must not.

    This is Control C's smoke test — the analysis is image-free, so it must work
    on the signature columns alone.
    """
    from aacr27.decomposition import label_site_variance

    frame, _ = make_cohort(n_patients=800, n_sites=16, seed=20)
    rng = np.random.default_rng(0)

    site_map = {s: rng.normal() for s in frame["tss"].unique()}
    frame["SIG_site_driven"] = (
        3.0 * frame["tss"].map(site_map) + rng.normal(scale=0.3, size=len(frame))
    )
    frame["SIG_clean_label"] = rng.normal(size=len(frame))

    out = label_site_variance(
        frame, ["SIG_site_driven", "SIG_clean_label"], min_site_n=10
    )
    got = out.set_index("signature")["r2_site_only"]

    assert got["SIG_site_driven"] > 0.7, f"site-driven label scored {got['SIG_site_driven']:.3f}"
    assert got["SIG_clean_label"] < 0.2, f"clean label scored {got['SIG_clean_label']:.3f}"


def test_label_site_variance_refuses_degenerate_cohort():
    from aacr27.decomposition import label_site_variance

    frame, _ = make_cohort(n_patients=60, n_sites=2, seed=21)
    with pytest.raises(ValueError, match="site"):
        label_site_variance(frame, ["SIG_test"], min_site_n=500)


# ------------------------------------------- global axis / corrected estimand ---

def _expr_with_global_axis(n=300, n_genes=400, seed=0):
    """Expression where a dominant latent axis drives ALL genes, plus a small
    immune-specific module that is independent of it."""
    rng = np.random.default_rng(seed)
    axis = rng.normal(size=n)                      # composition / purity axis
    immune = rng.normal(size=n)                    # immune-specific biology
    loadings = rng.uniform(0.6, 1.0, size=n_genes)  # every gene loads on the axis
    E = np.outer(axis, loadings) + rng.normal(scale=0.5, size=(n, n_genes))
    genes = [f"G{i}" for i in range(n_genes)]
    # 30 genes additionally carry immune-specific signal
    for j in range(30):
        E[:, j] += 1.5 * immune
    idx = [f"P{i:04d}" for i in range(n)]
    return pd.DataFrame(E, index=idx, columns=genes), pd.Series(axis, index=idx), pd.Series(immune, index=idx)


def test_global_axis_recovers_dominant_latent_factor():
    from aacr27 import globalaxis as ga
    expr, true_axis, _ = _expr_with_global_axis()
    axis = ga.compute_global_axis(expr, method="pc1")
    r = abs(np.corrcoef(axis.values.to_numpy(), true_axis.to_numpy())[0, 1])
    assert r > 0.9, f"PC1 correlated only {r:.3f} with the planted axis"
    assert axis.variance_explained > 0.3


def test_residualise_removes_the_axis():
    from aacr27 import globalaxis as ga
    expr, true_axis, _ = _expr_with_global_axis()
    axis = ga.compute_global_axis(expr, method="pc1")
    scores = pd.Series(expr.iloc[:, :30].mean(axis=1), index=expr.index)
    resid = ga.residualise(scores, axis)
    before = abs(np.corrcoef(scores, axis.values)[0, 1])
    after = abs(np.corrcoef(resid, axis.values)[0, 1])
    assert before > 0.5, "planted axis should dominate the raw score"
    assert after < 1e-8, f"residual still correlates {after:.4f} with the axis"


def test_three_way_report_separates_composition_from_specific_signal():
    """The arm that discriminates H1 from H2."""
    from aacr27 import globalaxis as ga
    expr, true_axis, immune = _expr_with_global_axis()
    axis = ga.compute_global_axis(expr, method="pc1")

    sig = pd.Series(expr.iloc[:, :30].mean(axis=1), index=expr.index)
    # An "image" that reads ONLY composition, not immune biology.
    rng = np.random.default_rng(1)
    image_axis_only = pd.Series(
        true_axis.to_numpy() + rng.normal(scale=0.3, size=len(expr)), index=expr.index
    )
    out = ga.three_way_report(image_axis_only, sig, axis)
    # The planted signature is roughly half composition and half immune, so an
    # axis-only image lands near 0.5 on the RAW score. The exact value is a
    # property of the synthetic mix; the load-bearing assertions are the two
    # residual ones below.
    assert out["r_raw"] > 0.45, "composition-reading image should predict the raw signature"
    assert abs(out["r_residual"]) < 0.25, (
        f"residual should collapse for a composition-only image, got {out['r_residual']:.3f}"
    )

    # An image that ALSO reads the immune module.
    image_both = pd.Series(
        true_axis.to_numpy() + immune.to_numpy() + rng.normal(scale=0.3, size=len(expr)),
        index=expr.index,
    )
    out2 = ga.three_way_report(image_both, sig, axis)
    assert out2["r_residual"] > out["r_residual"] + 0.2, (
        "an image carrying immune signal must show a larger residual correlation"
    )


def test_excess_over_null_is_stable_when_percentile_is_a_hair_trigger():
    """On a near-degenerate null a percentile flips wildly; the z-excess does not."""
    from aacr27 import globalaxis as ga
    from aacr27 import signatures as sg
    null = np.full(100, 0.60) + np.random.default_rng(0).normal(scale=0.002, size=100)

    hi = ga.excess_over_null(0.61, null, n=900)
    lo = ga.excess_over_null(0.59, null, n=900)
    assert hi.value > 0 and lo.value < 0
    # Percentile saturates; the z-excess stays proportionate and bounded.
    assert sg.null_percentile(0.61, null) > 99
    assert sg.null_percentile(0.59, null) < 1
    assert abs(hi.value) < 0.1, "z-excess should stay small for a small r difference"


def test_null_degeneracy_flags_a_concentrated_null():
    from aacr27 import globalaxis as ga
    tight = np.full(100, 0.6) + np.random.default_rng(0).normal(scale=0.001, size=100)
    wide = np.random.default_rng(0).uniform(0.1, 0.9, size=100)
    assert ga.null_degeneracy(tight)["degenerate"] == 1.0
    assert ga.null_degeneracy(wide)["degenerate"] == 0.0


def test_expected_random_set_correlation_matches_published_values():
    """k*rho/(1+(k-1)*rho) — the reason the null must be concentrated."""
    from aacr27 import globalaxis as ga
    assert ga.expected_random_set_correlation(29, 0.10) == pytest.approx(0.763, abs=0.01)
    assert ga.expected_random_set_correlation(100, 0.10) == pytest.approx(0.917, abs=0.01)
    assert ga.expected_random_set_correlation(50, 0.20) == pytest.approx(0.926, abs=0.01)


# ------------------------------------------------- outcome arm (PFI / Venet) ---

def _survival_from_score(score, *, seed=0, strength=1.0, censor=0.4):
    """Right-censored times whose hazard depends on `score` with known strength."""
    from aacr27 import outcome as oc
    rng = np.random.default_rng(seed)
    lp = strength * (np.asarray(score) - np.mean(score)) / (np.std(score) + 1e-9)
    t = rng.exponential(scale=np.exp(-lp))          # higher score -> shorter time
    c = rng.exponential(scale=np.quantile(t, 1 - censor) * 2)
    return oc.Survival(time=np.minimum(t, c) + 1e-6, event=(t <= c).astype(float), name="PFI")


def test_concordance_index_recovers_known_direction():
    from aacr27 import outcome as oc
    rng = np.random.default_rng(0)
    score = rng.normal(size=600)
    surv = _survival_from_score(score, strength=1.5)
    c = oc.concordance_index(score, surv)
    assert c > 0.65, f"C-index {c:.3f} too low for a strong planted effect"
    # Reversing the score must mirror C about 0.5.
    c_rev = oc.concordance_index(-score, surv)
    assert c_rev == pytest.approx(1 - c, abs=0.02)


def test_concordance_index_is_half_for_noise():
    from aacr27 import outcome as oc
    rng = np.random.default_rng(1)
    score = rng.normal(size=600)
    surv = _survival_from_score(rng.normal(size=600), strength=1.0, seed=2)
    assert oc.concordance_index(score, surv) == pytest.approx(0.5, abs=0.06)


def test_cox_score_recovers_hazard_direction_without_dichotomising():
    from aacr27 import outcome as oc
    rng = np.random.default_rng(3)
    score = rng.normal(size=800)
    surv = _survival_from_score(score, strength=1.0, seed=3)
    res = oc.cox_score(score, surv)
    assert res["hr"] > 1.2, f"HR {res['hr']:.3f} should exceed 1 per SD"
    assert res["n_events"] > 50
    assert np.isfinite(res["c_index"])


def test_survival_rejects_bad_endpoint():
    from aacr27 import outcome as oc
    with pytest.raises(ValueError, match="endpoint"):
        oc.Survival(np.array([1.0]), np.array([1.0]), name="NOPE")


def test_outcome_null_detects_a_genuinely_predictive_signature():
    """A signature driving the hazard must beat random sets that do not."""
    from aacr27 import outcome as oc
    rng = np.random.default_rng(4)
    n = 700
    real = rng.normal(size=n)
    surv = _survival_from_score(real, strength=1.2, seed=4)

    idx = [f"P{i:04d}" for i in range(n)]
    obs = pd.Series(real, index=idx)
    nulls = pd.DataFrame(
        {f"r{j}": rng.normal(size=n) for j in range(60)}, index=idx
    )
    out = oc.outcome_null(obs, nulls, surv)
    assert out["beats_null"], f"excess {out['excess_z']:.3f} [{out['excess_lo']:.3f}, {out['excess_hi']:.3f}]"
    assert out["observed"] > out["null_mean"]
    assert out["endpoint"] == "PFI"


def test_outcome_null_does_not_fire_on_a_null_signature():
    """Venet's actual finding: a signature no better than random must not pass."""
    from aacr27 import outcome as oc
    rng = np.random.default_rng(5)
    n = 700
    driver = rng.normal(size=n)
    surv = _survival_from_score(driver, strength=1.2, seed=5)

    idx = [f"P{i:04d}" for i in range(n)]
    obs = pd.Series(rng.normal(size=n), index=idx)          # unrelated to hazard
    nulls = pd.DataFrame({f"r{j}": rng.normal(size=n) for j in range(60)}, index=idx)
    out = oc.outcome_null(obs, nulls, surv)
    assert not out["beats_null"], (
        f"a signature unrelated to outcome should not beat random "
        f"(excess {out['excess_z']:.3f}, lo {out['excess_lo']:.3f})"
    )


def test_endpoint_warning_fires_for_os_not_pfi():
    from aacr27 import outcome as oc
    types = pd.Series(["LUAD"] * 10 + ["LUSC"] * 10)
    assert oc.endpoint_warning(types, "PFI") is None
    msg = oc.endpoint_warning(types, "OS")
    assert msg and "TCGA-CDR" in msg


def test_concordance_index_exact_hand_computed():
    """Exact value on a tiny case, computed by hand.

    times  [1, 2, 3]   events [1, 1, 0]   scores [3, 2, 1]
    Comparable pairs (the shorter-time member must have an event):
      (0,1): t 1<2, e0=1 -> comparable; s 3>2 -> concordant
      (0,2): t 1<3, e0=1 -> comparable; s 3>1 -> concordant
      (1,2): t 2<3, e1=1 -> comparable; s 2>1 -> concordant
    3 of 3 concordant -> C = 1.0
    """
    from aacr27 import outcome as oc
    # concordance_index guards at n>=10, so replicate the 3-row pattern 5x. Exact
    # ties introduced by replication score 0.5 each, so the ceiling is just under 1.
    big_t = np.array([1.0, 2.0, 3.0] * 5)
    big_e = np.array([1.0, 1.0, 0.0] * 5)
    big_s = np.array([3.0, 2.0, 1.0] * 5)
    # Replicating the pattern preserves perfect concordance except for exact ties,
    # which score 0.5 each; assert the direction is still perfect-or-near.
    c = oc.concordance_index(big_s, oc.Survival(big_t, big_e))
    assert c > 0.9, f"expected near-perfect concordance, got {c:.4f}"

    # And the anti-concordant ordering must mirror it.
    c_rev = oc.concordance_index(-big_s, oc.Survival(big_t, big_e))
    assert c_rev == pytest.approx(1 - c, abs=1e-9), "C must be symmetric under score reversal"


def test_concordance_index_returns_nan_when_no_events():
    from aacr27 import outcome as oc
    surv = oc.Survival(np.arange(1.0, 21.0), np.zeros(20))
    assert np.isnan(oc.concordance_index(np.arange(20.0), surv))


# ------------------------------------------- PRIMARY endpoint, end-to-end ---

def _cohort_with_expression(n=400, n_genes=2000, n_sites=14, dim=48, immune_visible=True, seed=0):
    """Cohort where the image may or may not carry immune-specific signal.

    A dominant global axis drives every gene. A 30-gene immune module carries
    extra, independent signal. `immune_visible` controls whether the embedding
    sees that module — which is exactly what the primary endpoint must detect.

    GENE-POOL SIZE IS PART OF THE FIXTURE'S CORRECTNESS, not a performance knob.
    At the original n_genes=300 the 30-gene signature was 10% of the pool, so a
    size-matched "random" draw contained 5.5 of the 30 real signature genes on
    average — each null draw was 18% real signal. That is not a null. It made the
    null artificially predictable and, once reliabilities were computed on the
    residualised scores, pushed every draw past the |r/sqrt(alpha)| >= 1 guard so
    the primary came out undefined. Real data is 163 genes in a pool of 41,046
    (0.40%); n_genes=2000 puts the fixture at 1.5%, the same regime.
    """
    rng = np.random.default_rng(seed)
    site_ids = [f"{i:02d}" for i in range(n_sites)]
    site_of = rng.choice(site_ids, size=n)
    patients = [f"TCGA-{site_of[i]}-{i:04d}" for i in range(n)]

    immune = rng.normal(size=n)

    # THE COVARIANCE STRUCTURE IS PART OF THE FIXTURE'S CORRECTNESS.
    #
    # The original generator was a SINGLE all-positive factor plus noise. Three
    # things follow from that, all of them wrong, and together they made the
    # corrected primary endpoint impossible to exercise:
    #
    #   PC1 explained 48-72% of variance   (real NSCLC, within type: 7.9%)
    #   PC1 WAS the gene-wise mean, r=1.00 (real: 0.81)
    #   mean pairwise gene-gene correlation after removing PC1: -0.0002
    #                                       (real: +0.020)
    #
    # The last one is fatal. If PC1 is exactly the grand mean, residualising it
    # forces the average residual covariance to zero by construction, so a random
    # gene set's residualised reliability is ~0 — below the floor where
    # Spearman's correction is defined. Every null draw was then dropped and the
    # disattenuated primary came out undefined, for a reason that was an artefact
    # of the generator rather than a property of the method or the data.
    #
    # Real transcriptomes are not one factor. They have a pervasive positive
    # component AND several heterogeneous, mixed-sign programmes of comparable
    # size, so PC1 is a blend that removes neither completely. Reproducing that
    # gives residualised random-set alpha 0.56-0.65 here against 0.67-0.74
    # measured on NSCLC — the same regime, which is what a fixture owes you.
    # WHAT THIS FIXTURE STILL DOES NOT REPRODUCE, stated so nobody reads more
    # into a passing test than it earns. Measured against real NSCLC:
    #
    #                      fixture        real NSCLC (within-type PC1)
    #   PC1 VE             0.083-0.088    0.079          matched
    #   r(axis, gene mean) 0.38-0.54      0.808          under
    #   residual rho       ~0.000         +0.017         NOT matched
    #   null alpha_resid   0.05-0.09      0.749          NOT matched
    #
    # Residual co-expression is the one that resists: with a handful of Gaussian
    # factors, PC1 is a true leading eigenvector and residualising it drives mean
    # residual covariance to zero. Real data keeps +0.017 because thousands of
    # genes carry structure no single component absorbs. So this fixture sits in
    # a MORE degenerate reliability regime than the data — roughly 60% of null
    # draws are dropped by the |r/sqrt(alpha)| >= 1 guard, against 0% on NSCLC.
    # That direction is fine for a test (it stresses the guard harder than
    # reality), but it means the fixture cannot validate the SIZE of a
    # disattenuation correction, only that the machinery runs and discriminates.
    #
    # SEVERAL COMPARABLE PROGRAMMES, none dominant. A single dominant factor does
    # not work no matter how its loadings are drawn: after z-scoring, one factor
    # is collinear with the gene-wise mean (measured r > 0.99 for every variant
    # tried), and residualising the mean forces average residual covariance to
    # zero. Real NSCLC has PC1 at r = 0.44 with the gene-wise mean, which is only
    # reachable with many factors of similar size.
    E = sum(
        np.outer(rng.normal(size=n), rng.normal(scale=0.35, size=n_genes))
        for _ in range(5)
    ) + rng.normal(scale=1.0, size=(n, n_genes))

    # The axis is whatever PC1 of the generated matrix turns out to be — it is
    # not one of the terms above and cannot be chosen in advance. The signature's
    # loading ON the axis is therefore injected AFTER the axis exists, which is
    # what makes the H2 branch testable: a signature that is partly composition.
    background = (E - E.mean(axis=0)) / E.std(axis=0)
    axis = np.linalg.svd(background - background.mean(axis=0), full_matrices=False)[0][:, 0]
    # A left singular vector has unit NORM (sd ~ 1/sqrt(n) = 0.05 here), which
    # the embedding's noise term below would swamp.
    axis = axis / axis.std()
    if np.corrcoef(axis, background.mean(axis=1))[0, 1] < 0:
        axis = -axis                       # match compute_global_axis's sign convention

    E[:, :30] += 1.2 * axis[:, None] + 1.0 * immune[:, None]

    genes = [f"G{i}" for i in range(n_genes)]
    expr = pd.DataFrame(E, index=patients, columns=genes)

    load_axis = rng.normal(size=dim)
    X = np.outer(axis, load_axis) + rng.normal(scale=0.8, size=(n, dim))
    if immune_visible:
        X = X + np.outer(immune, rng.normal(size=dim))

    frame = pd.DataFrame({
        "patient_id": patients,
        "tss": site_of,
        "cancer_type": rng.choice(["LUAD", "LUSC"], size=n),
        "purity": np.clip(0.5 + 0.15 * rng.normal(size=n), 0.05, 0.99),
        "stage": rng.choice(["I", "II"], size=n),
        "SIG_immune": expr.iloc[:, :30].mean(axis=1).to_numpy(),
    })
    sigset = signatures.SignatureSet(name="t", sets={"SIG_immune": genes[:30]})
    return frame, X, expr, sigset


def test_run_audit_computes_the_primary_when_axis_supplied():
    """The primary must exist, and must fire only when immune signal is visible."""
    from aacr27 import globalaxis as ga

    # n_null_sets must leave >= 5 draws AFTER disattenuation drops the ones with
    # |r/sqrt(alpha)| >= 1. This fixture drops roughly 60% (see the generator's
    # docstring), so 12 draws left exactly 5 — passing on the guard boundary.
    cfg = experiment.AuditConfig(
        n_folds=4, n_boot=200, n_null_sets=60, min_patients_per_signature=100,
        run_combat=False,
    )

    frame, X, expr, sigset = _cohort_with_expression(immune_visible=True, seed=30)
    axis = ga.compute_global_axis(expr, method="pc1")
    res = experiment.run_audit(
        frame, X, ["SIG_immune"], cfg, expression=expr, signature_set=sigset,
        global_axis=axis, verbose=False,
    )

    assert res.primary_excess is not None, "primary was not computed"
    assert res.three_way is not None and len(res.three_way) == 1
    row = res.three_way.iloc[0]
    assert row["r_axis"] > 0.3, "image should track the global axis"
    assert row["r_residual"] > 0.15, (
        f"immune-visible image should retain residual signal, got {row['r_residual']:.3f}"
    )
    assert "immune-specific excess" in res.headline()


def test_run_audit_primary_is_none_without_axis():
    """No axis -> no primary, and the omission must be recorded, not silent."""
    cfg = experiment.AuditConfig(n_boot=100, min_patients_per_signature=100, run_combat=False)
    frame, X, _, _ = _cohort_with_expression(seed=31)
    res = experiment.run_audit(frame, X, ["SIG_immune"], cfg, verbose=False)
    assert res.primary_excess is None
    assert any("PRIMARY NOT COMPUTED" in n for n in res.notes)
    assert "not computed" in res.headline()


def test_residual_collapses_when_image_cannot_see_immune_module():
    """The discriminating case: composition-only image -> residual near zero."""
    from aacr27 import globalaxis as ga
    cfg = experiment.AuditConfig(
        n_folds=4, n_boot=100, n_null_sets=60, min_patients_per_signature=100,
        run_combat=False,
    )
    frame, X, expr, sigset = _cohort_with_expression(immune_visible=False, seed=32)
    axis = ga.compute_global_axis(expr, method="pc1")
    res = experiment.run_audit(
        frame, X, ["SIG_immune"], cfg, expression=expr, signature_set=sigset,
        global_axis=axis, verbose=False,
    )
    row = res.three_way.iloc[0]
    assert row["r_raw"] > 0.3, "composition-only image should still predict the RAW signature"
    assert abs(row["r_residual"]) < 0.2, (
        f"residual should collapse without immune signal, got {row['r_residual']:.3f}"
    )


# ------------------------------------------ reliability / disattenuation ---

def _reliability_artefact(n=600, n_pool=600, k=9, load_hi=0.75, load_lo=0.25, seed=0):
    """A pool where every gene measures the SAME latent factor, at different precision.

    gene_j = L_j * latent + sqrt(1 - L_j^2) * noise_j

    Under this model the score of a k-gene set has corr(score, latent) = sqrt(alpha),
    so a HIGH-loading ("curated") set and a LOW-loading ("random") set are measuring
    the identical construct — the curated one simply measures it more precisely.

    The image sees only `latent`. So the TRUE construct-level correlation is the same
    for both sets, and any raw advantage of the curated set is pure reliability
    artefact. Disattenuation must equalise them.
    """
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(n_pool)]
    latent = rng.normal(size=n)

    loadings = np.full(n_pool, load_lo)
    loadings[:k] = load_hi                      # the "curated" module
    E = (
        loadings[None, :] * latent[:, None]
        + np.sqrt(1 - loadings[None, :] ** 2) * rng.normal(size=(n, n_pool))
    )

    idx = [f"P{i:04d}" for i in range(n)]
    expr = pd.DataFrame(E, index=idx, columns=genes)
    image = pd.Series(latent + rng.normal(scale=0.4, size=n), index=idx)
    return expr, genes[:k], image


def _coexpressed_vs_random(**kw):
    return _reliability_artefact(**kw)


def test_cronbach_alpha_separates_curated_from_random_sets():
    from aacr27 import signatures as sg
    expr, curated, _ = _reliability_artefact()
    rng = np.random.default_rng(1)
    random_set = list(rng.choice(expr.columns[9:], size=9, replace=False))

    a_cur = sg.cronbach_alpha(expr, curated)
    a_rnd = sg.cronbach_alpha(expr, random_set)
    assert a_cur > 0.6, f"co-expressed module alpha too low: {a_cur:.3f}"
    assert a_rnd < a_cur - 0.2, (
        f"random set alpha {a_rnd:.3f} should be well below curated {a_cur:.3f}"
    )


def test_mean_pairwise_correlation_matches_alpha_identity():
    """alpha == k*rho/(1+(k-1)*rho) — the identity the whole argument rests on."""
    from aacr27 import globalaxis as ga
    from aacr27 import signatures as sg
    expr, curated, _ = _reliability_artefact()
    rho = sg.mean_pairwise_correlation(expr, curated)
    alpha = sg.cronbach_alpha(expr, curated)
    assert alpha == pytest.approx(
        ga.expected_random_set_correlation(len(curated), rho), rel=1e-6
    )


def test_disattenuation_removes_the_reliability_false_positive():
    """THE load-bearing test for this correction.

    The image carries NO immune-specific signal. The curated set nonetheless
    correlates with the image prediction more strongly than random sets do,
    purely because it is a more reliable scale. Uncorrected, that reads as a
    positive result. Disattenuated, it must not.
    """
    from aacr27 import globalaxis as ga
    from aacr27 import signatures as sg

    expr, curated, image = _reliability_artefact(seed=7)
    rng = np.random.default_rng(7)

    def score(genes):
        z = (expr[genes] - expr[genes].mean()) / expr[genes].std()
        return z.mean(axis=1)

    r_obs = float(np.corrcoef(image, score(curated))[0, 1])
    a_obs = sg.cronbach_alpha(expr, curated)

    r_null, a_null = [], []
    for _ in range(60):
        rs = list(rng.choice(expr.columns[len(curated):], size=len(curated), replace=False))
        r_null.append(float(np.corrcoef(image, score(rs))[0, 1]))
        a_null.append(sg.cronbach_alpha(expr, rs))
    r_null, a_null = np.array(r_null), np.array(a_null)

    # 1. The uncorrected contrast fires — this is the false positive.
    naive = ga.excess_over_null(r_obs, r_null, n=len(image))
    assert naive.lo > 0, (
        "setup invalid: the reliability artefact should produce a spurious "
        f"positive, got excess {naive.value:.3f} [{naive.lo:.3f}, {naive.hi:.3f}]"
    )

    # 2. Disattenuated, the advantage should largely disappear.
    corrected = ga.excess_over_null(
        r_obs, r_null, n=len(image),
        reliability_observed=a_obs, reliability_null=a_null,
    )
    assert corrected.value < naive.value, (
        f"disattenuation must shrink the excess: {corrected.value:.3f} "
        f"vs naive {naive.value:.3f}"
    )
    assert abs(corrected.value) < 0.5 * abs(naive.value), (
        f"disattenuation should remove most of the artefact: corrected "
        f"{corrected.value:.3f} vs naive {naive.value:.3f}"
    )
    assert "disattenuated" in corrected.method


def test_disattenuate_guards_both_invalid_regimes():
    from aacr27 import globalaxis as ga
    # Unusable reliability -> NaN, not an exploded value.
    assert np.isnan(ga.disattenuate(0.5, 0.0))
    assert np.isnan(ga.disattenuate(0.5, -0.1))
    # Perfect reliability is a no-op.
    assert ga.disattenuate(0.5, 1.0) == pytest.approx(0.5)
    # Valid correction stays in range.
    assert ga.disattenuate(0.5, 0.5) == pytest.approx(0.7071, abs=1e-4)
    # r/sqrt(alpha) == 1 exactly is REJECTED: Fisher z of 1 is unbounded, so a
    # handful of such draws would dominate any null mean.
    assert np.isnan(ga.disattenuate(0.5, 0.25))
    assert np.isnan(ga.disattenuate(0.9, 0.5))


def test_excess_reports_dropped_null_draws_rather_than_hiding_them():
    """No silent caps: dropped draws must be visible in the method string."""
    from aacr27 import globalaxis as ga
    r_null = np.array([0.50, 0.52, 0.48, 0.51, 0.49, 0.95])
    a_null = np.array([0.60, 0.60, 0.60, 0.60, 0.60, 0.10])  # last -> |r/sqrt(a)|>1
    est = ga.excess_over_null(0.7, r_null, n=400,
                              reliability_observed=0.8, reliability_null=a_null)
    assert "dropped" in est.method, est.method
    assert "1/6" in est.method, est.method


# ================= REGRESSION TESTS FOR THE rev.3 PIPELINE FIXES =================

def test_fix21_effective_tests_identity_returns_k():
    """#21: I(|lambda| >= 1), not (lambda > 1). Identity returned 0.0 before."""
    # Independent variables: Meff == k.
    assert stats.effective_tests(np.eye(29)) == pytest.approx(29.0)
    assert stats.effective_tests(np.eye(5)) == pytest.approx(5.0)
    # Perfectly correlated block: eigenvalues are [k, 0, ..., 0], so Meff == 1.
    assert stats.effective_tests(np.ones((5, 5))) == pytest.approx(1.0, abs=1e-9)
    # Partially correlated sits strictly between the two extremes.
    C = np.full((10, 10), 0.5)
    np.fill_diagonal(C, 1.0)
    assert 1.0 < stats.effective_tests(C) < 10.0


def test_fix21_effective_tests_raises_on_nan():
    """#21: silently coercing a missing correlation to 0 asserts independence."""
    M = np.eye(4)
    M[0, 1] = M[1, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        stats.effective_tests(M)


def test_fix17_impute_median_does_not_mutate_caller():
    """#17: np.asarray does not copy float64, so the old code rewrote the input."""
    X = np.array([[1.0, 2.0], [np.nan, 4.0]])
    before = X.copy()
    out = models.impute_median(X)
    assert np.array_equal(X, before, equal_nan=True), "caller's array was mutated"
    assert np.isfinite(out).all()


def test_fix17_impute_accepts_training_fold_medians():
    """#17/#15: test folds must be imputed with TRAINING medians."""
    train = np.array([[1.0], [3.0], [5.0]])
    med = models.fit_impute_median(train)
    test = np.array([[np.nan]])
    assert models.impute_median(test, med)[0, 0] == pytest.approx(3.0)


def test_fix26_r2_ladder_shares_one_complete_case_mask():
    """#26: rungs fitted on different samples make the 'increment' meaningless."""
    from aacr27.decomposition import decompose_signature
    rng = np.random.default_rng(0)
    n = 400
    f = pd.DataFrame({
        "tss": rng.choice([f"{i:02d}" for i in range(10)], size=n),
        "cancer_type": rng.choice(["A", "B"], size=n),
        "stage": rng.choice(["I", "II"], size=n),
        "purity": rng.uniform(.2, .9, size=n),
        "SIG": rng.normal(size=n), "pred": rng.normal(size=n)})
    f.loc[f.index[:150], "purity"] = np.nan
    res = decompose_signature(f, "SIG", image_pred_col="pred")
    assert res.table["n_used"].nunique() == 1, (
        f"rungs used different n: {sorted(res.table['n_used'].unique())}"
    )
    assert "adj_r2" in res.table.columns


def test_fix30_degenerate_partition_is_refused_not_returned():
    """#30: an empty fold makes a CV estimate silently undefined."""
    # One site holds almost everything -> no balanced 5-way site-disjoint split.
    sites = pd.Series(["00"] * 480 + [f"{i:02d}" for i in range(1, 21)])
    pats = pd.Series([f"TCGA-{s}-{i:04d}" for i, s in enumerate(sites)])
    with pytest.raises(RuntimeError, match="degenerate"):
        splits.preserved_site_split(pats, sites, n_folds=5, seed=0)


def test_fix31_repeated_partitions_expose_partition_variance():
    """#31: a single partition hides how much the answer depends on the draw."""
    frame, _ = make_cohort(n_patients=700, n_sites=25, seed=1)
    reps = splits.repeated_preserved_site_splits(
        frame["patient_id"], frame["tss"], n_folds=5, repeats=5, seed=0)
    assert len(reps) >= 2
    assignments = {tuple(r.fold.tolist()) for r in reps}
    assert len(assignments) > 1, "repeats produced identical partitions"
    v = splits.partition_variance([float(r.balance["n"].std()) for r in reps])
    assert v["n_partitions"] == len(reps)


def test_fix15_combat_unestimable_for_held_out_sites():
    """#15: under preserved-site CV, ComBat cannot be applied to a new site."""
    frame, _ = make_cohort(n_patients=600, n_sites=20, seed=2)
    split = splits.preserved_site_split(frame["patient_id"], frame["tss"], n_folds=5, seed=0)
    diag = models.combat_estimable_fraction(frame["tss"], split)
    assert diag["estimable_fraction"] == 0.0, (
        "a held-out site has no training data, so ComBat must be unestimable"
    )
    # Under patient-level folds it is estimable, which is the contrast.
    rsplit = splits.random_patient_split(frame["patient_id"], n_folds=5, seed=0)
    assert models.combat_estimable_fraction(frame["tss"], rsplit)["estimable_fraction"] > 0.9


def test_fix32_embedding_layer_selection():
    """#32: the real parquet stores 14 layers x 768, not a flat vector."""
    from aacr27.data import _stack_embeddings
    col = pd.Series([np.arange(14 * 4, dtype=float).reshape(14, 4) for _ in range(6)])
    last = _stack_embeddings(col, layer=-1)
    first = _stack_embeddings(col, layer=0)
    assert last.shape == (6, 4) and first.shape == (6, 4)
    assert not np.allclose(last, first), "layer selection had no effect"
    with pytest.raises(IndexError):
        _stack_embeddings(col, layer=99)
    # A flat vector column must still work.
    flat = pd.Series([np.arange(4, dtype=float) for _ in range(3)])
    assert _stack_embeddings(flat).shape == (3, 4)


def test_fix19_constant_predictor_baselines_exist():
    """#19: the sanity floor a model must beat."""
    frame, _ = make_cohort(n_patients=200, seed=3)
    for kind in ("cohort_mean", "type_mean"):
        M = models.build_baseline_matrix(frame, kind)
        assert M.shape[0] == len(frame)
    with pytest.raises(ValueError, match="unknown baseline"):
        models.build_baseline_matrix(frame, "nonsense")


def test_fix9_null_draw_count_permits_bh_rejection():
    """#9: B=100 makes the declared BH family logically incapable of rejection."""
    cfg = experiment.AuditConfig()
    m, q = 29, 0.05
    smallest_p = 1.0 / (cfg.n_null_sets + 1)
    assert smallest_p <= q / m, (
        f"B={cfg.n_null_sets} gives min p={smallest_p:.4f}, but BH needs "
        f"<= {q/m:.6f} at m={m}"
    )


def test_fix10_gene_pool_ordering_is_deterministic():
    """#10: list(set) differs per process; the seed must reproduce the draw."""
    genes = {f"G{i}" for i in range(500)}
    a = signatures.random_gene_sets(sorted(genes), {"S": 20}, n_per_signature=3, seed=0)
    b = signatures.random_gene_sets(sorted(genes), {"S": 20}, n_per_signature=3, seed=0)
    assert a["S"].sets == b["S"].sets


# ===================== WEEK-1 VERIFICATION REGRESSION TESTS =====================

def test_week1_ancestry_asian_label_is_not_dropped():
    """The real file spells East Asian 'ASIAN', not 'EAS'.

    An earlier ANCESTRY_GROUPS listing 'EAS' routed all 633 Asian-ancestry
    patients into ADMIXED, silently destroying the subgroup the ancestry arm is
    powered on. This pins the mapping.
    """
    from aacr27.data import _broad_ancestry, ANCESTRY_VERIFIED_COUNTS
    assert _broad_ancestry("ASIAN") == "ASIAN"
    assert _broad_ancestry("asian") == "ASIAN"
    # Aliases must canonicalise, not fall through to ADMIXED.
    assert _broad_ancestry("EAS") == "ASIAN"
    assert _broad_ancestry("SAS") == "ASIAN"
    for g in ("EUR", "AFR", "AMR"):
        assert _broad_ancestry(g) == g
    assert _broad_ancestry("") == "UNKNOWN"
    assert _broad_ancestry("EUR,AFR") == "ADMIXED"
    # Verified counts are recorded so a future file change is detectable.
    assert ANCESTRY_VERIFIED_COUNTS["ASIAN"] == 633
    assert sum(ANCESTRY_VERIFIED_COUNTS.values()) == 10126


def test_week1_purity_column_named_plainly_is_found():
    """The ABSOLUTE table's column is literally `purity` — no column says ABSOLUTE."""
    import tempfile
    frame = pd.DataFrame({
        "array": ["TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01"],
        "sample": ["TCGA-OR-A5J1-01A-11D-A29H-01", "TCGA-OR-A5J2-01A-11D-A29H-01"],
        "call status": ["called", "called"],
        "purity": [0.90, 0.72],
        "ploidy": [2.0, 2.1],
    })
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "TCGA_ABSOLUTE_purity.csv"
        frame.to_csv(p, index=False)
        from aacr27.data import load_purity
        out, prov = load_purity(p)
    assert len(out) == 2
    assert out["patient_id"].tolist() == ["TCGA-OR-A5J1", "TCGA-OR-A5J2"]
    assert out["purity"].max() == pytest.approx(0.90)
    # ABSOLUTE is not expression-derived, so it must NOT be flagged circular.
    assert not any("CIRCULARITY" in n for n in prov.notes)


def test_week1_cpe_is_flagged_circular():
    """CPE embeds ESTIMATE, which is expression-derived — flag it."""
    import tempfile
    frame = pd.DataFrame({"Sample ID": ["TCGA-AA-1111-01"], "CPE": [0.6]})
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "aran_consensus.csv"
        frame.to_csv(p, index=False)
        from aacr27.data import load_purity
        _, prov = load_purity(p)
    assert any("CIRCULARITY" in n for n in prov.notes)


# --------------------------------------------------------------------------
# Regressions for the four defects found by the 2026-08-18 full-code audit.
# --------------------------------------------------------------------------

def test_fix_residualise_matrix_matches_column_by_column():
    """The vectorised residualiser must equal the statsmodels path exactly."""
    from aacr27 import globalaxis as ga

    rng = np.random.default_rng(3)
    n, k = 200, 25
    a = rng.normal(size=n)
    grp = rng.integers(0, 3, size=n)
    Y = np.outer(a, rng.uniform(0.3, 0.9, size=k)) + rng.normal(size=(n, k))
    frame = pd.DataFrame(Y, columns=[f"c{i}" for i in range(k)])

    fast = ga.residualise_matrix(Y, a, groups=grp)
    slow = ga.residualise(frame, pd.Series(a), within=pd.Series(grp)).to_numpy()
    assert np.nanmax(np.abs(fast - slow)) < 1e-10

    # A single NaN must not drop the whole column.
    Y2 = Y.copy()
    Y2[0, 0] = np.nan
    out = ga.residualise_matrix(Y2, a, groups=grp)
    assert np.isfinite(out[1:, 0]).all(), "one NaN wiped out an entire column"


def test_fix_general_alpha_equals_textbook_cronbach_on_residualised_items():
    """THE BUG: residualised correlations were disattenuated with RAW-score alpha.

    The general form must agree with the textbook definition computed directly on
    the residualised item matrix, otherwise the correction is being applied with
    the wrong reliability.
    """
    from aacr27 import globalaxis as ga

    rng = np.random.default_rng(11)
    n, k = 250, 30
    axis = rng.normal(size=n)
    Z = (np.outer(rng.normal(size=n), rng.uniform(0.3, 0.8, size=k))
         + np.outer(axis, rng.uniform(0.4, 0.9, size=k))
         + rng.normal(scale=0.7, size=(n, k)))

    R = ga.residualise_matrix(Z, axis)
    textbook = k / (k - 1) * (
        1 - R.var(axis=0, ddof=1).sum() / R.sum(axis=1).var(ddof=1)
    )
    mine = signatures.alpha_from_score_and_item_variances(
        R.mean(axis=1), R.var(axis=0, ddof=1)
    )
    assert abs(mine - textbook) < 1e-12, f"{mine} != {textbook}"

    # And on RAW unit-variance items it must agree with the existing closed forms.
    Zz = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)
    frame = pd.DataFrame(Zz, columns=[f"g{i}" for i in range(k)])
    assert abs(
        signatures.alpha_from_score_and_item_variances(Zz.mean(axis=1), Zz.var(axis=0, ddof=1))
        - signatures.cronbach_alpha(frame, list(frame.columns))
    ) < 1e-8


def test_fix_raw_and_residualised_reliability_actually_differ():
    """Guards the reason the fix matters: the two alphas are not interchangeable.

    On NSCLC the raw-score gap was +0.023 and the residualised gap +0.123. If a
    refactor ever makes these agree, the correction has silently reverted.
    """
    from aacr27 import globalaxis as ga

    rng = np.random.default_rng(5)
    n, k, pool = 300, 40, 400
    axis = rng.normal(size=n)
    E = np.outer(axis, rng.uniform(0.5, 1.0, size=pool)) + rng.normal(scale=0.6, size=(n, pool))
    Z = (E - E.mean(axis=0)) / E.std(axis=0, ddof=1)
    R = ga.residualise_matrix(Z, axis)
    iv = R.var(axis=0, ddof=1)

    idx = rng.choice(pool, size=k, replace=False)
    raw = signatures.alpha_from_scores(Z[:, idx].mean(axis=1), k)
    res = signatures.alpha_from_score_and_item_variances(R[:, idx].mean(axis=1), iv[idx])
    assert raw > 0.9, f"a size-{k} set on a one-axis pool should look reliable raw, got {raw}"
    assert res < raw - 0.3, (
        f"residualising must cost a random set its reliability: raw={raw:.3f} res={res:.3f}"
    )


def test_fix_per_fold_alpha_is_selected_on_training_rows_only():
    """`fixed_alpha` accepts a per-fold sequence, and rejects a wrong-length one."""
    rng = np.random.default_rng(9)
    n, d = 200, 12
    patients = pd.Series([f"TCGA-01-{i:04d}" for i in range(n)])
    X = rng.normal(size=(n, d))
    y = X[:, 0] * 2 + rng.normal(scale=0.5, size=n)
    split = splits.random_patient_split(patients, n_folds=4, seed=0)

    per_fold = experiment._select_alpha_per_fold(X, y, split)
    assert len(per_fold) == split.n_folds
    assert np.isfinite(per_fold).all()

    P = models.cross_val_predict_multi(X, y.reshape(-1, 1), split, patients.to_numpy(),
                                       fixed_alpha=per_fold)
    assert np.isfinite(P).all()

    with pytest.raises(ValueError, match="fixed_alpha has"):
        models.cross_val_predict_multi(X, y.reshape(-1, 1), split, patients.to_numpy(),
                                       fixed_alpha=np.array([1.0, 2.0]))


def test_fix_pooled_isi_ci_is_patient_clustered_not_signature_resampled():
    """The registration says patient-clustered; the code used to resample signatures.

    A signature-resampling bootstrap cannot produce an interval when there is
    only ONE signature (every draw is identical, so the CI has zero width). A
    patient-clustered one can, which is exactly the observable difference.
    """
    from aacr27 import globalaxis as ga

    cfg = experiment.AuditConfig(
        n_folds=4, n_boot=200, n_null_sets=60, min_patients_per_signature=100,
        run_combat=False,
    )
    frame, X, expr, sigset = _cohort_with_expression(immune_visible=True, seed=30)
    axis = ga.compute_global_axis(expr, method="pc1")
    res = experiment.run_audit(frame, X, ["SIG_immune"], cfg, expression=expr,
                               signature_set=sigset, global_axis=axis, verbose=False)

    est = res.primary_excess
    assert est is not None
    assert "patient-clustered-bootstrap" in est.method
    assert est.hi > est.lo, "one signature gave a zero-width CI — still resampling signatures"
    assert any("patient-clustered bootstrap" in n for n in res.notes)
    assert any("null subsampled" in n for n in res.notes), "subsampling must be disclosed"


def test_fix_observed_and_null_heads_are_trained_the_same_way():
    """Symmetry: the ISI's observed r must come from a head refit on the RESIDUAL.

    The raw-trained number is kept as `r_residual` for the descriptive three-way
    table, but the two must be reported separately and must not be conflated.
    """
    from aacr27 import globalaxis as ga

    cfg = experiment.AuditConfig(
        n_folds=4, n_boot=100, n_null_sets=60, min_patients_per_signature=100,
        run_combat=False,
    )
    frame, X, expr, sigset = _cohort_with_expression(immune_visible=True, seed=30)
    axis = ga.compute_global_axis(expr, method="pc1")
    res = experiment.run_audit(frame, X, ["SIG_immune"], cfg, expression=expr,
                               signature_set=sigset, global_axis=axis, verbose=False)

    row = res.immune_excess.iloc[0]
    assert "r_residual" in row and "r_residual_refit" in row
    assert np.isfinite(row["r_residual_refit"])
    assert row["r_residual"] != row["r_residual_refit"], (
        "refit and raw-trained residual correlations are identical — the "
        "symmetric head is not being used"
    )
    assert "r_residual_refit" in res.three_way.columns


def test_fix_undefined_disattenuation_is_loud_and_never_silently_uncorrected():
    """If the correction cannot be defined, the primary must be None and say why.

    Substituting the uncorrected excess would reintroduce the artefact that had a
    100% false-positive rate in simulation, so it must never be promoted.
    """
    from aacr27 import globalaxis as ga

    # Force the degenerate regime: a null so small it cannot survive the guard.
    cfg = experiment.AuditConfig(
        n_folds=4, n_boot=50, n_null_sets=6, min_patients_per_signature=100,
        run_combat=False,
    )
    frame, X, expr, sigset = _cohort_with_expression(immune_visible=True, seed=30)
    axis = ga.compute_global_axis(expr, method="pc1")
    res = experiment.run_audit(frame, X, ["SIG_immune"], cfg, expression=expr,
                               signature_set=sigset, global_axis=axis, verbose=False)

    if res.primary_excess is None:
        assert any("NOT COMPUTED" in n or "UNDEFINED" in n for n in res.notes)
        row = res.immune_excess.iloc[0]
        # The uncorrected value is visible, but under a name nobody can quote.
        assert "excess_z_UNCORRECTED_do_not_report" in row
        assert row["null_draws_dropped_by_disattenuation"] > 0


# --------------------------------------------------------------------------
# The ancestry arm. Supplementary, non-inferiority, and confounded by design —
# these tests exist to make sure the confounding is SURFACED, not smoothed over.
# --------------------------------------------------------------------------

def _ancestry_cohort(n=1600, delta=0.0, confound=False, seed=0):
    """Cohort with a known per-group signal, optionally confounded with type.

    `delta` shifts the noise level for group B so its true correlation differs
    from group A by a known amount. `confound` makes group B concentrate in a
    cancer type that is intrinsically easier to predict, so a CRUDE contrast sees
    a difference that is really disease mix — which is exactly the failure the
    stratified estimator has to catch.
    """
    rng = np.random.default_rng(seed)
    group = rng.choice(["EUR", "AFR"], size=n, p=[0.7, 0.3])

    if confound:
        # AFR is 80% "EASY", EUR is 20% — the disease mix differs sharply.
        p_easy = np.where(group == "AFR", 0.8, 0.2)
        ctype = np.where(rng.random(n) < p_easy, "EASY", "HARD")
    else:
        ctype = rng.choice(["EASY", "HARD"], size=n)

    truth = rng.normal(size=n)
    # EASY tumours are intrinsically better predicted than HARD ones.
    noise_scale = np.where(ctype == "EASY", 0.5, 1.4)
    # ...and group AFR optionally gets a genuine extra handicap on top.
    noise_scale = noise_scale * np.where(group == "AFR", 1.0 + delta, 1.0)

    pred = truth + rng.normal(scale=noise_scale, size=n)
    site = np.array([f"S{i % 40:02d}" for i in range(n)])
    return truth, pred, pd.Series(group), pd.Series(ctype), pd.Series(site)


def test_cramers_v_is_calibrated_against_a_permutation_null():
    """The analytic correction is not enough on a sparse table; permutation is.

    A single draw proves nothing here — with a very sparse table the statistic is
    so noisy that any one seed can land anywhere. The property that must hold is
    distributional: under independence the permutation p-value is not
    systematically small, and under a real association it is.
    """
    from aacr27 import ancestry as am

    small_p = 0
    for seed in range(10):
        rng = np.random.default_rng(seed)
        n = 1000
        a = pd.Series(rng.choice(["x", "y"], size=n))
        many = pd.Series(rng.choice([f"L{i}" for i in range(40)], size=n))  # independent
        out = am.cramers_v(a, many)
        assert out["v_uncorrected"] >= out["v"], "correction must not inflate"
        cal = am.cramers_v_calibrated(a, many, n_perm=60, seed=seed)
        small_p += int(cal["perm_p"] < 0.05)

    assert small_p <= 2, (
        f"independent variables flagged as associated in {small_p}/10 seeds — "
        "the permutation null is not calibrated"
    )

    # A real association must be detected, and must clear the null's own spread.
    rng = np.random.default_rng(0)
    a = pd.Series(rng.choice(["x", "y"], size=600))
    perfect = a.map({"x": "p", "y": "q"})
    assert am.cramers_v(a, perfect)["v"] > 0.95
    strong = am.cramers_v_calibrated(a, perfect, n_perm=60, seed=0)
    assert strong["perm_p"] < 0.05
    assert strong["v"] > strong["v_perm_p95"]


def test_group_correlations_recovers_planted_per_group_signal():
    from aacr27 import ancestry as am

    truth, pred, group, _, _ = _ancestry_cohort(delta=1.5, seed=1)
    out = am.group_correlations(truth, pred, group, which=("EUR", "AFR"))
    r = out.set_index("group")["r"]
    assert r["EUR"] > r["AFR"], "planted handicap for AFR was not recovered"
    assert out["n"].sum() == len(group)


def test_contrast_reports_mde_and_refuses_to_overclaim_a_null():
    from aacr27 import ancestry as am

    # No planted difference, and a small group -> must NOT be called detectable,
    # and the MDE must be reported so the null is interpretable.
    truth, pred, group, _, _ = _ancestry_cohort(n=400, delta=0.0, seed=2)
    per_group = am.group_correlations(truth, pred, group, which=("EUR", "AFR"))
    out = am.contrast_vs_reference(per_group)
    assert len(out) == 1
    row = out.iloc[0]
    assert not row["detectable"]
    assert np.isfinite(row["mde_80pct"]) and row["mde_80pct"] > 0
    assert row["delta_lo"] < 0 < row["delta_hi"]


def test_stratified_contrast_kills_a_confound_the_crude_one_reports():
    """THE test for this module: disease mix must not masquerade as ancestry."""
    from aacr27 import ancestry as am

    truth, pred, group, ctype, _ = _ancestry_cohort(n=2400, delta=0.0,
                                                    confound=True, seed=3)
    per_group = am.group_correlations(truth, pred, group, which=("EUR", "AFR"))
    crude = am.contrast_vs_reference(per_group).iloc[0]
    strat = am.stratified_contrast(truth, pred, group, ctype,
                                   which=("EUR", "AFR")).iloc[0]

    # AFR is concentrated in the EASY type, so crude sees a large spurious gain.
    assert crude["delta_z"] > 0.15, (
        f"fixture failed to plant a confound: crude delta_z={crude['delta_z']:.3f}"
    )
    # Within cancer type there is no real difference, so it must collapse.
    assert abs(strat["delta_z"]) < abs(crude["delta_z"]) / 2
    assert strat["delta_lo"] < 0 < strat["delta_hi"], (
        f"stratified contrast still reports a difference: {strat['delta_z']:.3f}"
    )
    assert strat["n_strata_used"] == 2


def test_stratified_contrast_still_finds_a_real_within_stratum_difference():
    from aacr27 import ancestry as am

    truth, pred, group, ctype, _ = _ancestry_cohort(n=2400, delta=1.5,
                                                    confound=False, seed=4)
    strat = am.stratified_contrast(truth, pred, group, ctype,
                                   which=("EUR", "AFR")).iloc[0]
    assert strat["delta_z"] < 0, "a real handicap must survive stratification"
    assert strat["delta_hi"] < 0, "and must be called detectable"


def test_site_identifiability_counts_shared_sites():
    from aacr27 import ancestry as am

    # Fully segregated: no site carries both groups -> contrast unidentifiable.
    n = 600
    group = pd.Series(["EUR"] * 400 + ["AFR"] * 200)
    site = pd.Series([f"S{i % 20:02d}" for i in range(400)]
                     + [f"T{i % 10:02d}" for i in range(200)])
    ctype = pd.Series(["A"] * n)
    out = am.site_identifiability(group, site, ctype, which=("EUR", "AFR"))
    assert out.iloc[0]["n_sites_supporting"] == 0

    # Fully shared: every site carries both.
    site_shared = pd.Series([f"S{i % 10:02d}" for i in range(n)])
    out2 = am.site_identifiability(group, site_shared, ctype, which=("EUR", "AFR"))
    assert out2.iloc[0]["n_sites_supporting"] == 10
    assert out2.iloc[0]["frac_sites_supporting"] == 1.0


def test_stratified_contrast_handles_nullable_string_strata():
    """Regression: nullable dtypes made the mask return pd.NA, not False."""
    from aacr27 import ancestry as am

    truth, pred, group, ctype, _ = _ancestry_cohort(n=800, seed=5)
    ctype = ctype.astype("string")
    ctype.iloc[:20] = pd.NA
    out = am.stratified_contrast(truth, pred, group, ctype, which=("EUR", "AFR"))
    assert len(out) == 1 and np.isfinite(out.iloc[0]["delta_z"])


def test_fix_pooled_bootstrap_survives_a_single_nan_patient():
    """REGRESSION: one non-finite patient silently zeroed the whole interval.

    `_corr_columns` centres an entire column at once, so a single NaN makes that
    column's correlation NaN on every draw. The point estimate masks non-finite
    pairs itself, so the two paths disagreed: pan-TCGA gave 16 finite
    per-signature excesses and 0 of 1000 usable bootstrap draws.
    """
    rng = np.random.default_rng(21)
    n, B = 300, 40

    boot_inputs = []
    for _ in range(3):
        y = rng.normal(size=n)
        p = 0.6 * y + rng.normal(scale=0.8, size=n)
        Yn = rng.normal(size=(n, B))
        Pn = 0.05 * Yn + rng.normal(scale=1.0, size=(n, B))
        boot_inputs.append({
            "signature": "s", "y": y, "p": p, "Yn": Yn, "Pn": Pn,
            "rel_obs": 0.95, "rel_null": np.full(B, 0.90),
        })

    cfg = experiment.AuditConfig(n_boot=100, seed=0)
    clean = experiment._pooled_isi_bootstrap(
        np.array([0.3, 0.3, 0.3]), boot_inputs, cfg, [], max_null_for_boot=B)
    assert np.isfinite(clean.lo) and clean.hi > clean.lo

    # Now poison ONE patient in ONE signature, as a missing global axis would.
    boot_inputs[1]["y"][7] = np.nan
    boot_inputs[2]["Pn"][11, 3] = np.nan
    notes = []
    dirty = experiment._pooled_isi_bootstrap(
        np.array([0.3, 0.3, 0.3]), boot_inputs, cfg, notes, max_null_for_boot=B)

    assert np.isfinite(dirty.lo) and np.isfinite(dirty.hi), (
        "a single NaN patient still destroys the interval"
    )
    assert any("dropped as" in x for x in notes), "the dropped patients must be reported"
    assert any("2 of 300" in x for x in notes), f"expected 2 dropped, notes: {notes}"


# ---------------------------------------------------------------------------
# The number checker (scripts/19_check_numbers.py) and the scorer's --arms flag.
#
# These two guard PROCESS rather than statistics, and both exist because
# something went wrong that reading the code did not reveal:
#
#   * v1 of the number checker reported "62 claims OK" while checking nothing,
#     because every regex matched only the correct literal. It passed by
#     construction. The only thing that found it was injecting a wrong number.
#   * `--arms` was added, smoke-tested by hand, and never covered by the suite.
# ---------------------------------------------------------------------------

import importlib.util  # noqa: E402
import os  # noqa: E402
import subprocess  # noqa: E402

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def _load_script(name):
    """Import a `NN_name.py` script, whose module name is not a valid identifier."""
    path = SCRIPTS / name
    spec = importlib.util.spec_from_file_location(path.stem.lstrip("0123456789_"), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_number_checker_is_not_vacuous():
    """Injecting a wrong number must produce a MISMATCH.

    This is the regression guard for the defect that shipped: a checker whose
    patterns match only the right answer is worse than no checker, because it
    reports a coverage count that reads as reassurance.

    It runs the checker's own `--self-test`, which substitutes a deliberately
    wrong value for each family of claim and asserts each is caught. `--fast`
    is NOT used: the self-test needs the full authority table.
    """
    out = subprocess.run(
        [sys.executable, str(SCRIPTS / "19_check_numbers.py"), "--self-test"],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=900)
    assert out.returncode == 0, (
        "the number checker's self-test failed -- at least one pattern is "
        f"vacuous and passes by construction:\n{out.stdout[-4000:]}")
    assert "SELF-TEST PASSED" in out.stdout
    # A self-test that injected nothing would also "pass". Require that it
    # actually exercised a non-trivial number of claim families.
    assert out.stdout.count("CAUGHT") >= 12, (
        f"too few claim families exercised:\n{out.stdout[-2000:]}")


def test_number_checker_agrees_with_the_frozen_results():
    """The documents must currently agree with the frozen artefacts."""
    out = subprocess.run(
        [sys.executable, str(SCRIPTS / "19_check_numbers.py")],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=900)
    assert out.returncode == 0, f"a document has drifted:\n{out.stdout[-4000:]}"


def test_scorer_arms_flag_offers_exactly_the_four_arms():
    """`--arms` must accept every arm and reject anything else.

    Covers the flag's contract only. The partial-run BEHAVIOUR -- that a partial
    summary carries `arms_requested`, `arms_complete: false`, and
    `NOT RUN (excluded by --arms)` rows rather than raising KeyError -- is
    covered by `test_partial_arm_run_reports_status_instead_of_raising`, which
    became possible on 2026-09-05 when that logic was extracted into the pure
    `plan_arms` / `arm_status_rows`. Before the extraction this docstring said
    the behaviour was "NOT covered here ... smoke-tested by hand", because
    reaching it needed the real cohort and hours of compute.
    """
    mod = _load_script("13_scorer_sensitivity.py")
    assert set(mod.ARMS) == {
        "mean_z__uncorrected", "mean_z__disattenuated",
        "ssgsea__uncorrected", "ssgsea__disattenuated",
    }
    # Canonical order matters: `main` re-sorts the user's `--arms` through
    # `[a for a in ARMS if a in set(args.arms)]`, so ARMS is the single source
    # of the order that reaches `arms_requested` in summary.json.
    assert list(mod.ARMS) == sorted(mod.ARMS, key=lambda a: (
        mod.SCORERS.index(a.split("__")[0]),
        ["uncorrected", "disattenuated"].index(a.split("__")[1])))

    out = subprocess.run(
        [sys.executable, str(SCRIPTS / "13_scorer_sensitivity.py"),
         "--arms", "not_an_arm"],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=120)
    assert out.returncode != 0, "an unknown arm must be rejected, not ignored"
    assert "invalid choice" in out.stderr


def test_scorer_arm_selection_preserves_canonical_order():
    """A user's `--arms` order must not leak into the reported order."""
    mod = _load_script("13_scorer_sensitivity.py")
    requested = ["ssgsea__uncorrected", "mean_z__disattenuated"]
    plan = mod.plan_arms(requested)
    assert plan["selected"] == ["mean_z__disattenuated", "ssgsea__uncorrected"]
    assert not plan["complete"], "a partial run must not look complete"


def test_partial_arm_run_reports_status_instead_of_raising():
    """The partial-run path, covered without the cohort or the compute.

    This is the oldest untested behaviour in the project. Until `plan_arms` and
    `arm_status_rows` were extracted (2026-09-05) the only way to reach this
    code was a real NSCLC run of several hours, so the flag whose entire purpose
    is to make a partial run affordable could only be exercised by paying for a
    full one. The extraction is what makes the assertion below cost nothing.

    The contract, stated as the three things a reader of a partial
    `summary.json` must be able to do: see that the run was partial, see which
    arms were asked for, and tell an EXCLUDED arm from a REQUESTED-BUT-MISSING
    one -- without a KeyError and without inferring anything from an absent key.
    """
    mod = _load_script("13_scorer_sensitivity.py")
    ssgsea_only = ["ssgsea__uncorrected", "ssgsea__disattenuated"]
    plan = mod.plan_arms(ssgsea_only)

    assert plan["complete"] is False
    assert plan["selected"] == ssgsea_only
    assert plan["skipped"] == ["mean_z__uncorrected", "mean_z__disattenuated"]
    assert plan["scorers_run"] == ["ssgsea"]
    assert plan["recon_scorers"] == ["ssgsea"]
    # One scorer cannot be compared with another, and the report must say so
    # rather than emit an empty agreement section that reads as "no difference".
    assert plan["agreement_possible"] is False

    # The arm that RAN, the arm that was requested and produced nothing, and the
    # two that were never asked for, must be three distinguishable states.
    results = {"ssgsea__uncorrected": {"excess": 0.18}}
    rows = mod.arm_status_rows(plan, results)
    assert [r["arm"] for r in rows] == list(mod.ARMS), "every arm must appear"
    status = {r["arm"]: r["status"] for r in rows}
    assert status["ssgsea__uncorrected"] == "ok"
    assert status["ssgsea__disattenuated"] == "REQUESTED BUT MISSING"
    assert status["mean_z__uncorrected"] == "NOT RUN (excluded by --arms)"
    assert status["mean_z__disattenuated"] == "NOT RUN (excluded by --arms)"

    # A full run must still look complete, and agreement must be computable.
    full = mod.plan_arms(list(mod.ARMS))
    assert full["complete"] is True
    assert full["skipped"] == []
    assert full["agreement_possible"] is True
    # With no results at all, a COMPLETE plan reports every arm as requested and
    # missing -- never as excluded, which would misdescribe a failed full run.
    assert {r["status"] for r in mod.arm_status_rows(full, {})} == {
        "REQUESTED BUT MISSING"}

    # An unknown arm is rejected in the pure layer too, not only by argparse.
    with pytest.raises(ValueError, match="unknown arm"):
        mod.plan_arms(["ssgsea__uncorrected", "not_an_arm"])


def test_partial_arm_status_check_can_fail():
    """The status check above must be capable of failing.

    A test that only ever sees correct inputs proves nothing about its own
    assertions. This feeds `arm_status_rows` a plan whose selection disagrees
    with the results it is given and confirms the status actually changes --
    i.e. that the function reads its arguments rather than returning a fixed
    table that happens to match.
    """
    mod = _load_script("13_scorer_sensitivity.py")
    plan = mod.plan_arms(list(mod.ARMS))
    all_present = {a: {"excess": 0.0} for a in mod.ARMS}
    assert {r["status"] for r in mod.arm_status_rows(plan, all_present)} == {"ok"}

    # Drop one result and exactly one status must move, and only that one.
    missing_one = dict(all_present)
    del missing_one["mean_z__disattenuated"]
    rows = {r["arm"]: r["status"]
            for r in mod.arm_status_rows(plan, missing_one)}
    assert rows["mean_z__disattenuated"] == "REQUESTED BUT MISSING"
    assert [a for a, s in rows.items() if s != "ok"] == ["mean_z__disattenuated"]


def test_version_check_matches_the_lockfile():
    """The pinned versions must be enforced, not merely printed.

    Until 2026-09-03 every script PRINTED numpy/pandas/sklearn/scipy versions
    and nothing compared them to `requirements-lock.txt`. A silent version
    drift produces exactly the symptom A9 is investigating -- a number that
    moved while the code, seeds and inputs did not -- so it must be ruled out
    mechanically rather than by reading a banner.
    """
    out = subprocess.run(
        [sys.executable, str(SCRIPTS / "20_check_versions.py"), "--strict"],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, (
        "installed packages do not match requirements-lock.txt; numbers "
        f"produced here are not comparable to the frozen results:\n{out.stdout}")


def test_version_check_actually_fails_on_a_mismatch():
    """And the check must be able to go red. A green check nobody has seen fail
    is not evidence -- the same reasoning as the number checker's --self-test."""
    mod = _load_script("20_check_versions.py")
    lock = Path(__import__("tempfile").mkdtemp()) / "requirements-lock.txt"
    lock.write_text("numpy==2.9.9\npandas==2.2.3\nscipy==1.15.3\n")
    real, argv = mod.LOCK, sys.argv
    try:
        mod.LOCK = lock
        sys.argv = ["x", "--strict"]
        assert mod.main() == 1, "a wrong pin did not fail --strict"
        sys.argv = ["x"]
        assert mod.main() == 0, "a wrong pin must only warn without --strict"
    finally:
        mod.LOCK, sys.argv = real, argv


def test_frozen_artefacts_still_hash_to_their_recorded_values():
    """No frozen result may change without someone deciding that it should.

    This is the complement to the number checker, not a duplicate of it.
    `19_check_numbers.py` verifies that the DOCUMENTS agree with the artefacts,
    reading the artefacts as ground truth -- so an edited frozen artefact makes
    it re-verify every document against the edited value.

    That is not hypothetical. Measured 2026-09-04: setting
    `nsclc_v3/summary.json`'s `primary_excess.value` to 0.99 and rerunning
    `19_check_numbers.py` printed "OK: every checked claim agrees with the
    frozen results" and exited 0, because the real triples merely became
    "unmatched", which that script reports rather than fails on. The one thing
    the guardrails most want to prevent was the one thing nothing detected.
    """
    out = subprocess.run(
        [sys.executable, str(SCRIPTS / "21_provenance_manifest.py")],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, (
        "a frozen artefact no longer hashes to its recorded value. If the "
        "change was deliberate, rerun with --write and commit the manifest "
        f"diff on its own, with the reason:\n{out.stdout}")


def test_provenance_manifest_actually_fails_on_an_edited_artefact():
    """And it must be able to go red, on a copy -- never on the real results."""
    import shutil
    import tempfile

    mod = _load_script("21_provenance_manifest.py")
    tmp = Path(tempfile.mkdtemp())
    real_results, real_manifest = mod.RESULTS, mod.MANIFEST
    try:
        # A miniature stand-in for a frozen run directory, so the real
        # `results/` tree is never written to by a test.
        (tmp / "nsclc_v3").mkdir(parents=True)
        target = tmp / "nsclc_v3" / "summary.json"
        target.write_text('{"primary_excess": {"value": 0.318240180054566}}')
        mod.RESULTS, mod.MANIFEST = tmp, tmp / "PROVENANCE.json"

        assert mod.write() == 0
        assert mod.verify() == 0, "a freshly written manifest must verify"

        target.write_text('{"primary_excess": {"value": 0.99}}')
        assert mod.verify() == 1, "an edited frozen artefact was NOT detected"

        target.unlink()
        assert mod.verify() == 1, "a deleted frozen artefact was NOT detected"
    finally:
        mod.RESULTS, mod.MANIFEST = real_results, real_manifest
        shutil.rmtree(tmp, ignore_errors=True)


# ------------------------------------------- the map must match the tree ---
#
# `results/` is a mix of frozen outputs, superseded runs, diagnostics and dead
# stubs, and `results/README.md` is the only thing that says which is which.
# Mistaking a stub for a result has cost this project real time more than once.
# Nothing verified that the map covered the tree: rows were added by hand, three
# in one session, and a fourth could simply be forgotten.

def test_results_readme_classifies_every_artefact_that_exists():
    """Every directory and log in results/ must be named in results/README.md.

    The failure this prevents is silent and specific: a run finishes, writes a
    directory, nobody classifies it, and six sessions later someone reads it as
    a result. The check is deliberately crude -- it asks only whether the NAME
    appears -- because anything cleverer would need to parse prose and would be
    wrong often enough to be switched off.
    """
    results = SCRIPTS.parent / "results"
    readme = results / "README.md"
    assert readme.exists(), "results/README.md is the map; it must exist"
    text = readme.read_text(errors="replace")

    unclassified = sorted(
        p.name for p in results.iterdir()
        # Dotfiles are tooling, not artefacts. PROVENANCE.json is the manifest
        # this file's sibling test already guards, not a run output.
        if not p.name.startswith(".")
        and p.name not in ("README.md", "PROVENANCE.json")
        and p.name not in text)

    assert not unclassified, (
        "results/README.md does not mention: " + ", ".join(unclassified) +
        "\nAdd a row classifying each as FROZEN / SUPERSEDED / DIAGNOSTIC / "
        "DEAD. Rule 8 forbids deleting them, so the fix is a line in the map, "
        "never an rm.")


def test_results_readme_check_can_fail():
    """The map check must be able to go red, or it is decoration."""
    import tempfile

    tmp = Path(tempfile.mkdtemp())
    try:
        (tmp / "README.md").write_text("# map\n\n`known_run/` is FROZEN.\n")
        (tmp / "known_run").mkdir()
        (tmp / "unclassified_run").mkdir()
        text = (tmp / "README.md").read_text()
        unclassified = sorted(
            p.name for p in tmp.iterdir()
            if not p.name.startswith(".") and p.name != "README.md"
            and p.name not in text)
        assert unclassified == ["unclassified_run"], (
            f"the map check did not flag an unclassified directory: "
            f"{unclassified}")
    finally:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)


# ------------------------------------------------ ENVIRONMENT.md is live ---

_ENV_GENERATED_END = "<!-- END GENERATED"


def _strip_volatile(text: str) -> list[str]:
    """ENVIRONMENT.md's GENERATED part, minus the one row that cannot be stable.

    `render()` records the git commit it was generated at, which is by
    construction one commit behind the tree that contains it -- committing the
    file changes the hash the file claims. So an exact-equality test would be
    red on every commit, permanently, and would be switched off within a day.
    Everything else in the generated part IS stable on a given platform, and
    that is what is compared. Measured 2026-09-05: the git-commit row is the
    ONLY line that differed between the committed file and a fresh render.

    Truncated at `<!-- END GENERATED`, added 2026-09-05. ENVIRONMENT.md now
    carries a hand-authored appendix -- the second platform's table and the
    sort-stability record -- which no machine can introspect and `render()`
    therefore cannot produce. Without the marker, keeping that appendix and
    keeping this test green were mutually exclusive, and a test that forbids the
    document from telling the truth is the one that gets deleted. The marker
    itself is emitted BY `render()`, so it cannot go missing on one side only.

    `PYTHONHASHSEED` excluded 2026-09-07, and this one was a REAL defect rather
    than a nuisance. `reproduce.sh` line 50 exports `PYTHONHASHSEED=0` for
    determinism, so a `render()` called from inside `reproduce.sh --check` sees
    `0` while the committed file -- generated from an interactive shell --
    records `(unset)`. The rows are both correct; the variable describes how the
    interpreter was LAUNCHED, not the machine the file documents, which is why
    ENVIRONMENT.md's own prose already says "exported by `reproduce.sh`. NOT set
    by default in an interactive shell". Comparing it meant `reproduce.sh
    --check` could **never** exit 0 on the very platform it validates, while the
    same test passed standalone -- so the repository's one-command reproduction
    entry point, cited in the paper's Reproducibility section and in
    `pipeline/README.md`, failed by construction. It went unnoticed because
    `--check` gained the suite step in session 28 and was not run end to end
    again until 2026-09-07. The variable is still WRITTEN to the file; it is
    only excluded from the equality comparison.
    """
    head = text.split(_ENV_GENERATED_END)[0]
    return [ln for ln in head.splitlines()
            if "| git commit |" not in ln and "`PYTHONHASHSEED`" not in ln]


def test_environment_write_preserves_hand_authored_tail():
    """`--write` must not eat the sections below the END GENERATED marker.

    REGRESSION GUARD for a bug measured 2026-09-05. `--write` was
    `out.write_text(render())`, which silently deleted every hand-authored
    section below the marker -- in practice the whole "Sort stability" record of
    A9 and A10, recovered only because the file was in git. The marker itself
    survived, because `render()` emits it, so any check that looked for the
    marker reported success while the content it protects was gone. This asserts
    on the CONTENT, not the marker.
    """
    import tempfile

    mod = _load_script("23_export_environment.py")
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "ENVIRONMENT.md"
        sentinel = "SENTINEL-hand-authored-line-that-no-render-produces"
        out.write_text(
            "stale generated part\n"
            f"{mod.MARKER} -- everything below this line is hand-authored. -->\n"
            f"\n## A hand-written section\n\n{sentinel}\n")

        spliced = mod._splice(mod.render(), out)
        assert sentinel in spliced, "the hand-authored tail was destroyed"
        assert "stale generated part" not in spliced, (
            "the generated head must be refreshed, not preserved")
        # And the check must be able to fail: with no marker there is nothing
        # to preserve, so the tail is legitimately dropped.
        out.write_text(f"no marker here\n{sentinel}\n")
        assert sentinel not in mod._splice(mod.render(), out)


def test_environment_md_matches_a_fresh_render_on_this_platform():
    """The committed environment export must not have gone stale.

    SKIPPED off-platform, deliberately. The file documents the machine the
    frozen results were produced on; on any other machine a mismatch is the
    correct and expected answer, not a defect, so asserting there would make
    the suite fail for every collaborator and on HPC4.
    """
    mod = _load_script("23_export_environment.py")
    committed = SCRIPTS.parent / "ENVIRONMENT.md"
    assert committed.exists(), (
        "ENVIRONMENT.md is referenced by the paper's Reproducibility section "
        "and by results/README.md; it must exist")

    fresh = mod.render()
    text = committed.read_text(errors="replace")

    # The file writes the interpreter as "| python | 3.13.9 (CPython) |" and
    # the architecture as "| machine | arm64 |". Both must name THIS machine
    # before an equality assertion is meaningful.
    import platform as _platform

    if (f"| python | {_platform.python_version()} " not in text
            or f"| machine | {_platform.machine()} |" not in text):
        pytest.skip(
            "ENVIRONMENT.md describes a different platform than this one "
            f"(here: python {_platform.python_version()}, "
            f"{_platform.machine()}); a mismatch off-platform is expected, "
            "not a defect")

    assert _strip_volatile(text) == _strip_volatile(fresh), (
        "ENVIRONMENT.md is stale on the platform it claims to describe. "
        "Regenerate it:\n"
        "    python scripts/23_export_environment.py --write\n"
        "and commit the diff. Everything in that file is read from the "
        "running interpreter, so a difference means the environment moved.")


# ---------------------------------------------------------------------------
# The ISI regression test, and a millisecond reproducer for A9.
# Added 2026-09-06 at the author's instruction to make it "the most rigorous"
# option available, which here meant two tiers and a MEASURED tolerance rather
# than a chosen one.
#
# WHAT THIS DELIBERATELY DOES NOT DO.
#
# The obvious regression test -- assert results/nsclc_v3/summary.json still
# reads 0.318240180054566 -- would be circular. That file is one of the 90
# artefacts `21_provenance_manifest.py` already hashes with SHA-256, so a test
# re-reading it is a second and weaker check of something already checked
# exactly. Worse, it would pass forever while the ESTIMATOR rotted underneath
# it, because nothing in it runs the estimator.
#
# So this pins the estimator: a small deterministic synthetic cohort driven
# through the REAL functions the audit uses -- preserved_site_split,
# _select_alpha_per_fold, cross_val_predict_multi, corr_ci, fisher_z. It takes
# milliseconds, so it runs on every invocation of the suite.
#
# AND IT TURNED OUT TO REPRODUCE A9. Measured 2026-09-06 on both machines:
#
#     macOS 15 / arm64,  Python 3.13.9, openblas 0.3.21     1.036173699004058
#     HPC4 / x86_64,     Python 3.12.14, scipy-openblas 0.3.27  0.9752094760755513
#     relative gap                                          5.89%
#
# The real NSCLC ISI diverges by 6.87% between the same two machines. So this
# 200x40 synthetic probe reproduces the platform defect at essentially the same
# relative magnitude as the full 944x768 pipeline -- which cost 8,525 s per side
# to observe with 18_bisect_platform.py. Anyone chasing A9 further should start
# here: pin `Ridge(solver=...)` explicitly, rerun this function on both
# machines, and see which pinning makes the two numbers agree. That experiment
# is now seconds rather than hours.
#
# THE TWO TIERS.
#
#   * On the reference platform the assertion is BIT-EXACT. Seeds are fixed, the
#     platform is fixed, the arithmetic is deterministic; any tolerance there
#     would be slack for no reason, and any estimator change fails.
#
#   * Off-platform the assertion is a RELATIVE bound of 12%. The first draft of
#     this test used an absolute 0.03, reasoned from the real ISI's scale
#     (gap 0.0218, CI half-width ~0.057). That was wrong twice over: this probe's
#     value is ~1.04 rather than ~0.318, so an absolute bound means something
#     different here, and the probe's own cross-platform spread had not been
#     measured. It was then measured, and it is 5.89% -- an absolute 0.03 bound
#     would have FAILED on HPC4 for a reason that is not a regression.
#
#     12% is roughly twice the largest spread measured on either the probe
#     (5.89%) or the real estimand (6.87%), leaving headroom for a third
#     platform, and it sits below the NSCLC bootstrap CI half-width of ~17.9%
#     relative, so it cannot tolerate a change large enough to move a reported
#     conclusion.
#
# Skipping off-platform, the way the ENVIRONMENT.md test does, was the other
# option and is strictly weaker: it asserts nothing at all on the machines most
# people reproducing this work will actually use.
#
# ---------------------------------------------------------------------------
# RE-DERIVED 2026-09-05, because the reasoning above was calibrated against a
# BUG.
#
# The 5.89% "cross-platform spread" that set the 12% bound was not a platform
# difference at all. It was A9: `splits.py:177` used a non-stable sort, so the
# two architectures built different cross-validation partitions from identical
# inputs. With `kind="stable"` applied there and at `signatures.py:212`, the
# probe returns
#
#     macOS arm64   0.9752094760755506
#     HPC4  x86_64  0.9752094760755513
#     relative      6.83e-16
#
# That is the pipeline's GENUINE cross-platform floating-point noise, and the
# old bound was about 1.8e14 times looser than reality -- wide enough to admit
# essentially any regression the probe was built to catch.
#
# The new bound is 1e-12: still ~1,400x the measured 6.83e-16, so it tolerates
# a third platform, a different BLAS, and reduction-order variation without
# tolerating anything that could move a reported number. If this assertion ever
# fires at, say, 1e-3, that is a real change in the estimator and NOT a platform
# difference -- the whole point of tightening it is that the distinction is now
# meaningful.
#
# The frozen constant moves with it. 1.036173699004058 was the pre-fix macOS
# value; it is kept below as `_SYNTHETIC_ISI_PREFIX_MACOS` because it is the
# number quoted throughout 14-SCIENCE-AUDIT.md and 09-PAPER-DRAFT.md as evidence
# of the defect, and a reader tracing those numbers must be able to find it.
# It differs from the fixed value by 6.25%.
_SYNTHETIC_ISI_FROZEN = 0.9752094760755506
_SYNTHETIC_ISI_HPC4_MEASURED = 0.9752094760755513
_SYNTHETIC_ISI_PREFIX_MACOS = 1.036173699004058   # A9, before kind="stable"
_SYNTHETIC_ISI_REL_TOL_OFF_PLATFORM = 1e-12


def _synthetic_isi() -> float:
    """One pooled z-excess from the real estimator chain on fixed synthetic data.

    Built with a seeded Generator rather than loaded from disk, so this test has
    no data dependency and cannot be quietly weakened by an input file moving.
    """
    rng = np.random.default_rng(0)
    n, p, n_null = 200, 40, 8

    X = rng.standard_normal((n, p))
    # The target is genuinely predictable from X, so the excess is a real
    # positive number. A probe pinned near zero would still pass if the ridge
    # stage were deleted outright, which is the opposite of a regression test.
    beta = rng.standard_normal(p)
    obs = X @ beta + 2.0 * rng.standard_normal(n)
    nulls = rng.standard_normal((n, n_null)) + 0.15 * (X @ beta)[:, None]

    patients = pd.Series([f"P{i:04d}" for i in range(n)])
    sites = pd.Series([f"S{i % 20:02d}" for i in range(n)], dtype="string")
    split = splits.preserved_site_split(patients, sites, n_folds=5, seed=0)

    alpha = experiment._select_alpha_per_fold(X, obs, split)
    assert np.all(np.isfinite(np.asarray(alpha, dtype=float))), (
        "synthetic fold alphas are not finite -- the probe is broken, not the "
        "estimator")

    pred = models.cross_val_predict_multi(
        X, np.column_stack([obs, nulls]), split, patients.to_numpy(),
        fixed_alpha=alpha)
    r_obs = stats.corr_ci(obs, pred[:, 0]).value
    r_null = np.asarray(
        [stats.corr_ci(nulls[:, j], pred[:, j + 1]).value for j in range(n_null)],
        dtype=float)
    return float(stats.fisher_z(r_obs) - np.nanmean(stats.fisher_z(r_null)))


def _on_reference_platform() -> bool:
    """True only on the machine ENVIRONMENT.md describes.

    Reads the committed export rather than hard-coding "darwin"/"arm64", so the
    repo names its reference platform in exactly one place and this test cannot
    drift away from it.
    """
    import platform as _platform

    env = SCRIPTS.parent / "ENVIRONMENT.md"
    if not env.exists():
        return False
    text = env.read_text(errors="replace")
    return (f"| python | {_platform.python_version()} " in text
            and f"| machine | {_platform.machine()} |" in text)


def test_isi_estimator_regression():
    """Pin the estimator chain. Bit-exact on-platform, 12% relative elsewhere."""
    got = _synthetic_isi()

    if _on_reference_platform():
        assert got == _SYNTHETIC_ISI_FROZEN, (
            "The estimator chain changed on the platform the frozen results "
            f"were produced on.\n  got    {got!r}\n  frozen "
            f"{_SYNTHETIC_ISI_FROZEN!r}\n"
            "This assertion is exact by design: same seeds, same platform, "
            "deterministic arithmetic. If the change was intended, re-measure "
            "on BOTH machines and update the two constants together -- the "
            "off-platform value is what calibrates the tolerance.")
        return

    rel = abs(got - _SYNTHETIC_ISI_FROZEN) / abs(_SYNTHETIC_ISI_FROZEN)
    assert rel <= _SYNTHETIC_ISI_REL_TOL_OFF_PLATFORM, (
        f"The estimator chain moved by {rel:.2%} against the frozen reference "
        f"value, which exceeds the {_SYNTHETIC_ISI_REL_TOL_OFF_PLATFORM:.0%} "
        "off-platform bound.\n"
        f"  got            {got!r}\n"
        f"  macOS frozen   {_SYNTHETIC_ISI_FROZEN!r}\n"
        f"  HPC4 measured  {_SYNTHETIC_ISI_HPC4_MEASURED!r}\n"
        f"  pre-fix macOS  {_SYNTHETIC_ISI_PREFIX_MACOS!r} (A9, 6.25% out)\n"
        "Since the A9/A10 sort fixes the two platforms agree to 6.83e-16, so "
        "this bound is 1e-12 rather than the old 12%. A deviation at 1e-3 or "
        "above is a real change in the estimator, not a platform difference. "
        "If it lands near 6.25% instead, suspect a reverted kind=\"stable\".")


def test_isi_regression_probe_can_fail():
    """The pin must be able to go red. Otherwise it is decoration.

    Mirrors how 19_check_numbers.py --self-test and the results/README.md test
    prove themselves: inject a wrong value and confirm the comparison rejects
    it, rather than trusting that a passing assertion means anything.
    """
    got = _synthetic_isi()
    wrong = _SYNTHETIC_ISI_FROZEN * 1.5

    assert got != wrong, "sanity: the probe cannot equal a deliberately wrong value"
    rel = abs(got - wrong) / abs(wrong)
    assert rel > _SYNTHETIC_ISI_REL_TOL_OFF_PLATFORM, (
        "A 50% injected error did not exceed the off-platform tolerance, so "
        "that tolerance is too loose to catch a regression on any platform.")


# --------------------------------------------------------------------------
# A9: the cross-validation partition is not platform-invariant (CLOSED)
# --------------------------------------------------------------------------

# Measured by running scripts/24_split_determinism.py on BOTH machines:
# macOS arm64 py3.13.9 and Einstein HPC4 x86_64 py3.12.14, numpy 2.1.3.
_STABLE_FOLD_HASH = "9f31e33c2caaf1e3"


def _fold_hash(assign_fn) -> str:
    """Hash the NSCLC site->fold assignment produced by `assign_fn`."""
    import hashlib

    from aacr27 import splits as _splits

    frame = pd.read_parquet(SCRIPTS.parent / "data" / "interim"
                            / "cohort_nsclc.parquet")
    orig = _splits._assign_sites_greedy
    _splits._assign_sites_greedy = assign_fn
    try:
        fold = _splits.preserved_site_split(
            frame["patient_id"], frame["tss"].astype("string"),
            n_folds=5, seed=0).fold.to_numpy().astype(np.int64)
    finally:
        _splits._assign_sites_greedy = orig
    return hashlib.sha256(np.ascontiguousarray(fold).tobytes()).hexdigest()[:16]


@pytest.mark.skipif(
    not (Path(__file__).resolve().parents[1] / "data" / "interim"
         / "cohort_nsclc.parquet").exists(),
    reason="needs the staged NSCLC cohort")
def test_site_size_ties_exist():
    """The precondition for A9. If this ever goes false, A9 cannot recur.

    `_assign_sites_greedy` orders sites with a NON-STABLE sort, so the partition
    is only platform-dependent while sites are tied on size. Asserting the tie
    structure explicitly means a future cohort change that removes the ties is
    reported, rather than silently making the next test vacuous.
    """
    frame = pd.read_parquet(SCRIPTS.parent / "data" / "interim"
                            / "cohort_nsclc.parquet")
    per_site = frame.groupby(frame["tss"].astype("string"))["patient_id"].nunique()
    sizes = per_site.to_numpy(dtype=float)
    _, counts = np.unique(sizes, return_counts=True)
    tied = int(counts[counts > 1].sum())
    assert len(per_site) == 68, f"expected 68 NSCLC sites, got {len(per_site)}"
    assert tied == 51, (
        f"{tied} of {len(per_site)} sites are tied on size, not 51. The tie "
        "structure is what makes np.argsort's default quicksort produce a "
        "different partition per architecture (A9). If the cohort changed, "
        "re-measure on both machines.")


@pytest.mark.skipif(
    not (Path(__file__).resolve().parents[1] / "data" / "interim"
         / "cohort_nsclc.parquet").exists(),
    reason="needs the staged NSCLC cohort")
def test_stable_sort_partition_is_platform_invariant():
    """A stable tie-break makes the partition identical on every platform.

    This is A9's fix, pinned. It deliberately does NOT assert the SHIPPED
    partition hash, because that value is platform-dependent -- which is the
    finding. `24_split_determinism.py` carries the full two-machine table.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "_sd24", SCRIPTS / "24_split_determinism.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    assert _fold_hash(mod._stable_greedy) == _STABLE_FOLD_HASH, (
        "The stable-sort site-to-fold assignment no longer matches the value "
        f"{_STABLE_FOLD_HASH} measured on macOS arm64 AND HPC4 x86_64. Either "
        "the assignment logic or the cohort changed; re-measure on both "
        "machines before updating this constant.")


@pytest.mark.skipif(
    not (Path(__file__).resolve().parents[1] / "data" / "interim"
         / "cohort_nsclc.parquet").exists(),
    reason="needs the staged NSCLC cohort")
def test_partition_invariance_check_can_fail():
    """The invariance pin must be able to go red, or it is decoration.

    Feeds the assignment a deliberately different tie-break (ascending rather
    than descending size) and asserts the hash comparison rejects it.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "_sd24", SCRIPTS / "24_split_determinism.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    def _wrong_greedy(per_site, n_folds, *, seed=0, jitter=0.0):
        flipped = per_site.iloc[::-1]
        return mod._stable_greedy(flipped, n_folds, seed=seed, jitter=jitter)

    assert _fold_hash(_wrong_greedy) != _STABLE_FOLD_HASH, (
        "A deliberately different site ordering produced the SAME fold hash, "
        "so the invariance test above cannot detect a changed partition.")


def test_cohort_scope_restricts_to_the_run_and_records_the_mismatch():
    """Scope 'run' must narrow to the run; scope 'file' must expose the gap.

    The defect this pins, found 2026-09-06: `11_close_science_gaps.py` takes
    `--cohort`/`--expr` independently of `--results`, both defaulting to the
    whole pan-TCGA tables. Stages A4, A4b and A6 therefore scored all 7,168
    pan-TCGA patients regardless of the cohort named by `--results`, and a run
    launched as `--results results/nsclc_v3_stablesort` wrote 202-site
    pan-cancer numbers into an NSCLC directory under an NSCLC heading. Nothing
    in the output named the cohort, so no check could see it.

    Asserted in BOTH directions: 'run' matches, and 'file' does NOT -- a test
    that only checked 'run' would still pass if the scope were ignored.
    """
    mod = _load_script("11_close_science_gaps.py")
    cohort = [f"P{i}" for i in range(100)]        # the wide table
    run = [f"P{i}" for i in range(10)]            # what the run actually used

    run_scope = mod.cohort_scope_status(cohort, run, "run")
    assert run_scope["n_selected"] == 10, (
        "scope 'run' did not restrict the cohort to the run's patients")
    assert run_scope["matches_run"] is True
    assert run_scope["n_missing_from_cohort"] == 0

    file_scope = mod.cohort_scope_status(cohort, run, "file")
    assert file_scope["n_selected"] == 100, (
        "scope 'file' must keep the whole table -- it is the pre-fix behaviour "
        "kept only to reproduce pancancer_v3's A6 block")
    assert file_scope["matches_run"] is False, (
        "scoping to the whole 100-patient table while the run used 10 was NOT "
        "reported as a mismatch. This is the exact condition that put "
        "pan-cancer numbers in an NSCLC directory, and it must be visible.")

    # A patient the run used but the cohort file lacks must be counted, not
    # silently dropped: a shrinking cohort is a result, but never a quiet one.
    short = mod.cohort_scope_status(cohort[:5], run, "run")
    assert short["n_missing_from_cohort"] == 5
    assert short["matches_run"] is False

    with pytest.raises(ValueError):
        mod.cohort_scope_status(cohort, run, "everything")


def test_score_join_keeps_the_patient_key_when_the_indexes_differ():
    """`join_scores_to_meta` must return `patient_id` as a column.

    The defect this pins, found 2026-09-06 by running A4 for the first time
    since the cohort-scope fix shipped: `DataFrame.join(how="inner")` preserves
    the LEFT index's name only when the two indexes are IDENTICAL. Measured on
    pandas 2.2.3, the name is dropped to None for a subset, a superset, a
    partial overlap, and even for the same set in a different order. `scored` is
    built by reassigning a bare `pd.Index`, which is unnamed, so the joined
    index came back anonymous, `reset_index()` produced a column called `index`,
    and `frame["patient_id"]` raised KeyError. A4 died, A4b cascaded off A4's
    unset `cols`, and every stage swallowing its own exception meant the script
    still exited 0 and wrote a `science_gaps.json` of nulls.

    Asserted in BOTH directions, because the identical-index case CANNOT fail:
    the fixed helper keeps the key, and the naive join on the same inputs does
    NOT. Without the second assertion this test would pass against the broken
    code, which is how it got shipped.
    """
    mod = _load_script("11_close_science_gaps.py")
    # A genuine subset, which is the real pipeline's situation: not every
    # patient in the clinical table has expression. Toy fixtures with matching
    # indexes are exactly what let this through.
    meta = pd.DataFrame(
        {"tss": ["A", "B", "C"]},
        index=pd.Index(["P1", "P2", "P3"], name="patient_id"))
    scored = pd.DataFrame(
        {"SIG_X": [0.1, 0.2]}, index=pd.Index(["P1", "P2"]))  # UNNAMED
    assert scored.index.name is None, "fixture must reproduce the unnamed index"

    joined = mod.join_scores_to_meta(meta, scored)
    assert "patient_id" in joined.columns, (
        "the patient key did not survive the score join")
    assert list(joined["patient_id"]) == ["P1", "P2"]

    naive = meta.join(scored, how="inner").reset_index()
    assert "patient_id" not in naive.columns, (
        "the naive join now preserves the index name, so this test can no "
        "longer fail and is no longer evidence of anything. Re-derive the "
        "trigger against the installed pandas before trusting the helper.")

    # A duplicated key silently multiplies rows in a per-patient control, so the
    # helper asserts on content, not only on the column's presence.
    dupe = pd.DataFrame(
        {"tss": ["A", "B"]},
        index=pd.Index(["P1", "P1"], name="patient_id"))
    with pytest.raises(AssertionError):
        mod.join_scores_to_meta(dupe, scored)


def test_science_gaps_fingerprint_matches_the_writer():
    """The reader's digest must equal the writer's, computed on one input.

    `11_close_science_gaps.py` compares a digest it computes itself against one
    `AuditResult.cohort_fingerprint` wrote. Two independent implementations of
    "sha256 of the sorted unique ids" is precisely the shape that drifts -- a
    stray strip(), a different separator, sorting after stringifying rather
    than before -- and the drift would present as a permanent FATAL on every
    correct run. Asserting they agree in prose is not a measurement; this runs
    both.
    """
    from aacr27.experiment import AuditResult

    mod = _load_script("11_close_science_gaps.py")
    # Deliberately awkward: out of order, duplicated, and non-string ids, so a
    # difference in ordering or de-duplication between the two would show.
    patients = ["TCGA-B", "TCGA-A", "TCGA-B", "TCGA-C", 7, "10", "9"]
    stub = types.SimpleNamespace(
        predictions=pd.DataFrame({"patient": patients}))

    written = AuditResult.cohort_fingerprint(stub)
    assert mod.fingerprint_patients(patients) == written["patient_set_sha256"], (
        "the science-gaps reader and the summary.json writer disagree about "
        "the cohort digest, so every run carrying a fingerprint would be "
        "declared a mismatch")
    assert written["n_patients"] == 6, (
        f"fixture expected 6 unique ids, writer saw {written['n_patients']}")


def test_cohort_fingerprint_mismatch_is_detected_and_absence_is_tolerated():
    """The fingerprint check must fire on a mismatch AND stay quiet on absence.

    Both directions are required. As of 2026-09-07 no results directory on disk
    carries a fingerprint, so every real invocation takes the "absent" branch:
    a test that ran only against the real runs would exercise the case where
    the check is a guaranteed no-op and prove nothing. That is exactly how the
    A11 cohort-scope fix shipped broken -- verified only against pan-cancer,
    where the patient filter changes nothing.

    The overloaded `cohort` key is pinned too. `12_partition_variance.py` and
    `13_scorer_sensitivity.py` already write `"cohort": "nsclc"` as a bare
    string, so a check that merely tested truthiness of `summary["cohort"]`
    would crash or misread on runs that exist today.
    """
    mod = _load_script("11_close_science_gaps.py")
    patients = {"P1", "P2", "P3"}
    digest = mod.fingerprint_patients(patients)

    match = mod.cohort_fingerprint_status(
        {"cohort": {"n_patients": 3, "patient_set_sha256": digest}}, patients)
    assert match["status"] == "match", match

    # The case the check exists for: summary and predictions describe different
    # cohorts. This is A11, mechanised.
    wrong = mod.cohort_fingerprint_status(
        {"cohort": {"n_patients": 7168,
                    "patient_set_sha256": mod.fingerprint_patients(
                        {"OTHER1", "OTHER2"})}},
        patients)
    assert wrong["status"] == "mismatch", (
        "a summary describing a different cohort was not detected, so the "
        f"fingerprint is decorative: {wrong}")

    # An identical digest but a contradictory count is still incoherent.
    bad_n = mod.cohort_fingerprint_status(
        {"cohort": {"n_patients": 944, "patient_set_sha256": digest}}, patients)
    assert bad_n["status"] == "mismatch", (
        f"n_patients disagreeing with the patient set was accepted: {bad_n}")

    # Absence must be tolerated: every frozen run predates the writer.
    for summary in ({}, None, {"cohort": None},
                    {"cohort": "nsclc"},  # the overloaded string key
                    {"cohort": {"n_patients": 3}}):  # dict, but no digest
        got = mod.cohort_fingerprint_status(summary, patients)
        assert got["status"] == "absent", (
            f"{summary!r} should read as 'no fingerprint recorded', got {got}")
        assert got["observed_n"] == 3


def test_a_run_records_which_patients_it_used():
    """`summary.json` must identify the cohort, not just the settings.

    The defect this pins is A11's root cause: a results directory recorded the
    estimator's configuration and its outputs and never said whose data
    produced them, so `11_close_science_gaps.py` could write 202-site pan-TCGA
    numbers into an NSCLC directory with nothing in either file contradicting
    it. It was caught only because three values matched another cohort's to
    seventeen significant figures.

    Asserted on all three properties that make the fingerprint useful, because
    a hash that fails any one of them would still LOOK like a fingerprint:
    order-invariance (or every re-run reports a false mismatch), sensitivity to
    the patient SET (or it detects nothing), and sensitivity at CONSTANT n (or
    it is just a row count wearing a hash).
    """
    from aacr27.experiment import AuditResult

    def fp(patients):
        r = AuditResult.__new__(AuditResult)
        r.predictions = pd.DataFrame({"patient": patients})
        return r.cohort_fingerprint()

    base = fp(["P2", "P1", "P1", "P3"])          # duplicated rows, unsorted
    assert base["n_patients"] == 3, "rows were counted instead of patients"

    assert fp(["P3", "P1", "P2"]) == base, (
        "the fingerprint changed under a row reordering, so every honest "
        "re-run would report a cohort mismatch")

    # Same SIZE, different people. This is the case a bare n_patients cannot
    # see, and the one that makes the hash worth storing at all.
    other = fp(["P1", "P2", "P4"])
    assert other["n_patients"] == base["n_patients"]
    assert other["patient_set_sha256"] != base["patient_set_sha256"], (
        "two different 3-patient cohorts hashed identically -- the fingerprint "
        "is not a function of the patient set and cannot detect a swap")

    assert fp(["P1", "P2"])["patient_set_sha256"] != base["patient_set_sha256"]


def test_strict_mode_detects_a_pattern_that_matches_nothing():
    """`--strict` must go red on a vacuous pattern, and plain mode must not.

    Both halves matter. A pattern matching no line anywhere PASSES the default
    checker -- it contributes no comparison, so nothing can disagree -- and a
    green "OK: every checked claim agrees" is printed over it. Three
    consecutive sessions shipped such a pattern, and `--strict` found two more
    that had been passing silently for handoffs. Asserting only that `--strict`
    is green today would not prove it can see anything, so this injects an
    unmatchable pattern and requires the two modes to DISAGREE.
    """
    src = (SCRIPTS / "19_check_numbers.py").read_text()
    anchor = '    ("ancestry_cramers_v", r"Cram'
    assert anchor in src, "anchor for the injection moved; update this test"
    mutant = src.replace(
        anchor,
        '    ("ancestry_cramers_v", r"ZZZ_NEVER_APPEARS\\s*(\\d?\\.\\d{2,})", 5e-4,\n'
        '     "MUTANT deliberately unmatchable", None),\n' + anchor, 1)

    # The mutant must live BESIDE the real script: 19_check_numbers.py derives
    # REPO from `Path(__file__).parents[1]`, so a copy in a temp directory
    # cannot find results/ and exits 1 for the wrong reason -- which is exactly
    # what happened on the first attempt, and is why the NAME is asserted below
    # and not merely the exit code.
    alt = SCRIPTS / "_19_strict_mutant_tmp.py"
    try:
        alt.write_text(mutant)
        strict = subprocess.run(
            [sys.executable, str(alt), "--strict", "--fast"],
            cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=900)
        plain = subprocess.run(
            [sys.executable, str(alt), "--fast"],
            cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=900)
    finally:
        alt.unlink(missing_ok=True)

    assert strict.returncode == 1, (
        "--strict did not fail on a pattern that matches nothing, so it "
        f"cannot detect a vacuous claim.\n{strict.stdout[-1500:]}")
    assert "MUTANT deliberately unmatchable" in strict.stdout, (
        "--strict failed but did not NAME the vacuous pattern, which is the "
        "part that makes it actionable.")
    assert plain.returncode == 0, (
        "plain mode was expected to pass the mutant -- if it now fails, the "
        "vacuity is being caught some other way and this test is testing the "
        f"wrong thing.\n{plain.stdout[-1500:]}")


def test_no_claim_pattern_is_vacuous_on_the_real_tree():
    """`--strict` must be GREEN on the repository as it actually stands.

    The test above proves `--strict` CAN see a vacuous pattern; it proves
    nothing about this tree, because it runs against a mutant. Until this test
    existed, `--strict` was a flag someone had to remember to pass, and four
    consecutive sessions shipped or nearly shipped a pattern that matched
    nothing -- most recently a Results ssGSEA point estimate whose line-scoped
    pattern was defeated by a line wrap. Each time the suite stayed green,
    because `test_number_checker_agrees_with_the_frozen_results` runs the
    checker WITHOUT `--strict`. This closes that gap: a pattern that checks no
    text can no longer be committed under a green suite.

    `--fast` is deliberately NOT used. It trims the authority table to skip the
    live pytest collection and abstract recount, and while it was measured on
    2026-09-07 to report the same zero vacuous patterns as the full run, that
    measurement only covers a tree where nothing is vacuous. A pattern keyed to
    one of the authorities `--fast` omits is exactly the case where the two
    could diverge, so the full table is used.

    Agreement between documents and artefacts is NOT this test's business --
    that is `test_number_checker_agrees_with_the_frozen_results`. Keeping them
    apart keeps the failure message specific about which of the two broke.
    """
    out = subprocess.run(
        [sys.executable, str(SCRIPTS / "19_check_numbers.py"), "--strict"],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=900)

    # Assert on CONTENT before the exit code. A checker that died early, or one
    # whose scan silently matched nothing at all, would also exit 0 on some
    # refactors; requiring the scan line proves it really walked the documents.
    assert "checked " in out.stdout and " claims across " in out.stdout, (
        "the checker did not report a scan, so its exit status says nothing "
        f"about whether any pattern ran:\n{out.stdout[-2000:]}")
    assert "VACUOUS" not in out.stdout, (
        "at least one CLAIMS pattern matches no text anywhere in the "
        "repository. It passes by construction and checks nothing. Fix the "
        "pattern (WRAP: usually repairs a line guard broken by a line wrap) "
        f"or delete it:\n{out.stdout[-2500:]}")
    assert out.returncode == 0, (
        "19_check_numbers.py --strict failed on the real tree:\n"
        f"{out.stdout[-4000:]}")


def test_every_exempt_file_is_actually_scanned_and_the_check_can_fail():
    """An EXEMPT entry must name a file the checker would otherwise scan.

    THE DEFECT THIS PINS, found 2026-09-07. `19_check_numbers.py` applies
    `EXEMPT` inside `for rel in DOCS:` -- so an EXEMPT key that is not in DOCS
    is unreachable. Both keys (`05-PRE-REGISTRATION.md`, `02-PROJECT-DECISION
    .md`) were missing from DOCS, which made the entire table dead code. It
    survived because the OBSERVABLE BEHAVIOUR was identical: the files went
    unscanned either way. What differed was the reason, and the reason is what
    everyone was relying on -- the module docstring gives the exemption its own
    section, and handoffs carried "frozen ... and in the checker's EXEMPT table"
    forward as a live fact for eleven sessions. Nobody had run the checker and
    looked for the word `exempt` in its output; it never appeared.

    Why this matters beyond tidiness: the protection was accidental. Adding
    `05-PRE-REGISTRATION.md` to DOCS for any reason -- and it stated results, so
    it was a plausible addition -- would have silently started scanning a file
    whose test count is frozen at 67 on purpose, reported a MISMATCH against the
    live count, and invited someone to "fix" the pre-registration to match. The
    exemption existed precisely to prevent that and would not have fired.

    This is the project's recurring shape: a guard being present is not the
    guarded thing being safe.
    """
    mod = _load_script("19_check_numbers.py")

    assert mod.EXEMPT, "EXEMPT is empty; this test is guarding nothing"
    missing = sorted(set(mod.EXEMPT) - set(mod.DOCS))
    assert not missing, (
        "these EXEMPT entries name files the checker never scans, so the "
        "exemption is unreachable and the file is protected only by accident: "
        f"{missing}. Either add each to DOCS (so the exemption fires and is "
        "reported under `skipped:`) or delete the entry and say plainly in the "
        "docstring that the file is simply not scanned.")

    # CAN IT FAIL? Assert the property is discriminating rather than trivially
    # true of any two collections. Without this, a refactor that made EXEMPT a
    # subset of DOCS by construction -- or emptied one of them -- would leave a
    # test that passes while checking nothing.
    fake_exempt = dict(mod.EXEMPT)
    fake_exempt["99-NOT-IN-DOCS.md"] = "a file no DOCS entry names"
    assert sorted(set(fake_exempt) - set(mod.DOCS)) == ["99-NOT-IN-DOCS.md"], (
        "the subset check does not detect an EXEMPT key absent from DOCS, so "
        "it would pass no matter what EXEMPT contained")

    # And the exemption must be OBSERVABLE, not merely reachable: a reader of
    # the output has to be able to tell that a file was skipped on purpose
    # rather than silently dropped. This is the half that was missing -- the
    # dead table was invisible precisely because nothing printed.
    out = subprocess.run(
        [sys.executable, str(SCRIPTS / "19_check_numbers.py")],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=900)
    for name, reason in mod.EXEMPT.items():
        assert f"{name} (exempt:" in out.stdout, (
            f"{name} is in EXEMPT and in DOCS, but the checker never reports "
            f"it as exempt, so the skip is invisible to anyone reading the "
            f"output. Expected a `skipped:` line naming it. Reason on file: "
            f"{reason}\n{out.stdout[-2500:]}")


# ------------------------------------------------- the seed coupling ---
#
# Three library functions take their RNG seed from a DEFAULT ARGUMENT, and
# their only callers do not pass one:
#
#   decomposition.permutation_calibration  <- 11_close_science_gaps.py:403
#   ancestry.cramers_v_calibrated          <- 06_run_ancestry.py:93,94
#   ancestry.pooled_over_signatures        <- 06_run_ancestry.py:189,190,199
#
# They are reproducible today for one reason only: each default happens to be
# 0, and `AuditConfig.seed` also happens to be 0. Nothing enforces that. Change
# the config seed to study seed sensitivity and these three arms would silently
# keep using 0 -- not crash, not warn, just quietly fail to move, which is the
# worst of the three outcomes because the run would still look successful.
#
# Two other functions were on the suspect list carried in the handoff chain and
# are NOT defects: `signatures.random_gene_sets` is called with
# `seed=config.seed` from experiment.py, and `signatures.split_half_reliability`
# is called with `seed=seed` from 13_scorer_sensitivity.py. Their defaults are
# never reached. Checked 2026-09-07 by reading every call site.
#
# The coupling is left in place rather than rewired: rewiring means threading a
# seed through two standalone scripts, and re-running them would change
# artefacts that PROVENANCE.json hashes. So it is made LOUD instead. This test
# is what makes the coupling safe, and it is the reason it may be left alone.

_SEED_COUPLED = (
    ("aacr27.decomposition", "permutation_calibration", "11_close_science_gaps.py"),
    ("aacr27.ancestry", "cramers_v_calibrated", "06_run_ancestry.py"),
    ("aacr27.ancestry", "pooled_over_signatures", "06_run_ancestry.py"),
)


def seed_default_mismatches(config_seed, entries):
    """Pure: the entries whose default seed differs from `config_seed`.

    Pure so that BOTH directions are testable. A check that has only ever been
    run on the passing case is not a check -- this project shipped an A11 fix
    verified only where it was a no-op, and that is the mistake being avoided.
    """
    import importlib
    import inspect as _inspect

    out = []
    for mod_name, func_name, caller in entries:
        mod = importlib.import_module(mod_name)
        default = _inspect.signature(getattr(mod, func_name)).parameters["seed"].default
        if default != config_seed:
            out.append((f"{mod_name}.{func_name}", default, config_seed, caller))
    return out


def test_seeded_defaults_agree_with_the_audit_config_seed():
    """The three unwired seeds must equal `AuditConfig.seed`."""
    config_seed = experiment.AuditConfig().seed
    bad = seed_default_mismatches(config_seed, _SEED_COUPLED)
    assert not bad, (
        "a function whose caller does not pass a seed now defaults to a "
        "DIFFERENT seed than AuditConfig, so that arm would silently use the "
        "wrong RNG stream while the run still reported success:\n"
        + "\n".join(
            f"  {name}: default {d!r} != AuditConfig.seed {c!r} (called by {caller})"
            for name, d, c, caller in bad)
        + "\nFix by passing seed= explicitly at the call site, which is the "
          "real repair; changing the default only moves the coupling.")


def test_seed_default_check_can_fail():
    """Prove the check above is not vacuous, in both directions.

    Without this, `test_seeded_defaults_agree_with_the_audit_config_seed` would
    pass just as happily if `seed_default_mismatches` returned [] for every
    input -- which is exactly how a guard becomes decorative.
    """
    # Direction 1: a config seed that does NOT match the defaults must be caught.
    caught = seed_default_mismatches(4242, _SEED_COUPLED)
    assert len(caught) == len(_SEED_COUPLED), (
        "every one of the three coupled functions defaults to 0, so a config "
        f"seed of 4242 must flag all {len(_SEED_COUPLED)}; got {caught}")
    assert all(d == 0 and c == 4242 for _, d, c, _ in caught)

    # Direction 2: a matching seed must be reported clean, so the check is not
    # simply always-positive.
    assert seed_default_mismatches(0, _SEED_COUPLED) == []

    # Direction 3: a function whose default genuinely differs is detected on
    # its own, not merely as part of a sweep.
    assert seed_default_mismatches(
        0, (("aacr27.ancestry", "cramers_v_calibrated", "x"),)) == []


# ------------------------------------------------- US spelling in prose ---

# The four artefacts that are actually submitted. Internal working documents
# (14-SCIENCE-AUDIT.md, 15-CHECKLIST.md, 16-YOUR-TASKS.md, the pipeline and
# results READMEs) are deliberately NOT here -- they were excluded from the
# conversion on purpose, because they sit next to British code identifiers
# like `residualise_matrix` and mixing conventions there reads worse than
# leaving them.
_SUBMITTED_DOCS = (
    "09-PAPER-DRAFT.md",
    "07-ABSTRACT-DRAFT.md",
    "submission/COVER-LETTER.md",
    "poster/poster.html",
)

# British forms that have a US counterpart. `analyse` is here; `analysis` and
# `analyses` are NOT, because they are identical in both conventions -- a naive
# list that included them would fire on correct prose and get itself deleted,
# which is how the previous WORD_CLAIMS attempt died.
_BRITISH = (
    "tumour", "colour", "favour", "behaviour", "labour", "centre", "licence",
    "defence", "analyse", "catalogue", "programme", "artefact", "generalis",
    "specialis", "normalis", "randomis", "optimis", "summaris", "categoris",
    "characteris", "dichotomis", "residualis", "standardis", "harmonis",
    "utilis", "organis", "recognis", "minimis", "maximis", "regularis",
    "penalis", "modelling", "labelled", "labelling", "signalling",
)

# Exact tokens the stem list above matches but which are CORRECT US English.
# `analyses` contains the stem `analyse` yet is the ordinary plural of
# `analysis` and identical in both conventions; flagging it is precisely the
# false-positive-on-correct-prose that killed the previous word-claim attempt.
_NOT_BRITISH = frozenset({"analyses", "analysis"})


def british_spellings_in_prose(text, *, is_markdown=True):
    """Pure: [(line_no, word)] for British spellings in PROSE only.

    Three exclusions, each one a trap this project already hit:

    * fenced blocks and `inline code`, because the library really does contain
      `residualise_matrix` and `globalaxis.residualisation`; renaming those
      would break the frozen pipeline and the provenance manifest;
    * everything from `## References` onward, because reference titles are
      quoted verbatim and Schmauch's *Nat Commun* title really is spelled
      "tumours";
    * HTML comments, which carry drafting notes rather than submitted prose.

    Pure so both directions are testable. The 2026-09-06 conversion was driven
    by a fixed find-and-replace list and silently missed every British word not
    on that list -- four of them survived in submitted prose for a day, and one
    (*signalling*) was sitting in the very checklist line that cited it as
    evidence the manuscript had not been converted. A word list cannot notice
    its own omissions; a scan over the finished text can.
    """
    import re as _re

    pattern = _re.compile(r"\b\w*(?:" + "|".join(_BRITISH) + r")\w*\b", _re.I)
    hits, in_fence, past_refs = [], False, False
    for lineno, line in enumerate(text.split("\n"), 1):
        if line.strip().startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        if is_markdown and _re.match(r"^##\s+References", line):
            past_refs = True
        if past_refs:
            continue
        clean = _re.sub(r"`[^`]*`", "", line)
        clean = _re.sub(r"<!--.*?-->", "", clean)
        hits.extend((lineno, m.group(0)) for m in pattern.finditer(clean)
                    if m.group(0).lower() not in _NOT_BRITISH)
    return hits


def test_submitted_documents_use_american_spelling():
    """The submission artefacts must stay in US spelling.

    Decided and executed in commit `3ee560a`; see the record in
    16-YOUR-TASKS.md. This test exists because the decision was executed by a
    word list, and a word list is not a check -- `dichotomisation`,
    `signalling`, `generalises` and `generalised` all survived it.
    """
    repo = Path(__file__).resolve().parents[2]
    bad = {}
    for rel in _SUBMITTED_DOCS:
        path = repo / rel
        assert path.exists(), f"submitted document is missing: {rel}"
        hits = british_spellings_in_prose(
            path.read_text(encoding="utf-8"), is_markdown=rel.endswith(".md"))
        if hits:
            bad[rel] = hits
    assert not bad, (
        "British spelling survives in prose that gets submitted (US spelling "
        "was decided and converted in commit 3ee560a):\n"
        + "\n".join(f"  {rel}:{ln}  {w}"
                    for rel, hits in bad.items() for ln, w in hits)
        + "\nIf a hit is inside a reference title or a code identifier it "
          "should be in backticks or below `## References`, which this scan "
          "already skips -- fix the markup rather than the word list.")


def test_american_spelling_check_can_fail():
    """Prove the scan above is not vacuous, and that its exclusions work.

    Without this the test would pass just as happily if the regex matched
    nothing at all.
    """
    # Direction 1: plain British prose is caught.
    assert british_spellings_in_prose("the tumour microenvironment\n") == [
        (1, "tumour")]

    # Direction 2: correct US prose is clean, and the two words that are
    # identical in both conventions do NOT fire.
    assert british_spellings_in_prose(
        "the tumor microenvironment\nAnalyses are exploratory; one analysis.\n"
    ) == []

    # Direction 3: each documented exclusion actually excludes.
    assert british_spellings_in_prose("call `residualise_matrix` here\n") == []
    assert british_spellings_in_prose(
        "```\nresidualise_matrix(x)\n```\n") == []
    assert british_spellings_in_prose("<!-- tumour note -->\n") == []
    assert british_spellings_in_prose(
        "## References\n1. RNA-Seq expression of tumours from slides\n") == []

    # Direction 4: the exclusions are scoped, not blanket -- prose on a line
    # that ALSO contains a code span is still scanned.
    assert british_spellings_in_prose(
        "the tumour is scored by `residualise_matrix`\n") == [(1, "tumour")]

    # Direction 5: the reference cutoff is markdown-only, so the HTML poster
    # cannot silently drop everything after a line that happens to match.
    assert british_spellings_in_prose(
        "## References\ntumour\n", is_markdown=False) == [(2, "tumour")]

    # Direction 6: the `analyses` exclusion is a whole-token exemption, not a
    # licence to ignore the `analyse` family. Catching this distinction is the
    # whole reason the exclusion is an exact-token set and not another stem.
    assert british_spellings_in_prose("we analysed the cohort\n") == [
        (1, "analysed")]
    assert british_spellings_in_prose("we analyse the cohort\n") == [
        (1, "analyse")]
    assert british_spellings_in_prose("two analyses; one analysis\n") == []


# --------------------------------- cohort fingerprint, on a real directory ---

def _make_results_dir(tmp_path, patients, summary=None, name="run_v9"):
    """A minimal results directory: predictions.csv.gz, optionally summary.json.

    Deliberately built from files rather than mocks. The point of the test
    below is the plumbing -- which file is opened and which column is read --
    and a mock would assume exactly the thing under test.
    """
    import gzip
    import json as _json

    res = tmp_path / name
    res.mkdir()
    frame = pd.DataFrame({
        "patient": list(patients) * 2,
        "model": ["ridge_embedding"] * len(patients) + ["ridge_site"] * len(patients),
        "target": ["SIG_A"] * (2 * len(patients)),
        "y_true": np.arange(2 * len(patients), dtype=float),
    })
    with gzip.open(res / "predictions.csv.gz", "wt") as fh:
        frame.to_csv(fh, index=False)
    if summary is not None:
        (res / "summary.json").write_text(_json.dumps(summary))
    return res


def test_science_gaps_reads_a_real_fingerprint_from_a_real_directory(tmp_path):
    """Exercise the fingerprint WIRING, not just the pure comparison.

    `test_cohort_fingerprint_mismatch_is_detected_and_absence_is_tolerated`
    covers all four branches of `cohort_fingerprint_status`, which is pure. It
    cannot cover the part that actually broke things in A11: reading the right
    file, taking the patient set from the right column, and letting a mismatch
    stop the run. As of 2026-09-07 no results directory on disk carries a
    fingerprint, so the real script has only ever taken the "absent" branch --
    the match and mismatch paths have never executed outside a unit test.

    Backfilling fingerprints into the frozen runs was considered and DECLINED:
    it would mean teaching the reader to look in side files, i.e. adding a
    production code path that exists only for testing. A fixture directory
    gives the same coverage and adds nothing to the shipped reader.
    """
    mod = _load_script("11_close_science_gaps.py")
    patients = ["TCGA-A1-0001", "TCGA-A2-0002", "TCGA-A3-0003"]
    digest = mod.fingerprint_patients(patients)

    # --- the branch that has never run for real: a genuine match -------------
    res = _make_results_dir(
        tmp_path, patients,
        summary={"cohort": {"n_patients": 3, "patient_set_sha256": digest}})
    preds, run_patients, fp = mod.read_and_check_cohort(res)
    assert fp["status"] == "match", fp
    assert run_patients == set(patients)
    # The digest must be computed over the PREDICTIONS' patients, not copied
    # out of the summary -- otherwise the check compares a value to itself.
    assert fp["observed"] == digest
    assert fp["observed_n"] == 3
    # `preds` must come back usable: main() pivots it immediately afterwards,
    # so a reader that returned the right verdict and the wrong frame would
    # still break the stage it guards.
    assert (preds["model"] == "ridge_embedding").sum() == 3

    # --- the branch that exists for A11: predictions from another cohort -----
    other = _make_results_dir(
        tmp_path, patients,
        summary={"cohort": {"n_patients": 7168,
                            "patient_set_sha256": mod.fingerprint_patients(
                                ["TCGA-ZZ-9999"])}},
        name="run_mismatch")
    with pytest.raises(SystemExit) as excinfo:
        mod.read_and_check_cohort(other)
    message = str(excinfo.value)
    assert "FATAL" in message and "run_mismatch" in message, message
    assert "7168" in message and "3" in message, (
        "the refusal must name both cohort sizes, or the person reading it "
        f"cannot tell which directory is wrong: {message}")

    # --- and absence must still be tolerated on a real directory -------------
    legacy = _make_results_dir(tmp_path, patients, summary=None, name="run_old")
    _, _, absent = mod.read_and_check_cohort(legacy)
    assert absent["status"] == "absent", absent
    assert absent["observed_n"] == 3

    # A summary that exists but predates the writer is absent, not a mismatch.
    stringy = _make_results_dir(tmp_path, patients, summary={"cohort": "nsclc"},
                                name="run_stringkey")
    _, _, got = mod.read_and_check_cohort(stringy)
    assert got["status"] == "absent", got


def test_science_gaps_fingerprint_reader_can_fail(tmp_path):
    """Prove the fixture above would notice a reader that stopped checking.

    Direction: if the wiring took its "observed" patients from the summary
    instead of from predictions.csv.gz, every directory would match itself and
    the test above would still pass on the match case. Pin the discriminating
    property -- a directory whose predictions were changed while the summary
    was not MUST become a mismatch.
    """
    mod = _load_script("11_close_science_gaps.py")
    patients = ["TCGA-A1-0001", "TCGA-A2-0002", "TCGA-A3-0003"]
    summary = {"cohort": {"n_patients": 3,
                          "patient_set_sha256": mod.fingerprint_patients(patients)}}

    # Same summary, one patient swapped in the predictions: the count still
    # agrees, so only the digest can catch this.
    swapped = _make_results_dir(
        tmp_path, ["TCGA-A1-0001", "TCGA-A2-0002", "TCGA-B9-0009"],
        summary=summary, name="run_swapped")
    with pytest.raises(SystemExit) as excinfo:
        mod.read_and_check_cohort(swapped)
    assert "FATAL" in str(excinfo.value)

    # And the unswapped control passes, so the assertion above is not simply
    # always-true.
    clean = _make_results_dir(tmp_path, patients, summary=summary,
                              name="run_clean")
    assert mod.read_and_check_cohort(clean)[2]["status"] == "match"


def test_the_public_snapshot_can_actually_render_the_figures(tmp_path):
    """The Data Availability Statement, made mechanical instead of asserted.

    The paper promises that "derived intermediates sufficient to reproduce every
    reported figure and table are deposited with the analysis code". That
    sentence was FALSE twice, and both times the prose was checked by reading
    rather than by running:

      * until 2026-09-07 `null_draws.npz` was excluded by suffix, so figure 1B
        silently drew a mean +/- SD band -- and the script printed that it had
        drawn the real distribution;
      * until later the same day `data/interim/expr_nsclc.parquet` and the GMT
        were the only route to figure 2's panel sizes, and `pipeline/data/` is
        excluded from the snapshot on purpose. `10_make_figures.py` exited 1 at
        figure 2 for every reader reproducing from the deposit.

    Reading the snapshot builder could not have caught either one; running the
    figure script against a staged tree catches both. So that is what this does:
    stage a real snapshot into a temp dir and render from it, with NO `data/`
    reachable -- which is exactly a reader's situation.
    """
    out = tmp_path / "snap"
    staged = subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert staged.returncode == 0, f"staging failed:\n{staged.stdout[-2000:]}"
    assert not (out / "pipeline" / "data").exists(), (
        "this test is only meaningful while the snapshot excludes data/; if that "
        "changed deliberately, the reasoning here needs rewriting, not the path")

    env = dict(os.environ, MPLCONFIGDIR=str(tmp_path / "mpl"))
    rendered = subprocess.run(
        [sys.executable, str(out / "pipeline" / "scripts" / "10_make_figures.py"),
         "--outdir", str(tmp_path / "figs")],
        cwd=out / "pipeline", capture_output=True, text=True, timeout=900, env=env)
    assert rendered.returncode == 0, (
        "the public snapshot cannot render the manuscript figures, so the Data "
        f"Availability Statement is false:\n{rendered.stdout[-3000:]}"
        f"\n{rendered.stderr[-3000:]}")

    # Presence of an exit code is not enough -- figure 1B has a fallback that
    # renders a DIFFERENT panel and keeps going. Require the real null, and
    # require figure 2 to have actually plotted rather than been skipped.
    assert "panel B null: violins" in rendered.stdout, (
        "panel B fell back to a mean +/- SD band inside the snapshot, which is "
        f"the 2026-09-07 defect returning:\n{rendered.stdout[-3000:]}")
    assert "figure2: Spearman" in rendered.stdout, (
        f"figure 2 did not render from the snapshot:\n{rendered.stdout[-3000:]}")
    # Every manuscript figure, by its real stem. Named explicitly rather than
    # globbed: a glob would pass on a tree that rendered only figure 0.
    for name in ("figure0_schematic", "figure1_isi", "figure2_reliability",
                 "figure3_outcome", "figure4_controls"):
        assert (tmp_path / "figs" / f"{name}.png").exists(), f"{name}.png missing"


def test_the_snapshot_figure_check_can_fail(tmp_path):
    """Prove the test above is discriminating, not merely green.

    Direction: remove the deposited panel sizes from a staged tree -- which is
    precisely the state the snapshot was in before 2026-09-07 -- and rendering
    MUST fail loudly rather than skipping figure 2 or plotting NaN. Without this,
    a future change that made `gene_set_sizes()` return `{}` on a missing deposit
    would leave the test above passing while the defect returned.
    """
    out = tmp_path / "snap"
    subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300, check=True)
    deposit = out / "pipeline" / "results" / "gene_set_sizes.csv"
    assert deposit.exists(), "the deposit is not staged -- figure 2 cannot render"
    deposit.unlink()

    env = dict(os.environ, MPLCONFIGDIR=str(tmp_path / "mpl"))
    rendered = subprocess.run(
        [sys.executable, str(out / "pipeline" / "scripts" / "10_make_figures.py"),
         "--outdir", str(tmp_path / "figs")],
        cwd=out / "pipeline", capture_output=True, text=True, timeout=900, env=env)
    assert rendered.returncode != 0, (
        "removing the panel-size deposit did not fail the render, so the check "
        "above cannot detect the defect it exists for")
    assert "FATAL" in rendered.stdout + rendered.stderr


def test_the_public_snapshot_can_verify_its_own_provenance(tmp_path):
    """The deposit's integrity checker must not report corruption that isn't there.

    Found 2026-09-07 the same way both Data Availability defects were: by running
    a consumer inside a staged tree. `21_provenance_manifest.py` printed
    **FAILED: 0 changed, 4 missing** to any reader who tried to verify the
    deposit. Nothing was wrong with it -- the four are `.npz` probes the snapshot
    excludes on purpose. A checker that cries corruption over its own deliberate
    exclusions is worse than no checker: the next reader who sees FAILED has no
    way to tell it from a real tampering, so they learn to ignore it.

    Note what the earlier audit missed. The excluded `.npz` files WERE examined
    on 2026-09-07 and correctly cleared -- but only against the question "is this
    a figure or table input?". Nobody asked whether anything else in the deposit
    read them. Clearing a file for one consumer says nothing about the others.
    """
    out = tmp_path / "snap"
    staged = subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert staged.returncode == 0, f"staging failed:\n{staged.stdout[-2000:]}"

    script = out / "pipeline" / "scripts" / "21_provenance_manifest.py"
    assert not (out / "pipeline" / "scripts" / "22_build_public_snapshot.py").exists(), (
        "the builder excludes itself by design, and that absence is how the "
        "manifest detects it is running inside a deposit; if it is now staged, "
        "the detection needs rewriting, not this assertion")

    checked = subprocess.run(
        [sys.executable, str(script)],
        cwd=out / "pipeline", capture_output=True, text=True, timeout=300)
    assert checked.returncode == 0, (
        "the public deposit cannot verify its own provenance, so a reader is "
        f"told the artefacts are corrupt when they are not:\n{checked.stdout[-3000:]}")
    # Exit 0 alone is not enough: it would also pass if the check silently
    # tolerated every absence. Require it to NAME the excluded artefacts.
    assert "NOT DEPOSITED" in checked.stdout, (
        f"the excluded artefacts were not reported at all:\n{checked.stdout[-3000:]}")
    for name in ("bisect_macos.npz", "bisect_hpc4.npz",
                 "alpha_margin_macos.npz", "alpha_margin_fast_macos.npz"):
        assert name in checked.stdout, f"{name} not accounted for in the deposit"


def test_the_snapshot_provenance_check_can_fail(tmp_path):
    """Prove the test above discriminates rather than tolerating any absence.

    Direction: delete an artefact the snapshot DOES stage. That is real
    corruption of the deposit, and it must still fail inside the deposit, where
    the tolerance for excluded files lives.
    """
    out = tmp_path / "snap"
    subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300, check=True)

    victim = out / "pipeline" / "results" / "nsclc_v3" / "summary.json"
    assert victim.exists(), "expected a deposited artefact to delete"
    victim.unlink()

    checked = subprocess.run(
        [sys.executable, str(out / "pipeline" / "scripts" / "21_provenance_manifest.py")],
        cwd=out / "pipeline", capture_output=True, text=True, timeout=300)
    assert checked.returncode == 1, (
        "deleting a DEPOSITED artefact did not fail the deposit's provenance "
        "check, so the tolerance added for excluded files is too broad")
    assert "MISSING" in checked.stdout


def test_an_undeposited_artefact_is_still_required_on_the_working_tree():
    """The tolerance is scoped to the deposit, and that scoping must hold.

    Inside a snapshot an artefact flagged `deposited: false` may be absent. On
    the working repository it may NOT -- that is where the frozen results live,
    and `bisect_macos.npz` vanishing there is exactly the event the manifest
    exists to catch. The two cases are separated by whether the snapshot builder
    is present, so this pins the working-tree half.
    """
    import json as _json
    import shutil
    import tempfile

    mod = _load_script("21_provenance_manifest.py")
    tmp = Path(tempfile.mkdtemp())
    real_results, real_manifest = mod.RESULTS, mod.MANIFEST
    real_tracked = list(mod.TRACKED_FILES)
    try:
        (tmp / "nsclc_v3").mkdir(parents=True)
        probe = tmp / "bisect_macos.npz"          # a suffix the snapshot excludes
        probe.write_bytes(b"not really an npz")
        mod.RESULTS, mod.MANIFEST = tmp, tmp / "PROVENANCE.json"
        mod.TRACKED_FILES = ["bisect_macos.npz"]

        assert mod.write() == 0
        recorded = _json.loads((tmp / "PROVENANCE.json").read_text())["files"]
        assert recorded["bisect_macos.npz"]["deposited"] is False, (
            "a .npz that is not null_draws.npz is not staged by the snapshot; "
            "if that changed, this test's premise changed with it")

        assert not mod._is_public_deposit(), (
            "the real tree has the builder in it, so this must not read as a deposit")
        probe.unlink()
        assert mod.verify() == 1, (
            "an excluded-from-snapshot artefact vanished from the WORKING tree "
            "and was tolerated; the deposit tolerance has leaked out of scope")
    finally:
        mod.RESULTS, mod.MANIFEST = real_results, real_manifest
        mod.TRACKED_FILES = real_tracked
        shutil.rmtree(tmp, ignore_errors=True)


def test_pattern_reach_is_recorded_and_can_detect_a_lost_document():
    """The reach report must record real files, not an empty dict that prints.

    `--strict` asks only whether a pattern matched ANYTHING, so a pattern that
    used to guard three documents and now guards one still passes. That is not
    hypothetical: on 2026-09-07, rewording a checklist item turned "rebuild
    stages **173 files**" into "rebuild now stages **173**", one of two guarded
    documents silently lost its guard, and `--strict` stayed green.

    So this pins the machinery the report is built on: reach must agree with the
    match COUNT it is derived from, and a pattern known to span two documents
    must be seen spanning two. A report that always printed "1 file" would be as
    useless as no report, and nothing else would notice.
    """
    mod = _load_script("19_check_numbers.py")
    table = mod.authority(slow=False)
    counts, reach = {}, {}
    mod.scan(mod.REPO, table, mod.every_stored_run(), quiet=True,
             match_counts=counts, match_files=reach)

    assert counts, "no patterns were exercised at all"
    for desc, n in counts.items():
        files = reach.get(desc, set())
        if n:
            assert files, (
                f"{desc!r} compared {n} line(s) but recorded no document; the "
                "reach bookkeeping has drifted from what was actually checked")
        else:
            assert not files, f"{desc!r} matched nothing yet recorded {files}"

    spanning = {d for d, f in reach.items() if len(f) >= 2}
    assert spanning, (
        "not one pattern reached two documents, which cannot be true of this "
        "repository -- the per-file recording is collapsing everything into one")
    # The cohort sizes are quoted in the abstract, the manuscript, the poster and
    # the checklists. If this ever reaches one document, either the reach
    # recording broke or the repository really did lose those guards; both are
    # worth failing on.
    wide = [d for d in spanning if "cohort size" in d]
    assert wide, f"no cohort-size pattern spans documents; reach={sorted(spanning)[:5]}"


def test_compare_live_detects_a_deposit_that_is_missing_a_staged_file():
    """The staged tree and the published tree are DIFFERENT ARTIFACTS.

    Every check in this repository measured the staged one. The Data
    Availability Statement was fixed in three consecutive sessions and verified
    each time against a fresh build; nobody asked the published repository what
    it actually contains, because `git ls-remote` answers only what refs exist
    and a ref cannot reveal a missing file. The deposit sat three sessions
    stale, and the files absent from it were exactly the ones added to make that
    statement true.

    `--compare-live` closes that gap, but it needs network, so the fetch cannot
    be a pytest. This is the part that decides the verdict, and it is pure.
    """
    snap = _load_script("22_build_public_snapshot.py")

    same = {"README.md": "aaa", "pipeline/reproduce.sh": "bbb"}
    clean = snap.compare_trees(same, dict(same))
    assert clean == {"missing_from_live": [], "only_live": [],
                     "content_differs": []}, (
        "identical trees must compare clean, or every real difference is noise")

    # The real 2026-09-07 shape: four files staged, none of them published.
    staged = dict(same, **{"README.md": "aaa",
                           "pipeline/results/gene_set_sizes.csv": "ccc"})
    live = {"pipeline/reproduce.sh": "bbb"}
    diff = snap.compare_trees(staged, live)
    assert diff["missing_from_live"] == ["README.md",
                                         "pipeline/results/gene_set_sizes.csv"]
    assert diff["only_live"] == []
    assert diff["content_differs"] == []

    # A file published that a rebuild would no longer stage is the other
    # direction, and it is the dangerous one: it means the deposit carries
    # something this repository has decided not to publish.
    diff = snap.compare_trees({"a": "1"}, {"a": "1", "leaked.txt": "2"})
    assert diff["only_live"] == ["leaked.txt"]
    assert diff["missing_from_live"] == []

    # Content drift is REPORTED and must never be confused with a set
    # difference: it moves on every prose edit, so gating on it would make this
    # red almost always -- the friction that made the byte total ungated too.
    diff = snap.compare_trees({"a": "1"}, {"a": "2"})
    assert diff["content_differs"] == ["a"]
    assert diff["missing_from_live"] == [] and diff["only_live"] == []


def test_compare_live_gates_on_the_file_set_and_not_on_content():
    """The verdict itself must fire, in both directions.

    A guard nobody has watched fail is indistinguishable from one that works, so
    this drives `report_comparison` rather than re-asserting `compare_trees`.
    """
    snap = _load_script("22_build_public_snapshot.py")

    ok = snap.report_comparison(
        {"missing_from_live": [], "only_live": [], "content_differs": []},
        "deadbee", 173, 173)
    assert ok == 0, "an identical file set must pass"

    drifted = snap.report_comparison(
        {"missing_from_live": [], "only_live": [],
         "content_differs": ["09-PAPER-DRAFT.md"] * 14},
        "deadbee", 173, 173)
    assert drifted == 0, (
        "content drift must NOT gate -- 14 files differed on 2026-09-07 purely "
        "because the deposit was three sessions behind, and a check that goes "
        "red on ordinary prose edits gets ignored")

    stale = snap.report_comparison(
        {"missing_from_live": ["README.md"], "only_live": [],
         "content_differs": []},
        "8b7cce1", 173, 169)
    assert stale == 1, "a file staged but not published must FAIL"


def test_the_number_checker_names_curation_rather_than_tampering_in_a_deposit():
    """`FROZEN ARTEFACT MISSING` reads as corruption. In a deposit it is not.

    Same liability `21_provenance_manifest.py` carried until 2026-09-07: a check
    that prints a claim it cannot support teaches the next reader to ignore it,
    and then a real failure looks the same as this one. Measured inside a
    freshly staged snapshot, `19_check_numbers.py` refused with the bare message
    over two `.npz` probes the snapshot excludes ON PURPOSE.

    The refusal itself is correct and stays -- the authorities really are
    unreadable there. Only the diagnosis changes.
    """
    checker = _load_script("19_check_numbers.py")
    missing = ["pipeline/results/bisect_macos.npz"]

    text = checker._missing_artefact_message(missing)
    assert "bisect_macos.npz" in text
    # This IS the working repository, so the builder is present and the deposit
    # wording must not appear.
    assert checker.BUILDER.is_file(), (
        "the snapshot builder is missing from the working tree, which would "
        "make every deposit-detection in this repository wrong")
    assert "NOTHING IS WRONG" not in text, (
        "the deposit wording leaked into the working repository, where a "
        "missing frozen artefact really does mean something is wrong")


def test_the_deposit_wording_fires_when_the_builder_is_absent(tmp_path):
    """Both directions, with the builder's presence as the ONLY variable.

    Same tree, same missing artefact, opposite diagnosis. Proven by moving one
    thing, because two runs that differ in several ways prove nothing about
    which one mattered -- the lesson that cost this project six sessions.
    """
    checker = _load_script("19_check_numbers.py")
    missing = ["pipeline/results/bisect_hpc4.npz"]
    real_builder = checker.BUILDER
    assert real_builder.is_file()

    try:
        checker.BUILDER = tmp_path / "22_build_public_snapshot.py"
        assert not checker.BUILDER.exists()
        deposit_text = checker._missing_artefact_message(missing)
        assert "NOTHING IS WRONG WITH THIS DEPOSIT" in deposit_text, (
            "a deposit reader is still told an intact deposit is corrupt")
        assert "bisect_hpc4.npz" in deposit_text

        checker.BUILDER.write_text("# the builder is present again\n")
        working_text = checker._missing_artefact_message(missing)
    finally:
        checker.BUILDER = real_builder

    assert "NOTHING IS WRONG" not in working_text, (
        "restoring the builder must restore the strict diagnosis; otherwise "
        "the tolerance leaks into the working repository")
    assert working_text != deposit_text


def test_the_coverage_floor_can_fail_and_gates_on_counts_not_percentages():
    """Coverage must be able to REGRESS loudly, and on the right quantity.

    The obvious gate is a floor on each document's checked FRACTION. It is the
    wrong one, and this is the measurement that says so: across one session of
    ordinary documentation work on 2026-09-07, two of thirteen documents lost a
    percentage point while every per-document checked COUNT held exactly.
    Writing prose about numbers moves the denominator and removes no guard, so a
    fraction floor would have fired twice for no defect -- and a gate that goes
    red on ordinary work gets ignored.
    """
    checker = _load_script("19_check_numbers.py")
    floors = checker.COVERAGE_FLOOR
    assert floors, "the floor table is empty, so this gate guards nothing"

    # cov maps rel -> (seen, checked, unchecked, in_style); only [1] is gated.
    at_floor = {rel: (10_000, n, [], 0) for rel, n in floors.items()}
    assert checker._coverage_floor_verdict(at_floor) == 0

    # A denominator explosion -- exactly what adding a table of exit codes does
    # -- must NOT fire, because nothing stopped being checked.
    diluted = {rel: (1_000_000, n, [], 0) for rel, n in floors.items()}
    assert checker._coverage_floor_verdict(diluted) == 0, (
        "the floor is reading a fraction; ordinary prose would make it red")

    # One fewer literal checked in one document must fail.
    victim = sorted(floors)[0]
    regressed = dict(at_floor)
    regressed[victim] = (10_000, floors[victim] - 1, [], 0)
    assert checker._coverage_floor_verdict(regressed) == 1, (
        f"losing a checked literal in {victim} did not fail the gate")

    # A document vanishing from the scan entirely is the loudest version of the
    # same failure and must not read as "nothing to check, therefore fine".
    dropped = {k: v for k, v in at_floor.items() if k != victim}
    assert checker._coverage_floor_verdict(dropped) == 1, (
        f"{victim} disappearing from the scan passed silently")
