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


def _d1_pinned_groupkfold():
    """scikit-learn 1.6.1's non-shuffle GroupKFold with `kind="stable"` -- ledger D1's
    in-process patch (results/session49_diagnostics/d1_pinned_groupkfold/), inlined
    so the deposit can run it. The reference `stable_group_kfold` must reproduce."""
    from sklearn.model_selection import GroupKFold
    from sklearn.utils import check_array

    class _Pinned(GroupKFold):
        def _iter_test_indices(self, X, y, groups):
            groups = check_array(groups, input_name="groups", ensure_2d=False, dtype=None)
            unique_groups, group_idx = np.unique(groups, return_inverse=True)
            n_samples_per_group = np.bincount(group_idx)
            indices = np.argsort(n_samples_per_group, kind="stable")[::-1]
            n_samples_per_group = n_samples_per_group[indices]
            n_samples_per_fold = np.zeros(self.n_splits)
            group_to_fold = np.zeros(len(unique_groups))
            for group_index, weight in enumerate(n_samples_per_group):
                lightest_fold = np.argmin(n_samples_per_fold)
                n_samples_per_fold[lightest_fold] += weight
                group_to_fold[indices[group_index]] = lightest_fold
            indices = group_to_fold[group_idx]
            for f in range(self.n_splits):
                yield np.where(indices == f)[0]

    return _Pinned


def test_stable_group_kfold_is_the_d1_patch_exactly():
    """A12, pinned 2026-09-24: the shipped helper IS what D1 measured, split for split.

    The pinned NSCLC and pan-cancer runs are reported as the D1 patch's numbers
    made permanent, so any difference here -- train or test indices, fold order,
    number of folds -- would make them a third, unmeasured variant. Inputs cover
    the case that matters (one row per group, every key tied), mixed group sizes
    with ties, and string labels.
    """
    pinned = _d1_pinned_groupkfold()
    rng = np.random.default_rng(1)
    compared = 0
    for t in range(150):
        n, k = int(rng.integers(10, 400)), int(rng.integers(2, 7))
        if t % 3 == 0:
            g = rng.permutation(n)
        elif t % 3 == 1:
            g = rng.integers(0, max(k, n // 4), size=n)
        else:
            g = np.array([f"P{int(x):04d}" for x in rng.integers(0, max(k, n // 3), size=n)],
                         dtype=object)
        if len(np.unique(g)) < k:
            continue
        ref = list(pinned(n_splits=k).split(np.zeros((n, 1)), groups=g))
        got = list(splits.stable_group_kfold(g, k))
        assert len(got) == len(ref) == k
        for (a_tr, a_te), (b_tr, b_te) in zip(ref, got):
            assert np.array_equal(a_tr, b_tr) and np.array_equal(a_te, b_te), f"case {t}"
        compared += 1
    assert compared >= 120, f"only {compared} cases compared"


def test_stable_group_kfold_known_answer_and_what_it_leaves_alone():
    """The pinned order is a fixed answer, and it moves nothing but ties.

    `caf92a3f83ddb5b9` is the fold hash of a 944-patient permutation under the
    stable order that macOS arm64 and HPC4 were measured to AGREE on
    (`26_sort_audit.py --library-order`, A12). A known answer, not a
    self-comparison: it would catch a change of tie rule on any platform.
    With every group a different size there are no ties, so the pinned helper
    must equal scikit-learn's own GroupKFold -- the pin changes tie order only.
    """
    import hashlib
    import inspect
    import re

    from sklearn.model_selection import GroupKFold

    g = np.random.default_rng(0).permutation(944)
    fold = np.full(944, -1)
    for f, (_, te) in enumerate(splits.stable_group_kfold(g, 5)):
        fold[te] = f
    assert hashlib.sha256(fold.astype(np.int64).tobytes()).hexdigest()[:16] == "caf92a3f83ddb5b9"
    assert sorted(np.bincount(fold)) == [188, 189, 189, 189, 189]

    sizes = np.repeat(np.arange(20), np.arange(1, 21))      # group g has g+1 rows
    lib = list(GroupKFold(n_splits=5).split(np.zeros((len(sizes), 1)), groups=sizes))
    got = list(splits.stable_group_kfold(sizes, 5))
    assert all(np.array_equal(a[1], b[1]) and np.array_equal(a[0], b[0])
               for a, b in zip(lib, got))
    # ...and the library's own sort is NOT what ships: on all-tied groups the two
    # disagree on this platform or the other (A12), so assert that none of the
    # three call sites reaches it. (`StratifiedGroupKFold`, the unused stratified
    # branch, shuffles with a seed and is a different class.)
    library = re.compile(r"(?<!Stratified)GroupKFold\(")
    for fn in (splits.random_patient_split, models.site_prediction_control,
               models.site_prediction_control_oof_combat):
        src = inspect.getsource(fn)
        assert not library.search(src), fn.__name__
        assert "stable_group_kfold(" in src, fn.__name__
    with pytest.raises(ValueError, match="number of groups"):
        list(splits.stable_group_kfold(np.array([0, 0, 1, 1]), 3))


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


def _shared_mask_case():
    """Synthetic stand-in for pan-cancer's 19 NaN-axis patients: a block of rows
    NaN in EVERY target column, plus one column with its own extra NaN so the
    grouping has two masks to keep apart."""
    rng = np.random.default_rng(48)
    n, d, m = 300, 40, 25
    patients = pd.Series([f"TCGA-02-{i:04d}" for i in range(n)])
    X = rng.normal(size=(n, d))
    Y = X[:, :1] * rng.normal(size=(1, m)) + rng.normal(size=(n, m))
    Y[rng.choice(n, 7, replace=False), :] = np.nan
    Y[rng.choice(n, 3, replace=False), 4] = np.nan
    split = splits.random_patient_split(patients, n_folds=5, seed=0)
    alphas = np.array([0.5, 3.0, 37.92690190732246, 1e3, 1e4])
    return X, Y, split, patients.to_numpy(), alphas


def test_shared_mask_fast_path_is_bitwise_equal_to_per_column_ridge():
    """The fast path is a speed change ONLY, and "only" means bit for bit.

    `_ridge_shared_mask_predict` computes X'X once per fold and mask instead of
    once per column. If it moved a single ULP, every pan-cancer number produced
    through it would stop reproducing the frozen runs, and the change would be
    an estimator-path change requiring its own decision. So the comparison is
    `np.array_equal`, not a tolerance, against the sklearn loop it replaces.
    """
    X, Y, split, patients, alphas = _shared_mask_case()
    for k in range(split.n_folds):
        tr, _ = split.indices(k)
        assert not np.isfinite(Y[tr]).all(axis=0).any(), (
            "the case must drive every column into the fallback, or it tests nothing")
    try:
        models.SHARED_MASK_FAST_PATH = False
        ref = models.cross_val_predict_multi(X, Y, split, patients, fixed_alpha=alphas)
        models.SHARED_MASK_FAST_PATH = True
        fast = models.cross_val_predict_multi(X, Y, split, patients, fixed_alpha=alphas)
    finally:
        models.SHARED_MASK_FAST_PATH = True
    assert np.isfinite(ref).sum() > 0.9 * ref.size
    assert np.array_equal(np.isnan(ref), np.isnan(fast))
    assert np.array_equal(ref[np.isfinite(ref)], fast[np.isfinite(fast)]), (
        f"max |diff| {np.nanmax(np.abs(ref - fast)):.3e}: the fast path is no longer "
        "bitwise equal to per-column sklearn Ridge")


def test_the_shared_mask_equality_check_can_fail(monkeypatch):
    """A one-ULP perturbation inside the fast path must be visible to the check."""
    X, Y, split, patients, alphas = _shared_mask_case()
    original = models._ridge_shared_mask_predict

    def nudged(*args):
        out = original(*args)
        return np.nextafter(out, np.inf)

    models.SHARED_MASK_FAST_PATH = False
    ref = models.cross_val_predict_multi(X, Y, split, patients, fixed_alpha=alphas)
    models.SHARED_MASK_FAST_PATH = True
    monkeypatch.setattr(models, "_ridge_shared_mask_predict", nudged)
    fast = models.cross_val_predict_multi(X, Y, split, patients, fixed_alpha=alphas)
    assert not np.array_equal(ref[np.isfinite(ref)], fast[np.isfinite(fast)])


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
import re  # noqa: E402
import subprocess  # noqa: E402

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


# ---------------------------------------------------------------------------
# TRIAGE: the twelve scripts/ with no mention in this file (measured 2026-09-08)
#
# All TWELVE src/aacr27/ modules are exercised by name here; twelve of the
# twenty-six scripts/ are not. That asymmetry is not automatically a defect --
# the library is where the estimand lives -- but it had never been triaged, so a
# decision is recorded per script rather than "add tests" or "leave it".
#
#   (a) NEEDS NETWORK OR data/  -- cannot run in CI or from a deposit, and a
#       mock would test the mock. Their real failure mode is a changed upstream
#       URL, which no unit test can see.
#         01_fetch_data.py, 03_fetch_signatures.py, 05_run_nsclc.py,
#         09_run_pancancer.py, 02_run_audit.py (--real path only; --demo is
#         synthetic and IS covered through the library)
#       DECISION: no pytest. The `--check` path of reproduce.sh exercises them
#       where inputs exist, and 20_check_versions.py guards the environment.
#
#   (b) MULTI-HOUR -- 05 and 09 are the cohort runs (227 s to 4.75 h); 15/16/17
#       are alpha probes, one of which ran >5 h against a docstring promising
#       "seconds".
#         05_run_nsclc.py, 09_run_pancancer.py, 15_alpha_probe.py,
#         16_alpha_margin.py, 17_alpha_margin_fast.py
#       DECISION: no pytest. Their outputs are frozen artefacts and PROVENANCE
#       hashes all 114, which is a stronger guard than re-running them would be.
#
#   (c) CLOSED NEGATIVE, deliberately frozen -- the CPTAC arm was dropped after
#       Gate 1 failed on four measured routes and is reported as Limitations 10.
#         04_budget_cptac.py, 07_triage_cptac.py
#       DECISION: no pytest, ever. Testing them would assert that a closed
#       negative still computes, which is not a property anyone needs.
#
#   (d) PURE AND CHEAP -- the only bucket that warranted new tests.
#         08_count_abstract.py  -- CLOSED 2026-09-08, four tests, all mutation-
#                                  proven, including the exit-2 path session 39
#                                  found dying with a raw traceback.
#         14_repro_probe.py     -- CLOSED 2026-09-16 (session 49): run on a
#                                  synthetic interim tree; deterministic, and
#                                  a changed input is flagged first at stage 1.
#
#   ADDED LATER, outside the 2026-09-08 count: 12_partition_variance.py's shard
#   and failure-kind behaviour, 27_partition_sweep.py and 00_build_interim.py
#   (new 2026-09-16) are tested at the end of this file (session 49), with
#   mutation checks. 28_type_only_concordance.py (new 2026-09-17) needs
#   data/interim and belongs to bucket (a); its one output is an authority in
#   19_check_numbers.py, which is what guards it.
#
#   NOT UNTESTED, despite appearing in the list: 26_sort_audit.py is tested
#   OUTSIDE pytest by `--self-test` in BOTH directions and is in the
#   checks-that-can-fail list; so is poster/make_qr.py --check. A test that lives
#   outside pytest still counts. Do not add a redundant pytest for either.
# ---------------------------------------------------------------------------


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
        # 2400 s, raised from 900 on 2026-09-17 (session 50): the self-test grew
        # to ~276 injected errors and took 698 s at 252 on a quiet macOS machine,
        # so 900 s left no room for load -- the closing suite of that session
        # timed out at 900 s under a foreign job at load1 27-40 with every other
        # test green (results/full_suite_20260917_run50.log).
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=2400)
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

    scanned = [p for p in results.iterdir() if not p.name.startswith(".")]
    # A globbing test that finds nothing reports "all clear" while proving
    # nothing -- how 19_check_numbers.py v1 announced "62 claims OK" over an
    # empty scan. `readme.exists()` above only catches results/ vanishing
    # WHOLESALE; it says nothing about an exclusion filter quietly broadening
    # until everything is skipped. results/ held 101 entries on 2026-09-08.
    assert len(scanned) >= 60, (
        f"only {len(scanned)} entries were scanned in results/; this check has "
        "gone vacuous and would pass over an unclassified run directory")

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


# ------------------------------- a row must not say `pending` for ever ---
#
# Since session 46 a run's row is written BEFORE the run, with the literal
# word `pending` where the result goes, and filled from the log afterwards.
# Nothing checked the second half. Session 54 found two rows that still read
# `pending` for runs that had long finished: `full_suite_20260920_run53d.log`,
# a RED run whose two failures were the reason for the next run, and
# `wrapper_20260916_run47e.log`. The classification check above passes over
# both, because it asks only whether the NAME appears.
#
# Decidable only where the log itself says it finished: a pytest terminal
# summary, or the timing wrapper's SUITE_EXIT / TOTAL_ELAPSED_SECONDS lines.
# A row whose log is absent (a run in a clean clone writes its log into the
# working tree, not the clone) or still being written is skipped, which is
# exactly the state a row is in while its own run is in flight.

_PENDING_MARKER = re.compile(r"RESULT: pending|literal word `pending`|with `pending`")
_RESULT_FILLED = re.compile(r"RESULT\b(?!: pending)")
_LOG_FINISHED = re.compile(
    r"^(=+ )?\d+ (passed|failed)\b.* in [\d.]+s|^SUITE_EXIT=|^TOTAL_ELAPSED_SECONDS=",
    re.M)


def _rows_left_pending(text, results):
    """(rows decided, names of finished logs whose row was never filled)."""
    decided, stale = 0, []
    for line in text.splitlines():
        if not line.startswith("| `") or not _PENDING_MARKER.search(line):
            continue
        names = re.findall(r"`([^`]+\.log)`", line.split("|")[1])
        logs = [results / n for n in names]
        if not names or not all(
                p.exists() and _LOG_FINISHED.search(p.read_text(errors="replace"))
                for p in logs):
            continue
        decided += 1
        if not _RESULT_FILLED.search(line):
            stale.append(", ".join(names))
    return decided, stale


def test_no_results_row_still_reads_pending_after_its_log_finished():
    """A row written before its run must be filled once the log has finished."""
    results = SCRIPTS.parent / "results"
    decided, stale = _rows_left_pending(
        (results / "README.md").read_text(errors="replace"), results)
    # 30 rows were decidable on 2026-09-21. Far fewer means the marker or the
    # finished-log pattern stopped matching, and an empty scan reports clean.
    assert decided >= 25, (
        f"only {decided} pending-convention rows could be decided; the check "
        "has gone vacuous")
    assert not stale, (
        "these logs have finished but their results/README.md row still reads "
        "`pending`: " + "; ".join(stale) + "\nFill the RESULT from the log; "
        "keep the original row text, since it records the expectation.")


def test_the_pending_row_check_can_fail(tmp_path):
    """Flags a finished-but-unfilled row, and nothing else."""
    (tmp_path / "done.log").write_text("....\n3 passed in 1.23s\n")
    (tmp_path / "running.log").write_text("....\n")
    row = "| `{}` | **ROW WRITTEN BEFORE THE RUN, with the literal word `pending`.** {} |"
    text = "\n".join([
        row.format("done.log", "**RESULT: pending.**"),
        row.format("running.log", "**RESULT: pending.**"),
        row.format("absent.log", "**RESULT: pending.**"),
    ])
    assert _rows_left_pending(text, tmp_path) == (1, ["done.log"]), (
        "the check did not flag exactly the one finished, unfilled row")
    filled = text.replace(
        "`done.log` | **ROW WRITTEN BEFORE THE RUN, with the literal word "
        "`pending`.** **RESULT: pending.**",
        "`done.log` | **ROW WRITTEN BEFORE THE RUN, with the literal word "
        "`pending`.** **RESULT: pending.** **RESULT: GREEN, 3 passed.**")
    assert filled != text
    assert _rows_left_pending(filled, tmp_path) == (1, []), (
        "a filled row was still flagged")


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


# `ENVIRONMENT.md` writes IMPORT names; `requirements-lock.txt` pins
# DISTRIBUTION names. The two differ for exactly one package here, and that is
# precisely why the lockfile is an independent authority rather than a second
# copy of the same thing.
_IMPORT_TO_DISTRIBUTION = {"sklearn": "scikit-learn"}

# Packages `ENVIRONMENT.md` reports that the lockfile deliberately does not pin,
# each with the reason. Kept as a named set rather than an omission so the
# reconciliation below can prove the two lists cover the table.
_ENV_PACKAGES_NOT_PINNED: dict[str, str] = {}


def _environment_md_package_versions(text: str) -> dict[str, str]:
    """Parse the `## Packages` table into {package: version}."""
    import re as _re
    section = text.split("## Packages", 1)
    if len(section) < 2:
        return {}
    body = section[1].split("\n## ", 1)[0]
    rows = _re.findall(r"^\|\s*([A-Za-z][\w.-]*)\s*\|\s*([0-9][^|]*?)\s*\|", body, _re.M)
    return {name: version for name, version in rows}


def _lockfile_pins(text: str) -> dict[str, str]:
    import re as _re
    return {m.group(1).lower(): m.group(2)
            for m in _re.finditer(r"^([A-Za-z][\w.-]*)==([^\s;#]+)", text, _re.M)}


def test_environment_md_versions_agree_with_the_lockfile():
    """`ENVIRONMENT.md`'s versions, checked against an authority it never reads.

    `test_environment_md_matches_a_fresh_render_on_this_platform` compares the
    committed file against `23_export_environment.py`'s own `render()`. Both
    sides come from that one script, so a renderer that reported the wrong
    version -- reading the wrong distribution, truncating, or caching a stale
    import -- would produce a committed file and a fresh render that agreed
    perfectly. It also SKIPS entirely off-platform, so on a collaborator's
    machine and on HPC4 it asserts nothing at all.

    `requirements-lock.txt` is written by a different process for a different
    purpose and the exporter never reads it, which makes it the independent
    side. It is also platform-independent, so unlike its sibling this check
    runs everywhere rather than skipping where it is most needed.
    """
    committed = SCRIPTS.parent / "ENVIRONMENT.md"
    lock = SCRIPTS.parent / "requirements-lock.txt"
    assert committed.exists() and lock.exists()

    reported = _environment_md_package_versions(
        committed.read_text(errors="replace"))
    pins = _lockfile_pins(lock.read_text(errors="replace"))

    assert len(reported) >= 5, (
        f"parsed only {len(reported)} package row(s) out of ENVIRONMENT.md's "
        "`## Packages` table; the parser and the table have drifted apart and "
        "this check is near-vacuous")
    assert len(pins) >= 5, f"parsed only {len(pins)} pin(s) from the lockfile"

    # UNCLASSIFIED-MEMBER, both directions: every package the table reports is
    # either pinned or documented as deliberately unpinned, and the documented
    # set may not name a package that has left the table.
    unclassified = sorted(
        name for name in reported
        if _IMPORT_TO_DISTRIBUTION.get(name, name).lower() not in pins
        and name not in _ENV_PACKAGES_NOT_PINNED)
    assert not unclassified, (
        f"ENVIRONMENT.md reports these packages and the lockfile pins none of "
        f"them: {unclassified}. Add each to the lockfile or to "
        "`_ENV_PACKAGES_NOT_PINNED` with a reason -- an undocumented omission "
        "is how a set goes quietly stale.")
    dead = sorted(set(_ENV_PACKAGES_NOT_PINNED) - set(reported))
    assert not dead, f"these entries name packages ENVIRONMENT.md no longer reports: {dead}"

    compared, wrong = 0, []
    for name, version in sorted(reported.items()):
        dist = _IMPORT_TO_DISTRIBUTION.get(name, name).lower()
        if dist not in pins:
            continue
        compared += 1
        if pins[dist] != version:
            wrong.append(f"{name}: ENVIRONMENT.md says {version}, "
                         f"requirements-lock.txt pins {dist}=={pins[dist]}")
    assert compared >= 5, (
        f"only {compared} package(s) were actually compared; the import-to-"
        "distribution mapping has drifted and this check is near-vacuous")
    assert not wrong, (
        "ENVIRONMENT.md disagrees with the lockfile:\n  " + "\n  ".join(wrong))


def test_the_environment_lockfile_agreement_check_can_fail():
    """Prove it discriminates, in the direction a renderer defect would take.

    Mutate the PARSED table in memory -- never the repository -- the way a
    renderer reporting a stale version would, and require the comparison to
    name that package.
    """
    committed = (SCRIPTS.parent / "ENVIRONMENT.md").read_text(errors="replace")
    lock = (SCRIPTS.parent / "requirements-lock.txt").read_text(errors="replace")
    reported = _environment_md_package_versions(committed)
    pins = _lockfile_pins(lock)

    victims = [n for n in reported
               if _IMPORT_TO_DISTRIBUTION.get(n, n).lower() in pins]
    assert victims, "expected at least one package present in both"
    name = sorted(victims)[0]
    mutated = dict(reported)
    mutated[name] = "0.0.1"

    wrong = [n for n, v in mutated.items()
             if _IMPORT_TO_DISTRIBUTION.get(n, n).lower() in pins
             and pins[_IMPORT_TO_DISTRIBUTION.get(n, n).lower()] != v]
    assert wrong == [name], (
        f"reporting {name} as 0.0.1 while the lockfile pins "
        f"{pins[_IMPORT_TO_DISTRIBUTION.get(name, name).lower()]} did not "
        "register, so this check cannot detect the defect it exists for")
    # Positive control: unmutated, the same comparison must be silent.
    assert not [n for n, v in reported.items()
                if _IMPORT_TO_DISTRIBUTION.get(n, n).lower() in pins
                and pins[_IMPORT_TO_DISTRIBUTION.get(n, n).lower()] != v]


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
#     macOS 26.5.2 / arm64, Python 3.13.9, openblas 0.3.21   1.036173699004058
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


def test_cohort_fingerprint_matches_its_own_stated_definition():
    """The independent companion to the reader/writer sync test above (F6.4).

    That test runs the reader and the writer on one input and asks only whether
    they agree, so a defect both share -- the same wrong separator, the same
    missing de-duplication -- passes it. This one holds the writer to the
    definition it writes into every `summary.json` beside the digest, and to a
    KNOWN ANSWER computed once from that sentence (session 50, 2026-09-17), not
    from either implementation at test time.
    """
    from aacr27.experiment import AuditResult

    patients = ["TCGA-B", "TCGA-A", "TCGA-B", "TCGA-C", 7, "10", "9"]
    written = AuditResult.cohort_fingerprint(
        types.SimpleNamespace(predictions=pd.DataFrame({"patient": patients})))
    assert written["patient_set_sha256_note"] == (
        "sha256 of the sorted unique patient ids, newline-joined"), (
        "the stated definition changed; recompute the known answer from the new one")
    # sha256(b"10\n7\n9\nTCGA-A\nTCGA-B\nTCGA-C"): the six unique ids as
    # strings, sorted as strings, joined by newlines.
    assert written["patient_set_sha256"] == (
        "51aca02c361083470492421b4abf5d720fc62a14979a07b7d63757dd9ccb44e3")


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
#   decomposition.permutation_calibration  <- 11_close_science_gaps.py, main()
#   ancestry.cramers_v_calibrated          <- 06_run_ancestry.py, main(), twice
#   ancestry.pooled_over_signatures        <- 06_run_ancestry.py, main(), 3 times
#
# [Named by function since 2026-09-22. The first line read
# `11_close_science_gaps.py:403` until then, which by that date was the
# `label_site_variance` call; `permutation_calibration` is called 24 lines
# further down.]
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


def test_audit_config_seed_is_the_seed_the_frozen_runs_recorded():
    """The independent companion to the seed sync test above (F6.4, session 50).

    That test asserts every unwired default equals `AuditConfig.seed`, so a change
    to `AuditConfig.seed` itself moves both sides together and passes it while
    every rerun stops reproducing the frozen results. The frozen runs recorded
    the seed they used in `config.json`; that record, not the code, is the
    external anchor.
    """
    import json as _json

    results = Path(__file__).resolve().parents[1] / "results"
    for run in ("nsclc_v3", "pancancer_v3"):
        recorded = _json.loads((results / run / "config.json").read_text())["seed"]
        assert experiment.AuditConfig().seed == recorded, (
            f"AuditConfig.seed is {experiment.AuditConfig().seed} but {run} was "
            f"frozen with seed {recorded}; a rerun would not reproduce it")


def test_audit_config_n_boot_is_the_count_the_frozen_runs_recorded():
    """E2, decided 2026-09-24 (session 59): `AuditConfig.n_boot` defaulted to
    2,000 while every frozen run recorded 1,000, so a bare `AuditConfig()`
    produced an interval the protocol never used. Anchored to the runs' own
    `config.json`, as the seed test above is, so a change to the default
    cannot move both sides together."""
    import json as _json

    results = Path(__file__).resolve().parents[1] / "results"
    for run in ("nsclc_v3", "pancancer_v3"):
        recorded = _json.loads((results / run / "config.json").read_text())["n_boot"]
        assert experiment.AuditConfig().n_boot == recorded, (
            f"AuditConfig.n_boot is {experiment.AuditConfig().n_boot} but {run} "
            f"was frozen with n_boot {recorded}")


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
    # Session 59 (B-21): Limitation 8's note, lifted out of the manuscript
    # 2026-09-24 and submitted as Supplementary Note 1.
    "SUPPLEMENTARY-NOTE-1.md",
    # B-21, 2026-09-24: the Results passages lifted to bring the body under
    # 5,000 words, submitted as Supplementary Note 2.
    "SUPPLEMENTARY-NOTE-2.md",
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
    bad, seen = {}, 0
    for rel in _US_SPELLING_DOCS:
        path = repo / rel
        if not path.exists():
            # A staged snapshot renames and drops documents by design; on the
            # working tree a missing submitted document is itself the defect.
            assert not _ON_WORKING_TREE, f"submitted document is missing: {rel}"
            continue
        seen += 1
        hits = british_spellings_in_prose(
            path.read_text(encoding="utf-8"), is_markdown=rel.endswith(".md"))
        if hits:
            bad[rel] = hits
    assert seen >= 1, "no document was actually read; the scan is vacuous"
    assert not bad, (
        "British spelling survives in prose that gets submitted (US spelling "
        "was decided and converted in commit 3ee560a):\n"
        + "\n".join(f"  {rel}:{ln}  {w}"
                    for rel, hits in bad.items() for ln, w in hits)
        + "\nIf a hit is inside a reference title or a code identifier it "
          "should be in backticks or below `## References`, which this scan "
          "already skips -- fix the markup rather than the word list.")


# ---------------------------------------------------------------------------
# The NINTH axis, named in session 45: UNCLASSIFIED MEMBER.
#
# Session 44 found SILENT-SCOPE -- a guard whose word list was complete and
# whose FILE SET was not. This is the sharper variant. `_SUBMITTED_DOCS` names
# four documents, and the comment above it names five more as deliberate
# exclusions. Nine documents accounted for, inclusions AND exclusions written
# down, which is exactly what a complete classification looks like. The tree
# holds THIRTY-FOUR tracked prose documents. Twenty-five belonged to neither
# set, and because the guard documented its own exclusions nobody asked whether
# the two lists together covered anything.
#
# A set that documents its exclusions looks finished. That appearance is the
# defect: it converts "nobody has classified these" into something that reads
# like "these were considered and rejected".
#
# The fix is not a longer list -- a longer list has the same hole one document
# later. It is the reconciliation below, which fails until every document in
# the tree sits in exactly one bucket. Adding a document now forces a decision.
# ---------------------------------------------------------------------------

# Deposited prose a reader meets as the project's own finished writing, held to
# the same US spelling as the submitted four. Both were MEASURED clean before
# being added (zero British hits each, confirmed by an independent raw grep, so
# the zero is not the scan failing to look), which makes this a free ratchet
# rather than a conversion.
_DEPOSITED_PROSE_DOCS = (
    "snapshot/README.md",                          # the public repo's root README
    "pipeline/results/supplementary/CAPTIONS.md",  # becomes the submitted captions
    # Session 59 (B-11, decided 2026-09-24): converted to US spelling and
    # enforced from then on -- CITATION.cff's title had said "tumour" while the
    # manuscript's says "tumor". `axis_residualised` stays, as a data value in
    # backticks. DATA-LICENSES.md is the data note split out of LICENSE (B-7).
    "CITATION.cff",
    "pipeline/results/supplementary/MANIFEST.md",
    "DATA-LICENSES.md",
)

_US_SPELLING_DOCS = _SUBMITTED_DOCS + _DEPOSITED_PROSE_DOCS

# Everything else, each with the reason it keeps its spelling. "Working
# document" is a real reason here: these sit beside British code identifiers
# like `residualise_matrix`, and mixing conventions reads worse than leaving
# them. A reason of "" would be a placeholder, and the test rejects it.
_SPELLING_NOT_ENFORCED = {
    "00-PROJECT-BRIEF.md": "historical brief; not deposited, records what was believed then",
    "01-LANDSCAPE-AUDIT.md": "working audit",
    "02-PROJECT-DECISION.md": "historical record; its wording states what was true when written",
    "03-WEEK1-VERIFICATION.md": "working record",
    "04-CODE-AUDIT.md": "working audit, quotes code identifiers throughout",
    "05-PRE-REGISTRATION.md": "protocol frozen at tag prereg-2026-08-18; retro-editing it would falsify the record",
    "06-ANCESTRY-ARM.md": "working analysis note",
    "08-READINESS-ASSESSMENT.md": "working record",
    "10-PROJECT-REFRAME.md": "historical record of the pivot",
    "11-SUBMISSION-PACK.md": "working logistics",
    "12-MENTOR-OUTREACH.md": "working logistics",
    "13-NEXTGEN-EXTENDED-ABSTRACT.md": "superseded draft for a path not taken",
    "14-SCIENCE-AUDIT.md": "deposited, but an audit record quoting British code identifiers",
    "15-CHECKLIST.md": "working checklist",
    "16-YOUR-TASKS.md": "private logistics, never deposited",
    "17-CPTAC-PLAN.md": "closed-negative working record",
    "18-SUPPLEMENTARY-INVENTORY.md": "working inventory",
    "bagaev.html": "a saved third-party web page, not project prose",
    "pipeline/ENVIRONMENT.md": "generated environment record",
    "pipeline/README.md": "sits beside British code identifiers, deliberately",
    "pipeline/cluster/a9-investigation-scratch/README.md": "scratch notes from the A9 investigation",
    "pipeline/results/README.md": "working classification of every run and log",
    "poster/README.md": "build notes for the poster",
    "submission/ABSTRACT-PORTAL-DRY-RUN.md": "operational script for the author, not submitted prose",
    "submission/PREPRINT-PLAN.md": "working plan",
    "submission/SUBMISSION-CHECKLIST.md": "working checklist",
}


def _tracked_prose(repo: Path) -> set[str]:
    """Every TRACKED prose document, DISCOVERED rather than listed.

    Discovery is the whole point: a hand-maintained list cannot notice a
    document nobody added to it. Handoffs are excluded because they are
    archives.

    TRACKED, and the word is load-bearing. The first version of this walked the
    filesystem, and within the hour it failed on `poster/poster_preview_tmp.html`
    -- a 228 KB stray that appeared in the working tree from outside this
    session and belongs to nobody's classification. Several Claude sessions and
    an external review run against this repository, and a check that ranges over
    whatever happens to be on disk reports their scratch files as project
    defects. `git ls-files` is the right universe: a document becomes the
    project's when it is staged, which is also the moment someone can be asked
    to classify it. Safe here because every caller is gated on
    `_ON_WORKING_TREE`, so a staged snapshot never reaches this.
    """
    out = subprocess.run(
        ["git", "ls-files", "-z", "*.md", "*.html", "*.cff"],
        cwd=repo, capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, (
        f"git ls-files failed in {repo}, so this check cannot know what the "
        f"project's documents are:\n{out.stderr[-500:]}")
    return {f for f in out.stdout.split("\0") if f and "HANDOFF" not in Path(f).name}


# `16-YOUR-TASKS.md` is never deposited (there is a test that enforces exactly
# that), so its presence is the cheapest proof we are on the working tree and
# not inside a staged snapshot, where most of these documents are absent by
# design and a completeness check would be meaningless.
_ON_WORKING_TREE = (Path(__file__).resolve().parents[2] / "16-YOUR-TASKS.md").exists()


@pytest.mark.skipif(not _ON_WORKING_TREE,
                    reason="classification is a working-tree property")
def test_every_prose_document_is_classified_for_spelling():
    """No prose document may sit outside BOTH spelling buckets.

    This is the guard that would have caught session 44's defect a session
    earlier, had the figures been documents: not by knowing about figures, but
    by refusing to let anything go unclassified.
    """
    repo = Path(__file__).resolve().parents[2]
    found = _tracked_prose(repo)
    assert len(found) >= 25, (
        f"only {len(found)} prose document(s) discovered; the scan is "
        "near-vacuous and would pass over an empty tree")

    classified = set(_US_SPELLING_DOCS) | set(_SPELLING_NOT_ENFORCED)

    unclassified = sorted(found - classified)
    assert not unclassified, (
        "these prose documents belong to neither spelling bucket, so nobody "
        "has decided whether a reader is meant to see US spelling in them:\n"
        + "\n".join(f"  {f}" for f in unclassified)
        + "\nAdd each to `_US_SPELLING_DOCS` (and make it pass) or to "
          "`_SPELLING_NOT_ENFORCED` with the reason it keeps its spelling.")

    dead = sorted(classified - found)
    assert not dead, (
        "these documents are classified but no longer exist, so the "
        "classification is drifting away from the tree:\n"
        + "\n".join(f"  {f}" for f in dead))

    blank = sorted(k for k, v in _SPELLING_NOT_ENFORCED.items() if not v.strip())
    assert not blank, (
        f"exclusions with no stated reason are placeholders, not decisions: {blank}")


@pytest.mark.skipif(not _ON_WORKING_TREE,
                    reason="classification is a working-tree property")
def test_every_prose_document_is_scanned():
    """The same reconciliation for the bracket and bold scans.

    `_prose_documents()` is a glob plus a hand-maintained `extra` list, and
    session 45 measured it covering 30 of 34 documents. One of the four it
    missed was `poster/poster.html` -- named in `_SUBMITTED_DOCS` in this very
    file, so two hand-maintained document sets disagreed with nothing to see it.
    """
    repo = Path(__file__).resolve().parents[2]
    found = _tracked_prose(repo)
    scanned = {p.relative_to(repo).as_posix() for p in _prose_documents(repo)}
    # VACUITY FLOOR (F10.4, session 50). An empty `found` -- a failed git call,
    # a changed suffix list -- makes `missed` empty and this test green over
    # nothing. Measured 2026-09-17: 34 tracked prose documents, 33 scanned.
    assert len(found) >= 30 and len(scanned) >= 30, (
        f"only {len(found)} tracked / {len(scanned)} scanned prose documents; "
        "the reconciliation below would pass while comparing nothing")

    missed = sorted(found - scanned - set(_PROSE_NOT_SCANNED))
    assert not missed, (
        "these prose documents are scanned by no delimiter or bold check, and "
        "are not recorded as deliberate exclusions:\n"
        + "\n".join(f"  {f}" for f in missed)
        + "\nMeasure the noise each one adds FIRST, then either add it to "
          "`_prose_documents()` or record it in `_PROSE_NOT_SCANNED`.")

    overlap = sorted(scanned & set(_PROSE_NOT_SCANNED))
    assert not overlap, (
        f"recorded as not scanned, yet scanned: {overlap}")


@pytest.mark.skipif(not _ON_WORKING_TREE,
                    reason="classification is a working-tree property")
def test_the_prose_classification_checks_can_fail():
    """Both reconciliations must notice a document that belongs to no bucket.

    Paired with a positive control, because a detector that is red on
    everything passes its own can-fail test. Mutation is applied to the REAL
    sets, in the direction the real defect takes: a document appears in the
    tree and nobody classifies it.
    """
    repo = Path(__file__).resolve().parents[2]
    found = _tracked_prose(repo)
    classified = set(_US_SPELLING_DOCS) | set(_SPELLING_NOT_ENFORCED)

    # Positive control: the real tree is fully classified right now.
    assert not (found - classified), "precondition: the real tree is classified"

    # Mutation: a new document lands and nobody classifies it.
    assert ({"19-NEW-DOCUMENT.md"} | found) - classified, (
        "the reconciliation cannot see an unclassified document, which is the "
        "only thing it exists to catch")

    # Mutation: a classified document is deleted and the entry goes stale.
    assert set(list(classified) + ["99-DELETED.md"]) - found, (
        "the reconciliation cannot see a dead classification entry")

    # And the scanned-set reconciliation, same two directions.
    scanned = {p.relative_to(repo).as_posix() for p in _prose_documents(repo)}
    assert not (found - scanned - set(_PROSE_NOT_SCANNED)), (
        "precondition: every document is scanned or recorded as excluded")
    assert ({"20-UNSCANNED.md"} | found) - scanned - set(_PROSE_NOT_SCANNED), (
        "the scanned-set reconciliation cannot see an unscanned document")


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


def test_the_public_snapshot_can_build_the_supplementary_tables(tmp_path):
    """The third consumer of the deposit, pinned instead of hand-checked.

    `10_make_figures.py` and `21_provenance_manifest.py` each earned a test only
    after they were found broken inside a staged tree. `25_build_supplementary.py`
    was verified by hand in sessions 37, 38 and 39 -- three times, always green,
    never pinned. A hand verification is one session from being forgotten, and
    the manuscript cites all nine S-tables by number, so a deposit that cannot
    render them is a deposit that cannot support the paper.

    Note the shape difference that makes this safe: 25 takes no `--outdir`, it
    writes into `results/supplementary/` of the tree it is run from. Run inside a
    temp snapshot that is exactly what we want, and the committed tables are
    untouched.
    """
    out = tmp_path / "snap"
    staged = subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert staged.returncode == 0, f"staging failed:\n{staged.stdout[-2000:]}"
    assert not (out / "pipeline" / "data").exists(), (
        "this test is only meaningful while the snapshot excludes data/; a "
        "reader reproducing from the deposit has no data/ to fall back on")

    built = subprocess.run(
        [sys.executable, str(out / "pipeline" / "scripts" / "25_build_supplementary.py")],
        cwd=out / "pipeline", capture_output=True, text=True, timeout=600)
    assert built.returncode == 0, (
        "the public deposit cannot render the supplementary tables the "
        f"manuscript cites by number:\n{built.stdout[-3000:]}\n{built.stderr[-3000:]}")

    # Exit 0 is not enough: a builder that skipped eight tables and wrote one
    # would also exit 0. Name every table, so a partial render fails here.
    supp = out / "pipeline" / "results" / "supplementary"
    for stem in ("S8_ancestry", "S3_small_panel_bracket", "S1_per_signature_excess",
                 "S2_reliability_alpha", "S7_negative_controls",
                 "S4_label_site_variance", "S9_decomposition",
                 "S6_partition_variance", "S10_platform_reproducibility"):
        table = supp / f"{stem}.tsv"
        assert table.exists(), f"{stem}.tsv was not rendered from the deposit"
        # A zero-row table is a silent failure the exit code cannot show.
        assert len(table.read_text().splitlines()) > 1, f"{stem}.tsv has no rows"
    for meta in ("MANIFEST.md", "CAPTIONS.md"):
        assert (supp / meta).exists(), f"{meta} was not rendered from the deposit"


def test_the_snapshot_supplementary_check_can_fail(tmp_path):
    """Prove the test above discriminates rather than merely being green.

    Direction: delete a source artefact the snapshot DOES stage. The builder must
    refuse loudly. Without this, a future change that made a missing source
    render an empty table would leave the test above passing while the deposit
    quietly shipped a supplement with a blank S3.
    """
    out = tmp_path / "snap"
    subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300, check=True)

    victim = out / "pipeline" / "results" / "pancancer_v3" / "small_panel_bracket.csv"
    assert victim.exists(), "expected S3's source to be staged in the deposit"
    victim.unlink()

    built = subprocess.run(
        [sys.executable, str(out / "pipeline" / "scripts" / "25_build_supplementary.py")],
        cwd=out / "pipeline", capture_output=True, text=True, timeout=600)
    assert built.returncode != 0, (
        "removing a supplementary source did not fail the build, so the check "
        "above cannot detect the defect it exists for")
    assert "MISSING SOURCE" in built.stdout + built.stderr


def test_every_library_module_imports_inside_the_deposit(tmp_path):
    """The last unswept axis of the "excluded but needed" class: the library.

    The class is *a file the snapshot excludes that something the deposit ships
    actually needs*. It has been swept twice. Session 38 swept the ARTIFACT axis
    (staged tree vs published tree); session 39 swept the CONSUMER axis by
    running every script inside a staged tree, and found three members including
    an advertised entry point that exits 2 ten seconds in. Both sweeps left the
    same gap: no module of `src/aacr27/` had ever been imported from a deposit in
    isolation. Every script that exercises the library also does a dozen other
    things, so a module-level dependency on an excluded file would surface as
    that script's failure and be diagnosed as that script's bug.

    The sweep came back clean, and the reason is structural rather than lucky: no
    module uses `__file__`, none has a path-valued default argument, and none has
    a module-level path constant -- every path enters as a caller-supplied
    argument. This test pins that property by exercising it, so a future module
    that resolves a repository path at import time fails here rather than in
    whichever script happens to import it first.
    """
    out = tmp_path / "snap"
    staged = subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert staged.returncode == 0, f"staging failed:\n{staged.stdout[-2000:]}"
    assert not (out / "pipeline" / "data").exists(), (
        "the deposit ships no data/; that absence is the point of this sweep")

    src = out / "pipeline" / "src"
    modules = sorted(p.stem for p in (src / "aacr27").glob("*.py"))
    # Discovered rather than hardcoded, so a new module is swept automatically --
    # but an empty or truncated glob must not pass vacuously.
    assert len(modules) >= 12, (
        f"only {len(modules)} module(s) found in the deposit's library; a sweep "
        "over an empty glob would pass while proving nothing")

    env = dict(os.environ, PYTHONPATH=str(src))
    failures = []
    for stem in modules:
        name = "aacr27" if stem == "__init__" else f"aacr27.{stem}"
        proc = subprocess.run(
            [sys.executable, "-c", f"import {name}"],
            cwd=out / "pipeline", capture_output=True, text=True, timeout=300,
            env=env)
        if proc.returncode != 0:
            failures.append(f"{name}: {proc.stderr.strip()[-400:]}")
    assert not failures, (
        "a library module cannot be imported from the public deposit, so the "
        "deposit ships code it cannot run:\n" + "\n".join(failures))


def _deposit_table_rows(readme: Path) -> list[tuple[str, int]]:
    """(command, recorded exit code) for each row of the deposit table.

    The table is `pipeline/README.md`'s "What runs from the deposit alone":
    rows of the form | `command` | **exit** | note |.
    """
    import re as _re
    text = readme.read_text(encoding="utf-8")
    start = text.index("### What runs from the deposit alone")
    end = text.index("\n### ", start + 10)
    rows = []
    for line in text[start:end].splitlines():
        m = _re.match(r"^\| `([^`]+)` \| \**(\d+)\** \|", line)
        if m:
            rows.append((m.group(1), int(m.group(2))))
    return rows


def _deposit_exit_mismatches(rows, got) -> list[str]:
    return [f"{cmd}: the table says exit {want}, the deposit gives {code}"
            for (cmd, want), code in zip(rows, got) if code != want]


def test_the_deposit_table_states_each_commands_real_exit_code(tmp_path):
    """The deposit table, run rather than carried (ledger F6.2, session 50).

    `pipeline/README.md` tells a deposit reader what exit code each command gives
    inside the published snapshot. Those codes were measured by hand on
    2026-09-08 and carried since; session 49 re-measured the table by hand, which
    is the argument for doing it mechanically. This stages a fresh snapshot with
    no `data/`, runs every listed command there, and compares.

    The pytest row is left out on purpose: running the suite inside a test of the
    suite recurses, and that row's three counts are re-measured by hand, as the
    README says. Measured cost when added: 21 s on macOS on AC power, staging
    included (2026-09-17).
    """
    rows = _deposit_table_rows(SCRIPTS.parent / "README.md")
    # VACUITY FLOOR: the table had eleven rows when this test was written.
    assert len(rows) >= 10, f"the deposit table parsed to only {len(rows)} rows"
    rows = [(cmd, code) for cmd, code in rows if "pytest" not in cmd]

    out = tmp_path / "snap"
    staged = subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert staged.returncode == 0, f"staging failed:\n{staged.stdout[-2000:]}"
    assert not (out / "pipeline" / "data").exists()

    env = dict(os.environ, MPLCONFIGDIR=str(tmp_path / "mpl"))
    got = []
    for cmd, _ in rows:
        argv = cmd.split()
        if argv[0] == "python3":
            argv[0] = sys.executable
        proc = subprocess.run(argv, cwd=out / "pipeline", capture_output=True,
                              text=True, timeout=600, env=env)
        got.append(proc.returncode)
    bad = _deposit_exit_mismatches(rows, got)
    assert not bad, (
        "pipeline/README.md's deposit table no longer states what the deposit "
        "does; re-measure the row and edit the table:\n  " + "\n  ".join(bad))

    # CAN FAIL: a table carrying one wrong code must be named.
    wrong = [(rows[0][0], rows[0][1] + 1)] + rows[1:]
    assert _deposit_exit_mismatches(wrong, got) == [
        f"{rows[0][0]}: the table says exit {rows[0][1] + 1}, "
        f"the deposit gives {got[0]}"]


# ---------------------------------------------------------------------------
# The REFERENCE axis: paths the deposit's own PROSE names.
#
# Fourth axis of the "excluded but needed" class. ARTIFACT (#38), CONSUMER (#39)
# and LIBRARY (#40) are swept. This one is the shape the class was FIRST found
# in and the only one never mechanised: the Data Availability Statement named
# files that were not in the deposit, and that was caught three times by reading
# and never by a tool. Generalised: every path a deposited file names must either
# resolve inside the deposit or say plainly that it does not.
#
# Scope is deliberately PROSE (`.md`, `.cff`) and not code. A script naming
# `07-ABSTRACT-DRAFT.md` is naming an input it scans on the working tree, and its
# behaviour in a deposit is covered by the CONSUMER sweep; a *sentence* naming it
# is an instruction to a reader who has only the deposit.
# ---------------------------------------------------------------------------

_DEPOSIT_PROSE_EXEMPT = {
    "05-PRE-REGISTRATION.md":
        "tagged protocol -- its links record what the protocol cited when it was "
        "frozen at prereg-2026-08-18, and retro-editing them would falsify the "
        "record rather than repair it",
    "14-SCIENCE-AUDIT.md":
        "chronological audit record -- its references are provenance for "
        "decisions taken on dated evidence, not instructions to a reader",
    "pipeline/results/README.md":
        "the results log classification -- each row records what a named run "
        "reported on a date, so a filename inside one is a mention in a record, "
        "not a pointer a reader is meant to follow, and the rows are protected "
        "from editing for that reason",
}
# Keyed on the snapshot-RELATIVE path, never the basename: three different
# README.md files are deposited, and a basename key would have exempted
# pipeline/README.md -- the document whose defect prompted this sweep.


def _normalise_prose(text: str) -> str:
    """Strip markdown emphasis and collapse whitespace.

    Without this the check is defeated by exactly the two things that defeated
    session 38's search: a line wrap in the middle of the phrase, and `**bold**`
    markers inside it.
    """
    import re as _re
    return " ".join(_re.sub(r"[*_`]", "", text).split())


_NOT_DEPOSITED_MARKER = (
    r"not (?:be )?part of this deposit|not in this deposit|not deposited"
    r"|private working repository|\| no \|"
)


def _undeposited_document_references(snap: Path, repo: Path) -> tuple[list[str], int]:
    """Every reference in deposited PROSE to a working document the deposit lacks.

    Returns (findings, number_of_prose_files_scanned).
    """
    import re as _re

    deposited = {p.name for p in snap.rglob("*") if p.is_file()}
    # Working documents = author-facing markdown that lives outside the pipeline.
    # Parenthesised deliberately: `-` binds tighter than `|`, so dropping these
    # brackets silently leaves every deposited root document in the set and the
    # sweep reports its own bug as nine defects. It did, on first run.
    working = sorted(
        ({p.name for p in repo.glob("*.md")}
         | {p.name for p in (repo / "submission").glob("*.md")})
        - deposited
    )
    assert len(working) >= 8, (
        f"only {len(working)} undeposited working document(s) found; this sweep "
        "would be near-vacuous -- did the document layout change?")

    pattern = _re.compile("|".join(_re.escape(n) for n in working))
    marker = _re.compile(_NOT_DEPOSITED_MARKER, _re.I)

    prose = [p for p in sorted(snap.rglob("*")) if p.is_file() and p.suffix in {".md", ".cff"}]
    findings = []
    for f in prose:
        rel = f.relative_to(snap).as_posix()
        if rel in _DEPOSIT_PROSE_EXEMPT:
            continue
        body = f.read_text(encoding="utf-8")
        lines = body.split("\n")
        # Scope the disclaimer to a WINDOW of lines around the mention, not to a
        # blank-line paragraph. Paragraph scope was the first draft and it had a
        # hole big enough to hide a real member: CITATION.cff has no blank lines,
        # so its "paragraph" was the whole 50-line header, and an unrelated
        # "private working repository" further up granted the file a blanket pass
        # -- the sweep stayed silent on a defect that was live in it. The window
        # is +/-2 lines, which is wide enough to survive the wrap that defeated
        # session 38's search and narrow enough that the disclaimer must really
        # be next to the thing it disclaims.
        for m in pattern.finditer(body):
            line = body[: m.start()].count("\n") + 1
            lo, hi = max(0, line - 1 - 2), min(len(lines), line + 2)
            if marker.search(_normalise_prose("\n".join(lines[lo:hi]))):
                continue
            findings.append(f"{rel}:{line} names {m.group(0)}")
    return findings, len(prose)


def test_no_deposited_prose_points_the_reader_at_an_undeposited_document(tmp_path):
    """A reader holding only the deposit must not be sent to a file it lacks.

    Found by sweeping this axis for the first time on 2026-09-08. Four deposited
    files -- the manuscript, `CITATION.cff`, `pipeline/README.md` and the
    generated supplementary `MANIFEST.md` -- told the reader to "see
    `16-YOUR-TASKS.md`", a private planning document that is deliberately not
    deposited and should not be. Each is now marked as undeposited where it is
    named, which is the only fix that keeps the pointer honest without publishing
    the author's task list.
    """
    out = tmp_path / "snap"
    staged = subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert staged.returncode == 0, f"staging failed:\n{staged.stdout[-2000:]}"

    repo = SCRIPTS.resolve().parents[1]
    findings, scanned = _undeposited_document_references(out, repo)
    assert scanned >= 9, (
        f"only {scanned} prose file(s) scanned in the deposit; a sweep over an "
        "empty glob would pass while proving nothing")

    # The exemptions must be REACHED, or they are dead entries that hide members.
    # Session 36 found the number checker's EXEMPT table had never executed once
    # in eleven sessions; an exemption nobody can see fire is indistinguishable
    # from a missing file.
    import re as _re
    for name, why in _DEPOSIT_PROSE_EXEMPT.items():
        staged_copy = out / name
        assert staged_copy.exists(), (
            f"{name} is exempt from this sweep but is not in the deposit at all, "
            f"so the exemption is dead. Reason on file: {why}")
        deposited = {p.name for p in out.rglob("*") if p.is_file()}
        working = ({p.name for p in repo.glob("*.md")}
                   | {p.name for p in (repo / "submission").glob("*.md")}) - deposited
        hits = _re.compile("|".join(_re.escape(n) for n in sorted(working)))
        assert hits.search(staged_copy.read_text(encoding="utf-8")), (
            f"{name} is exempt but no longer references any undeposited "
            "document, so the exemption should be removed rather than carried")

    assert not findings, (
        "deposited prose points a reader at a document the deposit does not "
        "ship, and does not say so:\n  " + "\n  ".join(findings))


def test_the_undeposited_reference_check_can_fail(tmp_path):
    """Prove the sweep discriminates.

    Direction: plant the exact sentence that was live for three sessions -- a
    bare "See `16-YOUR-TASKS.md`" with no disclaimer -- into a deposited prose
    file that is not exempt, and require the sweep to name that file.
    """
    out = tmp_path / "snap"
    subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300, check=True)

    repo = SCRIPTS.resolve().parents[1]
    victim = out / "pipeline" / "README.md"
    assert victim.exists(), "expected pipeline/README.md in the deposit"
    assert not (out / "16-YOUR-TASKS.md").exists(), (
        "16-YOUR-TASKS.md is now deposited, so it is the wrong probe for this "
        "test -- pick another undeposited working document")
    victim.write_text(
        victim.read_text(encoding="utf-8")
        + "\n\nSee `16-YOUR-TASKS.md` for the pending decision.\n",
        encoding="utf-8")

    findings, _ = _undeposited_document_references(out, repo)
    assert any("pipeline/README.md" in f for f in findings), (
        "planting an undisclosed pointer to an undeposited document did not "
        f"trip the sweep, so it cannot detect the defect it exists for: {findings}")


def _readme_document_table_rows(readme_text: str) -> list[tuple[str, bool]]:
    """Parse the `## Documents` table into (filename, claimed_deposited) rows."""
    import re as _re
    row = _re.compile(r"^\|\s*`\.\./([^`]+)`\s*\|\s*(yes|no)\s*\|", _re.M)
    return [(m.group(1), m.group(2) == "yes") for m in row.finditer(readme_text)]


def test_deposited_sibling_documents_match_the_readme_table(tmp_path):
    """`pipeline/README.md`'s Documents table must describe the real deposit.

    Until 2026-09-08 that section said all six sibling documents were
    "deliberately not part of the published pipeline" and that the paths "will
    not resolve" -- while `05-PRE-REGISTRATION.md`, `09-PAPER-DRAFT.md` and
    `14-SCIENCE-AUDIT.md` all resolved. It then told the reader to ask the author
    for the protocol, which was sitting beside that very file. Three false
    statements about the artifact a reader actually holds, in the paragraph whose
    only job is to describe it.

    The fix was a `deposited` column, and a column is only worth having if it is
    checked -- so this reads the table out of the STAGED tree and compares each
    row against what is really there.
    """
    out = tmp_path / "snap"
    staged = subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert staged.returncode == 0, f"staging failed:\n{staged.stdout[-2000:]}"

    rows = _readme_document_table_rows(
        (out / "pipeline" / "README.md").read_text(encoding="utf-8"))
    assert len(rows) >= 6, (
        f"parsed only {len(rows)} row(s) from the Documents table; the parser "
        "and the table have drifted apart and this check is near-vacuous")
    # A table that is all-yes or all-no proves nothing about the column.
    assert {claimed for _, claimed in rows} == {True, False}, (
        "the deposited column has only one value, so it cannot be shown to "
        "discriminate between a deposited and an undeposited sibling")

    wrong = [
        f"../{name}: table says deposited={claimed}, deposit says "
        f"{(out / name).exists()}"
        for name, claimed in rows if (out / name).exists() != claimed
    ]
    assert not wrong, (
        "pipeline/README.md's Documents table misdescribes the deposit:\n  "
        + "\n  ".join(wrong))


def test_the_readme_document_table_check_can_fail(tmp_path):
    """Prove the table check discriminates, in the direction the defect took.

    The real defect claimed a deposited document was absent. Reproduce that:
    flip a `yes` row to `no` in the staged copy and require a mismatch.
    """
    out = tmp_path / "snap"
    subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300, check=True)

    readme = out / "pipeline" / "README.md"
    body = readme.read_text(encoding="utf-8")
    rows = _readme_document_table_rows(body)
    truthful = [n for n, claimed in rows if claimed and (out / n).exists()]
    assert truthful, "expected at least one deposited sibling document to flip"

    name = truthful[0]
    body = body.replace(f"| `../{name}` | yes |", f"| `../{name}` | no |", 1)
    readme.write_text(body, encoding="utf-8")

    flipped = _readme_document_table_rows(body)
    wrong = [n for n, claimed in flipped if (out / n).exists() != claimed]
    assert name in wrong, (
        f"flipping {name} to `no` while it is deposited did not register as a "
        "mismatch, so the table check cannot detect the defect it exists for")


# The deposit's own root README is generated from `snapshot/README.md` and is
# the published repository's front page, not a sibling document of `pipeline/`.
# The Documents table describes `../NAME` siblings, so this one is excluded by
# kind and not by oversight.
_DEPOSIT_OWN_ROOT_README = "README.md"


def _deposited_root_documents(out: Path) -> set[str]:
    """Root-level sibling documents actually present in a staged deposit."""
    return {p.name for p in out.glob("*.md")} - {_DEPOSIT_OWN_ROOT_README}


def test_the_readme_document_table_covers_the_deposit_in_both_directions(tmp_path):
    """The Documents table must be COMPLETE, not merely accurate where it speaks.

    `test_deposited_sibling_documents_match_the_readme_table` walks the TABLE
    and checks each row against the deposit. That direction is structurally
    unable to see a document that is deposited and has no row at all: there is
    nothing to iterate over. Flipping a row is caught; never writing one is not.

    Session 45 read all 25 comparison tests and named seven whose two sides
    share an origin, so an identical defect on both passes. This is the
    independent companion for that one, and the independent property is the
    DEPOSIT'S OWN ROOT CONTENTS -- a set the table does not get a vote in.

    One end is anchored outside both sides as well. `_SUBMITTED_DOCS` is
    maintained a thousand lines above for US-spelling enforcement and knows
    nothing about staging; every root-level document it names must still appear
    in this table. Two hand-maintained sets that were written for unrelated
    reasons disagreeing is a signal; agreeing is weak evidence they are both
    right, which is the most an inexpensive check of this kind can offer.
    """
    out = tmp_path / "snap"
    staged = subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert staged.returncode == 0, f"staging failed:\n{staged.stdout[-2000:]}"

    deposited = _deposited_root_documents(out)
    # A floor: with an empty deposit every set comparison below is vacuously
    # true, which is exactly how a "found nothing" check goes quietly dead.
    assert len(deposited) >= 3, (
        f"only {len(deposited)} root-level sibling document(s) in the staged "
        "deposit; this check has gone near-vacuous and would pass over a "
        "deposit that had lost documents")

    rows = _readme_document_table_rows(
        (out / "pipeline" / "README.md").read_text(encoding="utf-8"))
    listed = {name for name, _ in rows}
    claimed_yes = {name for name, claimed in rows if claimed}

    # THE DIRECTION THE SIBLING TEST CANNOT WALK.
    unlisted = sorted(deposited - claimed_yes)
    assert not unlisted, (
        f"these documents are in the deposit but the Documents table does not "
        f"say so: {unlisted}. A reader of the published repository holds them "
        "and the table that describes the repository omits them.")

    # The other direction too, so this test stands alone if its sibling is ever
    # removed or narrowed.
    overclaimed = sorted(claimed_yes - deposited)
    assert not overclaimed, (
        f"the Documents table claims these are deposited and they are not: "
        f"{overclaimed}")

    missing_submitted = sorted(
        name for name in _SUBMITTED_DOCS if "/" not in name and name not in listed)
    assert not missing_submitted, (
        f"these submission-facing root documents have no row in the Documents "
        f"table at all: {missing_submitted}. `_SUBMITTED_DOCS` and this table "
        "are maintained separately and for different reasons; when they "
        "disagree, one of them is stale.")


def test_the_two_way_document_table_check_can_fail(tmp_path):
    """Prove the NEW direction fires, in the direction its defect would take.

    The sibling test's can-fail flips a row that exists. The defect this one
    exists for is an ABSENT row, so the mutation must delete one -- and it is
    applied to the STAGED copy in `tmp_path`, never to the repository.
    """
    out = tmp_path / "snap"
    subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300, check=True)

    readme = out / "pipeline" / "README.md"
    body = readme.read_text(encoding="utf-8")
    deposited = _deposited_root_documents(out)
    rows = _readme_document_table_rows(body)
    droppable = [n for n, claimed in rows if claimed and n in deposited]
    assert droppable, "expected a deposited row to delete"

    name = droppable[0]
    mutated = "\n".join(
        ln for ln in body.splitlines() if f"| `../{name}` |" not in ln)
    assert mutated != body, f"deleting the {name} row changed nothing"

    listed_after = {n for n, claimed in _readme_document_table_rows(mutated) if claimed}
    assert name not in listed_after, "the mutation did not remove the row"
    # The positive control: the real check must now report exactly that name.
    assert sorted(deposited - listed_after) == [name], (
        f"deleting the {name} row from the table did not register as an "
        "uncovered deposited document, so the two-way check cannot detect the "
        "defect it exists for")


# ---------------------------------------------------------------------------
# Delimiter balance: the cheapest detector for a HALF-APPLIED EDIT.
#
# Session 40 found a manuscript sentence reading "...which is the defect the
# CLOSED note records. Re-measure both / rebuild figures are themselves already
# superseded -- see below)" -- an edit had replaced the first half of a sentence
# and taken the rest of the clause with it, leaving an orphan `)`. It had been
# ungrammatical for a full session. The number checker cannot see it (no number
# changed), --coverage cannot see it, and a spell-checker cannot see it. An
# unmatched bracket can, and it costs milliseconds.
# ---------------------------------------------------------------------------

# `{` joined on 2026-09-17 (F8.2): measured at zero hits on the real tree first.
_BALANCE_PAIRS = {"(": ")", "[": "]", "{": "}"}


# Documents deliberately outside every prose scan, with the reason. Kept as a
# named set rather than an omission so `test_every_prose_document_is_scanned`
# can prove the scanned set and this one together cover the tree.
_PROSE_NOT_SCANNED = {
    "bagaev.html": "a saved third-party web page, not project prose",
}


def _prose_documents(repo: Path) -> list[Path]:
    """Author-facing prose. Handoffs are excluded: they are archives.

    The `extra` list is hand-maintained, which is the exposure session 45
    measured: it covered 30 of the 34 tracked prose documents, and one of the
    four it missed was `poster/poster.html` -- a SUBMITTED artefact, already
    named in `_SUBMITTED_DOCS` two thousand lines above. Two hand-maintained
    document sets in one file, disagreeing, with nothing to notice.
    `test_every_prose_document_is_scanned` now reconciles this list against the
    tree, so a new document has to be classified rather than silently skipped.
    Noise from the three additions was measured before they were added: zero
    hits from the bracket scan and zero from the bold scan, on all three.
    """
    docs = [p for p in sorted(repo.glob("*.md")) if "HANDOFF" not in p.name]
    docs += sorted((repo / "submission").glob("*.md"))
    for extra in ("CITATION.cff", "snapshot/README.md",
                  "pipeline/README.md", "pipeline/results/README.md",
                  "pipeline/results/supplementary/CAPTIONS.md",
                  "pipeline/results/supplementary/MANIFEST.md",
                  "pipeline/ENVIRONMENT.md",
                  "poster/poster.html", "poster/README.md",
                  "pipeline/cluster/a9-investigation-scratch/README.md"):
        p = repo / extra
        if p.exists():
            docs.append(p)
    return docs


def _unbalanced_delimiters(text: str) -> list[tuple[int, str]]:
    """[(line, why)] for brackets that never pair up, ignoring code and HTML.

    Code is skipped because `results/README.md` quotes shell and Python that is
    legitimately unbalanced when sliced by line, and HTML comments are skipped
    because four `numcheck: ignore` markers live in them.
    """
    import re as _re
    closers = {v: k for k, v in _BALANCE_PAIRS.items()}
    text = _re.sub(r"<!--.*?-->", "", text, flags=_re.S)
    fence = _re.compile(r"^\s*```")
    inline = _re.compile(r"`[^`]*`")

    out: list[tuple[int, str]] = []
    stack: list[tuple[str, int]] = []
    in_fence = False
    for lineno, raw in enumerate(text.split("\n"), start=1):
        if fence.match(raw):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        for ch in inline.sub("", raw):
            if ch in _BALANCE_PAIRS:
                stack.append((ch, lineno))
            elif ch in closers:
                if stack and stack[-1][0] == closers[ch]:
                    stack.pop()
                else:
                    out.append((lineno, f"closing {ch!r} with no opener"))
    out.extend((lineno, f"opening {ch!r} never closed") for ch, lineno in stack)
    return out


def _unpaired_marks(text: str) -> list[tuple[int, str]]:
    """[(line, why)] for paragraphs whose paired marks do not pair (F8.2).

    Backticks, straight double quotes and curly double quotes are checked per
    PARAGRAPH, not per line: an inline code span or a quotation legitimately
    wraps across a line break, and per line the backtick check alone produced 22
    false hits on the real tree where per paragraph it produces none. Measured
    2026-09-17 over the 33 scanned documents before adding it: zero hits for
    each of the three marks. Code fences and HTML comments are skipped, as in
    `_unbalanced_delimiters`.
    """
    import re as _re
    text = _re.sub(r"<!--.*?-->", lambda m: "\n" * m.group(0).count("\n"),
                   text, flags=_re.S)
    out: list[tuple[int, str]] = []
    para: list[str] = []
    start, in_fence = 1, False

    def flush() -> None:
        body = "\n".join(para)
        if body.count("`") % 2:
            out.append((start, "odd number of backticks"))
        prose = _re.sub(r"`[^`]*`", "", body)
        if prose.count('"') % 2:
            out.append((start, "odd number of straight double quotes"))
        if prose.count("\u201c") != prose.count("\u201d"):
            out.append((start, "curly double quotes do not pair"))

    for lineno, raw in enumerate(text.split("\n"), start=1):
        if _re.match(r"^\s*```", raw):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        if raw.strip():
            if not para:
                start = lineno
            para.append(raw)
        elif para:
            flush()
            para = []
    if para:
        flush()
    return out


def test_prose_documents_have_balanced_delimiters():
    """A half-applied edit usually orphans a bracket. Catch it mechanically."""
    repo = Path(__file__).resolve().parents[2]
    docs = _prose_documents(repo)
    assert len(docs) >= 20, (
        f"only {len(docs)} prose document(s) discovered; the sweep that found "
        "the session-40 defect covered twenty, so this is near-vacuous")

    bad = {}
    for d in docs:
        text = d.read_text(encoding="utf-8")
        hits = _unbalanced_delimiters(text) + _unpaired_marks(text)
        if hits:
            bad[str(d.relative_to(repo))] = hits[:8]
    assert not bad, (
        "unbalanced delimiters outside code -- usually a sentence half-replaced "
        "by an edit that took the rest of the clause with it:\n"
        + "\n".join(f"  {f}: {h}" for f, h in bad.items()))


def test_the_delimiter_balance_check_can_fail():
    """Prove the scan discriminates, using the real defect's shape.

    Verified against history as well as against this planted case: run the same
    function over `git show 7ab22cb^:09-PAPER-DRAFT.md` and it reports
    `line 906: closing ')' with no opener`, which is the session-40 defect. That
    is left out of the test itself because pinning a test to a commit hash makes
    it fail for reasons that have nothing to do with the property being checked.
    """
    broken = (
        "A sentence that is fine.\n\n"
        "Some prose which is the defect the CLOSED note records. Re-measure "
        "both / rebuild figures are themselves already superseded -- see "
        "below).\n"
    )
    hits = _unbalanced_delimiters(broken)
    assert hits == [(3, "closing ')' with no opener")], (
        f"the balance scan did not catch an orphan closing paren: {hits}")

    # And it must NOT fire on the things prose legitimately contains.
    for benign in (
        "A citation [1] and a CI [0.2692, 0.3133] and a link [text](url).\n",
        # A fence only opens at line start, which is why this fixture puts it
        # there -- the first draft inlined it and the scan correctly read the
        # following line as prose.
        "Inline code may be unbalanced: `f(x` here.\n\n```\nfoo(bar\n```\n",
        "An HTML comment <!-- numcheck: ignore: a ) inside --> is skipped.\n",
    ):
        assert not _unbalanced_delimiters(benign), (
            f"the balance scan fires on correct prose, which is how the "
            f"previous WORD_CLAIMS attempt died: {benign!r}")

    # F8.2's marks, both directions: each half-applied mark is caught, and the
    # shapes prose legitimately contains are not.
    for broken_mark, why in (
        ("An edit that dropped the close of `a code span.\n", "odd number of backticks"),
        ('He said "half a quotation and stopped.\n', "odd number of straight double quotes"),
        ("A \u201ccurly quotation that never closed.\n", "curly double quotes do not pair"),
    ):
        assert _unpaired_marks(broken_mark) == [(1, why)], (
            f"the mark scan missed {why!r}: {_unpaired_marks(broken_mark)}")
    assert _unbalanced_delimiters("A dict {k: v and nothing else.\n") == [
        (1, "opening '{' never closed")], "an orphan brace was not caught"
    for benign in (
        "A span `that wraps\nacross a line` is one span.\n",
        'A "quotation that\nwraps" is one quotation.\n',
        "A \u201ccurly one\u201d and a code span with a quote `\"x` inside.\n",
    ):
        assert not _unpaired_marks(benign), (
            f"the mark scan fires on correct prose: {benign!r}")


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


def _files_the_checker_reads(monkeypatch) -> tuple[object, set[str]]:
    """Every file under results/ that `19_check_numbers.authority()` opens.

    Measured, not declared: every way the checker reads a file (builtins.open,
    Path.read_text/read_bytes/open, np.load, pd.read_csv, pd.read_parquet) is
    wrapped for the one call, and each path that resolves inside results/ is
    recorded. A declared list would be a second definition of the authority set,
    free to disagree with the first, and every defect of this class so far has
    lived in exactly that gap.
    """
    import builtins

    checker = _load_script("19_check_numbers.py")
    res = checker.RESULTS.resolve()
    seen: set[str] = set()

    def record(p):
        try:
            q = Path(p).resolve()
        except TypeError:        # a file descriptor or an open handle
            return
        if res in q.parents:
            seen.add(q.relative_to(res).as_posix())

    def wrap(real, method=False):
        if method:
            return lambda self, *a, **k: (record(self), real(self, *a, **k))[1]
        return lambda f, *a, **k: (record(f), real(f, *a, **k))[1]

    monkeypatch.setattr(builtins, "open", wrap(builtins.open))
    for name in ("read_text", "read_bytes", "open"):
        monkeypatch.setattr(Path, name, wrap(getattr(Path, name), method=True))
    for owner, name in ((np, "load"), (pd, "read_csv"), (pd, "read_parquet")):
        monkeypatch.setattr(owner, name, wrap(getattr(owner, name)))
    try:
        checker.authority(slow=False)
    finally:
        monkeypatch.undo()
    seen.discard("PROVENANCE.json")
    return checker, seen


def test_every_file_the_checker_reads_as_an_authority_is_hashed(monkeypatch):
    """The checker's authority set must sit inside PROVENANCE.json (F3.4, session 50).

    `19_check_numbers.py` compares prose with artefacts and cannot see an edited
    artefact by construction; `21_provenance_manifest.py` exists to catch that,
    but only for the files it tracks. Nothing reconciled the two sets, and when
    session 50 measured them, six files the checker reads as authorities (E6's
    PLAGE summary, E8's omega summary, E16's four subsample-AUROC records) were
    not hashed: session 49 hashed its authorities before those three analyses
    landed. An authority that is not hashed can drift with the checker green.

    The one-way direction is the one that matters. PROVENANCE also hashes files
    no document quotes (figure inputs, the frozen primaries' secondary tables),
    by design. `every_stored_run()`'s reads are deliberately excluded: that is
    the whitelist of values any stored run produced, which includes the synthetic
    `demo/` run, and it authorises nothing on its own.
    """
    import json as _json

    checker, seen = _files_the_checker_reads(monkeypatch)
    # POSITIVE CONTROLS on the instrument: files the checker is known to read,
    # through three different read paths. If these vanish, the wrapping broke,
    # and an empty `seen` would make the assertion below pass on nothing.
    for known in ("nsclc_v3/summary.json", "pancancer_v3/per_signature.csv",
                  "sort_fix_fold_change.json"):
        if (checker.RESULTS / known).exists():
            assert known in seen, f"the read instrument missed {known}; it is blind"
    # VACUITY FLOOR. Measured 2026-09-17 on the working tree: 55 files, 53 of
    # them deposited (the two bisection .npz files are not). A public deposit
    # therefore reads about 53; the floor sits a little under each.
    floor = 50 if _ON_WORKING_TREE else 45
    assert len(seen) >= floor, (
        f"the checker read only {len(seen)} files under results/ (floor {floor}); "
        "either the authority table shrank or the instrument stopped seeing reads")

    hashed = set(_json.loads((checker.RESULTS / "PROVENANCE.json").read_text())["files"])
    unhashed = sorted(seen - hashed)
    assert not unhashed, (
        f"{len(unhashed)} file(s) are read by 19_check_numbers.py as an authority "
        "but not hashed by 21_provenance_manifest.py, so an edit to them would go "
        "unnoticed. Add each to TRACKED_DIRS or TRACKED_FILES and rerun --write "
        "in its own commit:\n  " + "\n  ".join(unhashed))

    # CAN FAIL: the same comparison against a manifest missing one real
    # authority must name it. Proves the set difference is not vacuous.
    victim = sorted(seen & hashed)[0]
    assert sorted(seen - (hashed - {victim})) == [victim]


def test_a_claim_whose_authority_vanished_is_reported_not_skipped():
    """`--strict` must see a claim row whose authority key is gone (F3.4, session 50).

    `scan` skips a claim with `if key not in table: continue`, and for CLAIMS
    that happens before the match count is recorded, so the vacuity check never
    heard of the row: an authority that vanished (a renamed JSON field read
    under an `if`, a loader block deleted) left its pattern checking nothing,
    silently. `_absent_claim_keys` is what `--strict` now fails on.
    """
    checker = _load_script("19_check_numbers.py")
    table = checker.authority(slow=False)
    # The fast table lacks exactly the live-measured keys, by design. That they
    # show up is the positive control: the helper really does see absences.
    live = {"abstract_chars", "abstract_headroom", "pancancer_n_sites",
            "snapshot_files", "test_count", "abstract_title_chars",
            "abstract_naive_len", "abstract_body_words",
            # Session 57: A16's extent, recomputed from the interim expression
            # index in the slow path beside the site count (`_a16_extent`).
            "a16_ns_not_tumour", "a16_ns_adjacent_normal", "a16_ns_recurrence",
            "a16_ns_normal_only", "a16_pan_averaged", "a16_pan_no_tumour",
            "a16_pan_total",
            # Session 58: the audit block's breakdown of those counts.
            "a16_ns_single_tumour", "a16_ns_tumour_normal",
            "a16_pan_single_tumour", "a16_pan_tumour_normal",
            "a16_pan_tumour_recurrence", "a16_pan_tumour_metastasis",
            "a16_pan_multi_other", "a16_pan_normal_only",
            "a16_pan_recurrence_only", "a16_pan_metastasis_only",
            "a16_pan_type_thca", "a16_pan_type_kirc", "a16_pan_type_prad",
            "a16_pan_type_luad", "a16_pan_type_lihc",
            # Session 58: the dry run's percentages, derived from the live
            # abstract recount beside the headroom.
            "abstract_pct_of_limit", "abstract_naive_excess_pct",
            # Session 58: live beside pancancer_n_sites since 2026-09-04, and
            # quoted by a pattern for the first time (the "of the 68 NSCLC
            # sites" sentences), so its absence from the fast table now shows.
            "nsclc_n_sites",
            # Session 60: E15's site attrition, beside the site counts.
            "e15_ns_sites_dropped", "e15_ns_patients_dropped",
            "e15_pan_sites_dropped", "e15_pan_patients_dropped", "e15_pan_retained"}
    absent = set(checker._absent_claim_keys(table))
    assert absent, "the fast table should lack the live keys; the helper is blind"
    assert absent <= live, (
        f"claim keys with no authority value: {sorted(absent - live)}. Their "
        "patterns are skipped without being counted, so they check nothing.")
    # CAN FAIL: remove a real authority and the helper must name it.
    victim = "sortaudit_isi"
    assert victim in table
    del table[victim]
    assert victim in checker._absent_claim_keys(table)


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


def _render_supplement_in_a_staged_tree(tmp_path):
    """Stage a deposit and re-render the supplement inside it.

    Returns the directory holding the freshly rendered files. Staging rather than
    copying the repository is deliberate: `22_build_public_snapshot.py` copies
    from the WORKING TREE, so an edited-but-not-re-rendered generator is visible
    here, which is precisely the defect this axis exists to catch. It also means
    the builder writes into the temp tree and never into `results/`.
    """
    out = tmp_path / "snap"
    staged = subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"),
         "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300)
    assert staged.returncode == 0, f"staging failed:\n{staged.stdout[-2000:]}"

    built = subprocess.run(
        [sys.executable, str(out / "pipeline" / "scripts" / "25_build_supplementary.py")],
        cwd=out / "pipeline", capture_output=True, text=True, timeout=600)
    assert built.returncode == 0, (
        f"the supplement did not render:\n{built.stdout[-3000:]}\n{built.stderr[-3000:]}")
    return out / "pipeline" / "results" / "supplementary"


def _generated_artefact_drift(rendered):
    """Names of committed supplementary files that a fresh render does not reproduce."""
    committed = SCRIPTS.parent / "results" / "supplementary"
    drifted = []
    compared = 0
    for path in sorted(committed.glob("*.tsv")) + sorted(committed.glob("*.md")):
        fresh = rendered / path.name
        assert fresh.exists(), f"{path.name} is committed but a fresh render did not produce it"
        compared += 1
        if fresh.read_bytes() != path.read_bytes():
            drifted.append(path.name)
    # A glob that matched nothing would report "no drift" while proving nothing --
    # exactly how 19_check_numbers.py v1 reported "62 claims OK" while checking
    # nothing. Nine tables plus MANIFEST.md and CAPTIONS.md is the floor.
    assert compared >= 11, (
        f"only {compared} generated artefact(s) were compared; the sweep has "
        "gone vacuous and would pass over a stale supplement")
    return drifted


def test_committed_generated_artefacts_match_a_fresh_render(tmp_path):
    """The sixth axis of the "producer and consumer disagree" class: GENERATED.

    Five axes are swept. Session 38 took ARTIFACT (staged tree vs published
    tree), 39 took CONSUMER (every script run inside a staged tree), 40 took
    LIBRARY (every module imported from a deposit), 41 took REFERENCE (every path
    the deposit's prose names) and sized POINTER (`file.py:NNN`). All five ask
    whether some consumer can still read what the repository ships. None asks the
    prior question: is the committed file even what its own generator produces?

    The gap is not hypothetical. Session 41 edited `25_build_supplementary.py`
    for a REFERENCE-axis fix and had to REMEMBER to re-render by hand; nothing in
    the suite would have caught a forgotten render, and the eleven files here are
    the supplement the manuscript cites by S-number. `ENVIRONMENT.md` has had
    such a test since session 35 and `PROVENANCE.json` is covered by `--write`;
    these eleven and `results/gene_set_sizes.csv` had nothing.

    `gene_set_sizes.csv` is deliberately NOT checked here: regenerating it needs
    `data/`, which the deposit does not ship and a clean checkout does not have.
    It was verified by hand against a fresh render on 2026-09-08 and is hashed in
    `PROVENANCE.json`, which detects an edit even though it cannot detect staleness.
    """
    drifted = _generated_artefact_drift(_render_supplement_in_a_staged_tree(tmp_path))
    assert not drifted, (
        "these committed files are NOT what their generator now produces: "
        f"{drifted}. Edit `25_build_supplementary.py` and re-render -- never "
        "hand-edit the generated file.")


def test_the_generated_artefact_render_check_can_fail(tmp_path):
    """Prove the sweep above discriminates, in the direction the real defect takes.

    The defect is never a corrupted output file; it is an edited GENERATOR whose
    committed artefact was not re-rendered. So the mutation edits the generator
    inside the staged tree and re-renders, rather than perturbing a rendered file
    -- which would only prove that `!=` works on bytes.
    """
    out = tmp_path / "snap"
    subprocess.run(
        [sys.executable, str(SCRIPTS / "22_build_public_snapshot.py"), "--out", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=300, check=True)

    builder = out / "pipeline" / "scripts" / "25_build_supplementary.py"
    original = builder.read_text()
    marker = "Genetic-ancestry composition of the pan-TCGA cohort and the "
    assert marker in original, (
        "the caption this mutation edits has been reworded; point it at a "
        "string that still exists, or the can-fail test is dead code")
    builder.write_text(original.replace(marker, "Ancestry composition (mutated) and the "))

    built = subprocess.run(
        [sys.executable, str(builder)],
        cwd=out / "pipeline", capture_output=True, text=True, timeout=600)
    assert built.returncode == 0, "the mutated generator should still render, just differently"

    drifted = _generated_artefact_drift(out / "pipeline" / "results" / "supplementary")
    assert "CAPTIONS.md" in drifted, (
        "editing a caption in the generator did not make the committed "
        "CAPTIONS.md read as stale, so the sweep cannot detect the defect it "
        f"exists for (drift reported: {drifted})")


# ---------------------------------------------------------------------------
# The supplement's SHAPES, checked against the data model rather than against
# the dataframe that produced them.
#
# Session 45 read all 25 comparison tests in this file and split them in two.
# Eight are PROTECTED: they compare against an external or independently derived
# authority -- the literature, a naive reference loop, algebra, a permutation
# null. Seven are STRUCTURALLY BLIND: both sides share an origin, so a defect
# present in both passes unnoticed. Session 43 gave one of the seven an
# independent companion: the committed figures are compared to a fresh render
# AND, separately, asserted to embed Type 42 fonts, which is a property no fresh
# render can launder. The other six were left without one.
#
# This is the companion for the supplement, and the blindness is exact and
# citable. `25_build_supplementary.py` prints the shape into MANIFEST.md and into
# CAPTIONS.md from `len(df)` and `len(df.columns)` of the SAME dataframe in the
# SAME loop iteration. The two files therefore agree by construction; they
# cannot disagree, whatever the builder does. And
# `test_committed_generated_artefacts_match_a_fresh_render` compares committed
# bytes against freshly rendered bytes, which shares that origin too. A builder
# that dropped a cohort, lost a signature or duplicated every row would produce
# a table, a manifest and a caption that all agreed with each other and with a
# fresh render -- and every existing test would pass.
#
# The independent authority is the data model the manuscript states: sixteen
# Hallmark signatures, two cohorts, the partition counts the variance analysis
# reports, two platforms. None of those numbers is read from a dataframe here.
#
# Three tables have no derivable shape and are recorded as such WITH a reason
# rather than omitted, because an omission is precisely what makes an incomplete
# set look finished -- the ninth axis, UNCLASSIFIED-MEMBER. The reconciliation
# in the test below is what stops this set acquiring one.
_SIGNATURES = 16                        # the Hallmark panel; "16/16 beat their null"
_COHORTS = ("NSCLC", "pan-TCGA")

_SUPPLEMENT_EXPECTED_ROWS = {
    "S1_per_signature_excess": (
        _SIGNATURES * len(_COHORTS), "16 Hallmark signatures x 2 cohorts"),
    "S2_reliability_alpha": (
        _SIGNATURES * len(_COHORTS), "16 Hallmark signatures x 2 cohorts"),
    "S9_decomposition": (
        _SIGNATURES * len(_COHORTS), "16 Hallmark signatures x 2 cohorts"),
    "S4_label_site_variance": (
        _SIGNATURES,
        "16 Hallmark signatures x ONE cohort -- Control C's label-side variance "
        "is a pan-TCGA analysis only, which is why this table is half the size "
        "of S1 rather than a sign that a cohort went missing"),
    "S6_partition_variance": (
        (6 + 5) + (24 + 24),
        "the FROZEN sweep: 6 NSCLC partitions, all usable, plus 6 pan-cancer "
        "partitions of which 5 are usable (split_seed=5 stopped on 'SVD did not "
        "converge' in the covariate-baseline regression -- session 48 traced "
        "it; it is not the degeneracy guard) -- and the STABLE24 sweep, 24 "
        "partitions per cohort, every one computed (ledger B3). The frozen run "
        "is what it is, so 12 frozen rows would mean the frozen table changed"),
    "S10_platform_reproducibility": (
        2, "2 platforms: macOS arm64 and Linux x86_64"),
    # Added 2026-09-16 with the table itself. S5 is S4's NSCLC counterpart, so
    # it is the SAME shape as S4 and for the same reason -- one cohort, not two.
    # If this ever equals 32 the two Control C tables have been merged and one
    # of these two entries is now wrong.
    "S5_label_site_variance_nsclc": (
        _SIGNATURES,
        "16 Hallmark signatures x ONE cohort -- the NSCLC counterpart to S4, "
        "which is pan-TCGA. Both are half the size of S1 because Control C is "
        "reported per cohort in its own table rather than stacked"),
}

_SUPPLEMENT_ROWS_NOT_DERIVABLE = {
    "S8_ancestry":
        "key-value scalars; the row count is however many quantities the "
        "ancestry summary chooses to report, not a product of anything",
    "S3_small_panel_bracket":
        "one row per panel size k, and the set of k values lives in the "
        "builder -- deriving it from the builder would be the same blindness "
        "this check exists to remove",
    "S7_negative_controls":
        "one row per (cohort, control, site); the per-cohort site counts are "
        "data-dependent (NSCLC contributes two controls, pan-TCGA one) and "
        "cannot be derived without reading the files being checked",
}

_COMMITTED_SUPPLEMENT = SCRIPTS.parent / "results" / "supplementary"


def _tsv_shape(path):
    """(rows excluding header, columns) for a committed supplementary table."""
    lines = path.read_text().splitlines()
    return len(lines) - 1, len(lines[0].split("\t"))


def _caption_shapes():
    """S-number -> (rows, columns) as CAPTIONS.md states them, in prose."""
    text = (_COMMITTED_SUPPLEMENT / "CAPTIONS.md").read_text()
    out, current = {}, None
    for line in text.splitlines():
        head = re.match(r"\*\*Supplementary Table (S\d+)\.\*\*", line)
        if head:
            current = head.group(1)
        shape = re.match(r"\*(\d+) rows x (\d+) columns\.", line)
        if shape and current:
            out[current] = (int(shape.group(1)), int(shape.group(2)))
            current = None
    return out


@pytest.mark.skipif(
    not (_COMMITTED_SUPPLEMENT / "S1_per_signature_excess.tsv").exists(),
    reason="the committed .tsv tables are not staged into the deposit "
           "(.tsv is not in RESULT_KEEP_SUFFIX); a deposit reader regenerates "
           "them, and that path is covered by the snapshot build test")
def test_supplementary_shapes_match_an_independently_derived_expectation():
    """Every caption's stated shape, checked against the data model.

    This is the independent companion the supplement's blind sync tests lack.
    See the comment above for why comparing CAPTIONS.md to MANIFEST.md, or
    either to a fresh render, cannot see a wrong shape.
    """
    committed = sorted(p.name[:-4] for p in _COMMITTED_SUPPLEMENT.glob("*.tsv"))
    assert len(committed) >= 10, (
        f"only {len(committed)} supplementary table(s) found; the sweep has "
        "gone vacuous and would pass over a supplement that had lost tables")

    # UNCLASSIFIED-MEMBER, both directions: the IN list and the OUT list must
    # together cover the universe, and neither may name a table that is gone.
    classified = set(_SUPPLEMENT_EXPECTED_ROWS) | set(_SUPPLEMENT_ROWS_NOT_DERIVABLE)
    unclassified = sorted(set(committed) - classified)
    assert not unclassified, (
        f"these supplementary tables are in neither the derivable set nor the "
        f"documented not-derivable set: {unclassified}. Add each to exactly one, "
        "with a reason -- a set that documents its exclusions looks finished, "
        "which is how 25 of 34 prose documents went unguarded until session 45.")
    dead = sorted(classified - set(committed))
    assert not dead, (
        f"these entries name supplementary tables that no longer exist: {dead}")
    overlap = sorted(set(_SUPPLEMENT_EXPECTED_ROWS) & set(_SUPPLEMENT_ROWS_NOT_DERIVABLE))
    assert not overlap, f"classified both ways: {overlap}"

    captions = _caption_shapes()
    wrong = []
    for name, (expected, why) in sorted(_SUPPLEMENT_EXPECTED_ROWS.items()):
        rows, cols = _tsv_shape(_COMMITTED_SUPPLEMENT / f"{name}.tsv")
        if rows != expected:
            wrong.append(f"{name}: {rows} rows, but the data model gives "
                         f"{expected} ({why})")
        # The caption is what a reviewer reads. A table that is right while its
        # caption describes a different one is the same defect to that reader.
        num = name.split("_")[0]
        assert num in captions, f"CAPTIONS.md states no shape for {num}"
        if captions[num] != (rows, cols):
            wrong.append(f"{name}: CAPTIONS.md says {captions[num]}, the file "
                         f"is ({rows}, {cols})")
    assert not wrong, "supplementary shapes disagree with the data model:\n  " + \
        "\n  ".join(wrong)

    # Row COUNT alone is weak: 32 rows could be 32 signatures in one cohort.
    # Assert the factorisation, which is what the caption actually claims.
    for name in ("S1_per_signature_excess", "S2_reliability_alpha"):
        lines = (_COMMITTED_SUPPLEMENT / f"{name}.tsv").read_text().splitlines()
        cohorts = {ln.split("\t")[0] for ln in lines[1:]}
        signatures = {ln.split("\t")[1] for ln in lines[1:]}
        assert cohorts == set(_COHORTS), (
            f"{name} covers cohorts {sorted(cohorts)}, not {list(_COHORTS)}")
        assert len(signatures) == _SIGNATURES, (
            f"{name} covers {len(signatures)} signatures, not {_SIGNATURES}")


@pytest.mark.skipif(
    not (_COMMITTED_SUPPLEMENT / "S1_per_signature_excess.tsv").exists(),
    reason="see the test above")
def test_the_supplementary_shape_check_can_fail():
    """Prove the check discriminates, in the direction the real defect takes.

    The real defect is a builder that emits the wrong number of rows while every
    file that describes it agrees, so the mutation moves the EXPECTATION -- the
    only side that is independent -- and the check must notice. It mutates the
    dict IN MEMORY and restores it: nothing is written to the repository, which
    is the pattern session 45 adopted after session 44's mutation tests wrote
    into `scripts/` and `results/` and relied on a `finally:` to undo it.
    """
    original = dict(_SUPPLEMENT_EXPECTED_ROWS)
    try:
        _SUPPLEMENT_EXPECTED_ROWS["S6_partition_variance"] = (
            12 + 48, "as if the frozen split_seed=5 partition had been computed")
        with pytest.raises(AssertionError, match="disagree with the data model"):
            test_supplementary_shapes_match_an_independently_derived_expectation()
    finally:
        _SUPPLEMENT_EXPECTED_ROWS.clear()
        _SUPPLEMENT_EXPECTED_ROWS.update(original)

    # And the other direction: an unclassified table must fail too.
    try:
        removed = _SUPPLEMENT_EXPECTED_ROWS.pop("S10_platform_reproducibility")
        with pytest.raises(AssertionError, match="neither the derivable set"):
            test_supplementary_shapes_match_an_independently_derived_expectation()
    finally:
        _SUPPLEMENT_EXPECTED_ROWS["S10_platform_reproducibility"] = removed


def test_abstract_counter_applies_the_aacr_counting_rules():
    """`08_count_abstract.py` gates the 2,600-character AACR limit and had no test.

    Session 41's triage of the thirteen unmentioned `scripts/` put this one in the
    "pure, cheap, and genuinely untested" bucket and called it the cheapest
    high-value test available: the abstract sits at 2,515 with 85 characters of
    headroom, and going over is HARD-BLOCKED at submission rather than trimmed.

    The rules are counter-intuitive in three ways, and each is asserted here
    rather than described: spaces are free, a table costs a flat 800 whatever it
    contains, and the NextGen variant uses the OPPOSITE space convention.
    """
    m = _load_script("08_count_abstract.py")

    # Spaces are free. This is the rule that makes a naive len() overcount by
    # ~17% and provoke needless cutting of an abstract that already fits.
    assert m.count("a b c") == 3
    assert m.count("  lots   of\n\tspace  ") == len("lotsofspace")
    assert m.count("") == 0

    # The limit and the flat table cost are the two constants a future edit is
    # most likely to get wrong; pin them to the Call for Abstracts.
    assert m.LIMIT == 2600
    assert m.TABLE_COST == 800
    assert m.NEXTGEN_LIMIT == 8000

    # extract() must take TITLE and BODY and leave working notes out -- the
    # draft file carries budget notes and TODOs in the same document.
    draft = ("## TITLE\nA Title Here\n\n## BODY\nBody text.\n\n"
             "## NOTES\nDo not count this sentence.\n")
    title, body = m.extract(draft)
    assert title == "A Title Here"
    assert "Body text." in body
    assert "Do not count" not in body, "working notes leaked into the counted body"

    # Bold markers are formatting for reading the draft, not submitted characters.
    _, bolded = m.extract("## TITLE\nT\n\n## BODY\n**Background:** x\n")
    assert "**" not in bolded

    # The draft closes its body with a `---` rule before the notes, and the rule
    # is never pasted. Until 2026-09-17 it was counted: every reported total was
    # 3 characters high (2,515 for a pasted 2,512). Dashes INSIDE the body stay.
    _, ruled = m.extract("## TITLE\nT\n\n## BODY\nBody text.\n\n---\n\n## NOTES\nx\n")
    assert ruled == "Body text.", f"a closing rule was counted as body: {ruled!r}"
    _, dashed = m.extract("## TITLE\nT\n\n## BODY\nA -- dash, and 3---4 inside.\n")
    assert dashed == "A -- dash, and 3---4 inside."


def test_abstract_counter_reports_a_missing_draft_instead_of_crashing(tmp_path, capsys):
    """The exact defect session 39 found by running the deposit's own entry point.

    `reproduce.sh --check` calls this script, and the public snapshot deliberately
    does NOT stage `07-ABSTRACT-DRAFT.md`. Every reader running the advertised
    command therefore reaches this path. Measured inside a staged snapshot on
    2026-09-07 it died with a raw `FileNotFoundError` traceback, which reads like
    a broken deposit and is not one. Exit 2 plus an explanation is the contract.
    """
    m = _load_script("08_count_abstract.py")
    missing = tmp_path / "07-ABSTRACT-DRAFT.md"
    assert not missing.exists()

    argv = sys.argv
    sys.argv = ["08_count_abstract.py", str(missing)]
    try:
        rc = m.main()
    finally:
        sys.argv = argv

    assert rc == 2, f"a missing draft must return 2, not {rc} (and never a traceback)"
    out = capsys.readouterr().out
    assert "MISSING" in out
    assert "not part of the public snapshot" in out, (
        "the message must tell a deposit reader that nothing is wrong with the "
        "deposit, or the exit code just relocates the confusion")


def test_abstract_counter_fails_when_the_abstract_is_over_the_limit(tmp_path, capsys):
    """Prove the gate can go red -- otherwise it is decoration on a hard block.

    Two directions: over by content, and over by TABLES ALONE. The second matters
    because tables cost a flat 800 regardless of size, so four tables overrun the
    limit on their own and a counter that ignored `--tables` would still look
    green on a fitting body.
    """
    m = _load_script("08_count_abstract.py")
    draft = tmp_path / "a.md"

    over = "x" * (m.LIMIT + 1)
    draft.write_text(f"## TITLE\nT\n\n## BODY\n{over}\n")
    argv = sys.argv
    try:
        sys.argv = ["08_count_abstract.py", str(draft)]
        assert m.main() == 1, "an over-length abstract did not fail"
        assert "OVER by" in capsys.readouterr().out

        # A body that fits, made over-length by four tables at 800 each.
        draft.write_text("## TITLE\nT\n\n## BODY\nshort body\n")
        sys.argv = ["08_count_abstract.py", str(draft)]
        assert m.main() == 0, "a short abstract with no tables should pass"
        sys.argv = ["08_count_abstract.py", str(draft), "--tables", "4"]
        assert m.main() == 1, (
            "four tables cost 3,200 against a 2,600 limit; --tables is not "
            "being charged at all")
    finally:
        sys.argv = argv


def _normalise_svg(text: str) -> str:
    """Strip only IDENTIFIERS and the render timestamp from an SVG.

    matplotlib names clip paths and marker defs with a per-render hash, so two
    byte-different SVGs are routinely identical drawings. Session 42 got this
    wrong TWICE before getting it right -- one regex reported four of five files
    clean and one stale, the next reported all five stale -- because `\\b` does
    not match after the underscore in `C0_0_b11c0c70ad`. Matching identifier
    ATTRIBUTES structurally is reliable where matching their spelling is not,
    and the can-fail companion proves this hides no geometry, colour or style.
    """
    text = re.sub(r'<dc:date>[^<]*</dc:date>', '<dc:date/>', text)
    text = re.sub(r'id="[^"]*"', 'id="X"', text)
    text = re.sub(r'xlink:href="#[^"]*"', 'xlink:href="#X"', text)
    return re.sub(r'url\(#[^)]*\)', 'url(#X)', text)


def _strip_pdf_date(blob: bytes) -> bytes:
    """Drop the PDF creation timestamp, which changes on every render."""
    return re.sub(rb'/CreationDate \([^)]*\)', b'/CreationDate ()', blob)


def _render_figures(tmp_path):
    """Render all five figures in all three formats into a temp directory."""
    out = tmp_path / "figs"
    env = dict(os.environ, MPLCONFIGDIR=str(tmp_path / "mpl"))
    r = subprocess.run(
        [sys.executable, str(SCRIPTS / "10_make_figures.py"),
         "--outdir", str(out), "--svg"],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=900, env=env)
    assert r.returncode == 0, f"figure render failed:\n{r.stdout[-3000:]}\n{r.stderr[-2000:]}"
    return out


# The FreeType that rendered the committed figures. Text extents, and with them
# the laid-out coordinates in the PDFs and the tight bounding box of the PNGs,
# depend on it: HPC4's matplotlib 3.10.0 wheel bundles FreeType 2.6.1, this conda
# build links 2.13.3, and there every PDF and PNG differs (text moved by
# fractions of a point, two PNGs a pixel narrower) while all five SVGs are
# identical after normalisation (render probe, HPC4 job 129926, 2026-09-16).
# So the SVGs are compared on every platform and the PDFs and PNGs only where
# the FreeType matches.
_FIGURES_FREETYPE = "2.13.3"


def _comparable_figure_suffixes() -> set[str]:
    import matplotlib.ft2font as _ft

    if _ft.__freetype_version__ == _FIGURES_FREETYPE:
        return {".png", ".pdf", ".svg"}
    return {".svg"}


def _figure_drift(fresh):
    """Committed figure files a fresh render does not reproduce.

    Compares only the formats `_comparable_figure_suffixes()` allows here.
    """
    committed = SCRIPTS.parent / "results" / "figures"
    suffixes = _comparable_figure_suffixes()
    drifted, compared = [], 0
    for path in sorted(committed.glob("figure*.*")):
        if path.suffix not in suffixes:
            continue
        other = fresh / path.name
        assert other.exists(), (
            f"{path.name} is committed but a fresh render did not produce it "
            "-- note SVG needs the --svg flag, which reproduce.sh does NOT pass")
        compared += 1
        if path.suffix == ".png":
            same = path.read_bytes() == other.read_bytes()
        elif path.suffix == ".pdf":
            same = _strip_pdf_date(path.read_bytes()) == _strip_pdf_date(other.read_bytes())
        else:
            same = _normalise_svg(path.read_text()) == _normalise_svg(other.read_text())
        if not same:
            drifted.append(path.name)
    # Five figures x the comparable formats. A glob that matched fewer has gone vacuous.
    assert compared >= 5 * len(suffixes), (
        f"only {compared} figure file(s) compared; the check has gone vacuous")
    return drifted


def test_committed_figures_match_a_fresh_render(tmp_path):
    """The GENERATED-ARTEFACT axis applied to the fifteen committed figures.

    Session 42 swept that axis over the supplement and `gene_set_sizes.csv` and
    found it clean, but could not cover the figures in the same test: rendering
    needs matplotlib and takes minutes, and the SVGs need a flag the pipeline's
    own `reproduce.sh` never passes.

    That last point is the live hazard and is why this test exists. `reproduce.sh`
    calls `10_make_figures.py` with NO `--svg`, so a session that re-rendered via
    the reproduce path would refresh five PDFs and five PNGs and leave five SVGs
    untouched. Nothing would have noticed a stale SVG, and the poster embeds them.

    Measured 2026-09-08: all fifteen are in sync. PNGs byte-identical; PDFs differ
    by exactly six bytes, all inside `/CreationDate`; SVGs identical once the
    per-render matplotlib identifiers are normalised.
    """
    drifted = _figure_drift(_render_figures(tmp_path))
    assert not drifted, (
        f"these committed figures are NOT what the script now renders: {drifted}. "
        "Re-render with `10_make_figures.py --outdir <tmp> --svg` and compare "
        "before committing; do not hand-edit a figure.")


def test_the_figure_render_check_can_fail(tmp_path):
    """Prove the normalisation above hides identifiers and NOTHING else.

    The risk here is not a comparison that fails to fire -- it is one that is too
    generous and reports green over a real change. So the mutations are content:
    a stroke width, a colour, a coordinate, and a PDF byte outside the date.
    """
    fresh = _render_figures(tmp_path)
    assert not _figure_drift(fresh), "baseline must be clean before mutating"

    svg = fresh / "figure2_reliability.svg"
    original = svg.read_text()
    for label, mutated in (
        ("stroke width", original.replace("stroke-width: 0.8", "stroke-width: 9.9", 1)),
        ("colour", original.replace("#0072b2", "#ff0000", 1)),
        # Re-targeted 2026-09-24 when H87 re-rendered the figures at print size
        # (the old target, x="208.420142", no longer exists); any coordinate
        # written once in the file will do.
        ("coordinate", original.replace('x="123.853408"', 'x="999.999999"', 1)),
    ):
        assert mutated != original, f"the {label} mutation matched nothing; re-target it"
        svg.write_text(mutated)
        assert "figure2_reliability.svg" in _figure_drift(fresh), (
            f"a changed {label} was NOT detected -- the SVG normalisation is too "
            "aggressive and would pass a genuinely stale figure")
    svg.write_text(original)

    # And a PDF change that is not the timestamp must still be caught -- where
    # PDFs are compared at all (see _comparable_figure_suffixes).
    if ".pdf" not in _comparable_figure_suffixes():
        return
    pdf = fresh / "figure3_outcome.pdf"
    blob = pdf.read_bytes()
    pdf.write_bytes(blob[:-40] + bytes(40))
    assert "figure3_outcome.pdf" in _figure_drift(fresh), (
        "a PDF content change was not detected; stripping /CreationDate has "
        "removed too much")


def _journal_print_sizes(tmp_path, **rc):
    """Draw the four journal figures and measure them, writing nothing.

    Runs the figure functions of `10_make_figures.py` with `save` replaced by the
    script's own `print_size`, so what is measured is exactly what a render
    measures. `rc` overrides rcParams AFTER the script sets them -- the can-fail
    below uses it to put the old 7 pt ticks back.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.environ.setdefault("MPLCONFIGDIR", str(tmp_path / "mpl"))
    mod = _load_script("10_make_figures.py")
    plt.rcParams.update(rc)
    mod.PRINT_SIZES.clear()

    def _measure(fig, outdir, stem):
        mod.print_size(fig, stem)
        plt.close(fig)

    mod.save = _measure
    try:
        for fn in (mod.figure1, mod.figure2, mod.figure3, mod.figure4):
            fn(tmp_path)
    finally:
        plt.rcdefaults()
    return mod, dict(mod.PRINT_SIZES)


def test_journal_figures_print_no_text_under_8_pt(tmp_path):
    """H87 (F11.6): at double-column width no figure prints text under 8 pt.

    Measured 2026-09-24 before the fix: the smallest text in figures 1-4 printed
    at 5.61 / 5.84 / 5.87 / 6.96 pt at 6.75 in -- legends and annotations set
    explicitly at 6-6.5 pt, below even the 7 pt ticks F11.6 had measured. Every
    drawn text item is counted, and none may sit off the canvas, where it would
    either be cropped or widen the figure past the column.
    """
    mod, sizes = _journal_print_sizes(tmp_path)
    assert set(sizes) == set(mod.JOURNAL_FIGURES), sizes.keys()
    for stem, r in sizes.items():
        assert r["n_text"] >= 20, f"{stem}: only {r['n_text']} text items measured"
        assert r["width_in"] <= mod.JOURNAL_WIDTH_IN + 1e-6, (stem, r["width_in"])
        assert not r["outside_canvas"], (stem, r["outside_canvas"])
        assert r["printed_smallest_pt"] >= mod.JOURNAL_MIN_PT - 1e-9, (
            f"{stem} prints its smallest text at {r['printed_smallest_pt']:.2f} pt")
    assert (mod.JOURNAL_WIDTH_IN, mod.JOURNAL_MIN_PT) == (6.75, 8.0)


def test_the_journal_print_size_check_can_fail(tmp_path):
    """Put the old 7 pt tick labels back and the measurement must see them."""
    _, sizes = _journal_print_sizes(tmp_path, **{"xtick.labelsize": 7,
                                                 "ytick.labelsize": 7})
    low = [s for s, r in sizes.items() if r["printed_smallest_pt"] < 8.0]
    assert len(low) == 4, f"7 pt ticks went unseen in {set(sizes) - set(low)}"


def test_the_poster_and_the_figure_script_share_one_okabe_ito_palette():
    """Two files, one colour list, hand-checked three times -- now pinned.

    The poster is printed at A0 and the figures are embedded in it, so a palette
    that drifts between them is a visible defect at print time and an
    accessibility regression that no existing check can see. Okabe-Ito is chosen
    because it is colourblind-safe; silently losing one member loses that.
    """
    script = (SCRIPTS / "10_make_figures.py").read_text()
    poster = (SCRIPTS.parents[1] / "poster" / "poster.html").read_text()

    okabe_ito = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]
    for colour in okabe_ito:
        assert colour.lower() in script.lower(), (
            f"{colour} is no longer in the figure script's palette")
        assert colour.lower() in poster.lower(), (
            f"{colour} is in the figure script but NOT in the poster; the two "
            "have drifted apart and the A0 print will not match the figures")

    # Guard the pairing itself: if someone swaps the figure script to a
    # different scheme, the loop above must not pass because the poster happens
    # to carry the old colours somewhere unrelated.
    assert "Okabe-Ito" in script, (
        "the figure script no longer names Okabe-Ito; if the scheme changed "
        "deliberately, update this test and the poster together")


# Poster figure-box geometry, measured in a browser on 2026-09-24 by evaluating
# getBoundingClientRect() over `.figbox` (1 CSS px = 25.4/96 mm; the poster is
# authored at true print size), each value ROUNDED DOWN. First measured
# 2026-09-08 as 1160x714 / 1587x687 / 1160x501 / 1595x409; the layout was
# changed on 2026-09-24 so that the two 16-row figures could print at 18 pt
# (TODO-PRINT-6, ledger H88). `10_make_figures.py`'s POSTER_BOXES_PX is the
# render's copy; the first test below holds the two equal.
_POSTER_FIGBOXES = [
    # (figure stem, box width px, box height px)
    ("figure0_schematic", 1160, 644),
    ("figure3_outcome", 1587, 616),
    ("figure2_reliability", 1160, 430),
    ("figure1_isi", 1594, 601),
]
_POSTER_MIN_PT = 18.0
_POSTER = SCRIPTS.parents[1] / "poster"


def _svg_size_px(path):
    """An SVG's intrinsic size in CSS px, from its pt width and height."""
    head = path.read_text()[:1500]
    w_pt = float(re.search(r'width="([\d.]+)pt"', head).group(1))
    h_pt = float(re.search(r'height="([\d.]+)pt"', head).group(1))
    return w_pt * 96 / 72, h_pt * 96 / 72


def _poster_smallest_text(tmp_path):
    """Draw the four poster figures in poster mode and measure, writing nothing."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.environ.setdefault("MPLCONFIGDIR", str(tmp_path / "mpl"))
    mod = _load_script("10_make_figures.py")
    mod.POSTER = True
    mod.PRINT_SIZES.clear()

    def _measure(fig, outdir, stem):
        mod.print_size(fig, stem)
        plt.close(fig)

    mod.save = _measure
    try:
        for fn in (mod.figure0, mod.figure1, mod.figure2, mod.figure3):
            fn(tmp_path)
    finally:
        plt.rcdefaults()
    return mod, dict(mod.PRINT_SIZES)


def _printed_smallest_pt(stem, box_w_px, box_h_px, smallest_pt):
    """What the poster prints: the render's smallest text times `contain`'s scale."""
    w_px, h_px = _svg_size_px(_POSTER / "figures" / f"{stem}.svg")
    return smallest_pt * min(box_w_px / w_px, box_h_px / h_px)


def test_a0_print_legibility(tmp_path):
    """Every poster figure prints its smallest text at 18 pt or more on A0.

    TODO-PRINT-6, decided 2026-09-24 (the most rigorous option). Until then the
    poster embedded the journal figures in boxes of another shape, and their
    7 pt tick labels printed at 13.5 / 12.2 / 10.5 / 7.2 pt (figure0 / 3 / 2 /
    1); `figure1_isi` sat height-constrained in a 422 x 108 mm box. Now each
    poster figure is drawn with its box's aspect at 1/2.25 its size, so the
    browser enlarges it by 2.25. Measured here from four independent pieces: the
    box table, the committed SVG's own size, the smallest text actually drawn,
    and the poster's own <img> sources -- a poster pointing back at the journal
    figures, or a render made for another box, fails.
    """
    mod, sizes = _poster_smallest_text(tmp_path)
    assert dict((s, (w, h)) for s, w, h in _POSTER_FIGBOXES) == mod.POSTER_BOXES_PX
    assert mod.POSTER_MIN_PT == _POSTER_MIN_PT

    html = (_POSTER / "poster.html").read_text()
    embedded = re.findall(r'<img src="([^"]+\.svg)"', html)
    assert sorted(embedded) == sorted(f"figures/{s}.svg" for s, _, _ in _POSTER_FIGBOXES), (
        f"the poster embeds {embedded}, not its four poster renders")

    for stem, w, h in _POSTER_FIGBOXES:
        r = sizes[stem]
        assert not r["outside_canvas"], (stem, r["outside_canvas"])
        assert not r["overlaps"], (stem, r["overlaps"])
        # Drawn for THIS box: the committed SVG has the box's aspect, so the
        # poster letterboxes it by under two pixels either way.
        w_px, h_px = _svg_size_px(_POSTER / "figures" / f"{stem}.svg")
        scale = min(w / w_px, h / h_px)
        assert w - w_px * scale < 2 and h - h_px * scale < 2, (stem, w_px, h_px)
        printed = _printed_smallest_pt(stem, w, h, r["smallest_pt"])
        assert printed >= _POSTER_MIN_PT - 1e-6, (
            f"{stem} prints its smallest text at {printed:.2f} pt on A0")
    assert len(_POSTER_FIGBOXES) == 4, "the poster carries four figure boxes"


def test_the_a0_legibility_check_can_fail():
    """Prove the check fires, in the direction the real defect takes.

    The defect is a figure that scales DOWN -- because its box got shorter than
    the one it was drawn for. Halving the box must take the printed size far
    below the floor.
    """
    stem, w, h = "figure1_isi", 1594, 601
    baseline = _printed_smallest_pt(stem, w, h, 8.0)
    assert baseline >= _POSTER_MIN_PT - 1e-6
    shrunk = _printed_smallest_pt(stem, w, h // 2, 8.0)
    assert shrunk < baseline * 0.6, (
        "halving the figure box did not materially shrink the printed size, so "
        "the arithmetic is not actually reading the box geometry")
    assert shrunk < _POSTER_MIN_PT, "the floor would not have caught a halved box"


def test_committed_poster_figures_match_a_fresh_render(tmp_path):
    """The GENERATED-ARTEFACT axis applied to the four committed poster renders."""
    out = tmp_path / "poster-figs"
    env = dict(os.environ, MPLCONFIGDIR=str(tmp_path / "mpl"))
    r = subprocess.run(
        [sys.executable, str(SCRIPTS / "10_make_figures.py"), "--poster",
         "--outdir", str(out)],
        cwd=SCRIPTS.parent, capture_output=True, text=True, timeout=900, env=env)
    assert r.returncode == 0, f"poster render failed:\n{r.stdout[-3000:]}\n{r.stderr[-2000:]}"
    fresh = sorted(p.name for p in out.glob("*.svg"))
    committed = sorted(p.name for p in (_POSTER / "figures").glob("*.svg"))
    assert fresh == committed == sorted(f"{s}.svg" for s, _, _ in _POSTER_FIGBOXES)
    drifted = [n for n in fresh
               if _normalise_svg((out / n).read_text())
               != _normalise_svg((_POSTER / "figures" / n).read_text())]
    assert not drifted, (
        f"these committed poster figures are NOT what the script now renders: "
        f"{drifted}. Re-render with `10_make_figures.py --poster`; do not hand-edit.")


# The number checker treats this shape as a LIVE claim about the current suite.
# `131-test` and `131 passed` are invisible to it; `131 tests` is not.
_LIVE_TEST_COUNT = re.compile(r'\**(\d{2,4})\**\s+(?:\w+\s+)?tests\b')


def _run_rows(text):
    """Table rows in results/README.md that classify a named .log file."""
    return [(i, line) for i, line in enumerate(text.splitlines(), 1)
            if line.startswith("| `") and ".log`" in line.split("|")[1]]


def test_historical_run_rows_state_no_live_test_count():
    """A record of the past must not need editing when the present changes.

    `results/README.md`'s run rows record what a NAMED run reported on a DATE.
    They are protected history. But the number checker cannot tell a historical
    sentence from a live one -- it only sees the shape -- so a row written as
    "**131 tests**, unchanged from _run3" reads as a claim about the CURRENT
    suite and goes red the moment the suite grows, demanding an edit to a record
    that was correct when written.

    This file's own convention is the hyphenated form (`131-test`), which the
    pattern cannot see. The convention was documented but nothing enforced it,
    and the cost is measured: session 41 introduced FOUR such phrasings in one
    new row plus four more of "an idle figure at 131 tests still does not exist";
    session 42 fixed all eight and then, writing three fresh rows an hour later,
    introduced TWO MORE of its own. Ten instances, two sessions, one convention
    everybody agreed with. That is what a guard is for.
    """
    readme = SCRIPTS.parent / "results" / "README.md"
    rows = _run_rows(readme.read_text())
    # Every timed run in this project's history has a row; a handful means the
    # parser stopped matching the table format rather than the rows disappearing.
    assert len(rows) >= 30, (
        f"only {len(rows)} run row(s) parsed out of results/README.md; the "
        "scan has gone vacuous and would pass over any number of live claims")

    # Honour the SAME exemption the number checker honours. A line carrying an
    # explicit `numcheck: ignore` is a decision already taken and recorded with a
    # reason -- re-flagging it here would be a second guard overruling the first.
    exempt = [(n, line) for n, line in rows if "numcheck: ignore" in line]
    assert exempt, (
        "no run row carries a `numcheck: ignore` marker any more, so this "
        "exemption branch is dead code -- delete it or find where it went")

    offenders = []
    for lineno, line in rows:
        if "numcheck: ignore" in line:
            continue
        for m in _LIVE_TEST_COUNT.finditer(line):
            offenders.append(f"  line {lineno}: ...{line[max(0, m.start()-55):m.end()+25]}...")
    assert not offenders, (
        "these historical run rows state a test count in the LIVE form, so they "
        "will go red every time the suite grows:\n" + "\n".join(offenders) +
        "\n\nRewrite as `NNN-test` (hyphenated) or \"the suite held NNN at the "
        "time\". Do NOT change the number -- the row records what that run reported.")


def test_the_historical_run_row_check_can_fail():
    """Prove the guard fires on the exact phrasings that were shipped twice."""
    real = [
        "| `full_suite_x.log` | **131 tests**, unchanged from `_run3`. | notes |",
        "| `full_suite_x.log` | An idle figure at 131 tests still does not exist. | n |",
        "| `full_suite_x.log` | — 136 tests against the 104-test run | n |",
    ]
    for line in real:
        rows = _run_rows(line)
        assert rows, f"the row parser did not recognise this as a run row: {line!r}"
        assert _LIVE_TEST_COUNT.search(rows[0][1]), (
            f"the guard missed a phrasing that was actually shipped: {line!r}")

    # And it must NOT fire on the correct forms, or it would forbid the fix.
    for line in [
        "| `full_suite_x.log` | a 136-test run against the 104-test run | n |",
        "| `full_suite_x.log` | **136 passed in 191.45s**, ZERO failures | n |",
        "| `full_suite_x.log` | the suite held 131 at the time | n |",
        "| `full_suite_x.log` | across 125-, 131- and 131-test runs | n |",
    ]:
        rows = _run_rows(line)
        assert not _LIVE_TEST_COUNT.search(rows[0][1]), (
            f"the guard fires on a CORRECT historical form, which would make the "
            f"documented fix impossible: {line!r}")


# Text a reader sees but no prose scan reaches. `_SUBMITTED_DOCS` covers four
# MARKDOWN/HTML documents; the five deposited figures are rendered from a python
# script, so every axis label, legend entry and schematic caption in them was
# outside the US-spelling guard entirely. It had been wrong the whole time:
# figure0's poster panel read "Residualise the SIGNATURE" while the manuscript,
# abstract, cover letter and poster prose all said "residualized". The word list
# already contained `residualis` -- only the SCOPE was missing.
#
# Data values are NOT prose and must keep their spelling: `axis_residualised` is
# a value inside the frozen, provenance-hashed outcome_arm.csv of six result
# directories, and "correcting" it would select nothing and empty a figure.
_FIGURE_TEXT_ALLOWED = {"axis_residualised"}


# Calls whose string arguments end up as pixels. `box` is the script's own
# schematic helper -- figure0 draws entirely through it, which is where the
# defect lived, so leaving it out would have made this guard miss the very
# instance that prompted it.
_DRAWS_TEXT = {
    "text", "annotate", "set_xlabel", "set_ylabel", "set_title", "suptitle",
    "set_xticklabels", "set_yticklabels", "legend", "colorbar", "set_label",
    "box",
}


def test_draws_text_covers_every_text_call_the_figure_script_makes():
    """UNCLASSIFIED-MEMBER, the reverse direction (F10.9).

    `_DRAWS_TEXT` decides what the spelling guard can see. Nothing checked that
    it still names every text-drawing call the script makes, so a new
    `ax.annotate(...)` — or a `set_xticklabels` in a figure added later — would
    silently fall outside the guard. This walks the script's AST for calls whose
    name is text-drawing by SHAPE (ends in `label` or `title`, or is one of
    matplotlib's text entry points) and asserts each is classified.
    """
    import ast

    shape = {"text", "annotate", "legend", "colorbar", "figtext", "bar_label", "table"}
    used = set()
    for node in ast.walk(ast.parse((SCRIPTS / "10_make_figures.py").read_text())):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            name = node.func.attr
            # A getter draws nothing: `get_legend_handles_labels` (2026-09-24,
            # the poster keys) ends in "labels" by shape only.
            if name.startswith("get_"):
                continue
            if name in shape or name.endswith(("label", "labels", "title")):
                used.add(name)
    assert used, "the AST walk found no text call at all; the check has gone vacuous"
    assert used <= _DRAWS_TEXT, (
        f"these text-drawing calls are outside _DRAWS_TEXT, so the figure spelling "
        f"guard cannot see their strings: {sorted(used - _DRAWS_TEXT)}")


def _figure_drawn_strings():
    """String literals reaching a drawing call, as (line, text).

    Deliberately NOT every literal in the file: `print()` messages and column
    names are not text a reader sees, and folding them in here would make a
    guard named for rendered text fire on things that never render.
    """
    import ast
    tree = ast.parse((SCRIPTS / "10_make_figures.py").read_text())
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = (fn.attr if isinstance(fn, ast.Attribute)
                else fn.id if isinstance(fn, ast.Name) else None)
        if name not in _DRAWS_TEXT:
            continue
        for sub in ast.walk(node):
            if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
                out.append((sub.lineno, sub.value))
    return out


def test_rendered_figure_text_uses_american_spelling():
    """The figures are deposited and printed on the poster; they are prose too."""
    lits = _figure_drawn_strings()
    # A parse that returned a handful of strings would pass over anything.
    assert len(lits) >= 40, (
        f"only {len(lits)} drawn string(s) parsed out of 10_make_figures.py; "
        "the scan has gone vacuous -- most likely a drawing call was renamed "
        "and _DRAWS_TEXT no longer names it")

    offenders = []
    for lineno, s in lits:
        if s in _FIGURE_TEXT_ALLOWED:
            continue
        low = s.lower()
        for brit in _BRITISH:
            if brit in low:
                offenders.append(f"  line {lineno}: {brit!r} in {s[:70]!r}")
    assert not offenders, (
        "these strings in 10_make_figures.py use British spelling and are "
        "RENDERED INTO THE DEPOSITED FIGURES:\n" + "\n".join(offenders) +
        "\n\nFix the generator and re-render (no --outdir; pass --svg, because "
        "the poster embeds the SVGs). If the string is a DATA VALUE rather than "
        "prose, add it to _FIGURE_TEXT_ALLOWED with the reason.")


def test_the_figure_spelling_check_can_fail():
    """Mutate toward the real defect, and keep a positive control.

    The first half replays the exact string that shipped. The second proves the
    check does not simply flag everything -- if it did, it would also condemn
    the data value the figures must keep, and the 'fix' would break Figure 3.
    """
    scan = lambda s: [b for b in _BRITISH if b in s.lower()]          # noqa: E731
    assert scan("1.  Residualise the SIGNATURE"), (
        "the guard misses the phrasing that was actually shipped in figure0")
    assert scan("Axis-residualised $r$ with image prediction"), (
        "the guard misses the axis label that was actually shipped")
    for ok in ("Axis-residualized $r$ with image prediction",
               "1.  Residualize the SIGNATURE",
               "Raw score", "excess Fisher z"):
        assert not scan(ok), f"the guard fires on correct US prose: {ok!r}"
    # The data value trips the word list on purpose -- the allowlist, not the
    # pattern, is what keeps it out of the offenders list.
    assert scan("axis_residualised"), (
        "if this stops matching, _FIGURE_TEXT_ALLOWED is dead weight and the "
        "comment explaining why it exists has become misleading")
    assert "axis_residualised" in _FIGURE_TEXT_ALLOWED


# ---------------------------------------------------------------------------
# Reading the ARTEFACT instead of the generator.
#
# `test_rendered_figure_text_uses_american_spelling` reads the literal strings
# out of `10_make_figures.py`. That is the source, not the artefact, and session
# 44's handoff named the consequence exactly: "a hand-edited figure would pass
# it, and figures are not provenance-hashed, so nothing else would notice
# either."
#
# The committed `.svg`s cannot close that gap -- matplotlib's `svg.fonttype`
# defaults to `path`, so they contain zero `<text>` elements and no tool can
# read a word out of one. The `.pdf`s CAN, because the Type-42 font pin session
# 43 added for an unrelated reason (Type 3 fonts are unacceptable to journals)
# embeds real fonts with real character codes. Decompress the content streams
# and read the operands of the text-showing operators.
#
# matplotlib splits a single text call across several operands at kerning pairs
# -- "residualized" arrives as "r" + "esidualized" -- so the strings must be
# JOINED before searching. Measured on the five committed PDFs: 679 drawn
# strings, and "residualized" is findable in three of them only after joining.
# ---------------------------------------------------------------------------


def _pdf_drawn_text(path: Path) -> list[str]:
    """Strings actually drawn into a PDF, in content-stream order."""
    import zlib
    raw = path.read_bytes()
    out: list[str] = []
    for chunk in re.findall(rb"stream\r?\n(.*?)endstream", raw, re.S):
        try:
            data = zlib.decompress(chunk)
        except Exception:
            continue
        for m in re.finditer(rb"\(((?:[^()\\]|\\.)*)\)", data):
            b = (m.group(1).replace(rb"\(", b"(")
                 .replace(rb"\)", b")").replace(rb"\\", b"\\"))
            # Type-42 subsets encode two bytes per glyph; the low byte carries
            # the character. Keep printable ASCII and drop positioning noise.
            s = "".join(chr(c) for c in b[1::2] if 32 <= c < 127)
            if s.strip():
                out.append(s)
    return out


def test_committed_figure_pdfs_are_free_of_british_spelling():
    """The deposited FIGURES, read as artefacts rather than as source.

    This is the independent companion for
    `test_rendered_figure_text_uses_american_spelling`. That test and the
    figures share one origin -- `10_make_figures.py` -- so a figure edited by
    hand, or rendered from a generator state that no longer exists, passes it
    silently. Nothing else would notice: the figures are deliberately NOT in
    `PROVENANCE.json`, which is what makes re-rendering them safe.
    """
    pdfs = sorted((SCRIPTS.parent / "results" / "figures").glob("*.pdf"))
    assert len(pdfs) >= 5, (
        f"only {len(pdfs)} committed figure PDF(s) found; this sweep has gone "
        "vacuous and would pass over a deposit that had lost figures")

    per_file = {p.name: _pdf_drawn_text(p) for p in pdfs}
    total = sum(len(v) for v in per_file.values())
    assert total >= 400, (
        f"only {total} drawn string(s) extracted from {len(pdfs)} PDFs; the "
        "extractor has stopped seeing text -- most likely the font type or the "
        "stream compression changed -- and a zero here would be meaningless")

    # POSITIVE CONTROL: the extractor must be able to see a real word, joined
    # across the kerning splits. If this fails, a clean result proves nothing.
    corpus = {name: "".join(strings).lower() for name, strings in per_file.items()}
    assert any("residualized" in c for c in corpus.values()), (
        "the extractor found no occurrence of 'residualized' in any committed "
        "figure, so it cannot be shown to read words at all and a British-"
        "spelling sweep over it would be vacuous")

    offenders = []
    for name, joined in sorted(corpus.items()):
        for brit in _BRITISH:
            if brit in joined and not any(
                    allowed.lower().replace("_", "") in joined
                    for allowed in _FIGURE_TEXT_ALLOWED if brit in allowed.lower()):
                offenders.append(f"  {name}: {brit!r} appears in the RENDERED text")
    assert not offenders, (
        "these COMMITTED FIGURES render British spelling, which the "
        "source-reading guard cannot see:\n" + "\n".join(offenders) +
        "\n\nRe-render from the generator (no --outdir; pass --svg, because the "
        "poster embeds the SVGs) rather than editing a figure by hand.")


def test_the_rendered_figure_spelling_check_can_fail():
    """Prove it fires on the artefact, in the direction the defect would take.

    The mutation is applied to the EXTRACTED strings in memory. Nothing is
    written to the repository -- session 45, 46 and 47's pattern.
    """
    pdfs = sorted((SCRIPTS.parent / "results" / "figures").glob("*.pdf"))
    assert pdfs, "no committed figure PDFs"
    strings = _pdf_drawn_text(pdfs[0])
    assert strings, f"no text extracted from {pdfs[0].name}"

    clean = "".join(strings).lower()
    hits = [b for b in _BRITISH if b in clean]
    assert not hits, (
        f"{pdfs[0].name} already renders British spelling {hits}; the control "
        "half of this test is meaningless until that is fixed")

    # The exact defect session 44 found printed on poster panel 1 and in three
    # deposited figures, reproduced here against the real extraction path.
    mutated = "".join(strings + ["Axis-r", "esidualised"]).lower()
    caught = [b for b in _BRITISH if b in mutated]
    assert "residualis" in caught, (
        "injecting 'residualised' across a kerning split did not register, so "
        "this check cannot detect the defect it exists for -- most likely the "
        "strings are being searched individually rather than joined")


# The same defect as `_LIVE_TEST_COUNT`, one axis over. A test count goes stale
# when the suite grows; a DURATION TALLY goes stale when the calendar moves, and
# nothing has to happen at all. Session 42 wrote "It has been blocked on
# available memory for EIGHT consecutive sessions" into the decline row; by
# session 44 it was ten, and the row had to be edited though nothing about that
# run had changed. The correct form is past-tense and dated -- "at that point it
# had been blocked for eight consecutive sessions" -- which is true forever.
_LIVE_DURATION_TALLY = re.compile(
    r'\b(?:has|have)\s+(?:now\s+)?been\b[^.|]{0,80}?'
    r'\b(?:for|across|over)\s+'
    r'(?:\d{1,3}|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|'
    r'thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|nineteen|twenty)\s+'
    r'(?:more\s+|further\s+|consecutive\s+)?'
    r'(?:sessions?|handoffs?|days?|weeks?)\b',
    re.I)


def test_historical_run_rows_state_no_live_duration_tally():
    """A record of the past must not need editing when the CALENDAR changes.

    `test_historical_run_rows_state_no_live_test_count` catches a count that
    goes stale when the suite grows. This catches the shape that goes stale
    when nothing happens at all: a running tally of elapsed sessions, days or
    handoffs, written in the present perfect inside a row that records a run
    that finished on a fixed date.
    """
    readme = SCRIPTS.parent / "results" / "README.md"
    rows = _run_rows(readme.read_text())
    # Same floor as its sibling: a handful of rows means the parser stopped
    # matching the table, not that the rows went away.
    assert len(rows) >= 30, (
        f"only {len(rows)} run row(s) parsed out of results/README.md; the "
        "scan has gone vacuous and would pass over any number of live tallies")

    offenders = []
    for lineno, line in rows:
        if "numcheck: ignore" in line:
            continue
        for m in _LIVE_DURATION_TALLY.finditer(line):
            offenders.append(
                f"  line {lineno}: ...{line[max(0, m.start() - 40):m.end() + 20]}...")
    assert not offenders, (
        "these historical run rows carry a duration tally in the LIVE form, so "
        "they go stale with the calendar rather than with any event:\n"
        + "\n".join(offenders) +
        "\n\nRewrite in the past tense against a date -- \"at that point it had "
        "been blocked for eight consecutive sessions\". Do NOT just bump the "
        "number; that is the edit this guard exists to make unnecessary.")


def test_the_duration_tally_check_can_fail():
    """Mutate in the direction the real defect took, both ways.

    A detector that fired on everything would pass the first half of this
    test while being worthless, so the second half is the positive control:
    the corrected phrasings must NOT trip it, or the documented fix would be
    forbidden and the guard would force the defect back in.
    """
    real = [
        "| `pancancer_v3_stablesort.log` | It has been blocked on available "
        "memory for EIGHT consecutive sessions; the floor was 3.59 GiB. | n |",
        "| `full_suite_x.log` | That advice has now been carried across five "
        "handoffs unexecuted. | notes |",
        "| `full_suite_x.log` | The sort audit has been deferred for 13 "
        "sessions. | notes |",
    ]
    for line in real:
        rows = _run_rows(line)
        assert rows, f"the row parser did not recognise this as a run row: {line!r}"
        assert _LIVE_DURATION_TALLY.search(rows[0][1]), (
            f"the guard missed a phrasing that was actually shipped: {line!r}")

    for line in [
        "| `pancancer_v3_stablesort.log` | At that point it had been blocked "
        "for eight consecutive sessions; the floor was 3.59 GiB. | n |",
        "| `full_suite_x.log` | Blocked on memory for eight sessions when "
        "this was decided on 2026-09-08. | notes |",
        "| `full_suite_x.log` | Twenty samples on 2026-09-09 gave a floor of "
        "7.69 GiB. | notes |",
        "| `full_suite_x.log` | **146 passed in 209.21s**, ZERO failures | n |",
    ]:
        rows = _run_rows(line)
        assert not _LIVE_DURATION_TALLY.search(rows[0][1]), (
            "the guard fires on a CORRECT dated form, which would make the "
            f"documented fix impossible: {line!r}")


def _strip_code_and_comments(text: str) -> str:
    """Remove HTML comments, fenced blocks and inline code before scanning prose."""
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    text = re.sub(r"^\s*```.*?^\s*```", "", text, flags=re.S | re.M)
    return re.sub(r"`[^`]*`", "", text)


def _odd_bold_paragraphs(text: str):
    """[(paragraph index, first line)] where `**` markers do not pair up."""
    out = []
    for i, para in enumerate(_strip_code_and_comments(text).split("\n\n")):
        if para.count("**") % 2:
            first = next((ln for ln in para.strip().splitlines() if ln.strip()), "")
            out.append((i, first[:110]))
    return out


def test_prose_documents_have_balanced_bold_markers():
    """The half-edit signature the bracket scan is blind to.

    `test_prose_documents_have_balanced_delimiters` catches an orphaned `(` or
    `[`, which is how session 40's defect showed up. It cannot see an orphaned
    `**`, and session 42 identified that as the likeliest REMAINING half-edit
    signature -- then found one on the first run.

    The real instance, in `submission/SUBMISSION-CHECKLIST.md`: the text read
    `**$4,200 (AACR member) / $5,000 (non-member)*` -- bold opened with `**` and
    closed with a single `*` -- while the parallel clause two lines later
    correctly read `**$2,550, reduced to $2,200**`. It renders wrong, and no
    number was involved, so nothing else in this repository could have seen it.

    Scoped per PARAGRAPH rather than per line, because bold legitimately spans a
    line wrap but never a blank line. Noise measured before shipping: exactly one
    hit across 29 documents, and it was the genuine defect.
    """
    repo = Path(__file__).resolve().parents[2]
    docs = _prose_documents(repo)
    assert len(docs) >= 20, (
        f"only {len(docs)} prose document(s) discovered; the scan is near-vacuous")

    offenders = []
    for d in docs:
        for _, first in _odd_bold_paragraphs(d.read_text(errors="replace")):
            offenders.append(f"  {d.relative_to(repo)}: {first!r}")
    assert not offenders, (
        "these paragraphs open a `**` they never close, so they render wrong:\n"
        + "\n".join(offenders))


def test_the_bold_balance_check_can_fail():
    """It must catch the real shape and stay quiet on correct markdown."""
    real = "Fees are **$4,200 (member) / $5,000\n(non-member)* and rising.\n"
    assert _odd_bold_paragraphs(real), "the guard missed the exact defect that shipped"

    for ok in (
        "A **bold** phrase and an *italic* one.\n",
        "Bold that **spans a\nline wrap** is fine.\n",
        "Two paragraphs.\n\nEach with **its own** bold.\n",
        "Code is exempt: `a ** b` and `x **kwargs`.\n",
        "<!-- a comment with ** in it -->\n",
    ):
        assert not _odd_bold_paragraphs(ok), f"false positive on correct markdown: {ok!r}"

    # A fenced block full of Python power operators must not fire.
    fenced = "text\n\n```python\nx = a ** b\ny = c ** d\n```\n\nmore text\n"
    assert not _odd_bold_paragraphs(fenced), "fenced code is not being stripped"


_BLOCK_START = re.compile(r"^(\s*$|[-*+] |\d+\. |#|\||>|<|```|\*\*|---|===)")


def _flush_left_continuations(text: str):
    """[(line number, text)] where flush-left prose follows an indented list line.

    Markdown allows it ("lazy continuation"), and it renders as part of the
    bullet above, so nothing looks wrong. Every list item in these documents
    indents its continuation lines, so a flush-left line straight after one is
    the scar of an insertion that split a sentence and carried its tail away.
    """
    lines = text.splitlines()
    out, fenced = [], False
    for i in range(1, len(lines)):
        if lines[i - 1].lstrip().startswith("```"):
            fenced = not fenced
        if fenced:
            continue
        prev, cur = lines[i - 1], lines[i]
        if (re.match(r"^ {2,}\S", prev) and re.match(r"^[A-Za-z(`\"']", cur)
                and not _BLOCK_START.match(cur)):
            out.append((i + 1, cur[:100]))
    return out


def test_no_list_item_carries_a_stranded_flush_left_line():
    """An insertion that splits a sentence strands its tail (session 51).

    The real instance, in `14-SCIENCE-AUDIT.md`: the F6.7 section ended "...it
    is now declared. The / report stops at `src/`...", and session 50's next
    five appends each landed between "The" and "report", so the sentence ended
    up flush-left at the bottom of an unrelated bullet (F7.4), reading as part
    of it. No number was wrong and every checker was green. A second instance
    of the class, in `15-CHECKLIST.md` the same week, had a different shape (a
    heading inserted above a paragraph's tail) that this scan cannot see; the
    class is named in the handoff, not closed by this test.

    Noise measured before shipping: zero hits across the scanned documents
    once the real instance was repaired, one (the real instance) before.
    """
    repo = Path(__file__).resolve().parents[2]
    docs = [d for d in _prose_documents(repo) if d.suffix == ".md"]
    assert len(docs) >= 20, (
        f"only {len(docs)} markdown document(s) discovered; the scan is near-vacuous")
    offenders = [f"  {d.relative_to(repo)}:{n}: {line!r}"
                 for d in docs
                 for n, line in _flush_left_continuations(d.read_text(errors="replace"))]
    assert not offenders, (
        "flush-left text directly under an indented list line -- usually a "
        "sentence tail an insertion carried away from its paragraph:\n"
        + "\n".join(offenders))


def test_the_stranded_line_scan_can_fail():
    """It must catch the real shape and stay quiet on ordinary markdown."""
    real = ("- **F7.4.** Seventeen statements now have authorities. The checker\n"
            "  reads 872 claims after this. The\n"
            "report stops at `src/`: scripts changed far more.\n")
    assert _flush_left_continuations(real), "the guard missed the shape that shipped"
    for ok in (
        "- a bullet\n  with an indented continuation\n\nA new paragraph.\n",
        "- one\n  two\n- three\n",
        "1. first\n   continued\n2. second\n",
        "- item\n  more\n## Heading\n",
        "- item\n  more\n| a | b |\n",
        "  ```\n  code\n```\nplain text after a fence\n",
    ):
        assert not _flush_left_continuations(ok), f"false positive: {ok!r}"


def test_abstract_counter_scores_the_nextgen_rule_with_the_opposite_convention():
    """`--nextgen` inverts the space rule, which is exactly what an edit gets wrong.

    The regular AACR abstract excludes spaces from its 2,600-character limit; the
    NextGen Stars extended abstract INCLUDES them against 8,000. Session 42 tested
    the regular path in three directions and left this one untested, noting that
    the asymmetry is the kind of thing a refactor silently unifies.
    """
    m = _load_script("08_count_abstract.py")
    assert m.NEXTGEN_LIMIT == 8000
    assert m.LIMIT == 2600, "the two limits must not be unified"

    draft = SCRIPTS.parents[1] / "13-NEXTGEN-EXTENDED-ABSTRACT.md"
    assert draft.is_file(), "the NextGen draft is gone; this test needs re-scoping"
    title, body = m.extract(draft.read_text())
    assert title and body, "the NextGen draft lost its '## TITLE' / '## BODY' sections"

    # The two conventions must give DIFFERENT answers on the same text, or one
    # of them has been changed to match the other.
    with_spaces = len(title) + len(body)
    without_spaces = m.count(title) + m.count(body)
    assert with_spaces > without_spaces, (
        "the NextGen count is not including spaces; the two rules have been "
        "collapsed into one and the 8,000 limit is now being scored wrongly")
    assert with_spaces <= m.NEXTGEN_LIMIT, (
        f"the NextGen abstract is over its limit: {with_spaces} / {m.NEXTGEN_LIMIT}")


def _pdf_font_defects(figures_dir):
    """PDFs that carry Type 3 glyphs or embed no font program at all.

    Returns (defects, examined) so a caller can floor the denominator: a scan
    that examined nothing must never read as a pass.
    """
    defects, examined = [], 0
    for pdf in sorted(Path(figures_dir).glob("figure*.pdf")):
        raw = pdf.read_bytes()
        examined += 1
        if re.search(rb"/Subtype\s*/Type3", raw):
            defects.append(f"{pdf.name}: /Subtype /Type3")
        elif not re.search(rb"/FontFile[23]?\b", raw):
            defects.append(f"{pdf.name}: no /FontFile -- no embedded font program")
    return defects, examined


def test_committed_pdf_figures_embed_a_real_font_and_are_not_type3():
    """Matplotlib's default `pdf.fonttype` is 3, and Type 3 is the format
    journal production systems most commonly reject at preflight.

    WHY THIS IS NOT COVERED BY THE FRESH-RENDER TEST, which is the whole reason
    it needs to exist separately. `test_committed_figures_match_a_fresh_render`
    compares the committed files against a fresh render of the SAME script; if
    someone drops the `pdf.fonttype` line from `10_make_figures.py`, BOTH sides
    become Type 3 and that test stays green while every figure regresses. A
    sync check can only see disagreement, never a defect both sides share.

    Measured 2026-09-08, before the fix: all five committed PDFs carried
    `/Subtype /Type3` and NO `/FontFile` of any kind. Type 3 glyphs are drawn as
    PDF content streams rather than backed by an embedded font program, so the
    figures rendered correctly everywhere -- which is exactly why nine sessions
    of figure work never noticed. `ps.fonttype` shares the default and is the
    sharper exposure: EPS is on AACR's accepted-source-format list for revised
    manuscripts and PDF is not.
    """
    figures = SCRIPTS.parent / "results" / "figures"
    defects, examined = _pdf_font_defects(figures)
    assert examined >= 5, (
        f"only {examined} committed PDF(s) examined; this scan has gone vacuous "
        "and would report a pass over an empty set")
    assert not defects, (
        "committed figure PDFs have a font defect:\n  " + "\n  ".join(defects)
        + "\nSet `pdf.fonttype`/`ps.fonttype` to 42 in 10_make_figures.py and "
          "re-render -- do NOT hand-edit the PDFs.")


def test_the_pdf_font_check_can_fail(tmp_path):
    """Mutate in the direction the real defect takes: the generator, not a byte.

    Perturbing a committed PDF would only prove that a regex works on bytes.
    The defect this guards against is someone deleting an rcParam, so the
    mutation deletes the rcParam, re-renders, and asserts the scan goes red.
    """
    src = (SCRIPTS / "10_make_figures.py").read_text()
    assert '"pdf.fonttype": 42' in src, (
        "the rcParam this test mutates is gone from 10_make_figures.py; either "
        "the fix was reverted or it moved, and this test is now vacuous")

    # The generator cannot simply be copied out of the tree and re-run: it
    # resolves both `aacr27` and `results/` relative to its own __file__, so a
    # copy in a temp directory fails on missing frozen inputs rather than on the
    # rcParam, and writing a mutated copy back INTO scripts/ would be a test
    # writing into the repository. So the mutation is applied where it actually
    # bites -- matplotlib's own PDF writer -- and the detector is run against
    # genuine Type 3 output rather than against a corrupted byte.
    out = tmp_path / "figs"
    out.mkdir()
    script = tmp_path / "emit.py"
    script.write_text(
        "import matplotlib\n"
        "matplotlib.use('Agg')\n"
        "import matplotlib.pyplot as plt\n"
        "import sys\n"
        "plt.rcParams['pdf.fonttype'] = int(sys.argv[2])\n"
        "fig, ax = plt.subplots()\n"
        "ax.plot([0, 1], [0, 1])\n"
        "ax.set_xlabel('a label with real glyphs')\n"
        "fig.savefig(sys.argv[1])\n")
    env = dict(os.environ, MPLCONFIGDIR=str(tmp_path / "mpl"))

    for fonttype, expect_defect in ((3, True), (42, False)):
        target = out / f"figure_fonttype{fonttype}.pdf"
        r = subprocess.run(
            [sys.executable, str(script), str(target), str(fonttype)],
            capture_output=True, text=True, timeout=300, env=env)
        assert r.returncode == 0, f"emitting fonttype {fonttype} failed:\n{r.stderr[-2000:]}"
        defects, examined = _pdf_font_defects(out)
        assert examined >= 1, "the detector examined nothing"
        hit = any(target.name in d for d in defects)
        if expect_defect:
            assert hit, (
                "a genuine matplotlib Type 3 PDF was NOT flagged -- the detector "
                "cannot see the defect it exists for")
            assert any("Type3" in d for d in defects), f"wrong defect kind: {defects}"
            target.unlink()
        else:
            # The positive control that matters: a detector that flags
            # everything would have passed the assertion above while being
            # useless. Type 42 output must come back clean.
            assert not hit, (
                f"a Type 42 PDF was flagged as defective ({defects}) -- the "
                "detector is red regardless of input and proves nothing")


# ---------------------------------------------------------------------------
# Session 49: the partition-sweep merge (27) and the shard's failure kinds (12)
# ---------------------------------------------------------------------------


def _write_shard(root: Path, seed: int, *, isi=None, kind=None):
    d = root / f"seed_{seed}"
    d.mkdir(parents=True, exist_ok=True)
    if isi is not None:
        pd.DataFrame([{"split_seed": seed, "isi": isi, "lo": isi - 0.02,
                       "hi": isi + 0.02, "n_signatures": 16, "runtime_s": 1.0}]
                     ).to_csv(d / "per_partition.csv", index=False)
    if kind is not None:
        pd.DataFrame([{"split_seed": seed, "kind": kind, "exception": "X",
                       "reason": "r"}]).to_csv(d / "refused_partitions.csv", index=False)


def test_partition_sweep_merge_pools_once_and_records_what_it_could_not(tmp_path):
    """Shards from two passes become ONE table and ONE spread.

    A seed present in both passes must agree bitwise (that is the evidence the
    passes compute the same ISI); a seed that failed in one pass and succeeded
    in the other counts as computed; a seed with no record at all is reported
    with the kind the caller names, because a timed-out shard cannot say so.
    """
    from scipy import stats as _st

    m = _load_script("27_partition_sweep.py")
    a, b, out = tmp_path / "a", tmp_path / "b", tmp_path / "out"
    _write_shard(a, 0, isi=0.30)
    _write_shard(a, 1, isi=0.31)
    _write_shard(a, 2, isi=0.29)
    _write_shard(a, 3, kind="linalg_error")
    _write_shard(b, 2, isi=0.29)          # duplicate, bitwise equal
    _write_shard(b, 3, isi=0.32)          # computed in the second pass
    _write_shard(b, 5, kind="linalg_error")
    rc = m.main(["merge", "--cohort", "pancancer", "--source", f"default={a}",
                 "--source", f"isi_only={b}", "--expect", "0-5",
                 "--not-computed-kind", "4=timeout", "--out", str(out)])
    assert rc == 0
    rows = pd.read_csv(out / "per_partition.csv", float_precision="round_trip")
    assert list(rows["split_seed"]) == [0, 1, 2, 3]
    assert list(rows["mode"]) == ["default", "default", "default", "isi_only"]
    nc = pd.read_csv(out / "not_computed.csv")
    assert list(zip(nc["split_seed"], nc["kind"])) == [(4, "timeout"), (5, "linalg_error")]

    import json as _json
    s = _json.loads((out / "summary.json").read_text())
    isi = np.array([0.30, 0.31, 0.29, 0.32])
    assert s["isi_sd_across_partitions"] == float(np.std(isi, ddof=1))
    lo = s["isi_sd_across_partitions"] * np.sqrt(3 / _st.chi2.ppf(0.975, 3))
    assert s["isi_sd_chi2_95ci"][0] == pytest.approx(lo, rel=1e-12)
    assert s["duplicate_seeds"] == [{"split_seed": 2, "modes": ["default", "isi_only"],
                                     "isi": ["0.29", "0.29"], "bitwise_equal": True}]
    assert s["rows_by_mode"] == {"default": 3, "isi_only": 1}
    # A deposited summary carries no path from the machine that merged it.
    assert str(tmp_path) not in (out / "summary.json").read_text()
    assert s["sources"]["default"] == {"read_from": "<outside the repository>/a", "origin": None}
    with pytest.raises(SystemExit, match="--origin names no --source"):
        m.main(["merge", "--cohort", "pancancer", "--source", f"default={a}",
                "--origin", "nosuch=x", "--expect", "0-5", "--out", str(tmp_path / "o3")])

    # It must refuse to overwrite, and refuse to pool a one-ULP disagreement.
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        m.main(["merge", "--cohort", "pancancer", "--source", f"default={a}",
                "--expect", "0-5", "--out", str(out)])
    _write_shard(b, 2, isi=float(np.nextafter(0.29, 1.0)))
    with pytest.raises(SystemExit, match="split_seed=2 differs"):
        m.main(["merge", "--cohort", "pancancer", "--source", f"default={a}",
                "--source", f"isi_only={b}", "--expect", "0-5",
                "--out", str(tmp_path / "out2")])


def test_partition_sweep_reads_every_float_exactly(tmp_path):
    """pandas' default CSV parser is not round-trip exact; the merge must be.

    Measured 2026-09-16: 16 of 44 cells of the committed S6 table sit one ULP
    from the per_partition.csv they were read from. A derived spread computed
    from mis-read inputs disagrees with its source in the last digit.
    """
    m = _load_script("27_partition_sweep.py")
    rng = np.random.default_rng(49)
    vals = rng.uniform(0.2, 0.4, 400)
    src = tmp_path / "s"
    (src).mkdir()
    pd.DataFrame({"split_seed": range(400), "isi": vals, "lo": vals - 0.05,
                  "hi": vals + 0.05, "n_signatures": 16, "runtime_s": 1.0}
                 ).to_csv(src / "per_partition.csv", index=False)
    rows, _ = m._read_source("x", src)
    exact = np.array([float(line.split(",")[1]) for line in
                      (src / "per_partition.csv").read_text().splitlines()[1:]])
    assert np.array_equal(exact, vals), "the fixture itself did not round-trip"
    default = pd.read_csv(src / "per_partition.csv")["isi"].to_numpy()
    assert np.array_equal(rows["isi"].to_numpy(), vals), (
        "the merge mis-read a float; pandas' default parser disagreed on "
        f"{int((default != vals).sum())} of 400 values on this platform")


def test_partition_variance_shard_names_each_failure_kind(tmp_path, monkeypatch):
    """A shard that cannot compute a partition records WHAT failed, and exits 0.

    Session 48 found every pan-cancer "refusal" was a LAPACK failure caught by
    an `except ValueError` and reported as geometry. The kinds are the fix; this
    is the test session 48 verified by hand only.
    """
    m = _load_script("12_partition_variance.py")
    n = 40
    cohort = pd.DataFrame({"patient_id": [f"P{i}" for i in range(n)],
                           "tss": [f"S{i % 8}" for i in range(n)],
                           "cancer_type": "LUAD", "SIG_A": np.arange(n, dtype=float)})
    monkeypatch.setattr(m, "load_nsclc", lambda: (cohort, np.zeros((n, 2)), None, None))
    monkeypatch.setattr(m.globalaxis, "compute_global_axis", lambda *a, **k: None)

    class _Est:
        value, lo, hi, n = 0.3, 0.25, 0.35, 1

    def fake_run_audit(cohort, X, sig_cols, cfg, **kw):
        if cfg.split_seed == 11:
            raise np.linalg.LinAlgError("SVD did not converge")
        if cfg.split_seed == 12:
            raise RuntimeError("degenerate site-disjoint partition: fold 3 is empty")
        if cfg.split_seed == 13:
            raise ValueError("some other value problem")
        return types.SimpleNamespace(primary_excess=_Est())

    monkeypatch.setattr(m.experiment, "run_audit", fake_run_audit)
    out = tmp_path / "shard"
    monkeypatch.setattr(sys, "argv", ["12_partition_variance.py", "--cohort", "nsclc",
                                      "--partitions", "4", "--start-seed", "10",
                                      "--shard", "--outdir", str(out)])
    assert m.main() == 0, "a shard whose every seed is accounted for must exit 0"
    rows = pd.read_csv(out / "per_partition.csv")
    assert list(rows["split_seed"]) == [10]
    refused = pd.read_csv(out / "refused_partitions.csv")
    assert list(zip(refused["split_seed"], refused["kind"])) == [
        (11, "linalg_error"), (12, "degeneracy_guard"), (13, "other_ValueError")]
    # The can-fail: the kind must come from the exception, not from the message.
    assert m._failure_kind(np.linalg.LinAlgError("degenerate site-disjoint partition")) \
        == "linalg_error"
    assert m._failure_kind(RuntimeError("SVD did not converge")) == "other_RuntimeError"


def _synthetic_raw(root: Path) -> None:
    """A raw tree with every quirk 00_build_interim.py must reproduce (A13, A16)."""
    import gzip as _gz

    import pyarrow as _pa
    import pyarrow.parquet as _pq

    rng = np.random.default_rng(13)

    def emb():
        return [list(map(float, rng.normal(size=4))) for _ in range(14)]

    rows = [  # filename, sample, type
        ("TCGA-05-0001-01Z-00-DX1.a.h5", "TCGA-05-0001-01", "LUAD"),
        ("TCGA-05-0001-01Z-00-DX1.a.h5", "TCGA-05-0001-11", "LUAD"),  # join artefact
        ("TCGA-05-0002-01Z-00-DX1.b.h5", None, "LUSC"),               # null sample
        ("TCGA-AA-0003-01Z-00-DX1.c.h5", "TCGA-AA-0003-01", "COAD"),
        ("TCGA-05-0004-01Z-00-DX1.d.h5", "TCGA-05-0004-06", "LUAD"),  # not primary
        ("TCGA-05-0005-01Z-00-DX1.e.h5", "TCGA-05-0005-01", "LUAD"),  # no expression
    ]
    table = _pa.table({
        "filename": [r[0] for r in rows],
        "embedding": [emb() for _ in rows],
        "_PATIENT": [r[0][:12] for r in rows],
        "sample": [r[1] for r in rows],
        "cancer type abbreviation": [r[2] for r in rows],
        "OS": [1, 1, 0, 1, 0, 1],
    })
    (root / "provgigapath").mkdir(parents=True)
    _pq.write_table(table, root / "provgigapath" / "embeddings.parquet")

    (root / "expression").mkdir()
    samples = ["TCGA-05-0001-11", "TCGA-05-0001-01", "TCGA-05-0002-01", "TCGA-AA-0003-01"]
    # ENSG0 duplicates G1's symbol (dropped: first wins); ENSG99 has no symbol.
    genes = ["ENSG1.1", "ENSG0.3"] + [f"ENSG{i}.2" for i in range(2, 11)] + ["ENSG99.1"]
    vals = rng.normal(size=(len(genes), len(samples))).round(4)
    lines = ["sample\t" + "\t".join(samples)]
    lines += [g + "\t" + "\t".join(f"{v:.4f}" for v in row) for g, row in zip(genes, vals)]
    with _gz.open(root / "expression" / "tcga_RSEM_gene_tpm.gz", "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    (root / "expression" / "ensembl_to_hugo.csv").write_text(
        "ensembl_gene_id,symbol\nENSG0,G1\n" + "".join(f"ENSG{i},G{i}\n" for i in range(1, 11)))
    (root / "ancestry").mkdir()
    (root / "ancestry" / "UCSF_Ancestry_Calls.csv").write_text(
        '"","Patient_ID","pam.ancestry.cluster"\n"1","TCGA-05-0001","EUR"\n')
    (root / "purity").mkdir()
    (root / "purity" / "TCGA_ABSOLUTE_purity.csv").write_text(
        "array,sample,purity\nTCGA-05-0002-01,x,0.5\n")
    (root / "signatures").mkdir()
    (root / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt").write_text(
        "HALLMARK_X\turl\tG1\tG2\tG3\tG4\tG5\n"
        "HALLMARK_Y\turl\tG6\tG7\tG8\tG9\tG10\n")


def test_interim_builder_reproduces_each_recorded_choice(tmp_path):
    """00_build_interim.py on a synthetic raw tree: the A13 choices, and A16's.

    Each assertion pins one choice that the frozen files were found to embody:
    the primary filter on the parquet's `sample`, the filename fallback, the LAST
    encoder layer, TOIL's sample order and scale, first-symbol-wins, patients
    kept for ANY TOIL sample, the FIRST sample per patient (an adjacent normal,
    here -- that is A16), and no DataFrame.attrs in the written parquet.
    """
    import pyarrow.parquet as _pq

    m = _load_script("00_build_interim.py")
    raw, out = tmp_path / "raw", tmp_path / "interim"
    _synthetic_raw(raw)
    assert m.main(["--raw", str(raw), "--out", str(out)]) == 0

    meta = pd.read_parquet(out / "meta.parquet")
    assert list(meta["filename"].str[:12]) == [
        "TCGA-05-0001", "TCGA-05-0002", "TCGA-AA-0003", "TCGA-05-0005"]
    assert meta.loc[1, "sample"] == "TCGA-05-0002-01Z-00-DX1"   # recovered from filename
    X = np.load(out / "X.npy")
    src = _pq.read_table(raw / "provgigapath" / "embeddings.parquet").column("embedding").to_pylist()
    assert np.array_equal(X[0], np.asarray(src[0][-1])), "not the last encoder layer"

    expr = pd.read_parquet(out / "expression_hugo.parquet")
    assert list(expr.index) == ["TCGA-05-0001-11", "TCGA-05-0001-01",
                                "TCGA-05-0002-01", "TCGA-AA-0003-01"]
    assert list(expr.columns) == [f"G{i}" for i in range(1, 11)]   # duplicate and unmapped dropped
    toil = pd.read_csv(raw / "expression" / "tcga_RSEM_gene_tpm.gz", sep="\t", index_col=0,
                       float_precision="round_trip")
    assert expr["G1"].tolist() == toil.loc["ENSG1.1"].tolist(), "the FIRST id for a symbol wins"

    cohort = pd.read_parquet(out / "cohort_nsclc.parquet")
    assert list(cohort["patient_id"]) == ["TCGA-05-0001", "TCGA-05-0002"]
    assert "cancer_type" in cohort.columns and {"SIG_HALLMARK_X", "SIG_HALLMARK_Y"} <= set(cohort)
    en = pd.read_parquet(out / "expr_nsclc.parquet")
    assert en.loc["TCGA-05-0001"].tolist() == expr.loc["TCGA-05-0001-11"].tolist(), (
        "the frozen NSCLC files take the FIRST TOIL sample per patient (A16)")
    assert b"PANDAS_ATTRS" not in (_pq.read_schema(out / "cohort_nsclc.parquet").metadata or {})

    # --tumour-only is the A16 correction: the tumour sample, into another directory.
    fixed = tmp_path / "fixed"
    assert m.main(["--raw", str(raw), "--out", str(fixed), "--tumour-only"]) == 0
    enf = pd.read_parquet(fixed / "expr_nsclc.parquet")
    assert enf.loc["TCGA-05-0001"].tolist() == expr.loc["TCGA-05-0001-01"].tolist()
    assert "TCGA-05-0001-11" not in pd.read_parquet(fixed / "expression_hugo.parquet").index

    with pytest.raises(SystemExit, match="refusing to overwrite"):
        m.main(["--raw", str(raw), "--out", str(out)])
    with pytest.raises(SystemExit, match="must not write into data/interim"):
        m.main(["--raw", str(raw), "--out", str(m.ROOT / "data" / "interim"), "--tumour-only"])


def test_the_outcome_category_test_is_exact_and_matches_its_record(tmp_path):
    """31_outcome_category_test.py (E11): an exact enumeration, reproducible.

    The count P is a hypergeometric identity and the relabelling count is
    C(16, 6); both are asserted from first principles, and the committed JSON
    must equal a fresh run. Reads only frozen results, so it runs in a deposit.
    """
    import json as _json
    import math as _math

    m = _load_script("31_outcome_category_test.py")
    out = tmp_path / "e11.json"
    assert m.main(["--out", str(out)]) == 0
    got = _json.loads(out.read_text())
    assert got["winners_inside_group"] == 5 and got["group_size"] == 6
    assert got["p_all_winners_inside_group"] == _math.comb(6, 5) / _math.comb(16, 5)
    assert got["relabellings_enumerated"] == _math.comb(16, 6)
    assert 0 < got["p_excess_difference_one_sided"] <= 1
    committed = SCRIPTS.parent / "results" / "outcome_category_test.json"
    if committed.exists():
        assert _json_close(_json.loads(committed.read_text()), got), "the committed record is stale"


def test_tiff_output_is_rgb_at_600_dpi_and_composites_translucency_on_white(tmp_path):
    """F11.2: `--tiff` writes RGB, not RGBA, and a translucent pixel is blended.

    Matplotlib stores a fully transparent pixel as transparent WHITE, so simply
    dropping alpha is harmless there. It is not harmless for a translucent pixel,
    whose stored colour is the unblended one. Half-transparent black over a
    transparent face must come out mid-grey, not black.
    """
    import matplotlib.pyplot as _plt
    from PIL import Image

    m = _load_script("10_make_figures.py")
    fig, ax = _plt.subplots(figsize=(1, 1))
    ax.plot([0, 1], [0, 1], color="#0072B2")
    m.save_tiff_rgb(fig, tmp_path / "opaque.tiff")
    fig.patch.set_facecolor((0, 0, 0, 0))
    ax.patch.set_facecolor((0, 0, 0, 0))
    ax.lines[0].remove()
    ax.set_axis_off()
    ax.axvspan(0, 1, color="black", alpha=0.5)
    with _plt.rc_context({"savefig.facecolor": (0, 0, 0, 0), "savefig.bbox": "standard"}):
        m.save_tiff_rgb(fig, tmp_path / "translucent.tiff")
    _plt.close(fig)
    for name in ("opaque.tiff", "translucent.tiff"):
        im = Image.open(tmp_path / name)
        assert im.mode == "RGB", (name, im.mode)
        assert tuple(round(v) for v in im.info["dpi"]) == (600, 600)
    im = Image.open(tmp_path / "translucent.tiff")
    centre = im.getpixel((im.size[0] // 2, im.size[1] // 2))
    assert all(100 <= c <= 160 for c in centre), f"not blended onto white: {centre}"


def test_plage_matches_an_independent_derivation_and_fixes_its_sign():
    """`score_plage` (E6) against the leading eigenvector of the gene Gram matrix.

    A different algorithm for the same quantity: if the SVD route and the
    eigendecomposition route agree, the scorer computes PLAGE and not something
    near it. Also: the sign follows the set's mean z-score; a per-gene affine
    rescaling changes nothing (genes are standardised); a zero-variance gene is
    dropped rather than poisoning the set; and the audit's dispatch reaches it.
    """
    from aacr27 import experiment as _exp
    from aacr27 import signatures as _sig

    rng = np.random.default_rng(3)
    n, genes = 60, [f"G{i}" for i in range(12)]
    latent = rng.normal(size=n)
    X = rng.normal(size=(n, 12)) + np.outer(latent, np.linspace(0.5, 2.0, 12))
    expr = pd.DataFrame(X, columns=genes, index=[f"P{i}" for i in range(n)])
    sset = _sig.SignatureSet(name="t", sets={"A": genes[:8], "B": genes[4:]})
    got = _sig.score_plage(expr, sset)

    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    for name, cols in (("A", slice(0, 8)), ("B", slice(4, 12))):
        block = Z[:, cols]
        w = np.linalg.eigh(block.T @ block)[1][:, -1]
        ref = block @ w
        ref /= np.linalg.norm(ref)
        ref *= np.sign(ref @ block.mean(1))
        assert np.allclose(got[name].to_numpy(), ref, atol=1e-10), name
        assert got[name].to_numpy() @ block.mean(1) > 0, "sign must follow the mean z"
        assert np.isclose(np.linalg.norm(got[name].to_numpy()), 1.0)

    scaled = expr * np.linspace(0.1, 10, 12) + 7.0
    assert np.allclose(_sig.score_plage(scaled, sset).to_numpy(), got.to_numpy(), atol=1e-10)

    flat = expr.copy()
    flat["G0"] = 1.0
    only = _sig.score_plage(flat, _sig.SignatureSet(name="t", sets={"A": genes[:8]}))["A"]
    Zr = (X[:, 1:8] - X[:, 1:8].mean(0)) / X[:, 1:8].std(0, ddof=1)
    w = np.linalg.eigh(Zr.T @ Zr)[1][:, -1]
    ref = Zr @ w
    ref = ref / np.linalg.norm(ref) * np.sign((Zr @ w) @ Zr.mean(1))
    assert np.allclose(only.to_numpy(), ref, atol=1e-10), "a flat gene must be dropped"

    via = _exp._score_sets(expr, sset, _exp.AuditConfig(scorer="plage"))
    assert np.array_equal(via.to_numpy(), got.to_numpy())
    with pytest.raises(ValueError, match="'plage'"):
        _exp._score_sets(expr, sset, _exp.AuditConfig(scorer="nope"))


def test_interim_verify_names_a_known_platform_build_and_nothing_else(tmp_path, monkeypatch, capsys):
    """`verify` has three answers without a reference: the frozen bytes, a
    recorded platform build (accepted and NAMED as that build), and anything
    else (a failure)."""
    m = _load_script("00_build_interim.py")
    frozen, linux, other = (tmp_path / n for n in ("a.parquet", "b.parquet", "c.parquet"))
    for f, body in ((frozen, b"frozen"), (linux, b"linux"), (other, b"other")):
        f.write_bytes(body)
    monkeypatch.setattr(m, "RECORDED", {"cohort_nsclc.parquet": m.sha16(frozen)})
    monkeypatch.setattr(m, "RECORDED_PLATFORM",
                        {"cohort_nsclc.parquet": {m.sha16(linux): "a test build"}})
    assert m.verify({"cohort_nsclc.parquet": frozen}, None) == 0
    assert m.verify({"cohort_nsclc.parquet": linux}, None) == 0
    assert "KNOWN BUILD" in capsys.readouterr().out
    assert m.verify({"cohort_nsclc.parquet": other}, None) == 1
    monkeypatch.setattr(m, "RECORDED_PLATFORM", {})
    assert m.verify({"cohort_nsclc.parquet": linux}, None) == 1, "the platform record is what accepted it"


def test_the_figure_comparison_narrows_only_on_a_different_freetype(monkeypatch):
    """The C2 guard must narrow to SVG off the committed FreeType, and only there."""
    import matplotlib.ft2font as _ft

    monkeypatch.setitem(globals(), "_FIGURES_FREETYPE", _ft.__freetype_version__)
    assert _comparable_figure_suffixes() == {".png", ".pdf", ".svg"}
    monkeypatch.setitem(globals(), "_FIGURES_FREETYPE", "0.0.0-not-a-freetype")
    assert _comparable_figure_suffixes() == {".svg"}, (
        "on a different FreeType the SVGs must still be compared, and nothing else")


def test_the_supplement_transcribes_its_source_csvs_exactly():
    """An INDEPENDENT companion to the render-sync test above (F6.4).

    That test asks whether the committed tables are what the generator writes; a
    generator that mis-reads its inputs passes it. Until 2026-09-16 this one did:
    pandas' default CSV parser put 16 of S6's 44 cells a ULP away from their
    source. Read both sides with the csv module and float(), which is exact.
    """
    import csv as _csv

    results = SCRIPTS.parent / "results"
    sup = results / "supplementary"
    if not (sup / "S6_partition_variance.tsv").exists():
        pytest.skip("the rendered tables are not staged in this tree")

    def rows(p, delim=","):
        with open(p, newline="") as fh:
            return list(_csv.DictReader(fh, delimiter=delim))

    compared, bad = 0, []
    src = {}
    for sweep, suffix in (("frozen", ""), ("stable24", "_stable24")):
        for co, lab in (("nsclc", "NSCLC"), ("pancancer", "pan-TCGA")):
            for r in rows(results / f"partition_variance_{co}{suffix}" / "per_partition.csv"):
                src[(sweep, lab, r["split_seed"])] = r
    s8 = rows(sup / "S6_partition_variance.tsv", "\t")
    assert len(s8) == len(src), "S6 and its sources differ in row count"
    for r in s8:
        for k in ("isi", "lo", "hi", "runtime_s"):
            compared += 1
            if float(r[k]) != float(src[(r["sweep"], r["cohort"], r["split_seed"])][k]):
                bad.append(f"S6 {r['sweep']} {r['cohort']} {r['split_seed']} {k}")
    s3 = rows(sup / "S1_per_signature_excess.tsv", "\t")
    for co, lab in (("pancancer_v3", "pan-TCGA"), ("nsclc_v3", "NSCLC")):
        stored = {r["signature"]: r for r in rows(results / co / "immune_excess.csv")}
        cols = [c for c in ("excess_z", "excess_lo", "excess_hi", "p_two_sided", "q_bh")
                if c in next(iter(stored.values()))]
        for r in (x for x in s3 if x["cohort"] == lab):
            for c in cols:
                compared += 1
                if float(r[c]) != float(stored[r["signature"]][c]):
                    bad.append(f"S1 {lab} {r['signature']} {c}")
    assert compared >= 59 * 4 + 32 * 3, f"only {compared} cells compared; the check has gone vacuous"
    assert not bad, f"{len(bad)} supplementary cell(s) differ from their source: {bad[:8]}"


@pytest.mark.skipif(not _ON_WORKING_TREE,
                    reason="handoff prompts are working documents and are never deposited")
def test_handoff_prompts_follow_the_name_two_prose_scans_rely_on():
    """Two prose scans skip any file with `HANDOFF` in its name (F8.12).

    That exclusion is only as good as the naming. A handoff saved under another
    name would be scanned (its superseded numbers would go red), and an
    unrelated document with `HANDOFF` in its name would silently escape every
    scan. Both directions are checked, against what the files actually say.
    """
    import re as _re

    repo = SCRIPTS.parents[1]
    named = sorted(p.name for p in repo.glob("*.md") if "HANDOFF" in p.name)
    assert len(named) >= 20, f"only {len(named)} handoff prompts found; the glob has gone vacuous"
    bad_names = [n for n in named if not _re.fullmatch(r"\d{2}-HANDOFF-PROMPT\.md", n)]
    assert not bad_names, f"files the scans skip without being handoffs by name: {bad_names}"
    def looks_like_a_handoff(p):
        return "Master prompt #" in p.read_text(errors="replace")[:200]

    # The positive control: the content signature must recognise every file that
    # IS a handoff, or the check below could never fire.
    unsigned = [n for n in named if not looks_like_a_handoff(repo / n)]
    assert not unsigned, f"handoff prompts the content signature does not recognise: {unsigned}"
    unnamed = [p.name for p in repo.glob("*.md")
               if "HANDOFF" not in p.name and looks_like_a_handoff(p)]
    assert not unnamed, f"handoff prompts saved under another name are being scanned: {unnamed}"


def test_repro_probe_is_deterministic_and_names_the_first_changed_stage(tmp_path, capsys):
    """14_repro_probe.py, the last open "pure and cheap" script (F10.1).

    Its whole purpose is that the FIRST differing line between two runs names the
    stage that moved. So: identical inputs must print identical output, and a
    one-ULP change to one expression value must first show at stage 1 (the
    input digest) and then at stage 3 (the mean expression), with the gene
    order (stage 2) untouched.
    """
    m = _load_script("14_repro_probe.py")
    rng = np.random.default_rng(14)
    n, g = 40, 300
    genes = [f"G{i}" for i in range(g)]
    patients = [f"TCGA-05-{i:04d}" for i in range(n)]
    interim = tmp_path / "interim"
    interim.mkdir()
    pd.DataFrame({"patient_id": patients}).to_parquet(interim / "cohort_nsclc.parquet")
    expr = pd.DataFrame(rng.normal(size=(n, g)), index=patients, columns=genes)
    expr.to_parquet(interim / "expr_nsclc.parquet")
    np.save(interim / "X_nsclc.npy", rng.normal(size=(n, 8)))
    gmt = tmp_path / "sets.gmt"
    gmt.write_text("A\turl\t" + "\t".join(genes[:12]) + "\nB\turl\t" + "\t".join(genes[50:62]) + "\n")
    m.INTERIM, m.GMT = interim, gmt

    def run():
        assert m.main() == 0
        return [ln for ln in capsys.readouterr().out.splitlines() if ln[:2] in {f"{k} " for k in "1234567"}]

    first, second = run(), run()
    assert first == second, "the probe is not deterministic on identical inputs"
    assert {ln[0] for ln in first} == set("1234567"), "a stage printed nothing"

    # One ULP in one of 40 values is absorbed by the column mean (measured: only
    # stage 1 moves), so the input digest is the only stage that is guaranteed
    # to see a change that small -- which is why it is printed first.
    expr.iloc[3, 7] = float(np.nextafter(expr.iloc[3, 7], np.inf))
    expr.to_parquet(interim / "expr_nsclc.parquet")
    third = run()
    changed = [a[0] for a, b in zip(first, third) if a != b]
    assert changed and changed[0] == "1", f"a changed input was not flagged first at stage 1: {changed}"

    expr.iloc[3, 7] += 1e-9
    expr.to_parquet(interim / "expr_nsclc.parquet")
    fourth = run()
    changed = [a[0] for a, b in zip(first, fourth) if a != b]
    assert changed[0] == "1" and "3" in changed and "2" not in changed, (
        f"a changed value must move the mean-expression stage and not the gene order: {changed}")


def test_the_figure_data_colours_are_the_published_okabe_ito_values():
    """F5.5: the palette test above checks two files AGAINST EACH OTHER.

    A colour mistyped identically in both would pass it. This one derives the
    palette from the published RGB triplets instead (Okabe M, Ito K. Color
    Universal Design: how to make figures and presentations that are friendly
    to colorblind people, 2008), and requires every DATA colour the figure
    script assigns -- the cohort colours and the category map -- to be one of
    them. C_NULL is grey on purpose (the null is not a category) and is the
    only exception.
    """
    import ast

    published_rgb = {
        "black": (0, 0, 0), "orange": (230, 159, 0), "sky blue": (86, 180, 233),
        "bluish green": (0, 158, 115), "yellow": (240, 228, 66), "blue": (0, 114, 178),
        "vermillion": (213, 94, 0), "reddish purple": (204, 121, 167),
    }
    palette = {"#{:02X}{:02X}{:02X}".format(*rgb) for rgb in published_rgb.values()}

    tree = ast.parse((SCRIPTS / "10_make_figures.py").read_text())
    data_colours = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        tuples = [t for t in node.targets if isinstance(t, ast.Tuple)]
        if tuples and isinstance(node.value, ast.Tuple):
            for t, v in zip(tuples[0].elts, node.value.elts):
                if isinstance(t, ast.Name) and t.id.startswith("C_") and isinstance(v, ast.Constant):
                    data_colours[t.id] = v.value
        if names == ["CAT_COLOUR"] and isinstance(node.value, ast.Dict):
            for k, v in zip(node.value.keys, node.value.values):
                data_colours[f"CAT_COLOUR[{k.value}]"] = v.value
    assert len(data_colours) >= 8, f"found only {data_colours}; the AST walk has gone vacuous"
    off = {k: v for k, v in data_colours.items()
           if v.upper() not in palette and k != "C_NULL"}
    assert not off, f"data colours outside the published Okabe-Ito set: {off}"
    assert data_colours["C_NULL"].upper() not in palette, (
        "C_NULL was meant to sit OUTSIDE the palette; if that changed, say why here")


def test_the_fetch_script_verifies_every_raw_input_and_derives_what_it_can(tmp_path, capsys):
    """01_fetch_data.py, rewritten 2026-09-17: every raw file is hash-verified.

    Offline: the network sources are never reached because the file that would
    be downloaded is already present. What is exercised is the part a reader
    depends on -- the copy from resources/, the tab-to-comma derivation, the
    MISSING and HASH DIFFERS verdicts, and the exit code.
    """
    import hashlib as _h

    m = _load_script("01_fetch_data.py")
    # The real table: exact paths the pipeline reads, full hashes, the purity
    # .tsv listed before the .csv derived from it.
    real = list(m.RAW_FILES)
    assert len(real) >= 8 and all(len(v[0]) == 64 for v in m.RAW_FILES.values())
    assert real.index("purity/TCGA_ABSOLUTE_purity.tsv") < real.index("purity/TCGA_ABSOLUTE_purity.csv")
    assert (SCRIPTS.parent / "resources" / "ensembl_to_hugo.csv").exists()
    assert "embeddings.parquet" in " ".join(real), "the builder reads embeddings.parquet"

    raw, res = tmp_path / "raw", tmp_path / "resources"
    res.mkdir()
    gene_map = b"ensembl_gene_id,symbol\nENSG1,A\n"
    (res / "ensembl_to_hugo.csv").write_bytes(gene_map)
    tsv = b"array\tpurity\nTCGA-01\t0.5\n"
    (raw / "purity").mkdir(parents=True)
    (raw / "purity" / "p.tsv").write_bytes(tsv)
    sha = lambda b: _h.sha256(b).hexdigest()  # noqa: E731
    m.RAW, m.RESOURCES = raw, res
    m.RAW_FILES = {
        "expression/ensembl_to_hugo.csv": (sha(gene_map), "resources"),
        "purity/p.tsv": (sha(tsv), "https://invalid.example/never-fetched"),
        "purity/p.csv": (sha(tsv.replace(b"\t", b",")), "tabs-to-commas"),
        "signatures/s.gmt": (sha(b"x"), "03_fetch_signatures.py"),
    }
    assert m.main(["--all"]) == 1
    out = capsys.readouterr().out
    assert (raw / "expression" / "ensembl_to_hugo.csv").read_bytes() == gene_map
    assert (raw / "purity" / "p.csv").read_bytes() == b"array,purity\nTCGA-01,0.5\n"
    assert "MISSING        signatures/s.gmt" in out and out.count("VERIFIED") == 3

    (raw / "signatures").mkdir(exist_ok=True)   # the failed fetch already made it
    (raw / "signatures" / "s.gmt").write_bytes(b"y")
    assert m.main(["--verify"]) == 1
    assert "HASH DIFFERS   signatures/s.gmt" in capsys.readouterr().out
    (raw / "signatures" / "s.gmt").write_bytes(b"x")
    assert m.main(["--verify"]) == 0

    hgnc = tmp_path / "hgnc.txt"
    hgnc.write_text("hgnc_id\tsymbol\tensembl_gene_id\nH:1\tA\tENSG1\nH:2\tB\t\n")
    assert m.derive_gene_map(hgnc, tmp_path / "derived.csv") == sha(gene_map), (
        "the map is HGNC's rows that carry an Ensembl id, two columns, in order")


def _json_close(a, b, rel=1e-12) -> bool:
    """Equal structure and strings; floats equal to `rel`."""
    import math as _math

    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_json_close(a[k], b[k], rel) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_json_close(x, y, rel) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        return _math.isclose(a, b, rel_tol=rel, abs_tol=0.0)
    return a == b


def test_json_close_is_not_vacuous():
    assert _json_close({"a": [1, 0.1 + 0.2]}, {"a": [1, 0.3]})
    assert not _json_close({"a": 0.3}, {"a": 0.3 + 1e-9})
    assert not _json_close({"a": 1}, {"a": 1, "b": 2})
    assert not _json_close({"a": "x"}, {"a": "y"})


def test_the_outcome_power_curve_reproduces_the_stated_mde(tmp_path):
    """30_outcome_power.py (E10): the curve must agree with the MDE it extends.

    The MDE the manuscript quotes (0.046 in C-index for NSCLC) and the power
    curve come from one normal approximation, so power at the MDE is 0.80 by
    construction, and the pan-cancer table's STORED MDE column must equal the
    script's recomputation. Reads only frozen results, so it runs in a deposit.
    """
    import json as _json

    m = _load_script("30_outcome_power.py")
    for cohort in ("nsclc", "pancancer"):
        out = tmp_path / f"{cohort}.json"
        assert m.main(["--cohort", cohort, "--out", str(out)]) == 0
        got = _json.loads(out.read_text())
        assert got["n_comparisons"] == 32
    nsclc = _json.loads((tmp_path / "nsclc.json").read_text())
    assert round(nsclc["mde_80pct_c_index_median"], 3) == 0.046, "the manuscript's MDE"
    curve = nsclc["power_by_c_index_advantage"]
    assert curve["0.02"]["median"] < curve["0.03"]["median"] < curve["0.04"]["median"] < 0.80
    committed = SCRIPTS.parent / "results" / "outcome_power_nsclc.json"
    if committed.exists():
        # Staleness, not bit identity: scipy's normal cdf and libm's atanh round
        # differently on Linux/x86_64 (a ULP in the power values, HPC4 job
        # 130188), so the leaves are compared to 1e-12.
        assert _json_close(_json.loads(committed.read_text()), nsclc), (
            "the committed curve is stale")
    # Can fail: a C-index/Fisher-z conversion that drops the factor of 2 must trip
    # the internal power-at-MDE assertion.
    orig = m.power
    m.power = lambda d, se: orig(d / 2, se)
    with pytest.raises(AssertionError):
        m.main(["--cohort", "nsclc", "--out", str(tmp_path / "bad.json")])
    m.power = orig


# ----------------------------------------------------- ORCID check digits ---
# Session 57. The ORCID iDs sit in the manuscript, the cover letter,
# CITATION.cff, the dry run and the checklists, and nothing checked them: a
# mutation sweep of the dry run moved every ORCID digit group and `--strict`
# passed each time -- correctly, since an identifier is not a result. But an
# iD with one digit wrong attributes a paper to a stranger, and ORCID's last
# character is an ISO 7064 MOD 11-2 check digit, which catches every
# single-digit substitution and every adjacent transposition of unequal
# digits (measured on the example below: 144 of 144 and 6 of 6). So the guard
# costs nothing, and it needs no list of iDs -- none are written here.

_ORCID = re.compile(r"\b(\d{4}-\d{4}-\d{4}-\d{3}[\dX])\b")


def _orcid_check_digit(base15: str) -> str:
    """ISO 7064 MOD 11-2, as ORCID specifies it for the final character."""
    total = 0
    for c in base15:
        total = (total + int(c)) * 2
    r = (12 - total % 11) % 11
    return "X" if r == 10 else str(r)


def _orcid_valid(orcid: str) -> bool:
    d = orcid.replace("-", "")
    return len(d) == 16 and _orcid_check_digit(d[:15]) == d[15]


@pytest.mark.skipif(not _ON_WORKING_TREE,
                    reason="scans the tracked documents of the working tree")
def test_every_orcid_in_the_documents_has_a_valid_check_digit():
    repo = Path(__file__).resolve().parents[2]
    found: dict[str, set[str]] = {}
    for rel in _tracked_prose(repo):
        text = (repo / rel).read_text(encoding="utf-8", errors="replace")
        for m in _ORCID.finditer(text):
            found.setdefault(m.group(1), set()).add(rel)
    n_docs = len(set().union(*found.values())) if found else 0
    # Measured 2026-09-22: 3 iDs in 5 tracked prose documents (CITATION.cff,
    # the cover letter, the dry run, the submission checklist and
    # 16-YOUR-TASKS.md; handoffs are archives and not scanned). The first
    # draft of this floor said 10 documents, counted with the handoffs in,
    # and the test refused it. 2026-09-24: the author sent three more
    # (S. Kodilkar, A. Raut, S. Raut), so six are on file and the floor
    # moved with them -- a floor left at three would let three vanish.
    assert len(found) >= 6 and n_docs >= 4, (
        f"only {len(found)} ORCID iD(s) across {n_docs} document(s) -- six "
        "are on file, in five documents, so the scan has lost some")
    bad = {o: sorted(fs) for o, fs in found.items() if not _orcid_valid(o)}
    assert not bad, f"ORCID iD(s) failing the MOD 11-2 check digit: {bad}"


def test_the_orcid_check_digit_rule_can_fail():
    # ORCID's own documentation example iD (not a real researcher's).
    good = "0000-0002-1825-0097"
    assert _orcid_valid(good)
    pos = [i for i, c in enumerate(good) if c.isdigit()]
    for i in pos:
        for d in "0123456789":
            if d != good[i]:
                bad = good[:i] + d + good[i + 1:]
                assert not _orcid_valid(bad), f"missed a one-digit error: {bad}"
    swaps = [(a, b) for a, b in zip(pos, pos[1:]) if b == a + 1 and good[a] != good[b]]
    assert swaps, "the example has no adjacent unequal digits to transpose"
    for a, b in swaps:
        bad = good[:a] + good[b] + good[a] + good[b + 1:]
        assert not _orcid_valid(bad), f"missed a transposition: {bad}"


# ------------------------------------------ markdown tables: dropped cells ---
# Session 58. GitHub-flavoured markdown splits a table row on EVERY unescaped
# pipe -- inside a code span too -- pads a row that has too few cells, and
# DROPS every cell past the header's count. Five rows of `results/README.md`
# had more cells than their header: two from a bare absolute-value bar in prose
# (`|r/sqrt(α)|`, `|Δ|`), three from a method cell pasted twice or a note
# appended as a new cell. 1,612 characters never rendered, among them E9's
# smallest null reliabilities and A10's score shift in standard deviations.
# `19_check_numbers.py` reads the raw text, so those numbers were checked all
# along; only a reader on GitHub lost them. Short rows are left alone: they
# lose nothing (three in `results/README.md` put class and notes in one cell).

_GFM_DELIM = re.compile(r"^\s*\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)*\|?\s*$")


def _gfm_cells(line: str) -> int:
    s = line.strip()
    parts = re.split(r"(?<!\\)\|", s)
    if s.startswith("|"):
        parts = parts[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        parts = parts[:-1]
    return len(parts)


def _gfm_overlong_rows(text: str) -> tuple[int, int, list[tuple[int, int, int]]]:
    """(tables, body rows, [(line, header cells, row cells)] for rows too long)."""
    lines = text.split("\n")
    tables = rows = 0
    bad: list[tuple[int, int, int]] = []
    fence = False
    i = 0
    while i < len(lines):
        if lines[i].lstrip().startswith(("```", "~~~")):
            fence = not fence
        elif (not fence and lines[i].lstrip().startswith("|")
              and i + 1 < len(lines) and _GFM_DELIM.match(lines[i + 1])):
            n = _gfm_cells(lines[i])
            tables += 1
            j = i + 2
            while j < len(lines) and lines[j].lstrip().startswith("|"):
                rows += 1
                if _gfm_cells(lines[j]) > n:
                    bad.append((j + 1, n, _gfm_cells(lines[j])))
                j += 1
            i = j
            continue
        i += 1
    return tables, rows, bad


@pytest.mark.skipif(not _ON_WORKING_TREE,
                    reason="scans the tracked documents of the working tree")
def test_no_markdown_table_row_has_more_cells_than_its_header():
    # CAN FAIL, first: a bare bar in a cell, and a cell pasted twice.
    planted = "| a | b |\n|---|---|\n| x | y |r| z |\n| x | y | y |\n| x \\| y | z |\n"
    assert _gfm_overlong_rows(planted) == (1, 3, [(3, 2, 4), (4, 2, 3)])
    repo = Path(__file__).resolve().parents[2]
    tables = rows = 0
    bad: dict[str, list] = {}
    for rel in sorted(f for f in _tracked_prose(repo) if f.endswith(".md")):
        t, r, b = _gfm_overlong_rows((repo / rel).read_text(encoding="utf-8"))
        tables, rows = tables + t, rows + r
        if b:
            bad[rel] = b
    # Measured 2026-09-23: 171 tables and 1,124 body rows across the tracked
    # markdown documents, handoffs excluded.
    assert tables >= 150 and rows >= 1000, (
        f"only {tables} table(s) and {rows} row(s) found -- the scan has gone "
        "near-vacuous")
    assert not bad, (
        "table rows with more cells than their header; GitHub drops the excess, "
        f"so that text never renders. Escape a bare bar as \\| or merge the cell: {bad}")


# ---------------------------------------------------------------------------
# Session 59 (B-19, decided by the author 2026-09-24): the three held-back code
# fixes, each with the behavior it replaces pinned down.
# ---------------------------------------------------------------------------

def test_site_control_with_no_qualifying_site_returns_an_empty_typed_frame():
    """B-19b: it raised `KeyError: 'auroc'` (E16, array 130207)."""
    rng = np.random.default_rng(0)
    X = rng.normal(size=(60, 4))
    sites = pd.Series([f"S{i % 12}" for i in range(60)])          # 5 per site
    patients = pd.Series([f"P{i}" for i in range(60)])
    out = models.site_prediction_control(X, sites, patients, min_site_n=20)
    assert out.empty and list(out.columns) == ["site", "n_patients", "auroc"]


def test_plate_within_site_is_nan_when_the_plate_carries_no_information():
    """B-19a: an all-missing plate used to give a structural 0 (E17)."""
    from aacr27 import decomposition

    rng = np.random.default_rng(1)
    n = 240
    frame = pd.DataFrame({
        "tss": [f"S{i % 6}" for i in range(n)],
        "cancer_type": ["A"] * n,
        "purity": rng.uniform(0.3, 0.9, n),
        "plate": [None] * n,
        "sig": rng.normal(size=n),
    })
    out = decomposition.label_site_variance(frame, ["sig"])
    assert np.isnan(out["r2_plate_within_site"]).all()
    assert np.isnan(out["r2_site_plus_plate"]).all()
    assert "not estimable" in out.attrs["plate_note"]
    # ...and with two or more real plate levels the rung is computed as before.
    frame["plate"] = [f"PL{i % 4}" for i in range(n)]
    out2 = decomposition.label_site_variance(frame, ["sig"])
    assert np.isfinite(out2["r2_plate_within_site"]).all()
    assert "plate_note" not in out2.attrs


def test_ols_fallback_matches_the_default_fit_and_is_loud(monkeypatch):
    """B-19c: when statsmodels' SVD fails, the pivoted-QR refit must give the
    same R^2 and adjusted R^2 as the default would have, and must warn."""
    from aacr27 import decomposition
    import statsmodels.api as sm

    rng = np.random.default_rng(2)
    X = np.column_stack([rng.normal(size=300), (np.arange(300) % 7 == 0).astype(float)])
    y = X @ np.array([0.5, 1.0]) + rng.normal(size=300)
    want = decomposition._r2(y, X)

    real_ols = sm.OLS

    class _Failing(real_ols):
        def fit(self, *a, **k):
            raise np.linalg.LinAlgError("SVD did not converge")

    monkeypatch.setattr(decomposition.sm, "OLS", _Failing)
    with pytest.warns(RuntimeWarning, match="pivoted QR"):
        got = decomposition._r2(y, X)
    assert got[2:] == want[2:]
    assert abs(got[0] - want[0]) < 1e-10 and abs(got[1] - want[1]) < 1e-10
