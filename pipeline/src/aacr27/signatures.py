"""TME signature scoring, and the random-gene-set null.

Two responsibilities:

1. Score curated gene-expression signatures from a bulk RNA matrix. These are
   the regression targets the image model tries to predict.

2. Generate size-matched (optionally expression-matched) RANDOM gene sets and
   score them identically — the Venet null.

   NOTE ON SCOPE: applied to EXPRESSION targets this null cannot by itself
   support a 'the model learned nothing immune-specific' conclusion, because a
   dominant global expression axis makes random-set predictability the expected
   result. See `globalaxis` for the corrected estimand and `outcome` for the
   version of this null run against a clinical endpoint, where Venet's logic
   actually applies.

Why the null matters
--------------------
Venet, Dumont & Detours, PLoS Comput Biol 2011;7:e1002240
(https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002240)
showed that >90% of random gene sets with >100 genes significantly predict
breast cancer outcome, and that 28 of 47 published breast signatures were no
better than same-size random signatures.

PRIOR ART — READ BEFORE CLAIMING NOVELTY.
HE2RNA (Schmauch et al., Nat Commun 2020;11:3877, doi:10.1038/s41467-020-17678-4)
ALREADY ran a size-matched random-gene-list null against image-predicted
expression: "the mean correlation coefficient R0 obtained for 10,000 random
lists of the same number of genes, for all 28 cancer types". Its published
answer favours the signatures — 75% (B-cell) and 86% (T-cell) of cancer types
had immune signatures significantly better predicted than random, though only
36% and 50% on a per-gene metric.

So a bare random-set null is NOT a contribution. What is open is that HE2RNA's
comparison is uncorrected for two things, both of which bias it toward the
curated set:
  1. the global expression axis (see `globalaxis`) — "beats random" is equally
     consistent with the image reading composition faithfully;
  2. reliability — curated sets are co-expressed modules and are therefore more
     predictable than random draws regardless of biology (see `cronbach_alpha`).

IMPORTANT — gene sets are NOT bundled
-------------------------------------
This module deliberately ships no curated gene lists other than one small,
unambiguously citable example. Signature definitions must be fetched from their
sources so that the provenance is auditable and the licence is respected. Run
`scripts/01_fetch_data.py --signatures` for retrieval instructions. The one
exception is CYT (below), a published two-gene metric included as a smoke test.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

# Rooney et al., Cell 2015;160:48-61 — cytolytic activity is the geometric mean
# of GZMA and PRF1 expression. Two genes, unambiguous, safe to hard-code.
CYT = ("GZMA", "PRF1")


@dataclass
class SignatureSet:
    """A named collection of gene sets."""

    name: str
    sets: dict[str, list[str]] = field(default_factory=dict)
    source: str = ""

    def __len__(self) -> int:
        return len(self.sets)

    @property
    def sizes(self) -> dict[str, int]:
        return {k: len(v) for k, v in self.sets.items()}

    def filter_to(self, available: set[str], *, min_genes: int = 5) -> "SignatureSet":
        """Restrict every set to genes present in the expression matrix.

        Sets falling below `min_genes` after filtering are dropped and reported,
        because scoring a 3-gene remnant of a 40-gene signature and calling it by
        the original name is a silent misattribution.
        """
        kept, dropped = {}, {}
        for name, genes in self.sets.items():
            present = [g for g in genes if g in available]
            if len(present) >= min_genes:
                kept[name] = present
            else:
                dropped[name] = (len(present), len(genes))
        out = SignatureSet(name=self.name, sets=kept, source=self.source)
        if dropped:
            out.coverage_dropped = dropped  # type: ignore[attr-defined]
        return out

    @classmethod
    def from_json(cls, path: str | Path, *, name: str | None = None) -> "SignatureSet":
        path = Path(path)
        payload = json.loads(path.read_text())
        sets = payload.get("sets", payload)
        return cls(
            name=name or payload.get("name", path.stem),
            sets={k: list(v) for k, v in sets.items()},
            source=payload.get("source", str(path)),
        )

    @classmethod
    def from_gmt(cls, path: str | Path, *, name: str | None = None) -> "SignatureSet":
        """Read a GMT file (MSigDB format): name<TAB>description<TAB>gene1<TAB>..."""
        path = Path(path)
        sets: dict[str, list[str]] = {}
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sets[parts[0]] = [g for g in parts[2:] if g]
        return cls(name=name or path.stem, sets=sets, source=str(path))

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps({"name": self.name, "source": self.source, "sets": self.sets}, indent=2)
        )


def example_signature_set() -> SignatureSet:
    """A one-entry set for smoke-testing the pipeline end to end."""
    return SignatureSet(
        name="example",
        sets={"CYT_Rooney2015": list(CYT)},
        source="Rooney et al., Cell 2015;160:48-61 (doi:10.1016/j.cell.2014.12.033)",
    )


# ---------------------------------------------------------------- scoring ---

def score_mean_z(expr: pd.DataFrame, sig: SignatureSet) -> pd.DataFrame:
    """Mean of per-gene z-scores. Fast, transparent, and the default here.

    `expr` is samples x genes, already log-transformed. Z-scoring is done across
    samples within the supplied matrix, which means scores are cohort-relative.
    If you score TCGA and CPTAC separately you have silently redefined the target
    between discovery and validation — fix the scaling policy once and apply it
    identically to both. See scripts/01_fetch_data.py --expression.
    """
    z = (expr - expr.mean(axis=0)) / (expr.std(axis=0).replace(0, np.nan))
    out = {}
    for name, genes in sig.sets.items():
        cols = [g for g in genes if g in z.columns]
        if not cols:
            continue
        out[name] = z[cols].mean(axis=1)
    return pd.DataFrame(out, index=expr.index)


# Cells per sample-block in `score_ssgsea`. The two dense per-sample tables are
# materialised one block at a time so peak memory does not scale with the cohort:
# pan-cancer is ~11k patients x 20k genes, which as one float64 table would be
# 1.8 GB. 8M cells is ~64 MB of weights plus ~32 MB of depths.
_SSGSEA_BLOCK_CELLS = 8_000_000


def _member_positions(columns: pd.Index, sig: SignatureSet) -> dict[str, np.ndarray]:
    """Column positions occupied by each gene set, resolved ONCE for all sets.

    `columns.isin(genes)` is an O(n_genes) hash sweep per gene set; at 16
    signatures x 1000 null draws that sweep alone costs more than the scoring it
    was preparing for. Inverting the column index once yields the same positions,
    including for duplicate column labels, where `isin` marks every position that
    shares the label.
    """
    by_label: dict[object, list[int]] = {}
    for pos, label in enumerate(columns):
        by_label.setdefault(label, []).append(pos)

    out: dict[str, np.ndarray] = {}
    for name, genes in sig.sets.items():
        hits = sorted({p for g in genes if g in by_label for p in by_label[g]})
        if hits:
            out[name] = np.asarray(hits, dtype=np.intp)
    return out


def _ssgsea_sample_tables(
    expr: pd.DataFrame, alpha: float
) -> tuple[np.ndarray, np.ndarray]:
    """Per-sample |rank|^alpha weights, and depth in the descending rank order.

    `depth[i, g] = n_genes - position of g` — the number of positions of the
    running enrichment sum that gene g contributes to in sample i. Both tables
    depend on the SAMPLE only and never on the gene set, which is the whole
    reason `score_ssgsea` became affordable.
    """
    ranks = expr.rank(axis=1, method="average").to_numpy(dtype=float)
    n_samples, n_genes = ranks.shape

    depth = np.empty((n_samples, n_genes), dtype=np.int32)
    countdown = np.arange(n_genes, 0, -1, dtype=np.int32)
    for i in range(n_samples):
        # `np.argsort(-row)` reproduced EXACTLY as the reference loop calls it,
        # `kind="stable"` included. Ties are common in log-expression (every gene
        # at the floor shares a rank), tied genes are ordered by argsort's
        # tie-breaking, and the score depends on that order — so this has to
        # reproduce the call, not merely its intent. Sorting the row ascending
        # and reversing, for instance, reverses tie groups and changes the answer.
        #
        # `kind="stable"` is LOAD-BEARING and was added 2026-09-05. The default
        # `kind="quicksort"` is not stable, so tied genes were ordered by an
        # implementation-defined rule that differs between numpy's arm64 and
        # x86_64 kernels — the same defect as A9 (`splits.py:177`), and on this
        # data a larger one. Measured on the real NSCLC matrix (944 x 41046):
        # every one of the 944 samples carries ties, 35,394,448 of 38,747,424
        # rank entries (91%) are tied, the depth table differs between the two
        # sort kinds, and the difference SURVIVES the O(k) reduction into the 16
        # signature scores — up to 1.282e-02 absolute against a score sd of
        # 2.27e-02, i.e. 56% of a standard deviation, on 15,060 of 15,104 cells
        # and all 16 signatures. Tied genes carry identical weights, so a swap
        # only moves the score when exactly ONE of the pair is a set member;
        # with 41,046 genes and ~200-gene sets that happens constantly.
        #
        # This never touched the frozen headline results: those run
        # `scorer="mean_z"` (see results/nsclc_v3/config.json), which does not
        # reach this function. It DID reach the ssGSEA arms of the A7
        # scorer-sensitivity analysis. Any tie convention is a valid ssGSEA —
        # the statistic is genuinely ambiguous under ties — so the fix is to
        # PIN one, not to find the "right" one. Stable pins tied genes to column
        # order, which is fixed by the input file.
        depth[i, np.argsort(-ranks[i], kind="stable")] = countdown

    # In place: `ranks` is not needed afterwards and it is the largest array here.
    np.abs(ranks, out=ranks)
    np.power(ranks, alpha, out=ranks)
    return ranks, depth


def score_ssgsea(
    expr: pd.DataFrame, sig: SignatureSet, *, alpha: float = 0.25
) -> pd.DataFrame:
    """Single-sample GSEA (Barbie et al., Nature 2009).

    Rank-based, so it is invariant to per-sample scaling and is what most of the
    TME-signature literature uses.

    WHY THIS IS NOT THE TEXTBOOK WALK
    ---------------------------------
    The naive form walks the whole ranked gene list once per sample PER GENE SET,
    recomputing the descending sort order and the |rank|^alpha weights every
    time. Those depend on the SAMPLE ONLY, so with 1000 null draws per signature
    the identical work was repeated 1000 times: measured at 6.3 ms per
    patient-set, one NSCLC ISI run (944 patients x 16 signatures x 1001 sets)
    would have taken about 24 hours, which is why the ssGSEA sensitivity analysis
    had never been run.

    This form computes the two per-sample tables once and then reduces each gene
    set in O(k) rather than O(n_genes), by an EXACT algebraic identity rather
    than an approximation. Writing the running enrichment score at position p as

        ES_p = hits_p / T_hit - misses_p / T_miss

    the statistic is sum_p ES_p / n_genes. An element at position q enters both
    cumulative sums at q and stays in them to the end, so it appears in exactly
    depth_q = n_genes - q of the terms, giving

        sum_p hits_p   = sum over MEMBERS of w_q * depth_q
        sum_p misses_p = n_genes*(n_genes+1)/2 - sum over MEMBERS of depth_q

    Both sums range over member genes alone, so a gene set costs work
    proportional to its own size and nothing else. Equality with the
    pre-optimisation loop is asserted to 1e-10 — on a tie-heavy matrix as well as
    a continuous one — in `test_ssgsea_optimised_matches_reference_loop`. Keep
    that test: a sensitivity analysis run with a subtly different ssGSEA answers
    a different question than the one it claims to. On the real NSCLC matrix,
    where 85% of each row's values are tied, the two agree to 1.8e-15.

    Measured on NSCLC (944 patients, 41046 genes, a 1000-draw null family):
    12.0 s, i.e. 0.0127 ms per patient-set against the old 6.26 ms — 494x, and
    the 24-hour ISI run becomes ~3 minutes of scoring. Most of what is left is
    the per-call `rank` + `argsort` (~9 s of the 12), which is amortised over
    however many sets are passed in ONE call — so pass the whole null family at
    once rather than looping draw by draw.
    """
    n_samples, n_genes = expr.shape
    members = _member_positions(expr.columns, sig)
    if not members:
        return pd.DataFrame({}, index=expr.index)

    results = {name: np.empty(n_samples) for name in members}
    all_depth = n_genes * (n_genes + 1) / 2.0

    block = max(1, _SSGSEA_BLOCK_CELLS // max(n_genes, 1))
    for start in range(0, n_samples, block):
        stop = min(start + block, n_samples)
        weights, depth = _ssgsea_sample_tables(expr.iloc[start:stop], alpha)

        for name, idx in members.items():
            n_miss = n_genes - len(idx)
            w, d = weights[:, idx], depth[:, idx]
            total_hit = w.sum(axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                scores = (
                    (w * d).sum(axis=1) / total_hit
                    - (all_depth - d.sum(axis=1, dtype=np.float64)) / n_miss
                ) / n_genes
            # The reference loop's degeneracy guard, unchanged: a set carrying no
            # weight, or one covering every gene, has no enrichment to measure.
            # NaN weights (a member gene with an undefined rank) still propagate
            # through the arithmetic, exactly as they did before.
            results[name][start:stop] = np.where(
                (total_hit == 0) | (n_miss == 0), np.nan, scores
            )

    return pd.DataFrame(results, index=expr.index)


# ------------------------------------------------------------- Venet null ---

def random_gene_sets(
    genes: list[str],
    reference: "SignatureSet | dict[str, int]",
    *,
    n_per_signature: int = 100,
    seed: int = 0,
    match_expression: pd.Series | None = None,
    n_bins: int = 10,
) -> dict[str, SignatureSet]:
    """Generate size-matched random gene sets, one family per real signature.

    Parameters
    ----------
    reference
        Either a `SignatureSet` (preferred — enables true expression matching,
        because the real gene lists are needed to compute the target per-bin
        histogram) or a plain {name: n_genes} dict for size-only matching.
        Pass the set AFTER `filter_to()` so sizes reflect scorable genes.
    match_expression
        Optional mean expression per gene. When supplied, random sets are drawn
        to match the real signature's *expression-level distribution*, not just
        its size. This is the stricter null: immune genes tend to be moderately
        expressed, and a purely size-matched draw can be dominated by low-count
        genes that are simply noisier. Report both if you can; the strict null
        is the one a reviewer will ask for.

    Returns
    -------
    {signature_name: SignatureSet with n_per_signature random sets}
    """
    rng = np.random.default_rng(seed)
    pool = np.asarray(genes)
    out: dict[str, SignatureSet] = {}

    # Normalise the reference into {name: gene_list or None} plus {name: size}.
    if isinstance(reference, SignatureSet):
        gene_lists = {k: list(v) for k, v in reference.sets.items()}
        sizes = {k: len(v) for k, v in gene_lists.items()}
    else:
        gene_lists = {k: None for k in reference}
        sizes = dict(reference)

    genes_by_bin: dict[int, np.ndarray] | None = None
    bin_of: dict[str, int] | None = None
    if match_expression is not None:
        shared = match_expression.reindex(pool).dropna()
        if len(shared) < len(pool) * 0.5:
            raise ValueError("match_expression covers <50% of the gene pool")
        bins = pd.qcut(shared.rank(method="first"), n_bins, labels=False).astype(int)
        bin_of = dict(zip(shared.index, bins))
        genes_by_bin = {
            b: shared.index[bins.to_numpy() == b].to_numpy() for b in range(n_bins)
        }

    for sig_name, size in sizes.items():
        real_genes = gene_lists.get(sig_name)
        matched = genes_by_bin is not None and real_genes is not None

        if genes_by_bin is not None and real_genes is None:
            import warnings

            warnings.warn(
                f"{sig_name}: match_expression was supplied but no gene list is "
                "available (a sizes dict was passed instead of a SignatureSet), so "
                "this family falls back to SIZE-ONLY matching.",
                stacklevel=2,
            )

        target = (
            _bin_histogram(real_genes, bin_of, len(genes_by_bin)) if matched else None
        )

        draws: dict[str, list[str]] = {}
        for i in range(n_per_signature):
            if target is None:
                picked = rng.choice(pool, size=min(size, len(pool)), replace=False)
            else:
                picked = _draw_expression_matched(target, genes_by_bin, rng)
            draws[f"{sig_name}__random_{i:03d}"] = list(map(str, picked))

        kind = "expression-matched" if matched else "size-matched"
        out[sig_name] = SignatureSet(
            name=f"{sig_name}__null",
            sets=draws,
            source=f"{kind} random null (n={size}, {n_per_signature} draws)",
        )
    return out


def _bin_histogram(
    genes: list[str], bin_of: dict[str, int], n_bins: int
) -> np.ndarray:
    """How many of a real signature's genes fall in each expression bin.

    This is the target the null must reproduce. Genes absent from the expression
    matrix contribute nothing, so the drawn sets can be slightly smaller than the
    nominal signature size — which is correct, since those genes are not scorable
    for the real signature either.
    """
    counts = np.zeros(n_bins, dtype=int)
    for g in genes:
        b = bin_of.get(g)
        if b is not None:
            counts[b] += 1
    return counts


def _draw_expression_matched(
    target: np.ndarray,
    genes_by_bin: dict[int, np.ndarray],
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw a random set reproducing `target`'s per-bin expression histogram.

    Previously this spread genes UNIFORMLY across bins, which matched a uniform
    expression distribution rather than the real signature's — a null that is
    systematically wrong for immune signatures, whose genes cluster in the
    moderate-expression range. It now reproduces the signature's own histogram
    bin by bin.

    If a bin holds fewer candidates than requested, the shortfall is redistributed
    to the nearest bins with capacity so the returned set still has the right size.
    """
    picked: list[str] = []
    shortfall = 0

    for b, want in enumerate(target):
        avail = genes_by_bin.get(b, np.asarray([]))
        take = min(int(want), len(avail))
        if take:
            picked.extend(rng.choice(avail, size=take, replace=False))
        shortfall += int(want) - take

    if shortfall > 0:
        already = set(picked)
        spare = np.asarray(
            [g for arr in genes_by_bin.values() for g in arr if g not in already]
        )
        if len(spare):
            picked.extend(rng.choice(spare, size=min(shortfall, len(spare)), replace=False))

    return np.asarray(picked)


def null_percentile(observed: float, null_scores: np.ndarray) -> float:
    """Where the real signature falls in its own random-set distribution.

    This is the headline number per signature. A curated signature whose
    image-predictability sits at the 55th percentile of random sets of the same
    size is not carrying specific information — which is exactly the Venet
    finding, transposed from outcome prediction to image prediction.
    """
    null_scores = np.asarray(null_scores, dtype=float)
    null_scores = null_scores[np.isfinite(null_scores)]
    if len(null_scores) == 0:
        return np.nan
    return float(sps.percentileofscore(null_scores, observed, kind="mean"))


def null_summary(
    observed: dict[str, float],
    nulls: dict[str, np.ndarray],
    *,
    expr: pd.DataFrame | None = None,
    sig: "SignatureSet | None" = None,
    null_families: dict[str, "SignatureSet"] | None = None,
    n_obs: int = 0,
) -> pd.DataFrame:
    """Per-signature comparison of real vs random predictability.

    The headline column is `excess_z` (excess over the null mean in Fisher z with
    a CI), NOT `percentile`. The percentile is retained as a labelled secondary
    because on a concentrated null it is a hair trigger: +/-0.01 in r moves it
    from the 1st to the 99th percentile. Its bias direction is stated in the
    column name so it cannot be quoted innocently.

    When `expr` and `sig` are supplied, `rho_bar` (the signature's mean pairwise
    gene-gene correlation) is reported next to `null_rho_bar`. That pair is the
    number that tells a reviewer which regime the analysis is in: comparable
    rho_bar means the null is fair, a large gap means reliability is doing the
    work and disattenuation is mandatory.
    """
    from .globalaxis import excess_over_null, null_degeneracy

    rows = []
    for name, obs in observed.items():
        draws = np.asarray(nulls.get(name, []), dtype=float)
        draws = draws[np.isfinite(draws)]
        row: dict[str, object] = {"signature": name, "observed": obs, "n_null": len(draws)}

        if len(draws) == 0:
            rows.append(row)
            continue

        deg = null_degeneracy(draws)
        est = excess_over_null(obs, draws, n=n_obs or len(draws))
        row.update({
            "excess_z": est.value, "excess_lo": est.lo, "excess_hi": est.hi,
            "beats_null": bool(est.lo > 0) if np.isfinite(est.lo) else False,
            "null_mean": deg.get("r_mean"), "null_sd": deg.get("r_sd"),
            "null_degenerate": deg.get("degenerate"),
            "percentile_SECONDARY_biased_on_tight_null": null_percentile(obs, draws),
        })

        if expr is not None and sig is not None and name in sig.sets:
            row["rho_bar"] = mean_pairwise_correlation(expr, sig.sets[name])
            row["alpha"] = cronbach_alpha(expr, sig.sets[name])
            if null_families and name in null_families:
                fam = null_families[name]
                rhos = [mean_pairwise_correlation(expr, g) for g in fam.sets.values()]
                alphas = [cronbach_alpha(expr, g) for g in fam.sets.values()]
                row["null_rho_bar"] = float(np.nanmean(rhos))
                row["null_alpha"] = float(np.nanmean(alphas))
                if np.isfinite(row["alpha"]) and np.isfinite(row["null_alpha"]):
                    row["reliability_gap"] = float(row["alpha"] - row["null_alpha"])

        rows.append(row)

    out = pd.DataFrame(rows)
    return out.sort_values("excess_z", ascending=False) if "excess_z" in out.columns else out


# ----------------------------------------------------------- reliability ---

def cronbach_alpha(expr: pd.DataFrame, genes: list[str]) -> float:
    """Internal consistency of a gene set, as a scale.

        alpha = k * rho_bar / (1 + (k - 1) * rho_bar)

    where rho_bar is the mean pairwise gene-gene correlation across samples.
    This is the SAME quantity as the expected correlation between two
    independent size-k random sets drawn from a pool with that mean correlation,
    which is why reliability and the null are not separable concerns.

    Why it must be corrected for
    ----------------------------
    Curated signatures are co-expressed modules: high rho_bar, high alpha. Random
    draws of the same size are near-independent: low rho_bar, low alpha. A score
    with higher reliability is measured with less error and is therefore more
    predictable from ANY predictor, entirely independently of whether it carries
    immune-specific biology. Comparing a curated set to a random set without
    disattenuating is biased in favour of the curated set.

    In simulation this bias alone produced a 100% false-positive rate.
    """
    cols = [g for g in genes if g in expr.columns]
    k = len(cols)
    if k < 2:
        return float("nan")
    sub = expr[cols].dropna(axis=0, how="any")
    if len(sub) < 5:
        return float("nan")
    corr = np.corrcoef(sub.to_numpy(dtype=float), rowvar=False)
    off = corr[~np.eye(k, dtype=bool)]
    rho_bar = float(np.nanmean(off))
    if not np.isfinite(rho_bar) or rho_bar <= 0:
        return 0.0
    return float(k * rho_bar / (1 + (k - 1) * rho_bar))


def mean_pairwise_correlation(expr: pd.DataFrame, genes: list[str]) -> float:
    """rho_bar for a gene set — measure this in week 2; it decides the regime.

    Small rho_bar with small k puts you in the ATTENUATION regime (random sets
    are noisy, curated sets win on precision alone). Large k with any rho_bar
    puts you in the DEGENERACY regime (every random draw measures the same axis).
    The correction differs, so measure rather than assume.
    """
    cols = [g for g in genes if g in expr.columns]
    if len(cols) < 2:
        return float("nan")
    sub = expr[cols].dropna(axis=0, how="any")
    if len(sub) < 5:
        return float("nan")
    corr = np.corrcoef(sub.to_numpy(dtype=float), rowvar=False)
    return float(np.nanmean(corr[~np.eye(len(cols), dtype=bool)]))


def set_reliabilities(expr: pd.DataFrame, sig: SignatureSet) -> pd.Series:
    """Cronbach's alpha for every set in a SignatureSet."""
    return pd.Series(
        {name: cronbach_alpha(expr, genes) for name, genes in sig.sets.items()},
        name="alpha",
    )


def split_half_reliability(
    expr: pd.DataFrame,
    genes: list[str],
    *,
    scorer=None,
    residualiser=None,
    n_splits: int = 50,
    seed: int = 0,
) -> dict[str, float]:
    """Reliability of a gene set AS SCORED, via Spearman-Brown split-half.

    Cronbach's alpha is the closed form for mean-z scoring only. It is the wrong
    statistic for ssGSEA (rank-based, non-linear in the genes) and for scores
    that have been residualised on the global axis — and both are used here. This
    is the general form: split the gene set in half, score each half with the
    SAME scorer, residualise each half on the SAME axis, correlate the halves,
    and step the correlation up to full length with Spearman-Brown:

        r_full = 2 * r_half / (1 + r_half)

    Parameters
    ----------
    scorer
        callable(expr, SignatureSet) -> DataFrame. Defaults to `score_mean_z`.
    residualiser
        callable(Series) -> Series applied to each half's score before
        correlating. Pass a global-axis residualiser to get the reliability of
        the quantity the ISI actually uses.
    """
    cols = [g for g in genes if g in expr.columns]
    if len(cols) < 4:
        return {"reliability": float("nan"), "n_genes": len(cols), "n_splits": 0}

    scorer = scorer or score_mean_z
    rng = np.random.default_rng(seed)
    halves = []

    for _ in range(n_splits):
        perm = rng.permutation(cols)
        a, b = list(perm[: len(perm) // 2]), list(perm[len(perm) // 2 :])
        sset = SignatureSet(name="half", sets={"a": a, "b": b})
        sc = scorer(expr, sset)
        if "a" not in sc.columns or "b" not in sc.columns:
            continue
        sa, sb = sc["a"], sc["b"]
        if residualiser is not None:
            sa, sb = residualiser(sa), residualiser(sb)
        ok = np.isfinite(sa) & np.isfinite(sb)
        if ok.sum() < 10:
            continue
        r = float(np.corrcoef(sa[ok], sb[ok])[0, 1])
        if np.isfinite(r):
            halves.append(r)

    if not halves:
        return {"reliability": float("nan"), "n_genes": len(cols), "n_splits": 0}

    r_half = float(np.median(halves))
    r_full = 2 * r_half / (1 + r_half) if r_half > -1 else float("nan")
    return {
        "reliability": float(np.clip(r_full, 0.0, 1.0)),
        "r_half_median": r_half,
        "r_half_sd": float(np.std(halves, ddof=1)) if len(halves) > 1 else 0.0,
        "n_genes": len(cols),
        "n_splits": len(halves),
    }


def alpha_from_scores(scores: np.ndarray, k: int) -> float:
    """Cronbach's alpha computed directly from a mean-z SCORE column. O(n), not O(n*k^2).

    Computing alpha the naive way needs a k x k correlation matrix per gene set,
    which is unaffordable for 1000 null draws x 16 signatures. But when the score
    is the mean of z-scored genes, alpha is recoverable from the score alone.

    Let z be the n x k z-scored submatrix (unit variance, ddof=1) and
    score_i = (1/k) * sum_j z_ij. Then sum_j z_ij = k * score_i, so

        sum(all entries of the correlation matrix)
            = (1/(n-1)) * sum_i (sum_j z_ij)^2
            = k^2 * sum_i score_i^2 / (n-1)

    Subtract the k diagonal ones to get the off-diagonal sum, divide by k(k-1)
    for rho_bar, then alpha = k*rho_bar / (1 + (k-1)*rho_bar).

    THIS IS THE FIX FOR A REAL BUG. The pipeline previously assigned the OBSERVED
    signature's alpha to every null draw, which made the disattenuation cancel
    exactly and do nothing — reintroducing the reliability artefact the whole
    correction exists to remove.

    SCOPE: this closed form is valid for RAW mean-z scores only, because it
    assumes every item has unit variance. For the axis-residualised scores the
    ISI actually disattenuates, use
    `alpha_from_score_and_item_variances` instead — the two differ by a factor of
    five on the null sets.
    """
    s = np.asarray(scores, dtype=float)
    s = s[np.isfinite(s)]
    if len(s) < 5 or k < 2:
        return float("nan")
    total = (k ** 2) * float(np.sum(s ** 2)) / (len(s) - 1)
    off_diag = total - k
    rho_bar = off_diag / (k * (k - 1))
    if not np.isfinite(rho_bar) or rho_bar <= 0:
        return 0.0
    return float(k * rho_bar / (1 + (k - 1) * rho_bar))


def alphas_from_score_frame(scores: pd.DataFrame, sizes: dict[str, int]) -> np.ndarray:
    """Vectorised `alpha_from_scores` over every column of a null score frame."""
    return np.array(
        [alpha_from_scores(scores[c].to_numpy(), sizes.get(c, 0)) for c in scores.columns],
        dtype=float,
    )


def alpha_from_score_and_item_variances(
    score: np.ndarray, item_variances: np.ndarray | pd.Series
) -> float:
    """General-form Cronbach alpha for a composite whose items are NOT unit-variance.

        alpha = k/(k-1) * (1 - sum_j var(item_j) / var(sum_j item_j))

    WHY THIS EXISTS — IT FIXES A REAL METHODOLOGICAL BUG
    ----------------------------------------------------
    `cronbach_alpha` and `alpha_from_scores` both assume the items are z-scored
    genes with unit variance. That holds for a RAW mean-z signature score. It does
    NOT hold once the score has been residualised on the global expression axis,
    which is exactly the quantity the ISI disattenuates. The pipeline previously
    corrected residualised correlations using raw-score reliabilities.

    That mismatch is not cosmetic. Measured on NSCLC
    (SIG_HALLMARK_INTERFERON_GAMMA_RESPONSE, n=944):

        observed alpha    raw 0.9861   residualised 0.9850
        null mean alpha   raw 0.9630   residualised 0.8620
        reliability GAP   raw +0.0232  residualised +0.1230

    The true gap is 5.3x the reported one. The reason is itself a finding: random
    gene sets are internally consistent ONLY because every gene loads on the
    global axis, so removing the axis makes them fall apart; a curated immune
    module has coherent structure beyond the axis and barely loses any.

    Because the composite is the MEAN of k items, sum_j item_j = k * score, so
    the denominator is k^2 * var(score) and only the per-item variances are
    needed — see `globalaxis.residual_item_variances`, which computes them for
    every gene at once.
    """
    v = np.asarray(item_variances, dtype=float)
    v = v[np.isfinite(v)]
    k = len(v)
    s = np.asarray(score, dtype=float)
    s = s[np.isfinite(s)]
    if k < 2 or len(s) < 5:
        return float("nan")

    denominator = (k**2) * float(np.var(s, ddof=1))
    if denominator <= 0:
        return float("nan")
    alpha = k / (k - 1) * (1.0 - float(v.sum()) / denominator)
    # Alpha is a variance ratio and can go negative when items are, on average,
    # negatively correlated. Report 0 rather than a negative "reliability", and
    # never above 1.
    return float(np.clip(alpha, 0.0, 1.0))
