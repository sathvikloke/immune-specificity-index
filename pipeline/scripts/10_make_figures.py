#!/usr/bin/env python3
"""Render the manuscript figures from the frozen result files.

    python scripts/10_make_figures.py --outdir results/figures

Reads `results/nsclc_v3` and `results/pancancer_v3`, and derives no statistic
that the paper also quotes, so a figure cannot disagree with the numbers in the
text. (It named `pancancer_v2` until 2026-09-07 and had done since before v3 was
the reported run; and "computes nothing" was never true -- figure 2 recomputes
two Spearman correlations and figure 3 an MDE. Both claims are corrected here
rather than left as flattering description.)

Figure 1B draws the null as a real distribution. Until 2026-08-19 only the null's
mean and SD reached disk, so that panel showed a mean +/- SD band and said so;
`AuditResult.save()` now persists the per-draw correlations to `null_draws.npz`
and the panel shows what was actually sampled.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as sps  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from aacr27 import signatures as sig_mod  # noqa: E402

NSCLC = ROOT / "results" / "nsclc_v3"
PAN = ROOT / "results" / "pancancer_v3"
GMT = ROOT / "data" / "raw" / "signatures" / "h.all.v2024.1.Hs.symbols.tme.gmt"
EXPR = ROOT / "data" / "interim" / "expr_nsclc.parquet"

# Figure 2's x-axis needs ONE derived quantity from the analysis inputs: how many
# genes of each signature survive expression filtering. That is 16 integers, and
# until 2026-09-07 the only way to obtain them was to read a 41,046-column
# parquet and the GMT -- neither of which is staged into the public snapshot,
# because `pipeline/data/` is excluded there on purpose (the inputs are public
# but large). So `10_make_figures.py` exited 1 at figure 2 for every reader
# reproducing from the deposit, which made the paper's Data Availability
# Statement false a SECOND way, after `null_draws.npz` was fixed the same day.
#
# The fix is the same shape as that one: deposit the derived quantity next to the
# results, where the snapshot's `.csv` rule already carries it. It is 16 rows.
# Recomputation from the inputs stays authoritative when they are present, and
# the two are compared rather than merely preferred -- a deposit nobody checks is
# a deposit that can drift.
GENE_SIZES = ROOT / "results" / "gene_set_sizes.csv"
N_SIGNATURES = 16

# Okabe-Ito: colourblind-safe.
C_PAN, C_NSCLC, C_NULL = "#0072B2", "#D55E00", "#999999"
CAT_COLOUR = {"immune": "#0072B2", "proliferation": "#D55E00",
              "development": "#009E73", "signaling": "#CC79A7",
              "pathway": "#E69F00"}

# Liberzon et al., Cell Systems 2015, Table 1 process categories. COAGULATION and
# COMPLEMENT are curated as IMMUNE — this is not our grouping, and it is the
# reason the immune-vs-other test is p=0.016 rather than the p=0.00008 a
# hand-made split produced.
CAT = {
    "ALLOGRAFT_REJECTION": "immune", "COAGULATION": "immune", "COMPLEMENT": "immune",
    "INFLAMMATORY_RESPONSE": "immune", "INTERFERON_ALPHA_RESPONSE": "immune",
    "INTERFERON_GAMMA_RESPONSE": "immune",
    # IL6_JAK_STAT3_SIGNALING is left OUT of the immune group. MSigDB's own gene
    # set page states no process category for it, so the assignment could not be
    # verified. Excluding it is the conservative choice: it has a negative excess,
    # so counting it as immune would only strengthen the immune-vs-other contrast
    # (p=0.016 -> 0.008).
    "IL6_JAK_STAT3_SIGNALING": "signaling",
    "IL2_STAT5_SIGNALING": "signaling", "TNFA_SIGNALING_VIA_NFKB": "signaling",
    "TGF_BETA_SIGNALING": "signaling", "E2F_TARGETS": "proliferation",
    "G2M_CHECKPOINT": "proliferation", "MYC_TARGETS_V1": "proliferation",
    "ANGIOGENESIS": "development", "EPITHELIAL_MESENCHYMAL_TRANSITION": "development",
    "HYPOXIA": "pathway",
}

# Extra output formats, off by default so the routine render stays fast.
#
# PNG and PDF are always written. The two opt-in formats exist for two different
# destinations, and neither is a matter of taste:
#
#   --tiff  AACR's accepted graphics formats for REVISED manuscripts are EPS,
#           TIFF, AI, PSD, PNG and PS. *PDF is not on that list*, and Editorial
#           Policies separately asks that images be captured in "an uncompressed
#           format such as TIFF". 600 dpi is the line-art resolution the
#           submission guides cite. LZW is lossless, so "uncompressed" in the
#           sense that matters -- no generation loss -- at a tenth the size of
#           raw (~34 MB per figure at this size).
#   --svg   For the A0 poster, where these figures are enlarged 1.5-1.7x beyond
#           their native 300 dpi size (figure3 prints 330 mm wide against a
#           190 mm native width, i.e. an effective 172 dpi). SVG is
#           resolution-independent, so the enlargement costs nothing.
#
# Neither is written by default, so neither lands in the repository unless asked.
EXTRA_FORMATS: list[str] = []


def save(fig, outdir: Path, stem: str) -> None:
    """Write one figure in every requested format, then close it."""
    fig.savefig(outdir / f"{stem}.png")
    fig.savefig(outdir / f"{stem}.pdf")
    if "tiff" in EXTRA_FORMATS:
        fig.savefig(outdir / f"{stem}.tiff", dpi=600,
                    pil_kwargs={"compression": "tiff_lzw"})
    if "svg" in EXTRA_FORMATS:
        fig.savefig(outdir / f"{stem}.svg")
    plt.close(fig)


plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.4,
})


def short(s: pd.Series) -> pd.Series:
    return s.str.replace("SIG_HALLMARK_", "", regex=False).str.replace("_", " ")


def load(path: Path, name: str) -> pd.DataFrame:
    d = pd.read_csv(path / f"{name}.csv")
    if "signature" in d.columns:
        d["sig"] = d["signature"].str.replace("SIG_HALLMARK_", "", regex=False)
        d["label"] = short(d["signature"])
    return d


def _rel(p: Path) -> str:
    """Display path, relative to the pipeline root when it is under it."""
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def _sizes_from_inputs() -> dict[str, int] | None:
    """Recompute the panel sizes from the analysis inputs, or None if absent.

    Returns None rather than raising when `data/` is missing, because that is
    the normal state of the public snapshot, not an error.
    """
    if not (EXPR.exists() and GMT.exists()):
        return None
    expr_cols = set(pd.read_parquet(EXPR).columns)
    sets = sig_mod.SignatureSet.from_gmt(GMT).sets
    # GMT keys are HALLMARK_*; the result files key on the bare name. Without
    # stripping the prefix the join silently yields all-NaN and the Spearman
    # correlation in figure 2 comes out NaN rather than raising.
    return {k.replace("HALLMARK_", ""): len([g for g in v if g in expr_cols])
            for k, v in sets.items()}


def _sizes_from_deposit() -> dict[str, int] | None:
    """Read the deposited panel sizes, or None if they are not deposited.

    An INCOMPLETE deposit is rejected here rather than tolerated. A file that
    exists but covers fewer than the 16 signatures would map to NaN in figure 2
    and yield a NaN Spearman rather than an error -- the same shape of silent
    hole that let an inadequate `null_draws.npz` render panel B with no null at
    all. Coverage is not presence.
    """
    if not GENE_SIZES.exists():
        return None
    d = pd.read_csv(GENE_SIZES)
    if list(d.columns) != ["signature", "n_genes"] or len(d) != N_SIGNATURES:
        raise SystemExit(
            f"FATAL: {_rel(GENE_SIZES)} is malformed -- expected {N_SIGNATURES} "
            f"rows of (signature, n_genes), found {len(d)} row(s) of "
            f"{list(d.columns)}. Regenerate it with --write-gene-sizes on a "
            f"machine that has data/, or delete it to fall back to recomputation.")
    return dict(zip(d["signature"], d["n_genes"].astype(int)))


def gene_set_sizes() -> dict[str, int]:
    """The 16 post-filtering panel sizes figure 2 plots against.

    Recomputation wins when `data/` is present; the deposit serves the snapshot.
    When both exist they are COMPARED, so a stale deposit fails loudly on the
    one machine able to notice.
    """
    computed, deposited = _sizes_from_inputs(), _sizes_from_deposit()
    if computed is not None and deposited is not None and computed != deposited:
        diff = {k: (computed.get(k), deposited.get(k))
                for k in set(computed) | set(deposited)
                if computed.get(k) != deposited.get(k)}
        raise SystemExit(
            f"FATAL: {_rel(GENE_SIZES)} disagrees with a fresh recomputation "
            f"from data/. computed vs deposited: {diff}\n"
            f"       The deposit is what the public snapshot ships, so this "
            f"must be resolved, not ignored. Re-run with --write-gene-sizes.")
    if computed is not None:
        return computed
    if deposited is not None:
        print(f"  figure2: panel sizes read from {_rel(GENE_SIZES)} "
              f"(data/ is absent -- this is the deposited-snapshot path)")
        return deposited
    raise SystemExit(
        f"FATAL: figure 2 needs the per-signature panel sizes and neither route "
        f"is available.\n"
        f"       recomputation needs {_rel(EXPR)} and {_rel(GMT)};\n"
        f"       the deposit would be at {_rel(GENE_SIZES)}.\n"
        f"       Obtain the inputs (see README) or the deposited results.")


# ------------------------------------------------------------- schematic ---

def figure0(outdir: Path) -> None:
    """The argument, as a diagram. Poster panel 1; not a manuscript figure.

    Draws no data. It exists because the central point — two hypotheses predict
    the SAME observed correlation, so the usual control cannot separate them — is
    the one thing a passer-by must grasp in five seconds, and prose does not do
    it at poster distance.
    """
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    fig, ax = plt.subplots(figsize=(7.8, 4.3))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.2)
    ax.axis("off")

    def box(x, y, w, h, text, fc, ec, size=7.5, weight="normal"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.10",
                                    facecolor=fc, edgecolor=ec, linewidth=1.2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=size, fontweight=weight)

    def arrow(x1, y1, x2, y2, colour="#444444", ls="-"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=11, color=colour,
                                     linewidth=1.1, linestyle=ls,
                                     shrinkA=3, shrinkB=3))

    ax.text(5.0, 5.9, "Why the published control is not enough",
            ha="center", fontsize=10.5, fontweight="bold")

    box(0.15, 4.35, 2.8, 1.05,
        "H1  The image carries\nimmune-specific signal", "#DCE9F5", C_PAN, 8)
    box(0.15, 2.30, 2.8, 1.05,
        "H2  The image reads\ncomposition faithfully", "#FBE3D5", C_NSCLC, 8)

    box(3.60, 3.32, 2.5, 1.05,
        "Identical observed\ncorrelation with a\ncurated signature",
        "#EFEFEF", "#666666", 8, "bold")
    arrow(2.95, 4.75, 3.55, 4.20)
    arrow(2.95, 2.90, 3.55, 3.50)

    box(6.75, 3.32, 3.10, 1.05,
        "A size-matched random-set\nnull CANNOT separate them\n"
        "-- both load on the same axis", "#FFF3CD", "#B8860B", 7.5)
    arrow(6.12, 3.85, 6.70, 3.85)

    # Caption sits in the clear band BETWEEN the arrowheads and the boxes; at
    # y=2.05 it ran through the H2 box and both dashed arrows.
    ax.text(5.0, 1.62, "Two corrections make the hypotheses distinguishable",
            ha="center", fontsize=8.8, fontweight="bold")

    box(0.55, 0.30, 3.9, 1.05,
        "1.  Residualise the SIGNATURE\non the global expression axis",
        "#E3F2E8", "#009E73", 7.8)
    box(5.10, 0.30, 4.3, 1.05,
        "2.  Disattenuate using the\nRESIDUALISED score's reliability",
        "#E3F2E8", "#009E73", 7.8)
    arrow(4.30, 3.25, 3.10, 2.05, ls="--")
    arrow(5.45, 3.25, 6.75, 2.05, ls="--")

    save(fig, outdir, "figure0_schematic")
    print("  figure0: schematic (poster panel 1, no data)")


# ------------------------------------------------------------------ fig 1 ---

def figure1(outdir: Path) -> None:
    """Per-signature immune-specific excess, both cohorts, plus observed vs null."""
    pan, nsc = load(PAN, "immune_excess"), load(NSCLC, "immune_excess")
    order = pan.sort_values("excess_z")["sig"].tolist()
    pan = pan.set_index("sig").loc[order].reset_index()
    nsc = nsc.set_index("sig").loc[order].reset_index()
    y = np.arange(len(order))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 4.4),
                                   gridspec_kw={"width_ratios": [1.15, 1]})

    ax1.axvline(0, color="k", lw=0.8, zorder=1)
    for d, c, off, lab in ((pan, C_PAN, 0.16, "Pan-TCGA (n=7,168)"),
                           (nsc, C_NSCLC, -0.16, "NSCLC (n=944)")):
        ax1.hlines(y + off, d["excess_lo"], d["excess_hi"], color=c, lw=1.6, alpha=.85)
        ax1.plot(d["excess_z"], y + off, "o", ms=3.6, color=c, label=lab, zorder=3)
    ax1.set_yticks(y)
    ax1.set_yticklabels([o.replace("_", " ") for o in order])
    ax1.set_xlabel("Immune-specific excess (Fisher $z$)")
    ax1.set_title("A  Per-signature ISI, disattenuated", loc="left", fontweight="bold")
    ax1.legend(loc="lower right", frameon=False)

    # Observed vs null. The null is the ACTUAL 1,000 draws per signature, not a
    # summary band — see the module docstring.
    draws_path = PAN / "null_draws.npz"
    _fig1b_null = "violins"
    data = []
    if draws_path.exists():
        z = np.load(draws_path)
        key = {k.replace("SIG_HALLMARK_", ""): k for k in z.files}
        data = [z[key[s]] for s in order if s in key]
        if len(data) != len(order):
            # THE SECOND SILENT HOLE, closed 2026-09-07. Previously an npz that
            # existed but did not cover all 16 signatures fell through BOTH
            # arms: no violins, and no band either, so panel B rendered with no
            # null at all and nothing said so. `else` guarded only a missing
            # file, not an inadequate one. Coverage is not the same as presence.
            print(f"\n  !! WARNING: {draws_path} covers {len(data)} of "
                  f"{len(order)} signatures.\n"
                  "     Falling back to the mean +/- SD band for ALL of them, "
                  "rather than\n     drawing a partial null that would look "
                  "complete.\n")
            _fig1b_null = "band"
    else:
        _fig1b_null = "band"

    if _fig1b_null == "violins":
        vp = ax2.violinplot(data, positions=y, vert=False, widths=0.85,
                            showextrema=False, showmedians=True)
        for b in vp["bodies"]:
            b.set_facecolor(C_NULL)
            b.set_alpha(0.55)
            b.set_edgecolor("none")
        vp["cmedians"].set_color("#444444")
        vp["cmedians"].set_linewidth(0.9)
        ax2.plot([], [], "s", color=C_NULL, alpha=.55, ms=6,
                 label="Random-set null (1,000 draws)")
    else:
        # LOUD, and it must stay loud. This branch draws a DIFFERENT FIGURE
        # from the published one -- a mean +/- SD band in place of the real
        # 1,000 draws -- and until 2026-09-07 it did so in silence while
        # `main()` unconditionally printed that panel B showed the real
        # distribution. Anyone reproducing from the public snapshot hit exactly
        # this path, because `null_draws.npz` was excluded from the snapshot by
        # suffix, and was told in the script's own output that they had got the
        # violins. The file is now staged (it is 124 KB), so this branch should
        # be unreachable for a reader who cloned the snapshot -- if you are
        # seeing it, something is missing that the paper says is deposited.
        print("\n  !! WARNING: figure 1B fell back to a mean +/- SD BAND.\n"
              f"     {draws_path} is absent or does not cover all "
              f"{len(order)} signatures,\n"
              "     so the 1,000 per-signature null draws could not be drawn.\n"
              "     THIS IS NOT THE PUBLISHED PANEL B, which shows the real\n"
              "     distribution as violins. Re-run the pan-cancer analysis to\n"
              "     regenerate null_draws.npz, or obtain it with the deposited\n"
              "     results.\n")
        ax2.hlines(y, pan["null_mean_r"] - pan["null_sd_r"],
                   pan["null_mean_r"] + pan["null_sd_r"], color=C_NULL, lw=4,
                   alpha=.55, label="Random-set null (mean $\\pm$ SD)")
        # Annotate the ARTEFACT ITSELF, not just the terminal. A figure that
        # leaves the machine carrying a claim it cannot support is the failure
        # mode; a warning scrolled past in a log is not a control.
        ax2.text(0.98, 0.02,
                 "null shown as mean $\\pm$ SD\n(per-draw nulls unavailable)",
                 transform=ax2.transAxes, ha="right", va="bottom",
                 fontsize=5.5, color="#B00020")
    ax2.plot(pan["r_residual_refit"], y, "o", ms=4, color=C_PAN,
             label="Curated signature", zorder=3)
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlabel("Axis-residualised $r$ with image prediction")
    ax2.set_title("B  Pan-TCGA: curated vs random", loc="left", fontweight="bold")
    ax2.legend(loc="upper left", frameon=False, fontsize=6.5)
    ax2.set_xlim(-0.02, 0.58)

    save(fig, outdir, "figure1_isi")
    print(f"  figure1: 16 signatures, pan excess "
          f"{pan.excess_z.min():.3f}-{pan.excess_z.max():.3f} "
          f"[panel B null: {_fig1b_null}]")
    return _fig1b_null


# ------------------------------------------------------------------ fig 2 ---

def figure2(outdir: Path) -> None:
    """Reliability gap vs panel size, and why the raw-score alpha misleads."""
    sizes = gene_set_sizes()
    pan, nsc = load(PAN, "immune_excess"), load(NSCLC, "immune_excess")
    for name, d in (("pancancer", pan), ("nsclc", nsc)):
        d["k"] = d["sig"].map(sizes)
        # ASSERT THE JOIN. The comment on the HALLMARK_ prefix has warned since
        # 2026-08-19 that an unmatched key yields all-NaN and a NaN Spearman
        # "rather than raising" -- and nothing raised. A documented hazard with
        # no guard is a hazard, so this is now the guard.
        missing = sorted(d.loc[d["k"].isna(), "sig"])
        if missing:
            raise SystemExit(
                f"FATAL: figure 2's panel-size join left {len(missing)} of "
                f"{len(d)} {name} signatures unmatched: {missing}\n"
                f"       Spearman would have been NaN and the panel would have "
                f"rendered anyway. Panel sizes cover: {sorted(sizes)}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 3.4))

    for d, c, m, lab in ((pan, C_PAN, "o", "Pan-TCGA"), (nsc, C_NSCLC, "s", "NSCLC")):
        ax1.scatter(d["k"], d["reliability_gap"], s=26, color=c, marker=m,
                    label=f"{lab}, residualised", zorder=3, edgecolor="none")
        ax1.scatter(d["k"], d["alpha_observed_raw"] - d["alpha_null_mean_raw"], s=22,
                    facecolor="none", edgecolor=c, marker=m, lw=.8,
                    label=f"{lab}, raw score")
    rho, p = sps.spearmanr(pan["k"], pan["reliability_gap"])
    rho_n, p_n = sps.spearmanr(nsc["k"], nsc["reliability_gap"])
    ax1.set_xlabel("Genes in signature (after filtering)")
    ax1.set_ylabel("Curated $-$ random Cronbach's $\\alpha$")
    ax1.set_title("A  Reliability gap scales with panel size",
                  loc="left", fontweight="bold")
    # The statistic goes INSIDE the axes: as a second title line it overprinted
    # panel B's title, and the k=36 point sat underneath the legend.
    ax1.text(0.97, 0.93,
             f"Spearman $\\rho$ = {rho:.2f}, p = {p:.3f}  (pan-TCGA)\n"
             f"$\\rho$ = {rho_n:.2f}, p = {p_n:.3f}  (NSCLC)",
             transform=ax1.transAxes, ha="right", va="top", fontsize=6.5)
    ax1.legend(frameon=False, ncol=2, fontsize=6, loc="upper right",
               bbox_to_anchor=(1.0, 0.80))
    ax1.axhline(0, color="k", lw=0.6)
    ax1.set_ylim(-0.04, 0.50)
    ax1.set_xlim(0, 225)

    # Panel B: where the alpha actually goes when the axis is removed.
    w = 0.35
    x = np.arange(2)
    for i, (d, lab, c) in enumerate(((pan, "Pan-TCGA", C_PAN), (nsc, "NSCLC", C_NSCLC))):
        vals = [d["alpha_null_mean_raw"].median(), d["alpha_null_mean"].median()]
        ax2.bar(x + (i - .5) * w, vals, w, color=c, alpha=.85, label=f"{lab}, random sets")
        obs = [d["alpha_observed_raw"].median(), d["alpha_observed"].median()]
        ax2.plot(x + (i - .5) * w, obs, "k_", ms=14, mew=1.6,
                 label="Curated signature" if i == 0 else None)
    ax2.set_xticks(x)
    ax2.set_xticklabels(["Raw score", "Axis-residualised"])
    ax2.set_ylabel("Median Cronbach's $\\alpha$")
    ax2.set_ylim(0, 1.05)
    ax2.set_title("B  Random sets are consistent only because\nthey load on the axis",
                  loc="left", fontweight="bold")
    ax2.legend(frameon=False, loc="lower left", fontsize=6, ncol=1)
    for i, (d, c) in enumerate(((pan, C_PAN), (nsc, C_NSCLC))):
        for j, v in enumerate([d["alpha_null_mean_raw"].median(), d["alpha_null_mean"].median()]):
            ax2.text(j + (i - .5) * w, v + .02, f"{v:.2f}", ha="center",
                     fontsize=6, color=c, fontweight="bold")

    save(fig, outdir, "figure2_reliability")
    print(f"  figure2: Spearman rho={rho:.3f} p={p:.4f}; "
          f"pan gap raw {(pan.alpha_observed_raw-pan.alpha_null_mean_raw).mean():+.3f} "
          f"-> residualised {pan.reliability_gap.mean():+.3f}")


# ------------------------------------------------------------------ fig 3 ---

def figure3(outdir: Path) -> None:
    """Outcome excess over null, pan-TCGA, by MSigDB process category."""
    o = load(PAN, "outcome_arm")
    o = o[o["adjustment"] == "axis_residualised"].copy()
    o["cat"] = o["sig"].map(CAT)
    o = o.sort_values("excess_z").reset_index(drop=True)
    y = np.arange(len(o))

    crit = sps.norm.ppf(0.975)
    se = (o["excess_hi"] - o["excess_lo"]) / (2 * crit)
    mde = float(((crit + sps.norm.ppf(0.80)) * se).median())

    n = load(NSCLC, "outcome_arm")
    n = n[n["adjustment"] == "axis_residualised"].set_index("sig")
    nse = (n["excess_hi"] - n["excess_lo"]) / (2 * crit)
    nmde = float(((crit + sps.norm.ppf(0.80)) * nse).median())

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.axvspan(-mde, mde, color="#000000", alpha=.055, zorder=0,
               label=f"Pan-TCGA MDE $\\pm${mde:.3f} (80% power)")
    ax.axvline(0, color="k", lw=0.8, zorder=1)

    for i, r in o.iterrows():
        c = CAT_COLOUR[r["cat"]]
        ax.hlines(i, r["excess_lo"], r["excess_hi"], color=c, lw=1.6, alpha=.85)
        ax.plot(r["excess_z"], i, "o", ms=4.5, color=c, zorder=3,
                markeredgecolor="k" if r["beats_null"] else "none", mew=.7)
    ax.plot(n.reindex(o["sig"])["excess_z"].to_numpy(), y, "x", ms=4,
            color="#666666", label=f"NSCLC (MDE $\\pm${nmde:.3f})", zorder=2)

    ax.set_yticks(y)
    ax.set_yticklabels(o["sig"].str.replace("_", " "))
    ax.set_xlabel("Outcome excess over random-set null (Fisher $z$), within cancer type")
    ax.set_title("Prognostic advantage over random gene sets is proliferative and stromal,\n"
                 "not immune (PFI, pan-TCGA n=7,168)", loc="left", fontweight="bold")

    handles = [plt.Line2D([], [], marker="o", ls="", color=CAT_COLOUR[k], label=k)
               for k in ["immune", "proliferation", "development", "pathway", "signaling"]]
    handles.append(plt.Line2D([], [], marker="o", ls="", color="w",
                              markeredgecolor="k", label="beats null"))
    leg1 = ax.legend(handles=handles, frameon=False, loc="lower right", fontsize=6.5,
                     title="MSigDB process category", title_fontsize=6.5)
    ax.add_artist(leg1)
    ax.legend(frameon=False, loc="upper left", fontsize=6.5)

    save(fig, outdir, "figure3_outcome")
    imm = o[o["cat"] == "immune"]
    print(f"  figure3: {int(o.beats_null.sum())}/16 beat null; "
          f"immune {int(imm.beats_null.sum())}/{len(imm)}; MDE {mde:.4f}")


# ------------------------------------------------------------------ fig 4 ---

def figure4(outdir: Path) -> None:
    """The control that decides how much of this to believe: site recoverability."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0))

    for path, c, lab in ((PAN, C_PAN, "Pan-TCGA"), (NSCLC, C_NSCLC, "NSCLC")):
        sc = pd.read_csv(path / "site_control.csv")
        ax1.hist(sc["auroc"], bins=np.linspace(0.5, 1.0, 26), color=c, alpha=.7,
                 label=f"{lab} ({len(sc)} sites, median {sc['auroc'].median():.3f})")
    ax1.axvline(0.5, color="k", lw=0.8, ls="--")
    ax1.set_xlabel("One-vs-rest AUROC for tissue source site")
    ax1.set_ylabel("Sites")
    ax1.set_title("A  Site is almost perfectly recoverable\nfrom the embeddings",
                  loc="left", fontweight="bold")
    ax1.legend(frameon=False, fontsize=6.5, loc="upper left")

    ps = pd.read_csv(PAN / "per_signature.csv")
    ns = pd.read_csv(NSCLC / "per_signature.csv")
    x = np.arange(2)
    w = 0.35
    for i, (d, lab, c) in enumerate(((ps, "Pan-TCGA", C_PAN), (ns, "NSCLC", C_NSCLC))):
        ax2.bar(x + (i - .5) * w,
                [d["r_covariates"].median(), d["r_preserved_site"].median()],
                w, color=c, alpha=.85, label=lab)
    ax2.set_xticks(x)
    ax2.set_xticklabels(["Covariates only\n(site+type+purity+stage)", "Image embedding"])
    ax2.set_ylabel("Median $r$ with signature")
    ax2.set_title("B  Pan-cancer the image adds nothing\nover covariates ($\\Delta r=-0.000$)",
                  loc="left", fontweight="bold")
    ax2.legend(frameon=False, fontsize=6.5)

    save(fig, outdir, "figure4_controls")
    print(f"  figure4: pan emb-over-cov {ps.embedding_over_covariates.median():+.4f}, "
          f"nsclc {ns.embedding_over_covariates.median():+.4f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "figures")
    ap.add_argument("--tiff", action="store_true",
                    help="also write 600 dpi LZW TIFF -- the format AACR accepts "
                         "for revised manuscripts (PDF is not on their list)")
    ap.add_argument("--svg", action="store_true",
                    help="also write SVG, for the A0 poster where these figures "
                         "are enlarged past their native 300 dpi size")
    ap.add_argument("--write-gene-sizes", action="store_true",
                    help=f"recompute the 16 panel sizes from data/ and write "
                         f"{_rel(GENE_SIZES)}, the deposit that lets figure 2 "
                         f"render from the public snapshot. Renders nothing.")
    args = ap.parse_args()

    if args.write_gene_sizes:
        sizes = _sizes_from_inputs()
        if sizes is None:
            raise SystemExit(
                f"FATAL: --write-gene-sizes needs the analysis inputs, and "
                f"{_rel(EXPR)} or {_rel(GMT)} is missing. The deposit can only "
                f"be regenerated on a machine that has data/.")
        d = (pd.DataFrame(sorted(sizes.items()), columns=["signature", "n_genes"])
             .sort_values("signature", kind="stable"))
        GENE_SIZES.parent.mkdir(parents=True, exist_ok=True)
        d.to_csv(GENE_SIZES, index=False)
        print(f"wrote {_rel(GENE_SIZES)}: {len(d)} signature(s), "
              f"{d['n_genes'].min()}-{d['n_genes'].max()} genes after filtering")
        return 0

    args.outdir.mkdir(parents=True, exist_ok=True)

    EXTRA_FORMATS.extend(f for f in ("tiff", "svg") if getattr(args, f))
    extra = f" (+ {', '.join(EXTRA_FORMATS)})" if EXTRA_FORMATS else ""

    print(f"rendering into {args.outdir}: png, pdf{extra}")
    figure0(args.outdir)
    fig1b_null = figure1(args.outdir)
    figure2(args.outdir)
    figure3(args.outdir)
    figure4(args.outdir)
    # CONDITIONAL, since 2026-09-07. This NOTE used to print unconditionally,
    # which made it an affirmative false statement in exactly the case where it
    # mattered: a reader whose `null_draws.npz` was missing got the band AND was
    # told they had got the real distribution. It was not a stale caveat, it was
    # a claim the script could not check -- and the one population guaranteed to
    # hit it was anyone reproducing from the public snapshot, which excluded the
    # file by suffix until it was staged.
    if fig1b_null == "violins":
        print("\nNOTE: figure 1B draws the null as a REAL distribution -- the "
              "1,000 per-signature draws persisted to null_draws.npz by "
              "AuditResult.save(), not a mean +/- SD band. The band was what "
              "this script drew before 2026-08-19; that is obsolete.")
    else:
        print("\nNOTE: figure 1B drew the null as a mean +/- SD BAND, NOT the "
              "real distribution. See the warning above. The rendered figure "
              "does not match the published panel B and is annotated in the "
              "artefact itself to say so.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
