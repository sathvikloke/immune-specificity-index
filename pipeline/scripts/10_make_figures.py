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
#   --svg   Written for the A0 poster until 2026-09-24, which enlarged these
#           figures 1.5-1.7x beyond their native 300 dpi size; SVG is
#           resolution-independent, so the enlargement cost nothing. The
#           poster now embeds its own renders (`--poster`, below), and the
#           committed journal SVGs remain as the render-drift record.
#
# Neither is written by default, so neither lands in the repository unless asked.
EXTRA_FORMATS: list[str] = []


def save_tiff_rgb(fig, path: Path, dpi: int = 600) -> None:
    """600 dpi LZW TIFF in RGB, not matplotlib's RGBA (ledger F11.2, 2026-09-17).

    Print production expects RGB or CMYK line art; an alpha channel is at best
    ignored and at worst composited against black. The figures are drawn on an
    opaque white face, so dropping alpha loses nothing. That is checked, and a
    figure with any transparent pixel is composited onto white instead.
    """
    import io

    from PIL import Image

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    buf.seek(0)
    rgba = Image.open(buf).convert("RGBA")
    if rgba.getextrema()[3][0] < 255:
        rgb = Image.new("RGB", rgba.size, "white")
        rgb.paste(rgba, mask=rgba.getchannel("A"))
    else:
        rgb = rgba.convert("RGB")
    rgb.save(path, compression="tiff_lzw", dpi=(dpi, dpi))


def save(fig, outdir: Path, stem: str) -> None:
    """Measure one figure's print size, write it in every requested format, close it."""
    print_size(fig, stem)
    if POSTER:
        # SVG only -- the poster embeds nothing else -- and UNCROPPED, so its
        # size is exactly the canvas and its aspect exactly the box's. NOT
        # `bbox_inches=None`: that means "use rcParams", which say "tight", and
        # the first poster renders came out cropped and padded -- 535 x 306 px
        # for a 515 x 287 px canvas -- so the browser scaled them 2.11x, not
        # 2.25x, and their 8 pt text printed at 16.9 pt. Measured in the
        # browser, not assumed from the canvas.
        with plt.rc_context({"savefig.bbox": "standard"}):
            fig.savefig(outdir / f"{stem}.svg")
        plt.close(fig)
        return
    fig.savefig(outdir / f"{stem}.png")
    fig.savefig(outdir / f"{stem}.pdf")
    if "tiff" in EXTRA_FORMATS:
        save_tiff_rgb(fig, outdir / f"{stem}.tiff")
    if "svg" in EXTRA_FORMATS:
        fig.savefig(outdir / f"{stem}.svg")
    plt.close(fig)


# JOURNAL PRINT SIZE (ledger H87, closing F11.6; 2026-09-24). Every figure is
# drawn at the size it prints: saved no wider than JOURNAL_WIDTH_IN (a double
# column) and no visible text under JOURNAL_MIN_PT once scaled to that width,
# MEASURED from the drawn artists by `print_size()` for every figure on every
# render. Until 2026-09-24 the type was 7 pt ticks and legends on canvases up to
# 7.8 in wide, plus explicit 6-6.5 pt legends and annotations: at 6.75 in the
# smallest text printed at 5.61 / 5.84 / 5.87 / 6.96 pt in figures 1-4 (F11.6
# measured the 7 pt ticks alone, 6.1-7.4 pt).
JOURNAL_WIDTH_IN = 6.75
JOURNAL_MIN_PT = 8.0
# The canvas width that saves at JOURNAL_WIDTH_IN once `savefig.pad_inches`
# (0.1 in, matplotlib's default, left alone) is added on both sides -- provided
# nothing is drawn outside the canvas, which constrained layout ensures and
# `print_size()` measures rather than assumes.
FIG_W = JOURNAL_WIDTH_IN - 2 * 0.1
JOURNAL_FIGURES = ("figure1_isi", "figure2_reliability", "figure3_outcome",
                   "figure4_controls")
PRINT_SIZES: dict[str, dict] = {}

# POSTER GEOMETRY (ledger H88, TODO-PRINT-6; 2026-09-24). A poster figure is
# drawn on a canvas with EXACTLY the aspect of its A0 box and 1/POSTER_SCALE its
# size, with the journal type, and saved without cropping -- so the poster's
# `object-fit: contain` enlarges it by exactly POSTER_SCALE and every glyph,
# line and marker prints POSTER_SCALE times its journal size: the journal's
# 8 pt floor becomes POSTER_MIN_PT. Until 2026-09-24 the poster embedded the
# journal figures in boxes of another shape, and their smallest text printed at
# 7.2-13.5 pt on A0 (a 7 pt tick measure; the legends were smaller still).
POSTER_MIN_PT = 18.0
POSTER_SCALE = POSTER_MIN_PT / JOURNAL_MIN_PT
# `.figbox` width and height in CSS px (1 px = 25.4/96 mm on the A0 sheet),
# measured in a browser with getBoundingClientRect after two layout changes made
# FOR these figures on 2026-09-24: the grid's rows went from 1.241fr / 1fr to
# 1.10fr / 1fr, and panel 4's pan-cancer partition paragraph moved to panel 6.
# Figures 1 and 3 each list 16 signatures one per row, and at an 18 pt floor
# that needs about 150 mm of box height apiece; the old layout gave them 182 mm
# and 108 mm. Each value is the measured size ROUNDED DOWN (1160.10 x 644.64
# -> 1160 x 644, and so on): a canvas drawn for a box a fraction of a pixel
# larger than the real one is scaled a fraction under POSTER_SCALE, and prints
# its 8 pt text at 17.99 pt. `test_a0_print_legibility` holds the same table;
# if the poster's layout moves, re-measure both.
POSTER_BOXES_PX = {
    "figure0_schematic": (1160, 644),
    "figure3_outcome": (1587, 616),
    "figure2_reliability": (1160, 430),
    "figure1_isi": (1594, 601),
}
POSTER = False          # set by --poster


def canvas(stem: str, journal: tuple[float, float]) -> tuple[float, float]:
    """The figure size to draw at: the journal size, or the poster box / POSTER_SCALE."""
    if not POSTER:
        return journal
    w_px, h_px = POSTER_BOXES_PX[stem]
    return (w_px / 96 / POSTER_SCALE, h_px / 96 / POSTER_SCALE)


def by_mode(journal, poster):
    """A layout choice that differs between the journal and the poster canvas."""
    return poster if POSTER else journal


def print_size(fig, stem: str) -> dict:
    """Measure one figure as it will print: saved width, and its smallest text.

    The saved width is the tight bounding box plus `savefig.pad_inches`, as
    `savefig(bbox_inches="tight")` writes it. The printed size of a glyph is its
    point size times JOURNAL_WIDTH_IN / saved width -- the scale at which the
    figure fills a double column. Every visible, non-empty Text artist with a
    drawn extent is counted: tick labels, legend entries and titles,
    annotations, panel titles.
    """
    import matplotlib.text as mtext

    # Count what a draw actually DRAWS. `fig.findobj(Text)` also returns the
    # labels of ticks outside the view limits, which are never drawn (figure 3's
    # x axis carries an undrawn 0.125), so the draw itself is intercepted.
    drawn: dict[int, mtext.Text] = {}
    original = mtext.Text.draw

    def _record(self, renderer):
        if self.get_visible() and self.get_text().strip():
            drawn[id(self)] = self
        return original(self, renderer)

    mtext.Text.draw = _record
    try:
        fig.canvas.draw()
    finally:
        mtext.Text.draw = original
    renderer = fig.canvas.get_renderer()
    if POSTER:          # saved uncropped: the SVG is exactly the canvas
        width, height = fig.get_size_inches()
    else:
        pad = plt.rcParams["savefig.pad_inches"]
        bbox = fig.get_tightbbox(renderer)
        width, height = bbox.width + 2 * pad, bbox.height + 2 * pad
    texts = [t for t in drawn.values() if t.get_window_extent(renderer).width > 0]
    if not texts:
        raise SystemExit(f"FATAL: {stem} has no visible text to measure")
    sizes = [t.get_fontsize() for t in texts]
    # A glyph past the canvas edge is either cropped or widens the saved figure
    # past the column, depending on the renderer's bbox; count it either way. A
    # half-pixel slack absorbs antialiasing.
    canvas = fig.bbox
    outside = [t.get_text()[:40] for t in texts
               if (e := t.get_window_extent(renderer)).x0 < canvas.x0 - 0.5
               or e.x1 > canvas.x1 + 0.5 or e.y0 < canvas.y0 - 0.5 or e.y1 > canvas.y1 + 0.5]
    # Two drawn texts whose boxes intersect -- tick labels packed tighter than
    # their own height, titles running into each other -- are unreadable
    # whatever their size. A 1-pixel slack absorbs boxes that merely touch.
    ext = [(t, t.get_window_extent(renderer)) for t in texts]
    overlaps = []
    for i, (ta, ea) in enumerate(ext):
        for tb, eb in ext[i + 1:]:
            if (min(ea.x1, eb.x1) - max(ea.x0, eb.x0) > 1
                    and min(ea.y1, eb.y1) - max(ea.y0, eb.y0) > 1):
                overlaps.append((ta.get_text()[:30], tb.get_text()[:30]))
    smallest = min(sizes)
    if POSTER:          # object-fit: contain, and the aspect is the box's
        w_px, h_px = POSTER_BOXES_PX[stem]
        scale = min(w_px / 96 / width, h_px / 96 / height)
    else:
        scale = JOURNAL_WIDTH_IN / width
    rec = {"width_in": width, "height_in": height, "smallest_pt": smallest,
           "n_text": len(sizes), "outside_canvas": outside, "overlaps": overlaps,
           "printed_smallest_pt": smallest * scale}
    PRINT_SIZES[stem] = rec
    return rec


plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "legend.title_fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.4,
    # TYPE 42, NOT MATPLOTLIB'S DEFAULT TYPE 3. Measured 2026-09-08: every
    # committed .pdf carried `/Subtype /Type3` and NO `/FontFile` of any kind,
    # because matplotlib defaults `pdf.fonttype` and `ps.fonttype` to 3 and this
    # block never overrode them. Type 3 glyphs are drawn as PDF content streams
    # rather than backed by an embedded font program, so the figures RENDER
    # correctly everywhere -- which is exactly why nobody noticed -- while being
    # the format journal production systems most commonly reject or flag at
    # preflight. `ps.fonttype` matters for the same reason and is the sharper
    # exposure of the two: EPS is on AACR's accepted-source-format list for
    # revised manuscripts (SUBMISSION-CHECKLIST.md section 3) and PDF is not, so
    # the path this project would actually be asked for inherits the same
    # default. 42 means TrueType, embedded as a subset with `/FontFile2`.
    #
    # This does not affect the .png (raster) or the .svg (`svg.fonttype` is
    # `path`, so SVG text is already converted to outlines and carries no font
    # dependency at all -- which is why the poster embeds SVGs safely).
    "pdf.fonttype": 42, "ps.fonttype": 42,
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

    fig, ax = plt.subplots(figsize=canvas("figure0_schematic", (7.8, 4.3)))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    if POSTER:          # saved uncropped, so the drawing must fill the canvas
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    def box(x, y, w, h, text, fc, ec, size=8, weight="normal"):
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
        "-- both load on the same axis", "#FFF3CD", "#B8860B", 8)
    arrow(6.12, 3.85, 6.70, 3.85)

    # Caption sits in the clear band BETWEEN the arrowheads and the boxes; at
    # y=2.05 it ran through the H2 box and both dashed arrows.
    ax.text(5.0, 1.62, "Two corrections make the hypotheses distinguishable",
            ha="center", fontsize=8.8, fontweight="bold")

    box(0.55, 0.30, 3.9, 1.05,
        "1.  Residualize the SIGNATURE\non the global axis, within cancer type",
        "#E3F2E8", "#009E73", 8)
    box(5.10, 0.30, 4.3, 1.05,
        "2.  Disattenuate using the\nRESIDUALIZED score's reliability",
        "#E3F2E8", "#009E73", 8)
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=canvas("figure1_isi", (FIG_W, 5.4)),
                                   layout="constrained",
                                   gridspec_kw={"width_ratios": [1.25, 1]})

    ax1.axvline(0, color="k", lw=0.8, zorder=1)
    for d, c, off, lab in ((pan, C_PAN, 0.16, "Pan-TCGA (n=7,168)"),
                           (nsc, C_NSCLC, -0.16, "NSCLC (n=944)")):
        ax1.hlines(y + off, d["excess_lo"], d["excess_hi"], color=c, lw=1.6, alpha=.85)
        ax1.plot(d["excess_z"], y + off, "o", ms=3.6, color=c, label=lab, zorder=3)
    ax1.set_yticks(y)
    ax1.set_yticklabels([o.replace("_", " ") for o in order])
    ax1.set_xlabel("Immune-specific excess\n(Fisher $z$)")
    ax1.set_title("A  Per-signature ISI,\ndisattenuated", loc="left", fontweight="bold")
    # Below the axes: inside, at 8 pt, the key sat on the lowest signature's CIs.
    # On the poster's short, wide canvas both keys go to the right of panel B,
    # as one legend (below).
    if not POSTER:
        ax1.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)

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
        vp = ax2.violinplot(data, positions=y, orientation="horizontal", widths=0.85,
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
                 fontsize=8, color="#B00020")
    ax2.plot(pan["r_residual_refit"], y, "o", ms=4, color=C_PAN,
             label="Curated signature", zorder=3)
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlabel("Axis-residualized $r$\nwith image prediction")
    ax2.set_title("B  Pan-TCGA: curated vs\nrandom", loc="left", fontweight="bold")
    if POSTER:
        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax2.legend(h1 + h2, l1 + l2, frameon=False, loc="upper left",
                   bbox_to_anchor=(1.02, 1.0), ncol=1)
    else:
        ax2.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)
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

    # THE POSTER DRAWS PANEL A ALONE (TODO-PRINT-6, 2026-09-24). Its box is
    # 307 x 114 mm, and two panels with their keys cannot hold 18 pt text in it
    # (measured: the plots shrank to ~15 mm tall and the keys ran off the
    # canvas). Panel B's content is already printed as numbers in poster panel
    # 5's own text -- the random sets' median alpha of 0.974, falling to 0.801
    # without the axis, while curated signatures hold at 0.976 -- which is the
    # rule TODO-PRINT-5 applied when figure 4 left the poster.
    if POSTER:
        fig, ax1 = plt.subplots(figsize=canvas("figure2_reliability", (FIG_W, 4.0)),
                                layout="constrained")
        ax2 = None
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIG_W, 4.0), layout="constrained")

    for d, c, m, lab in ((pan, C_PAN, "o", "Pan-TCGA"), (nsc, C_NSCLC, "s", "NSCLC")):
        ax1.scatter(d["k"], d["reliability_gap"], s=26, color=c, marker=m,
                    label=f"{lab}, residualized", zorder=3, edgecolor="none")
        ax1.scatter(d["k"], d["alpha_observed_raw"] - d["alpha_null_mean_raw"], s=22,
                    facecolor="none", edgecolor=c, marker=m, lw=.8,
                    label=f"{lab}, raw score")
    rho, p = sps.spearmanr(pan["k"], pan["reliability_gap"])
    rho_n, p_n = sps.spearmanr(nsc["k"], nsc["reliability_gap"])
    ax1.set_xlabel("Genes in signature (after filtering)")
    ax1.set_ylabel(by_mode("Curated $-$ random Cronbach's $\\alpha$",
                           "Curated $-$ random\nCronbach's $\\alpha$"))
    ax1.set_title(by_mode("A  Reliability gap scales\nwith panel size",
                          "Reliability gap scales with panel size"),
                  loc="left", fontweight="bold")
    # The statistic goes INSIDE the axes: as a second title line it overprinted
    # panel B's title, and the k=36 point sat underneath the legend.
    ax1.text(0.97, 0.93,
             f"Spearman $\\rho$ = {rho:.2f}, p = {p:.3f}  (pan-TCGA)\n"
             f"$\\rho$ = {rho_n:.2f}, p = {p_n:.3f}  (NSCLC)",
             transform=ax1.transAxes, ha="right", va="top", fontsize=8)
    # The key goes BELOW the axes: at 8 pt, inside, it covered two points and
    # the y-axis label. The y range leaves the band above 0.5 empty for the
    # statistic, which sits above the highest point (0.42). On the poster the
    # key goes to the right, where the wide box has room.
    ax1.legend(frameon=False, **by_mode(
        {"ncol": 2, "loc": "upper center", "bbox_to_anchor": (0.5, -0.2),
         "columnspacing": 0.8, "handletextpad": 0.3},
        {"ncol": 1, "loc": "upper left", "bbox_to_anchor": (1.02, 1.0)}))
    ax1.axhline(0, color="k", lw=0.6)
    ax1.set_ylim(-0.04, 0.62)
    ax1.set_xlim(0, 225)

    if ax2 is not None:
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
        ax2.set_xticklabels(["Raw score", "Axis-residualized"])
        ax2.set_ylabel("Median Cronbach's $\\alpha$")
        ax2.set_ylim(0, 1.05)
        ax2.set_title("B  Random sets are consistent\nonly because they load\non the axis",
                      loc="left", fontweight="bold")
        ax2.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1)
        # The value goes INSIDE its bar, in white: above the bar it ran into the
        # curated signature's mark, which sits just over the random-set median.
        for i, d in enumerate((pan, nsc)):
            for j, v in enumerate([d["alpha_null_mean_raw"].median(),
                                   d["alpha_null_mean"].median()]):
                ax2.text(j + (i - .5) * w, v - .03, f"{v:.2f}", ha="center", va="top",
                         fontsize=8, color="white", fontweight="bold")

    save(fig, outdir, "figure2_reliability")
    print(f"  figure2: Spearman rho={rho:.3f} p={p:.4f}; "
          f"pan gap raw {(pan.alpha_observed_raw-pan.alpha_null_mean_raw).mean():+.3f} "
          f"-> residualised {pan.reliability_gap.mean():+.3f}")


# ------------------------------------------------------------------ fig 3 ---

def figure3(outdir: Path) -> None:
    """Outcome excess over null, pan-TCGA, by MSigDB process category."""
    o = load(PAN, "outcome_arm")
    # `axis_residualised` is a DATA VALUE in the frozen, provenance-hashed
    # outcome_arm.csv of six result directories -- not prose. It keeps the
    # British spelling deliberately. US spelling is enforced on text this
    # script DRAWS; a find-and-replace here would select nothing and empty
    # the figure silently.
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

    fig, ax = plt.subplots(figsize=canvas("figure3_outcome", (FIG_W, 5.0)), layout="constrained")
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
    ax.set_title(by_mode("Prognostic advantage over random gene sets is\n"
                         "proliferative and stromal, not immune\n(PFI, pan-TCGA n=7,168)",
                         "Prognostic advantage over random gene sets is proliferative\n"
                         "and stromal, not immune (PFI, pan-TCGA n=7,168)"),
                 loc="left", fontweight="bold")

    handles = [plt.Line2D([], [], marker="o", ls="", color=CAT_COLOUR[k], label=k)
               for k in ["immune", "proliferation", "development", "pathway", "signaling"]]
    handles.append(plt.Line2D([], [], marker="o", ls="", color="w",
                              markeredgecolor="k", label="beats null"))
    # On the poster's wide canvas both keys sit to the right of the axes.
    leg1 = ax.legend(handles=handles, frameon=False, title="MSigDB process category",
                     **by_mode({"loc": "lower right"},
                               {"loc": "lower left", "bbox_to_anchor": (1.02, 0.0)}))
    ax.add_artist(leg1)
    # Below the axes: inside, at 8 pt, this key sat on the top row's NSCLC mark.
    ax.legend(frameon=False, **by_mode(
        {"loc": "upper center", "bbox_to_anchor": (0.5, -0.1), "ncol": 2},
        {"loc": "upper left", "bbox_to_anchor": (1.02, 1.0), "ncol": 1}))

    save(fig, outdir, "figure3_outcome")
    imm = o[o["cat"] == "immune"]
    print(f"  figure3: {int(o.beats_null.sum())}/16 beat null; "
          f"immune {int(imm.beats_null.sum())}/{len(imm)}; MDE {mde:.4f}")


# ------------------------------------------------------------------ fig 4 ---

def figure4(outdir: Path) -> None:
    """The control that decides how much of this to believe: site recoverability."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIG_W, 3.4), layout="constrained")

    for path, c, lab in ((PAN, C_PAN, "Pan-TCGA"), (NSCLC, C_NSCLC, "NSCLC")):
        sc = pd.read_csv(path / "site_control.csv")
        ax1.hist(sc["auroc"], bins=np.linspace(0.5, 1.0, 26), color=c, alpha=.7,
                 label=f"{lab} ({len(sc)} sites, median {sc['auroc'].median():.3f})")
    ax1.axvline(0.5, color="k", lw=0.8, ls="--")
    ax1.set_xlabel("One-vs-rest AUROC for tissue source site")
    ax1.set_ylabel("Sites")
    ax1.set_title("A  Site is almost perfectly recoverable\nfrom the embeddings",
                  loc="left", fontweight="bold")
    ax1.legend(frameon=False, loc="upper left")

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
    # Above the bars, in the headroom: at 8 pt in the corner it sat on a bar.
    ax2.set_ylim(0, 0.8)
    ax2.legend(frameon=False, loc="upper center", ncol=2)

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
                    help="also write SVG. The committed journal SVGs are compared "
                         "against a fresh render by the test suite; the poster no "
                         "longer embeds them (see --poster)")
    ap.add_argument("--poster", action="store_true",
                    help="render the four A0 poster figures (TODO-PRINT-6, "
                         "decided 2026-09-24): each drawn with the exact aspect "
                         "of its poster box at 1/POSTER_SCALE its size and saved "
                         "uncropped as SVG, so the poster enlarges it by "
                         "POSTER_SCALE and no text prints under POSTER_MIN_PT "
                         "(18 pt). Writes to ../poster/figures/, which the poster "
                         "embeds and which is committed; journal figures are "
                         "untouched.")
    ap.add_argument("--write-gene-sizes", action="store_true",
                    help=f"recompute the 16 panel sizes from data/ and write "
                         f"{_rel(GENE_SIZES)}, the deposit that lets figure 2 "
                         f"render from the public snapshot. Renders nothing.")
    args = ap.parse_args()

    if args.poster:
        global POSTER
        POSTER = True
        if args.outdir == ap.get_default("outdir"):
            args.outdir = ROOT.parent / "poster" / "figures"
        if args.tiff or args.svg:
            raise SystemExit("FATAL: --poster writes SVG only; drop --tiff/--svg")

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

    if POSTER:
        print(f"rendering the A0 poster figures into {_rel(args.outdir)}: svg, "
              f"uncropped, at 1/{POSTER_SCALE:g} of each box")
        figure0(args.outdir)
        fig1b_null = figure1(args.outdir)
        figure2(args.outdir)
        figure3(args.outdir)
        print("\nprinted on A0 (1 CSS px = 25.4/96 mm): smallest text")
        bad = []
        for stem, (w_px, h_px) in POSTER_BOXES_PX.items():
            r = PRINT_SIZES[stem]
            ok = (r["printed_smallest_pt"] >= POSTER_MIN_PT - 1e-9
                  and not r["outside_canvas"])
            print(f"  {stem:22s} {r['printed_smallest_pt']:5.2f} pt  (box "
                  f"{w_px * 25.4 / 96:.0f} x {h_px * 25.4 / 96:.0f} mm, "
                  f"{r['n_text']} text items)"
                  f"{'' if ok else '  <-- BELOW THE FLOOR OR OFF THE CANVAS'}")
            if not ok:
                bad.append(stem)
        if bad:
            raise SystemExit(
                f"FATAL: {bad} would print text under {POSTER_MIN_PT:g} pt on "
                f"A0, or draw text off the canvas "
                f"({ {s: PRINT_SIZES[s]['outside_canvas'] for s in bad} }).")
        return 0

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
    # THE PRINT-SIZE CHECK, on every journal render (H87). The committed
    # figures are the manuscript's, so a figure whose smallest text prints
    # under the floor at double-column width, or that draws text past its own
    # canvas, is refused rather than written silently. Poster renders are
    # judged against their A0 boxes instead, above.
    print(f"\nprinted at {JOURNAL_WIDTH_IN:g} in (double column): smallest text")
    too_small = []
    for stem in JOURNAL_FIGURES:
        r = PRINT_SIZES[stem]
        ok = (r["printed_smallest_pt"] >= JOURNAL_MIN_PT - 1e-9
              and not r["outside_canvas"])
        print(f"  {stem:22s} {r['printed_smallest_pt']:5.2f} pt  "
              f"(saved {r['width_in']:.3f} in wide, {r['n_text']} text items)"
              f"{'' if ok else '  <-- BELOW THE FLOOR OR OFF THE CANVAS'}")
        if not ok:
            too_small.append(stem)
    if too_small:
        raise SystemExit(
            f"FATAL: {too_small} would print text under {JOURNAL_MIN_PT:g} pt at "
            f"{JOURNAL_WIDTH_IN:g} in, or draw text off the canvas "
            f"({ {s: PRINT_SIZES[s]['outside_canvas'] for s in too_small} }). "
            "Raise the type or narrow the canvas; do not lower the floor.")

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
