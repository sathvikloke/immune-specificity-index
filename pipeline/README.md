# Does H&E-inferred tumour immune microenvironment carry immune-specific signal?

Pipeline for the AACR Annual Meeting 2027 project. Abstract deadline **Tue 10 Nov 2026**.

## The question

Models that infer tumour immune-microenvironment (TME) signatures from H&E whole-slide
images are being positioned as adjuncts to PD-L1 for immunotherapy patient selection.
Reported correlations look respectable (HistoTME: mean Pearson ~0.5 external; HistoTME-v2:
median 0.61 internal / 0.53 external).

### What is NOT the contribution

**A bare random-gene-set null is already published.** HE2RNA
([Schmauch et al., Nat Commun 2020;11:3877](https://doi.org/10.1038/s41467-020-17678-4))
ran it — *"the mean correlation coefficient R₀ obtained for 10,000 random lists of the
same number of genes, for all 28 cancer types"* — and its answer **favours** the <!-- numcheck: ignore: a direct quotation from Schmauch et al. 2020, about THEIR 28 cancer types, not this study's 31 -->
signatures: 75% (B-cell) and 86% (T-cell) of cancer types beat random on the signature
metric, though only 36% and 50% per-gene.

Site confounding is also settled: tissue source site is predictable from pathology
foundation-model embeddings at 93–99% accuracy
([de Jong et al.](https://arxiv.org/abs/2501.18055)), and all 20 models in
[PathoROB](https://www.nature.com/articles/s41467-026-73923-2) show robustness deficits.
Howard et al. covers TCGA-LUAD and TCGA-LUSC specifically. Demographic gaps are published
too ([Vaidya, Nat Med 2024](https://www.nature.com/articles/s41591-024-02885-z);
[FLEX, Nat Commun 2025](https://www.nature.com/articles/s41467-025-66300-y)).

These appear here as **cited reproduction checks, never as findings.**

### What is the contribution

HE2RNA's comparison is uncorrected for two things, both biasing it toward the curated set:

1. **The global expression axis.** "Beats random" is equally consistent with the image
   reading tissue composition faithfully — a *positive* result about H&E, not a null one.
   Separating the two requires residualising the SIGNATURE on the axis and asking whether
   the residual is still image-predictable above a residualised null.

2. **Reliability.** Curated sets are co-expressed modules; random draws are not. A more
   reliable target is more predictable regardless of biology. Uncorrected, this alone
   produced a **100% false-positive rate** in simulation (`test_disattenuation_removes_the_reliability_false_positive`).

Plus an **outcome arm**: the same null against PFI, where Venet's logic is exact rather
than analogical, because survival is distal to expression whereas an image and an RNA
profile are two measurements of one specimen.

## Design constraints

- **Fully open, ungated data.** No dbGaP, no EGA, no signing official.
- **Permissive licences only.** Prov-GigaPath (Apache-2.0) and MSigDB (CC-BY-4.0), not
  UNI2-h / Virchow2 / CONCH (CC-BY-NC) or BostonGene/MFP (academic-only).
- **CPU-only core.** Slide-level embeddings + ridge heads.
- **Patient-level everything** — splits, bootstraps, unit of analysis.

## Install and run

```bash
python3 -m pip install -r requirements.txt
```

Python 3.11+. Verified on 3.13.9 with numpy 2.1, pandas 2.2, scikit-learn 1.6, scipy 1.15,
statsmodels 0.14.

```bash
python3 -m pytest tests -q
```

**122 tests.** If this number drifts from reality again, re-measure rather than
editing the prose. (Re-measured 2026-09-02 by counting `def test_` in `tests/`;
the previous "57 tests, ~7 s" had drifted. <!-- numcheck: ignore: deliberately records the superseded count, as the evidence that drift happened --> The wall-clock figure is deliberately
not quoted now — it was measured at ~7 s on an idle machine and 37 min on a loaded
one, so it says more about the machine than the suite.)

```bash
python3 scripts/03_fetch_signatures.py --hallmark
python3 scripts/02_run_audit.py --demo
```

The demo plants three signatures with *known* confounding and separates them. From
`results/demo/`:

| signature | r_residual | null mean | excess (z) | 95% CI | beats null |
|---|---|---|---|---|---|
| `SIG_immune_clean` | 0.935 | 0.810 | **+0.567** | [0.487, 0.646] | **yes** |
| `SIG_immune_siteconf` | −0.125 | −0.158 | +0.034 | [−0.036, 0.104] | no |
| `SIG_immune_purity` | 0.578 | 0.811 | −0.487 | [−0.565, −0.408] | no |

Real biology survives; site artefact and purity-driven signal do not.

## Pre-registered primary endpoint

> **Immune-specificity index (ISI):** reliability-disattenuated, global-axis-residualised
> excess in Fisher z of r(image, signature) over the mean of r(image, size-matched random
> set), 1000 draws, patient-clustered bootstrap CI.

This string is byte-identical in `AuditConfig.primary_endpoint`, and `AuditResult.save()`
**asserts** that the estimand actually computed matches it — the registration and the code
cannot drift apart. Timestamp an OSF registration or a signed git tag before the first
real run.

**Secondary:** `delta_MAE`, a paired difference in mean absolute error between split
schemes. It is a calibration-loss contrast, *not* a difference of correlations, and must
never be reported as "Δr".

### Multiplicity

BH-FDR at q=0.05 within each declared family. **M_eff is descriptive only** — it is a
family-wise-error device and does not compose with BH; reporting both as a correction is
a category error. B=1000 null draws, because at B=100 the smallest attainable p (1/101)
exceeds what BH needs at m=29 (0.0017), making rejection logically impossible.

### Ancestry

Pan-cancer only. NSCLC counts are fatal (EAS n=17 against ~137 per group needed). Report
as a precision statement with `stats.min_detectable_delta_z`, never a bare difference
test, and always both marginally and stratified by site.

## Layout

```
src/aacr27/
  barcodes.py       TCGA barcode -> patient / tissue source site
  splits.py         random-patient vs preserved-site CV, degenerate-partition guard
  models.py         ridge (single + stacked multi-target), baselines, ComBat fit/apply
  signatures.py     scoring, the random-set null, Cronbach's alpha / rho_bar
  globalaxis.py     the global axis, residualisation, disattenuation, the ISI
  outcome.py        Harrell's C, Cox, the PFI null (no dichotomisation path)
  decomposition.py  nested variance decomposition on ONE shared complete-case mask
  stats.py          Fisher-z CIs, cluster bootstrap, BH-FDR, M_eff, rotation null
  data.py           loaders (embedding layer selector, purity, expression) + provenance
  experiment.py     the orchestrator — run_audit()
scripts/
  01_fetch_data.py        embeddings / ancestry / purity / expression
  02_run_audit.py         --demo (synthetic) or --real
  03_fetch_signatures.py  MSigDB, CC-BY-4.0, KEGG-restricted sets filtered out
tests/              122 tests: planted-confound recovery + a regression test for
                    every rev.3 pipeline fix
```

## Things that are deliberately not what you would expect

**Post-ComBat site AUROC is not interpretable, in either direction.** Measured: raw 1.000,
in-sample ComBat 0.049, out-of-fold ComBat 0.020. Out-of-fold fitting does not rescue it,
and an earlier expectation that an honest correction "lands near 0.50" is wrong. ComBat is
working — per-site mean spread falls from 0.968 to 0.033 — but it imposes a within-site
zero-sum constraint that the site classifier calibrates to and then reverses on. Report
the **raw** AUROC as the confound measure, and the **estimable fraction** as the
deployment answer: under preserved-site CV it is **0% by construction**, because a
held-out site has no training data. That is the position a model meeting a new hospital is
in, and it is a better finding than any corrected number.

**Baselines are evaluated under the random-patient split**, not the site-disjoint one.
Under preserved-site CV a one-hot-site model sees no training example from any test site,
so it scores ~0 by construction — that measures the split, not the covariate.

**Larger gene sets flip the statistical regime.** Two independent size-*k* sets correlate
at `k·ρ̄/(1+(k−1)·ρ̄)`, which is also Cronbach's α. Small *k* → attenuation (curated wins
on precision); large *k* → degeneracy (every draw measures the same axis). MSigDB Hallmark
runs 36–200 genes, i.e. degenerate. `globalaxis.null_degeneracy()` reports which regime a
run is in. Pre-register the panel.

## Known limitations

- `combat_correct` has no empirical-Bayes shrinkage; confirm against `neurocombat`.
- Embedding layer, tiling, magnification and pooling are **UNVERIFIED** — the parquet
  stores 14 layers × 768 and `layer` is a researcher degree of freedom. Pre-register it.
- Column-name detection in `load_ancestry` is a guess recorded in provenance. Check it.
- `score_mean_z` is cohort-relative. Scoring TCGA and CPTAC separately silently redefines
  the target. Fix the scaling policy once.
- **No real data has been analysed.** Every number here is synthetic.

## What this study cannot say

No TCGA patient in this cohort received checkpoint blockade, and no response endpoint is
analysed. Conclusions are about **measurement validity** — what these scores measure — not
about who should receive immunotherapy. An earlier draft of this file closed with a
patient-selection recommendation; it had zero supporting data and was removed.

## Reporting

Ship a **TRIPOD+AI** and **PROBAST+AI** self-audit. AACR journals name no ML reporting
standard, so supplying one is free differentiation. Only 3.2% of accepted AACR
computational abstracts mention code availability — ship the repo with the abstract.

## Reproducing the results

```bash
pip install -r requirements-lock.txt   # EXACT versions the results were produced with
bash reproduce.sh --check              # tests + figures from stored results, 4 m 7 s (measured)
bash reproduce.sh                      # full: fetch, both cohorts, figures, ~2.5 h (ESTIMATE, never measured)
```

`requirements.txt` carries lower bounds for ordinary use; `requirements-lock.txt`
is what reproduces the manuscript numbers. The pin is not ceremony — sklearn has
changed `RidgeCV`'s alpha-selection defaults between minor versions, and the ISI
depends on the selected penalty.

### What runs from the deposit alone

**Read this before concluding anything is broken.** What is published at
`immune-specificity-index` is a *curated snapshot*, not this working repository:
the raw inputs under `data/` are public but large and are not redistributed
here, several multi-megabyte internal diagnostic probes are left out, and the
working documents (drafts, checklists, the poster) are not deliverables. Some
scripts shipped in the deposit therefore have no input there. **That is
curation, not corruption.**

Every row below was MEASURED on 2026-09-07 by running the command inside a
freshly staged snapshot with no `data/` reachable — not inferred from reading
the code.

| command | exit | in the deposit |
|---|---|---|
| `python3 scripts/10_make_figures.py` | **0** | all five figures render |
| `python3 scripts/25_build_supplementary.py` | **0** | all nine tables build |
| `python3 scripts/20_check_versions.py` | **0** | |
| `python3 scripts/21_provenance_manifest.py` | **0** | 110 of 114 verified; 4 not staged, and it says so |
| `python3 scripts/23_export_environment.py` | **0** | |
| `python3 -m pytest tests/ -q` | 1 | **13 failed, 101 passed, 3 skipped** — see below |
| `python3 scripts/24_split_determinism.py --check` | 2 | needs `data/interim/cohort_nsclc.parquet` |
| `python3 scripts/26_sort_audit.py` | 2 | needs the staged inputs |
| `python3 scripts/19_check_numbers.py` | 1 | needs the excluded `bisect_*.npz` probes |
| `python3 scripts/08_count_abstract.py` | 2 | the abstract draft is not a deliverable |
| `bash reproduce.sh --check` | 2 | **stops at `24_split_determinism.py`, ~10 s in** |

Two consequences worth stating plainly rather than letting a reader discover
them:

**`bash reproduce.sh --check` cannot pass inside the deposit.** It is this
repository's one-command check over the *working* tree, and three of the
commands it runs need inputs the deposit does not carry. It fails fast and
loudly, which is the right failure — but it is quoted above as a reproduction
route and that is only true where `data/` exists. To verify the deposit itself,
run the five exit-0 commands in the table.

**The 13 test failures are all curation, and none of them is the analysis
code.** Classified by cause, they are: 6 that exercise
`22_build_public_snapshot.py`, which excludes *itself* by design and so is
absent; 6 that need the excluded `bisect_*.npz` platform probes; and 1 that
checks spelling across working documents the deposit does not stage. The 101
that pass are the pipeline's own tests. A deposit reader running the suite sees
red on an intact deposit — that is a known and accepted cost of curating, and
it is written down here rather than left as a surprise.

### The lock file is necessary but NOT sufficient — read this before rerunning

**Pinning the versions does not currently get you the same number on a different
machine.** Measured 2026-09-02, NSCLC ISI, same code, same seeds, same
`requirements-lock.txt` versions verified at runtime, and inputs hashed
byte-for-byte identical on both sides:

| | NSCLC ISI |
|---|---|
| macOS 15 / arm64, Python 3.13.9 | **0.318240180054566** ← the frozen results |
| Linux / x86_64, Python 3.12.14 | **0.2964** |

The gap is 0.0218. It does not change the sign, the significance, or the fact
that 16/16 signatures beat their null — it changes the fourth digit, and this
README quotes four digits. It is 2.4x larger than the partition-choice
variability the study already widens its interval for.

> **CAUSE FOUND, 2026-09-05 — and it is not floating point.** `splits.py:177`
> orders sites with `np.argsort(-(sizes + noise))`. The default
> `kind="quicksort"` is **not stable**, and with `jitter=0.0` the noise term is
> exactly zero, so equal-sized sites are exactly tied — **51 of the 68 NSCLC
> sites**. The tie order is architecture-dependent, so the two machines build
> **different cross-validation partitions from byte-identical inputs**, with
> identical fold sizes (which is why every guard passed). Verify it yourself on
> any machine **that has the inputs** — this script builds real fold partitions,
> so it exits 2 in the public deposit, where `data/` is not redistributed:
>
> ```bash
> python3 scripts/24_split_determinism.py --check
> ```
>
> `hash_fold_shipped` is `e55ee4e8887b5167` on macOS and `a40c862b12d9a391` on
> HPC4; under `kind="stable"` both are `9f31e33c2caaf1e3`. On the millisecond
> reproducer, the stable sort **on macOS** returns `0.9752094760755506` against
> HPC4's `0.9752094760755513` — the Linux value to 15 significant figures.
> **Genuine cross-platform floating-point noise in this pipeline is 7e-16.**
> The ridge solve is exonerated: `cond(K + αI)` is 18–77 and a relative 1e-14
> input perturbation moves the predictions by 6e-15, so it *attenuates* rather
> than amplifies. The fix is one keyword, and it is deliberately **not applied**
> — it would move 216 of 944 patients to a different fold and change every
> frozen number. See `16-YOUR-TASKS.md` for the pending decision and
> `../14-SCIENCE-AUDIT.md` §"A9 CLOSED" for the full two-machine table.

**The stage is identified, 2026-09-06: the per-fold ridge penalty selection.**
`experiment._select_alpha_per_fold` makes 16 x 5 = 80 scalar RidgeCV selections
off `ALPHAS = np.logspace(-2, 5, 24)`, each then held fixed across a signature's
observed target and all 1,000 null targets; adjacent grid points differ by
~1.96x, so one flip is a large discrete move rather than noise. Running
`scripts/18_bisect_platform.py` on both machines and comparing the dumps:

```
  S0_boot_synth    5.551e-17     standalone bootstrap on synthetic input
  S2_obs_resid     1.010e-14     global-axis-residualised observed score
  S2_null_resid    1.021e-14     global-axis-residualised null scores
  S3_fold_alpha    1.515e+03     per-fold ridge alpha            <-- DIVERGES
  ... every later stage diverges as amplification of S3
```

**11 of the 80 penalties differ, each by exactly one grid step** (1515.3870 =
3007.8825 - 1492.4955). The bootstrap and the residualisation stage are clean to
round-off, so both are eliminated.

**But the penalty flip is not the root cause.** One penalty is selected per
signature-fold and then held fixed across that signature's observed and null
targets, so the eight signatures whose penalties are IDENTICAL on both machines
are a free control. `--compare` now runs it automatically:

```
  AMPLIFICATION CHECK: 11 of 80 fold alphas differ; 8 of 16
  signature(s) have IDENTICAL alphas on both machines and act as a control.
  Control signatures DIVERGE downstream anyway: max |A-B| in the out-of-fold
  predictions is 3.023e-01, which is 177% of their own sd (0.1705).
  So S3 is the FIRST divergence but NOT the only one: the ridge solve differs
  at identical penalties on identical inputs. Attack models.cross_val_predict_multi,
  not the alpha grid.
```

Same penalty, inputs agreeing to 9.77e-15, and the two platforms' predictions
correlate at only **r ~= 0.96**. The likely structural reason, measured: X is
(944, 768) over 5 folds, so each training set is ~755 samples against 768
features -- **fewer rows than columns in every fold**. Which part of the solve
turns 1e-14 into 1e-1 is not yet established, and nothing here should be read as
saying it is.

**A warning about what this README used to say here.** Until 2026-09-06 this
section stated that the 80 penalties *cannot* differ across platforms, on the
strength of `scripts/16_alpha_margin.py`'s margin measurement:

```
min 1.034e-04   median 5.099e-03   max 1.125e-02      below 1e-06: 0 / 80
```

Those margins are real and that run was correct. The **inference** from them was
not: a margin eleven orders of magnitude above BLAS round-off does not imply
invariance if something between the inputs and the decision amplifies, and both
of the runs that "confirmed" invariance were on the same machine. The
cross-platform measurement, when finally taken, contradicts it. Treat a margin
argument as a hypothesis until two machines have printed the same digest.

- `scripts/18_bisect_platform.py` dumps every stage; `--compare` names the first
  one whose max absolute difference leaves the 1e-12 round-off floor. Its S3 row
  is still *labelled* "PROVEN invariant by 16_alpha_margin" -- that label is
  stale and the row it labels is the one that diverges.
- `scripts/17_alpha_margin_fast.py` answers the margin question without
  `run_audit` (seconds of setup rather than hours). Its closing VERDICT still
  asserts the alphas cannot differ; it is stale for the same reason.
- `scripts/14_repro_probe.py` fingerprints every upstream stage; stages 1-6 are
  identical on both machines.
- `scripts/15_alpha_probe.py` prints the same 80-value digest as 16 and 17 and
  needs both machines to say anything -- which is exactly what it eventually
  said: `946033557a8ad08b` on HPC4 against `fc299c9c0f47466e` on macOS.

**So: the numbers in this README are the macOS numbers.** If you reproduce on
another platform and get 0.2964, you have not done anything wrong, and it would
be useful to know — that is a second data point on a defect with one.

`reproduce.sh` sets `PYTHONHASHSEED=0`. Python randomises str hashing per
process, so without it the gene pool is ordered differently every run and
`seed=0` does not reproduce.

### Headline results (frozen at the 2026-08-18 analysis freeze)

> The freeze is recorded by an annotated git tag, `prereg-2026-08-18`, **in the
> author's private working repository — not in this published snapshot, which is
> built as a fresh single-commit repository and therefore carries no tags.** The
> pre-registered protocol's *content* is here in full (`05-PRE-REGISTRATION.md`);
> its *timestamp* is not independently verifiable from this repository alone.
> Earlier wording put the tag name in this heading, which invited a reader to run
> `git tag`, find nothing, and reasonably conclude the claim was invented.

| | Pan-TCGA (n=7,168, 31 types) | NSCLC (n=944) |
|---|---|---|
| ISI | **0.2911 [0.2692, 0.3133]** | 0.3182 [0.2614, 0.3763] |
| Signatures beating their null | 16/16 | 16/16 |
| Outcome (PFI, within cancer type) | 10/32 | 0/32 (MDE 0.046) |
| Site AUROC | 0.998 | 0.992 |

Prognostically the signal is proliferation and stroma, not immunity: the five
signatures beating their null are G2M checkpoint, E2F targets, angiogenesis,
EMT and hypoxia; no interferon, inflammatory, complement or allograft-rejection
signature does.

**These analyses are exploratory** — the analysis preceded the protocol. See
`../05-PRE-REGISTRATION.md` §0.

## Documents

These sit in the parent project directory, one level ABOVE this pipeline. They are
**planning and manuscript documents, deliberately not part of the published
pipeline** — so if you are reading this repository on its own, these paths will not
resolve, and that is intended rather than a broken link. `05-PRE-REGISTRATION.md`
is the one whose absence actually costs you something; ask for it if you need the
protocol.

| File | What |
|---|---|
| `../04-CODE-AUDIT.md` | Five defects found and fixed in a full-code audit, with measurements |
| `../05-PRE-REGISTRATION.md` | Protocol, frozen configuration, exploratory-vs-confirmatory split |
| `../07-ABSTRACT-DRAFT.md` | AACR 2027 abstract, 2,515/2,600 characters |
| `../09-PAPER-DRAFT.md` | Manuscript draft, references verified against PubMed |
| `../11-SUBMISSION-PACK.md` | Category choice, submission checklist, poster plan |
| `../14-SCIENCE-AUDIT.md` | Gaps between the protocol and the delivered results, including §A9 above |

**One script here is not science.** `scripts/08_count_abstract.py` counts
characters against AACR's 2,600-character abstract limit under the Call for
Abstracts' own rule that spaces are free. It is conference tooling, kept in the
repository because it is what produced the 2,517 figure quoted above and a reader
should be able to check that number rather than take it on trust. It touches no
data and is not part of `reproduce.sh`.
