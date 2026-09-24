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
- **Permissive licences only.** Not UNI2-h / Virchow2 / CONCH (CC-BY-NC) or
  BostonGene/MFP (academic-only). Two objects are involved and they carry
  DIFFERENT licences, which is worth stating plainly because the two are easy to
  conflate: the Prov-GigaPath **encoder model** is Apache-2.0, while what this
  pipeline actually consumes is the derived **TCGA slide-embeddings dataset**
  `seandavis/tcga_provgigapath_embeddings`, which is CC-BY-4.0. MSigDB is
  CC-BY-4.0. Both licences were verified at source on 2026-09-08 against the
  HuggingFace cards. Note also that the encoder model is *gated* on HuggingFace
  while the embeddings dataset is not — which costs this pipeline nothing,
  because it never downloads the model, only the embeddings.
- **CPU-only core.** Slide-level embeddings + ridge heads.
- **Patient-level everything** — splits, bootstraps, unit of analysis.

## Input provenance — the six external sources

Every source the manuscript's Data Availability Statement names, with the exact
accession that was downloaded and the date it was fetched. Recorded 2026-09-08;
before that, the manuscript named the sources and nothing recorded *which file*.
Sizes are as fetched, on disk, and are why the raw inputs are not redistributed.

| source | exact file and pin (what `scripts/01_fetch_data.py` fetches and hash-verifies) | licence | size | fetched |
|---|---|---|---|---|
| Slide embeddings | HuggingFace dataset `seandavis/tcga_provgigapath_embeddings`, revision `073115403c2fc5134ee8d1332c603edba591dddb`, file `provgigapath_embeddings_with_metadata.parquet`, saved as `data/raw/provgigapath/embeddings.parquet` (ungated) | **CC-BY-4.0**, verified on the dataset card 2026-09-08 | 467 MB | 2026-08-17 |
| Expression | Xena TOIL `https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz`, log2(TPM + 0.001) as shipped (floor −9.9658) | open, no login | 741 MB | 2026-08-17 |
| Gene symbols | HGNC's complete set as served on 2026-08-17, reduced to `ensembl_gene_id,symbol`. HGNC does not archive that version (its August 2026 monthly archives map 41,037 of TOIL's genes, not 41,046), so the reduced map is committed at `pipeline/resources/ensembl_to_hugo.csv` and copied from there | **CC0** (genenames.org, checked 2026-09-17) | 1.0 MB map | 2026-08-17 |
| Tumour purity | GDC PanCanAtlas ABSOLUTE, `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` (`api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5`) → `TCGA_ABSOLUTE_purity.tsv`; the `.csv` the pipeline reads is the same file with tabs replaced by commas | open (GDC) | 902 KB | 2026-08-17 |
| Outcome | TCGA-CDR (Liu et al. 2018). **The pipeline reads the CDR columns that ship pre-merged inside the embeddings parquet**; no separate CDR file is read. PFI is used. | open | — | 2026-08-17 |
| Ancestry | Carrot-Zhang et al. 2020, `UCSF_Ancestry_Calls.csv` from the GDC CCG-AIM-2020 page (`api.gdc.cancer.gov/data/fdfa536a-c3c8-405d-99d9-bc9375b5084c`) | open | 1.9 MB | 2026-08-17 |
| Signatures | MSigDB `h.all.v2024.1.Hs.symbols.gmt`, fetched by `scripts/03_fetch_signatures.py` — the version is **pinned in the fetch URL**, not "latest" — which also writes the 16-set `.tme.gmt` subset | **CC-BY-4.0**, attribution in `data/raw/signatures/ATTRIBUTION.txt` | 49 KB | 2026-08-12 |

Every one of these was re-downloaded from the source named here on 2026-09-16/17
and reproduced the sha256 recorded in `01_fetch_data.py` (the gene map is the
exception explained in its row). `python scripts/01_fetch_data.py --verify`
checks a local copy against the same hashes.

**On the Prov-GigaPath licence, which looks contradictory and is not.** Two
different objects: the *encoder model* `prov-gigapath/prov-gigapath` is
**Apache-2.0** and is *gated* on HuggingFace; the *derived TCGA slide-embeddings
dataset* used here is **CC-BY-4.0** and ungated. This pipeline never downloads the
model — only the embeddings — so the gate costs reproduction nothing. Both were
read off their HuggingFace cards on 2026-09-08.

**From `data/raw/` to `data/interim/`.** The cohort scripts read six derived
files from `data/interim/`. `scripts/00_build_interim.py` builds them from the
raw sources above; it was reconstructed on 2026-09-16 from the files the
frozen results were computed from, and on macOS it reproduces five of the six
byte for byte (`--verify` compares against their recorded hashes; the sixth,
the 10,535 × 41,046 expression matrix, needs about 5 GiB and is built on
HPC4). One property of those files is a known defect rather than a design
choice: the expression profile is not always the tumour (A16 in
`14-SCIENCE-AUDIT.md`: 51 of 944 NSCLC and 618 of 7,168 pan-TCGA patients).
`--tumour-only` builds corrected inputs into a separate directory for the
sensitivity analysis.

## Install and run

```bash
python3 -m pip install -r requirements.txt
```

Python 3.11+. Verified on 3.13.9 (macOS, where the frozen results were produced)
and 3.12.14 (Linux), with numpy 2.1, pandas 2.2, scikit-learn 1.6, scipy 1.15,
statsmodels 0.14.

```bash
python3 -m pytest tests -q
```

**203 tests.** If this number drifts from reality again, re-measure rather than
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
cannot drift apart. The protocol was tagged `prereg-2026-08-18` in the private
working repository before the frozen results were computed; two changes made the
day after filing are recorded in its Appendix A.

**Secondary:** `delta_MAE`, a paired difference in mean absolute error between split
schemes. It is a calibration-loss contrast, *not* a difference of correlations, and must
never be reported as "Δr".

### Multiplicity

BH-FDR at q=0.05 within each declared family. **M_eff is descriptive only** — it is a
family-wise-error device and does not compose with BH; reporting both as a correction is
a category error. B=1000 null draws, because at B=100 the smallest attainable p (1/101)
exceeds what BH needs at m=16 (0.05/16 = 0.0031), making rejection logically
impossible. (This read m=29 until 2026-09-17, the signature count of an earlier
design; the registered family is the 16 Hallmark signatures.)

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
  ancestry.py       the supplementary ancestry arm and its site-identifiability check
  experiment.py     the orchestrator — run_audit()
scripts/
  00_build_interim.py     data/raw -> the six data/interim files (A13)
  01_fetch_data.py        fetches every raw input and verifies each against the
                          sha256 of the copy the results used (--all, --verify)
  02_run_audit.py         --demo (synthetic) or --real
  03_fetch_signatures.py  MSigDB v2024.1.Hs (pinned in the URL, not "latest"),
                          CC-BY-4.0. A KEGG_/KEGG_MEDICUS_/BIOCARTA_ filter
                          exists but is a MEASURED NO-OP on the default path:
                          the pipeline fetches Hallmark only, and 0 of those
                          50 sets carry a restricted prefix. It bites only for
                          the optional --c2cp collection.
  ...                     04-34: each script's docstring says what it does and
                          which audit or ledger item it closes
  27_partition_sweep.py   pools sharded partition sweeps, computes the spread once
tests/              203 tests: planted-confound recovery + a regression test for
                    every rev.3 pipeline fix
```

## Things that are deliberately not what you would expect

**Post-ComBat site AUROC is not interpretable, in either direction.** On synthetic data
with a planted multi-dimensional site fingerprint (the measurement recorded in the docstring
of `models.site_prediction_control_oof_combat`): raw 1.000, in-sample ComBat 0.049,
out-of-fold ComBat 0.020. On the frozen NSCLC run the out-of-fold median is 0.005
(`results/nsclc_v3/site_control_combat.csv`). Out-of-fold fitting does not rescue it,
and an earlier expectation that an honest correction "lands near 0.50" is wrong. ComBat is
working — on that synthetic data the per-site mean spread falls from 0.968 to 0.033 — but
it imposes a within-site zero-sum constraint that the site classifier calibrates to and
then reverses on. *(Scoped 2026-09-22: until then this passage said only "Measured:",
so the synthetic figures read as this pipeline's own results. No stored artefact carries
them; the stored runs give 0.005.)* Report
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
run is in. The panel, the 16 Hallmark TME sets, is fixed in the protocol.

**The library's defaults are not the whole registered configuration** (checked
2026-09-16 against all ten stored `config.json` files). `AuditConfig.n_boot`
defaults to 1,000 — every stored run used 1,000 (300 for the synthetic demo), and
the cohort scripts pass it explicitly; the default was 2,000 until 2026-09-24, and
a test now anchors it to the frozen runs' recorded value.
The other differences are per-run choices each script states: ComBat is off and
`stage_col` is empty for the four pan-cancer runs, and the demo uses 40 null sets
and a 100-patient floor. Call `run_audit` directly with a bare `AuditConfig()`
and, since 2026-09-24, you get the registered 1,000-draw bootstrap; the other
per-run choices above still have to be passed.

## Known limitations

- `combat_correct` has no empirical-Bayes shrinkage; confirm against `neurocombat`.
- Tiling, magnification and pooling inside the embeddings release are the
  release's own and are not re-derived here. The parquet stores 14 layers × 768;
  the final layer (`layer=-1`) is used, a choice fixed before the frozen runs and
  not varied.
- Column-name detection in `load_ancestry` is a guess recorded in provenance. Check it.
- `score_mean_z` is cohort-relative. Scoring TCGA and CPTAC separately silently redefines
  the target. Fix the scaling policy once.
- The demo table above is synthetic. Every other number in this repository comes
  from TCGA data. (Until 2026-09-17 this line said no real data had been
  analysed, a leftover from before the first real run.)

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
bash reproduce.sh --check              # tests + figures from stored results, 8 m 29 s (measured 2026-09-17, needs data/)
bash reproduce.sh                      # full: fetch, build interim, both cohorts, figures -- 8,751 s
                                       # (2 h 26 m) from an empty data/ on HPC4, 8 CPUs, measured
                                       # 2026-09-17 (job 130192); both cohorts reproduced bit for bit
```

**On a cluster, fetch on a login node first.** On Einstein HPC4 on 2026-09-17 a
compute node could not reach `data.broadinstitute.org` (TLS handshake timeout,
job 130180) although the login node could. Running
`python scripts/03_fetch_signatures.py && python scripts/01_fetch_data.py --all`
on the login node downloaded and verified every raw input in 37 s; the batch
job's `bash reproduce.sh` then finds them present, re-verifies them and needs no
network.

`requirements.txt` carries lower bounds for ordinary use; `requirements-lock.txt`
is what reproduces the manuscript numbers. The pin is not ceremony — sklearn has
changed `RidgeCV`'s alpha-selection defaults between minor versions, and the ISI
depends on the selected penalty.

**There is no continuous integration.** No workflow runs on push, and the
repository has no CI configuration. Every check is run by hand, or by
`reproduce.sh --check`, which chains them: the version check, provenance, the
split-determinism check, the sort audit and its self-test, the suite, the
number checker and its self-test, and coverage. The platform-dependence checks in
particular (`24_split_determinism.py --check`, `26_sort_audit.py`,
`--library-order`) mean something only when run on each machine that produces
numbers. Their per-machine records are in `results/` and are asserted against
each other by `19_check_numbers.py`.

### What runs from the deposit alone

**Read this before concluding anything is broken.** What is published at
`immune-specificity-index` is a *curated snapshot*, not this working repository:
the raw inputs under `data/` are public but large and are not redistributed
here, several multi-megabyte internal diagnostic probes are left out, and the
working documents (drafts, checklists, the poster) are not deliverables. Some
scripts shipped in the deposit therefore have no input there. **That is
curation, not corruption.**

Every row below was MEASURED by running the command inside a freshly staged
snapshot with no `data/` reachable — not inferred from reading the code. The
pytest and provenance rows were re-measured on 2026-09-17 in session 50 at the
187-test size. Every other row's exit code is now asserted by
`test_the_deposit_table_states_each_commands_real_exit_code`, which stages a
snapshot and runs each command; only the pytest row stays hand-measured,
because running the suite inside a test of the suite would recurse.
**Re-measure the pytest row whenever the suite changes**: its three numbers
must sum to the suite's size, and that arithmetic is what caught the previous
figures — which summed to 117 against a 122-test suite, and then to 125 against
a 152-test suite four sessions later, undetected for exactly as long as nobody
did the addition.

| command | exit | in the deposit |
|---|---|---|
| `python3 scripts/10_make_figures.py` | **0** | all five figures render |
| `python3 scripts/25_build_supplementary.py` | **0** | all ten tables build |
| `python3 scripts/20_check_versions.py` | **0** | |
| `python3 scripts/21_provenance_manifest.py` | **0** | 216 of 220 verified; 4 not staged, and it says so (re-measured 2026-09-17, session 50) |
| `python3 scripts/23_export_environment.py` | **0** | |
| `python3 -m pytest tests/ -q` | 1 | **42 failed, 135 passed, 10 skipped** — see below (re-measured 2026-09-17 at the 187-test size, 219 s) |
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

**The 42 test failures are all curation, and none of them is the analysis
code.** Classified by cause — mechanically, by pairing each failing test name
with its failure text, not by eye — they are:

| n | cause |
|---|---|
| 22 | exercise `22_build_public_snapshot.py`, which excludes *itself* by design and so is absent (staging fails, or the file is not found) |
| 8 | need the excluded `bisect_*.npz` platform probes, directly or because `19_check_numbers.py` refuses without them and says so in as many words |
| 8 | check the committed figures or `poster/poster.html`, none of which is deposited |
| 3 | are vacuity floors that cannot be met: they require more prose documents, or more `results/` entries, than the deposit stages |
| 1 | reads the NextGen abstract draft, which is not deposited |

The five buckets sum to 42, and that arithmetic is measured rather than
asserted: each failing test name is paired with the first error line of its own
failure text, read one by one (session 50), and the grouping above is what that
pairing returns. The split differs from session 49's 20/8/6/4/1 in how
checker-subprocess failures are assigned: read by their text, they are the
checker refusing without the probes, not unmet floors. The 135 that pass are the
pipeline's own tests. Two of them are new in session 47 and are worth noting here because they move the *passing* side rather than the failing one: the `ENVIRONMENT.md`-against-lockfile check and its can-fail read two files that ARE deposited, so unlike the blind sibling they compare — rather than skip — inside the deposit and on any other platform.

**Two of those causes did not exist when this section was first written, and it
said so in as many words: "There is no fourth cause."** It was true when
measured and false four sessions later, because the sentence recorded a
*measurement* in the grammar of a *law*. There are six causes today. A count
here is a dated observation; a claim that no further cause can arise is not
something a single run can support.

This row is re-measured rather than carried, and the record of it drifting is
the argument for doing so: the builder-dependent group was 6 at the
**117-test** size, 13 at the **125-test** size, 15 at the **155-test** size,
15 again at the **157-test** size, 20 at the **181-test** size and 22 at the
**187-test** size. (Those six
are hyphenated deliberately, so
the number checker does not read a past measurement as a current claim; the live
count is stated twice above.) Between the **125-test** and **155-test**
measurements the row went unre-measured for four sessions and drifted by 15
failures. The **157-test** re-measurement, taken one session later, moved the
skip count from 6 to 8 and left the failures and passes untouched — the two
tests added that session read the committed `.tsv` tables, which are not staged,
and skip on that condition rather than failing. The **181-test**
re-measurement (2026-09-17, session 49, after two dozen more test cases) moved the
failures from 39 to 39 — the same total over different causes, because the
builder group grew by five while the probe group shrank by five: the checker
now names the two excluded probes and refuses once rather than failing test by
test. Passes moved 116 to 132 and skips 8 to 10. **Re-measuring promptly is what
makes a two-line change legible as a two-line change.** A deposit
reader running the suite sees red on an intact deposit — that is a known and
accepted cost of curating, and it is written down here rather than left as a
surprise.

### The lock file is necessary but NOT sufficient — read this before rerunning

**Pinning the versions does not currently get you the same number on a different
machine.** Measured 2026-09-02, NSCLC ISI, same code, same seeds, same
`requirements-lock.txt` versions verified at runtime, and inputs hashed
byte-for-byte identical on both sides:

| | NSCLC ISI |
|---|---|
| macOS 26.5.2 / arm64, Python 3.13.9 | **0.318240180054566** ← the frozen results |
| Linux / x86_64, Python 3.12.14 | **0.2964** |

The gap is 0.0218. It does not change the sign, the significance, or the fact
that 16/16 signatures beat their null — it changes the second decimal place,
and this README quotes four. It is 2.5 times the standard deviation of partition
choice over six partitions (0.0089), a variability the study already widens its
interval for; 2.457 to full precision (`results/sort_fix_fold_change.json`).

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
> than amplifies. The fix is one keyword and **it was applied on 2026-09-05**:
> `splits.py` now passes `kind="stable"`, and the `--check` above reports the
> *shipped* partition as identical to the stable one. Applying it moved 216 of
> 944 patients to a different fold and would change every frozen number, so the
> author's decision was to fix the sort **and** keep the frozen results as the
> reported primary, reporting both numbers. See `../14-SCIENCE-AUDIT.md`
> §"A9 CLOSED" for the full two-machine table.

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
saying it is. *(Answered since, and recorded here rather than rewritten: no part
of the solve does. The "control" signatures diverge because the two machines
trained on different partitions -- the non-stable sort described above -- and
the solve itself attenuates, cond 18-77.)*

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
  once carried the label "PROVEN invariant by 16_alpha_margin"; it now says 11 of
  80 alphas differ across platforms, which is the measurement.
- `scripts/17_alpha_margin_fast.py` answers the margin question without
  `run_audit` (seconds of setup rather than hours). Its closing VERDICT asserted
  the alphas cannot differ until 2026-09-17; it now says why a wide margin does
  not imply that.
- `scripts/14_repro_probe.py` fingerprints every upstream stage; stages 1-6 are
  identical on both machines.
- `scripts/15_alpha_probe.py` prints the same 80-value digest as 16 and 17 and
  needs both machines to say anything -- which is exactly what it eventually
  said: `946033557a8ad08b` on HPC4 against `fc299c9c0f47466e` on macOS.

**So: the numbers in this README are the frozen macOS numbers, computed under
the old sort.** The code as shipped produces the corrected-ordering values
instead -- NSCLC 0.2922 and pan-TCGA 0.2968 -- and for NSCLC those agree between
macOS/arm64 and Linux/x86_64 to about 1e-16 (the pan-TCGA value was computed on
Linux). If you get them, you have reproduced the current code; 0.3182 needs the
old sort on macOS, and 0.2964 was the old sort on Linux. `reproduce.sh` compares
against the corrected-ordering runs for this reason.

`reproduce.sh` sets `PYTHONHASHSEED=0`. Python randomises str hashing per
process, so without it the gene pool is ordered differently every run and
`seed=0` does not reproduce.

### Headline results (the frozen primary runs, computed on 2026-08-19 and 2026-08-20 under the protocol tagged 2026-08-18)

> The protocol freeze is recorded by an annotated git tag, `prereg-2026-08-18`, **in the
> author's private working repository — not in this published snapshot, which is
> built as a fresh single-commit repository and therefore carries no tags.** The
> pre-registered protocol's *content* is here in full (`05-PRE-REGISTRATION.md`);
> its *timestamp* is not independently verifiable from this repository alone.
> Earlier wording put the tag name in this heading, which invited a reader to run
> `git tag`, find nothing, and reasonably conclude the claim was invented. The
> tag froze the protocol, not these numbers: the NSCLC and pan-TCGA runs below
> were computed the following two days, after the two changes the protocol's
> Appendix A records.

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

These sit one level ABOVE this pipeline, in the parent project directory. **Some
are part of this deposit and some are not**, and the `deposited` column below
says which — so a `../` link that fails to resolve can be told apart from a
broken one. The column is not prose: `test_deposited_sibling_documents_match_the_readme_table`
stages a snapshot and checks every row against it, so a document that gains or
loses deposit status cannot leave this table behind.

| File | deposited | What |
|---|---|---|
| `../04-CODE-AUDIT.md` | no | Five defects found and fixed in a full-code audit, with measurements |
| `../05-PRE-REGISTRATION.md` | yes | Protocol, frozen configuration, exploratory-vs-confirmatory split |
| `../07-ABSTRACT-DRAFT.md` | no | AACR 2027 abstract, 2,555/2,600 characters |
| `../09-PAPER-DRAFT.md` | yes | Manuscript draft, references verified against PubMed |
| `../SUPPLEMENTARY-NOTE-1.md` | yes | Supplementary Note 1: Limitation 8's forensic account of the ordering defect, lifted out of the manuscript 2026-09-24 |
| `../SUPPLEMENTARY-NOTE-2.md` | yes | Supplementary Note 2: the sensitivity analyses in full (partition, other analysis choices, two scorer checks, the NSCLC label-side control, the site control at smaller cohort sizes, ancestry), lifted out of the manuscript 2026-09-24 to bring its body under 5,000 words |
| `../11-SUBMISSION-PACK.md` | no | Category choice, submission checklist, poster plan |
| `../14-SCIENCE-AUDIT.md` | yes | Gaps between the protocol and the delivered results, including §A9 above |
| `../DATA-LICENSES.md` | yes | Where each input dataset comes from, under its own license; the code's MIT license is `../LICENSE` |

The rows marked `no` are planning documents kept only in the private working
repository. **`05-PRE-REGISTRATION.md` is deposited** — until 2026-09-08 this
section asserted the opposite of all three `yes` rows and told the reader to ask
the author for the protocol, which had been sitting beside this file since the
deposit was first pushed.

**One script here is not science.** `scripts/08_count_abstract.py` counts
characters against AACR's 2,600-character abstract limit under the Call for
Abstracts' own rule that spaces are free. It is conference tooling, kept in the
repository because it is what produced the 2,517 figure quoted above and a reader
should be able to check that number rather than take it on trust. It touches no
data and is not part of `reproduce.sh`.
