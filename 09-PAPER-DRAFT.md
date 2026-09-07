# Paper draft v1

**Target, revised 2026-08-19:** submit to ***Cancer Research*** (flagship,
IF 22.6, 2025 JCR) **first**, then transfer to *Cancer Research Communications*
(IF 3.3) via AACR's Manuscript Transfer Service if declined — the transfer carries
the reviews and reviewer identities, so aiming high costs weeks rather than a
review cycle. CRC accepts only the "transfer requires edits" pathway.

**The flagship's written gate**, which this paper must clear: technology/method
papers "must demonstrate the potential of the approach by providing new biological
insights to be considered". The biological insight here is *the prognostic content
of TME signatures is proliferative, not immune* — that must lead the framing, with
the disattenuation correction supporting it rather than heading it. **The observed
bar from the 2025-26 corpus is an independent external cohort**, which this work
does not yet have; see §7 of the pre-registration.

**Status:** complete draft against the 2026-08-18 analysis freeze
([`05-PRE-REGISTRATION.md`](05-PRE-REGISTRATION.md), git tag `prereg-2026-08-18`).
All numbers trace to `pipeline/results/nsclc_v3` and `pipeline/results/pancancer_v3`.
`pancancer_v2` is **identical on every immune-specificity number** — mean excess
0.291144, range 0.165–0.390, checked 2026-09-02 — and differs only in carrying
`null_draws.npz`, the persisted per-draw null correlations that make the
family-level rotation null and Figure 1B's distribution possible. An earlier
version of this line said v2; it was already reporting v3 results, since it quotes
the rotation null.
Figures are rendered by `pipeline/scripts/10_make_figures.py` into
`pipeline/results/figures/` as PNG at 300 dpi and vector PDF, with `--svg` and
`--tiff` (600 dpi) for the poster and for journal revision respectively.
**References verified against PubMed 2026-08-19** with PMIDs recorded.

---

# Proliferation, not immunity: H&E-inferred tumor microenvironment signatures are immune-specific but carry no prognostic advantage over random gene sets

## Abstract

*Journal format.* ***Cancer Research* abstracts are unstructured** — verified
across 40 published articles, which parse as an unlabeled body plus a
**Significance** statement, with no Background/Methods/Results/Conclusions
headings. Observed lengths: 131–264 words including Significance, with
Significance itself 18–40 words. The version below is written to that format and
measures **234 + 28 = 262 words** — inside the observed ceiling, though the
*official* cap is unverified because `aacrjournals.org` refuses automated
requests. The structured 461-word version it replaces is folded below and is the
right shape for the preprint, which has no such convention.

Deep learning models predict tumor microenvironment (TME) gene signatures from
routine H&E as low-cost surrogates for transcriptomic profiling,
but published controls cannot separate immune biology from two artifacts: a
dominant global expression axis onto which every large gene set loads, and the
higher reliability of curated modules relative to random gene sets. In 7,168 TCGA
patients across 31 cancer types and 944 non–small cell lung cancers, using
slide-level Prov-GigaPath embeddings, 16 MSigDB Hallmark TME signatures, and
preserved-site cross-validated ridge heads, we defined an
immune-specificity index: the excess image–signature correlation over 1,000 size-
and expression-matched random gene sets, after residualizing each score on the
within–cancer type first principal component of expression and disattenuating for
the residualized score's reliability. The index was 0.291 (95% CI, 0.269–0.313)
pan-cancer and 0.318 (0.261–0.376) in NSCLC, all 16 signatures exceeding their
null. Reliability dominated: removing the global axis reduced random
sets' Cronbach's α from 0.97 to 0.80 while curated signatures held above 0.97,
making the curated-versus-random gap 5.3- to 13-fold larger than on raw scores
and scaling it inversely with panel size — 0.084 at 160 genes to 0.453 at 10.
Against progression-free interval within cancer type, five signatures beat their
null — G2M checkpoint, E2F targets, angiogenesis, epithelial–mesenchymal
transition, and hypoxia — whereas no interferon, inflammatory, complement, or
allograft-rejection signature did (0 of 6 MSigDB-immune vs. 5 of 10 other
signatures; *P* = 0.031). H&E thus carries immune-specific information that is
not prognostic.

**Significance:** Image-based tumor microenvironment inference is advancing
toward clinical proposals faster than its controls. The immune signal is real,
but the prognostic value is carried by proliferation, not immunity.

<details>
<summary>Superseded structured abstract (461 words) — kept for the preprint</summary>

**Background.** Deep learning models predict tumor microenvironment (TME) gene
signatures from routine H&E and are increasingly proposed as low-cost surrogates
for transcriptomic profiling. Whether these predictions carry immune-*specific*
information, and whether that information is clinically useful, has not been
established, because the published controls cannot separate immune biology from
two artifacts: a dominant global expression axis that every large gene set loads
onto, and the higher measurement reliability of curated co-expressed modules
relative to random gene sets.

**Methods.** We analyzed 7,168 TCGA patients across 31 cancer types and, as a
second cohort, 944 non-small cell lung cancer patients, using public slide-level
Prov-GigaPath embeddings and 16 MSigDB Hallmark TME signatures. Ridge heads were
fitted under preserved-site cross-validation. We defined an immune-specificity
index (ISI): the excess, in Fisher z, of the image–signature correlation over
that of 1,000 size- and expression-matched random gene sets, after residualizing
every score on the within-cancer-type first principal component of expression and
disattenuating for reliability computed on the residualized score. Observed and
null heads were refit identically. The same null was run against progression-free
interval with concordance computed within cancer type.

**Results.** The ISI was 0.291 (95% CI 0.269–0.313) pan-cancer and 0.318
(0.261–0.376) in NSCLC, with all 16 signatures exceeding their null in both
cohorts. Residualizing the global axis cost little. Reliability mattered
substantially more than previously appreciated: random gene sets are internally
consistent largely because their genes load on the global axis, and removing that
axis reduced their Cronbach's α from 0.97 to 0.80 while curated signatures held
above 0.97. The resulting curated-versus-random reliability gap was 5.3-fold
(NSCLC) to 13-fold (pan-cancer) larger than the gap computed on raw scores, and
it scaled inversely with panel size; subsampling signatures to fixed size, the
gap rose monotonically from 0.084 at 160 genes to 0.453 at 10 (Spearman
ρ = −1.000, p < 0.0001). Against
outcome, five signatures beat their random-set null — G2M checkpoint, E2F
targets, angiogenesis, epithelial–mesenchymal transition and hypoxia — while no
interferon, inflammatory, complement or allograft-rejection signature did (0 of 6
MSigDB-immune versus 5 of 10 other signatures, Mann–Whitney p = 0.031). Tissue
source site was recoverable from the embeddings at median AUROC 0.998.

**Conclusions.** H&E carries genuinely immune-specific TME information beyond
tissue composition. That information is not prognostic: the outcome signal
carried by TME signatures is proliferative and stromal, and image-inferred immune
scores should not be used as prognostic biomarkers. Reliability disattenuation is
mandatory when comparing curated gene sets against random ones, and its magnitude
depends on panel width, leaving small clinical panels most exposed.

**Significance.** Image-based TME inference is advancing toward clinical
proposals faster than its controls. We show the signal is real but that its
prognostic value is carried by proliferation rather than immunity, and we give a
correction that any gene-set null requires.

</details>

---

## Introduction

Whole-slide H&E images are the cheapest, most universally available data in
oncology, and a substantial literature now predicts bulk transcriptomic
signatures from them. HE2RNA (1) predicted gene
expression across 28 TCGA cancer types; HistoTME and successors predict curated
TME signatures directly; foundation models such as Prov-GigaPath (2),
UNI (3) and Virchow (4) have raised slide-level representation quality substantially. The
implicit promise is a low-cost surrogate for transcriptomic immune profiling, and
by extension for immunotherapy stratification.

Two objections stand between that promise and the evidence, and neither is
addressed by the controls in current use.

**The global-axis objection.** Bulk expression has a dominant axis — tumor
purity, stromal fraction, proliferation, tissue of origin — that essentially
every gene set of appreciable size loads onto. A model that reads tissue
composition faithfully from an image will therefore predict an immune signature
well, without carrying any immune-specific information. Crucially, this
hypothesis and the "genuine immune signal" hypothesis predict the *same* raw
correlation. Comparing a curated signature against a random gene set does not
separate them, because both load on the same axis. HE2RNA ran precisely such a
random-list null — 10,000 size-matched lists per pathway per cancer type — and
reported that 75% and 86% of cancer types had B-cell and T-cell signatures
significantly better predicted than random. That result is consistent with either
hypothesis.

**The reliability objection.** Curated signatures are, by construction,
co-expressed modules: high mean inter-gene correlation, high internal
consistency. Random gene sets of the same size are not. A more reliable target is
measured with less error and is therefore more predictable from *any* predictor,
entirely independently of biology. Comparing the two without disattenuation is
biased in favor of the curated set.

We address both. The signature — not the image prediction — is residualized on
the global axis, and the contrast is disattenuated using the reliability of the
residualized score. We then ask the question that matters clinically: does any of
this predict outcome better than a random gene set does? That is Venet's question
(5), which showed that most random
gene sets of >100 genes significantly predict breast cancer outcome and that 28
of 47 published signatures did no better than same-size random sets. It has never
been asked of the TME signature family, and never of image-inferred scores.

## Materials and Methods

### Cohorts and data

All data are public and permissively licensed; no controlled-access data was
used. Slide-level Prov-GigaPath embeddings (2) were
obtained from
`seandavis/tcga_provgigapath_embeddings` (HuggingFace, CC-BY-4.0), which ships
TCGA-CDR clinical annotation merged. The parquet stores 14 slide-encoder layers
of 768 dimensions per slide; the **final layer** was used and this choice was
fixed in advance. Only primary solid tumor slides (sample code 01) were
retained, of which 99.5% are diagnostic (DX) formalin-fixed slides. Slides were
averaged to one vector per patient before any splitting.

Expression is Xena TOIL RSEM log2(TPM+1), mapped from Ensembl to HUGO symbols via
HGNC (42,356 pairs), giving 41,046 genes. Tumor purity is PanCanAtlas ABSOLUTE;
CPE and ESTIMATE were deliberately avoided because they are expression-derived
and would make the purity control circular. Outcome is progression-free interval
from TCGA-CDR, which Liu and colleagues (6) recommend for all but 4 of 33 TCGA
types.

The **pan-cancer cohort** is 7,168 patients across 31 cancer types and 619 tissue
source sites (2,442 PFI events). The **NSCLC cohort** is 944 patients (LUAD 471,
LUSC 473; 328 events).

Signatures are the 16 tumor-microenvironment-relevant sets of MSigDB Hallmark
v2024.1 (CC-BY-4.0). BostonGene/MFP definitions were considered and rejected:
their license restricts use to academic and non-profit purposes and defines
"Source Form" to include data, which is incompatible with this project's
permissive-license constraint.

### Scoring, splits and models

Signatures were scored as the mean of per-gene z-scores across the cohort, and
the identical scorer was applied to observed and null sets — mixing scorers is a
scale mismatch rather than a null (measured sd ratio ≈ 3× between mean-z and
ssGSEA).

Cross-validation is **preserved-site** (7): whole
tissue source sites are assigned to folds so that no site appears in more than
one fold, with site disjointness asserted at runtime. Sites are assigned
largest-first with joint size and outcome balancing, and a partition producing an
empty or badly undersized fold is refused rather than returned. Models are ridge
heads with the penalty tuned by generalized cross-validation **inside each
training fold**; imputation and standardization are likewise fitted on training
folds only.

Two orderings in this pipeline are made explicit because a seed alone does not
determine them. Sites of equal size are common — 51 of the 68 NSCLC sites share
a size with another site — and the largest-first ordering must therefore break
ties. The rank table underlying the ssGSEA scorer must likewise break ties among
genes of equal expression, which on log expression is the normal case rather
than the exception, since every gene resting on the floor shares a rank. Both
are now settled by a **stable** sort, so ties retain input order and the result
is fixed by the data rather than by the sorting implementation. This is a
deliberate pin rather than a correctness claim: any tie convention yields a
valid partition and a valid ssGSEA, but only a pinned one yields the same answer
on two machines. The results reported as primary were produced before this pin
and are reported unchanged, with the re-run under the pinned convention given
alongside them in Limitations.

### The immune-specificity index

For signature *s*:

    ISI_s = z[ r(pred_s, signature_s | axis) / sqrt(α_s) ]
            − mean_j z[ r(pred_sj, random_sj | axis) / sqrt(α_sj) ]

where *z* is Fisher's transform, the axis is the first principal component of
log expression computed **within cancer type**, and α is Cronbach's alpha.

Three details are load-bearing.

**Reliability is computed on the residualized score, not the raw score.**
Cronbach's closed form α = k·ρ̄/(1+(k−1)ρ̄) assumes unit-variance items and is
therefore valid only for a raw mean-z composite. Because the ISI disattenuates a
*residualized* correlation, we use the general form
α = k/(k−1)·(1 − Σvar(itemⱼ)/var(Σitemⱼ)). Since the composite is the mean of k
items, Σitemⱼ = k·score, so only per-gene residual variances are required; these
are computed once for all 41,046 genes.

**Observed and null heads are refit identically.** Each null draw's head is
refitted on its residualized random set; the observed head is likewise refitted
on the residualized signature. Correlating a residualized target against a head
trained on the raw signature would not be a like-for-like comparison.

**Draws whose disattenuation is undefined are dropped, not clipped.** When
|r/√α| ≥ 1, sampling error in α has pushed the corrected value outside the valid
range; Fisher z of 1 is unbounded and a handful of such draws would dominate a
null mean. The count of dropped draws is reported (0 of 16,000 in both cohorts).

Nulls are 1,000 gene sets per signature, matched on both size and the
signature's expression-level histogram. Intervals are patient-clustered
bootstraps: one patient set is resampled per draw and shared across signatures,
preserving between-signature correlation; complete cases are used and the number
of incomplete patients is reported.

### Outcome arm

Signature and null scores enter Harrell's C directly as continuous risk scores.
There is no dichotomisation path in the code, and the reduction from model to
scalar was fixed in advance, both to avoid documented inflation of type I error
(8) and C-hacking (9).

**Concordance is computed within cancer type.** Progression-free interval differs
greatly between tumor types, and cancer type alone reaches C = 0.676 pan-cancer;
an unstratified concordance therefore largely measures which cancer a sample is,
and any tissue-correlated score inherits that. Unstratified, 94% of signatures
"beat" their null; stratified, 31% do. Only the stratified figures are reported.

Every null result is reported with its **minimum detectable effect** at 80%
power, because a null outcome result is otherwise uninterpretable.

**Comparing immune against non-immune signatures.** The paper's central
biological claim — that the prognostic content of the TME family is
proliferative rather than immune — rests on a comparison between two groups of
signatures, and that comparison needs its own stated method.

Signatures are assigned to process categories using **MSigDB's own
categories** (10), not a grouping of our making; the assignment is a literal
table in `pipeline/scripts/10_make_figures.py` and is the single source of
truth that the numerical checker also reads. Six signatures are MSigDB-immune
(allograft rejection, coagulation, complement, inflammatory response, and
interferon α and γ). Note that **coagulation and complement are MSigDB-immune
and are counted as such**, which works against the contrast rather than for it.
IL6/JAK/STAT3 signalling is **excluded** from the immune group because MSigDB's
gene-set page states no process category for it; since its excess is negative,
including it would strengthen the contrast, so excluding it is the conservative
choice.

Two tests are reported on the axis-residualized outcome excess: a
**Mann–Whitney U test** on the continuous excess, and **Fisher's exact test** on
the counts beating their null.

**Both are reported two-sided.** The Mann–Whitney *P* is 0.031 (0.031219) and
Fisher's exact *P* is 0.093 (0.093407). Neither test was pre-registered as
directional, so no one-sided value is reported anywhere in this paper, including
in figure captions.

> **Provenance of this convention, recorded because it was not always so.**
> Until 2026-09-06 the paper quoted "*P* = 0.016; *P* = 0.093" in a single
> parenthesis. Recomputed from `outcome_arm.csv`, the first was the **one-sided**
> Mann–Whitney (0.015609) and the second the **two-sided** Fisher (0.093407) —
> each correct for the test it named, but a mixed convention, stated nowhere, in
> one parenthesis. It was found by recomputing the numbers in order to write this
> Methods paragraph, which did not previously exist. All four values (both tests,
> both tails, with and without IL6/JAK/STAT3 counted as immune) are now
> mechanically checked by `19_check_numbers.py`, so the convention cannot drift
> back silently. The contrast survives the change: *P* = 0.031.

### Reproducibility

104 automated tests; the analysis is deterministic given seed 0 **on a fixed
platform**. The configuration is written to `config.json` on every run and the
code asserts at save time that the computed estimand matches the pre-registered
string. Each run additionally records the cohort it was fitted on — the patient
count and a SHA-256 of the sorted patient identifiers — so that a downstream
analysis claiming to describe a run can be checked against the patients the run
actually used rather than trusted; the frozen runs predate this and do not carry
it. Package versions are pinned exactly in `requirements-lock.txt`.

Determinism does not currently extend across hardware, and this is stated rather
than discovered by a reader. Re-running the NSCLC primary on Linux/x86_64 with the
same seeds, the same pinned versions and inputs verified byte-identical returns
ISI 0.2964 against 0.3182 on the macOS/arm64 platform the results were frozen on.
All values reported here are the latter. The cause is now identified: the
fold-assignment routine sorts tied sites with a non-stable sort, so the two
architectures build different cross-validation partitions from identical inputs.
`scripts/24_split_determinism.py` reports the partition hash on any machine and,
with `--check`, asserts the stable-sort partition against the value measured on
both. Under a stable sort the partition is identical across platforms and the
residual numerical difference is 7 × 10⁻¹⁶. See Limitations 8.

## Results

### The image predicts TME signatures, and preserved-site splitting costs little

Under preserved-site cross-validation the
median image–signature correlation was 0.630 pan-cancer and 0.400 in NSCLC. The difference from a site-agnostic
patient-level split was small (median Δr = +0.020 and +0.023; paired calibration
loss Δ-MAE = 0.0088 [0.0078, 0.0098] and 0.0054 [0.0024, 0.0085]), indicating
that the models are not merely site detectors even though site is highly
recoverable (below).

**A necessary caveat.** Pan-cancer, a covariate-only baseline of site, cancer
type, purity and stage reached median r = 0.655, and the embedding added
**−0.000** over it. In NSCLC, where cancer type carries far less information, the
embedding added +0.031. The pan-cancer raw correlation is therefore substantially
tissue of origin, and must not be read as image performance. This is precisely
why the ISI is defined on within-type residualized scores.

### The immune-specific excess is positive and replicates

The ISI was **0.291 (95% CI 0.269–0.313)** pan-cancer and **0.318 (0.261–0.376)**
in NSCLC. All 16 signatures exceeded their own null in both cohorts (per-signature
excess 0.165–0.390 and 0.148–0.440; per signature in Supplementary Table
S3). Median residualized correlations were 0.394
and 0.372 against null means of 0.106 and 0.051.

The larger cohort gives a tighter interval and an overlapping estimate: this is a
replication across 31 diseases, not a single-disease finding.

Residualizing the global axis was cheap, because the image barely tracks it:
median image–axis correlation was 0.067 pan-cancer and **−0.023** in NSCLC. The
composition hypothesis in its strong form is therefore not what is happening —
the image is not simply reading the dominant expression axis.

> **Figure 1.** (A) Per-signature immune-specific excess with 95% CIs in both
> cohorts, ordered by pan-cancer excess. (B) Pan-TCGA observed axis-residualized
> correlation against the random-set null, drawn as **the actual null
> distribution** — the 1,000 per-signature draws that `AuditResult.save()`
> persists to `null_draws.npz`, shown as violins, not a summary band.

### Reliability, not composition, is the correction that matters

Random gene sets are internally consistent — median α 0.974 pan-cancer — but
almost entirely because their genes load on the global axis. Removing that axis
dropped the null's α to 0.801. Curated signatures barely moved:
their α fell from 0.984 to 0.976. Per-signature α, residualized and raw, with
the curated-minus-random gap, is Supplementary Table S4.

The consequence is that the reliability gap computed the conventional way, on raw
scores, is badly understated: +0.014 pan-cancer and +0.032 NSCLC, against
**+0.184 and +0.167** when computed on the residualized score that the estimand
actually disattenuates — a 13-fold and 5.3-fold difference.

The gap scales inversely with panel size (Spearman ρ = −0.50, p = 0.049
pan-cancer; ρ = −0.59, p = 0.015 in NSCLC): it averages +0.141 across the
200-gene sets and reaches +0.420 for the smallest, 36-gene angiogenesis set.

Because real clinical panels are smaller than any Hallmark set, we tested the
scaling *inside* the measured range rather than extrapolating below it, by
subsampling each signature to fixed k and comparing against random sets of the
same size. The gap rises monotonically: 0.084 at k=160, 0.141 at 80, 0.214 at 40,
0.290 at 20 and **0.453 at k=10** (Spearman ρ = −1.000, p < 0.0001) — 5.4-fold
across the range; the full bracket is Supplementary Table S2. These are raw-score alphas and therefore a conservative lower
bound on the residualized gap the index actually corrects with.

Hallmark sets are large, so
this cohort sits in the regime where the correction is mild. **Clinical panels of
5–20 genes sit at the opposite end**, where an uncorrected curated-versus-random
contrast would be dominated by measurement reliability rather than biology.

> **Figure 2.** (A) Curated-minus-random Cronbach's α against gene-set size,
> both cohorts, with raw-score α (open markers) and residualized α (filled)
> overlaid. (B) Median α of random sets and of curated signatures, before and
> after removing the global axis — the random sets collapse, the curated ones do
> not.

### The index depends on how the signature is scored

The reliability correction above is applied to a mean-of-z-scores composite. To
test whether the index survives a different scoring rule, we recomputed the whole
NSCLC arm under single-sample GSEA — a rank-walk statistic that shares no
arithmetic with the mean-z composite — holding the null sets, the splits, the
bootstrap draws and every seed fixed. All four arms (two scorers × corrected and
uncorrected) ran on one machine. As an acceptance check, the `mean_z`
disattenuated arm reproduces the frozen primary to full double precision, CI
bounds included (+0.3182 [0.2614, 0.3763], 16/16 beating null), confirming the
harness perturbs nothing but the scorer.

Three findings follow, in decreasing order of authority.

**First, on the like-for-like uncorrected contrast the effect survives the change
of scorer at roughly 56% of its magnitude**: +0.180 [0.139, 0.216] under ssGSEA
against +0.319 [0.266, 0.373] under mean-z, with 15 of 16 signatures still
beating their null (16/16 under mean-z). Rank-order agreement between the two
scorers is moderate across signatures (Pearson r = 0.633, Spearman ρ = 0.582),
and the residualized observed correlations themselves agree closely (r = 0.766,
mean |difference| 0.049). These ssGSEA figures were computed before the
tie-ordering defect described in Limitation 8 was found, and the rank table they
rest on is one of the two places it occurred. That arm has since been recomputed
under the pinned ordering, and the defect does not propagate to this contrast:
+0.174 [0.133, 0.210], a shift of 0.006, or 16% of the pre-fix interval's own
half-width, with the same 15 of 16 signatures beating their null and
essentially the same fraction of the mean-z magnitude (54% against 56%). The
movement the defect produces in the underlying scores — up to 0.61 of a score
standard deviation — therefore averages out of the index rather than
accumulating in it. Both values are reported, as for the primary. The mean-z
figures beside them are unaffected in any case, as that scorer never touches the
rank table.
This arm answers "does the scorer change the answer"
and it does not change it qualitatively — but it is the *uncorrected* contrast,
which had a 100% false-positive rate in our own simulation, and it therefore
cannot be quoted as an immune-specificity index.

**Second, the registered estimand is undefined under ssGSEA.** All 1,000 null
draws hit the disattenuation guard for all 16 signatures, so no corrected ssGSEA
index exists to report. This is a property of the reliability estimator, not of
ssGSEA: Cronbach's α is defined here for a composite that is the mean of k
z-scored genes, and an ssGSEA score is a rank-walk statistic with roughly 20-fold
smaller dispersion, so α collapses toward zero and every draw is dropped. The
observed signatures remain highly reliable under ssGSEA (split-half 0.934); it is
the *random* sets whose reliability collapses, from 0.799 under mean-z to 0.258.

**Third, and least authoritative, a scorer-agnostic reconstruction reverses the
sign** (−0.1333 [−0.1831, −0.0835], 0/16 beating null, against +0.3177 [0.2822,
0.3533] for the same estimator under mean-z — quoted to four decimals because at
three it is not distinguishable from the frozen primary, which it is not: the
reconstruction replaces the patient-clustered bootstrap with an across-signature
standard error). The mechanism is visible in the reliability
numbers: with null reliability 0.258 against 0.934 observed, disattenuation
multiplies the null by 1/√0.258 ≈ 1.97 and the observed by only 1.03. Dividing by
the square root of a near-zero reliability is precisely where disattenuation
becomes numerically violent, and this estimate also replaces the registered
patient-clustered bootstrap with an across-signature standard error. We report it
for completeness and do not treat it as an estimate of the index.

The defensible conclusion is not that the index is negative under ssGSEA, but
that **its magnitude — and under the scorer-agnostic estimator even its sign —
depends on the scoring rule, because rank-walk scores of random gene sets are
barely reliable enough to disattenuate against.**

### Prognostically, the signal is proliferation and stroma — not immunity

Pan-cancer, with concordance computed within cancer type, five signatures beat
their random-set null: **G2M checkpoint** (excess 0.069), **E2F targets**
(0.062), **angiogenesis** (0.051), **epithelial–mesenchymal transition** (0.038)
and **hypoxia** (0.030). No interferon, inflammatory, complement or
allograft-rejection signature did. Using MSigDB's own process categories
(10), 0 of 6 immune signatures beat their null
against 5 of 10 others (Mann–Whitney p = 0.031; Fisher p = 0.093).

In NSCLC no signature beat its null. **That is an underpowered null, not a
negative result**: the NSCLC minimum detectable effect is 0.046 in C-index, while
the real pan-cancer advantage is approximately 0.03 (pan-cancer MDE 0.017). The
two cohorts are consistent; the smaller one simply cannot see an effect of this
size. We report the NSCLC result only with its MDE attached.

This recapitulates Venet's decisive control (5) in a new setting. Venet showed that
adjusting for a proliferation metagene abrogated nearly all published
signature–outcome associations. Here, the signatures that retain prognostic value
over a random-set null *are* the proliferation and stromal programs, and the
immune programs do not.

> **Figure 3.** Outcome excess over the random-set null per signature,
> pan-TCGA, within cancer type, colored by MSigDB process category, with the
> minimum detectable effect shaded and NSCLC overlaid as crosses. Black-edged
> points beat their null. IL6/JAK/STAT3 is plotted as *signaling* rather than
> immune because MSigDB's gene-set page states no process category for it;
> counting it as immune would strengthen the contrast (p = 0.031 → 0.016), so
> excluding it is the conservative choice.

### Multiplicity, panel redundancy and the label-side control

A family-level permutation test that respects between-signature correlation —
every signature is scored against the same images and the same axis, so a naive
sigma/sqrt(m) band around the family mean is far too narrow — gives an observed
family mean correlation of 0.386 against a null family mean of 0.110
(95% range 0.098–0.122), p = 0.001 at 1,000 draws.

Benjamini–Hochberg correction at q = 0.05 across the 16-signature family changed
neither headline: 16 of 16 signatures still exceeded their null on the index, and
10 of 32 comparisons still did so on outcome — including hypoxia, which had the
narrowest margin (excess_lo = +0.0069). The effective number of independent tests
(11) is 5.00 for the 16 signatures, reported descriptively; dividing alpha by
it *and* applying Benjamini–Hochberg would be a category error and is not done.

Regressing each signature's ground truth on tissue source site **from RNA alone,
with no image**, site appears to explain a median 44% of the label pan-cancer.
That figure is almost entirely an artifact of two things: source site largely
determines which diseases a center contributes, and dummy-coding 202 sites
inflates R² by construction. Conditioning on cancer type reduces it to 5.5%, and
calibrating against site labels permuted *within* cancer type gives 4.1%
(p = 0.010). Plate within site — batch with no plausible biological reading —
explains 0.000 (per signature in Supplementary Table S6). So the label-side
site artifact is small once disease is
accounted for, and is not technical; the site information the embeddings carry is
not mirrored by a comparable artifact in the targets.

### Sensitivity to the choice of site-to-fold partition

Preserved-site cross-validation admits many valid assignments of sites to folds,
and the analysis uses one. We re-ran the primary endpoint under six assignments
in each cohort, varying only the partition with the null gene sets and bootstrap
draws held fixed (per partition, both cohorts, in Supplementary Table S8).

In NSCLC the between-partition sd is 0.0089, against a within-partition
bootstrap standard error of 0.0294. Combining the two as independent sources
widens the interval from [0.2542, 0.3695] to [0.2516, 0.3720]: in NSCLC the
reported interval is 4.5% too narrow, and partition choice accounts for 8.4%
of total variance in NSCLC.

In pan-cancer the between-partition sd is 0.0065, against a
bootstrap standard error of 0.0111, over the five assignments of six that the
analysis could evaluate. The same combination widens [0.2717, 0.3152] to
[0.2683, 0.3187]: the reported interval is 16.0% too narrow, and
partition choice accounts for 25.7% of total variance. The conclusion is
unchanged in both cohorts — the honest interval still excludes zero by a
wide margin.

The comparison between the two cohorts is the more informative result, and it
runs against the intuition that a larger cohort is a safer one. Partition choice
is a *larger* share of total uncertainty pan-cancer (25.7%) than in NSCLC
(8.4%), even though the pan-cancer cohort is 7.6 times larger. The reason is
that the two variance components scale differently: the bootstrap standard error
falls with sample size, from 0.0294 to 0.0111, while the between-partition
standard deviation barely moves, from 0.0089 to 0.0065. Partition variance
behaves as a floor that resampling patients cannot lower, because it comes from
the design rather than from sampling noise. Reporting a bootstrap interval alone
therefore becomes *more* misleading, not less, as cohorts grow — which is the
opposite of how such intervals are usually read.

Two further observations. The sixth pan-cancer assignment was refused by the
degeneracy guard when the residualization SVD failed to converge, so the
estimand is not computable for every valid partition, and a run that silently
substituted another partition would have hidden that. And this sensitivity is
specific to the index: for the split-robustness contrast Δr, which is a
difference *between* two split schemes rather than a contrast computed within
one, partition choice matters considerably more.

### Negative controls

Tissue source site was recoverable from the embeddings at median one-vs-rest
AUROC **0.998** pan-cancer and 0.992 in NSCLC, with 100% of evaluable sites above
0.9 in both (per site, both cohorts, plus the ComBat arm, in Supplementary
Table S5). We report this raw figure as the confounding measure and deliberately
do **not** report a post-ComBat AUROC: per-site centring imposes a within-site
zero-sum constraint that makes the site classifier systematically anti-predictive
(measured 0.005), so the corrected number is uninterpretable in either direction.
The deployment-relevant statement is that under preserved-site folds ComBat is
estimable for **0%** of held-out samples by construction — a model meeting a new
hospital is in exactly that position.

> **Figure 4.** (A) Distribution of one-vs-rest site AUROC across evaluable
> sites in both cohorts. (B) Median correlation achieved by a covariate-only
> baseline against the image embedding — pan-cancer the two are
> indistinguishable.

The image prediction's partial correlation with the signature, controlling for
ABSOLUTE purity, was 0.547 pan-cancer (against 0.627 unadjusted), so purity
explains part but not most of the association. Incremental adjusted R² of the
image over purity, type, site and stage was 0.036 and 0.054.

### Ancestry (supplementary)

A pan-cancer ancestry arm (6,580 patients; EUR 5,312, AFR 594, ASIAN 502, AMR
172) is reported as Supplementary Table S1 with a negative headline: **the
contrast is not identifiable from tissue source site**. Only 0.6–1.0% of sites
carry at least 10 patients of both the reference and any comparison group
(ancestry × site Cramér's V = 0.537 against a permutation null with 95th
percentile 0.088). Ancestry calls are the published UCSF consensus set
(12). Any TCGA study reporting ancestry
subgroup differences in an
image model faces this constraint; most do not check it.

## Discussion

The headline is two-sided, and both sides matter.

**H&E carries real immune-specific information.** After conditioning the
signature on the dominant expression axis and correcting for the measurement
reliability of curated modules, curated TME signatures remain substantially more
image-predictable than size- and expression-matched random gene sets, in 31
cancer types. This is a stronger result for the modality than the composition
hypothesis predicts, and it is obtained under site-disjoint validation.

**That information is not prognostic.** The prognostic content of the Hallmark
TME family, measured against the only fair comparator, is proliferative and
stromal. Immune programs — interferon α and γ, inflammatory response,
complement, allograft rejection — do not beat random gene sets on
progression-free interval within cancer type. Since these are exactly the
programs that image-based TME inference is proposed to deliver, the practical
implication is direct: **image-inferred immune scores should not be used as
prognostic biomarkers**, not because the image reads them poorly, but because the
underlying transcriptomic scores are not prognostic once compared fairly.

**A methodological point that generalises.** Any comparison of a curated gene set
against a random one is confounded by reliability, and the size of that
confounding depends on how much of the sets' internal consistency comes from a
shared global axis. Because that dependence is inverse in panel size, the
smallest panels — which are the ones nearest clinical deployment — are the most
exposed. Reporting a raw-score α understates the problem by up to an order of
magnitude.

### Limitations

1. **Exploratory.** The analysis preceded the protocol. Analysis choices were
   revised after seeing results, each to fix an identified defect, and the
   revision history is public. External validation is pre-registered but not
   run; see item 10 for why.
2. **One consortium, one embedding model.** All results use TCGA and
   Prov-GigaPath. Whether they hold for other foundation models is untested.
3. **Site is almost perfectly recoverable** from the embeddings (AUROC 0.998).
   Preserved-site cross-validation makes the predictions honest — no test site is
   seen in training — but it does not remove site information from the features,
   and no available harmonization is estimable for an unseen site.
4. **Pan-cancer the embedding adds nothing over covariates** (Δr = −0.000 over
   site + type + purity + stage; the full nested decomposition is
   Supplementary Table S7). The pan-cancer ISI is meaningful because it is
   defined within cancer type, but the raw pan-cancer correlation is not a
   measure of image performance.
5. **Slide-level embeddings, linear heads.** No attention-based multiple-instance
   learning was fitted. This is deliberate — it removes a tuning surface and
   matches standard linear-probe practice — but a stronger head might change the
   magnitude, though not the null comparison, which is symmetric by construction.
6. **The outcome arm maps C-index onto Somers' D** to reuse the Fisher-z
   machinery, including a 1/√(n−3) standard error derived for a correlation
   rather than a concordance statistic. The outcome intervals should be read as
   indicative.
7. **Per-signature nulls share a seed**, so two signatures of identical size and
   expression profile draw identical random sets. This is intentional for the
   family-level statistic but means the 16 nulls are not independent.
8. **The ISI is platform-dependent in its third decimal, and the cause is a
   non-deterministic fold assignment rather than arithmetic.** With seed,
   configuration and package versions pinned and the inputs verified
   byte-identical by hashing, the NSCLC ISI is 0.3182 on macOS/arm64 and 0.2964
   on Linux/x86_64. The fold-assignment routine orders tissue source sites
   largest-first with a sort that is not stable, and
   51 of the 68 NSCLC sites share a size with at least one other site; the order
   among tied sites is settled by the sort's internal tie-breaking, which
   differs between processor architectures. The two platforms therefore build
   *different cross-validation partitions* from identical inputs — with
   identical fold sizes, which is why the site-disjointness and fold-balance
   guards did not catch it. Hashing confirms it directly: the site labels, site
   sizes and embedding matrix hash identically on both machines while the fold
   assignment does not, and 216 of
   944 patients change fold. Because the training partitions differ, the
   per-fold ridge penalties differ too — 11 of 80, each by one grid step — so
   those flips are a consequence and not the cause. Substituting a stable sort
   makes the partition hash identical on both platforms and, on a synthetic
   cohort driven through the same estimator chain, reproduces the Linux value on
   macOS to 15 significant figures, a residual of 7 × 10⁻¹⁶. Genuine
   cross-platform floating-point noise in this pipeline is therefore six orders
   of magnitude smaller than the discrepancy long attributed to it, and the
   ridge solve is not implicated: the regularized system is well conditioned
   (condition number 18–77) and a relative 10⁻¹⁴ input perturbation moves the
   predictions by 6 × 10⁻¹⁵, so it attenuates rather than amplifies. Two earlier
   accounts of this limitation — one arguing from the width of the penalty
   selection margins that the penalties *could not* differ across platforms, the
   other attributing the divergence to amplification through an ill-conditioned
   ridge system — are contradicted by these measurements and are corrected here.
   The discrepancy is thus an instance of the partition-choice sensitivity
   already quantified above rather than a separate defect: 0.0218 is 2.4× the
   partition standard deviation of 0.0089. It does not alter the sign, the
   significance, or the finding that 16/16 signatures exceed their nulls, so no
   conclusion here turns on it; it does mean the fourth digit should not be
   over-read. **The stable sort is now the shipped default**, and a regression
   check fails if it is reverted; the frozen results retain the original sort so
   that they remain the runs that were pre-registered and reported.

   That separately reported re-run now exists. Re-running the pre-registered
   NSCLC configuration unchanged except for the two corrected orderings gives an
   ISI of 0.2922 [0.2355, 0.3527], against the pre-registered
   0.3182 [0.2614, 0.3763]. The point estimate moves by 0.0260, which is 2.9
   times the between-partition standard deviation reported above and falls
   inside the honest interval [0.2516, 0.3720] that the partition-variance
   analysis derived — so the correction moves the estimate by about as much as
   changing the partition does, which is what a fold-assignment defect should
   do, and no more. All 16 of 16 signatures still exceed their nulls and the
   interval still excludes zero. We report both numbers rather than replacing
   one with the other: the pre-registered value is what was registered and run,
   and the corrected value is what the fixed code produces, and a reader is
   entitled to see the difference rather than a single number chosen after the
   fact.

   A systematic audit of every other non-stable ordering in the pipeline, run
   after this cause was identified, found **one further instance and it is on
   the scoring path**: the single-sample GSEA rank table ordered tied genes by
   the same architecture-dependent rule. Ties there are not incidental — log
   expression has a floor, every gene resting on it shares an average rank, and
   all 944 samples carry ties, 35,394,448 of 38,747,424 rank entries in total.
   The two orderings produce different scores for 15,060 of 15,104
   patient–signature cells across all 16 signatures. The largest score
   difference between the two orderings is 0.0128, against a mean
   per-signature score standard deviation of 0.0211, and the
   ordering is confirmed architecture-dependent by hashing on both machines. The
   per-machine hashes for the fold assignment are Supplementary Table S9.
   This does not touch any result above that is computed from the
   mean-of-z-scores composite, which does not use the rank table. It does reach
   the ssGSEA sensitivity arm, which has now been recomputed under the pinned
   ordering: the uncorrected contrast moves from +0.180 to
   +0.174, 15 of 16 signatures beat their null either way, and the
   estimand stays undefined under ssGSEA for the reason given in Limitation 9.
   A 0.61-standard-deviation movement in the per-cell scores thus produces a
   0.006 movement in the index they aggregate to, so the defect is real at the
   score level and immaterial at the level of anything reported.
   Both orderings are now pinned. The remaining orderings were
   checked and are safe for stated reasons rather than by assumption: the
   Benjamini–Hochberg step-up assigns tied *p*-values the same *q* whichever
   order they take, and the site-AUROC table is consumed only through a median.
   One further order-dependent routine, a per-site summary, is not called
   anywhere in the pipeline; it has nonetheless been pinned, on the view that a
   latent instance of a defect this analysis has already been bitten by twice
   should not be left for whoever calls it next. With that, the audit reports
   no unpinned order-dependent operation anywhere in the pipeline.
9. **The index is defined for a mean-of-z-scores composite, and is undefined for
   at least one common alternative.** Under single-sample GSEA the reliability
   estimator that the index disattenuates with breaks down: Cronbach's α is
   defined here for a mean of k z-scored genes, an ssGSEA score is a rank-walk
   statistic with roughly 20-fold smaller dispersion, and α consequently collapses
   toward zero for random sets — 1,000 of 1,000 null draws were dropped for all 16
   signatures, so no corrected ssGSEA index exists. This is a limitation of the
   reliability estimator rather than a property of ssGSEA, but its practical
   consequence is real: the registered estimand is not scorer-portable as written.
   On the uncorrected contrast, which *is* comparable across scorers, the effect
   survives at ~56% of its magnitude (+0.180 vs +0.319, 15/16 vs 16/16 beating
   null), and at +0.174 under the pinned tie ordering of item 8, which
   leaves that reading unchanged; under a scorer-agnostic reconstruction it
   reverses sign. A
   scorer-portable reliability estimator is the natural next step and we have not
   built one.
10. **No external validation, because the intended external cohort does not
    carry the variable the estimand is defined on.** The index is defined under
    preserved-site cross-validation, so it requires a tissue source site. We
    attempted validation in CPTAC-3 and found no site variable in either of the
    two repositories that serve it: all 1,866 CPTAC-3 cases in the Genomic Data
    Commons return a missing tissue source site (the same query returns a named
    site for TCGA controls), and the Proteomic Data Commons exposes no
    institutional field at all — its site-like fields are anatomical
    (`biospecimen_anatomic_site`) or technical (`preservation_method`), and
    `tissue_collection_type` was null for every sample we inspected. The
    case-identifier prefix is not a usable substitute: 81% of CPTAC-3 Discovery
    cases fall into just two prefixes. Substituting some other grouping variable
    would redefine the estimand rather than validate it, so we report the
    attempt and its outcome instead. This is a constraint on any preserved-site
    analysis proposing CPTAC as a validation cohort, and we are not aware of it
    being noted elsewhere.

## Data Availability Statement

All primary data analyzed in this study are publicly available and were obtained
from existing repositories: Prov-GigaPath whole-slide-image embeddings for TCGA
(HuggingFace, CC-BY-4.0); Xena TOIL RSEM gene expression; PanCanAtlas ABSOLUTE
tumor purity calls; the TCGA Clinical Data Resource (TCGA-CDR) outcome table;
MSigDB Hallmark gene sets v2024.1 (CC-BY-4.0); and published UCSF genetic
ancestry calls for TCGA (12). No data were generated in this study. Derived
intermediates sufficient to reproduce every reported figure and table are
deposited with the analysis code below.

## Code Availability Statement

> ✅ **DECIDED 2026-09-04: Variant A, published as a CURATED public snapshot.**
> *Cancer Research* requires a Code Availability Statement separate from the
> Data Availability Statement, and a statement pointing at a private repository
> does not satisfy it.
>
> **What is being published is a curated snapshot, not this working
> repository.** A scan on 2026-09-04 found no credentials but did find, in
> tracked files, the author's cluster hostname and username, cluster and local
> home paths, a personal email address, and six internal working handoff
> documents. Flipping the working repository would publish all of those and
> leave them in its git history permanently. `pipeline/scripts/22_build_public_snapshot.py`
> stages the deliverables — `src/`, `scripts/`, `tests/`, `tools/`, the frozen
> result JSON/CSV, the manuscript, the pre-registration, the lockfile and
> `reproduce.sh` — and then **greps the staged tree for every one of those
> patterns and refuses to proceed if any survives**. Verified clean, and
> verified able to fail, on 2026-09-04: 141 files, 1,024,458 bytes.
>
> **PUSHED 2026-09-07.** The snapshot is live at
> `https://github.com/sathvikloke/immune-specificity-index`. It is a DIFFERENT
> repository from the private working one (`hne-tme-immune-specificity-audit`),
> which stays private and whose name must never be reused here — an earlier
> draft of this statement pointed at that private name. `CITATION.cff`'s
> `repository-code` is updated to match; the poster's QR (TODO-PRINT-3) is now
> unblocked and takes the same URL.
>
> Re-verified at push time: 169 staged files, 1,447,396 bytes, and the audit
> CLEAN — plus an independent `grep` for each forbidden pattern (cluster host,
> cluster username, cluster and local home paths, email) over the staged tree,
> returning nothing, because the tool passing its own check is not the same as
> the tree being clean. The patterns are not spelled out here: quoting them
> literally makes this file fail the very audit it is describing, which is why
> `22_build_public_snapshot.py` excludes itself. The rebuild caught a REAL leak
> that the 2026-09-04 build could not have seen: `science_gaps.json` had gained
> an absolute `cohort_path` under the author's home directory. It was fixed at
> source and the artefact regenerated, not hand-edited.

**Variant A — public GitHub repository. ← CHOSEN**

> All analysis code is publicly available at
> `https://github.com/sathvikloke/immune-specificity-index` under the
> MIT license. The repository contains the complete pipeline, the frozen
> configuration used for every reported result, an exact dependency lockfile,
> and a single-command reproduction script. The pre-registered protocol is
> timestamped at git tag `prereg-2026-08-18`, created before any result in this
> paper was computed. Reported values were produced on macOS 15 / arm64 with
> Python 3.13.9; see Limitations for the observed cross-platform spread.

**Variant B — archived Zenodo snapshot with a DOI.**

> All analysis code is archived at Zenodo under DOI
> [10.5281/zenodo.XXXXXXX — ⚠ to be minted before submission] and released under
> the MIT license. The deposit contains the complete pipeline, the frozen
> configuration used for every reported result, an exact dependency lockfile,
> and a single-command reproduction script, captured at the commit corresponding
> to this manuscript. The pre-registered protocol is included as git tag
> `prereg-2026-08-18`, created before any result in this paper was computed.
> Reported values were produced on macOS 15 / arm64 with Python 3.13.9; see
> Limitations for the observed cross-platform spread.

Variant B satisfies the requirement without making the repository public, and
Zenodo appears on AACR's own list of acceptable repositories. Variant A is
simpler and is what the 3.2% of accepted computational abstracts claiming code
availability generally do. Both name the platform, because Limitations item 8
reports that the estimate is platform-dependent in its fourth digit and a
reproduction claim that does not say where it was measured is not a claim.

## Authors' Contributions

**S. Loke:** Conceptualization, data curation, software, formal analysis,
validation, investigation, visualization, methodology, writing–original draft,
writing–review and editing.

*Sole-author submission. AACR requires a CRediT-style contributions statement for
every author; the roles above are the CRediT terms that apply.*

## Authors' Disclosures

No disclosures were reported.

*Wording taken from what* Cancer Research *actually prints; see
`submission/SUBMISSION-CHECKLIST.md` §8, which verified this phrasing across
published articles rather than from a third-party guide.*

## Acknowledgments

This study used only publicly available data; no funding was received for this
work, and no third party contributed to its design, analysis or interpretation.

**Use of generative AI.** Generative AI (Anthropic Claude) was used as a coding
assistant and as a drafting and editing aid throughout: writing and refactoring
the analysis pipeline, the verification tooling and the test suite, and drafting
and revising this manuscript. All analyses were specified, executed and checked
by the author, who takes full responsibility for the content, the correctness of
the reported numbers and the conclusions drawn. No AI system is listed as an
author, in accordance with AACR policy, which requires such use to be declared
both in the cover letter and in the Acknowledgments.

> **THIS SECTION EXISTED NOWHERE UNTIL 2026-09-04.** The cover letter stated that
> AI use "is declared in the Acknowledgments" while the manuscript had no
> Acknowledgments section, no Authors' Contributions and no Disclosures — an
> assertion about the manuscript that the manuscript did not support. The cover
> letter's own §7 had noted the gap without it being closed. Verify the CRediT
> role list and the disclosure wording against the *Cancer Research* author
> instructions before submission.

## References

**Bibliographic metadata for all twelve re-verified against PubMed on
2026-09-03** — authors, journal, year, volume, issue, page range and DOI
confirmed against the indexed record. This note covers the METADATA only; see
the cross-reference note below it. Two had to be
re-identified: the first search returned a paediatrics paper for Altman and a
yoga trial for Sonabend, so PMIDs are recorded here to prevent that recurring.

1. Schmauch B, Romagnoni A, Pronier E, et al. A deep learning model to predict
   RNA-Seq expression of tumours from whole slide images. *Nat Commun*
   2020;11(1):3877. PMID 32747659. [DOI](https://doi.org/10.1038/s41467-020-17678-4)
2. Xu H, Usuyama N, Bagga J, et al. A whole-slide foundation model for digital
   pathology from real-world data. *Nature* 2024;630(8015):181-188.
   PMID 38778098. [DOI](https://doi.org/10.1038/s41586-024-07441-w)
3. Chen RJ, Ding T, Lu MY, et al. Towards a general-purpose foundation model
   for computational pathology. *Nat Med* 2024;30(3):850-862. PMID 38504018.
   [DOI](https://doi.org/10.1038/s41591-024-02857-3)
4. Vorontsov E, Bozkurt A, Casson A, et al. A foundation model for
   clinical-grade computational pathology and rare cancers detection. *Nat Med*
   2024;30(10):2924-2935. PMID 39039250.
   [DOI](https://doi.org/10.1038/s41591-024-03141-0)
5. Venet D, Dumont JE, Detours V. Most random gene expression signatures are
   significantly associated with breast cancer outcome. *PLoS Comput Biol*
   2011;7(10):e1002240. PMID 22028643. [DOI](https://doi.org/10.1371/journal.pcbi.1002240)
6. Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer
   clinical data resource to drive high-quality survival outcome analytics.
   *Cell* 2018;173(2):400-416.e11. PMID 29625055.
   [DOI](https://doi.org/10.1016/j.cell.2018.02.052)
7. Howard FM, Dolezal J, Kochanny S, et al. The impact of site-specific digital
   histology signatures on deep learning model accuracy and bias. *Nat Commun*
   2021;12(1):4423. PMID 34285218. [DOI](https://doi.org/10.1038/s41467-021-24698-1)
8. Altman DG, Lausen B, Sauerbrei W, Schumacher M. Dangers of using "optimal"
   cutpoints in the evaluation of prognostic factors. *J Natl Cancer Inst*
   1994;86(11):829-35. PMID 8182763.
   [DOI](https://doi.org/10.1093/jnci/86.11.829)
9. Sonabend R, Bender A, Vollmer S. Avoiding C-hacking when evaluating survival
   distribution predictions with discrimination measures. *Bioinformatics*
   2022;38(17):4178-4184. PMID 35818973.
   [DOI](https://doi.org/10.1093/bioinformatics/btac451)
10. Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The
    Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell
    Syst* 2015;1(6):417-425. PMID 26771021.
    [DOI](https://doi.org/10.1016/j.cels.2015.12.004)
11. Li J, Ji L. Adjusting multiple testing in multilocus analyses using the
    eigenvalues of a correlation matrix. *Heredity (Edinb)* 2005;95(3):221-7.
    PMID 16077740. [DOI](https://doi.org/10.1038/sj.hdy.6800717)
12. Carrot-Zhang J, Chambwe N, Damrauer JS, et al. Comprehensive analysis of
    genetic ancestry and its molecular correlates in cancer. *Cancer Cell*
    2020;37(5):639-654.e6. PMID 32396860.
    [DOI](https://doi.org/10.1016/j.ccell.2020.04.012)

**Cross-reference check, 2026-09-04.** **Every entry above is cited in the
body.** The previous reference 9 (Rooney et al., immune cytolytic activity) was
listed but cited nowhere, and the manuscript contains no cytolytic-activity
discussion for it to attach to -- a leftover from an earlier draft. It was
**dropped by authorial decision on 2026-09-04** rather than given an invented
citation site. The list is
now **twelve references, all cited, all verified against PubMed**.

**Citation style converted 2026-09-04.** In-text citations were author--year
("Carrot-Zhang et al., *Cancer Cell* 2020") and are now **numbered in
parentheses**, with the reference list renumbered into **citation-sequence
order** (order of first appearance in the body), which is what a numbered style
requires. The style was **verified empirically rather than assumed**: across 24
*Cancer Research* articles in PubMed Central (2023-2025) whose full text encodes
its in-text citations, **18 render them parenthetically and 4 as superscripts**,
so parenthetical is the house style. `aacrjournals.org` could not be read to
confirm this from the author instructions directly -- it is behind a bot
challenge (HTTP 403) and blocked by browsing policy -- so this survey is the
evidence, and the instructions themselves remain an open item.

*Bibliographic metadata retrieved from PubMed.*
