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
([`05-PRE-REGISTRATION.md`](05-PRE-REGISTRATION.md); recorded by the annotated tag
`prereg-2026-08-18` in the **private** working repository, which the published
snapshot does not carry — see the Code Availability Statement).
All numbers trace to `pipeline/results/nsclc_v3` and `pipeline/results/pancancer_v3`.
`pancancer_v2` is **identical on every immune-specificity number** — mean excess
0.291144, range 0.165–0.390, checked 2026-09-02 — and differs only in **lacking**
`null_draws.npz`, the persisted per-draw null correlations that make the
family-level rotation null and Figure 1B's distribution possible.
(This sentence said v2 *carried* the file until 2026-09-07. It is the other way
round, verified by `ls`: `pancancer_v3/null_draws.npz` exists, `pancancer_v2/`
has no `.npz` at all. The direction is the whole content of the sentence, since
v3 is the reported run precisely because it has the draws.) An earlier
version of this line said v2; it was already reporting v3 results, since it quotes
the rotation null.
Figures are rendered by `pipeline/scripts/10_make_figures.py` into
`pipeline/results/figures/` as PNG at 300 dpi and vector PDF, drawn at their
printed size (no text under 8 pt at a 6.75 in double column), with `--tiff`
(600 dpi) for journal revision; the A0 poster's figures are separate renders
(`--poster`) drawn for their boxes.
**References verified against PubMed 2026-08-19** with PMIDs recorded.

---

# Proliferation, not immunity: H&E-inferred tumor microenvironment signatures are immune-specific but carry no prognostic advantage over random gene sets

## Abstract

*Journal format.* ***Cancer Research* abstracts are unstructured** — verified
across 40 published articles, which parse as an unlabeled body plus a
**Significance** statement, with no Background/Methods/Results/Conclusions
headings. Observed lengths: 131–264 words including Significance, with
Significance itself 18–40 words. The version below is written to that format and
measures **236 words of body and 28 of Significance, 264 in total** — at the top
of the observed range (262 until 2026-09-23, when the reliability clause was scoped
to the median curated signature), though the *official* cap is unverified because
`aacrjournals.org` refuses automated requests. All three counts are re-derived
by the checker on every run rather than maintained by hand (session 53), and the
measure is whitespace-delimited tokens containing a word character: it excludes
three spaced em-dashes and the `=` of a reported statistic, which a naive split
would count as four more words and put the total just past the observed ceiling.
The claim is therefore true under a word count and not under a token count, and
the margin against the ceiling is two words. The structured 456-word version it
replaces is folded below and is the right shape for the preprint, which has no
such convention.

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
sets' Cronbach's α from 0.97 to 0.80 while the median curated signature held above 0.97,
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
<summary>Superseded structured abstract (456 words) — kept for the preprint</summary>

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
retained, and all of them are diagnostic (DX) formalin-fixed slides. Slides were
averaged to one vector per patient before any splitting.

Expression is Xena TOIL RSEM on TOIL's log2(TPM + 0.001) scale, mapped from Ensembl to HUGO symbols via
HGNC (42,356 pairs), giving 41,046 genes. Expression samples were not
restricted by sample type. Pan-cancer, a patient's TOIL samples are averaged; in
NSCLC, the first sample in TOIL's order is used. So a minority of profiles
include or consist of non-tumor tissue; the effect is measured in Results
(tumor-only expression). Tumor purity is PanCanAtlas ABSOLUTE;
CPE and ESTIMATE were deliberately avoided because they are expression-derived
and would make the purity control circular. Outcome is progression-free interval
from TCGA-CDR, which Liu and colleagues (6) recommend for all but 4 of 33 TCGA
types.

The protocol registered the NSCLC analysis and a CPTAC validation (item 10);
the pan-cancer cohort was added after filing, on 2026-08-19, as a replication
under the same registered estimand and fixed settings, and the protocol's
deviation appendix records that too.
The **pan-cancer cohort** is 7,168 patients across 31 cancer types and 619 tissue
source sites (2,442 PFI events). Of these, 19 patients (16 from a site whose
other 64 patients all have testicular germ cell tumors) carry no clinical record
in the embeddings release and therefore no cancer type. They enter the
unresidualized analyses, but their within-type residualized scores are
undefined, so the index and its bootstrap use the 7,149 complete cases. The
bootstrap's complete-case mask was added the day after the protocol was filed,
when the pan-cancer interval first came back undefined; the protocol's deviation
appendix records it, and it does not touch NSCLC, where no patient is missing.
The **NSCLC cohort** is 944 patients (LUAD 471,
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

**Two split schemes are run, and the comparison between them is itself a
reported result.** Alongside the preserved-site partition, patients are assigned
to folds at random with no regard to site, using the same fold count and the
same seed; every model is fitted under both. Two quantities are reported from
the pair. Δ*r* is the per-signature difference in out-of-fold correlation
between the schemes. Δ-MAE is a **paired difference in mean absolute error**:
out-of-fold predictions from the same ridge head are matched on (patient,
signature) across the two schemes, and the difference is taken as preserved-site
minus random-patient, so a **positive value means preserved-site splitting costs
accuracy**. Its interval is a patient-clustered bootstrap on the paired
differences — the same resampling used for the index, not a normal
approximation.

**Covariate baselines are fitted under the random-patient split, and that scope
limits what they can be read as.** A covariate-only design matrix of one-hot
tissue source site, one-hot cancer type, numeric ABSOLUTE purity and one-hot
stage is fitted with the same ridge head and the same folds as the embedding,
and the embedding's advantage is its out-of-fold correlation minus that
baseline's, per signature, on that split. Simpler baselines — cohort mean,
per-type mean, site only, type only, purity only — are computed the same way.
Because the split is random over patients, a held-out patient's site is present
in training, so these baselines measure how much of a signature a covariate set
*explains within this cohort*; they are not estimates of how a covariate model
would transfer to an unseen site, which is the question the preserved-site
scheme asks.

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
alongside them in Limitations. A third ordering sat inside a library:
scikit-learn's grouped k-fold splitter orders equally sized patient groups
(after collapsing to patients, all of them) by the same kind of non-stable
sort, so the random-patient partition behind the secondary quantities
depended on the platform (Limitation 8), as did the patient-disjoint folds of
the site-classification control; the index uses neither. That ordering is now
pinned as well: the pipeline assigns those folds with the splitter's own
procedure under a stable sort, which changes nothing where no group sizes tie.
The frozen runs' secondary quantities are reported as run, and the pinned
re-run's beside them in Results.

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
There is no dichotomization path in the code, and the reduction from model to
scalar was fixed in advance, both to avoid documented inflation of type I error
(8) and C-hacking (9).

**Concordance is computed within cancer type.** Progression-free interval differs
greatly between tumor types, and cancer type alone reaches C = 0.676 pan-cancer
(Harrell's C with each type's PFI event proportion as the score, in sample;
0.678 and 0.681 under two other type-only scores);
an unstratified concordance therefore largely measures which cancer a sample is,
and any tissue-correlated score inherits that. Unstratified, 94% of signatures
"beat" their null; stratified, 31% do. Only the stratified figures are reported.
Stratification is not in the filed protocol, which fixed only Harrell's C: it was
introduced on 2026-08-19, the day after filing and after the unstratified
pan-cancer outcome had been seen, and the protocol's deviation appendix records
it. It leaves NSCLC's count unchanged.

Each concordance is folded to 0.5 + |C − 0.5|, so an anti-prognostic score is
not rewarded. It is then mapped onto Somers' D (D = 2C − 1) and compared with
its random sets in Fisher z, as the index is. The observed value's standard
error is taken as 1/√(n − 3), the formula for a correlation, with n the patients
who have a usable endpoint. The null mean's sampling error is added in
quadrature. For a concordance statistic that standard error is an approximation
(Limitation 6), and the intervals and minimum detectable effects below inherit
it.

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
IL6/JAK/STAT3 signaling is **excluded** from the immune group because MSigDB's
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

203 automated tests; the analysis is deterministic given seed 0 **on a fixed
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

Run end to end on Linux/x86_64 from the raw files — download and hash-verify,
rebuild the derived inputs, both cohorts, the ancestry supplement, the figures —
the pipeline reproduces the corrected-ordering values for both cohorts bit for
bit, in 2 h 26 min on eight CPUs. That covers the index and the preserved-site
partition it is computed on. The frozen runs' random-patient secondary
quantities are platform-specific, as Limitation 8 sets out; with that split
pinned (above), a re-run gives the same secondary quantities on either
platform.

**Software.** All analyses were performed in Python 3.13.9 (RRID:SCR_008394) on
macOS 26.5.2 (Darwin 25.5.0; arm64). Ridge regression with generalized cross-validation used
scikit-learn 1.6.1 (RRID:SCR_002577). Array and dataframe handling used NumPy
2.1.3 (RRID:SCR_008633) and pandas 2.2.3 (RRID:SCR_018214), with pyarrow 24.0.0
as the parquet engine. Rank and product-moment correlations, the normal
quantiles behind the Fisher *z* transform, permutation percentiles and the χ²
tests used SciPy 1.15.3 (RRID:SCR_008058). Ordinary least squares for the
covariate baselines and the Cox proportional-hazards fits of the outcome arm
used statsmodels 0.14.4 (RRID:SCR_016074); the Benjamini–Hochberg correction is
implemented in this project's own code rather than called from a library.
Figures were rendered with matplotlib 3.10.0 (RRID:SCR_008624), and the
automated tests run under pytest 8.3.4. Gene sets are the Hallmark collection of
the Molecular Signatures Database v2024.1.Hs (RRID:SCR_016863), pinned by
version in the fetch URL rather than resolved as "latest". The linear-algebra
backend loaded at run time is OpenBLAS 0.3.29 (numpy's build record names
0.3.21), reported because the cross-platform difference
described above turns on it. **scikit-learn's version is load-bearing rather
than incidental:** RidgeCV's alpha-selection defaults changed between minor
releases and the index depends on the selected penalty, so a lower-bound
specification is not a reproduction guarantee. Exact versions for every package
are pinned in `requirements-lock.txt`; `scripts/20_check_versions.py` checks the
four whose arithmetic reaches the estimand — NumPy, pandas, scikit-learn and
SciPy — against that lockfile and exits non-zero on drift.

**Where the remaining analyses are specified.** Twelve secondary analyses are
defined at first use in Results rather than here, because each is short and
reads better beside the number it produces:
- the family-level permutation test;
- the label-side control's regression and its within-type permutation;
- the site-classification AUROC;
- the ComBat-adjusted arm;
- the purity partial correlation;
- incremental adjusted R²;
- Cramér's V for ancestry against site;
- the Benjamini–Hochberg correction;
- the repeated site-to-fold partitions;
- the outcome arm's power curve;
- the enumeration over signature categories;
- the one-at-a-time sensitivity analyses.

They are listed here so a reviewer looking for them in Methods is not left
concluding they are missing. A plate-within-site control named in earlier
drafts is withdrawn, because the inputs carry no plate identifier (see the
label-side control in Results).

## Results

### The image predicts TME signatures, and preserved-site splitting costs little

Under preserved-site cross-validation the
median image–signature correlation was 0.630 pan-cancer and 0.400 in NSCLC. The difference from a site-agnostic
(random-patient, in the terms of Methods) patient-level split was small (median Δr = +0.020 and +0.023; paired calibration
loss Δ-MAE = 0.0088 [0.0078, 0.0098] and 0.0054 [0.0024, 0.0085]), indicating
that the models are not merely site detectors even though site is highly
recoverable (below).

**A necessary caveat.** Pan-cancer, a covariate-only baseline of site, cancer
type, purity and stage reached median r = 0.655, and the embedding added
**−0.000** over it. In NSCLC, where cancer type carries far less information, the
embedding added +0.031. The pan-cancer raw correlation is therefore substantially
tissue of origin, and must not be read as image performance. This is precisely
why the ISI is defined on within-type residualized scores.

**Pinned random-patient split.** The figures above come from the frozen runs,
whose random-patient folds depended on the platform (Methods). With that split
pinned, which does not touch the index, NSCLC gives a median Δr of +0.013,
Δ-MAE 0.0043 [0.0014, 0.0070], an embedding advantage over covariates of
+0.026 and a median site AUROC of 0.991; pan-cancer gives a median Δr of
+0.023, Δ-MAE 0.0093 [0.0083, 0.0102], an embedding advantage of −0.001 and a
median site AUROC of 0.999. No conclusion changes.

### The immune-specific excess is positive and replicates

The ISI was **0.291 (95% CI 0.269–0.313)** pan-cancer and **0.318 (0.261–0.376)**
in NSCLC. All 16 signatures exceeded their own null in both cohorts (per-signature
excess 0.165–0.390 and 0.148–0.440; Fig. 1A; per signature in Supplementary
Table S1). Median residualized correlations were 0.394
and 0.372 against null means of 0.106 and 0.051 (pan-cancer, against the full
null distribution, in Fig. 1B).

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
their α fell from 0.984 to 0.976 (Fig. 2B). Per-signature α, residualized and
raw, with the curated-minus-random gap, is Supplementary Table S2.

The consequence is that the reliability gap computed the conventional way, on raw
scores, is badly understated: +0.014 pan-cancer and +0.032 NSCLC, against
**+0.184 and +0.167** when computed on the residualized score that the estimand
actually disattenuates — a 13-fold and 5.3-fold difference.

The gap scales inversely with panel size (Fig. 2A; Spearman ρ = −0.50,
p = 0.048 pan-cancer; ρ = −0.59, p = 0.015 in NSCLC): it averages +0.141 across the
eleven 200-gene Hallmark sets (198–200 genes after filtering) and
reaches +0.420 for the smallest, 36-gene angiogenesis set.

Because real clinical panels are smaller than any Hallmark set, we tested the
scaling *inside* the measured range rather than extrapolating below it, by
subsampling each signature to fixed k and comparing against random sets of the
same size. The gap rises monotonically: 0.084 at k=160, 0.141 at 80, 0.214 at 40,
0.290 at 20 and **0.453 at k=10** (Spearman ρ = −1.000, p < 0.0001) — 5.4-fold
across the range; the full bracket is Supplementary Table S3. These are raw-score alphas and therefore a conservative lower
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

**Two further checks (Supplementary Note 2).** A third scorer, PLAGE,
attenuates the uncorrected contrast as ssGSEA does (+0.1353 [0.0889, 0.1834], 13
of 16 signatures above their null) but stays positive under the scorer-agnostic
reconstruction (+0.0772), so the reversal is specific to ssGSEA. Replacing
Cronbach's α with McDonald's ω lowers the reconstructed NSCLC index by about a
quarter, to 0.2147 [0.1739, 0.2554], with 15 of 16 signatures still above their
null; read it as a bound, not a better estimate.

### Prognostically, the signal is proliferation and stroma — not immunity

Pan-cancer, with concordance computed within cancer type, five signatures beat
their random-set null (Fig. 3): **G2M checkpoint** (excess 0.069), **E2F targets**
(0.062), **angiogenesis** (0.051), **epithelial–mesenchymal transition** (0.038)
and **hypoxia** (0.030). No interferon, inflammatory, complement or
allograft-rejection signature did. Using MSigDB's own process categories
(10), 0 of 6 immune signatures beat their null
against 5 of 10 others (Mann–Whitney p = 0.031; Fisher p = 0.093). All five that
do fall within the six signatures MSigDB files under proliferation, development
or pathway. Five winners drawn at random from the 16 would all land there with
probability 0.0014, and no relabeling of the 16 gives six signatures a larger
mean excess over the other ten than these six have (1 of 8,008). That grouping
was named after the result, and signatures that share genes are not
exchangeable, so both figures size the pattern rather than test a hypothesis.

In NSCLC no signature beat its null. **That is an underpowered null, not a
negative result**: the NSCLC minimum detectable effect is 0.046 in C-index (Harrell's C), while
the real pan-cancer advantage is approximately 0.03 (pan-cancer MDE 0.017). The
two cohorts are consistent; the smaller one simply cannot see an effect of this
size: at a C-index advantage of 0.03 the NSCLC arm's power is 45% (23% at 0.02,
68% at 0.04), by the same normal approximation that gives the MDE. We report
the NSCLC result only with its MDE attached.

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

### Multiplicity, panel redundancy and the label-side control (Control C)

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
The regression uses the 202 sites that contribute at least 10 patients (5,775
patients, after the 17 among them without a recorded cancer type drop out).
That figure is almost entirely an artifact of two things: source site largely
determines which diseases a center contributes, and dummy-coding 202 sites
inflates R² by construction. Conditioning on cancer type reduces it to 5.5%.
Calibrating against site labels permuted *within* cancer type, over every site
with 200 permutations per signature, leaves a median of 4.7% across the 16
signatures (range 2.8% to 6.0%; 6 of 16 at p ≤ 0.05). The first calibration
we ran, on allograft rejection alone, gave 4.1% (p = 0.010). So the
label-side site artifact is small once disease is accounted for (per signature
in Supplementary Table S4), and the site information the embeddings carry is
not mirrored by a comparable artifact in the targets. Whether the remainder is
technical cannot be told from these data. The plate identifier that would
isolate processing batch is absent from the inputs: neither the slide barcodes
nor the expression sample identifiers carry it. The plate-within-site column of
Supplementary Tables S4 and S5 is therefore zero by construction and is
evidence of nothing.

On the NSCLC scope the same control attributes 14.5% of the label to site at
the median, 9.9% once cancer type is conditioned on, and 7.0% under the
within-type permutation calibration (per signature in Supplementary Table S5;
Supplementary Note 2), so the two cohorts are tabulated separately
rather than stacked.

### Sensitivity to the choice of site-to-fold partition

Preserved-site cross-validation admits many valid assignments of sites to folds,
and the analysis uses one. We re-ran the primary endpoint under six assignments
in each cohort, varying only the partition with the null gene sets and bootstrap
draws held fixed, and later under 24 assignments in each cohort (per partition,
both cohorts and both sweeps, in Supplementary Table S6; the full analysis is
Supplementary Note 2). The 24-partition runs use the pinned tie ordering of
Limitation 8, so they sit beside the six-partition figures rather than replacing
them.

In NSCLC the between-partition sd is 0.0089 over six partitions and 0.0077 over
24, against a within-partition bootstrap standard error of 0.0294. Combining the
two as independent sources widens the interval from [0.2542, 0.3695] to
[0.2516, 0.3720]: the reported interval is 4.5% too narrow (3.3% at 24), and
partition choice accounts for 8.4% of total variance (6.4% at 24).

Pan-cancer the between-partition sd is 0.0065 over the five of six assignments
that could be evaluated and 0.0072 over 24, against a bootstrap standard error of
0.0111. The interval widens from [0.2717, 0.3152] to [0.2683, 0.3187]: the
reported interval is 16.0% too narrow (19.4% at 24), and
partition choice accounts for 25.7% of total variance (29.9% at 24). The
conclusion is unchanged in both cohorts — the honest interval still excludes
zero by a wide margin.

With 24 draws per cohort, the bootstrap component falls with cohort size as
sampling theory says it should, by 2.652 against the 2.756 predicted by the
square roots of the sample sizes, but the partition component does not fall: its
ratio between the cohorts is 1.060, with an F-based 95% interval of
[0.697, 1.611] that contains 1. Over these two cohorts the spread across
partitions behaves as a floor set by the partition design rather than by cohort
size, which is why it is the larger share of uncertainty in the larger cohort.
Two cohorts cannot show how the floor depends on the design, and we do not claim
that they do. For the split-robustness contrast Δr, partition choice matters
considerably more.

### Sensitivity to other analysis choices

None of the analyses in this section was pre-registered; each changes one choice
and recomputes only the index, against the corrected-ordering NSCLC value of
0.2922 [0.2355, 0.3527] (Supplementary Note 2). Removing a single
expression axis computed across the two NSCLC types, rather than within each,
collapses the NSCLC index to −0.0096 [−0.0630, 0.0466]; pan-cancer the same change
raises it to 0.3512 [0.2692, 0.4353]. The index is therefore a statement about
expression variation within a cancer type, which is how it was registered. With
3 folds the NSCLC index is 0.2188 and with 10 it is 0.3154, so its value is
specific to the registered five folds, though its sign and the excess of curated
over random sets hold at every number of folds tried. Restricting expression to
the primary-tumor sample raises the NSCLC index to 0.3448 [0.2731, 0.4079] and
moves pan-cancer to 0.2937 [0.2701, 0.3149]. Matching the null on size alone,
rather than on size and mean expression, gives 0.3041, so expression matching
lowers the NSCLC index by 0.012; the null subsample inside the bootstrap leaves
the estimate unchanged. All values elsewhere in this paper use the original
inputs.

### Negative controls

Tissue source site was recoverable from the embeddings at median one-vs-rest
AUROC **0.998** pan-cancer and 0.992 in NSCLC, with 100% of evaluable sites above
0.9 in both (Fig. 4A; per site, both cohorts, plus the ComBat arm, in
Supplementary Table S7). We report this raw figure as the confounding measure and deliberately
do **not** report a post-ComBat AUROC: per-site centring imposes a within-site
zero-sum constraint that makes the site classifier systematically anti-predictive
(measured 0.005), so the corrected number is uninterpretable in either direction.
The deployment-relevant statement is that under preserved-site folds ComBat is
estimable for **0%** of held-out samples by construction — a model meeting a new
hospital is in exactly that position.

The control does not depend on cohort size in any range we can evaluate:
pan-cancer the median stays at 0.998 at half the cohort and 0.994 at a quarter,
and sites contributing as few as five patients are still identifiable
(Supplementary Note 2).

> **Figure 4.** (A) Distribution of one-vs-rest site AUROC across evaluable
> sites in both cohorts. (B) Median correlation achieved by a covariate-only
> baseline against the image embedding — pan-cancer the two are
> indistinguishable.

The image prediction's partial correlation with the signature, controlling for
ABSOLUTE purity, was 0.547 pan-cancer (against 0.627 unadjusted), so purity
explains part but not most of the association. Incremental adjusted R² of the
image over purity, type, site and stage was 0.036 and 0.054. The covariate-only
baseline of the first Results section is drawn beside the image embedding in
Fig. 4B.

### Ancestry (supplementary)

A pan-cancer ancestry arm (6,580 patients) is reported as Supplementary Table S8
with a negative headline — **the contrast is not identifiable from tissue source
site**: only 0.6–1.0% of sites carry at least 10 patients of both the reference
and any comparison group (Supplementary Note 2). Ancestry calls are the
published UCSF consensus set (12).

## Discussion

The headline is two-sided, and both sides matter.

**H&E carries real immune-specific information.** After conditioning the
signature on the dominant expression axis and correcting for the measurement
reliability of curated modules, curated TME signatures remain substantially more
image-predictable than size- and expression-matched random gene sets, in 31
cancer types. This is a stronger result for the modality than the composition
hypothesis predicts, and it is obtained under site-disjoint (preserved-site)
validation. The
claim is made for the mean-of-z-scores composite the index is defined on, and
does not carry over to an arbitrary scoring rule: under single-sample GSEA the
uncorrected contrast survives at about half its magnitude, but the corrected
estimand does not exist at all, for the reason given below.

**That information is not prognostic.** The prognostic content of the Hallmark
TME family, measured against the only fair comparator, is proliferative and
stromal. Immune programs — interferon α and γ, inflammatory response,
complement, allograft rejection — do not beat random gene sets on
progression-free interval within cancer type. Since these are exactly the
programs that image-based TME inference is proposed to deliver, the practical
implication is direct: **image-inferred immune scores should not be used as
prognostic biomarkers**, not because the image reads them poorly, but because the
underlying transcriptomic scores are not prognostic once compared fairly.

**A methodological point that generalizes.** Any comparison of a curated gene set
against a random one is confounded by reliability, and the size of that
confounding depends on how much of the sets' internal consistency comes from a
shared global axis. Because that dependence is inverse in panel size, the
smallest panels — which are the ones nearest clinical deployment — are the most
exposed. Reporting a raw-score α understates the problem by up to an order of
magnitude.

The correction that repairs this confounding is itself not scorer-free, and that
is the sharper half of the point. Disattenuation divides by the square root of
the comparator's reliability, so a scoring rule that drives the reliability of
*random* sets toward zero — as a rank-walk statistic does, from 0.799 to 0.258
here — turns a modest correction into a numerically violent one and can remove
the estimand entirely. Any index built on a curated-versus-random comparison
therefore has to name its scorer as part of its definition rather than treat
scoring as an implementation detail, and has to report what its comparator's
reliability actually is. Ours does not survive that substitution, and we would
expect most such indices not to.

<!--
Deliberate omission, decided 2026-09-07: the partition-variance result is NOT
mentioned in the Discussion. Per the A5 addendum in 14-SCIENCE-AUDIT.md the
cross-cohort contrast is a contrast of point estimates whose chi-square
intervals overlap, so it is not a contribution and does not belong in the
paper's argument. It stays a Results paragraph plus Limitation 8. (Since
2026-09-24, B-21, the Results carry a summary and the full analysis is
Supplementary Note 2.) Do not
re-find this as a gap -- it is a decision, not an oversight. The scorer /
ssGSEA half of the same gap WAS acted on, in the two paragraphs above.
-->

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
   Supplementary Table S9). The pan-cancer ISI is meaningful because it is
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
8. **The ISI is platform-dependent in its second decimal, and the cause is a
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
   assignment does not (per-machine hashes in Supplementary Table S10).
   Substituting a stable sort
   makes the partition hash identical on both platforms and, on a synthetic
   cohort driven through the same estimator chain, reproduces the Linux value on
   macOS to 15 significant figures, a residual of 7 × 10⁻¹⁶. Genuine
   cross-platform floating-point noise in this pipeline is therefore six orders
   of magnitude smaller than the discrepancy long attributed to it.
   The discrepancy is thus an instance of the partition-choice sensitivity
   already quantified above rather than a separate defect: 0.0218 is 2.5 times
   the six-partition standard deviation of 0.0089 (2.9 times the 24-partition
   value). It does not alter the sign, the
   significance, or the finding that 16/16 signatures exceed their nulls, so no
   conclusion here turns on it; it does mean the estimate should not be read
   more finely than its second decimal. **The stable sort is now the shipped default**, and a regression
   check fails if it is reverted; the frozen results retain the original sort so
   that they remain the runs that were reported (and, for NSCLC, the one the
   protocol recorded at filing).

   The full forensic account of the ordering defect — the two corrected re-runs,
   the cross-platform validation to sixteen significant figures, the scope of the
   library ordering, and the pipeline-wide audit that found one further instance —
   is Supplementary Note 1; no claim here depends on reading it.
9. **The index is defined for a mean-of-z-scores composite, and is undefined for
   at least one common alternative.** Under single-sample GSEA, Cronbach's α —
   defined here for a mean of k z-scored genes — collapses toward zero for random
   sets: 1,000 of 1,000 null draws were dropped for all 16 signatures, so no
   corrected ssGSEA index exists. This is a limitation of the reliability
   estimator rather than of ssGSEA (PLAGE falls outside the definition for the
   same reason), but its consequence is real: the registered estimand is not
   scorer-portable as written. The uncorrected contrast survives at ~56% of its
   magnitude (+0.180 vs +0.319, 15/16 vs 16/16 beating null,
   and at +0.174 under the pinned tie ordering of item 8), a scorer-agnostic
   reconstruction reverses its sign, and the choice of reliability estimator
   moves the reconstructed index by about a quarter,
   with the same 15 or 16 of 16 signatures above their null (Results, "The
   index depends on how the signature is scored"). A scorer-portable
   reliability estimator is the natural next step and we have not built one.
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
11. **The combined partition-and-bootstrap interval treats the two sources of
    variance as independent, and they are not.** The partition analysis adds the
    between-partition variance to the bootstrap variance. Cross-validation
    variance decomposes into components that are correlated through the overlap
    between training and test sets, and no estimator of that variance is
    unbiased for every distribution (13). The combined interval is therefore an
    approximation whose direction of error we have not established. More
    partitions sharpen the estimate of the between-partition standard deviation;
    they do not remove this assumption.
12. **The index's value belongs to its registered choices.** Three choices were
    fixed in advance: removing the global axis within cancer type, five folds,
    and each cohort's expression samples as the frozen inputs hold them.
    Varying each one moves the index, and the axis and three folds move it by
    more than the half-width of its interval (Results, "Sensitivity to other
    analysis choices").
    - Defining the axis across NSCLC's two types removes the NSCLC index
      entirely; pan-cancer, the same change raises it.
    - The index rises with the number of folds.
    - Restricting NSCLC expression to the tumor sample raises it.

    Its sign and the curated-over-random excess held under every variant except
    the NSCLC global axis. The reported values are therefore the registered
    estimand, not a quantity expected to transfer unchanged across these
    choices.

## Data Availability Statement

All primary data analyzed in this study are publicly available and were obtained
from existing repositories: Prov-GigaPath (2) whole-slide-image embeddings for
TCGA (HuggingFace dataset `seandavis/tcga_provgigapath_embeddings`, revision
`073115403c2fc5134ee8d1332c603edba591dddb`, licensed CC-BY-4.0); Xena TOIL RSEM
gene expression; PanCanAtlas ABSOLUTE
tumor purity calls; the TCGA Clinical Data Resource (TCGA-CDR) outcome table;
MSigDB Hallmark gene sets v2024.1 (CC-BY-4.0); and published UCSF genetic
ancestry calls for TCGA (12). No new primary data were generated in this
study; everything deposited is derived from these sources. The raw inputs are
public but large — the embeddings parquet is 467 MB and the expression download
741 MB, and the mapped expression matrix built from it is 792 MB on disk,
expanding to 3.2 GiB as a dense in-memory copy — and are not redistributed with
the code. The six
derived input files the cohort analyses read are rebuilt from those raw
inputs by `scripts/00_build_interim.py`, which is deposited;
what is deposited alongside it are the derived intermediates the reported
figures and supplementary tables are built from, namely the frozen per-cohort
result tables and summaries, the per-signature null draws (`null_draws.npz`,
without which figure 1B's null distributions cannot be drawn), and the
per-signature gene counts surviving expression filtering (`gene_set_sizes.csv`,
figure 2's x-axis). Two tests hold this to account rather than leaving it as
prose: one renders all five figures inside a staged deposit with no input data
reachable, and one builds all ten supplementary tables there. Naming the
intermediates is deliberate: a reader can confirm in seconds that the deposit
holds what this sentence promises, which a claim of sufficiency alone does not
permit.

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
> Re-verified at the FIRST push, 2026-09-04: 169 staged files, 1,447,396 bytes,
> and the audit CLEAN. ⚠ **Those two figures described what was LIVE then, and
> both are now superseded — see the CLOSED note below for the current pair.** A
> rebuild on 2026-09-07 produced **172 files, 1,756,296 bytes**, because a root
> `README.md` was added — the landing page for this URL was a bare list of five
> files until then. That README was staged and not pushed for three sessions,
> which is the defect the CLOSED note records. Re-measure both figures at the
> next push rather than carrying either pair forward. (These
> rebuild figures are themselves already superseded — see below. Every edit to a
> staged file moves the byte count, so the pair is only true as of the rebuild
> that produced it, and it is quoted with its date for that reason.)
> The rebuild figures moved again later on 2026-09-07, 170 → 172 files and
> 1,493,616 → 1,756,296 bytes, for a reason that belongs in the record: the two
> `null_draws.npz` files (249,170 bytes) are now staged, because **without them
> the Data Availability Statement above was false.** Figure 1B's null cannot be
> drawn from anything else in the deposit, and `10_make_figures.py` silently
> substituted a mean ± SD band while printing that it had drawn the real
> distribution. The other 13,510 bytes are this session's edits to files that
> are themselves staged.
>
> ⚠ **The statement was false a SECOND way, and the first fix did not reveal
> it — running the figure script against a staged tree did.** With
> `null_draws.npz` in place, `10_make_figures.py` still **exited 1** inside the
> snapshot, at figure 2: its x-axis needs the count of genes per signature
> surviving expression filtering, whose only route was
> `data/interim/expr_nsclc.parquet` and the Hallmark GMT — and `pipeline/data/`
> is excluded from the snapshot deliberately, because the inputs are public but
> large. Every reader reproducing from the deposit hit it. The quantity is
> **16 integers**, so it is now deposited as `results/gene_set_sizes.csv`, which
> the existing `.csv` rule stages; recomputation from `data/` stays authoritative
> where those inputs exist and the two are **compared**, so a stale deposit fails
> loudly. Measured after the fix: the snapshot renders **all five figures,
> exit 0**, with figure 2 reproducing ρ = −0.500, p = 0.0485 identically and
> figure 1B drawing violins. All nine supplementary tables already built from the
> snapshot alone. **This is now a test, not a claim** — the suite stages a real
> snapshot into a temp directory, renders from it with no `data/` reachable, and
> a companion test deletes the deposit to prove the check can fail.
> After these changes a rebuild staged 173 files on 2026-09-07; as of
> 2026-09-24 a rebuild stages **322 files**, and the deposit was republished at
> that count the same day (commit `80b5535`); it had held the 173 since
> 2026-09-08. The 149 added in between are sessions 47 to 50's sensitivity runs
> (the corrected-ordering pan-cancer run, both 24-partition sweeps, the
> tumor-only runs, the global-axis runs, the seven one-at-a-time NSCLC variants,
> the third scorer, the reliability comparison), their derived records
> (partition scaling, Control C calibration, null reliability, subsample AUROC,
> outcome power and composition, the per-machine ordering records, the sort
> fix's fold-change census) and the
> scripts that produce them, plus the committed Ensembl-to-HUGO map and, since
> 2026-09-24, `DATA-LICENSES.md`, the data note split out of `LICENSE`,
> Supplementary Notes 1 and 2, and the two pinned-split runs' files.
> `pipeline/results/README.md` classifies every one of them, so this paragraph
> states the count and the kinds rather than a list that goes stale each time a
> run lands. The staged tree measured
> **1,780,474 bytes as of the rebuild at 19:46 CDT on 2026-09-07**, and
> **1,806,703 bytes as of the rebuild at 21:14 CDT on 2026-09-07** — the count
> unmoved across both, the byte total moved by edits to staged files, which is
> precisely why the two are treated differently — and the difference is
> deliberate. The **file count is a
> checked authority** in `19_check_numbers.py`, derived from a live rebuild
> rather than stored, because it moves only when a file is added or removed —
> which is the drift that went unnoticed four times (169 → 170 → 172 → 173). The
> **byte total is recorded as a dated observation and is not gated**, because
> `pipeline/scripts/` and this manuscript are both staged, so it moves on every
> character of every edit; a gate on it would be red almost always and would be
> ignored within a session. A stale dated observation reads as history; a stale
> bare figure reads as a fact. Every historical pair above is likewise phrased
> so the checker does not read it as a current claim.
> The audit was CLEAN on the rebuild too — plus an independent `grep` for each
> forbidden pattern (cluster host,
> cluster username, cluster and local home paths, email) over the staged tree,
> returning nothing, because the tool passing its own check is not the same as
> the tree being clean. The patterns are not spelled out here: quoting them
> literally makes this file fail the very audit it is describing, which is why
> `22_build_public_snapshot.py` excludes itself. The rebuild caught a REAL leak
> that the 2026-09-04 build could not have seen: `science_gaps.json` had gained
> an absolute `cohort_path` under the author's home directory. It was fixed at
> source and the artifact regenerated, not hand-edited.
>
> ✅ **CLOSED 2026-09-08 by a push, and verified from outside rather than
> asserted.** For three sessions the Data Availability Statement was true of this
> working tree and FALSE of the artifact a reader could reach. Established
> 2026-09-07 by querying the live repository's *contents* for the first time —
> every earlier check asked only for its refs, and a ref cannot reveal a missing
> file. At commit `8b7cce1` the published tree held **169 files, 1,448,379
> bytes** against a rebuild's **173**, and the four absent were `README.md`,
> `pipeline/results/gene_set_sizes.csv` and both `null_draws.npz` — **exactly the
> three files added to make this statement true**, plus the landing page.
> Measured against a tree reconstructed to match that commit, `10_make_figures.py`
> **exited 1** at figure 2 with figure 1B silently degraded to a mean ± SD band
> before it, while the nine supplementary tables built (exit 0) — so the table
> half was true and the figure half false.
>
> The author pushed at 18:08 CDT on 2026-09-08. The deposit was then commit
> `2b3bb75`, **173 files, 1,827,428 bytes**, and matched a fresh build with **zero
> file-set difference and zero content drift**. What closes this is not the push
> but the check on it: `22_build_public_snapshot.py --compare-live` fetches the
> published tree and diffs it against a fresh build, and it **exits 0**. Then the
> consumers were run *inside a clone of the pushed tree*, because running a
> producer is not the same as reading it: `10_make_figures.py` exits 0 and renders
> all five figures, reporting `panel B null: violins` and
> `figure2: panel sizes read from results/gene_set_sizes.csv (data/ is absent —
> this is the deposited-snapshot path)`, reproducing ρ = −0.500, p = 0.0485;
> `21_provenance_manifest.py` exits 0 with *120 of 124 verified; 4 not staged*;
> `25_build_supplementary.py` exits 0 with all ten tables. **Both defects this
> section describes as fixed are now genuinely fixed for a reader.**
>
> The general shape outlives the fix, which is why it stays written down: **the
> staged tree and the published tree are different artifacts, and until
> 2026-09-07 every check in this repository measured the staged one.** That gap
> is now a single command rather than an idea someone has to have, and it is the
> check to run before believing any statement in this section again.

**Variant A — public GitHub repository. ← CHOSEN**

> All analysis code is publicly available at
> `https://github.com/sathvikloke/immune-specificity-index` under the
> MIT license. The repository contains the complete pipeline, the frozen
> configuration used for every reported result, an exact dependency lockfile,
> a single-command reproduction script, and the pre-registered protocol. The
> reproduction script requires the public input data listed above to be
> downloaded first; those inputs are large and are not redistributed, so the
> script's verification mode does not run from the deposit alone.
> The frozen primary results were produced on macOS 26.5.2 / arm64 with
> Python 3.13.9; the corrected-ordering pan-cancer run and several later
> sensitivity analyses were produced on Linux/x86_64 with Python 3.12.14. See
> Limitations for the observed cross-platform spread.

⚠ **A tag claim was REMOVED here on 2026-09-07 and the gap is an OPEN AUTHOR
DECISION.** The paragraph above previously ended "The pre-registered protocol is
timestamped at git tag `prereg-2026-08-18`" — untrue of the repository it points
at: `git ls-remote --heads --tags` returns one ref and **zero tags**. The cause
is structural, not an oversight, so it cannot be patched here: the snapshot
builder publishes a *fresh* repository so no private string survives in history,
and a repository with no history has no tags. **Tagging the public repo is
ruled out** — its only commit is dated 2026-09-07. The protocol's *content* is
public; only its *timestamp* is private. **The three options, the trade-offs and
the recommendation (Zenodo, i.e. Variant B below) live in the project's working
task list (`16-YOUR-TASKS.md`, a planning document that is deliberately **not**
part of this deposit) — this manuscript is not the place to hold an unresolved
decision, and duplicating it invites the two copies to drift.** Nothing above
this line needs editing when it is settled
except, under Variant B, adding the DOI.

**Variant B — archived Zenodo snapshot with a DOI.**

> All analysis code is archived at Zenodo under DOI
> [10.5281/zenodo.XXXXXXX — ⚠ to be minted before submission] and released under
> the MIT license. The deposit contains the complete pipeline, the frozen
> configuration used for every reported result, an exact dependency lockfile,
> and a single-command reproduction script, captured at the commit corresponding
> to this manuscript. As in Variant A, that script needs the public inputs
> downloaded first and its verification mode does not run from the archive alone.
> The pre-registered protocol is included in the deposit, and
> the deposit's Zenodo publication date establishes its timestamp independently.
> The frozen primary results were produced on macOS 26.5.2 / arm64 with Python
> 3.13.9; the corrected-ordering pan-cancer run and several later sensitivity
> analyses were produced on Linux/x86_64 with Python 3.12.14. See Limitations
> for the observed cross-platform spread.

⚠ **Variant B carried the SAME false tag claim until 2026-09-07 and it is fixed
above.** It read "The pre-registered protocol is included as git tag
`prereg-2026-08-18`, created before any result in this paper was computed" — the
identical sentence removed from Variant A, and false for the identical reason
the moment the deposit is built from the public snapshot, which has no tags.
Recorded because of how it was found: the Variant A fix was verified against the
live repository and stopped there, leaving the unchosen variant to reintroduce
the defect if it were ever chosen. **A claim removed from one place is not a
claim removed.** Note also what Variant B actually buys: Zenodo's own
publication date is a third-party timestamp, but it is the date of *deposit*, not
of pre-registration — so it evidences "fixed before these results were
published", not "fixed before they were computed". That is weaker than the
private tag and stronger than nothing, and the Limitations wording should say so
rather than implying a tag exists.

Variant B satisfies the requirement without making the repository public, and
Zenodo appears on AACR's own list of acceptable repositories. Variant A is
simpler and is what the 3.2% of accepted computational abstracts claiming code
availability generally do. Both name the platform, because Limitations item 8
reports that the estimate was platform-dependent in its second decimal and a
reproduction claim that does not say where it was measured is not a claim.

## Authors' Contributions

**S. Loke:** Conceptualization, data curation, software, formal analysis,
validation, investigation, visualization, methodology, writing–original draft,
writing–review and editing.

**S. Kodilkar:** [NEEDS AUTHOR — CRediT roles].

**N. Movva:** [NEEDS AUTHOR — CRediT roles].

**M. Hota:** [NEEDS AUTHOR — CRediT roles].

**A. Raut:** [NEEDS AUTHOR — CRediT roles].

**S. Raut:** [NEEDS AUTHOR — CRediT roles].

**S. Bestavemula:** [NEEDS AUTHOR — CRediT roles].

*AACR requires a CRediT-style contributions statement for every author; the roles
above are the CRediT terms that apply. The six placeholders are deliberate and
must be filled before submission — a contributions statement is a factual claim
about what each named person did, and roles were not assigned on anyone's behalf.
Author order is the final byline the corresponding author gave on 2026-09-23.
(Until session 59 this sentence said the order followed the list of 2026-09-09
and was not yet confirmed; the list and its order both changed on 2026-09-23.)*

> **THIS SECTION WAS A SOLE-AUTHOR BLOCK UNTIL 2026-09-09.** It named S. Loke
> alone and stated "*Sole-author submission*". Six authors were confirmed on that
> date, with a seventh possible *[added 2026-09-23: the seventh, S. Bestavemula,
> joined the final byline that day]*. Every other place the sole-author framing
> survives is tracked in `16-YOUR-TASKS.md`, which is not part of this deposit;
> `CITATION.cff`, the poster byline and the portal author fields were still
> single-author when this note was written. **ORCID is required for every author
> before the copyright and COI forms can be completed, and three of the six iDs
> were still uncollected.**

## Authors' Disclosures

No disclosures were reported. [NEEDS AUTHOR — this sentence speaks for all seven
authors; confirm it with each co-author, whose disclosure forms are not yet
collected.]

*Wording taken from what* Cancer Research *actually prints, verified against
published articles rather than from a third-party guide. The verification record
is in the project's submission checklist, which is not part of this deposit.*

## Acknowledgments

This study used only publicly available data; no funding was received for this
work, and no third party contributed to its design, analysis or interpretation.

**Use of generative AI.** AI (Claude) was used for code editing and manuscript
editing. The corresponding author takes full responsibility for the content of
this work. No AI system is listed as an author, in accordance with AACR policy,
which requires such use to be declared both in the cover letter and in the
Acknowledgments.

> **THIS SECTION EXISTED NOWHERE UNTIL 2026-09-04.** The cover letter stated that
> AI use "is declared in the Acknowledgments" while the manuscript had no
> Acknowledgments section, no Authors' Contributions and no Disclosures — an
> assertion about the manuscript that the manuscript did not support. The cover
> letter's own §7 had noted the gap without it being closed. Verify the CRediT
> role list and the disclosure wording against the *Cancer Research* author
> instructions before submission.

## References

**Bibliographic metadata for references 1–12 re-verified against PubMed on
2026-09-03** — authors, journal, year, volume, issue, page range and DOI
confirmed against the indexed record. Reference 13 is not indexed in PubMed (the
query `Bengio Y[Author] AND Grandvalet Y[Author]` returned no record on
2026-09-16); its title, authors, volume, pages and year were confirmed against
the publisher's page the same day. This note covers the METADATA only; see
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
13. Bengio Y, Grandvalet Y. No unbiased estimator of the variance of K-fold
    cross-validation. *J Mach Learn Res* 2004;5:1089-1105. Not indexed in
    PubMed. [Publisher](https://www.jmlr.org/papers/v5/grandvalet04a.html)

**Cross-reference check, 2026-09-04.** **Every entry above is cited in the
body.** The previous reference 9 (Rooney et al., immune cytolytic activity) was
listed but cited nowhere, and the manuscript contains no cytolytic-activity
discussion for it to attach to -- a leftover from an earlier draft. It was
**dropped by authorial decision on 2026-09-04** rather than given an invented
citation site. That left twelve entries, every one cited and every one
PubMed-verified. On 2026-09-16 Bengio and Grandvalet was added as reference 13,
first cited in Limitations item 11, so the list is now
**thirteen references, all cited**, and all thirteen verified at source —
twelve against PubMed, one against its publisher.

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
