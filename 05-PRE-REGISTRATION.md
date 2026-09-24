# Pre-registration — H&E immune-specificity audit

**Filed 2026-08-18.** Target venue: AACR Annual Meeting 2027 (abstract deadline
Tue 2026-11-10) plus a full paper.

---

## 0. HONEST STATUS — read this first

**The discovery analysis has already been run. This document does not, and must
not, claim to pre-register it.**

By 2026-08-18 the NSCLC cohort (n=944, TCGA LUAD+LUSC) had been analysed
end-to-end five times and the pan-TCGA cohort (n=7,168) twice, and analysis
choices — the axis method, the disattenuation form, the head-refitting symmetry,
the bootstrap unit, and outcome stratification by cancer type — were changed
*after* seeing results, as documented in [04-CODE-AUDIT.md](04-CODE-AUDIT.md). Every one
of those changes was made to fix a defect rather than to move a number, and each
is recorded with the measurement that motivated it. But that is a
defect-driven revision history, not a blind protocol.

Therefore this document does two separate things, and the manuscript must keep
them separate:

| | status |
|---|---|
| **Part A — TCGA discovery result** (NSCLC n=944 and pan-TCGA n=7,168) | **EXPLORATORY.** Timestamped here so the analysis state is fixed from now on, and so any later change is visible as a change. Must be labelled exploratory in the abstract and the paper. |
| **Part B — CPTAC external validation** | **CONFIRMATORY.** The data has not been touched. No CPTAC slide has been downloaded, no feature extracted. This part is a genuine pre-registration and is the one that carries inferential weight. |

Anyone reviewing this project should hold Part A to the standard of a
well-controlled exploratory analysis and Part B to the standard of a
pre-registered confirmatory test.

---

## 1. Question

Do H&E foundation-model embeddings carry **immune-specific** information about
the tumour microenvironment, or do they read a dominant global expression axis
(purity, stromal fraction, proliferation) that essentially every gene set of
appreciable size loads onto?

The two hypotheses predict the *same* raw correlation, which is why a bare
random-set null cannot separate them — and why HE2RNA's published null
(Schmauch et al., Nat Commun 2020: 10,000 random lists, 75%/86% of cancer types
favouring immune signatures) does not settle the question. Two corrections are
required and are the contribution:

1. condition the **signature** on the global axis (not, as in the usual purity
   control, condition the image *prediction* on purity);
2. disattenuate for **reliability** — curated sets are co-expressed modules and
   are therefore more predictable than random draws regardless of biology.

## 2. Primary endpoint (byte-identical to `AuditConfig.primary_endpoint`)

> immune-specificity index (ISI): reliability-disattenuated,
> global-axis-residualised excess in Fisher z of r(image, signature) over the
> mean of r(image, size-matched random set), 1000 draws, patient-clustered
> bootstrap CI

Operationally, per signature *s*:

```
ISI_s = z[ r(pred_s , signature_s | axis) / sqrt(alpha_s) ]
        - mean_j z[ r(pred_sj, random_sj | axis) / sqrt(alpha_sj) ]
```

- `pred` is out-of-fold under **preserved-site** CV (Howard et al.,
  Nat Commun 2021;12:4423). The observed head and every null head are refit
  identically on the residualised target.
- `alpha` is Cronbach's α **of the residualised score**, general form
  α = k/(k−1)·(1 − Σvar(itemⱼ)/var(Σitemⱼ)). *Not* the raw-score closed form —
  the two differ 5.3-fold in NSCLC and 13-fold pan-cancer, see
  [04-CODE-AUDIT.md](04-CODE-AUDIT.md) §A.
- Pooled across signatures as a plain mean; CI by bootstrap over **patients**
  (944 resampled per draw, shared across signatures), null subsampled to 200
  draws per signature and that subsampling reported.
- Declared positive iff the pooled CI excludes 0.

**Fixed in advance and not to be changed:** `axis_method="pc1"`,
`axis_within_type=True`, `scorer="mean_z"`, `n_null_sets=1000`, `n_folds=5`,
`seed=0`, `match_null_expression=True`, `disattenuate=True`. The axis method in
particular is **not** a free sensitivity analysis — `meanz` collapses random-set
reliability to 0.025 and makes the estimand undefined ([04-CODE-AUDIT.md](04-CODE-AUDIT.md) §G).

## 3. Secondary endpoints

- **Split robustness (`delta_MAE`)** — paired difference in mean absolute error
  between random-patient and preserved-site CV, patient-clustered bootstrap.
  This is a calibration-loss contrast and must **not** be reported as "Δr".
- **Outcome arm** — Venet's original question against PFI (TCGA-CDR; Liu et al.,
  Cell 2018;173:400-416). Continuous scores into Cox; **no dichotomisation path
  exists in the code**; reduction to scalar fixed in advance as Harrell's C.
- **Three-way report** — r_raw / r_axis / r_residual per signature, the arm that
  separates "no immune signal" from "the image reads composition faithfully".

## 4. Multiplicity

Benjamini–Hochberg FDR at q=0.05 across the 16-signature family. `effective_tests`
(Li & Ji) is reported **descriptively only** — dividing α by M_eff *and* applying
BH is a category error and is not done.

## 5. Negative controls, all mandatory

Raw site-prediction AUROC of the embeddings; baselines (cohort_mean, type_mean,
site_only, type_only, purity_only, covariates); purity decomposition on one shared
complete-case sample with adjusted R². Post-ComBat site AUROC is **not**
interpretable in either direction and is never reported as success; the deployment
statement is the estimable fraction, which under preserved-site folds is 0% by
construction.

## 6. Part A — exploratory result as of this filing

Two cohorts, Prov-GigaPath slide embeddings (CC-BY-4.0), MSigDB Hallmark TME
sets (CC-BY-4.0), B=1000. Outcome concordance is computed WITHIN cancer type.

```
                        NSCLC n=944            PAN-TCGA n=7,168 (31 types)
PRIMARY   ISI     0.3182 [0.2614, 0.3763]   0.2911 [0.2692, 0.3133]
          beat null      16/16                     16/16
3-WAY     r_raw/axis/resid  0.400/-0.023/0.338     0.630/0.067/0.242
OUTCOME   beat null       0/32                     10/32
          MDE (C-index)   0.046                    0.017
CONTROL   site AUROC      0.992                    0.998
RELIAB.   mean gap       +0.1672                  +0.1843
```

**The two outcome results are consistent, not contradictory.** NSCLC could not
detect a C-index advantage below 0.046; the real pan-cancer advantage is ~0.03.
The NSCLC null was underpowered, and it must not be reported as a demonstration
that no signature is prognostic.

**What is prognostic is proliferation and stroma, not immunity.** Pan-cancer, the
five signatures beating their null are G2M checkpoint, E2F targets, angiogenesis,
EMT and hypoxia. No interferon, inflammatory, complement or allograft-rejection
signature does (0 of 6 MSigDB-immune vs 5 of 10 others, Mann-Whitney p=0.016).
An earlier pass reported p=0.00008 using a hand-made immune grouping; MSigDB
curates COAGULATION and COMPLEMENT as *immune* (Liberzon 2015 Table 1) and the
corrected test is the one above.

Interpretation fixed now, before validation, so it cannot be re-fit to the
answer:

- The image predicts curated TME signatures **specifically**, and this replicates
  across 31 cancer types with a tighter interval than the single-disease cohort.
  The gap is not reliability (the correction is applied and the finding survives)
  and not the global axis (removing it costs little).
- **The prognostic content of TME signatures is proliferation and stroma, not
  immunity.** Image-inferred immune scores should not be used as prognostic
  biomarkers — the image reads them accurately, they are simply not prognostic.
- The NSCLC 0-of-32 is an **underpowered null** and is reported as such, never as
  evidence that no signature is prognostic.
- Both halves are publishable; neither requires a clinical claim.

## 7. Part B — CONFIRMATORY validation (data untouched)

**Cohort.** CPTAC-LUAD + CPTAC-LSCC, 455 subjects, 2,218 WSI, 845.5 GB
(TCIA, verified 2026-08-12). No TCGA patient overlap.

> ### Resolved 2026-08-18 against the actual TCIA manifest
>
> `02-PROJECT-DECISION.md` §6 dropped CPTAC as "qualification-workflow
> specimen-segment" slides; `00-PROJECT-BRIEF.md` §8a re-budgeted it treating
> 4.9 slides/subject purely as a cost lever. Neither had read the manifest.
> `scripts/07_triage_cptac.py` now does. TCIA's "qualification workflow" wording
> is confirmed verbatim, and the discovery side is **99.48% DX** (10,117 of
> 10,170 slides, parsed from the embeddings parquet). *[Corrected 2026-09-17,
> in this addendum only; the tagged protocol never carried the figure. All
> 10,170 are diagnostic slides from primary-tumour samples: the 53 the parse
> missed carry lettered suffixes, DXA to DXU, which a digit-only pattern does not
> match. They belong to 7 patients (SARC, THYM, LGG, TGCT).]* What the manifest
> adds:
>
> | | slides |
> |---|---|
> | all slides in the cohort table | 2,138 |
> | − normal tissue | −765 (**35.8%**) |
> | − tumour, segment failed CPTAC's own QC | −235 |
> | **= usable** | **1,138** from **439 cases** (LUAD 615/229, LSCC 523/210) |
>
> **The decisive objection is not slide type — it is selection on the outcome.**
> Accepted segments were chosen for **median 70% tumour nuclei (IQR 60–75)**.
> This project's target is stromal/immune signature score, so CPTAC's QC selects
> directly on the compartment being predicted, depleting stroma and immune
> infiltrate and compressing the range of the target. That attenuates the
> correlation by construction.
>
> **Registered consequence, fixed now:** a surviving ISI in CPTAC is a
> **conservative** result and counts as replication. A **null** ISI in CPTAC does
> **not** refute Part A, because attenuation is predicted in advance — it will be
> reported as inconclusive, and the observed range restriction (sd of the
> signature scores in CPTAC vs TCGA) will be reported beside it. Committing to
> that asymmetry *before* seeing the data is the point of writing it here.
>
> ### RESOLVED 2026-08-20: the slides are FFPE. Part B is live.
>
> The decisive source is the STAR Methods of the CPTAC pan-cancer histopathology
> paper, which used **these exact TCIA slides** (2,755 slides / 657 patients
> across six types including LUAD and LSCC):
>
> > "Images consisted of control diagnostic H&E-stained slides sectioned from
> > **formalin-fixed paraffin-embedded (FFPE)** tissue blocks obtained from the
> > Clinical Proteomic Tumor Analysis Consortium (CPTAC), representing
> > mirror-image sections matched to the cognate OCT embedded tissue samples used
> > for molecular analyses."
>
> Wang, Hong, Demicco, et al. *Deep learning integrates histopathology and
> proteogenomics at a pan-cancer level.* Cell Rep Med 2023;4(9):101173.
> PMID 37582371, [DOI](https://doi.org/10.1016/j.xcrm.2023.101173). That paper
> separately partitions TCGA into FFPE vs frozen and reports the performance
> penalty for frozen, so its authors demonstrably track the distinction rather
> than using "FFPE" loosely.
>
> Corroborated by the CPTAC Biospecimen Core Resource SOP (Van Andel;
> [DOI](https://doi.org/10.1186/s12014-024-09450-3)), which states the CPTAC-wide
> rule that the larger piece is snap-frozen while "the corresponding mirror-image
> sections were formalin-fixed paraffin-embedded (FFPE) and H&E stained" — and
> whose named QC fields (viable tumour nuclei, total cellularity, necrosis) are
> exactly the columns in the TCIA manifest. The manifest **is** that FFPE
> qualification review.
>
> **There is no FFPE-vs-frozen choice to make.** Verified directly against the
> manifest: 2,138 rows, 2,138 distinct Slide_IDs, 2,138 distinct Specimen_IDs,
> **zero specimens carrying more than one slide.** Slide number = specimen number
> + 20. Filenames carry no `-DX`/`-TS`/`-BS` analogue of the TCGA convention
> because there is only one kind of slide.
>
> **Contradicting evidence, recorded rather than suppressed.** A 2022 *Diagnostics*
> paper asserts CPTAC slides are frozen — but for **UCEC, not lung**, from a group
> outside CPTAC, and contradicted by the CPTAC-affiliated group that assembled
> that image set. Separately, GDC lists CPTAC-3 preservation as snap-frozen with
> zero FFPE; that registers the **molecular aliquots**, and GDC hosts no CPTAC-3
> slide images at all. Neither displaces the primary source.
>
> **Residual risk, ~5%.** 19 slides were inspected visually (via the public NCI
> Imaging Data Commons, no full WSI downloaded) and showed paraffin ribbons and no
> freeze artefact. All 1,138 usable slides were not individually checked. LUAD's
> methods mention OCT embedding for a subset where LSCC's do not, so if a cheap
> further check is wanted, screen LUAD thumbnails first.
>
> **Also drop:** 4 of 2,138 slides (0.19%) are legacy CPTAC-2 cases with UUID
> filenames and a `_D1_D1` suffix (11LU013, 11LU016, 11LU022, 11LU035). Excluding
> them keeps the set homogeneous.
>
> **Revised budget:** 1,138 usable slides ≈ 434 GB ≈ 78 GPU-hours ≈ **$27–54** on
> a 4090. Filtering to usable slides first roughly halves both download and
> compute against the §8a figure.

**Committed before any data is touched:**

1. **Every USABLE slide is extracted** — `Specimen_Type == "tumor_tissue"` AND
   `Tumor_Segment_Acceptable == "Yes"`, i.e. 1,138 slides from 439 cases — and
   slides are averaged to patient level by `collapse_to_patient()`, identical to
   the TCGA arm. Normal-tissue slides and QC-failed segments are excluded by this
   rule, written down before extraction. Selecting one slide per patient is
   cheaper but introduces a post-hoc choice of *which*; that degree of freedom is
   closed here. (`scripts/07_triage_cptac.py`, `scripts/04_budget_cptac.py`.)
2. Prov-GigaPath, `layer=-1`, tile_size 256 at level 0, sample code `01` only.
   Both sides of the comparison use settings we control, which resolves the
   week-1 magnification/pooling ambiguity by construction.
3. Signatures are z-scored **within cohort**. This is a real choice: scoring TCGA
   and CPTAC on a pooled scale would leak; scoring separately redefines the
   target. Within-cohort is fixed now.
4. Success criterion, single and stated in advance: **the pooled ISI CI in CPTAC
   excludes 0, and its point estimate lies within the TCGA CI [0.2614, 0.3763].**
   A positive ISI outside that interval is reported as replication of *direction*
   but not of *magnitude*.
5. The outcome arm is re-run identically. The pre-registered expectation is that
   it again fails to beat the null; a positive outcome result in CPTAC would be
   the surprise and would be reported as such.
7.6 **Fallback.** If extraction is not funded, not completed before the deadline, or
   blocked by the slide-type question above, the abstract
   reports Part A **as exploratory, without validation**, and says so. It does not
   silently drop Part B.

## 8. Deviations

Any change to this protocol after filing goes in a dated appendix to this file,
with the reason, before the affected result is reported. `AuditResult.save()`
asserts that the estimand actually computed matches
`config.primary_endpoint`, so the code and this document cannot drift apart
silently.

## 9. Reproduction

```bash
python scripts/05_run_nsclc.py --outdir results/nsclc_v2
python -m pytest tests/ -q          # 67 tests
```

Every run writes `config.json` (the full frozen configuration) and
`summary.json` (estimand, CI, method string, and all run notes) into its output
directory.

## Appendix A — deviations after filing, recorded 2026-09-17

§8 requires each change to go here, dated and with its reason, **before** the
affected result is reported. The two changes below were made on 2026-08-19, the
day after filing, and the results they affect were reported without this
appendix. It was written on 2026-09-17, when the tagged protocol
(`prereg-2026-08-18`, commit `2d3f6fa`) was compared line by line with the code
that produced the frozen results; both changes are in commit `39664aa`.
Recording them late is itself a deviation from §8.

1. **The outcome arm's concordance is stratified by cancer type.** §3 fixed the
   reduction to a scalar as Harrell's C and said nothing about pooling across
   types. Run pan-TCGA, a C-index pooled over 31 cancer types mostly measures
   which cancer a sample is: cancer type alone reaches C = 0.676 on PFI, and
   unstratified, 94% of signatures "beat" their null. Only within-type pairs are
   now comparable whenever a cohort has more than one type, so NSCLC is
   stratified by LUAD and LUSC too. The change was made after the unstratified
   pan-cancer outcome had been seen. Every reported outcome figure, in both
   cohorts, is the stratified one; NSCLC's count, 0 of 32, is the same under
   both.
2. **The pooled bootstrap uses one complete-case mask across signatures.**
   Pan-TCGA, 19 of the 7,168 patients carry no cancer type and therefore no
   within-type axis value. One non-finite value made every bootstrap draw's
   correlation undefined while the point estimate masked it, so the interval
   came back empty. The bootstrap now drops those patients from every signature
   together, which keeps the patient-cluster structure, and reports the count.
   NSCLC has no such patient: its primary result, 0.3182 [0.2614, 0.3763], is
   the one §6 recorded at filing.

3. **The pan-TCGA cohort was added.** This protocol registers the NSCLC
   analysis (Part A) and a CPTAC validation (Part B); it does not mention a
   pan-cancer cohort. Pan-TCGA (7,168 patients, 31 cancer types) was added on
   2026-08-19 as a replication under the same estimand and the settings §2
   fixed, and it is where both changes above were found. Its results are
   therefore not pre-registered in the sense Part B is, and the manuscript no
   longer calls them so.

Neither change alters the primary estimand in §2: the point estimate, the null,
the disattenuation and the settings fixed in advance are as filed. The
manuscript's Methods state all three items and their dates.
