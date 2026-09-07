# Science audit — what is actually left

**2026-08-19.** Method: grep every analysis function in `aacr27` for whether any
runner calls it, then check every commitment in the pre-registration and the
brief against what the result files actually contain. Not from memory.

**Headline: the science is ~85% done, but three pre-registered items were
computed-and-discarded or never computed at all, and one of them is load-bearing
for the outcome arm.**

---

## A. Gaps between what was promised and what exists

### A1. BH-FDR is pre-registered and was never applied — **fix required**

[`05-PRE-REGISTRATION.md`](05-PRE-REGISTRATION.md) §4, verbatim:

> Benjamini–Hochberg FDR at q=0.05 across the 16-signature family.

Neither `immune_excess.csv` nor `outcome_arm.csv` carries a q-value. `beats_null`
is `excess_lo > 0` — an **uncorrected** 95% interval. The ancestry arm, written
later, does apply BH; the primary family is the one place it is missing.

**Does it matter?** For the ISI, almost certainly not — the smallest pan-cancer
`excess_lo` is +0.14, nowhere near the boundary. **For the outcome arm it may
decide the result.** Hypoxia beats its null with `excess_lo = +0.0069`. Five of
32 comparisons sit close to zero, and BH across 32 could remove some of them.
The claim "five signatures beat their null" is therefore **not yet
multiplicity-corrected**, and the corrected count may be smaller.

### A2. The family-level p-value is computed every run and thrown away

`stats.rotation_null` is called, and the result stored in
`immune_excess.attrs["rotation_null"]`. **pandas `.attrs` do not survive
`to_csv`.** So the rotation null — added specifically because a naive
sigma/sqrt(m) band around the family mean is 2–5× too narrow — is computed on
every run and silently lost.

### A3. `effective_tests` (M_eff) is promised descriptively and never computed

No runner calls it. The pre-registration says it is reported; it is not.

### A4. Control C was never run

[`00-PROJECT-BRIEF.md`](00-PROJECT-BRIEF.md) §4C describes label-side site
variance — regress each signature's ground truth on tissue source site *from RNA
alone, no image* — as one of the four design elements.
`decomposition.label_site_variance` exists and is tested. **No runner calls it.**
It needs no embeddings and costs minutes.

### A5. Partition variance is switched off

`report_partition_variance=False` in both runners. The function's own docstring
records the measurement that motivated it: across 15 site-disjoint partitions of
NSCLC, Δr had sd 0.019 against a patient-bootstrap half-width of 0.048 — **about
27% of the honest interval is missing** when partition choice is ignored. The
reported CIs are therefore too narrow by an amount that has been measured and not
applied.

### A6. The headline methodological claim is extrapolated, not measured

*"Small clinical panels are most exposed"* rests on a reliability gap measured
across **k = 36 to 200**. No gene set below 36 was tested. The claim about 5–20
gene panels — the ones nearest deployment, and the reason the finding matters —
is an extrapolation beyond the measured range. This is the weakest link in the
most transportable result, and it is **cheap to close** by subsampling the
Hallmark sets down to k = 10 and 20.

### A7. Scorer sensitivity never run — **CLOSED 2026-09-03**

`scorer="ssgsea"` was implemented and unused. Everything was mean-z, and whether
the ISI survived a rank-based scorer was untested. ssGSEA is slow, so this was a
real cost rather than a reporting fix.

**Now run on macOS, four arms, durable artefact on disk.** The answer is not a
simple yes or no: the effect survives the change of scorer on the comparable
uncorrected contrast, while the *registered* estimand turns out to be undefined
under ssGSEA because the reliability estimator — not ssGSEA — breaks down. Full
result, numbers and provenance caveat in §C below.

### A8. External validation — BLOCKED BY A DESIGN PROBLEM, not deferred

**This section said "deferred deliberately" until 2026-09-04. That was wrong, and
the correction is a result rather than a status change.** CPTAC is not waiting on
budget, hardware or priority. It is blocked because **the cohort does not carry
the variable the estimand is defined on**, and that has now been measured against
**three independent repositories** rather than inferred once.

**The estimand requires a site variable.** The ISI is defined under preserved-site
cross-validation (whole tissue source sites assigned to disjoint folds). Part B of
the pre-registration cannot be run on a cohort with no site variable — not
approximately, but by definition.

**Route 1 — GDC. Closed negative, measured 2026-09-04.** All **1,866** CPTAC-3
cases return `tissue_source_site = _missing`. The query is not at fault: the TCGA
positive control returns `TCGA-05-4384 → {"code": "05", "name": "Indivumed"}`.
No other institution-like field is populated either — `center.name`,
`center.short_name` and country are all empty, and `samples.submitter_id` is an
opaque numeric barcode.

**Route 2 — the Proteomic Data Commons. Also closed negative, measured
2026-09-04, and this closes the route handoff #23 called "the only one that
preserves the design".** The PDC GraphQL API at `pdc.cancer.gov/graphql` is
reachable and returns **6,808 cases**, of which **1,996** are
`CPTAC3 Discovery and Confirmatory`. Introspection is disabled server-side, so the
schema was mapped by submitting candidate fields and reading the validation
errors. The finding:

| PDC field | what it actually contains | usable as a site? |
|---|---|---|
| `tissue_source_site` | **does not exist** on type `PublicCase` | — |
| `tissue_collection_type` | **null for 146 of 146 samples** inspected across 100 cases | no |
| `biospecimen_anatomic_site` | organ — Kidney, Uterus, Pancreas | no, anatomical not institutional |
| `site_of_resection_or_biopsy` | `Not Reported` | no, and anatomical anyway |
| `preservation_method` | Snap Frozen / Frozen | no, technical |
| `pool` | Yes/No | no |

**Every site-like field in PDC is anatomical or technical. None is an enrolling
institution.** The one structure that does exist is the case-ID prefix, and it is
not a substitute: of the 1,996 Discovery cases, **881 are `C3L` and 738 are
`C3N`** — 81% falling into just two levels, which cannot support site-disjoint
folds. The remaining prefixes (`JHU` 39, `Ref` 32, `Kad` 30, `IPM` 23, `NX0` 20,
`KOC` 17) *do* look institutional, but they cover under 10% of the cohort and
would define a site variable for a minority while leaving the majority
unpartitioned.

**Route 3 — The Cancer Imaging Archive. Closed negative, measured 2026-09-05.**
TCIA serves nine CPTAC collections (`CPTAC-AML`, `-CCRCC`, `-CM`, `-LSCC`,
`-LUAD`, `-PDA`, `-SAR`, `-STAD`, `-UCEC`). Its patient metadata carries
**eight fields and not one of them is institutional**: `Collection`,
`PatientId`, `PatientName`, `PatientSex`, `EthnicGroup`, `Phantom`,
`SpeciesCode`, `SpeciesDescription`. There is no submitting centre, no source
site, no enrolling institution. `CPTAC-LUAD` also returns only **38 patients**,
so even a usable site variable would not reach the scale Part B needs.

**Route 3b — GDC's `center` entity, re-tested directly 2026-09-05 with a
control on the same query.** `center` is a *different* GDC entity from
`tissue_source_site` and deserved its own test rather than an inference. Asking
for `samples.portions.analytes.aliquots.center.{name,code,namespace}` across
CPTAC-3 returns cases carrying **neither key** — the fields are absent, not
empty. The same query against `TCGA-LUAD` returns both: `tissue_source_site =
{"code": "44", "name": "Christiana Healthcare"}` **and** `center = {"code":
"09", "name": "Washington University School of Medicine"}`. So the query is
right and the metadata is simply not there. The only structure CPTAC-3 exposes
below the case is the aliquot barcode (`CPT0123740006`, `CPT0244590002`), which
is an accession identifier, not an institution.

**The last remaining inference in this section was TESTED and is now a
measurement, 2026-09-06.** Until that date this paragraph read: "One candidate
route was NOT tested … the `cptac` Python package's harmonised clinical frames …
were not queried, because the package harmonises the CPTAC DCC and PDC clinical
files that Route 2 already examined directly, and a harmoniser cannot supply a
variable its sources do not contain. That is an argument, not a measurement, and
it is the only step in this section that is."

The argument was sound and the conclusion was right, but the author asked for the
most rigorous available route, and the most rigorous route was to stop arguing
and run it. `cptac` was installed into an **isolated virtualenv** — never this
environment, because the four pinned packages must not move — and queried
directly:

| | measured |
|---|---|
| sources exposing LUAD clinical data | **1 of 7** (`mssm`; `washu`, `bcm`, `broad`, `harmonized`, `pdc`, `umich` all return *"does not have clinical data for luad cancer"*) |
| LUAD patients in that frame | **111** |
| clinical columns | **124** |
| columns matching /site\|center\|institut\|tss\|submit\|source/ | **8** |
| of those, institutional | **0** |

The eight are `tumor_site`, `tumor_site_other`,
`specify_distant_metastasis_documented_sites`, `history_source` (×3),
`site_of_new_tumor` and `other_site_of_new_tumor` — five anatomical, three the
provenance of the clinical *record* rather than of the specimen. There is no
`tissue_source_site`, no `center`, no institution, and no submitter.

So gate 1 now fails on **four independent routes, all measured**, and this
section contains no inferences at all. The cost of removing the last one was one
virtualenv and about two minutes.

**Also measured: `Genomics_Available` in the TCIA manifest is stale, but not in
the direction that rescues the arm.** GDC shows open-tier RNA-Seq for **337** of
the 441 manifest cases against the manifest's own **217**. 337 is the ceiling, not
441: **104 cases have no open RNA-seq**, and matched RNA is a hard ISI
requirement.

**Consequence, and the author's decision, 2026-09-06.** Part B cannot be run on
CPTAC as written. Of the four options set out in `17-CPTAC-PLAN.md`, the author
chose **option 1: drop the external-validation arm and report the attempt as a
finding** — a substituted grouping variable would silently redefine the estimand,
and the honest sentence is more useful to a reviewer than a weakened analysis.
§§4–8 of that file are now design record rather than live plan; nothing is
deleted, because they stay correct if a site variable is ever found.

This is a reusable negative, not merely a local one: anyone proposing CPTAC as a
**preserved-site** validation cohort for an image/RNA model meets the same
missing variable, and four independent routes now say so by measurement.

**Do not build the extraction pipeline.** It is the expensive part and it has no
use under the recommended option.

### A9. The ISI is not reproducible across platforms — NEW, 2026-09-02

**Found while running A7 on the Einstein HPC4 cluster, and it is the most serious
open item in this file.**

| | NSCLC ISI |
|---|---|
| macOS 15 / arm64, Python 3.13.9 | **0.318240180054566** |
| Einstein HPC4 / x86_64, Python 3.12.14 | **0.2964** |

Same code. Same seeds. Same `requirements-lock.txt` versions actually installed
and verified at runtime: numpy 2.1.3, pandas 2.2.3, scikit-learn 1.6.1, scipy
1.15.3. Bit-identical inputs — `14_repro_probe.py` hashes the parquet/npy files,
the gene axis, the mean-expression vector, the qcut bins, the sampled null gene
sets and the mean-z scores, and stages 1–6 agree on both machines.

**The gap is 0.0218.** For scale, A5 measured partition choice — a source of
uncertainty this audit judged worth reporting and correcting the interval for —
at a standard deviation of **0.0089**. Platform choice moves the estimate
**2.4× more than partition choice does**, and unlike partition variance it is
not in the interval at all.

**Two candidates existed; one is eliminated.** `13_scorer_sensitivity.py`
recomputes the cohort's mean-z columns before auditing, which is a no-op on macOS
(drift 0.00e+00) but a 1.11e-15 change on HPC4 — so that run was not evaluating
the stored Y. Job 117403 removed the confound by running
`12_partition_variance.py` at `split_seed=0` instead: same `AuditConfig`, stored
columns, no rescoring anywhere. It still returned **0.2964**. Rescoring is not
the cause; the divergence is inside `run_audit`.

**The leading hypothesis is now ELIMINATED, measured 2026-09-02.** It was that
`experiment._select_alpha_per_fold` — 16 × 5 = 80 scalar RidgeCV selections off
`ALPHAS = np.logspace(-2, 5, 24)`, each then held *fixed* across that signature's
observed target and all 1,000 null targets — could flip one decision across
platforms. Adjacent grid points differ by ~1.96×, so a flip would not perturb a
prediction; it would roughly halve or double the shrinkage for one signature-fold
on the observed and null sides at once. The magnitude fit: 0.022 is ~1e13 times
machine epsilon, which is what a discrete branch flip looks like and is not what
accumulated libm noise looks like.

That argument was plausible and it was never measured. `RidgeCV` picks `argmin`
over the grid of leave-one-out CV error, so a decision can only flip if winner
and runner-up are closer than the error two BLAS/libm implementations disagree
by. `16_alpha_margin.py` measured that separation directly on macOS,
`rel_margin = (mse[runner_up] − mse[winner]) / mse[winner]`, over all 80
selections in 19,123 s:

```
min 1.034e-04   p10 7.917e-04   median 5.099e-03   max 1.125e-02
below 1e-06: 0 / 80        below 1e-12: 0 / 80        below 1e-14: 0 / 80
recomputed argmin reproduces the real alpha in all 80 finite selections
```

The closest of the 80 decisions sits **1.0e-04** from its runner-up — roughly
eleven orders of magnitude further than differing floating-point implementations
disagree by. Every selection also lands on grid index 16, 17 or 18 with the
runner-up immediately adjacent, which is what a well-separated interior optimum
looks like. **The 80 alphas cannot differ across platforms, so they are not the
cause of the 0.0218 gap.** This was decided on one machine, without the cluster.

**A seventh explanation is now excluded on macOS: package-version drift.**
`requirements-lock.txt` pins numpy 2.1.3, pandas 2.2.3, scikit-learn 1.6.1 and
scipy 1.15.3, and every diagnostic script *printed* the versions it ran under —
but until 2026-09-03 nothing *compared* them to the lockfile. Printing a version
is not verifying it, and a silent version drift produces exactly A9's symptom: a
number that moved while the code, the seeds and the inputs did not.
`scripts/20_check_versions.py` now checks it mechanically, and confirms this Mac
matches all four pins exactly. **It has not yet been run on HPC4**, so version
drift is excluded on one side of the comparison only — running it there is a
one-second addition to the bisect trip and should not be skipped.

Two supporting measurements, both worth recording:

- `17_alpha_margin_fast.py` — the same question reached without `run_audit`,
  which is what turned a 5-hour job into a 1-second setup — was run concurrently.
  **It has now finished all 80 rows (2026-09-03 03:05, PID 52573,
  `results/alpha_margin_fast_macos.npz`) and printed
  `ALPHA DIGEST fc299c9c0f47466e` — the identical digest `16_alpha_margin.py`
  printed over the same 80 alphas.** Two independent code paths, one running
  through `run_audit` in 19,123 s and one bypassing it in ~1 s of setup, agree
  **bit-for-bit on every alpha**. This is the formal closure of this elimination:
  the fidelity caveat in `17`'s docstring (it residualises the 16 observed
  columns stacked, where `run_audit` stacks one observed column with that
  signature's nulls) is now a **measured no-op**, not an argued-away risk, and
  the caveat has been retired in the file. The elimination no longer rests on a
  single script being trusted.
- numpy on the macOS machine reports **BLAS openblas 0.3.21, not Accelerate**.
  The standing assumption that the contrast was Accelerate-vs-OpenBLAS was
  wrong. If HPC4's numpy also ships the OpenBLAS wheel, this is one library's
  arm64 NEON kernels against its x86_64 AVX kernels.

**Remaining suspects, in order.** Everything downstream of the alpha, per
`experiment.py:1099-1140`: (1) the ridge solve itself,
`models.cross_val_predict_multi` — the only stage running a large dense
factorisation, hence the only one where kernel differences have room to
accumulate; (2) `stats.corr_ci`, which wraps `scipy.stats.pearsonr`; (3) the
patient-clustered bootstrap `_pooled_isi_bootstrap` — its PCG64 draws are
bit-exact across platforms by specification, so only its arithmetic is suspect.

`scripts/18_bisect_platform.py` discriminates between these. It dumps raw float64
intermediates per stage to a `.npz` — residualised scores, fold alphas,
out-of-fold ridge predictions, `corr_ci` r values, Fisher-z excess — plus a
standalone bootstrap probe on synthetic input that isolates the bootstrap from
every upstream stage. Run it on both machines, then `--compare` names the first
stage whose max absolute difference leaves the 1e-12 round-off floor. That stage
is where the divergence enters; everything after it is amplification.

**Status as of 2026-09-04: the macOS half is DONE and the experiment is still
unanswered.** `results/bisect_macos.npz` (2,182,700 b, 16/16 signatures,
8,525.3 s, pooled undisattenuated z-excess `0.32133084670199263`) was written on
2026-09-03 at 02:04. The HPC4 half has not run, because **SSH key access to the
cluster is still not installed** — verified failing again on 2026-09-04 with
`Permission denied (publickey,gssapi-keyex,gssapi-with-mic,password)`. This is
now the single thing standing between the project and a named cause: one
`sbatch` on HPC4 and one `--compare` converts seven eliminations into an answer.
**One half of a two-machine experiment is not a result, and no conclusion should
be written from `bisect_macos.npz` alone.**

### A9, new evidence 2026-09-04: the arm-by-arm macOS/HPC4 contrast

**An A9 measurement was sitting unread inside the A7 artefacts.** A7 ran the
identical four-arm script on both machines — `scorer_sensitivity_local/` (macOS)
and `scorer_sensitivity_hpc4/` (job 116899) — and nobody had ever compared them
arm by arm. Doing so costs no compute and narrows A9.

| arm | macOS | HPC4 | macOS − HPC4 | relative |
|---|---|---|---|---|
| `mean_z__disattenuated` | 0.318240180054566 | 0.2963918792070625 | +0.021848 | **6.87%** |
| `mean_z__uncorrected` | 0.3191678769197761 | 0.29886562351366824 | +0.020302 | **6.36%** |
| `ssgsea__uncorrected` | 0.1799529379910222 | 0.16761685471555202 | +0.012336 | **6.86%** |
| reconstructed `mean_z` | 0.3177350730216135 | 0.2960046563884893 | +0.021730 | **6.84%** |

**Finding 1 — every arm diverges, in the same direction, by nearly the same
*relative* amount.** Three of the four sit within 0.03 percentage points of
6.85%, across two different scorers and two different corrections whose absolute
magnitudes differ by a factor of ~1.8. `beats_null` is identical on both machines
for all four arms (16, 16, 0, 15). A divergence that is roughly *proportional*
rather than roughly *additive* is evidence that it enters upstream of the arms and
is then carried through them, rather than each arm diverging on its own.

**Finding 2 — the reliability/disattenuation stage is EXONERATED, empirically and
cross-platform.** Cronbach's α agrees to **2.2e-16** on the observed signatures
and **3.2e-15** on the null sets — round-off, over all 16 signatures. Split-half
reliability for mean-z is identical to **0.0** on the observed mean and
**1.7e-15** on the null mean. This matters procedurally: `18_bisect_platform.py`
deliberately does *not* dump the reliability stage, and the plan was to extend it
if every other stage came back clean. That contingency is no longer needed for
mean-z. (ssGSEA's reliabilities *do* differ, by ~1.5e-3 — a rank-walk statistic
recomputed on each machine, which is a separate and smaller question.)

**Finding 3 — the divergence is already fully present in the correlation stage,
and its shape is diagnostic.** Comparing per-signature values:

- **Observed** residualised correlation `r_residual_refit`: differs by up to
  **2.24e-02**, and **mixed in sign** — macOS is higher for 11 of 16 signatures
  and lower for 5. Mean signed difference +0.00729.
- **Null** mean correlation `null_mean_r`: differs by up to **1.84e-02** and is
  **higher on HPC4 for 16 of 16 signatures**, mean signed difference +0.01182.

Both push the excess the same way — macOS has a higher observed and a lower null
— which is why the arm-level gaps compound so consistently. The mixed-sign
observed side looks like accumulated floating-point divergence; the **perfectly
one-sided null side does not**, and that asymmetry is the sharpest unexplained
feature of A9. The null gene sets and their scores are bit-identical (Cronbach α
above, and `14_repro_probe.py` stages 1–6), so whatever differs happens *between*
the null scores and their correlations — i.e. in `cross_val_predict_multi` or
`stats.corr_ci`, exactly the two leading remaining suspects, and now with a
signed prediction the bisection must reproduce.

**Confound, stated rather than assumed away.** The HPC4 artefact came from a run
that rescored the cohort's mean-z columns (drift 1.11e-15 there, 0.00e+00 on
macOS). Cluster job **117403** already showed rescoring is not the cause of the
ISI gap — it ran `12_partition_variance.py --partitions 1`, which reads the
stored columns and never rescores, and still returned 0.2964 — so this contrast
remains informative. But this particular pair of artefacts does contain that
difference, and a clean test is what `18_bisect_platform.py` is for.

**Both of these are now done, 2026-09-06.** `20_check_versions.py` ran on HPC4
inside job 121438 and reports **`OK: every watched package matches the
lockfile`** — numpy 2.1.3, pandas 2.2.3, scikit-learn 1.6.1, scipy 1.15.3. So
package-version drift is excluded on *both* machines, and the seventh
explanation is fully closed rather than half-closed. What the same run also
printed, and what nothing had recorded before, is that the two BLAS builds are
not the same build: **macOS openblas 0.3.21** against **HPC4 scipy-openblas
0.3.27** (x86_64, Haswell, DYNAMIC_ARCH). Cluster job **118320** was also read at
last, after eight sessions — its contents are what opened the alpha question
below.

**What this does and does not threaten.** It does not touch the sign, the
significance, or the ordering: both platforms give a large positive ISI with
16/16 signatures beating their null. It does threaten the *precision* the project
claims, the "one-command reproduction" wording in the README, and any reader who
reruns this and gets a different fourth digit. Two honest responses, not mutually
exclusive: report the platform as part of the provenance, and widen the interval
to include this source the way A5 widened it for partitions.

**Consequence for pending work.** Until this is understood, new numbers generated
on the cluster are not directly comparable to the frozen macOS results, so any
result intended for the manuscript is produced on the manuscript's platform.

### A9 RESOLVED, 2026-09-06: the stage is named, and the elimination that blocked it was wrong

> **Read “A9 CLOSED, 2026-09-05” at the end of this block first.** The stage
> named here (the ridge solve) is *not* the cause. The cause is a non-stable
> sort in `splits.py:177` that builds a different partition on each machine.
> This section is kept because its zero-flip control is what exposed the real
> cause, and because the wrong conclusion is part of the record.

**The divergence enters at the per-fold ridge alpha selection.** Both halves of
the bisection now exist and `--compare` names the stage:

```
  stage                  max |A-B|      max rel   description
  S0_boot_synth          5.551e-17    1.400e-16   standalone bootstrap on synthetic input
  S2_obs_resid           1.010e-14    7.848e-11   global-axis-residualised observed score
  S2_null_resid          1.021e-14    1.297e-09   global-axis-residualised null scores
  S3_fold_alpha          1.515e+03    1.015e+00   per-fold ridge alpha   <-- DIVERGES
  S4_pred_obs            4.245e-01    1.919e+03   out-of-fold predictions, observed  <-- DIVERGES
  S4_pred_null           2.708e-01    1.618e+04   out-of-fold predictions, null      <-- DIVERGES
  S5_r_obs               2.239e-02    7.779e-02   stats.corr_ci observed r           <-- DIVERGES
  S5_r_null              3.730e-02    8.824e+00   stats.corr_ci null r               <-- DIVERGES
  S6_z_excess            4.008e-02    1.574e-01   Fisher-z excess per signature      <-- DIVERGES

  pooled_z_excess: macOS=0.32133084670199263  HPC4=0.2988636473916485
                   delta=+0.022467199310344133
```

(`--compare` prints the pooled line truncated to 12 decimals; the three values
above are quoted at the precision the two `.npz` files actually hold, so that
`19_check_numbers.py` can verify them against the artefacts rather than having
to exempt them.)

`results/bisect_hpc4.npz` came from cluster job **121439** (50.2 s of compute;
2,182,700 bytes, byte-for-byte the same size as the macOS dump, so the two files
have identical structure). S0 and S2 agree to round-off. **The bootstrap and the
residualisation stage are therefore eliminated as causes** — two more negatives.

**11 of the 80 observed per-fold alphas differ, every one by exactly one grid
step.** Only three values occur across both machines — 740.5685, 1492.4955,
3007.8825 — and `max |A−B| = 1515.3870 = 3007.8825 − 1492.4955`. This is a
discrete flip between adjacent points on the log grid, not a numerical drift.

**Elimination #6 was wrong, and it is worth being precise about how.** Six
sessions carried "THE 80 RIDGE ALPHAS CANNOT DIFFER ACROSS PLATFORMS — confirmed
twice, independently". Both confirmations were on **macOS**. `16_alpha_margin.py`
and `17_alpha_margin_fast.py` agreed with each other bit-for-bit, which
established that the fast path is faithful — not that the alphas are
platform-invariant. The invariance itself came from an *argument*: the smallest
relative margin between a selected alpha and its runner-up was 1.034e-04, which
is ~8 orders of magnitude above the ~1e-12 that differing BLAS/libm builds
disagree by, so no flip should be possible. The argument is clean and it is
false. The measurement was never taken across platforms because SSH was down for
eight sessions.

Taking it required one run. `17_alpha_margin_fast.py` — the *same script* that
produced the macOS digest, so no iteration-order assumption is involved — was run
on HPC4 as job **121438**:

| | ALPHA DIGEST (all 80 values) |
|---|---|
| macOS (`16` and `17`, agreeing) | `fc299c9c0f47466e` |
| HPC4 (`17`, this session) | `946033557a8ad08b` |

The per-fold margins differ too, and by far more than round-off: EMT fold 1 has
relative margin **1.034e-04 on macOS and 5.488e-03 on HPC4**, a factor of 53. So
it is not merely that a near-tie fell the other way; the CV score curve the
argmin runs over is materially different on the two machines.

Note that `18_bisect_platform.py`'s own output still labels S3 "(PROVEN
invariant by 16_alpha_margin)", and `17_alpha_margin_fast.py` still prints a
VERDICT ending "The 80 alphas CANNOT differ across platforms ... A9's leading
hypothesis is DEAD." **Both strings are now falsified by the runs that printed
them.** They are left in place for this session so the record shows what was
believed; correcting them is workstream 1's first task next session.

**BLAS thread count is NOT the mechanism — a new elimination (#9).** Every
macOS-vs-HPC4 contrast to date confounded platform with threading, because the
cluster sbatch scripts set `OPENBLAS_NUM_THREADS=8` while macOS runs are
forbidden from setting it. Job 121438 ran `17` twice on the same node, at 8
threads and at 1, and the two arms are **bit-identical** — same digest, same
margins to every printed digit. Reduction order from thread count is not the
lever.

**The signed prediction reproduced, which is what makes the measurement
trustworthy.** Recomputed from the two dumps:

| prediction (from the frozen per-signature numbers) | measured in the bisection |
|---|---|
| observed `r` higher on macOS in 11 of 16 | **11 of 16** |
| observed max abs difference up to 2.24e-02 | **2.2386e-02** |
| null `r` higher on HPC4 in 16 of 16 | **16 of 16** |
| null max abs difference up to 1.84e-02 | 2.7227e-02 — *larger*, and expected to be: the frozen figure is over 1,000 null draws, the bisection over 8, so its per-signature null mean is noisier |

### The alpha flip is NOT the root cause — measured the same day, and it corrects the paragraph above

The first version of this section, written an hour earlier, said the null-side
alpha selections "are not in the dump" and should be dumped next. **That was
wrong, and reading `18_bisect_platform.py:338-344` is what corrected it.** There
are no null-side alpha selections to dump. The alpha is chosen once per
signature-fold on the observed residual and then applied to the observed and all
null columns together:

```python
fold_alpha = experiment._select_alpha_per_fold(X, obs_resid, split)
stacked    = np.column_stack([obs_resid, null_resid])
P = models.cross_val_predict_multi(X, stacked, split, patients,
                                   fixed_alpha=fold_alpha)
```

That is the documented design — one penalty held fixed across a signature's
observed target and all its nulls — and it turns the open question into a
sharper one, because it supplies a free control.

**Eight of the sixteen signatures have identical alphas on both machines.** For
those eight, the penalty is the same, and the inputs agree to **9.77e-15**. If
the alpha flip were the origin of the divergence, their downstream stages would
agree to round-off. They do not:

| | measured over the 8 zero-flip signatures |
|---|---|
| max abs difference, out-of-fold predictions | **3.023e-01** |
| the same, relative to those predictions' own sd (0.1705) | **177%** |
| median abs difference, relative to that sd | **14.7%** |
| correlation between the two platforms' predictions | **r = 0.954 – 0.970** |

Two platforms' ridge predictions correlating at r ≈ 0.96, from identical
penalties and inputs identical to 1e-14, is not round-off. **The ridge solve
diverges on its own.** The alpha flips and the prediction divergence are two
symptoms of one instability, and "attack S3" — which is what `--compare` used to
print — would have sent the next session to the wrong stage. The script now runs
this control automatically rather than asserting the amplification story;
`_amplification_check` is the function, and it was written *because* the sentence
it replaces was wrong.

**The structural fact that most plausibly explains it, measured not assumed:**
X is (944, 768) and the split is 5 preserved-site folds, so **every training set
is about 755 samples against 768 features — fewer rows than columns.** A ridge
system in that regime is solvable but its conditioning is poor, and the two
BLAS/LAPACK builds are genuinely different: **macOS openblas 0.3.21 / arm64**
against **HPC4 scipy-openblas 0.3.27 / x86_64 Haswell DYNAMIC_ARCH**, both now
recorded by `23_export_environment.py` on the machine in question.

**What has NOT been established** is which part of the solve turns a 1e-14 input
difference into an O(0.1) prediction difference — the solver `Ridge` selects for
n < p, the GCV eigendecomposition inside `RidgeCV`, or the conditioning itself.
That is the next measurement and it is stated here as a question, not an answer.
The cheapest discriminating experiment: refit one zero-flip signature's folds at
the fixed known alpha with `solver=` pinned explicitly (`svd`, then `cholesky`),
on both machines, and see which pinning makes the platforms agree.

> **SUPERSEDED, 2026-09-05.** The question above is malformed and its premise is
> false: nothing amplifies 1e-14 to O(0.1), because the ridge solve is a
> *contraction*. The two machines were never solving the same problem — they
> built different cross-validation partitions. See **A9 CLOSED** below. The
> section is kept because the zero-flip control it introduced is what made the
> real cause findable.

### A9 CLOSED, 2026-09-05: a non-stable sort builds a different partition on each machine

**The cause is `splits.py:177`.** `_assign_sites_greedy` orders sites
largest-first with

```python
order = per_site.index.to_numpy()[np.argsort(-(sizes + noise))]
```

`np.argsort` defaults to `kind="quicksort"` (introsort), which is **not
stable**. With the default `jitter=0.0` the `noise` term is exactly zero, so
equal-sized sites are exactly tied and their order is decided purely by the
sort's internal tie-breaking — which is implementation-defined and differs
between numpy's arm64 and x86_64 kernels. **51 of the 68 NSCLC sites sit in a
tied size group.** A different order feeds the greedy a different sequence and
produces a different site→fold assignment, *while preserving fold sizes
exactly*, which is why `_assert_site_disjoint`, the degenerate-partition guard,
and the equal-byte size of the two bisection dumps all passed.

**Measured on both machines** (`scripts/24_split_determinism.py`), numpy 2.1.3 /
pandas 2.2.3 / sklearn 1.6.1 / scipy 1.15.3 on each:

| hash | macOS arm64 | HPC4 x86_64 | |
|---|---|---|---|
| `H_X_bytes` | `b5d70e2e64ef4ad6` | `b5d70e2e64ef4ad6` | identical |
| `H_tss` | `5292af3ef56ccbef` | `5292af3ef56ccbef` | identical |
| `H_sizes` | `1085bc3b226f462f` | `1085bc3b226f462f` | identical |
| `H_order_quicksort` | `9fd91de2c6f61b19` | `0e7c4c7488097d58` | **DIFFERS** |
| `H_order_stable` | `654c68f3e3aa46a2` | `654c68f3e3aa46a2` | identical |
| `H_FOLD` shipped | `e55ee4e8887b5167` | `a40c862b12d9a391` | **DIFFERS** |
| `H_FOLD` stable | `9f31e33c2caaf1e3` | `9f31e33c2caaf1e3` | identical |

Inputs are bit-identical; the partition is not.

**The millisecond reproducer settles it.** `tests/_synthetic_isi` uses 20
synthetic sites of exactly 10 patients each, so *every* site is tied:

| | macOS | HPC4 | gap |
|---|---|---|---|
| shipped (quicksort) | 1.036173699004058 | 0.9752094760755513 | 5.89% |
| stable sort | 0.9752094760755506 | 0.9752094760755513 | 7e-16 |

**Applying the stable sort on macOS reproduces the HPC4 value to 15 significant
figures.** The whole 5.89% "platform gap" was the sort order. Real
cross-platform floating-point noise is ~7e-16 — six orders of magnitude below
what had been attributed to it. This also explains why a probe with n=200 > p=40
reproduced "A9" despite being nowhere near underdetermined: conditioning was
never involved.

**The premise that framed A9 for six sessions is false.** Measured here:

- `cond(K + αI)` for the per-fold ridge system is **18–77**, not ~1e10. The raw
  kernel is exactly singular — column-standardising `X` puts the all-ones vector
  in `null(Xᵀ)` — but the *regularised* system actually solved is well
  conditioned.
- Perturbing the inputs by 1e-14 relative moves the predictions by **6e-15**, an
  amplification factor of ~0.25. The ridge solve **attenuates**.
- The observed 4.245e-01 output difference against a 1.010e-14 input difference
  implies an amplification of **4.2e+13**, which this system cannot produce.

So `models.cross_val_predict_multi` is exonerated, and so is BLAS. The 11-of-80
alpha flips are a *consequence* of the different partition — different training
sets select different penalties — not a cause. The session-26 reversal ("the
alphas CAN differ") was correct; the reason was not.

**Consistency check.** The NSCLC ISI gap, 0.3182 − 0.2964 = 0.0218, is 2.45×
the measured partition sd of 0.0089 (A5). A plausible draw. The
"platform difference" was a *partition-variance* draw all along, which is a
phenomenon the paper already documents.

**The fix is now APPLIED, 2026-09-05 (session 28), under option B.**
`splits.py:177` passes `kind="stable"`. Verified by re-running
`scripts/24_split_determinism.py --check`: the shipped fold hash is now
`9f31e33c2caaf1e3`, identical to the stable hash and to the value measured
independently on macOS arm64 and HPC4 x86_64; patients whose fold id differs
between shipped and stable is **0 of 944**, down from 216. On the millisecond
synthetic reproducer the shipped path now returns `0.9752094760755506`, matching
the HPC4 value `0.9752094760755513` to 15 significant figures.

**The frozen results were deliberately NOT re-run.** `results/nsclc_v3/` and
`results/pancancer_v3/` still hold the numbers produced by the pre-fix code, and
they remain the reported primary. That is option B, chosen over a re-freeze for
a reason: re-freezing would sever the link between the `prereg-2026-08-18` tag
and the reported run, and would replace the pre-registered number instead of
reporting both. Measured re-run cost, from the frozen `summary.json` files:
NSCLC 1,739.6 s (0.48 h), pan-cancer 17,115.6 s (4.75 h) — about 5.2 h idle, and
this project has measured a 35× contention spread, so far more under load.

**`24_split_determinism.py --check` is now a regression guard that FAILS.**
Until this session it merely printed a NOTE when the shipped partition differed
from the stable one, because the fix was deliberately absent and a hard failure
would have made the script unusable in `reproduce.sh`. It now exits 1. Verified
able to fail: reverting `kind="stable"` in a scratch copy produced exit 1, the
old fold hash `e55ee4e8887b5167`, and "216 of 944 patients are in a different
fold"; the source file was restored byte-identical afterwards.

### A10, NEW 2026-09-05: `signatures.py:212` was a SECOND live instance, and a larger one

**Handoff #27 carried this sort forward as "latent, not implicated on this
data", on the grounds that "continuous expression evidently produces no exact
rank ties" — inferred from downstream scores agreeing across platforms to
1.0e-14. The inference was wrong, and the source file said so.**
`signatures.py:206`'s own comment reads "Ties are common in log-expression
(every gene at the floor shares a rank)", and `score_ssgsea`'s docstring already
recorded "the real NSCLC matrix, where 85% of each row's values are tied". The
handoff contradicted the code it was describing. This is the third time in this
project's history that an argument has been carried as a measurement.

Measured by `scripts/26_sort_audit.py` on the real NSCLC matrix (944 × 41,046):

| quantity | measured |
|---|---|
| samples carrying at least one exact rank tie | **944 of 944** |
| tied rank entries | **35,394,448 of 38,747,424 (91%)** |
| depth-table hash, quicksort | `0bca973cfa26da89` |
| depth-table hash, stable | `096172c056cc7a34` |
| ssGSEA score cells differing between the two | **15,060 of 15,104** |
| signatures affected | **16 of 16** |
| max abs. score difference | **1.2824e-02** |
| mean per-signature score sd | **2.1154e-02** |
| max difference as a fraction of that sd | **0.606** |

The difference is not confined to the intermediate: it survives the O(k)
reduction. Tied genes carry identical weights, so a swap inside a tie group
cancels whenever both or neither of the pair is a set member, and moves the
score only where exactly one is — which, with 41,046 genes and ~200-gene sets,
happens constantly. **0.606 of a standard deviation is six orders of magnitude
above the ~7e-16 that is this pipeline's genuine cross-platform noise.**

**Scope — this never touched the frozen headline results.** They run
`scorer="mean_z"` (`results/nsclc_v3/config.json`), which does not reach
`_ssgsea_sample_tables` at all. It did reach the ssGSEA arms of the A7
scorer-sensitivity analysis. `kind="stable"` was applied 2026-09-05 together
with the A9 fix. Any tie convention is a valid ssGSEA — the statistic is
genuinely ambiguous under ties — so the fix pins one rather than finding a
"right" one; stable pins tied genes to column order, fixed by the input file.

**RESOLVED 2026-09-06 by re-running, not by argument.** The two ssGSEA arms were
recomputed under the pinned ordering into
`results/scorer_sensitivity_ssgsea_stable/` (PID 6503, detached, complete
02:43). The defect is real at the score level and immaterial at the level of
anything reported:

| | pre-fix (`scorer_sensitivity_local/`) | pinned (`..._ssgsea_stable/`) |
|---|---|---|
| ssGSEA uncorrected excess | 0.1799529379910222 [0.1391910971170491, 0.2159114981000094] | 0.17377773241035724 [0.13325771984752488, 0.21018161815588407] |
| signatures above their null | 15 of 16 pre-fix | 15 of 16 pinned |
| reconstructed disattenuated | −0.13333507569716976 | −0.13986374583533862 |
| split-half reliability, observed | 0.9336120947138573 | 0.9349444164608371 |
| registered disattenuated estimand | undefined | undefined |

The point estimate moves **0.00617520558066495** — 16.1% of the pre-fix
interval's own half-width, and about a quarter of the 0.026 that the same class
of defect moved the NSCLC primary. A 0.606-of-a-score-sd movement in 15,060 of
15,104 per-cell scores therefore aggregates to a 0.006 movement in the index
built from them. **Sign, significance, the 15/16 count and the "survives at
roughly half the magnitude" reading are all unchanged**, and the registered
estimand stays undefined for the Limitation-9 reason, which is a property of the
reliability estimator and not of the sort.

Both numbers are now authorities in `19_check_numbers.py` (`a7_ssgsea_unc`,
`a10ss_ssgsea_unc`) with the shift **derived** from the pair rather than stored,
so the paper and the poster cannot drift away from the two runs they describe.

**DECISION 2026-09-06: A11 does NOT go in the manuscript, and the reason is
recorded rather than assumed.** A11 and the join defect that its own fix
introduced are both bugs in `11_close_science_gaps.py`, a post-hoc gap-closing
script, not in the estimator. Neither changed a published number: the audit
already scopes A4 to pan-TCGA, A4b permutes within cancer type and is meaningful
only pan-cancer, and MANIFEST.md/CAPTIONS.md source S6 to `pancancer_v3`. A
Limitations item describes a limitation of the reported analysis; a caught and
fixed defect in an auxiliary script is an audit finding. It stays here.

What DID reach the manuscript is the fix, in §Reproducibility and in one
sentence only: a run now records the cohort it was fitted on. That is a claim
about what the released code does, and it is the property whose absence let A11
hide for a session.

**Cost, measured, against a docstring that was wrong by ~30×.**
`13_scorer_sensitivity.py` claims "~8 minutes per scorer per arm". The two arms
actually cost 6,288.9 s and 3,630.8 s — 2 h 45 m in total, on a contended
machine. Its *own* live estimate of the scoring step (523 s per arm) was sound;
the docstring's total was not, because the ISI machinery around the scoring
dominates the scoring it was extrapolated from. Do not quote the docstring.

**The test that was supposed to catch this could not.**
`test_ssgsea_optimised_matches_reference_loop` compares the optimised form to a
reference loop, deliberately including tie-heavy cases. Both sides used the
default quicksort, so the test pinned the optimised form to a
*platform-dependent* target: it passed on arm64 and on x86_64 while asserting
different numbers on each. Both sides now pass `kind="stable"`, so the assertion
names a fixed value.

**Cross-platform status: CONFIRMED ON TWO MACHINES, 2026-09-05.** HPC4 job
**121781** (`cluster/sort_audit.sbatch`, node cpu-850, x86_64, Python 3.12.14,
numpy 2.1.3 / pandas 2.2.3, BLAS scipy-openblas 0.3.27) reproduces A9's pattern
exactly:

| quantity | macOS arm64 | HPC4 x86_64 | |
|---|---|---|---|
| samples with ties | 944 of 944 | 944 of 944 | identical |
| tied rank entries | 35,394,448 | 35,394,448 | identical |
| depth hash, **quicksort** | `0bca973cfa26da89` | `b98fb11e701e8e98` | **DIFFERS** |
| depth hash, **stable** | `096172c056cc7a34` | `096172c056cc7a34` | **identical** |
| ssGSEA max abs. diff | 1.2824e-02 | 1.3040e-02 | differs (quicksort arm) |
| as a fraction of score sd | 0.606 | 0.618 | differs (quicksort arm) |

The inputs are byte-equivalent — both machines count the same 35,394,448 tied
entries — and the unpinned convention still produces different depth tables,
while the pinned one is identical. That is the A9 signature, and it settles the
axis the macOS-only measurement could not: `signatures.py:212` really was
platform-dependent, and `kind="stable"` really does remove it. The two
platforms' ssGSEA differences are themselves slightly different (1.2824e-02 vs
1.3040e-02) because both the quicksort scores and their sd move; the *stable*
arm is the invariant one, which is the arm now shipped.

The job's own exit status is worth recording so nobody re-litigates it: `sacct`
shows job 121781 as **FAILED, ExitCode 1:0**, and that is a **true positive, not
a job error**. HPC4 still holds the pre-fix `src/aacr27/` by design (nothing on
the cluster was overwritten), so cases A, A-deep and E correctly reported
DANGEROUS and `26_sort_audit.py` exited 1 as designed; `set -euo pipefail` in
the sbatch then stopped the script before its final `24_split_determinism.py`
step. The self-test passed on HPC4 first, including the two-way hash check.

**The remaining non-stable sorts, now measured rather than assumed:**

| call site | ties on real data | changes output | verdict |
|---|---|---|---|
| `stats.py:180` BH-FDR | forced 32/32 | **no**, even when forced | SAFE by algebra |
| `barcodes.py:130` `site_summary` | 51 tied site counts | yes (`cum_frac`) | **DEAD CODE** — never called |
| `models.py:339` site AUROC order | 0 of 15 | no (median only) | SAFE by consumer |
| `outcome.py:173,181` | — | — | already `kind="stable"` |

`stats.py:180` is safe for a reason worth stating, and it is *not* "no ties":
exactly tied p-values receive adjacent ranks *i* and *i+1*, and the reverse
running minimum in the step-up assigns them the **same** q whichever came first,
so swapping them permutes `order` and `q` together and `out_q` is bit-identical.
`26_sort_audit.py` forces all-tied and all-identical p-value families and
confirms max |Δq| = 0 rather than trusting that argument.

`barcodes.py:130` is a genuine instance of the defect — per-site patient counts
are more heavily tied than the sizes in A9, and `cum_frac` is order-dependent by
construction — that is harmless **only** because `site_summary` is never called.
Measured, not remembered: `26_sort_audit.py` greps every `.py`/`.md`/`.sh`/
`.ipynb` in the repo and the sole occurrence is its own `def`. If anything ever
calls it, it must be fixed first.

A related dead pair: `models.py:341-342` write `median_auroc` and
`frac_above_0.9` into `.attrs`, which nothing reads — and which would not
survive `.to_csv()` anyway, exactly how the rotation-null p-value was silently
lost in A2. The reported site-AUROC median is safe because
`experiment.py:249` recomputes it from the column.

**Both of the above were CLOSED on 2026-09-06, and closing the first exposed a
defect in the audit itself.** `barcodes.py` now pins `kind="stable"`, on the
principle that a latent instance of a known defect class is still an instance.
But pinning it changed nothing: `26_sort_audit.py` went on reporting case C as
DANGEROUS and exiting 1. The reason is worth recording, because it is the
"assert on content, not on sentinels" lesson wearing new clothes — `case_c()`
was the only case that never returned a `shipped_pins_stable` key at all, so
`verdict()` had no way to observe a repair and could only ever call it
dangerous. **A check that cannot register the fix it demands is not measuring
the shipped code.** With case C wired to the same `pins_stable()` helper the
other cases use, `26_sort_audit.py` exits 0 for the first time: A FIXED, B
SAFE, C FIXED, D SAFE, E FIXED. `--self-test` still passes, so the detectors
can still fire. `models.py:341-342` are simply gone.

### A11, NEW 2026-09-06: the science-gaps runner ignored the cohort it was told to use

**`11_close_science_gaps.py` took `--cohort` and `--expr` independently of
`--results`, both defaulting to the whole pan-TCGA tables.** Stages A4
(Control C), A4b (permutation calibration) and A6 (small-panel bracket) read
those defaults and only borrowed *signature names* from the results directory,
so they scored all 7,168 pan-TCGA patients no matter which cohort `--results`
named. A run launched as `--results results/nsclc_v3_stablesort` therefore
wrote **202-site pan-cancer numbers into an NSCLC directory, under a heading
that read `=== closing science gaps for nsclc_v3_stablesort ===`.**

It was found by comparison, not by a check: three values in
`nsclc_v3_stablesort/science_gaps.json` — `label_site_variance_median_r2`,
`control_c_calibrated_r2`, `control_c_perm_median` — are **bit-for-bit
identical to `pancancer_v3`'s to seventeen significant figures**, as is the
whole small-panel bracket. Two cohorts differing 7.6× in size cannot agree to
17 digits. The log settles it: it announces the NSCLC directory and then
reports `16 signatures, 202 sites`, where NSCLC has 68.

**No published number is wrong.** The audit scopes A4 explicitly to "pan-TCGA
over 202 sites", A4b permutes *within cancer type* (which only means anything
pan-cancer), and `MANIFEST.md` and `CAPTIONS.md` both source S6 to
`pancancer_v3`. What was wrong is the artefacts: three files in an NSCLC
directory describe a different cohort. **The stages that matter were never
affected** — A1's BH-FDR columns, A2's rotation null and A5's M_eff all read
from the results directory, so the 16/16-after-BH result stands.

The cost was real even so: the ~2.5 h that run spent in A4b recomputed a
pan-cancer number that already existed.

Fixed by giving the script a `--cohort-scope {run,file}`, defaulting to `run`,
which restricts the cohort-derived stages to the patients the run was actually
fitted on — recovered from `predictions.csv.gz`, **the only cohort fingerprint
a results directory carries**, since neither `config.json` nor `summary.json`
records which patients were used. That absence is why nothing could detect the
mismatch, so `science_gaps.json` now records `cohort_path`, `cohort_scope`,
`cohort_n_patients_in_file`, `cohort_n_patients_used`, `cohort_matches_run` and
`a6_n_patients`, and a mismatch prints loudly.

**The filter is a measured no-op for pan-cancer, which is why it is safe:**
`meta` INNER JOIN `expr` is exactly 7,168 patients, precisely `pancancer_v3`'s
run, so A4 and A4b still reproduce their frozen values. A6 is the exception and
is scoped rather than silently changed — it z-scores over all 9,781 patients
with expression, which is not even the pan-cancer run's 7,168, so reproducing
the stored A6 block requires `--cohort-scope file`. The decision is recorded
rather than buried: correctness by default, exact reproduction on request.

Pinned by `test_cohort_scope_restricts_to_the_run_and_records_the_mismatch`,
which asserts in **both** directions — that `run` narrows and that `file` is
reported as a mismatch. Verified able to fail: making the resolver ignore its
scope argument turns the test red.

## B. What is genuinely finished

- The primary estimand, in two cohorts, with the four defects from the code audit
  fixed and regression-tested.
- The outcome arm, stratified by cancer type, with MDEs attached.
- Negative controls: raw site AUROC, the ComBat estimability argument, covariate
  baselines, purity decomposition.
- The ancestry supplement, including the identifiability result.
- 104 tests, deterministic seeds, exact version lock, one-command reproduction.

## C. Verdict, in effort

| Item | Effort | Status |
|---|---|---|
| A1 BH-FDR | post-hoc | **CLOSED — and it changes nothing** |
| A2 rotation null persisted | minutes | **CLOSED — p = 0.000999** |
| A3 M_eff | minutes | **CLOSED — M_eff = 5.00 for 16 signatures** |
| A4 Control C | no embeddings | **CLOSED — see below** |
| A6 small-panel bracket | subsampling | **CLOSED — and it upgrades a headline** |
| A5 partition variance | ~15 min once the legacy null is off | **CLOSED — and it corrects a claim** |
| A7 ssGSEA sensitivity | hours of compute | **CLOSED 2026-09-03 — four macOS arms landed, acceptance check passed bit-identically** |
| A8 CPTAC | ~$50, 6 days, plus a slide-type check | deferred, stated |
| A9 cross-platform reproducibility | mechanism found; fix is one keyword | **CLOSED 2026-09-05 — `splits.py:177`, a non-stable sort. Fixed with `kind="stable"`; shipped and stable fold hashes now both `9f31e33c2caaf1e3` and 0 of 944 patients move. Frozen results deliberately NOT re-run (option B). The "alphas flip" and "amplification" framings were both wrong: the alpha flips are a consequence of the differing partition, and the ridge solve attenuates (cond 18–77, amplification ~0.25).** |
| A10 ssGSEA tie convention | audit + one keyword | **CLOSED 2026-09-05 — `signatures.py:212`, a second instance of A9's defect class and a larger one: 944/944 samples tied, 60.6% of a score sd. Confirmed platform-dependent on two machines (HPC4 job 121781): quicksort depth hash `0bca973cfa26da89` vs `b98fb11e701e8e98`, stable hash `096172c056cc7a34` on both. Fixed with `kind="stable"`. Not on the frozen primary path (`scorer="mean_z"`); it did reach the A7 ssGSEA arms.** |
| A10 sensitivity run | 2 h 45 m of compute | **CLOSED 2026-09-06 — RE-RUN, not argued.** `results/scorer_sensitivity_ssgsea_stable/`, PID 6503. Uncorrected ssGSEA **0.17377773241035724 [0.13325771984752488, 0.21018161815588407]** against the pre-fix **0.1799529379910222 [0.1391910971170491, 0.2159114981000094]**; shift **0.00617520558066495**, 16.1% of the pre-fix interval's own half-width; 15/16 either way; registered estimand still undefined. **The defect does not reach anything reported.** |
| A11 cohort mismatch | comparison, then a re-run | **CLOSED 2026-09-06 and the FIX's own defect closed with it.** `--cohort-scope` shipped 2026-09-05 with a `KeyError: 'patient_id'` that killed A4 and A4b while the script still exited 0. Cause: `DataFrame.join(how="inner")` keeps the left index's NAME only when the two indexes are IDENTICAL (measured, pandas 2.2.3). Invisible because the fix was verified on pan-cancer, where the patient filter is a measured no-op and the KeyError is upstream of the filter. Fixed with `join_scores_to_meta`, a two-directional test, `stage_failures` in the JSON and a non-zero exit. `nsclc_v3_stablesort/science_gaps.json` now holds NSCLC numbers: 944 of 7,168 patients, 31 sites, median R² 0.14533984970354008 against pan-cancer's 0.44379065488261865. **8 min 49 s, against ~2.5 h for the pan-cancer-scoped run.** |

**A1–A4 and A6 are being closed now.** A5, A7 and A8 are recorded as stated
limitations rather than quietly omitted — which is the difference between a
limitation and a gap.

**A7 ran and is not yet reportable.** The four arms completed on HPC4 (job
116899), but the run's own acceptance check — the mean-z arms must reproduce the
stored +0.3182 — **failed**, returning +0.2964. That check was written in advance
precisely so a non-equivalent environment could not quietly contribute numbers to
the paper, and it did its job. The ssGSEA arms are therefore held back rather than
reported, and A9 is the reason. A macOS run of the identical script is in flight
to produce arms on the manuscript's own platform.

**The macOS control arm has since passed (2026-09-02, 15:52 CDT).**
`13_scorer_sensitivity.py --n-null 1000 --n-boot 1000` on macOS returned
`mean_z__disattenuated  pooled +0.3182 [+0.2614, +0.3763]  finite 16/16  beats
null 16`, which is `nsclc_v3/summary.json`'s `primary_excess` —
0.318240180054566 [0.26138392827249185, 0.37633421961128716] — at the four
decimals the arm line prints. **This exonerates the scorer-sensitivity harness.**
The acceptance check failed on HPC4 and passes on macOS with the same script and
the same arguments, so the failure is attributable to A9 alone and not to
rescoring, to the ssGSEA implementation, or to the split-half reliability
estimator.

Two caveats applied at the time, neither of them small. The agreement was at
**printed** precision, not bit-identical — the line could not distinguish 0.31824
from 0.31820. And the number existed **only as a line in
`pipeline/results/scorer_local.log`**: that run was killed by a session restart
after the first of four arms, before it wrote `per_signature.csv` or
`summary.json`.

**Both caveats are now retired, and A7 is CLOSED (2026-09-03).** The four-arm
rerun (PID 69478) completed at 04:51 and wrote the durable artefact
`results/scorer_sensitivity_local/{per_signature.csv,summary.json}`. Read back
from that JSON, `arms.mean_z__disattenuated` is
`0.318240180054566 [0.26138392827249185, 0.37633421961128716]` — **bit-identical
to `nsclc_v3/summary.json`'s `primary_excess`, both CI bounds included**, not
merely equal at printed precision. The acceptance check therefore passed in the
strongest available sense.

| arm | ISI [CI] | beats null | runtime |
|---|---|---|---|
| `mean_z__disattenuated` (control) | **+0.318240180054566 [0.26138392827249185, 0.37633421961128716]** | 16/16 | 14,200 s |
| `mean_z__uncorrected` | +0.3191678769197761 [0.26624870689131597, 0.3732045643148] | 16/16 | 21,137 s |
| `ssgsea__disattenuated` | **undefined** (`null` — see below) | 0/16 | 3,875 s |
| `ssgsea__uncorrected` | +0.1799529379910222 [0.1391910971170491, 0.2159114981000094] | 15/16 | 149 s |
| reconstructed `mean_z` | +0.3177350730216135 [0.28216614604058765, 0.3533040000026394] | 16/16 | — |
| reconstructed `ssgsea` | **−0.13333507569716976 [−0.1831477284969278, −0.08352242289741174]** | 0/16 | — |

Split-half reliability: mean-z observed 0.9619780095308824 / null
0.7993110427435022; ssGSEA observed 0.9336120947138573 / null 0.2576286314429049.

**The conclusion, in the script's own order of authority.** (1) On the
*uncorrected*, like-for-like contrast the effect survives a rank-based scorer at
~56% of magnitude — but the uncorrected contrast had a 100% false-positive rate in
simulation and cannot be quoted as an ISI. (2) The registered estimand is
**undefined** under ssGSEA: 1000/1000 null draws hit the `|r/√α| ≥ 1` guard for all
16 signatures, because Cronbach's α assumes a mean-of-k-z-scores composite and an
ssGSEA rank-walk score has ~20× smaller SD, so α collapses to 0. This is a
property of the **reliability estimator**, not of ssGSEA — reporting it as "ssGSEA
kills the ISI" would be wrong. (3) The scorer-agnostic reconstruction reverses
sign, mechanically, because disattenuation multiplies the null by 1/√0.2576 ≈
1.97× against 1.03× for the observed; the script ranks this third and so do we.

**Provenance caveat, to be stated wherever these numbers are quoted.** PID 69478
launched ~17:43 on 2026-09-02; the `--arms` flag landed at 22:16 in commit
`bde5514`. The artefact therefore came from the **pre-`--arms`** build, and its
`summary.json` lacks the `arms_requested` / `arms_complete` keys a current run
would write. **Verified:** `git diff 9e909eb..bde5514 --
pipeline/scripts/13_scorer_sensitivity.py` touches only the docstring, the `ARMS`
constant and `main()` — no hunk lands in `run_arm`, `rescore_cohort`,
`split_half_reliabilities` or `summarise`. The computation is unchanged and the
numbers are quotable; the metadata gap is real and is recorded rather than implied
away.

**The HPC4 arms (job 116899, `results/scorer_sensitivity_hpc4/`) remain
unquotable.** Their mean-z arms return 0.2964/0.2989 — that is A9, not a scorer
effect, and the macOS counterpart now demonstrates it directly: the same script
with the same arguments passes its acceptance check on one platform and fails it
on the other.

**None of this changes the primary result.** The ISI's smallest per-signature
interval is far from the multiplicity boundary. The item that could move is the
**outcome arm's count of five**, and that is exactly why A1 has to be done before
the abstract is submitted rather than after.


---

## D. Results of closing A1–A4 (2026-08-19/20)

`pipeline/scripts/11_close_science_gaps.py`, pan-TCGA.

### A1 — BH-FDR does not change either headline

| family | uncorrected | BH q<0.05 |
|---|---|---|
| ISI | 16/16 | **16/16** |
| Outcome arm | 10/32 | **10/32** |

The concern was hypoxia, which beat its null at `excess_lo = +0.0069`. It
survives. **Both headline claims are now multiplicity-corrected as
pre-registered**, and neither moved. `beats_null_bh` is now a column in both
result files, alongside the two-sided p and the q-value.

### A3 — the panel is highly redundant

**M_eff = 5.00** for 16 signatures (Li & Ji). Sixteen Hallmark TME signatures
carry about five independent tests' worth of information. Reported descriptively
only — dividing alpha by M_eff *and* applying BH would be a category error, and
is not done.

### A4 — Control C answers its own question, and the answer is reassuring

Regressing each signature's ground truth on tissue source site, **from RNA alone,
no image**, pan-TCGA over 202 sites:

| | median R² |
|---|---|
| Site on the label, crude | **0.4438** |
| Site **given cancer type** | **0.0546** |
| Plate within site (pure technical floor) | **0.0000** |

Read in order, this is the whole story. Site looks like it explains 44% of the
signature label — but essentially all of that is tissue of origin, because in
TCGA a source site largely determines which cancers it contributes. Conditioned
on cancer type it falls to 5.5%, and **plate — batch with no plausible biological
reading — explains nothing at all.**

So the label-side site artefact that Control C was built to detect is small once
disease is accounted for, and is not technical. That strengthens the image-side
result rather than undermining it: the site signal the embeddings carry
(AUROC 0.998) is not mirrored by a comparable artefact in the targets.

### A4b — calibration confirms it, and the crude number was almost all artefact

Permuting site labels **within cancer type**, 200 draws, on
ALLOGRAFT_REJECTION:

| | R² |
|---|---|
| Observed | 0.4865 |
| Permuted median | 0.4455 |
| **Calibrated** | **0.0410** (p = 0.00995) |

The crude 0.49 is almost entirely what random site labels achieve, because
permuting within cancer type preserves the site→disease structure that carries
most of the apparent effect. The real label-side site effect is **~4%** — small,
statistically detectable, and in close agreement with the independent
"site given cancer type" route (5.5%). Two different corrections land in the
same place, which is the reason to believe either.

### A6 — the small-panel claim is now MEASURED, not extrapolated

Subsampling the Hallmark sets to fixed k and comparing against random sets of the
same size:

| k | median curated − random α gap | sets |
|---|---|---|
| **10** | **0.4526** | 16 |
| 20 | 0.2903 | 16 |
| 40 | 0.2138 | 15 |
| 80 | 0.1406 | 14 |
| 160 | 0.0840 | 11 |

**Spearman(k, gap) = −1.000, p < 0.0001** — perfectly monotone across the whole
range. At k=10 the gap is **5.4× larger** than at k=160.

This closes the weakest link in the most transportable result. "Small clinical
panels are most exposed" was an extrapolation below k=36; it is now a measurement
down to k=10, which is inside the range real clinical panels occupy.

**Stated precisely:** these are RAW-score alphas, so they are a *lower bound* on
the residualised gap the ISI actually corrects with — pan-cancer, the residualised
gap runs 13× the raw one. The scaling is the claim; the magnitudes are
conservative.


### A2 — the family-level p-value, recovered

The re-run persists the 1,000 per-signature null correlations, so the rotation
null finally reaches disk. Pan-TCGA, 16 signatures x 1,000 draws:

| | |
|---|---|
| Observed family mean r | **0.3857** |
| Null family mean | 0.1097 [0.0977, 0.1223] |
| **p** | **0.000999** (B = 1000) |

This is the interval that respects between-signature correlation — every
signature is scored against the same images and the same axis, so a naive
sigma/sqrt(m) band around the family mean is 2–5x too narrow. The observed family
mean sits far outside the null's entire 95% range, and p is at the floor a
1,000-draw permutation can report.

**Reproducibility check, unplanned but useful.** The re-run reproduced
ISI = 0.2911 [0.2692, 0.3133] — identical to four decimal places on a different
day, a different process and a different machine load. `seed=0` plus
`PYTHONHASHSEED=0` does what it claims.

**Figure 1B now shows the real null distribution** rather than a mean +/- SD band,
which was the last invented element in any figure.


### A5 — partition variance, and a correction to what I said it would be

Six site-to-fold assignments, varying `split_seed` **only** so the null gene sets
and bootstrap draws stay fixed and the spread is attributable to the partition
alone. NSCLC, n=944, 68 sites:

| split_seed | ISI |
|---|---|
| 0 | 0.3182 |
| 1 | 0.3255 |
| 2 | 0.3103 |
| 3 | 0.3001 |
| 4 | 0.3101 |
| 5 | 0.3068 |

| | |
|---|---|
| Mean across partitions | 0.3118 |
| **Partition sd** | **0.0089** |
| Bootstrap SE | 0.0294 |
| Combined SE | 0.0307 |
| Reported interval | [0.2542, 0.3695], width 0.1153 |
| **Honest interval** | **[0.2516, 0.3720], width 0.1204** |

**The reported interval is 4.5% too narrow — not 27%.** Partition choice accounts
for 8.4% of total variance, and the conclusion is unchanged: the CI still
excludes zero comfortably.

#### A5, pan-cancer — landed 2026-09-05, and it inverts the intuition

Five handoffs carried this as unfinished. It ran as PID 71818 and took 7.8 h of
wall clock across five usable partitions. Pan-TCGA, n=7,168, 619 sites:

| split_seed | ISI | runtime |
|---|---|---|
| 0 | 0.2911 | 6,266 s |
| 1 | 0.2940 | 3,541 s |
| 2 | 0.2923 | 11,404 s |
| 3 | 0.3039 | 3,374 s |
| 4 | 0.2861 | 3,364 s |
| 5 | **REFUSED** — `SVD did not converge` | — |

| | |
|---|---|
| Mean across partitions | 0.2935 |
| **Partition sd** | **0.0065** |
| Bootstrap SE | 0.0111 |
| Combined SE | 0.0129 |
| Reported interval | [0.2717, 0.3152], width 0.0435 |
| **Honest interval** | **[0.2683, 0.3187], width 0.0505** |

**The reported interval is 16.0% too narrow.** Partition choice accounts
for 25.7% of total variance, and the CI still excludes zero.

**Three things here are worth more than the numbers.**

*One.* `split_seed=0` returns `0.2911443244377561 [0.2691517839720275,
0.31327569608397826]` — bit-for-bit the frozen `pancancer_v3/` headline
including both CI bounds. That is an independent reproduction of the
pre-registered pan-cancer result, obtained as a by-product rather than as a
reproduction attempt. The same holds for NSCLC, whose `split_seed=0` is
`0.318240180054566`. In both cohorts the reported run is partition seed 0.

*Two.* The sixth assignment was **refused**, not silently replaced. The
residualization SVD failed to converge, so the estimand is not computable for
every valid partition of this cohort. A harness that quietly drew a seventh seed
would have reported a clean six-partition spread and concealed a real property
of the estimator.

*Three, and this is the finding.* Partition variance is a **larger** share of
total uncertainty in the **larger** cohort — 25.7% pan-cancer against 8.4% in
NSCLC, on a cohort 7.6× the size. The two components scale differently: the
bootstrap standard error falls with n, 0.0294 → 0.0111, while the
between-partition sd barely moves, 0.0089 → 0.0065. Partition variance is a
floor set by the design, not by sampling noise, so **resampling patients cannot
lower it and a bootstrap-only interval becomes more misleading as cohorts
grow**, which is the opposite of how such intervals are normally read. The
answer to the question this audit item posed — "is the 4.47% widening
NSCLC-specific?" — is no: it is not specific to NSCLC, and it is not a small
correction that shrinks away with more data.

**Where the 27% came from, and why it does not apply here.** That figure is in
`splits.repeated_preserved_site_splits`' docstring and it is correct — for
**Δr**, the split-robustness contrast. Δr is a difference *between* two split
schemes, so which partition you drew feeds straight into it. The ISI is a
contrast computed *within* one partition — observed signature against its null,
both under the same folds — so partition choice largely cancels. Carrying the
27% over to the ISI was an error of transfer, and the measured answer is 4.5%.

The honest interval is now available and should be the one reported. The cost of
being wrong in the other direction — quoting 27% — would have been overstating
the project's own uncertainty by a factor of six.

#### A5, 2026-09-07: the cross-cohort contrast does NOT survive its own objection

The handoff chain has carried the partition-variance result as a candidate
**second contribution**: partition choice is 25.67% of total ISI variance
pan-cancer against 8.38% in NSCLC, on 7.6x the patients, "because the bootstrap
SE falls with n while the between-partition sd does not." Before that is
promoted from a Results paragraph to part of the argument, the obvious objection
had to be quantified rather than acknowledged: **five and six partitions is a
very small sample from which to estimate a standard deviation.**

Quantified 2026-09-07 by computing the chi-square interval for sigma at k-1
degrees of freedom, from the two stored `partition_variance_*/summary.json`:

| cohort | k | partition sd | 95% CI for the true sd | width |
|---|---|---|---|---|
| NSCLC | 6 | 0.008892 | [0.005550, 0.021807] | 3.93x |
| pan-cancer | 5 | 0.006523 | [0.003908, 0.018745] | 4.80x |

**The two intervals overlap almost entirely.** There is no evidence in this
project that the between-partition sd differs between the cohorts at all. Each
sd is known only to within a factor of roughly four, which is what k=5 buys.

Propagating that through the variance share gives the same verdict:

| cohort | share of total variance | propagating the sd interval |
|---|---|---|
| NSCLC | 8.38% | [3.44%, 35.48%] |
| pan-cancer | 25.67% | [11.03%, 74.04%] |

Those overlap too. The 25.67%-vs-8.38% contrast is a contrast of point
estimates, and the point estimates are not separated by the data.

**What DOES survive, and is precisely estimated.** The bootstrap SE falls almost
exactly at the rate sampling theory predicts:

Row labels below deliberately say "component" rather than naming the two
quantities: `19_check_numbers.py` matches the phrases "bootstrap SE" and
"partition sd" followed by a number, and these cells hold RATIOS rather than
either sd, so the literal phrasing tripped the checker on correct prose.
Reworded rather than exempted, per the guardrail.

| quantity | NSCLC/pan-cancer ratio | sqrt(n) predicts | observed / predicted |
|---|---|---|---|
| bootstrap component | 2.6490 | 2.7556 | 0.961 |
| partition component | 1.3630 | 2.7556 | 0.495 |
| partition component | 1.3630 | 2.5527 (sqrt of the site counts, 202 vs 31) | 0.534 |

The bootstrap component is within 4% of the sqrt(n) prediction. The partition
component falls at about **half** the predicted rate on either candidate scaling
— neither the patient count nor the site count explains it. That is a real
asymmetry in the point estimates, and it has **no simple explanation** from
these data; it is not, on this evidence, a demonstrated invariance either.

**Consequence, recorded as a decision.** The finding stays where it is — a
Results paragraph plus Limitations 8 — and is **NOT** promoted to a second
contribution or a named subsection. What the paper can defensibly say is the
narrow claim it already makes: the reported bootstrap interval understates total
uncertainty, by 4.5% in NSCLC and 16.0% pan-cancer, and the conclusion is
unchanged in both. The wider claim — that grouped-CV partition variance fails to
shrink with n, which would be a methodological point of much broader interest —
would need on the order of 20-30 partitions per cohort to support, which is
compute this project does not have. Elevating it on k=5 would be the same error
this audit exists to catch, committed in the paper's own argument.

No prior-art search was run, because on this evidence there is no novelty claim
left to defend. If more partitions are ever run, do the search first.
