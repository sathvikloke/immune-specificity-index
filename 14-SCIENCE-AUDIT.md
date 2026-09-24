# Science audit — what is actually left

> **SUPPLEMENTARY TABLES WERE RENUMBERED ON 2026-09-24** (B-20, session 59), in the
> order the manuscript first cites them. Every table number and table file name
> written in this file before that date is an OLD number. Old → new: S3 → S1,
> S4 → S2, S2 → S3, S6 → S4, S10 → S5, S8 → S6, S5 → S7, S1 → S8, S7 → S9,
> S9 → S10. The tables' contents did not change (each file was compared byte for
> byte). Bisection STAGES (`S0`–`S9`, e.g. `S3_fold_alpha`) are a different
> namespace and were never renumbered.

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

*[Corrected 2026-09-21: the macOS row's "macOS 15" was never verified and is
wrong. By the laptop's own install record it ran macOS 26.5.2 (Darwin 25.5.0)
from 2026-07-19 until macOS 27.0 on 2026-09-20; see the session 54 block. The
table is kept as written.]*

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
reruns this and gets a different fourth digit *[2026-09-23: the gap, 0.3182
against 0.2964, is in the SECOND decimal — the correction session 50 made
elsewhere (F4.10) never reached this record]*. Two honest responses, not mutually
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
wrong, and reading `18_bisect_platform.py:338-344` is what corrected it.** [2026-09-22:
the lines meant are the per-signature loop's `_select_alpha_per_fold` call and the
`cross_val_predict_multi(..., fixed_alpha=fold_alpha)` that follows it. They were
338-344 in the file as read (`8981c25`), and the commit that wrote this sentence
(`0596a81`) also added lines above them, so the pointer was already stale in the
commit that created it.] There
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
reporting both. (That tag lives in the **private** working repository; this file
is published inside the snapshot, which carries no tags. The link the sentence
describes is real, but a reader of the snapshot cannot check it there.) Measured re-run cost, from the frozen `summary.json` files:
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
`signatures.py:206`'s own comment [in `signatures._ssgsea_sample_tables`; line
239 by 2026-09-22] reads "Ties are common in log-expression
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

[2026-09-22: the call sites in this table and the paragraphs below are line
numbers of 2026-09-05, and two had moved by 2026-09-22 — `models.py:339` into
`impute_median` (the sort it meant is the last line of
`models.site_prediction_control`) and `experiment.py:249` into a docstring (the
median is recomputed in `AuditResult.headline`). The live copies —
`ENVIRONMENT.md`'s sort census and `26_sort_audit.py`'s case labels — now name
the function: `stats.bh_fdr`, `barcodes.site_summary`,
`models.site_prediction_control`, `outcome.concordance_index`. This record keeps
its coordinates.]

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

[2026-09-22: the "C FIXED" above is what the audit printed on 2026-09-06, and
it printed it for the wrong reason. `verdict()` already returned DEAD CODE before
looking at the pin; case C reached FIXED only because its call detector matched
the bare name `site_summary`, so `ENVIRONMENT.md`'s sentence saying the function
is never called counted as a call — the self-defeating match fixed on
2026-09-16. Since then the audit reports C as DEAD CODE, which is the correct
verdict; the pin is still in place and still registered through
`shipped_pins_stable`. Verified 2026-09-22 against `ef4de1c`'s own source.]

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

### A12, NEW 2026-09-16: the RANDOM-patient split is platform-dependent too — a library ordering the A10 audit could not see

`splits.random_patient_split` shuffles patient order with a seeded generator and
hands the result to scikit-learn's `GroupKFold` (1.6.1), whose
`_iter_test_indices` orders groups by `np.argsort(n_samples_per_group)[::-1]` —
the default, unstable sort. After `collapse_to_patient` every group has size 1,
so every key ties and the order is the kernel's tie rule. Measured on an
all-ones int64 array: macOS arm64 returns hash `9085a1f89f19e541` (n=7,168;
`cb54ec6dc2893b62` at n=944), the HPC4 login node `db52d4be6eabc0f2`
(`84f9a4114100edfa`). Neither is the identity permutation.

Consequence, measured on the corrected-ordering NSCLC run (macOS
`results/nsclc_v3_stablesort/` against HPC4 job 129151's output):

| scheme | fold agreement | max abs prediction difference |
|---|---|---|
| preserved-site (the ISI) | 1.000 | 1.260e-14 |
| random-patient (embedding, all four baselines) | 0.193 | up to 1.340 |

0.193 is chance (0.2). Per-signature `r_random_patient` and `delta_r` differ by
up to 0.034 (0/16 bitwise equal); `secondary_delta_mae` is
0.005202722334422926 [0.0022131927198345775, 0.008001492984137] on macOS
against 0.00481 [0.00199, 0.00748] on HPC4 (full precision in the HPC4
checkout's `results/nsclc_v3_stablesort/summary.json`, which is not stored in
this repository — hence the rounding) *[2026-09-24: it IS stored now, as
`results/session48_diagnostics/job129151_nsclc_stablesort_hpc4/summary.json`
(copied by A6), hashed in PROVENANCE since 2026-09-24, and the checker reads
the full-precision interval from it]*. The ISI agrees to one ULP, as job 129151 reported. So the validation gate
is valid **for the preserved-site path only**, and every secondary quantity in
the paper is a macOS number. Manuscript Limitations 8 now says so.

Why A10's sweep missed it: `26_sort_audit.py` sweeps orderings in this
repository's source. This one is in a dependency. Not fixed: pinning it changes
the random split on macOS as well, so the frozen secondary numbers would no
longer reproduce on either machine — the A9 trade-off again, and it needs the
same explicit decision (Option B) before code changes. *[2026-09-24: decided by
the author — option (b) of the memo below — and done: see "Session 60: A12
pinned" at the end of this file.]*

**Now visible to the audit (ledger C9, 2026-09-17).** `26_sort_audit.py` has a
case F for this library ordering. *[2026-09-24: it now reports `FIXED`, and its
self-test's first F check asserts the shipped source reads as pinned.]* It reports the verdict `DOCUMENTED, NOT PINNED`
rather than `DANGEROUS`, so the audit's exit status is unchanged while the
decision is open, and it changes no split. `--self-test` fires in three
directions:
- on the shipped, unpinned source;
- on an injected pinned source, which must read as pinned;
- on groups of all-different sizes, which must NOT change the folds.

`--library-order OUT` records case F alone, per machine:
`pipeline/results/groupkfold_order_darwin_arm64.json` and
`groupkfold_order_linux_x86_64_login.json`, plus a compute-node record from HPC4
job 130190. The hashes reproduce the ones above, and the checker now asserts
both halves of the claim. The stable GroupKFold fold assignment of a 944-patient
permutation agrees across machines (`caf92a3f83ddb5b9`), and the library one
does not (`8fdee48c0ac8a457` against `8f4427a6cc4567e6`). If the author pins it,
the pinned helper must be named `stable_group_kfold` in `splits.py` for case F
to see it.

#### A12 decision memo (ledger D1), measured 2026-09-16/17 — for the author

The question: should the random-patient split's tie order be pinned? Measured
by re-running the shipped NSCLC configuration with scikit-learn's
`GroupKFold._iter_test_indices` patched to `np.argsort(..., kind="stable")` in
that process only (no source change), on both machines. Records:
`pipeline/results/session49_diagnostics/d1_pinned_groupkfold/`.

| comparison | ISI | per-signature Δr, max / mean abs diff | Δ-MAE |
|---|---|---|---|
| macOS unpinned control vs `nsclc_v3_stablesort/` | identical | 0 / 0 (16/16 bitwise) | identical |
| macOS pinned vs macOS unpinned | identical | 0.043 / 0.016 | 0.00428 [0.00138, 0.00703] vs 0.00520 [0.00221, 0.00800] |
| HPC4 unpinned vs macOS unpinned (A12 above) | one ULP | 0.034 / 0.014 | 0.00481 vs 0.00520 |
| **HPC4 pinned vs macOS pinned** | one ULP | **8.9e-16** (every secondary within 2.8e-14) | agree to 9e-17 |

So: (1) the tie order moves every random-split quantity by about as much as
changing platforms does, and by about as much as any other valid random
split would; mean Δr is 0.0205 unpinned, 0.0168 pinned. (2) Pinned, the two
machines agree on every secondary quantity to round-off, so pinning WOULD make
the secondary numbers portable. (3) No conclusion changes under any of the
four: the preserved-site ISI is untouched, Δ-MAE's interval excludes zero in
all of them, and the per-signature Δr values stay small.

Options, unchanged in kind from A9: **(a)** leave the split as is and keep
the documented platform dependence (the current state); **(b)** pin it in
`splits.random_patient_split` and report the secondary numbers both ways,
re-running NSCLC (about 4 minutes on an idle Mac) and pan-cancer (HPC4); **(c)**
pin it behind a configuration flag for new runs only. Recommendation: (b),
for the same reason A9 took it: the pinned numbers are measured to be
portable and the frozen ones stay reported beside them. **BLOCKED-AUTHOR on
the choice.** *[2026-09-24: the author chose (b); done — "Session 60: A12
pinned" at the end of this file.]*

#### A5 correction — every pan-cancer "refusal" was a LAPACK failure in a secondary analysis

`12_partition_variance.py` caught `(RuntimeError, ValueError)` and called the
result a refusal "by the degeneracy guard". `numpy.linalg.LinAlgError`
subclasses `ValueError`. HPC4 job 129294 captured the failing inputs with a full
traceback:

| seed (stable sort, HPC4) | failing call | matrix | gesdd retry | gesvd | numpy svd |
|---|---|---|---|---|---|
| 6, 11, 12 | covariate baseline — `models.cross_val_predict` → `RidgeCV` → `scipy.linalg.svd` | 5734 × 650, rank 590–595, 29–34 all-zero dummy columns | fails (copy and Fortran order) | converges | converges |
| 15 | variance decomposition — statsmodels OLS → `pinv` → `numpy.linalg.svd` | site dummies again | — | — | — |

The covariate baseline runs under the RANDOM-patient split, so which partition
fails depends on A12's platform-specific split as well as the LAPACK build. On
macOS, repeating only that regression for split_seed 0–23 fails for exactly one
of 120 seed-folds — **split_seed=5, fold 4, the frozen run's refused partition**
— and gesvd converges on it. Seeds 0–4, which completed, all pass.

Neither the baseline nor the decomposition feeds the ISI. The manuscript,
15-CHECKLIST.md, the inventory and the S8 caption said "the estimand is not
computable for every valid partition"; that is withdrawn. `--isi-only` switches
both analyses off for the partition sweep, and every partition not computed is
now recorded with its KIND (`degeneracy_guard`, `linalg_error`, other).

### A13, NEW 2026-09-16: nothing in the repository writes `data/interim/`

`reproduce.sh`'s full path runs `01_fetch_data.py --all`, then the cohort
scripts, which read six files from `data/interim/` (`X.npy`, `meta.parquet`,
`expression_hugo.parquet`, `cohort_nsclc.parquet`, `X_nsclc.npy`,
`expr_nsclc.parquet`). `01_fetch_data.py` downloads the embeddings and PRINTS
instructions for the other four sources; it writes no interim file. A search of
`scripts/` and `src/` finds no `to_parquet` or `np.save` for any of the six, and
no commit in the history adds or removes such a writer. The files are dated
2026-08-17 10:38–10:41, before the first commit that mentions them. `data.py`
carries the loaders (`load_embeddings`, `assemble_cohort`, `load_expression`,
`load_purity`), so the logic exists; the step that runs it and writes the
result does not. A reader with the raw sources cannot start the full path. To
close: a builder script whose outputs hash-match the six recorded prefixes
(`86305b9f7c97950f`, `6048e461aab97f87`, `2a36f4a3d7f36a43`, and
`0b396a1e2780b0b0`, `7a1408eac2f64f7c`, `8532e31cf0176acd` for the NSCLC three),
run where 3.2 GiB dense plus parsing fits — HPC4, not this laptop.

**Progress, 2026-09-16 (session 49).** `scripts/00_build_interim.py` now
exists. It was reconstructed from the files themselves: `meta.parquet` is
`data.load_embeddings(..., barcode_col="sample")`, whose primary-tumour
filter also removes the parquet's repeated slide rows; `X.npy` is the last
of 14 encoder layers; the NSCLC cohort is the LUAD and LUSC slides collapsed
to patients, joined to ancestry and ABSOLUTE purity, kept if the patient has
any TOIL sample, and scored on the re-read expression frame. On macOS it
reproduces `meta.parquet`, `X.npy`, `X_nsclc.npy`, `expr_nsclc.parquet` and
`cohort_nsclc.parquet` byte for byte (pyarrow 24.0.0); the cohort file
matched only after clearing the `join_coverage` attribute that pandas
otherwise writes into the file. Writing it surfaced A16.

**CLOSED, 2026-09-16.** On HPC4 (job 129998, raw inputs downloaded there from
their public sources and hash-identical to the laptop's; the embeddings
pinned to dataset revision `073115403c2fc5134ee8d1332c603edba591dddb`) the
builder reproduced `expression_hugo.parquet` byte for byte, which the laptop
cannot hold, together with `meta.parquet`, `X.npy`, `X_nsclc.npy` and
`expr_nsclc.parquet`. Linux's `cohort_nsclc.parquet` differs from the
macOS-built file only in its 16 `SIG_` columns, by at most 1.1e-15: the
mean-z reductions round differently on x86_64. All six files are therefore
reproduced from the raw sources; `--verify` reports that residual rather
than calling it a content difference, and `reproduce.sh` builds the files
when they are absent.

### A14, NEW 2026-09-16: the shared-mask fast path is bitwise identical on real data

Why it exists: 19 pan-cancer patients have no cancer type, so their
within-type axis is NaN and their residualised scores are NaN in every null
column. Every preserved-site training fold contains some of them — no clean
fold in any of 24 partitions (job 129663) — so every null target was refitted
one column at a time. On seed 7 (job 129664) that per-column refit was 708.2 of
850.37 s for the first signature's null.

`models.SHARED_MASK_FAST_PATH` shares the centring and `X'X` across columns
with the same missingness mask and keeps every per-column step of
scikit-learn 1.6.1's cholesky Ridge path. The claim is equality, not closeness,
so it was tested for equality:

| check | where | result |
|---|---|---|
| synthetic, every column forced through the fallback, two masks, per-fold alphas | macOS, `test_shared_mask_fast_path_is_bitwise_equal_to_per_column_ridge` | `np.array_equal`, with a one-ULP can-fail |
| real pan-cancer folds, split_seed=7, the first two signatures' null calls (each 7,168 × 1,001) | HPC4 Linux, job 129737 | bitwise equal, max abs diff 0.0, same NaN pattern; `ALL_BITWISE True` |

Speed on that job: reference 561.3 s and 550.83 s, fast path 290.87 s and
279.76 s per null call, about 1.9×. Record:
`pipeline/results/session48_diagnostics/job129737_verify_fast_path/`
(private repository; excluded from the deposit because the logs carry cluster
paths).

Two companion gates for the partition sweep, both run 2026-09-16:

- **`--isi-only` changes nothing the ISI reads.** Job 129666 ran NSCLC
  split_seed=0 with and without it: ISI 0.2922124012879652 and both interval
  ends identical to each other and to job 129201's row (`GATE_BITWISE True`).
  So `--isi-only` shards pool with default-mode shards.
- **Seed 4 is not dearer than seed 7.** Job 129665 profiled seed 4's first
  two signatures: 840.72 s and 571.22 s per signature's null, against seed 7's
  850.37 s and 854.32 s. The seven timeouts in array 129207 are therefore a
  matter of placement and node load at a cost that sits near the 4-hour cap
  (sixteen signatures at 570–850 s), not of which seeds they were. Seed 21, the
  last default-mode shard, also timed out with nothing written.

**NSCLC never takes the fast path**, so no NSCLC run is evidence about it.
Measured on macOS by wrapping the two functions during a full NSCLC ISI run
(split_seed=0, 1,000 null sets): 16 null calls of shape 944 × 1,001, zero
non-usable target columns in any training fold, zero fast-path invocations.
NSCLC has no patient without a cancer type and no missing value on its global
axis (0 of 944). The same run returned 0.2922124012879651
[0.23552139004939532, 0.3527272668764091], identical to
`results/nsclc_v3_stablesort/`.

**B2, 2026-09-17: the partition sweep reproduces the primary script bit for
bit on the real cohort.** Array 129874's seed 0 (fast path and `--isi-only`,
`12_partition_variance.py`) returned 0.29684911595367525 [0.274920163364846,
0.31893289976112255], exactly `09_run_pancancer.py`'s job 129159; its seed 19
equals the default-mode, reference-path value. The seven re-run shards took
90–112 min each, against 123–236 min for the default-mode shards that
finished and more than 240 min for those that did not.

**Memory (ledger F1.5), from `sacct` MaxRSS:** every pan-cancer partition job,
in any mode and with or without the fast path, peaked at 22.3–23.4 GiB; the
full pan-cancer run (job 129159) at 23.4 GiB; the NSCLC run (job 129151) at
3.2 GiB. The laptop has 24 GiB in total, so a pan-cancer run cannot fit there
whatever its free-memory floor on the day; that, not the sampled floors, is
why those runs moved to HPC4.

The pan-cancer partitions that timed out are re-running with the fast path and
`--isi-only` (array 129874), with seed 19 as a control on both changes
together (its default-mode, reference-path ISI is known).

### A15, NEW 2026-09-16: pandas' default CSV float parser is not exact, and not the same on every platform

`pd.read_csv` without `float_precision="round_trip"` returns some doubles one
or more ULPs away from the decimal string in the file. Measured on macOS: 209
of 400 random doubles written with `repr`, and 16 of the 44 cells of the
committed S8 table, including the NSCLC headline's upper bound (the table said
0.3763342196112871 against the frozen 0.37633421961128716). Linux reads some
cells differently again: rendering the supplement on HPC4 (job 129926) gave
27 differing cells in S3 and one in S8. That is why one of the three tests
`reproduce.sh`'s pre-flight suite failed there did so (job 129695); the other
two are the figure checks, whose PDFs and PNGs depend on the FreeType
matplotlib links (2.13.3 here, 2.6.1 in HPC4's wheel). All five SVGs match on
both, so the SVGs are now compared everywhere and the PDFs and PNGs only where
the FreeType matches.

With exact parsing, every supplementary table is byte-identical on both
machines except S3's NSCLC `p_two_sided` and `q_bh`. The builder derives those
from normal tail probabilities as small as 1e-41, and scipy returns them with
platform-dependent last digits (22 cells, all beyond the 15th significant
digit). `25_build_supplementary.py` now parses exactly and gives those two
derived columns to twelve significant digits. Ten committed files changed at
the ULP level on re-rendering; every transcribed cell now equals its source
exactly, which a new test checks independently of the generator.

The class is wider than the supplement. The same default parser reads the
purity table in `data.load_purity` (read only by the covariate baselines and
the decomposition, never by the ISI), the stored predictions in
`11_close_science_gaps.py` and the stored per-signature tables in
`10_make_figures.py`. None of those feeds a headline number, and none has been
changed.

### A16, NEW 2026-09-16: the expression profile is not always the tumour

Found while writing the interim builder (A13). Neither cohort restricts
expression to primary tumour (TCGA sample code 01), although the slides are
all primary-tumour diagnostic slides:

- **Pan-cancer.** `09_run_pancancer.py` (and `06_run_ancestry.py`,
  `11_close_science_gaps.py`, `12_partition_variance.py`) truncates every
  TOIL sample id to the patient and takes the mean. Of the 7,168 patients,
  6,550 have exactly one TOIL sample and it is the tumour. 578 are averages of
  several samples (520 tumour plus adjacent normal, 27 tumour plus recurrence,
  20 tumour plus metastasis, 11 others), and 40 have no tumour sample at all
  (28 adjacent normal only, 8 recurrence only, 4 metastasis only). The most
  affected types are THCA (67), KIRC (62), PRAD (50), LUAD and LIHC (48 each).
- **NSCLC.** `expr_nsclc.parquet` holds each patient's FIRST TOIL sample in
  TOIL's column order. For 851 of 944 patients that is the only sample and it
  is the tumour. Of the 90 with a tumour and a normal sample, 48 received the
  normal profile. Two with a tumour and a recurrence received the recurrence,
  and one patient has only a normal sample. So 51 of 944 NSCLC signature
  scores describe tissue other than the slide's tumour.

The values are TOIL's own scale, log2(TPM + 0.001) (floor −9.9658), not
log2(TPM + 1) as `05_run_nsclc.py`'s docstring says.

**Measured for NSCLC, 2026-09-17** (HPC4 job 130002, the unchanged
`05_run_nsclc.py` on `--tumour-only` inputs; 943 patients, the normal-only
patient dropped): **ISI 0.3448 [0.2731, 0.4079]** against the
corrected-ordering 0.2922 [0.2355, 0.3527] on the same machine, a shift of
+0.0526 (about seven 24-partition standard deviations, and inside neither
interval's complement: the intervals overlap). 16/16 signatures still beat
their null, and 14 of 16 per-signature excesses rise (the largest rise is
TGF-β, +0.131; angiogenesis falls by 0.041). Adjacent-normal profiles were
diluting the image–expression association. Δ-MAE barely moves (0.00484
against 0.00481). The outcome arm gains one signature (hypoxia,
axis-residualised, excess 0.067 [0.003, 0.132]): 1/32 against 0/32.
`results/nsclc_v3_tumouronly/`.

**Pan-cancer, first attempt failed, 2026-09-16** (job 130003): the unchanged
`09_run_pancancer.py` died after 103 s in the covariate baseline's gesdd SVD,
the D8 exposure, now on a primary-configuration run: the tumour-only cohort
(7,128 patients) gives a different random-patient split. The index is being
computed through `12_partition_variance.py --isi-only` at split_seed 0 (job
130155), which reproduces `09_run_pancancer.py` bit for bit on the original
inputs (ledger B2, below).

**Measured for pan-cancer, 2026-09-17** (job 130155, 6,293 s, fast path;
`pipeline/results/pancancer_tumouronly/`). Over 7,128 patients and 619 sites,
**ISI 0.2937 [0.2701, 0.3149]** against the corrected-ordering 0.2968
[0.2749, 0.3189]. The shift is −0.0031, less than half the 24-partition
standard deviation (0.0072). In NSCLC the correction matters, because an
adjacent-normal profile replaced the tumour for 48 patients. Pan-cancer it
does not, because mixed profiles were averaged with the tumour rather than
substituted for it. Both values are now in Results ("Sensitivity to other
analysis choices"). The frozen values stay the reported ones. Whether the
tumour-only values replace them or sit beside them is **BLOCKED-AUTHOR**, as
for A9.

Otherwise the frozen results are unchanged. `00_build_interim.py --tumour-only` builds corrected inputs
(sample code 01 only; patients without one are dropped), and the sensitivity
runs use them. Whether the correction replaces or sits beside the frozen
numbers is the author's decision, as it was for A9.

### D5, 2026-09-16: which OpenBLAS the results ran on

Every statement of "the" BLAS in this project (this file's A9 sections,
`ENVIRONMENT.md`, `20_check_versions.py`, the Software paragraph) quoted numpy's
build record, OpenBLAS 0.3.21. That is what numpy was compiled against, not
what runs. `threadpoolctl` finds one OpenBLAS in the macOS process, conda's
`libopenblasp-r0.3.29.dylib`, shared by numpy and scipy (scipy's own build
record says 0.3.29), and conda's install history shows that library added on
2025-08-21 with no BLAS, numpy or scipy change since, so 0.3.29 is what
produced the frozen results. On HPC4 the two wheels each load their own
library in the same process: 0.3.27 for numpy, 0.3.28 for scipy (job
129294). The earlier sections are left as written, since they record what
was printed; the exporter, the version check, `ENVIRONMENT.md` and the
manuscript now report the loaded library beside the build records. No
conclusion changes: A9 was a sort, not the BLAS.

### E1, NEW 2026-09-17: the index depends on removing the dominant axis WITHIN cancer type

The pre-registration fixed `axis_within_type=True` "in advance and not to be
changed", so this is an unregistered sensitivity analysis, reported because it
changes how the index should be read. NSCLC, macOS,
`scripts/29_nsclc_sensitivity.py`, ISI-only configuration whose reference
variant reproduces the corrected-ordering primary:

| axis | ISI | signatures above null | median r, curated (refit) | median r, random-set null | median null α | median uncorrected excess |
|---|---|---|---|---|---|---|
| within type (registered) | 0.2922 [0.2355, 0.3527] | 16/16 | 0.373 | 0.079 | 0.835 | 0.300 |
| global, both types together | −0.0096 [−0.0630, 0.0466] | 1/16 | 0.351 | 0.361 | 0.944 | −0.004 |

The curated signatures are about as predictable from the image either way. What
changes is the null: once the axis is removed only globally, the image predicts
a random gene set almost exactly as well as a curated signature, and the excess
vanishes before any reliability correction. The pipeline flags that null
itself: under the global axis, 10 of the 16 signatures have a null whose sd
across draws is below 0.02 (`degenerate` in
`results/nsclc_sensitivity/global_axis/immune_excess.csv`; 0 of 16 in the
reference). Their null means are r = 0.34–0.37 with sd 0.018–0.020, the
near-degenerate regime `globalaxis.py` describes, where a small excess is a
hair trigger. Which shared component the random sets are carrying was not
measured. So the positive index is a statement
about expression variation WITHIN a cancer type: the component every gene set
shares is removed there, and what remains is specific to the curated sets. It is
not a statement that survives an axis defined across types.

**Pan-cancer (HPC4 job 130189, 2026-09-17, 1,217 s, fast path;
`pipeline/results/pancancer_sensitivity/`).** The registered row is the
corrected-ordering primary, which the same ISI-only configuration reproduces bit
for bit (job 130185_1):

| axis | ISI | signatures above null | median r, curated (refit) | median r, random-set null | median null sd | median null α | median uncorrected excess |
|---|---|---|---|---|---|---|---|
| within type (registered) | 0.2968 [0.2749, 0.3189] | 16/16 | 0.391 | 0.100 | 0.025 | 0.801 | 0.306 |
| global, all 31 types together | 0.3512 [0.2692, 0.4353] | 15/16 | 0.725 | 0.452 | 0.069 | 0.805 | 0.424 |

The two cohorts move in OPPOSITE directions. In pan-cancer, removing only the
global first component (10.6% of expression variance) leaves much more
image-predictable variation in every gene set. Random sets rise from r = 0.100
to 0.452, but the curated signatures rise further, from 0.391 to 0.725. The
excess grows by 0.054, its interval widens from 0.044 to 0.166, one signature
drops below its null, and no null is degenerate. This fits the handoff's
reading that global PC1 is tissue of origin, so an axis removed only globally
leaves between-type structure in the scores. That reading was not measured
here. What was measured: the sign and size of the index depend on where the
axis is removed, in both cohorts, and in different directions. The registered
within-type choice is therefore part of the estimand's definition, not a
neutral default. Unregistered; reported as a sensitivity analysis.

### D3, D8, 2026-09-17: the two LAPACK exposure sites, census of 24 partitions

`pipeline/cluster/d3_d8_linalg_census.py` runs the audit's own baseline and
decomposition code for each split seed, with the index and nulls skipped, and
records every LinAlgError by site. HPC4 arrays 130202 and 130211;
`pipeline/results/linalg_census/`.

| cohort | seeds | covariate baseline fails | decomposition fails |
|---|---|---|---|
| NSCLC | 0–23 | none | none |
| pan-cancer | 0–23 | 6 (fold 2), 11 (fold 2), 12 (fold 1), 20 (fold 2) | 15 (MYC targets, rung `Mfull_+image`) |

- This reproduces session 48's five failing seeds exactly, on the same
  platform, and finds no others.
- **D8, the primary partition:** seed 0 passes both sites in both cohorts on
  Linux. The frozen macOS primaries also completed both analyses at seed 0.
  The primary runs were not near the failure in the partition they used.
- The macOS pan-cancer covariate census (session 48) failed only at seed 5
  fold 4; the random split differs across platforms (A12), so the failing
  seeds differ too.
- **The macOS half of D3** (the decomposition census) was not run. A local
  smoke test of two NSCLC seeds was stopped after 73 minutes at load 70–100
  from a foreign process.
- **D3 decision.** Document the exposure; do not change the code this
  session. A failed decomposition removes one secondary table, never the
  index. A gesvd or numpy fallback would change only the failing calls, but it
  is a change to a reported analysis's code path.
  **BLOCKED-AUTHOR:** add the fallback (and check bitwise on seed 0 first), or
  keep the documented exposure. *[2026-09-24, decided by the author: the
  decomposition's OLS now retries a failed SVD fit with a pivoted-QR solver
  (scipy `gelsy`) and warns. Only failing calls take that path, so seed 0 and
  every stored run are unchanged by construction; a test proves the refit
  gives the default fit's R² and adjusted R² to 1e-10. The covariate
  baseline's failures (sklearn `RidgeCV`) are left documented: they occur only
  in partition sweeps, never at the primary partition (D8).]* `09_run_pancancer.py` on the tumour-only inputs
  (job 130003) is the one primary-configuration run the exposure has actually
  stopped.

### F6.3, 2026-09-17: the eleventh axis, LIBRARY-ORDERING, swept

These are tie-dependent orderings inside dependencies, the class A12 opened.
Every sorting or grouping call in `src/aacr27/` was enumerated (grep for
GroupKFold, StratifiedGroupKFold, np.unique, get_dummies, groupby,
pivot_table, value_counts, sort_values, argsort, rank), and each member was
checked against the library source or by a hash on this machine.

| member | call sites | verdict |
|---|---|---|
| `GroupKFold._iter_test_indices`, `np.argsort(n_samples_per_group)[::-1]` | `splits.random_patient_split`; `models.site_prediction_control`; `models.site_prediction_control_oof_combat` | *[2026-09-24: no longer called — all three use `splits.stable_group_kfold`.]* **platform-dependent (A12).** Groups are all size 1 after collapsing to patients. Moves the random-patient secondaries and the site-control AUROC: NSCLC median 0.9916 (macOS) against 0.9894 (Linux, E16 at the full cohort); pan-cancer 0.99847747 against 0.99848372. The index uses none of them. Case F of `26_sort_audit.py` covers the first call site; the other two share its cause. |
| `StratifiedGroupKFold` | `splits.random_patient_split(stratify=...)` | safe: sklearn 1.6.1 sorts with `kind="mergesort"`, and no caller passes `stratify` |
| `np.unique` | `stats` (bootstrap groups) | safe: the values are distinct, so the order is defined |
| `pd.get_dummies`, `groupby(sort=True)`, `pivot_table` | decomposition, baselines, data, experiment | safe: sorted by distinct keys |
| `Series.rank(method="first")` | `signatures` (expression-matched null bins) | safe: pandas breaks ties by order of appearance (hash of an all-ties rank equals 1..n) |
| `DataFrame.sort_values` on a column with ties (default quicksort) | `site_control.csv` (AUROC, many ties at 1.0), `decomposition.csv`, `label_site_variance.csv`, `immune_excess.csv` | **row order only.** Measured on macOS: an all-ties sort is not the identity. The output CSVs' row order can differ across platforms while every value is identical. Every consumer (checker, figures) takes a median or keys by name. A byte comparison of these files across machines would fail for this reason alone. |
| `np.argsort(p)` in BH | `stats.bh_fdr` | safe, by derivation (A10 audit, case B) |

The Methods sentence on the library ordering now names the site-control folds
too.

### E3, E4, E5, 2026-09-17: three more NSCLC choices, one at a time

`pipeline/results/nsclc_sensitivity/`, macOS, the reference gate held bit for
bit. Now in Results with checker authorities.

| variant | ISI [95% CI] | above null | reading |
|---|---|---|---|
| E3, null matched on size only | 0.3041 [0.2434, 0.3536] | 16/16 | expression matching lowers the index by 0.012 (what matching buys) |
| E4, 400 null draws in the bootstrap | unchanged; CI [0.2336, 0.3512] | 16/16 | interval moves by at most 0.0019 |
| E4, all 1,000 draws | unchanged; CI [0.2333, 0.3511] | 16/16 | at most 0.0022; the registered 200-draw subsample is not a material choice |
| E5, 3 folds | 0.2188 [0.1673, 0.2722] | 15/16 | −0.073 |
| E5, 10 folds | 0.3154 [0.2607, 0.3711] | 16/16 | +0.023 |

The index rises with K: 0.219, 0.292 and 0.315 at 3, 5 and 10 folds. More
training data per fold predicts the curated signatures better than their
nulls. The registered K = 5 is therefore part of the estimand's definition.
The sign holds at every K, but the value does not transfer across K. No K
was refused: 68 NSCLC sites support 10 site-disjoint folds.

### E9, 2026-09-17: the mean-z null is far from the reliability cliff (NSCLC)

HPC4 array 130185, split_seed 0, ISI-only configuration, fast path; each
task's index equals its cohort's corrected-ordering primary bit for bit
(`pipeline/results/null_reliability/`). A draw is dropped when its reliability
is at or below 0.01 or when |r / sqrt(α)| reaches 1.

| cohort | null draws | smallest α | smallest 1st percentile of α, any signature | largest \|r / sqrt(α)\| | draws with ratio ≥ 0.9 | dropped |
|---|---|---|---|---|---|---|
| NSCLC (task 0, 258 s) | 16,000 | 0.0586 | 0.226 | 0.553 | 0 | 0 |
| pan-cancer (task 1, 4,148 s) | 16,000 | 0.172 | 0.243 | 0.390 | 0 | 0 |

The nearest draw is 5.9 times the reliability floor in NSCLC and 17 times it
pan-cancer. The ratio guard is never approached. The mean-z arm is far from
the cliff that removes every ssGSEA draw.

Timing, recorded for D2: task 1 is the first whole pan-cancer index timed with
the fast path, 4,129 s on cpu-848 against 11,022 s for the full run before it
(different configurations: ISI-only here). The global-axis run took 1,217 s:
with the axis removed globally, the 19 untyped patients no longer have a NaN
axis, so no null column needs refitting on its own.

### E11, E12, E19–E22, 2026-09-17: statements checked where they are made

- **E11.** `scripts/31_outcome_category_test.py` →
  `pipeline/results/outcome_category_test.json`. All five pan-cancer outcome
  winners lie inside the six signatures MSigDB files under proliferation,
  development or pathway.
  - By exact enumeration, P(all five inside) = 6/4,368 = 0.0014.
  - The six have a mean axis-residualised excess 0.0672 above the other ten.
    No other of the C(16, 6) = 8,008 relabellings does better (one-sided
    P = 1/8,008).
  - Reported in Results as descriptive, with both caveats: the grouping was
    named after the result, and G2M/E2F share genes, so the signatures are not
    exchangeable.
  - The pre-specified contrast (immune vs other: Mann–Whitney 0.031, Fisher
    0.093) is unchanged and remains the test.
- **E12.** Methods' Outcome arm now states the Somers' D mapping and the
  1/√(n − 3) standard error where they are used (n = patients with a usable
  endpoint; the null mean's sampling error added in quadrature,
  `globalaxis.excess_over_null`). Limitation 6 keeps the caveat.
- **E19.** The ssGSEA depth-table hashes are now quoted in Limitations 8: default
  sort `0bca973cfa26da89` (macOS) and `b98fb11e701e8e98` (Linux); pinned
  `096172c056cc7a34` on both. The checker's two-sided A10 assertion guards the
  artefacts, not the prose. S9 carries the fold hashes only, as its caption
  says.
- **E20.** *[2026-09-24: the macOS 0.0211 below is a misrounding of the stored
  0.021154, which prints 0.0212; its pattern's 0.005 floor hid it until session
  60's sweep. Supplementary Note 1 now says 0.0212; this record keeps its words.]*
  Limitations 8 now says its 0.0128 / 0.0211 / 0.61 are macOS values
  and gives Linux's beside them (0.0130 / 0.0211 / 0.62). All three Linux
  values are checker authorities read from `sort_audit_linux_x86_64.json`,
  with self-test rows.
- **E21.** Checked on 2026-09-17 while writing C9: `26_sort_audit.py --write`
  still records that the Linux artefact must not be rewritten, and why. C9
  writes its per-machine record through `--library-order`, which refuses the
  `sort_audit_*` names.
- **E22.** Limitation 10 says "we are not aware of it being noted elsewhere"
  (CPTAC-3 carries no tissue source site).
  - **Searched 2026-09-17.** PubMed: five queries pairing CPTAC with
    source-site, site-preserved, batch and external-validation terms returned
    no relevant record. As a positive control, the same style of query does
    not reach Howard et al. 2021 either (found only by author and title), so
    abstract-level search is weak here. Europe PMC full text:
    `"CPTAC" AND "tissue source site"` returned 20 records, 10 of them with
    histopathology.
  - **Read in full text:** the passages near CPTAC and "site" in five
    candidates:
    - FLEX, Nat Commun 2025, doi:10.1038/s41467-025-66300-y;
    - Murchan et al., J Pathol Inform 2024, doi:10.1016/j.jpi.2024.100396;
    - Lin et al., Cell Rep Med 2025, doi:10.1016/j.xcrm.2025.102527;
    - Murchan et al., Diagnostics 2024, doi:10.3390/diagnostics14050462;
    - Howard et al., Sci Adv 2024, doi:10.1126/sciadv.adq0856.
  - **What they do:** each pairs TCGA site-preserved or site-aware training
    with CPTAC used WHOLE, as one external cohort, which needs no site
    variable. None runs a preserved-site analysis inside CPTAC or notes the
    missing variable. The one preprint hit (doi:10.1101/2025.01.06.631471) is
    about genomic biomarkers.
  - **Limit:** bioRxiv has no keyword search through the tools available; only
    Europe PMC's preprint index was covered.
  - **Decision:** the sentence is kept as written. It claims only unawareness,
    and the search supports that.

### E14, E15, E17, 2026-09-17: Control C, re-examined

**E14, the NSCLC/pan-cancer gap (14.5% against 44%).** The manuscript said
"most of that gap is the dummy-coding term". That is contradicted by the
stored per-signature tables (S6, S10):

| median across 16 signatures | pan-cancer | NSCLC |
|---|---|---|
| site R², crude | 0.444 | 0.145 |
| site R², adjusted for its number of terms | 0.424 | 0.112 |
| cancer type alone | 0.385 | 0.036 |
| site given cancer type | 0.055 | 0.099 |
| site, calibrated within type (all sites, 200 permutations) | 0.047 (6/16 at p ≤ 0.05) | 0.070 (15/16) |

- Adjustment removes about 0.02–0.03 in both cohorts; that is the
  dummy-coding term, and it cannot explain a gap of 0.30.
- The gap is cancer type: pan-cancer, site stands in for disease.
- Beyond disease, NSCLC carries at least as much label-side site structure as
  pan-cancer.
- The Results paragraph is corrected, and each number is a checker authority
  (`lsv_*`, `e14_*`).
- The calibration records are
  `pipeline/results/control_c_calibration/{pancancer,nsclc}.json`, from HPC4
  array 130200 (23 min and 1 min). Each reproduces script 11's single-probe
  value bit for bit, and each covers all sites rather than Control C's sites
  of at least 10 patients: median crude R² 0.485 and 0.177 on that wider frame.

**E15, site attrition.** `label_site_variance` keeps sites with at least 10
patients (`min_site_n = 10`):

| cohort | sites kept | sites dropped | patients in dropped sites | notes |
|---|---|---|---|---|
| pan-cancer | 202 of 619 | 417 (largest dropped site: 9 patients) | 1,376 | 17 of the 5,792 retained patients have no cancer type and fall out of the complete-case mask, giving 5,775 |
| NSCLC | 31 of 68 | 37 | 148 | 796 |

Purity is not a column of the frame script 11 builds, so Control C never used
it, although `label_site_variance` accepts a purity column. Now stated in the
Results paragraph and the S6/S10 captions.

**E17, plate within site.** The finding is that the zero is not a
measurement.
- `barcodes.annotate` parses `plate` from the barcode's optional plate field (the sixth dash-separated component, as in `TCGA-XX-XXXX-01A-11R-A12B-07`).
  Neither barcode form the inputs use has one: the slide barcode is parsed
  from `sample` (`TCGA-06-0138-01`), and TOIL ids are 15 characters. So
  `plate` is None in all 10,170 rows of `meta.parquet`.
- `_design` then dummy-codes a single NaN level, drops the constant column,
  and returns the site design unchanged. `r2_plate_within_site` is exactly 0
  by construction, in every signature and both cohorts.
- The manuscript, the poster and this audit (A4) used that zero to conclude
  the label-side site artefact "is not technical". The claim is withdrawn in
  all three, the S6/S10 captions say the plate columns are evidence of
  nothing, and the checker's plate pattern is removed.
- **Not changed in code, deliberately.** `11_close_science_gaps.py` writes
  `label_site_variance.csv` into the frozen result directories, so making
  `label_site_variance` return NaN when the plate is absent would change
  frozen artefacts on the next run.
- **BLOCKED-AUTHOR:** fix the function (NaN plus a note when a categorical
  covariate is entirely missing), and decide whether to source plate
  identifiers from GDC aliquot metadata, which would be a new input.
  *[2026-09-24, decided by the author: fixed — the plate rung is NaN with a
  note whenever the plate column has fewer than two non-missing levels, and
  the frozen `label_site_variance.csv` files are not rewritten. Plate
  identifiers are NOT sourced from GDC: that would be a new, unregistered
  input for a control the manuscript has already withdrawn.]*

### C3 and D2, 2026-09-17: the whole path, measured, from an empty data directory

HPC4 job 130192 (`aacr27-s49repro3`, a clean clone of 95a0067, `data/raw`
symlinked to the inputs attempt 4 had already downloaded and hash-verified),
`bash reproduce.sh`, **exit 0 in 8,751 s = 2 h 25 m 51 s**. The six-item science
list's item 5 is closed. Log:
`pipeline/results/session49_diagnostics/job_logs/repro49c_130192.out`.

| stage | elapsed | result |
|---|---|---|
| versions, provenance | 15 s | lockfile OK; 211 artefacts hash to their frozen values |
| fetch, hash-verify | 15 s | all 8 raw inputs VERIFIED |
| build `data/interim` | 163 s | 5 of 6 byte-identical; `cohort_nsclc.parquet` the recorded Linux build (A13); exit 0 |
| determinism, sort audits | 42 s | fold hash matches both machines; no dangerous ordering |
| test suite | 885 s | 177 passed, 1 skipped (178 at that commit) |
| NSCLC (05) | 617 s | |
| pan-TCGA (09) | 6,944 s | with the shared-mask fast path |
| ancestry (06) | 57 s | |
| figures (10) | 26 s | five figures from the frozen results |

**Both cohorts reproduced the corrected-ordering values bit for bit**: NSCLC
0.292212 [0.2355, 0.3527] and pan-TCGA 0.296849 [0.2749, 0.3189], each with
delta +0.000000 against `results/*_stablesort/`. The NSCLC reference was
computed on macOS and the pan-cancer one on HPC4, so this is also a second
cross-platform agreement for NSCLC's index, now through the whole path from raw
files rather than through one script.

**D2, the fast path's cost, measured.** `09_run_pancancer.py` took 6,944 s here
against 11,022 s for the same script without the fast path (job 129159): a
1.59x speed-up on the whole run, and the first end-to-end figure with the fast
path on. The ISI-only shard of the same cohort takes 4,129 s (E9's capture,
job 130185_1). macOS ran 17,118 s before the fast path and has not been
re-timed. Every runtime sentence in `reproduce.sh` and `pipeline/README.md` now
quotes the measured figures with their platform and date.

### E6, 2026-09-17: a third scorer, and what it says about ssGSEA's reversal

`signatures.score_plage` adds PLAGE (Tomfohr et al. 2005; the leading singular
vector of a set's gene-standardised matrix) as an opt-in scorer, and
`scripts/32_third_scorer.py` runs it through `13_scorer_sensitivity.py`'s own
functions. HPC4 job 130215, NSCLC, n_null 1000, current code
(`pipeline/results/scorer_sensitivity_plage/`); the smoke stage at n_null 10 ran
first so a crash would cost minutes.

| arm (NSCLC, current code) | value | above null |
|---|---|---|
| mean-z, uncorrected (control, same run) | +0.2963 [0.2442, 0.3514] | 16/16 |
| PLAGE, uncorrected | +0.1353 [0.0889, 0.1834] | 13/16 |
| ssGSEA, uncorrected (stored, pinned ordering) | +0.1738 [0.1333, 0.2102] | 15/16 |
| PLAGE, registered disattenuated | NOT COMPUTED | — |
| ssGSEA, registered disattenuated | NOT COMPUTED | — |
| mean-z, split-half reconstruction | +0.2919 | 16/16 |
| PLAGE, split-half reconstruction | +0.0772 | 11/16 |
| ssGSEA, split-half reconstruction | −0.1333 | 0/16 |

Two readings change:
- The registered estimand is undefined under BOTH alternative scorers, for one
  reason: the reliability estimator assumes a mean composite. Limitation 9's
  claim is now about that assumption, not about ssGSEA.
- ssGSEA's sign reversal does NOT generalise. PLAGE, which is also not a mean
  composite, stays positive under the same scorer-agnostic estimator. So the
  reversal is a property of rank-walk scores of random sets, not of "any
  alternative to mean-z". The attenuation of magnitude does generalise: both
  alternatives sit well below mean-z.

### E8, 2026-09-17: McDonald's omega beside Cronbach's alpha

`scripts/33_omega_reliability.py`, macOS, 586 s, on the registered residualised
items of `results/nsclc_v3_stablesort` (`pipeline/results/omega_reliability/`).
Gate: the script's own alpha reproduces the run's stored `alpha_observed` to
1.1e-16.

| median over 16 signatures | Cronbach alpha | McDonald omega |
|---|---|---|
| curated sets | 0.9773 | 0.9773 |
| random sets (25 draws each) | 0.8358 | 0.2270 |
| curated − random gap | 0.1264 | 0.7333 |
| reconstructed index | 0.2917 [0.2584, 0.3251] | 0.2147 [0.1739, 0.2554] |

15 of 16 signatures still exceed their null under omega. The two estimators
agree on direction and on 15 or 16 of 16, and differ in magnitude by about a
quarter.

**Read omega as a bound, not an improvement.** Omega-total from a one-factor
model is interpretable where one factor fits. For curated modules it does
(omega equals alpha to four decimals). For random sets after the axis is
removed it does not, and that is precisely why their omega is so low: what
remains is multidimensional, so a single factor captures little of the
covariance the unweighted mean actually sums. Alpha, which does not model a
factor, is the more defensible choice for the null side, and the registered
estimand keeps it. The honest statement is that the reliability ESTIMATOR moves
the index by about a quarter — more than switching scorers within the
mean-composite family does.

### E16, 2026-09-17: the site AUROC does not stop being ~1 anywhere we can look

`cluster/e16_site_auroc_stress.py` recomputes the registered
`models.site_prediction_control` on random patient subsamples (HPC4 arrays
130207 and 130212; `pipeline/results/site_auroc_stress/`). Fractions below 1
are drawn three times; the medians below are over those repeats.

| cohort, threshold | full cohort | half | quarter |
|---|---|---|---|
| pan-cancer, sites ≥ 20 patients | 0.9985 (98 sites) | 0.9976 (40–42) | 0.9940 (10–14) |
| NSCLC, sites ≥ 20 patients | 0.9894 (15 sites) | 0.9805 (4–5) | 0.894–0.972 (1–2) |
| pan-cancer, sites ≥ 5 patients | 0.9981 (320 sites, 97% above 0.9) | — | 0.9963–0.9969 (114–123) |
| NSCLC, sites ≥ 5 patients | 0.9867 (47 sites, 94% above 0.9) | — | 0.9623–0.9783 (17–20) |

- Below a quarter of either cohort, NO site keeps the 20 patients the control
  requires, so the honest answer to "where does it stop being ~1" is that the
  binding limit is site eligibility, not discriminability.
- At a 5-patient threshold the median is still 0.987 and 0.998, so sites
  contributing five patients are identifiable.
- A LATENT BUG, found by this run and not fixed:
  `models.site_prediction_control` sorts its result by `auroc`, so with no
  qualifying site it raises `KeyError: 'auroc'` instead of returning an empty
  table (array 130207, tasks 0 and 1). No registered run reaches it. The
  diagnostic catches it and records "no site evaluated".
  **BLOCKED-AUTHOR:** return an empty frame with its columns, or leave it.
  *[2026-09-24, decided by the author: it returns an empty frame with its
  columns; a test pins it.]*
- The full-cohort NSCLC median here is 0.9894 on Linux against the frozen
  0.9916 on macOS. That is A12 again: this control uses `GroupKFold` too
  (F6.3).

### F-workstream decisions, 2026-09-17 (session 49)

Closed by measurement:
- **F1.4.** The cohort-fingerprint MATCH path had never read a real artefact.
  It has now: `read_and_check_cohort` on `results/nsclc_sensitivity/reference`
  and `/k3` returns `status: match`, 944 patients, sha `54c7688a…`. Those runs
  are the first written by code that records the fingerprint.
- **F3.1, F3.2 — measured 2026-09-17; a point-in-time reading, not a standing
  state.** `--coverage --strict` exits 0: **753 of 10,157 literals
  checked (7%), 817 claim comparisons, 0 mismatches**, and no document has
  fewer checked literals than its floor. PATTERN REACH is **205 patterns: 0
  matching nothing, 172 reaching one document, 33 reaching two or more**
  (session 48: 117 / 0 / 86 / 31). Nothing shrank, and no pattern is vacuous.
- **F5.1.** The poster's figure boxes measured again in a browser after panels 4
  and 6 changed: 714 / 687 / 501 / 409 px, unchanged, and no section overflows
  (`scrollHeight == clientHeight` for all seven).
- **F6.11.** `pipeline/README.md` now states that there is no CI, that every
  check is run by hand or by `reproduce.sh --check`, and that the
  platform-dependence checks only mean anything when run on each machine.
- **F8.11.** No tracked text file has trailing whitespace, and none has CRLF
  line endings.
- **F11.2.** `--tiff` now writes RGB at 600 dpi (`save_tiff_rgb`): print
  production expects RGB, and matplotlib's TIFF path was RGBA. Translucent
  pixels are composited onto white; a fully opaque figure loses nothing. Tested,
  with the composite branch proven able to fail.
- **F11.6.** The journal-size arithmetic, for the first time. Tick labels are
  7 pt in `rcParams`; the figures are 6.4 to 7.8 in wide natively, so scaled
  into a 6.75 in double column they print at **6.1 to 7.4 pt** — below the 8 pt
  floor most guidance gives, with the widest figure worst. Raising
  `xtick.labelsize` would re-render every committed figure, so it is the
  author's call, alongside TODO-PRINT-6 (poster legibility, measured
  2026-09-08). **BLOCKED-AUTHOR.**
- **F11.8.** The deposit figure test names all five stems explicitly
  (`figure0_schematic` … `figure4_controls`), so a tree that rendered only
  figure 0 fails it.
- **F12.1.** Every source without an upstream version string is pinned by
  content instead: the TOIL download and the gene map by sha256, ABSOLUTE and
  the ancestry calls by GDC file UUID and sha256, TCGA-CDR inside the
  embeddings parquet by HuggingFace revision `073115403c…` (`01_fetch_data.py`).
- **F12.4.** The embeddings dataset card was read: CC-BY-4.0 with no attribution
  text of its own. The Data Availability Statement now names the dataset id, the
  revision and the licence beside reference 2.
- **F12.5.** Reference 12 resolves (PubMed 32396860, doi
  10.1016/j.ccell.2020.04.012, *Cancer Cell* 37(5):639-654.e6), and its UCSF
  group is the source of the ancestry calls used.
- **F12.10.** The "no gated model weights" sweep, re-run after the builder
  landed: the only download in the repository is the HuggingFace DATASET
  (`repo_type="dataset"`). No model weights, no `timm`, no torch.
- **F14.2.** The submission checklist's statistical-software gap is closed by
  the Software paragraph (SciPy for the correlations, Fisher z quantiles,
  permutation percentiles and χ²; statsmodels for OLS and Cox; BH in this
  project's own code), and §6's other two items are updated.
- **F14.3.** §1's live claims about the public repository re-verified from
  GitHub's API: 0 tags, licence `NOASSERTION`, last push 2026-09-08T23:08Z.
- **F7.9.** S5's caption now states its composition: 98 pan-TCGA sites plus the
  15 NSCLC sites twice (embedding and ComBat-corrected), which is 128, and says
  pan-TCGA has no ComBat arm.
- **F6.5.** Already satisfied: a test asserts every package `ENVIRONMENT.md`
  reports is pinned in the lockfile or listed in `_ENV_PACKAGES_NOT_PINNED`
  with a reason, in both directions.

- **F5.6, F11.3.** Figure 3 draws no p-value; the drawn one is figure 2's
  Spearman (`rho = -0.50, p ... 0.048`). It is already checked, by the route
  F11.3 called (B): the committed SVG carries that text, and the render-sync
  test compares committed SVGs against a fresh render on every platform, so a
  change in the drawn p fails the suite. The checker keeps its
  dependency-freedom — no scipy — and route (A) is not needed.
- **F10.9.** `_DRAWS_TEXT` was proven complete against the script: it names
  twelve calls, the script makes seven, and all seven are named. That reading
  is now a test: `test_draws_text_covers_every_text_call_the_figure_script_makes`
  walks the AST for calls that are text-drawing by shape (ending in `label`,
  `labels` or `title`, or one of matplotlib's text entry points) and asserts
  each is classified, with a vacuity floor. Proven able to fail by removing
  three members.
- **F13.4.** The go/no-go countdowns in `16-YOUR-TASKS.md` were recomputed from
  the machine's date on 2026-09-17: 54 days to the abstract deadline, 26 to the
  membership cutoff, 81 to early registration. `15-CHECKLIST.md`'s line is
  dated 2026-09-16 and stays as a dated statement.
- **F14.6.** The placeholder list, re-located by content rather than by memory:
  `09-PAPER-DRAFT.md`'s Zenodo DOI (`10.5281/zenodo.XXXXXXX`, Part B item 4);
  `poster/poster.html`'s TODO-PRINT-1 (affiliation) and TODO-PRINT-2 (`[TK]`
  control number, assigned at submission); TODO-PRINT-6 (A0 legibility, now
  with the journal-size figure beside it, F11.6); five authors' CRediT roles
  and every affiliation as `[NEEDS AUTHOR]`; `submission/COVER-LETTER.md`'s
  repository DOI, bioRxiv DOI and one ORCID. Every one is Part B.

- **F8.8.** `cluster/` is still excluded from the snapshot wholesale: the
  2026-09-17 rebuild's 296 files contain no path under `cluster/`, checked by
  searching the staged tree rather than by reading the builder.
- **F11.4.** Confirmed: the figure-spelling guard has two independent routes,
  `test_rendered_figure_text_uses_american_spelling` (string literals reaching a
  drawing call, from the AST) and
  `test_committed_figure_pdfs_are_free_of_british_spelling` (text extracted from
  the committed PDFs). Neither can pass by sharing the other's defect.
- **F13.5.** Swept: every claim of cross-platform agreement in the manuscript is
  scoped to the index and the preserved-site partition, and Methods now names
  the site-classification control as a second thing the library ordering
  reaches. Reproducibility gains the whole-path statement C3 earned, with the
  same scope attached.

Declined in writing, with the reason:
- **F8.10.** segno is deliberately absent from the lockfile (recorded in
  `poster/make_qr.py`), and `make_qr.py --check` verifies the ENCODED CONTENT —
  URL, error level, quiet zone — so the library's version cannot affect the
  claim. Pinning it would add a dependency to satisfy a check that does not
  depend on it.
- **F12.3.** An assertion that the KEGG prefix filter "bites" would fire by
  design on Hallmark, where it is a measured no-op (0 of 50 sets restricted).
  The no-op is recorded instead.
- **F10.8.** Pytest markers for the deposit-staging tests are declined:
  `reproduce.sh --check` runs the whole suite, and a marker is a way for those
  tests to stop running without anyone noticing. The deposit's own expected
  failures are recorded in `pipeline/README.md`'s table instead.
- **F10.10.** `_ON_WORKING_TREE` gates three tests; two others run everywhere
  and fail inside a deposit BY DESIGN, because they assert that the working
  tree's staged snapshot matches the deposit. That trade-off is recorded in the
  deposit table: a deposit reader sees those failures and the table says which
  they are.
- **F8.9.** `results/PROVENANCE.json` stays unclassified in
  `results/README.md`, tolerated explicitly by the classification test: it is
  the manifest OF the results, not a result. Classifying it would make the map
  describe itself.

### F1.1, F4.2, F9.2, 2026-09-17 (session 50): Limitation 8 read cold, and an argument it carried was wrong

Limitation 8 had been edited about eighteen times in thirty days and never read
whole. Read cold, it is 150 lines against 2–20 for every other limitation, and
it carried one argument (item 1) and five statements (items 2 to 6) that the
stored artefacts contradict; a sixth statement sat in Results (item 7). Each was
checked against a file, not against memory.

1. **The pan-cancer scaling argument is contradicted by the folds.** It said
   the sort fix moved pan-cancer by fewer partition sds than NSCLC because
   pan-TCGA's 619 sites mean any tie resolved differently "perturbs a
   correspondingly smaller share of the partition". Nothing had measured the
   share. `scripts/34_sort_fix_fold_change.py` reads the per-patient
   preserved-site folds that both cohorts' frozen and corrected runs store:
   the fix moved **53.2%** of pan-cancer patients (3,812 of 7,168; adjusted
   Rand index 0.131, a nearly unrelated partition) against **22.9%** of NSCLC
   patients (216 of 944; ARI 0.537). The mechanism is also measured: the fix
   reorders tied sites only (no patient in an untied site moved, and every
   fold kept its size, in both cohorts), and pan-TCGA has 600 of 619 sites
   tied, holding 79.5% of patients, against 51 of 68 and 38.1%. So the larger
   cohort had MORE of its partition changed and moved LESS. The paragraph now
   states the measurement.
2. **"About as much as changing the partition does, and no more"** (the NSCLC
   move of 0.0260) overstated the comparison. Of the 276 pairs among the 24
   NSCLC partitions, 3 differ by at least that much; the median pair differs by
   0.0080. The pan-cancer move (0.0057) is typical: 153 of 276 pairs differ by
   more. The NSCLC sentence now says the move is large for a change of
   partition and still inside the honest interval.
3. **Two yardstick ratios were wrong at one decimal**, and none of the six had
   an authority: the platform gap is 2.457 standard deviations of the
   six-partition sweep, where the manuscript wrote 2.4, and the
   NSCLC-over-pan-cancer contrast is 3.347, where it wrote 3.4. The other four —
   the NSCLC move at 2.9 and 3.4 sds, the pan-cancer move at 0.87 (also written
   0.9) and 0.79 — hold. All now have `sffc_*`
   authorities with self-test rows, two of which inject exactly the old values.
   The "2.4×" also appeared four times in three other documents. In
   pipeline/README.md it was a present-tense statement in the deposited "start
   here" file, and is corrected there. The other three (this file's A9
   section, 15-CHECKLIST.md twice) are dated records of the session-26/27
   arithmetic, which divided the four-decimal gap by the four-decimal sd
   (2.449), and are left as written.
4. **"216 of 944 patients change fold" sat in a sentence about the two
   machines**, where it read as a macOS-against-Linux count. It is macOS's
   shipped-against-stable count (`split_determinism_darwin_arm64.json`); Linux's
   own is 180, and no machine-against-machine count exists. The sentence now
   says which comparison each number is.
5. **"Every secondary quantity in this paper is the macOS value, and none is
   taken from the Linux runs"** was false after session 49. E16's subsample
   stress test ran on HPC4, and its full-cohort NSCLC median (0.989) differs
   from the frozen macOS 0.992 because the site control's folds come from the
   same library routine (F6.3). The sentence is now scoped to Δr, Δ-MAE and the
   covariate baselines, and both it and Negative controls say why the two
   medians differ.
6. **"No unpinned order-dependent operation anywhere in the pipeline"** was a
   law over a set it did not cover: the library ordering (case F) is
   order-dependent and deliberately documented rather than pinned. The
   sentence now says the pipeline's own code and names the library call sites.
7. **Results, "the same machine".** The sentence that one pan-cancer partition
   in 24 fails the covariate regression was true on macOS and unscoped; the
   24-partition paragraph above it reports five Linux failures. It now names
   the machine and says the five Linux failures are that census's.

Also checked and kept: "validated" (now "validated for the index under the
corrected ordering", with the next paragraph's scope unchanged); the
platform-justification paragraph (the memory measurement is what licenses a
Linux number, so it stays, but it is the limitation's longest aside); F9.2's
question — the pan-cancer sd is used twice as a yardstick, once per sweep and
labelled, and never counted inside the interval it is compared with. Whether
Limitation 8 should move most of its account to a supplementary note is an
editorial choice for the author, recorded in the ledger.

### F1.6, 2026-09-17 (session 50): swap and compressor occupancy track nothing we report

The question carried since session 46: macOS's swap total read 7,168 to
30,720 MB across sessions, and compressor occupancy varied by several GiB. Does
either track anything the project reports? The only quantity the memory samples
can be set against is suite wall time, so each of the eight committed
`memsample_*.log` files was paired with the `full_suite_*.log` of the same tag
and summarised over its own samples.

| tag | suite (s) | mean available GiB | mean compressor GiB | swap total MB | mean swap used MB | mean load1 | power |
|---|---|---|---|---|---|---|---|
| 20260915_run47 | 255.00 | 4.66 | 9.09 | 18432 | 16970 | 5.1 | not recorded |
| 20260915_run47b | 257.49 | 4.61 | 9.06 | 18432 | 18413 | 5.7 | not recorded |
| 20260916_run48 | 259.83 | 8.12 | 6.38 | 30720 | 28567 | 9.7 | AC |
| 20260916_run48b | 465.81 | 6.95 | 6.64 | 30720 | 29970 | 25.0 | battery |
| 20260916_run49 | 280.97 | 5.95 | 5.84 | 11264 | 11094 | 6.7 | AC |
| 20260916_run49b | 287.98 | 6.36 | 5.65 | 11264 | 9857 | 19.8 | AC |
| 20260917_run49c | 280.02 | 5.68 | 6.61 | 6144 | 4694 | 8.2 | AC |
| 20260917_run49d | 277.52 | 4.63 | 5.96 | 6144 | 5883 | 9.3 | AC |

Spearman correlation with suite wall time, n = 8: compressor −0.571 (p 0.139),
swap used −0.048 (0.911), swap total −0.098 (0.818), available memory +0.548
(0.160), load1 +0.762 (0.028). Among the five AC runs alone, no covariate
reaches |ρ| 0.7 or p below 0.18.

**Answer: no.** Compressor occupancy moves the wrong way for a pressure effect
(the slower runs had less compressed memory) and swap does not move with wall
time at all. The swap TOTAL is the size of the swap files macOS has allocated
since boot, set by whatever else the machine ran, not by this project. Only
load correlates, and that rests on the single battery run, where load and power
are confounded, which is F2.2's unresolved question and stays declined there.
The memory fields remain recorded because they cost nothing and would show a
real pressure episode (a pan-cancer run, not a suite); as a covariate for suite
timing they are closed.

### F3.4, F3.5, F4.7, F4.10, F4.12, F11.5, F12.7, 2026-09-17 (session 50): reading items, and what they found

- **F3.4, the authority set against the provenance set.** Measured by wrapping
  every read the checker makes: `authority()` opens 55 files under `results/`,
  and six of them were not hashed (E6's PLAGE summary, E8's omega summary, E16's
  four subsample records). Session 49 hashed its authorities before those three
  analyses landed. PROVENANCE 212 → 220 (+8, 0 removed, 0 changed), and a test
  now keeps the two sets reconciled; it failed on the real tree, naming exactly
  those six, before the manifest was rewritten. A second gap: a claim row whose
  authority key is absent was skipped before `--strict` counted it, so a
  vanished authority left its pattern checking nothing, silently. Zero such keys
  today; `--strict` now fails on any, and a test proves the helper can see one.
  The other sets: EXEMPT is reconciled two ways by existing tests; TRACKED_DIRS
  mirrors the results/README.md classification by measurement (every untracked
  directory is DIAGNOSTIC, a written exclusion — figures F6.9, supplementary
  F7.5 — the synthetic `demo/`, or empty), not by a gate; the snapshot keep-sets'
  complement is the four `deposited: false` artefacts and the regenerable TSVs,
  both decided.
- **F3.5, "as of" claims outside run rows.** Eleven read. Two fixed: the
  supplementary inventory said every manuscript value from its artefacts was
  mechanically checked, which `--coverage` contradicts; the preprint plan's
  "123 staged files unpublished as of 2026-09-17" was already false on that
  date, and is now a record of session 49's check. The rest are dated and true,
  or checked live.
- **F4.7, the ten S-number citations.** Each points at what its caption
  describes. S9's pointer sat in the ssGSEA rank-table paragraph while the
  fold-assignment hashing it tabulates is claimed a page earlier; moved there.
  Open for the author: first citation runs S3, S4, S2, S6, S10, S8, S5, S1, S7,
  S9, and most journals number supplementary items in order of first citation.
  *[2026-09-24: renumbered in exactly this order, so the manuscript now cites
  S1 to S10 in sequence; see the note at the top of this file.]*
- **F11.5, Figure 4's asymmetry.** Not an asymmetry: no figure was cited in the
  running text at all; all four existed only as caption blocks. Callouts added in
  order (1A, 1B, 2B, 2A, 3, 4A, 4B). Figure 4B's content is first discussed in
  Results' opening subsection, before Figure 1, so it is called out where its
  caption sits. Its absence from the poster is deliberate (TODO-PRINT-5).
- **F4.10, the Code Availability blockquote.** Both variants said "Reported
  values were produced on macOS 15 / arm64" *[corrected 2026-09-21: 26.5.2, not
  15; see the session 54 block]*; the corrected-ordering pan-cancer
  run, the pan-cancer 24-partition sweep and several sensitivity analyses were
  produced on Linux. Both now say which. The closing paragraph said the estimate
  is "platform-dependent in its fourth digit" and Limitation 8's heading said
  "third decimal": the 0.0218 gap between 0.3182 and 0.2964 is in the SECOND
  decimal, and the wording had been wrong since the commit that introduced it
  with those numbers (`f95b96a`). Corrected in the manuscript (four places), the
  snapshot README, pipeline/README.md and the abstract notes; this file's A9
  section keeps its dated wording.
- **F4.12, Contributions, Disclosures, Acknowledgments.** The Acknowledgments'
  AI sentence said analyses were "checked by the author"; with six authors that
  now reads "the corresponding author", which is true whatever the co-authors'
  roles turn out to be. "No disclosures were reported." speaks for all six while
  five disclosure forms are uncollected, so it carries a `[NEEDS AUTHOR]` like
  the five contribution placeholders. The submission checklist still said
  "Single-author, so this is mechanical"; fixed. The cover letter's "copyright
  retained by the author" is the author's wording and is listed in
  16-YOUR-TASKS.md.
- **F12.7, the DX-slide figure.** The pre-registration addendum and
  `07_triage_cptac.py` said the discovery side is 99.48% diagnostic slides
  (10,117 of 10,170), which contradicted the manuscript's "every slide is a
  primary-tumor diagnostic slide". Measured from the slide filenames: all 10,170
  are diagnostic slides and all are sample type 01. The 53 the parse missed carry
  lettered suffixes (DXA to DXU) across 7 patients in SARC, THYM, LGG and TGCT,
  which a pattern requiring a digit after "DX" does not match. The manuscript
  was right; the addendum (never part of the tagged protocol) and the script now
  say so. Also from the same read: Limitation 8 gave the expression matrix as
  10,535 × 41,047, counting the parquet file's stored index as a gene column; it
  is 10,535 × 41,046, as the Methods say.

### F6.7, 2026-09-17 (session 50): what changed in the estimator since `prereg-2026-08-18`

239 commits separate the tag (`2d3f6fa`, 2026-08-18 22:13 CDT) from session
50's start. The question a reviewer asks is narrower: which changes to
`pipeline/src/` could have moved a reported number, and is each one declared?

**Between the tag and the frozen runs** (NSCLC v3 committed 2026-08-19 in
`2d4f200`, pan-TCGA v3 on 2026-08-20 in `3c9c512`), four commits touched `src/`:

| commit | change | effect on a reported number | declared |
|---|---|---|---|
| `2408bca` | new `ancestry.py` | the supplementary ancestry arm only | Results, S1 |
| `39664aa` | outcome arm stratified by cancer type; bootstrap complete-case mask | both reported outcome counts; the pan-cancer interval (NaN without it) | **not until 2026-09-17**: protocol Appendix A and two Methods sentences, written this session |
| `2d4f200` | MDE computed from the outcome interval | adds the MDE, a new descriptive quantity | Results report it as such |
| `ed04f4a` | null draws persisted to `null_draws.npz` | none (storage) | — |

The NSCLC primary the tag recorded (0.3182 [0.2614, 0.3763]) is the frozen
value to the digit, and NSCLC's outcome count is 0 of 32 under both the tagged
and the stratified arm.

**After the frozen runs**, eight commits touched `src/` (`39d618b`, `507b550`,
`64cac37`, `ef4de1c`, `db7b136`, `20a1748`, `80660da`, `87a76a6`). None can
alter a frozen result: those are hashed in PROVENANCE.json and the suite fails if
one changes. Where one alters what the current code computes, both numbers are
reported (A9/A10's corrected ordering, Option B), or the change is shown to be
bitwise-identical before it ships (the ssGSEA identity, the null-refit fast
path), or it is opt-in for a sensitivity analysis (`split_seed`, PLAGE).

So the one undeclared change was `39664aa`, and it is now declared. The
report stops at `src/`: scripts, tests and documents changed far more, but the
frozen numbers depend on them only through `src/` and the stored inputs.
*(This sentence was restored here on 2026-09-17 by session 51. Five of session 50's
later appends (`7d29447`, `c4dacb7`, `bb76a05`, `6aee608`, `77df3fe`) each
landed between its first word and the rest, which left it stranded at the
end of the F7.4 bullet below.)* The same
read showed the tagged protocol never mentions a pan-cancer cohort; Limitation
8's "pre-registered" pan-cancer value and configuration now say "frozen", and
Appendix A records the cohort's addition as its third item.

### F7.8, F10.3, F10.4, 2026-09-17 (session 50)

- **F7.8, DECLINED.** A `supplementary/index.json` would be a third generated
  copy of what `MANIFEST.md` (S-number, file, shape, sources) and `CAPTIONS.md`
  (caption) already carry, both written from one dict in
  `25_build_supplementary.py`. Nothing in the repository or the submission plan
  consumes one. If a scripted consumer appears, it is a short addition to the
  same generator reading the same dict.
- **F10.3, LISTED.** A mechanical scan of `tests/test_pipeline.py` for test
  functions that discover a set (glob, rglob, iterdir, git ls-files, DOCS and
  the document-set constants) and assert it empty finds 14. The definition is
  the scan's, written here so it can be rerun: a body matching both a discovery
  pattern and an emptiness assertion (`assert not x`, `== []`, `== set()`).
- **F10.4.** Of the 14, the scan flagged five without a floor. Read one by one:
  two are can-fail companions themselves, one runs on a fixed temporary tree,
  one has a floor the scan's pattern missed (`scanned >= 9`), and one,
  `test_every_prose_document_is_scanned`, had none: an empty `_tracked_prose`
  would have left it green over nothing. It now asserts at least 30 tracked
  and 30 scanned documents (34 and 33 measured).
- **F3.9, MEASURED, no change.** The FRACTIONS net (every N/M on a line saying
  "signature(s) beat/exceed/beating", "ISI |" or "beats null" must be 16/16)
  checks 7 lines across the scanned documents. Nine further lines carry a
  genuine non-16/16 fraction (15/16, 13/16, 11/16, 1/16, 0/16) beside null
  wording and escape the net only because the fraction is not written next to
  "signature(s) beat": two in the manuscript ("0/16 beating null", "15/16 vs
  16/16 beating"), one in this file, six in results/README.md. That is the
  false-positive exposure: one rewording away, and a trip is loud (a MISMATCH
  naming the line) and cheap to resolve ("15 of 16"), which is how #31 resolved
  it. Narrowing the guard would trade that for a blind spot.
- **F3.6, DECLINED with the measured count.** The five law-grammar phrases the
  ledger named occur 93 times in the checked documents ("the only" 58, "there
  is no" 20, "nothing else" 15, "in every case" 0, "never needs editing again"
  0). None is attached to a number an authority could check, so a WORD_CLAIMS
  row cannot test them, and a lint would flag all 93 for a human. The sweep was
  done by reading instead, over the deposited documents: the protocol's "no
  FFPE-vs-frozen choice" (dated, verified at the time), the manuscript's "no
  dichotomization path" (true: `outcome.py` has none), pipeline/README.md's "no
  continuous integration" (true: no CI configuration is tracked), and Results'
  "for no other" (true, now that it names macOS, F1.1 item 7).
- **F6.10, audited.** `20_check_versions.py` compares numpy, pandas,
  scikit-learn and scipy with the lockfile and prints every BLAS the process
  loads, which is how D5 saw numpy's OpenBLAS 0.3.27 and scipy's 0.3.28 side by
  side on HPC4; BLAS stays informational because it legitimately differs by
  machine and the corrected-ordering index agrees across two machines anyway
  (C3). The gap was the watch list: statsmodels, pinned at 0.14.4 and used in
  `globalaxis`, `decomposition`, `outcome` and `experiment`, was never compared.
  It is now watched; macOS and HPC4's environment both run 0.14.4, as do the
  other four packages on both.
- **F8.2, extended after measuring.** Over the 33 scanned prose documents the
  bracket-balance scan's four candidate extensions were measured before any was
  added: `{}` 0 hits, curly double quotes 0, straight double quotes 0, and
  backticks 22 per line but 0 per paragraph (inline code wraps). All four now
  run, the three marks per paragraph, each with a planted-break can-fail and a
  benign-shape check.
- **F3.8, and a wrong number it found.** WRAP- and PARA-scoped patterns match
  whitespace-normalised paragraphs, so wrapping cannot hide their claims. The
  exposure is the 100 LINE-scoped rows: an occurrence whose anchor and number
  sit on different lines is simply not seen. Measured by matching every
  line-scoped row per line and per normalised paragraph over the scanned
  documents: 12 row-and-document pairs have occurrences only the paragraph
  view finds. Three are artefacts of the paragraph view (a slow-only key, an
  NSCLC value in an NSCLC paragraph, a Fisher p beside a Mann-Whitney anchor);
  five are correct values, unchecked (the rotation p in 15-CHECKLIST.md, the
  small-panel Spearman and a tied-site count in this file, cohort sizes in the
  cover letter, the abstract draft and the NextGen abstract). **One was wrong:**
  the cover letter paired the counts 0 of 6 immune and 5 of 10 other
  signatures with a Mann–Whitney p of 0.016. Those counts give 0.031, the
  manuscript's value;
  0.016 is the p with IL6/JAK/STAT3 counted as immune (0 of 7 against 5 of 9),
  which the manuscript reports only as the Figure 3 caption's contrast. The line
  broke between "Mann–Whitney" and "*P*", so the row whose guard needs
  "MSigDB-immune" on the same line never saw it. Fixed to 0.031; a WRAP row
  anchored on the counts now reaches it, and a self-test row injects exactly
  0.016 (verified: 0 failures on the real tree, 1 with the old value). The five
  correct-but-unchecked occurrences now have rows of their own too (WRAP rows
  for the four that wrap, a line-scoped row for the audit heading whose guard
  sat on another line), each with a self-test row injected in the raw wrapped
  form; all six caught.
- **F7.4.** CAPTIONS.md had 1 of its 58 numeric literals checked. Seventeen
  caption statements now have authorities, read from the same files the tables
  are built from: S5's NSCLC median, its 98 and 15 evaluable sites and the
  post-ComBat median; S6's and S10's site counts, patient counts and both
  medians; S7's two EMT gains (ordinary and adjusted R-squared); S2's Spearman.
  Each has a self-test row; all 18 new rows (with F3.8's) were run in isolation
  and caught. The checker reads 872 claims after this, up from 817 at session
  49's close.

### Session 51, 2026-09-17: two claims whose scope is wider than their data

The ledger was closed at session 50's handoff; these came from the cold reads
the handoff asked for. Neither changes a number. Both change what a sentence
claims, so the wording is the author's (Part B B-22, B-30), and the measured
facts are recorded here. *[2026-09-23: both answered by the author and applied
that day — the abstract carries the two AI sentences with cut 4a, and B-30's
option (a), "the median curated set" (or "signature"), is in the four live
summaries. The superseded structured abstract folded in the manuscript keeps the
old wording as a record.]*

- **The curated-reliability clause states a median as a floor.** The
  conference abstract, the journal abstract, the cover letter and the NextGen
  abstract all say the curated signatures stayed above 0.97 once the global
  axis was removed. Read from `immune_excess.csv` (`alpha_observed`), that is
  true of the median curated signature: 0.976 pan-cancer, 0.977 NSCLC. It is not
  true of the curated signatures as a set. In each cohort 7 of the 16 fall below
  0.97. The lowest are angiogenesis (0.832 pan-cancer, 0.862 NSCLC) and TGF-β
  (0.846 pan-cancer, 0.903 NSCLC); coagulation, hypoxia, IL2/STAT5,
  IL6/JAK/STAT3 and complement sit at 0.927–0.959 pan-cancer and 0.937–0.966
  in NSCLC. The random-set figures in the same
  sentence are medians (Results says so). The Results paragraph and poster
  panel 4 scope the curated figure as a median and are correct; the four
  summary sentences do not. Scoping it costs 7 to 13 characters in the
  conference abstract, measured with `08_count_abstract.py` on scratch copies:
  "the median curated set" +8, "curated sets' median" +7, both clauses made
  explicit medians +11 or +13. With the compact AI pair and cut 4a (below) the
  abstract stays inside the limit with any of them. The archived four-section
  abstract in the manuscript's `<details>` block says the same and is left as
  written.
- **"External validation is pre-registered" in a conclusion.** The last
  sentence of the conference abstract, and the poster's conclusion line, say
  the analyses are exploratory and that external validation is pre-registered.
  The protocol did register a CPTAC validation, but the arm was dropped on
  2026-09-06 because CPTAC-3 has no tissue source site (Limitation 10), so in a
  conclusion the clause reads as a validation to come. Limitation 1 states it
  correctly ("pre-registered but not run; see item 10 for why"). Dropping the
  clause and keeping "Analyses are exploratory." saves 35 characters (cut 4a in
  `submission/ABSTRACT-PORTAL-DRY-RUN.md` §4). That also lets the compact AI
  pair fit without cutting EMT.

Also found and repaired in the same reads, none of them a number:

- **Pan-TCGA was still called pre-registered in three records** (this file's
  A5 block and two rows of `pipeline/results/README.md`) after session 50
  settled that it is "frozen", not registered. Each now carries a dated
  bracketed correction pointing at Appendix A, item 3; the records are not
  rewritten. `15-CHECKLIST.md`'s A9 item and its short version were corrected
  the same way.
- **Two sentence tails stranded by insertions.** The `src/` sentence restored in
  the F6.7 section above, and five wall-clock figures of the 2026-09-05
  five-partition run that `15-CHECKLIST.md` had shown under the 24-partition
  item since session 49 inserted that item's heading. The first shape (flush-left
  text under an indented bullet) is now caught by
  `test_no_list_item_carries_a_stranded_flush_left_line`, with zero noise
  measured first. The second shape is not mechanically detectable and is
  named in the handoff as a class.
- **Coverage.** The checker now reads every figure the abstract counter prints
  (title characters, length including spaces, body words, and the dry run's
  form of the headroom): four claims, 884 in all, each with a live self-test
  row. The dry run's length-including-spaces figure had been stale since the
  counter's closing-rule fix; the self-test re-plants that exact value and
  catches it.

## B. What is genuinely finished

- The primary estimand, in two cohorts, with the four defects from the code audit
  fixed and regression-tested.
- The outcome arm, stratified by cancer type, with MDEs attached.
- Negative controls: raw site AUROC, the ComBat estimability argument, covariate
  baselines, purity decomposition.
- The ancestry supplement, including the identifiability result.
- 203 tests, deterministic seeds, exact version lock, one-command reproduction.

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
| A8 CPTAC | ~$50, 6 days, plus a slide-type check | **CLOSED 2026-09-06 — NOT RUNNABLE AT ANY PRICE.** *(This cell read "deferred, stated" until 2026-09-20; it was the last row in this table still describing A8 as a scheduling choice. The cost and the slide-type check in the middle column are what the arm WOULD have needed and are kept as the record. The arm was dropped because no repository serving CPTAC-3 records a tissue source site, and the index is defined under preserved-site cross-validation — `17-CPTAC-PLAN.md`'s findings table, `09-PAPER-DRAFT.md` Limitations 10. Session 49 searched PubMed and Europe PMC for any report of the missing variable and found none, ledger E22.)* |
| A9 cross-platform reproducibility | mechanism found; fix is one keyword | **CLOSED 2026-09-05 — `splits.py:177`, a non-stable sort. Fixed with `kind="stable"`; shipped and stable fold hashes now both `9f31e33c2caaf1e3` and 0 of 944 patients move. Frozen results deliberately NOT re-run (option B). The "alphas flip" and "amplification" framings were both wrong: the alpha flips are a consequence of the differing partition, and the ridge solve attenuates (cond 18–77, amplification ~0.25).** |
| A10 ssGSEA tie convention | audit + one keyword | **CLOSED 2026-09-05 — `signatures.py:212`, a second instance of A9's defect class and a larger one: 944/944 samples tied, 60.6% of a score sd. Confirmed platform-dependent on two machines (HPC4 job 121781): quicksort depth hash `0bca973cfa26da89` vs `b98fb11e701e8e98`, stable hash `096172c056cc7a34` on both. Fixed with `kind="stable"`. Not on the frozen primary path (`scorer="mean_z"`); it did reach the A7 ssGSEA arms.** |
| A10 sensitivity run | 2 h 45 m of compute | **CLOSED 2026-09-06 — RE-RUN, not argued.** `results/scorer_sensitivity_ssgsea_stable/`, PID 6503. Uncorrected ssGSEA **0.17377773241035724 [0.13325771984752488, 0.21018161815588407]** against the pre-fix **0.1799529379910222 [0.1391910971170491, 0.2159114981000094]**; shift **0.00617520558066495**, 16.1% of the pre-fix interval's own half-width; 15/16 either way; registered estimand still undefined. **The defect does not reach anything reported.** |
| A11 cohort mismatch | comparison, then a re-run | **CLOSED 2026-09-06 and the FIX's own defect closed with it.** `--cohort-scope` shipped 2026-09-05 with a `KeyError: 'patient_id'` that killed A4 and A4b while the script still exited 0. Cause: `DataFrame.join(how="inner")` keeps the left index's NAME only when the two indexes are IDENTICAL (measured, pandas 2.2.3). Invisible because the fix was verified on pan-cancer, where the patient filter is a measured no-op and the KeyError is upstream of the filter. Fixed with `join_scores_to_meta`, a two-directional test, `stage_failures` in the JSON and a non-zero exit. `nsclc_v3_stablesort/science_gaps.json` now holds NSCLC numbers: 944 of 7,168 patients, 31 sites, median R² 0.14533984970354008 against pan-cancer's 0.44379065488261865. **8 min 49 s, against ~2.5 h for the pan-cancer-scoped run.** |
| A5 correction: the pan-cancer "refusals" | diagnosis, one HPC4 job | **CLOSED 2026-09-16 — LAPACK, not geometry.** Seeds 6, 11, 12 and 20 fail in the covariate baseline's gesdd SVD and seed 15 in the decomposition's OLS; the ISI reads neither. `--isi-only` removes both and is bitwise equal on NSCLC (job 129666) and on pan-cancer seed 19 (array 129667). |
| A12 random-patient split is platform-dependent | a library sort | **PINNED 2026-09-24** (the author chose the memo's option b; `splits.stable_group_kfold`; see "Session 60: A12 pinned"). *Until then:* **OPEN — a decision for the author.** The ISI is unaffected. The login and compute nodes give the same Linux tie order (job 129875). A macOS measurement of what pinning it would change is in this session's record (D1). |
| A13 nothing writes `data/interim/` | a builder script | **CLOSED 2026-09-16 — `scripts/00_build_interim.py`.** All six rebuilt from raw: five byte-identical on both platforms, `cohort_nsclc.parquet` byte-identical on macOS and within 1.1e-15 (its mean-z scores) on Linux (job 129998). |
| A14 the shared-mask fast path | one HPC4 job | **CLOSED 2026-09-16 — bitwise equal on real folds (job 129737); NSCLC never takes it.** |
| A15 pandas' CSV parser | one keyword per read | **CLOSED for the supplement 2026-09-16 — exact parsing, tables identical on both platforms except S3's derived p/q, now at 12 significant digits.** Other reads recorded, none on a headline path. |
| A16 expression is not always the tumour | a sensitivity run per cohort | **MEASURED 2026-09-17, in extent (51 of 944 NSCLC, 618 of 7,168 pan-TCGA) and in effect: NSCLC 0.3448 [0.2731, 0.4079] (+0.0526, job 130002), pan-cancer 0.2937 [0.2701, 0.3149] (−0.0031, job 130155).** Reported in Results beside the frozen values. Replace or keep beside: **DECIDED 2026-09-24 — kept beside, on the author's instruction to take the most rigorous course.** The frozen values are the registered configuration run on the inputs as built; the tumour-only profiles are a correction found after the results were known, and it RAISES the NSCLC index by about seven partition sds, so promoting it to primary would be a post hoc choice that improves the headline — the forking-paths move the pre-registration exists to rule out. Both stay reported, and A16's extent and effect are checked claims. |

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
disease is accounted for, and is not technical. **[WITHDRAWN 2026-09-17, ledger
E17: "is not technical" rested on the plate row above, and that row is zero by
construction. The `plate` field is empty in all 10,170 input rows, so the plate
term never entered the design. See "E14, E15, E17" below.]** That strengthens the image-side
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

**Scope, recorded 2026-09-17 (ledger E14).** This table is ONE signature, the
first column script 11 happens to see (A4b calibrates `cols[0]`), not a
summary of the 16. The manuscript quoted its 4.1% beside two 16-signature
medians. Rerun for every signature (HPC4 array 130200,
`pipeline/results/control_c_calibration/`): pan-cancer median **4.7%** (range
2.8–6.0%, 6 of 16 at p ≤ 0.05), NSCLC median **7.0%** (range 1.8–10.5%, 15 of
16). Both reruns reproduce this probe bit for bit. The "~4%" reading survives
as a median; the manuscript now gives both numbers with their scope. The
"agreement with 5.5%" argument still holds for pan-cancer.

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

**Correction, 2026-09-16 (session 48): row 5 was not a partition refusal.** See
"A5 correction" under A12 below. The index does not use the regression that
failed, so split_seed=5 was a computable partition lost to a LAPACK driver.

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
reproduction attempt. *[Corrected 2026-09-17: the pan-cancer result is the
frozen one, not a registered one. The tagged protocol never mentions a
pan-cancer cohort, and `05-PRE-REGISTRATION.md` Appendix A, item 3, records
its addition the day after filing.]* The same holds for NSCLC, whose `split_seed=0` is
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
| partition component | 1.3630 | 3.0171 (sqrt of the site counts, 619 vs 68) | 0.452 |

The bootstrap component is within 4% of the sqrt(n) prediction. The partition
component falls at about **half** the predicted rate on either candidate scaling
— neither the patient count nor the site count explains it. That is a real
asymmetry in the point estimates, and it has **no simple explanation** from
these data; it is not, on this evidence, a demonstrated invariance either.

**Correction, 2026-09-07 (later the same day) — the right numbers, for a
reason worth stating carefully.** The site-count row above originally read
`2.5527 (sqrt of the site counts, 202 vs 31)` with an observed/predicted of
`0.534`. Those two counts are **real and correctly transcribed**: `202` and
`31` are printed by `11_close_science_gaps.py` as `16 signatures, N sites` for
the pan-cancer-scoped and NSCLC-scoped runs respectively (see
`results/science_gaps_stablesort.log` and `..._run3.log`). The error is not
fabrication — it is that they are **the wrong quantity for this argument**.

`202` and `31` are Control C's *usable* site counts: the sites the label-side
variance regression retains after its own filtering. Partition variance does
not arise there. It arises from **assigning sites to folds**, and that
assignment ranges over every tissue source site in the cohort, filtered or not.
The counts that govern it are therefore the full TCGA tissue-source-site counts,
which are **68** for NSCLC and **619** for pan-TCGA. Both were counted directly
from the barcodes in each run's `predictions.csv.gz` and both are corroborated
independently inside the library: `decomposition.py`'s docstring records "641
distinct TSS codes over 8,535 patients" for the full parquet before cohort
filtering, and "n=956 with 68" for NSCLC, and Limitation 8 says "51 of the 68
NSCLC sites". So the prediction is sqrt(619/68) = 3.0171 and observed/predicted
is 0.452.

**The verdict is unchanged, and slightly strengthened:** site count explains the
partition component *less* well than patient count does, not more, so neither
candidate scaling rescues it. Recorded rather than silently edited, because the
failure mode here is a subtle one worth naming — a number can be real,
correctly copied, and still wrong, because it answers a different question than
the one being asked. A first pass at this correction asserted the counts "appear
nowhere else in the pipeline"; that assertion came from grepping only `.md`,
`.py` and `.html` and not the logs, and it was itself wrong. Both errors have
the same shape.

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

#### A5, 2026-09-17: 24 partitions per cohort (ledger B3, B8, B11)

Both sweeps now have 24 partitions (`split_seed` 0–23) under the pinned
ordering: NSCLC from HPC4 job 129201, pan-cancer merged from four HPC4 arrays
(`pipeline/results/partition_variance_pancancer_stable24/summary.json` records
each seed's mode and history). All 48 were computed. The three seeds computed
more than once agree bit for bit, and pan-cancer seed 0 equals the
corrected-ordering primary bit for bit. `27_partition_sweep.py derive` /
`compare` produced every number below.

| cohort | k | partition sd | 95% CI for the true sd | width | share of total variance [chi-square propagated] | interval understated by |
|---|---|---|---|---|---|---|
| NSCLC | 24 | 0.007658 | [0.005952, 0.010743] | 1.80x | 6.37% [3.95%, 11.80%] | 3.34% |
| pan-cancer | 24 | 0.007226 | [0.005616, 0.010136] | 1.80x | 29.86% [20.46%, 45.59%] | 19.41% |

| quantity | NSCLC/pan-cancer ratio | 95% interval | sqrt(patients) predicts | sqrt(sites) predicts |
|---|---|---|---|---|
| bootstrap component | 2.6520 | — | 2.7556 (observed/predicted 0.962) | — |
| partition component | 1.0599 | [0.6971, 1.6115] (F, 23 and 23 df) | 2.7556 — outside | 3.0171 — outside |

The 6-and-5-partition point estimates above (ratio 1.3630, "about half the
predicted rate") are superseded as estimates. They stay as the record of what
k = 5 and 6 showed. With k = 24 the partition component does not fall with
cohort size at all: its ratio's interval contains 1 and excludes both
candidate scalings. The between-partition sd is about 0.007 to 0.008 in both
cohorts. The shares differ because the bootstrap component shrinks and the
partition component does not. The "floor set by the design" reading that
2026-09-07 declined to build on is now what the data show, over these two
cohorts. Two cohorts cannot say how the floor depends on the design (site
count, site-size distribution, K), and the Results paragraph says so.

**B11, the promotion question re-asked with data.** The two sd intervals still
overlap, which is the expected result if the sds are equal, and the two share
intervals no longer do. So the data now support the contrast that was
declined on 2026-09-07. **The decline stands anyway, on the prior-art
ground:** partition-choice variance in grouped cross-validation is known
(Bengio and Grandvalet 2004, reference 13; Howard et al. 2021, reference 7).
The finding is reported in Results, beside the frozen values, not promoted to a
contribution. The narrow claim is unchanged in kind and larger in size: the
reported interval understates total uncertainty by 3.3% in NSCLC and 19.4%
pan-cancer at k = 24.

**Ledger E23.** Beyond 24 partitions? Not needed. The question E23 set was
whether B8's asymmetry is still undecided at 24, and it is decided: the ratio's
interval excludes both scalings.

### Session 52, 2026-09-20: two submission documents read cold, and a number that was right when written

Ledger items H13–H17. Session 51's handoff named `submission/SUBMISSION-CHECKLIST.md`
and `submission/COVER-LETTER.md` as the two documents nobody had read cold. Both
were read in full. The reads found one defect class in six places, and one number
whose drift nothing could have caught because nothing was checking it.

**H13. The cover letter's body length was right when written and 30% stale when
read.** The letter states a body length and a fallback — delete paragraph 4 and
*N* remain — and the author uses both at submission to decide whether it sets on
one page. Measured by the letter's own stated method at `713d0be` (2026-09-02),
the commit that wrote them: **609 words as pasted**, so the stated 614 was right
to within the separator lines, and 614 − 104 = 510 was exact arithmetic. Measured
at this session's inherited commit: **794**. The letter grew by about thirty per
cent while the number stood still, because **every dated correction note since
has been appended INSIDE the measured span** — including session 50's and this
session's own. Paragraph 8 alone is now 265 words and holds most of the drift.
The consequence was not cosmetic: the letter told the author it would set at one
page "but only just", and that cutting paragraph 4 would bring it to 510. At the
true length it will not set on one page, and paragraph 4 alone leaves 697.

The count, the paragraph-4 length and **the subtraction between them** are now
checked claims re-derived on every run (`_cover_letter_words`), each with its own
self-test row, each proven to fail on an injected drift. The subtraction gets its
own row deliberately: that is the exact shape the abstract-headroom defect took
in September — both counts right, the difference wrong — and here the difference
is the number the author acts on.

**A NEW INSTANCE OF AN OLD CLASS, AND THE REASON IT SURVIVED.** This is the
CONSUMER axis in a form the project had not met: a number that is correct at
birth and is invalidated by the *maintenance of its own document*. Every
correction note this project writes is an improvement; appended inside a measured
span, each one silently falsified a count. No guard could see it because no guard
was pointed at it — 614 and 510 were matched by a generic pattern and never
re-derived. Coverage is not verification.

**H14. `_snapshot_size` refused the tree, and it was right.** The first draft of
the new authority anchored the letter's closing boundary on the corresponding
author's email address, as a literal in `19_check_numbers.py`. That script is
published in the public snapshot, so the deposit's leak audit refused the staged
tree — one occurrence, correctly located, exit 1 — and `--strict` then failed
with `NO AUTHORITY snapshot_files`, which is session 50's F3.4 guard doing its
job downstream of it. Re-anchored on the first `orcid.org/` line after `[Date]`,
which identifies the same line without naming the address; the audit returns to
CLEAN and the measured span is byte-for-byte the same. **Recorded as a failure of
this session, not as a feature:** two guards fired in series and neither was
consulted by choice.

**H15. `--strict` exit 0 read through a pipe is not `--strict` exit 0.** Two
readings of the checker's status in this session came back green from a pipeline
whose last command was `grep`, so the 0 was grep's. The real status at that
moment was 1. zsh's `$?` after a pipeline is the LAST command's; `$pipestatus[1]`
is the checker's. §8 has carried this warning since session 48 and it still cost
two wrong readings here. Redirect to a file and test `$?`, or read
`$pipestatus[1]` — never `$?` after a pipe.

**H16. CPTAC was still on sale in three documents.** Session 51 removed a stale
CPTAC go/no-go from `16-YOUR-TASKS.md` §5. Fixing one instance says nothing about
the class, and the class had three more members, found by sweeping every scanned
document for CPTAC beside cost, schedule or decision language:

- `submission/SUBMISSION-CHECKLIST.md` §12 — the most expensive copy, because it
  quoted a price and a timeline: "CPTAC Part B is pre-registered, the FFPE gate
  is resolved, and the cost is ~$27–54 and ~6 days … only CPTAC makes it
  *competitive* at the flagship." Corrected in place, dated.
- `12-MENTOR-OUTREACH.md` — **outward-facing, and the worst of the three.** An
  email template the author would SEND told a prospective mentor that "CPTAC lung
  is budgeted at roughly $50 and six days of GPU time; the open question is
  whether its qualification-workflow slides are FFPE or frozen". The FFPE gate was
  never the blocker and was resolved; the blocker is the missing tissue source
  site. A lung or proteomics reviewer may well know CPTAC's metadata, so the old
  wording proposed an analysis that demonstrably cannot be run. Replaced with the
  closed-negative framing, which is also the stronger pitch.
- `14-SCIENCE-AUDIT.md`'s own A1–A16 verdict table — the A8 row read "deferred,
  stated", the last cell in that table still calling A8 a scheduling choice.

The same sweep found the same wording in `00-PROJECT-BRIEF.md`,
`03-WEEK1-VERIFICATION.md` and `08-READINESS-ASSESSMENT.md`. The first is exempt
as a historical record written 2026-08-11; the other two are dated planning
snapshots that read as history. They are left as written, deliberately, and named
here so a successor does not re-find them as new.

**H17. One outward-facing "timestamped" corrected.** The same mentor template said
the work was "timestamped 18 August 2026". The project has no third-party
timestamp — that is Part B item B-4, still open — so the fixing date is
self-attested. Changed to "protocol fixed 18 August 2026, before the reported
results were computed", which is what `05-PRE-REGISTRATION.md` Appendix A
supports.

**What the cold reads did NOT find.** The cover letter's Mann–Whitney *P* = 0.031
is correct (session 50 fixed it from 0.016). Its "was pre-registered … and it
cannot be run as registered" paragraph is accurate and scoped. Every number in
the letter still traces to the manuscript or the pre-registration. §5 of the
checklist is accurate and dated in its body — only its heading still advertised
the gap as open, and that is now fixed.

**H18. The heading-above-a-tail shape: swept, and not found.** Session 51 named
this shape (a heading inserted above a paragraph adopts the paragraph's tail),
repaired the one known instance in `15-CHECKLIST.md`, and recorded that no test
can see it. It is still not mechanically decidable, but it is mechanically
*narrowable*: the shape leaves a heading whose first following line reads as a
continuation. Scanning every heading in every scanned document for a next line
that begins lowercase or with a connective surfaced **18 candidates out of 408
headings across 28 documents**, and all 18 were read. Every one is either a
lowercase proper noun starting a sentence (`sklearn`, `bioRxiv`) or an ordinary
forward reference ("This is the crux…"). **No stranded tail of this shape
exists in the tree today.**

That negative is only worth the floor and the control behind it. The floor is
408 headings, not a handful; the control is a deliberately planted instance —
a heading inserted between the first sentence of this file's own B11 paragraph
and its tail — which the scan caught, raising the count from 18 to 19 and naming
the planted line. Nothing was written to the repository to do it; the planted
copy existed only in memory. The scan is recorded here rather than shipped as a
test: at 18 false positives per run it would be a guard nobody could read, and
the project's rule is to measure friction before adding a gate.

### Session 53, 2026-09-20: the self-counting class swept, and one hypothesis that did not survive its own check

Session 52 named a defect class it had found once: **a number that is correct
at birth and is invalidated by the maintenance of its own document** (S52-1,
S52-2). It fixed the instance — the cover letter's body length — and recorded
the axis. The project's standing rule is that fixing one instance says nothing
about the class, so this session swept for the class: every scanned document
that states a count **of itself**.

**The sweep found two more members, and they are not the same defect.**

**(a) `16-YOUR-TASKS.md` still described the cover letter as "one page, 614
words".** This is a genuine second instance. Session 52 corrected the letter's
own statement of its length and made it a checked claim; this copy, in the
document the *author* reads to decide what is done, survived untouched. It sits
under a `**Finished 2026-09-02:**` heading and was accurate on that date, so
under the project's convention it keeps its wording and gains a dated
correction rather than a rewrite. The consequence, had it stood: the author's
task list says the cover letter is finished and fits one page, when it is 801
words and does not, and the fallback of deleting paragraph 4 leaves 697 rather
than the 510 the old arithmetic promised. Corrected in place, pointing at B-24.

**(b) `09-PAPER-DRAFT.md` states that its abstract "measures 234 + 28 = 262
words" and is therefore "inside the observed ceiling" of 131–264. This looked
like the same drift and IS NOT.** The hypothesis was that the body had grown
since the count was written. It had not. `git log -S` located `509be54`
(2026-09-02), the commit that wrote the numbers; measured at that commit the
body is **238** by a naive whitespace split, and it is **238 today** — the
abstract text has not changed at all. The stated 234 is reproduced exactly, and
28 with it, by one specific measure: whitespace-delimited tokens **containing a
word character**, which drops four tokens — three spaced em-dashes and the `=`
of `*P* = 0.031`. So the number was right when written, is right now, and the
document had simply never said which of two defensible counts it meant.

That distinction matters enough to record as a near-miss. Six counting variants
were tried before concluding anything; had the check stopped at the naive split,
this session would have reported a 4-word drift that does not exist and
"corrected" a correct number. **A number that disagrees with your recount is not
yet a defect — first ask whether your recount is the document's method.** The
cover letter states its method; the manuscript did not, and that was the real
gap.

**What was actually wrong, and what it costs.** None of 234, 28 or 262 was bound
to any authority — `grep` over `19_check_numbers.py` returned nothing for all
three. They were literals matched by a generic pattern and never re-derived,
which is session 52's other rule, COVERAGE IS NOT VERIFICATION, in its second
instance. And the exposure is immediate rather than theoretical: **B-30 edits a
line inside the measured span.** `09-PAPER-DRAFT.md`'s abstract carries the
"curated signatures held above 0.97" clause, so whichever wording the author
chooses changes the body count, the sum, and possibly the ceiling claim — and
nothing would have said so.

Measured, by applying each candidate to the real span:

| Wording in the manuscript abstract | body | + Significance | total | vs the 131–264 ceiling |
|---|---|---|---|---|
| current ("curated signatures held above 0.97") | 234 | 28 | **262** | 2 words of margin |
| B-30 (a) "the median curated set held above 0.97" | 236 | 28 | **264** | **exactly on the ceiling, zero margin** |
| B-30 (b) "curated sets' median held above 0.97" | 235 | 28 | **263** | 1 word of margin |

So B-30 is not cost-free in the manuscript even though it is nearly free in the
conference abstract: option (a) spends the entire remaining margin against the
observed ceiling, option (b) spends half of it. Neither breaks the claim, and
the ceiling is an *observed* range across 40 articles rather than an enforced
cap — the official cap remains unverified because `aacrjournals.org` refuses
automated requests — so this is a consideration for the author, not a blocker.
It is recorded because the author is choosing that wording now and had no way to
see this cost. Under the naive token count the total is 266 today and already
past 264, which is precisely why the paragraph now states its measure.

**Closed by:** `_paper_abstract_words()` re-derives all three on every run, with
structural anchors that return nothing rather than measure the wrong span if
either end moves. Three patterns, three self-test rows, and all three proven to
fail against the real file — body +5, Significance −3, and the **sum pushed +4,
deliberately in the direction that carries it past the ceiling**, so a MISSED on
that row would mean the ceiling claim had gone unguarded. The file was restored
byte-identically after each (sha256 compared, not assumed). The sum has its own
row for the reason session 52 gave the cover letter's subtraction one: September's
abstract-headroom defect was both parts right and the arithmetic between them
wrong, and here the arithmetic is what the conclusion rests on.

### Session 53, 2026-09-20: an outward-facing ask that cannot be granted, and a hard deadline three days late in two documents

Session 52 swept the CPTAC class into `12-MENTOR-OUTREACH.md` and called that
instance the worst of three because the document is **outward-facing** — it
holds emails the author would actually send to named senior researchers. It
corrected the CPTAC paragraph and did not re-read the rest of the file. The
class rule applies to session 52's own fix: **fixing one instance says nothing
about the others**, and this document had a second, larger one.

**The document's headline instruction was "Ask A first", and Ask A is
impossible.** Ask A is a NextGen Stars recommendation letter, deadline **Wed
2026-09-23 — three days from this session**. Two independent blocks, either one
sufficient:

1. `13-NEXTGEN-EXTENDED-ABSTRACT.md` records that NextGen requires **AACR
   membership in good standing before the application can be opened**.
2. Membership is **not started** (Part B item B-1) and `08-READINESS-ASSESSMENT.md`
   puts processing at roughly 30 days. Three days is not thirty.

Two other documents already record the path as closed — `16-YOUR-TASKS.md` ("closed
to you on eligibility", and the September deadline pressure "is gone") and
`08-READINESS-ASSESSMENT.md` (struck through, "<5%, effectively closed"). So the
project knew. The knowledge simply never reached the file that would have acted
on it, and that file was the only one pointed outward.

**Why this outranks a wasted email.** Ask A and Ask B target the *same people* —
the document says so itself, and the whole strategy is that A warms up B. Ask B,
a senior author on the paper, has **no deadline** and is the thing actually worth
having. Opening with a request that cannot be granted spends the single first
impression B depends on, in front of precisely the reader most likely to notice:
the email names the recipient's own 2021 method and asks them to read a
manuscript in three days for an application that cannot be opened. Email A is
now marked DO NOT SEND, the table row struck through, and the reader is sent
straight to Email B. The template body is kept — its two-finding paragraph is
the best short statement of the result in the repository and is raw material
for B.

**Separately: the membership cutoff is stated three days LATE in two documents.**
`~13 Oct` was a 2026-cycle projection, superseded when the 2027 Call was read end
to end on 2026-09-17 — the Call requires the application "at least one month
prior" to the 10 Nov deadline, so the cutoff is **Sat 2026-10-10**. That
correction reached `16-YOUR-TASKS.md`, `11-SUBMISSION-PACK.md` and
`15-CHECKLIST.md`, each with a dated note. It did **not** reach
`08-READINESS-ASSESSMENT.md` (five places) or `10-PROJECT-REFRAME.md` (two).
Both now carry dated corrections.

That decision needs its reasoning on the record, because session 52 deliberately
left this same file's CPTAC wording alone and §3 of the handoff passed the
question forward. **The line drawn here is not "planning documents are exempt"
but "what is the reader going to do with it".** A stale CPTAC cost is harmless
history: acting on it is impossible, because the arm is dropped. A stale
membership deadline is not history at all — it is live, it is 20 days out, it is
the one hard gate in the project (no membership, no abstract), and it is wrong in
the **permissive** direction, which is the only direction that can cost anything.
A reader who trusts it believes they have until the 13th. The CPTAC wording in
those two files is therefore still left as session 52 left it, deliberately, and
only the dates are flagged.

**Not a defect, recorded so it is not re-found:** the 2027 Call's 23 Sep NextGen
deadline passes in three days, after which `13-NEXTGEN-EXTENDED-ABSTRACT.md` and
`01-LANDSCAPE-AUDIT.md` will state a future date that is past. Both read as dated
records of a closed path and are left as written.

### Session 53, 2026-09-20: the poster's drawn result was the least-checked thing in the project

`§3` of the session-52 handoff named the two biggest remaining coverage gaps —
`poster/poster.html` at 86 unchecked literals and `pipeline/README.md` at 138 —
and called them the only ones left worth naming. Opening the poster's list
settles which of the 86 matter: **the outcome panel's entire drawn result was
unchecked.** Five per-signature excesses, the multiplicity margin beside them,
and nothing binding any of them to an artefact.

That is the worst place in the project for an unchecked number. A poster is read
by people standing in front of it who cannot re-run anything, the numbers are
printed at A0 and cannot be corrected once the board is up, and these six are
exactly the figures a reviewer challenges — they are the paper's negative
result, the claim that the prognostic signal is proliferative rather than
immune.

**Bound to `pancancer_v3/outcome_arm.csv`, and the adjustment is named.** Each
signature has two rows there, `axis_residualised` and `raw`. The poster prints
the residualised arm. Taking the first matching row instead would not fail
loudly — under `raw` the *same five* signatures still beat their null, so the
shape of the claim survives, but the **order changes**: EMT leads at 0.0623 and
G2M falls to 0.0520, against 0.0687 and 0.0382 residualised. A row-order-
dependent read would therefore check one arm's numbers against the other arm's
printed values and could agree by luck on some of them. The selector names both
the signature and the adjustment, and asserts exactly one row matches.

**The manuscript's copy is bound too, and it needed WRAP.** `09-PAPER-DRAFT.md`
states the same five values, and they agree — verified before anything was
bound, since a poster and a paper disagreeing would have been the real finding.
They wrap: `**E2F targets**` ends one line and `(0.062)` opens the next, which
is precisely the shape session 50's F3.8 found a line-scoped pattern cannot see.
All five manuscript patterns are WRAP-scoped and use `\s+` between tokens, so
re-flowing the paragraph cannot blind them.

**Both copies got self-test rows, and that is not redundant.** A row that only
rewrites the poster would still report CAUGHT, because both documents read the
same authority and the poster's mismatch alone fails the run — and the
manuscript pattern would never once be proven to fire. A green row that proves
nothing is the vacuity defect this file has shipped before. Eleven rows were
added in total and each was verified to match **exactly once** across the
scanned tree before being trusted, then every one was injected against the real
file and caught, with the file restored byte-identically each time (sha256
compared, not assumed).

Coverage moved 861 → 872 checked and 928 → 939 claims at this commit; the poster
is no longer the least-covered document by the measure that matters, which is
not the percentage but whether the numbers a reader acts on are bound to
anything. `pipeline/README.md`'s 138 remain, and are mostly prose about the
pipeline rather than results a reviewer would challenge — a successor should
spend effort there only after checking that the same is still true.

### Session 53, 2026-09-20: the eight RRIDs resolved, and the blocker was the request, not the registry

F4.1/F12.8 had been BLOCKED across several sessions: "SciCrunch returns 403 to
scripts, to WebFetch and in the browser pane", so none of the manuscript's eight
RRIDs had ever been verified to resolve. That is now closed, and the diagnosis is
worth recording because the item was blocked on the wrong thing.

**The registry was never refusing us. The HTML pages are what 403s.** Requesting
the JSON resolver endpoint — `https://scicrunch.org/resolver/<RRID>.json` — with
an ordinary browser `User-Agent` returns **HTTP 200 for every one of the nine
identifiers tried**, including the two scikit-learn records. Previous sessions
recorded a service-level block and stopped; it was a content-negotiation and
user-agent problem, one flag wide.

**All eight resolve, and all eight are correct:**

| RRID | resolves to | cited as |
|---|---|---|
| SCR_008394 | Python Programming Language | Python 3.13.9 |
| SCR_008633 | NumPy | NumPy 2.1.3 |
| SCR_018214 | Pandas | pandas 2.2.3 |
| SCR_008058 | SciPy | SciPy 1.15.3 |
| SCR_016074 | statsmodel | statsmodels 0.14.4 |
| SCR_008624 | Matplotlib | matplotlib 3.10.0 |
| SCR_016863 | Molecular Signatures Database | MSigDB v2024.1.Hs |
| SCR_002577 | scikit-learn | scikit-learn 1.6.1 — **changed, see below** |

**The scikit-learn reconciliation, which F4.1 also asked for.** Both records are
live and neither is marked superseded, so this is a choice between duplicates
rather than a correction of a dead identifier. The manuscript cited
**SCR_019053**, whose registry name is "Sklearn", whose description reads
"Software Python package part of nonnegative matrix factorization NMF" — the
signature of a record generated from one paper's methods section — and which
carries "Scikit-learn" only as a synonym. **SCR_002577** is named exactly
`scikit-learn`, its description is the project's own tagline, its abbreviation
field is `scikit-learn`, and it is typed as a software application rather than a
toolkit. The manuscript now cites SCR_002577. This is a one-token change and
trivially reversible if the author prefers the record they originally chose;
nothing else in the paper depends on it.

**Method recorded so the next session does not re-find the 403:** fetch
`scicrunch.org/resolver/<RRID>.json` with a browser User-Agent and read
`hits.hits[0]._source.item.name`. The HTML page at the same path is what refuses
automated clients.

### Session 53, 2026-09-20: Limitation 8 folded, and the cover letter cut as far as mechanics can take it

**B-21 — Limitation 8's length.** Measured rather than assumed, and the carried
premise was slightly wrong: the handoff said "~150 lines against 2–20 for the
others", but the siblings actually run **2–23** lines (Limitation 9 is 23) and
Limitation 8 stood at **183**. A first measurement also reported Limitation 12 at
397 lines, which would have made it the real outlier — that was a scripting
error, not a finding: the Limitations section ends at a `##` heading and the scan
looked for `###`, so it ran to the end of the file. Limitation 12 is 16 lines.
Recorded because it was caught only by asking whether a 397-line limitation was
plausible.

The restructure keeps **every word**. The forensic material — the mechanism, the
two corrected re-runs, the sixteen-significant-figure cross-platform validation,
the scope of the library ordering, and the pipeline-wide ordering audit — is
folded into a `<details>` block titled as Supplementary Note 1, the same device
the manuscript already uses for its superseded structured abstract. The open
limitation is now **31 lines**: the defect, the two platform values, the cause,
the magnitude against the partition standard deviation, the consequence, and
that the stable sort ships. Total length went **up** by three lines, because
nothing was deleted and the summary text was added.

Folding rather than moving to a separate file is deliberate and is the
conservative choice. Every number in Limitation 8 stays in the same document, so
no document-scoped checker pattern breaks and the per-document coverage floor
cannot drop — both verified, 939 claims and 0 vacuous patterns before and after.
It also leaves the author's actual B-21 decision open and makes it mechanical: if
a real supplementary note is wanted at submission, the folded block lifts out
whole. If the author disagrees with the fold at all, it is two lines to revert.

**B-24 — the cover letter.** Three defects fixed outright, because none was a
judgment call. (a) Paragraph 8 carried an editorial aside addressed to the
author, not the editors — the italic "(True as of 2026-09-04: …)" note — which
would have been pasted into the cover-letter field verbatim; removed. (b)
"copyright retained by the author" said *author* singular for a six-author
submission; corrected. (c) The letter presented the cross-platform difference
with no cause, the pre-A9 framing, which reads as an unexplained defect to an
editor; it now names the non-stable sort and says the corrected values were
reproduced from raw files on both platforms.

The length is where mechanics run out. The body went **801 → 751** across two
compression passes over the boilerplate paragraph, and the second pass saved
**three words** — the point at which further squeezing stops being mechanical.
The calibration is the letter's own history: it measured 609 when it claimed to
set at one page "but only just", so the target is about 610 and it is roughly
140 over. Deleting paragraph 4, the designated fallback, leaves 647 and still
leaves ~40 to find. The header enumerates what is left, paragraph by paragraph,
so any further cut means choosing which part of the argument to lose — a
judgment about how to pitch the paper, which is the author's. The header now
states all of this, with the arithmetic, rather than only the word count.

*[Reworded 2026-09-22 by session 55. This paragraph restated the header's
per-paragraph enumeration in the live present tense. Every figure in it was
correct both when written and when re-measured, so nothing here was wrong —
but B-24 is an open decision to DELETE one of those paragraphs, and the moment
the author acts on it this record would have become false and `--strict` would
have gone red pointing at a record of the past. That is the failure session 54
recorded as H37, and the lesson session 53 recorded as H30: where a
cross-document count cannot be made a checked claim, stop asserting it and
point at the thing that re-derives it. The enumeration is now a set of checked
claims in the letter's own header, re-derived on every run.]*

**One guard fired on this session's own work, and it was the right one.** The
rewritten header first closed the paragraph-4 sentence with the word "leaving"
in place of the "and N remain" phrasing that the `cover_letter_body_minus_para4`
pattern session 52 added anchors on, which broke it. (Described rather than
quoted, deliberately, and corrected here after the first draft of this very
paragraph quoted the header's bolded fragment verbatim — which made a historical
record match the live pattern, so a future change to paragraph 4 would have
forced an edit to the record of the past. The coverage diff caught it: one
pattern silently went from reaching one document to two.) The
claim count fell 939 → 938 and `--strict` failed with `VACUOUS`, not with a
mismatch: the number was still correct, the pattern had simply stopped seeing it.
That is the failure mode F3.4 exists for, and it is the second time in two
sessions that a guard has caught the session that wrote it.

### Session 53, 2026-09-20: the last named coverage gap held a deposited results table, not prose

Session 52's handoff named `pipeline/README.md`'s 138 unchecked literals as the
largest remaining coverage gap, and this session's own earlier note predicted
they would turn out to be "mostly prose about the pipeline rather than results a
reviewer would challenge", advising a successor to check before spending effort
there. **That prediction was wrong, and checking is what found it.**

Most of the file's literals are indeed uninteresting — DOI fragments, arXiv
identifiers, journal article numbers, package versions, file sizes and measured
runtimes. But the file also carries a **full headline results table** restating
both cohorts' frozen numbers: the two ISIs with their confidence intervals, the
signature counts, the outcome arm, the NSCLC minimum detectable effect and both
site AUROCs. Twenty-one of the file's 459 literals were checked and **not one of
these was among them.**

That is a worse exposure than the poster's, for one reason: **`pipeline/README.md`
is deposited.** It ships in the public snapshot, and it is the first file anyone
opening the repository reads. A drift there does not merely misstate a number, it
publishes an artefact that contradicts the paper it exists to support — and
`--compare-live` would happily report the deposit in sync while it did so,
because that check compares the deposit against the working tree, not against
the frozen results.

Eleven patterns added, all anchored on the row label or the cohort name and
never on the digits, each bound to the authority the manuscript already uses, so
the README cannot disagree with the paper without failing the run. One of them,
the pan-TCGA cohort size, reaches three documents rather than one. Ten
injections were applied to the real file — one digit of the reported precision
each, which is the shape a transcription slip actually takes — and all ten were
caught, with the file restored byte-identically every time. Six self-test rows
make it permanent rather than a one-off check. Claims 939 → 952, patterns
315 → 326, still 0 matching nothing.

**The lesson is the one this project keeps relearning in new costumes: a guess
about where the risk is, written down confidently, is not a measurement.** The
guess here was written by the same session that then disproved it, one message
apart.

### Session 53, 2026-09-20: a sweep for the session's own wake, and the axis bit the session hunting it

Every change this session made to a number was followed by a sweep for other
documents restating it — the discipline the session had just spent its time
arguing for. It found two live contradictions, and the more instructive one was
self-inflicted.

**(a) The checklist contradicted the manuscript about an RRID.**
`submission/SUBMISSION-CHECKLIST.md` recorded that SciCrunch "gives the proper
citation as `Sklearn (RRID:SCR_019053)`, and that is what the manuscript uses",
and deferred the duplicate question to a future SciScore report. After H27 the
manuscript cited `SCR_002577`, so the checklist asserted something false about
its own manuscript and pointed at a third-party service for an answer already in
hand. Corrected, with the evidence, and the 2026-09-08 identifier list above it
left as the record of what was written that day.

**(b) A correction note written EARLIER IN THIS SESSION was stale within the
hour.** The dated note added to `16-YOUR-TASKS.md` to fix the cover letter's
"one page, 614 words" said the letter measures 801. Two hours later the same
session cut the letter to 751, and the correction note became exactly the kind of
artefact it had been written to fix: a number, correct when written, invalidated
by later work on the thing it described — session 52's CONSUMER axis, committed
by the session that swept for it.

That is worth more than the fix. The rule "a document stating a count of itself
gets a checked claim, not a note" was written for documents counting
*themselves*; this is a document counting *another* document, where the same
failure applies and a checked claim is not available, because the note is prose
about history rather than a live statement. The note now says so explicitly: it
records both figures, names its own staleness as the illustration, and points the
reader at the cover letter's own header, which re-derives the live numbers every
run. **Where a cross-document count cannot be made a checked claim, the durable
fix is to stop asserting it and point at the thing that re-derives it.**

Handoff files 24 and 43–52 also restate these numbers and are left untouched:
they are excluded from prose scans by name and are records of what was believed
on their dates.

### Session 53, 2026-09-20: TODO-PRINT-6 option (a) built, so the print decision can be made by looking

B-5 / F11.1 has stood open since 2026-09-08 with three options and a
measurement: once each figure is scaled into its poster box, the four A0 figures'
7 pt tick labels print at 13.5, 12.2, 10.5 and 7.2 pt, against poster guidance of
about 18 pt for axis labels and an absolute floor of 14. All four are below all
three. Option (a), the cheapest, was already specified in the poster's own
comment — "give the figure script a `--poster` flag that re-renders at ~2x font
sizes for these four figures only, leaving journal figures at 7 pt".

It is now built. `10_make_figures.py --poster [SCALE]` (default 2.0) multiplies
only the six type-size rcParams and writes to `results/figures-poster/`; dpi,
spines, grid and the Type 42 font setting are deliberately untouched, so the only
difference between a poster render and a journal render is glyph size. The
committed figures were verified **byte-identical across a `--poster` run** (md5
before and after), which is the property that makes the flag safe to leave in.

**Building it produced a number the author needs, and it is not the whole fix.**
At SCALE = 2 the four figures print at 27.0, 24.4, 21.0 and **14.4** pt. Three
clear the ~18 pt guidance comfortably; `figure1_isi` clears only the 14 pt floor,
because it is height-constrained in a 422 × 108 mm box and scales 1.03× into it
rather than 2.12×. **So option (a) alone does not bring every figure to guidance
— it brings three of four there and the fourth barely over the floor.** Reaching
18 pt on `figure1_isi` needs option (b), the taller panel, or about SCALE = 2.5
for that figure alone. That is a real constraint on the choice and it was not
visible before the flag existed.

The flag does not choose between (a), (b) and (c); it makes (a) renderable, so
the choice can be made by looking at a proof rather than at arithmetic. The
output directory is gitignored — regenerable renders, and the project does not
stage figure files — and classified in `results/README.md`, because the
classification guard scans the filesystem rather than the index and would
otherwise fail the moment anyone ran it.

*[2026-09-24, session 60 (H87, H88): the author chose the most rigorous option
and both halves are done. The "about SCALE = 2.5" above was an argument, not a
measurement: at any scale the 16 rows of figures 1 and 3 did not fit their A0
boxes at 18 pt without overlapping. What shipped instead: the journal figures
are drawn at print size (no text under 8 pt at 6.75 in; the smallest had been
5.6 pt, below the 7 pt ticks F11.6 measured); `--poster` now draws each poster
figure for its box at 1/2.25 of its size into the committed `poster/figures/`;
the grid's rows went 1.241 → 1.10 over 1, panel 4's pan-cancer partition
paragraph moved to panel 6, and the poster's figure 2 is its panel A. Measured
in a browser: every poster figure prints its smallest text at 18.00–18.01 pt.
`results/figures-poster/` is retired.]*

### Session 53, 2026-09-20: the cover letter's two required counts computed, and B-21 turns out to be about length

`submission/COVER-LETTER.md` opened with "The manuscript runs **[X] words** …
with **[N] figures and no tables**" and had carried those placeholders since it
was written. They are not decorative: the submission pack records that AACR asks
the cover letter to state both. Nothing in the repository had ever computed
either.

Both are now derived and bound. **8,082 words and 4 figures.** The span is the
letter's own definition — Introduction + Results + Discussion, excluding the
cover page, abstract, Materials and Methods, tables and references, with the
back matter (Data and Code Availability, the authorship sections,
Acknowledgments) excluded too. HTML comments and tags are stripped, so the
`<details>` markup contributes nothing, and the token measure is the same one
`_paper_abstract_words()` uses, so the letter's two numbers are commensurable
with the abstract's. The figure count is taken from the manuscript's own
citations, which reach Fig. 1–4; `figure0_schematic` is a poster panel with no
data and is correctly not among them.

**And computing it turned B-21 from a style question into a length question.**
The folded Supplementary Note 1 sits inside the counted span, and it is
**1,593 words — 19.7% of the body**. Lifting it into a real supplementary file,
which is exactly what B-21 asks, takes the manuscript from 8,082 to **6,489**
*[2026-09-23: B-34's clarification added six words to Results — 8,088 and
6,495 then. 2026-09-24, B-21 decided by the author: the note was LIFTED into
`SUPPLEMENTARY-NOTE-1.md` word for word, with a four-line pointer left in
Limitation 8, and the manuscript's body is 6,540 words by the same measure.
Both counts are re-derived by the checker on every run]*.
That is no longer a judgment about how a limitation reads; it is potentially the
difference between a compliant submission and an over-length one, and the author
should decide it with that number in view. `paper_suppnote_words` is computed on
every run so the figure stays live.

**One measurement error, caught and worth recording.** The first attempt located
the folded note by the first `<details>` in the file and reported **310 words**.
That block is the manuscript's *superseded structured abstract*, which has been
folded since 2026-09-02; the Limitation 8 fold is the second one. The error was
caught by noticing that 310 words could not be a 150-line block. The authority
now anchors on the first `<details>` **after the Discussion heading** and says
so in its docstring, because the file has two and will plausibly gain more.

### Session 53, 2026-09-20: the OS changed underneath the project, and nothing moved

The laptop was upgraded from **Darwin 25.5.0 to 27.0.0** (macOS 27, arm64) and
rebooted at 17:49, in the middle of this session. That is not cosmetic for this
project. The frozen headline results were computed on macOS 15 / arm64 *[corrected
2026-09-21: macOS 26.5.2 / Darwin 25.5.0; see the session 54 block]*, the
entire A9 investigation turned on a platform difference of 0.0218 in the second
decimal, and Limitation 8 exists because two platforms disagreed. A major OS
version is exactly the kind of change that could move a library's behaviour.

**Measured immediately after the reboot, before anything else:**

* **Provenance 220/220.** Every one of the 220 frozen artefacts hashes to its
  recorded value on the new OS. Not one byte moved.
* **The stable-sort partition still hashes to `9f31e33c2caaf1e3`**, the value
  measured on macOS arm64 and on HPC4 x86_64, and the shipped partition is
  identical to it.
* **`--strict` passes at 954 claims across 24 files**, so every checked claim
  still agrees with the artefacts it is derived from.
* The pinned toolchain survived intact: Python 3.13.9, numpy 2.1.3, scipy
  1.15.3, pandas 2.2.3, scikit-learn 1.6.1.

So the OS upgrade is a free cross-platform datapoint the project did not plan
and would not have paid for: **the frozen results are stable across a macOS
major version on the same hardware.** Recorded because the opposite result
would have been a serious finding, and because a successor comparing figures to
`_run53b` needs to know the platform string changed between them.

### Session 53, 2026-09-20: a duration measured across a sleeping machine is not a duration

`_run53c` reported the self-test step at **2,253.39 s against its 2,400 s
limit** — a margin of 146.6 s, or 1.065×. Taken at face value that would
overturn three prior measurements (2.5× at 314 rows, 2.1× at 321, 2.05× at 335)
and mean the guard is a few seconds of bad luck away from failing a green tree.
It was not taken at face value, and it should not be quoted.

**The run was not 40 minutes long.** The log was born 11:52:18 and last written
17:27:13 — **5 h 35 m of wall clock for a run pytest measured at 2,430 s.** The
machine slept for roughly 4.9 hours in the middle of it, and the reported
durations, which come from a monotonic clock that does not advance across a
Darwin system sleep, exclude that gap while the work was still interleaved with
whatever the machine did on either side of it.

**Why it slept, and the new form of an old rule.** `caffeinate -i -w <pid>` was
running, and it was correct, and it still failed — because it was a child of the
*shell session* that launched the suite, and that session ended. The project's
rule has been "hold `caffeinate -i -w <pid>` for every detached run" since
session 50. It is now sharper: **the holder must outlive the session, not merely
the command.** The fix is to detach the caffeinate itself, so it reparents to
PID 1 exactly as the run does; `_run53d` was launched that way and both
processes were verified at PPID 1 before anything else happened.

**The margin question is therefore still open, and tonight cannot close it.**
The re-run is also on a busy machine, and for two reasons that were measured
rather than assumed: the foreign THEIA job is back at 322% CPU, and two
`Metadata.framework` processes are reindexing at 87.7% and 77.4% after the OS
upgrade. A successor wanting the true per-row cost and the true margin needs a
genuinely quiet machine, which this one will not be until the reindex settles.
What `_run53d` can still establish — and what matters more — is that the tree is
green on macOS 27.

### Session 53, 2026-09-20: the first suite on macOS 27 went RED, and both failures were real

`_run53d` — the first run of this tree on the upgraded OS — reported **2 failed,
187 passed**. Kept as the record, because both failures were true positives and
one of them is a finding about the manuscript.

**First, the timing question the run existed to answer: the margin is fine.**
The self-test step took **1,318.64 s at 343 rows**, leaving **1,081 s against
the 2,400 s limit, 1.82×**. So `_run53c`'s alarming 2,253.39 s was the sleep
contamination diagnosed above and nothing else. **The limit is adequate and is
NOT being raised**, which was the change under consideration before this
measurement existed. That is the whole value of re-measuring rather than acting
on the first number: the "fix" would have been a permanent loosening of a guard
in response to an artefact.

**Failure 1 — `test_historical_run_rows_state_no_live_test_count`, and the
guard's own docstring predicted it.** A results row written earlier in this
session gave the count in the live form (the number, a space, then the word
"tests") instead of the hyphenated `189-test` convention *[reworded 2026-09-21:
this sentence first QUOTED the phrase, so the checker read the quotation as a
live count, and it went red the moment the suite grew, which is the failure it
describes]*, so it read as a claim about the current suite and would go red the
moment the suite grew. The docstring records that session 41 introduced four
such phrasings plus four more, session 42 fixed all eight and then introduced
two of its own within the hour — "ten instances, two sessions, one convention
everybody agreed with. That is what a guard is for." Session 53 is the third
session to do it, in a row written by the agent that had read that docstring
earlier the same day. The convention is not hard; remembering it while writing
prose about something else is. Fixed.

**Failure 2 — `test_environment_md_matches_a_fresh_render_on_this_platform`,
which is the OS upgrade arriving where it actually bites.** `ENVIRONMENT.md`
recorded `Darwin 25.5.0`; the machine is now `Darwin 27.0.0`. The test's
instruction is to regenerate — correct for its normal case, where the machine
drifted under you — but **wrong to follow blindly here**, because the file was
titled "Exact environment that produced the frozen results" and regenerating it
would have made that title false while looking like a routine fix.

Resolved by making the file say what is true rather than what is convenient: the
title is now "Exact environment, read from the running interpreter", and the
header records that the table describes the machine it was last generated on,
that the two coincided until the 2026-09-20 upgrade, and that the frozen results
were re-verified on the new OS with nothing moving. Renamed rather than quietly
regenerated into a contradiction.

**And it exposed an unverified claim in the manuscript.** Chasing the Darwin
string showed that **nothing in this project has ever recorded the macOS
marketing version** — the generator reads `platform.release()` only. The
manuscript states in three places that the frozen results were produced on
"macOS 15 / arm64", and that string has never been checked by any tooling. On
this machine today `platform.mac_ver()` and `platform.release()` agree at 27.0,
so a Darwin 25.5.0 machine was not macOS 15. What the August 2026 runs actually
ran on **cannot now be established from this repository**, because the earliest
generated copy of `ENVIRONMENT.md` is dated 2026-09-04 and postdates them. The
generator now records `product version` so the question can never be ambiguous
again, and the manuscript's "macOS 15" is flagged for the author rather than
rewritten to a value this session cannot establish — see Part B. *[Answered
2026-09-21, session 54: it CAN be established, from the machine rather than the
repository — the laptop's install record gives macOS 26.5.2 throughout. And the
reasoning above has a flaw worth keeping visible: "`mac_ver()` and `release()`
agree at 27.0" is true of macOS 27 and false of the release that mattered,
where Darwin 25.5.0 is macOS 26.5.2. See the session 54 block.]*

### Session 54, 2026-09-21: an OS name the machine could answer, two rows that never closed, and a poster lead bound to nothing

**B-33 answered from the machine, not the repository.** Session 53 concluded
that the OS behind the frozen results "cannot now be established from this
repository", and that was true: the repository never recorded it. The laptop
did. `softwareupdate --history` and `/Library/Receipts/InstallHistory.plist`,
two independent records, agree: macOS 26.5.1 installed 2026-06-28, macOS 26.5.2
on 2026-07-19, then nothing until macOS 27.0 on 2026-09-20. The Claude Code
transcript that covers 2026-08-11 to 2026-09-02 shows `05_run_nsclc.py
--outdir results/nsclc_v3` launched at 2026-08-19T16:24Z from this laptop's own
working directory. So the frozen results were produced on **macOS 26.5.2
(Darwin 25.5.0) / arm64**, which agrees with the Darwin string `ENVIRONMENT.md`
carried from 2026-09-04. "macOS 15" was wrong, and it had spread: nine live
statements (the manuscript three times, the deposited landing page, the
deposited README's A9 table, `reproduce.sh`'s timing header, the docstring of
`18_bisect_platform.py`, a comment in the test suite and the generator behind
`ENVIRONMENT.md`) now say 26.5.2, and five records (four in this file, one
results row) carry dated brackets instead of being rewritten. Every timing the
live statements vouch for was measured between 2026-09-04 and 2026-09-17,
inside the 26.5.2 window. The sweep was every OS version string in every
tracked text file, not a grep for the one phrase.

**The argument that preceded this had a flaw worth recording.** Session 53
reasoned that because `platform.mac_ver()` and `platform.release()` agree at
27.0 on the upgraded machine, a Darwin 25.5.0 machine was not macOS 15. The
conclusion was right and the premise does not transfer: on the release that
mattered the two numbers differ, Darwin 25.5.0 being macOS 26.5.2. Carried one
step further, the premise would have produced "macOS 25.5", a release that
never existed. An argument is not a measurement, even when it reaches the right
answer.

**Two results rows still read `pending` for runs that had long finished.**
`full_suite_20260920_run53d.log` — session 53's first suite on macOS 27, which
went RED and whose two failures were the reason for the next run — and
`wrapper_20260916_run47e.log`, from session 47. Since session 46 a row is
written before its run with the literal word `pending` and filled from the log
afterwards; nothing checked the second half, and the classification test passes
over both because it asks only whether a name appears. Measured before gating:
41 rows carry the placeholder, 30 of them name logs that exist and say they
finished (a pytest summary, or the wrapper's `SUITE_EXIT` line), and against
the map as committed at `9da6137` the detector flags exactly those two and
nothing else. A first draft that keyed on "row written before the run" flagged
nine older rows that predate the placeholder convention and record their
result as COMPLETE; it was narrowed to the placeholder rather than shipped
noisy. Both rows are filled from their logs, and the check is now a test pair.

**The suite grew by that pair, and a record went red the moment it did.** One
of the ten lines stating the live test count was session 53's own account of a
results row that had stated a count in the live form. It QUOTED the phrase, so
the checker read the quotation as a live claim, and the growth of the suite
made it false — the exact failure that paragraph describes. It now describes
the form instead of quoting it, which is this project's older rule arriving
late, and the claim total fell by one for that reason alone.

**The poster's lead sentence was bound to nothing.** Session 53 bound the
outcome panel and called the poster's drawn result the least-checked thing in
the project. It was not the only unchecked thing on the board. A mutation sweep
moved each literal the coverage report listed as unchecked on the poster by one
step in its last printed place, against the real file in a scratch clone, ran
`--strict`, and restored the file with its sha256 compared every time. The
first pass ran with `--fast`, which skips the live-measured authorities, and so
reported the tissue-source-site count as unguarded when it is not; the second
pass ran without it. What was genuinely unguarded: both cohorts' ISI and 95%
intervals in the lead — the largest type on the board — and every number in
the reliability panel and its small-panel table, the per-signature ranges, the
rotation null, both partition-widened intervals and their percentages, the
image-axis correlations, the NSCLC cohort and event counts, the gene count, and
most of the limitations panel. None of them was wrong. On an A0 board that
cannot be corrected once it is up, unchecked and wrong carry the same risk.
The poster writes HTML entities where the manuscript writes characters, so the
manuscript's patterns never reached it.

Fifty-two claims were added, forty-nine on the poster and three on the
deposited landing page, each anchored on the poster's own prose and each
verified to match exactly one line. Two authorities were derived rather than
typed — the small-panel ratio the poster prints in bold and the label-side
site R² as a fraction — and `_num` learnt `&minus;`, so the three negative
values are compared with their sign. Each claim has its own self-test row; two
rows first rewrote the abstract and the manuscript as well as the poster, which
would have let their patterns report CAUGHT while the poster's was never proven
to fire, and were narrowed to poster-only markup. The re-sweep of what remains
unchecked on the poster found layout sizes inside HTML comments, a licence
version, design constants (the null-draw count, the q threshold, the panel
sizes), an inequality, one deliberately approximate "about", a published figure
from another group's paper, and HPC4's NSCLC value, whose only artefact is a
directory classified as not quotable. Each is left on purpose. The landing
page's remaining literals are dated runtime measurements with no stored
artefact behind them.

**And the coverage report had been misplacing wrapped claims.** A `WRAP:`
pattern matches a paragraph flattened to one line, and every such match was
recorded against the paragraph's FIRST line. So MISMATCH reports for wrapped
claims pointed at the wrong line, and `--coverage` listed a guarded literal as
unchecked on its real line: measured on this tree, 123 literals in nine
documents, 105 of them in the manuscript and 11 on the poster, were in the
unchecked list while a pattern guarded them. The mutation sweep is what exposed
it, by catching five poster literals the report said nothing covered. Matches
are now attributed to the line the literal sits on. The report's checked total
moved from 933 to 930 under identical patterns, because a few spurious pairs no
longer count, and its unchecked lists are 123 shorter and now true. The
poster's coverage floor rose from 19 to 79 and the landing page's from 13 to
16, so the new guards cannot disappear silently.

### Session 54, 2026-09-21, part 2: the self-test made cheap, the triple rule made strict, and the manuscript's Results bound

**The self-test was 80–88% of every suite, and it did not need to be.** Each of
its rows injected one wrong value and then rewrote and rescanned every scanned
document, so its cost was rows × documents and grew with every claim this
project added — which is the main thing this project's sessions do. `scan`
keeps no state across documents: each is read from its own path and compared
alone, so a document the row did not touch contributes exactly its baseline.
The loop now rescans only the documents a row changed. That was measured, not
argued: on the same tree and the same 395 rows, the old loop and the new one
produced identical output for every row — verdict, mismatch count, description
and substitution count — in 28 min 9 s and 3 min 19 s respectively, both under
the same concurrent load. The new loop also re-proves the equivalence on its
first three changing rows on every run, and asserts that the per-document
baselines sum to the whole-tree baseline, so the day `scan` acquires
cross-document state it fails loudly instead of drifting.

**A bracketed triple's point estimate could drift freely.** The triple rule
identifies WHICH quantity a triple quotes by its point estimate and then checks
the two bounds against it. A point that identifies nothing was reported and
passed, on the grounds, written into the rule, that HPC4 diagnostics and ssGSEA
arms were absent from the table; they have since been added. So the bounds
were guarded and the headline number beside them was not: a mutation sweep
moved seven Results points (the ssGSEA, mean-z, PLAGE, reconstruction and
global-axis arms, and the A7 acceptance arm) and `--strict` stayed green for
every one. It also ran per line, so a triple the typesetting had wrapped was
neither checked nor reported. The rule now matches each paragraph flattened to
one line, and an unidentified point is a failure. Measured before gating: 14
unidentified triples across three documents. Three were real quantities with
authorities and no triple row (the pinned ssGSEA arm, the 24-partition sds with
their chi-square intervals, the tumor-only NSCLC run's one outcome gain);
three were stored runs' secondary triples, now recognised through the same
stored-run whitelist the full-precision rule uses, and only when all three
numbers match; and three were measurements this repository never stored (two
HPC4 values, one diagnostic run), now listed exactly as written with where each
lives, so a changed digit still fails.

**The manuscript's Results, swept the way the poster was.** 121 literals the
coverage report listed as unchecked in the abstract, Results and Discussion
were each moved one step in their last printed place against the real file,
with `--strict` run and the file restored by sha256 each time. Before: 4 of
121 caught. The unguarded included the whole partition-variance paragraph —
both widened intervals, the chi-square intervals for the sds and the variance
shares they imply, the ratios 2.6490, 2.7556 and 1.3630 the argument turns on —
the Control C sample sizes, the size-gap Spearman statistics, the rotation
null's family means, the multiplicity margin and M_eff, the NSCLC site AUROC,
the post-ComBat AUROC, the tumor-only and subsample counts, and the live
abstract's cohort size and fold range, which wrap. Sixty-nine patterns were
added, each anchored on its own sentence, and every new authority was derived
from stored values rather than typed (the chi-square arithmetic with
`27_partition_sweep.py`'s own formula, the Spearman statistics the way figure 2
computes them, the cohort and SE ratios, the disattenuation inflations). Every one matched its intended line;
a few also reached agreeing copies in this file, the checklist and the cover
letter. After: 86 of 121 caught, and what remains is design constants (the
null-draw count, permutation counts, panel sizes, thresholds), an observed
range from 40 published articles, two inequality bounds, two deliberately
approximate figures, B-30's clause, and four counts no stored artefact carries.

**Three numbers checked against their method before being called anything.**
(1) The structured abstract folded below the manuscript's abstract was "461
words" in the same paragraph that defines a word as a whitespace token
containing a word character. 461 is a naive whitespace count; by the stated
measure it is 456. The manuscript now says 456 and the checker re-derives it;
the two dated records of the old figure carry brackets. (2) "It averages
+0.141 across the 200-gene sets" is right, but only under a reading the
sentence does not state: the eleven Hallmark sets of nominal size 200, which
keep 198–200 genes after filtering. The six that keep exactly 200 average
0.152 pan-cancer, which is what a reader recomputing it would likely get. The
method is recorded in the checker; the wording is flagged for the author
(B-34). *[2026-09-23: applied on the author's instruction — the manuscript now
says "the eleven 200-gene Hallmark sets (198–200 genes after filtering)", and the
checker binds the eleven's size range as well as the mean.]* (3) "Only 0.6–1.0% of sites" is right under the stored denominator,
the sites where each comparison group appears; over all 619 sites it would read
0.5–1.0%. Recorded, not changed.

**Two of the new global patterns first matched the wrong sentence.** One
captured 0.2576 as 0.257 in this file's A7 block, a truncated capture rather
than a discrepancy; the other matched this file's 24-partition sentence, which
has the Results sentence's shape and different numbers. `--strict` failed on
both and both were tightened before anything was committed. A pattern is a
claim about every scanned document, not just the one it was written for.

**And the coverage report's bookkeeping had a second fault.** A claim captured
with its sign ("−0.023") was recorded under a string no literal has, because
literals are extracted unsigned; 16 more guarded literals in five documents
were being listed as unchecked. Fixed with part 1's line attribution.

### Session 55, 2026-09-22: the cover letter's science bound, and an inequality that is not a measurement

**The outward-facing document had bound its own word counts and almost none of
its science.** Sessions 52 and 53 made the cover letter's length, its
paragraph-4 fallback and the subtraction between them checked claims, and
session 53 added the two counts AACR requires it to state about the
manuscript. Every one of those is a count of a document. The letter also makes
nine claims about *results* — the cohort size, the family-level permutation
test and its draw count, both random-set reliabilities, both reliability-gap
folds, the small-panel gap at each end of the sweep, and the frozen NSCLC
index an editor is explicitly invited to reproduce — and a mutation sweep
found all nine unguarded.

**Measured, not read off the coverage list.** Every literal `--coverage`
reported unchecked was moved one step in its last printed place against the
real file, `--strict` run WITHOUT `--fast`, and the file restored and its
sha256 compared each time. A positive control ran first — the body-length
claim, which is known to be guarded — and was CAUGHT, so an UNCAUGHT verdict
means something. Result: **40 distinct (line, literal) pairs, 2 CAUGHT, 38
UNCAUGHT**, the file restored byte-identically at the end. The two CAUGHT were
the pan-cancer confidence bounds, caught by the global triple rule session 54
made strict; everything else the letter says about the science was free to
drift. Eleven patterns were added and two widened, taking the letter from 18
checked literals to 32 and the checker from 1,134 claims to 1,159.

**A rule that checks an inequality as a point estimate reports a defect that
does not exist.** The reliability sentence reads "α falls from 0.97 to 0.80
while curated signatures hold above 0.97". The first two are measured
quantities and are now bound to `pancancer_alpha_null_raw` (0.9740) and
`pancancer_alpha_null_resid` (0.8012). The third was written as a claim about
the curated median and `--strict` refused it: that median is **0.9764**, which
rounds to 0.98, not 0.97. The number is not a rounded median — it is a
THRESHOLD the author chose, and 0.9764 > 0.97 is true. The pattern was
withdrawn rather than the tolerance loosened, which is the distinction the
project's own rule demands: *establish the method before declaring a defect*.
That clause is B-30, and it is the author's to reword. *[2026-09-23: reworded on
the author's instruction to "the median curated set held above 0.97", option
(a); the threshold stays unbound for the reason above.]*

**The two reliability figures turned out to be a four-document claim.** The
same sentence appears in the conference abstract, the NextGen abstract, the
manuscript (twice) and the letter, with a different verb each time — "falls",
"dropped", "reduced". A pattern anchored on the letter's verb saw one of the
five; anchored on "from X to Y while curated" it sees all five. Checked first
that no NSCLC variant of the sentence exists, because the NSCLC pair is
0.9565 → 0.8346 and a widened pattern would have compared it against
pan-cancer authorities and failed for the wrong reason.

**The paragraph enumeration was right, and its scope was not stated.** The
header tells the author that deleting paragraph 4 leaves 647 words and gives
the four remaining argument paragraphs and the boilerplate paragraph one count
each *[the counts are no longer quoted here: they are re-derived live, and the
letter was trimmed on 2026-09-24]*. Those sum to 608, not 647. Before calling that a
defect the document's own method was established — the words between `[Date]`
and the ORCID line, HTML comments removed — and under it every member is
exactly right and the paragraphs sum to the body exactly, 751. The missing 39
words are the date line, the address block, the salutation, the sign-off and
the byline, which no cut reaches and which the sentence never said it
excluded. This is the near-miss shape session 53 recorded: a recount that
disagrees is not yet a defect. The sentence now states its scope, and all
seven numbers — five members, their sum, and the scaffolding derived by
subtraction — are checked claims.

**0.2964 is deliberately unbound, and that is a decision rather than a gap.**
The letter's reproducibility promise names two numbers. 0.3182 is
`nsclc_isi` and is now bound. Its partner 0.2964 is the pre-fix index the
HPC4 run returned; it lives in `results/scorer_sensitivity_hpc4/run.log`, and
both `results/README.md` and the checker's own header already state that this
directory's numbers are diagnostic and not quotable. Binding it would mean
deriving an authority from a `.log`, which F6.8 decided against. It is
recorded here so a successor neither "fixes" it nor lets it drift unnoticed.

**A record that restated the enumeration was reworded before it could go
red.** Session 53's block in this file restated the per-paragraph list in the
live present tense. Every figure in it was correct — but B-24 is an open
decision to DELETE one of those paragraphs, and the moment the author acts on
it that record becomes false and `--strict` goes red pointing at a record of
the past. That is session 54's H37 exactly. Applying session 53's own H30
lesson — where a cross-document count cannot be made a checked claim, stop
asserting it and point at what re-derives it — the paragraph now points at the
header, with a dated bracket saying why.

**The baseline suite went red, and the red was the harness's own doing.** The
run was launched with its log inside the clone it was auditing, so
`test_results_readme_classifies_every_artefact_that_exists` correctly refused
a file the committed map could not mention. Diagnosed by measurement rather
than argument: copying this session's `results/README.md` into the clone makes
that test pass in 1.20 s with nothing else unclassified. The inherited commit
`f04367d` is green. The fix for successors is one word of path — write a
detached suite's log OUTSIDE the tree under test.

**A covariate this project had never isolated.** Every previous run at the
463-row self-test size was on AC. This one ran on battery on a quiet machine:
**288.22 s against 276.68 s**, +4.2%, with the suite at 466.69 s against
481.10 s. Battery is not a material covariate at this size, and a quiet
machine on battery is worth about as much as a loaded one on AC — the two
effects very nearly cancel.

### Session 60, 2026-09-24: A12 pinned (ledger H86, the author's option b)

**The code.** `splits.stable_group_kfold` is scikit-learn 1.6.1's non-shuffle
`GroupKFold` procedure with `np.argsort(..., kind="stable")` and nothing else
changed — ledger D1's in-process patch made permanent — and all three call
sites use it: `random_patient_split`, `models.site_prediction_control` and
`models.site_prediction_control_oof_combat`. Two tests hold it: it equals the
D1-patched splitter split for split on 150 random inputs (all-tied groups,
mixed sizes, string labels), and a 944-patient permutation gives the fold hash
`caf92a3f83ddb5b9` — the stable hash macOS and HPC4 were recorded agreeing on —
while tie-free groups reproduce the library's own folds. Both fail when the pin
is reverted. `26_sort_audit.py` case F reads `FIXED`.

**NSCLC, re-run with the shipped code** into `results/nsclc_v3_pinned/` (hashed;
PROVENANCE 222 → 233). The expectation was written first and held exactly: the
ISI is `nsclc_v3_stablesort/`'s bit for bit; Δ-MAE is 0.004284906287348778
[0.0013828565937186866, 0.0070256075863230905], bit for bit D1's macOS pinned
run, as are its 14 per-signature columns. New, because D1 stored no site table:
the site-control median AUROC over the same 15 sites is 0.9906433445573468
(0.9915909090909092 unpinned). Median Δr is +0.0132 (frozen +0.0232),
embedding over covariates +0.0258 (+0.0307). No conclusion changes. The
manuscript reports them beside the frozen values (Results, "Pinned
random-patient split"), bound to the stored run.

**Pan-cancer, re-run on macOS** into `results/pancancer_v3_pinned/` (hashed;
PROVENANCE 236 → 246). The ledger named HPC4, but the cluster resolves only on
the campus VPN, which was down all session (GlobalProtect appeared to wait for
a sign-in); the pinned split is the one D1 measured to give the same
secondaries on both platforms, so the laptop could stand in. A memory sampler
ran beside it: peak resident 6.2 GB, swap at most 6.5 GB, no kill — HPC4's
23.4 GiB MaxRSS had counted the whole address space. 5,825 s. The expectation
written first was half right: the ISI agrees with HPC4's corrected-ordering run
to 5e-16 (up to nine units in the last place, not the one expected), and the
per-signature excesses to 9e-15. The pinned secondaries: Δ-MAE
0.009265838536153103 [0.008271740109471818, 0.01022153374323296] (frozen 0.0088
[0.0078, 0.0098], corrected ordering 0.0096 [0.0086, 0.0106]); median Δr +0.0230
(+0.0199 frozen); embedding over covariates −0.0012 (−0.0000); median site AUROC
0.998794 over the same 98 sites (0.998477). The null-beating and outcome-arm
counts are unchanged. The manuscript reports both cohorts' pinned values beside
the frozen ones.

### Session 60, 2026-09-24: the 1,072 unchecked numbers in this file (ledger H89)

**The sweep.** Every numeric literal in this file that `--coverage` listed as
unchecked and not structural — 1,072 of them at `a79fcb0` — was moved one step
in its last printed digit, one at a time, and the file rescanned alone (the
self-test's exactness argument: the scan keeps no state across documents),
behind a positive control that was caught first (0.3857 → 0.3858; a first
control, a 16-digit ISI, was NOT caught, because a one-ULP move sits inside the
full-precision rule's tolerance — so it could not serve). Every restore was
sha256-compared. Result: 7 CAUGHT (other rules already held them), 1,065 not.

**What was bound: 288.** Each is a current result that an authority holds,
bound at printed precision with a self-test row (282 in `AUDIT_SWEEP_SPEC` in
`19_check_numbers.py` — four of them E15's site-attrition counts, recomputed
live from the inputs beside the site counts — and A5's six NSCLC partition
rows, scoped to their table's paragraph; every row caught). The shapes are each literal's
own line context, generated and checked to match exactly one line across every
scanned document. Some needed authorities that did not exist and now do — the
E1 tables' medians, E4's interval moves, E5's shifts, E9's null-reliability
floor (its records were hashed for this, PROVENANCE 233 → 236), E11's excess
difference, E14's label-side medians, E15's wider frame, Control C's A4b probe,
A5's widths and observed/predicted ratios, the A12 memo's corrected-ordering
Δ-MAE and mean Δr against the pinned run's — each read from a hashed run, and
derived quantities computed from the same stored values the tables print. The
file's coverage floor went 103 → 546 checked literals.

**What was not, and why** — per literal, with its class, in
`pipeline/results/session60_diagnostics/audit_sweep_20260924.tsv`:
- 564 sit in sections that record an investigation rather than state a result
  the paper reports (A9's platform hunt, A13–A15, D3/D5/D8, the F-workstream
  and session blocks). Their numbers that ARE current results were searched for
  separately and 28 of them bound (the fold-change census, the median curated
  α, two outcome excesses); the rest are the investigation's own record.
- 145 are other records: dated run records with units, PIDs, job ids, pointers,
  versions, external API counts (A8's GDC/PDC figures of 2026-09-04), and
  closed or superseded values (A9's 0.2964 stays unbound by decision, H49).
- 37 are design constants (thresholds, draw counts, set sizes).
- 10 are diagnostics whose sources are not hashed (the HPC4 per-prediction and
  per-signature comparisons behind A12's 0.193 fold agreement).
- 14 are results left unbound, each for a stated reason: 11 are rows of A5's
  two summary tables, which share a header and every row label, so neither a
  line nor a paragraph can tell them apart (the quantities they print — sd,
  bootstrap and combined SE, both intervals — are bound where the manuscript
  states them); 3 repeat a cohort count already bound elsewhere in the same
  sentence's sense (41,046 genes, 944 patients, 619 sites).

**A class this found, measured across the whole checker.** The sweep turned up
one number passing only because its pattern's tolerance FLOOR (0.005) was a
hundred times its printed precision: Supplementary Note 1's macOS ssGSEA score
sd, printed 0.0211 against a stored 0.021154 (which prints 0.0212). Every one of
the checker's 861 matched claims was then tested for the same looseness
(printed value against the authority rounded at the printed precision): 18
candidates, 16 of them integer percentages read with a spurious decimal, and 2
real — that one, and the results README's 24-partition bootstrap SE (0.011073),
which the SIX-partition pattern claimed and passed only by its 5e-5 floor.
Both fixed: the note says 0.0212, the four patterns' floors are 5e-7 so the
printed precision governs, and the 24-partition SEs have their own authorities
and bindings.

