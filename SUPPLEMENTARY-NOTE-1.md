# Supplementary Note 1 — the full forensic account of the ordering defect

*Supplementary material for “Proliferation, not immunity: H&E-inferred tumor
microenvironment signatures are immune-specific but carry no prognostic advantage
over random gene sets.”*

The two corrected re-runs, the cross-platform validation to sixteen significant
figures, the scope of the library ordering, and the pipeline-wide audit that found
one further instance. It is evidence for Limitation 8 of the manuscript rather
than a further limitation, and no claim in the manuscript depends on reading it.
Where the text below says “this limitation” it means Limitation 8, and results
“reported above” are the manuscript's Results and Limitation 8.

*(Until 2026-09-24 this note was a folded block inside Limitation 8. It was
lifted out word for word, so that the manuscript's counted body excludes it.)*

*The mechanism, and two accounts of it that these measurements overturn.*
Against the partition a stable sort gives on both machines, 216 of 944
patients change fold on macOS and 180 on Linux. Because the training
partitions differ, the per-fold ridge penalties differ too — 11 of 80, each
by one grid step — so those flips are a consequence and not the cause. The
ridge solve is not implicated either: the regularized system is well
conditioned (condition number 18–77) and a relative 10⁻¹⁴ input perturbation
moves the predictions by 6 × 10⁻¹⁵, so it attenuates rather than amplifies.
Two earlier accounts of this limitation — one arguing from the width of the
penalty selection margins that the penalties *could not* differ across
platforms, the other attributing the divergence to amplification through an
ill-conditioned ridge system — are contradicted by these measurements and
are corrected here.

That separately reported re-run now exists. Re-running the pre-registered
NSCLC configuration unchanged except for the two corrected orderings gives an
ISI of 0.2922 [0.2355, 0.3527], against the pre-registered
0.3182 [0.2614, 0.3763]. The point estimate moves by 0.0260, which is 2.9
times the between-partition standard deviation reported above (3.4 times the
24-partition value measured under the pinned ordering) and falls
inside the honest interval [0.2516, 0.3720] that the partition-variance
analysis derived. It is nonetheless a large move for a change of partition:
of the 276 pairs among the 24 partitions, 3 differ by as much. All 16 of 16
signatures still exceed their nulls and the
interval still excludes zero. We report both numbers rather than replacing
one with the other: the pre-registered value is what was registered and run,
and the corrected value is what the fixed code produces, and a reader is
entitled to see the difference rather than a single number chosen after the
fact.

The equivalent re-run for the pan-cancer cohort has now been performed, and
the correction is smaller there than in NSCLC. Re-running the frozen
pan-cancer configuration unchanged except for the two corrected orderings
gives an ISI of 0.2968 [0.2749, 0.3189], against the frozen
0.2911 [0.2692, 0.3133]. The point estimate moves by 0.0057, which is
**0.87 times the pan-cancer between-partition standard deviation of 0.0065**
(0.79 times the 24-partition value) — that is, by less than changing the
cross-validation partition typically does, since 153 of the 276 pairs among
the 24 partitions differ by more — and it
falls inside the honest interval [0.2683, 0.3187] that the partition-variance
analysis derived. All 16 of 16 signatures still exceed their nulls and the
interval still excludes zero. In NSCLC the same fix moved the estimate by
2.9 partition standard deviations and here by 0.9, a factor of 3.3 (4.3 in
24-partition units), and the share of the partition it changed does not
explain the contrast. The fix reorders only sites of equal size: comparing
each cohort's frozen and corrected fold assignments, no patient in an
untied site changed fold and every fold kept its size. Pan-TCGA has far
more tied sites, 600 of 619 holding 79.5% of its patients against 51 of 68
holding 38.1% in NSCLC, and the fix moved 3,812 of 7,168 pan-cancer
patients (53.2%) to a different fold against 216 of 944 (22.9%) in NSCLC.
The adjusted Rand index between each cohort's two partitions is 0.13
pan-cancer and 0.54 in NSCLC. In pan-cancer the fix therefore amounted to a
nearly unrelated partition and moved the index by a typical amount for
one; in NSCLC it changed less of the partition and moved the index by more
than almost any pair of the 24 partitions differs.

This run was executed on Linux/x86_64 rather than on the macOS/arm64 machine
that produced the frozen results, because scoring the pan-cancer cohort
requires a dense copy of a 10,535 × 41,046 expression matrix — 3.2 GiB before
the scorer makes a second — held for the hours the run takes, and the
available laptop could not guarantee that: sampled every 30 s during test
runs on 2026-09-15 and 2026-09-16, its free memory fell as low as 2.6 GiB
and moved by up to 8.8 GiB between adjacent samples. **Reporting a
cross-platform number requires justifying it, and the justification is a
measurement rather than an argument.** The concern this item exists to
document is precisely that two platforms disagreed; under the *shipped*
ordering they did, and that disagreement is this limitation. Under the
corrected ordering they do not. Before the pan-cancer run, the pre-registered
NSCLC configuration was re-run on the same Linux machine under the corrected
ordering and compared against the macOS result reported above: the two agree
to sixteen significant figures, differing by one unit in the last place on
the point estimate and the upper bound and not at all on the lower
(0.2922124012879651 against 0.2922124012879652). That is six orders of
magnitude below the discrepancy this item documents, on the same cohort,
the same 68 sites, the same 1,000 null draws and the same patient-clustered
bootstrap. The platform is therefore validated for the index under the
corrected ordering, and the pan-cancer value above is reported on that
basis.

That validation covers the index and the preserved-site partition it is
computed on, and nothing wider. The random-patient partition behind Δ*r*,
Δ-MAE and every covariate baseline is assigned by a library routine
(scikit-learn's GroupKFold), which orders equally sized patient groups by the
same architecture-dependent rule. Under the corrected ordering the two
machines place only 19.3% of NSCLC patients in the same random-patient fold,
which is what chance assignment gives. Per-signature Δ*r* differs between
them by up to 0.034, and Δ-MAE is 0.0052 [0.0022, 0.0080] on macOS against
0.0048 [0.0020, 0.0075] on Linux, while the preserved-site predictions agree
to 1.3 × 10⁻¹⁴. Neither conclusion drawn from these secondary quantities
changes, but their digits are platform-specific. Every Δ*r*, Δ-MAE and
covariate-baseline value in this paper is the macOS value. The same library
routine assigns the folds of the site-classification control, so its
medians are platform-specific too: the full-cohort values under Negative
controls are the macOS ones, and the subsample stress test that follows
them was run on Linux, which is why its full-cohort NSCLC median reads 0.989
where the macOS run gives 0.992.

That library ordering is now pinned (2026-09-24): the pipeline assigns these
folds with the library's own procedure under a stable sort, which is identical
to the library wherever group sizes do not tie. Measured before the change by
patching the library in-process on both machines, the pinned partition gives
the same secondary quantities on macOS and on Linux to round-off, and the
shipped code reproduces that patched run bit for bit. The values quoted
elsewhere in this paper are the frozen runs' and stay as run; the pinned
re-run's are reported beside them in Results.

As in NSCLC, we report both numbers rather than replacing one with the other:
the frozen value is what was run under the registered estimand, and the
corrected value is what the fixed code produces. Nothing reported for pan-TCGA
depends on the re-run, because the frozen estimate remains the reported one and
the two cohorts are still reported on the same footing; what the re-run adds
is that the ordering's effect on the estimate is now measured in both
cohorts rather than measured in one and bounded in the other.

A systematic audit of every other non-stable ordering in the pipeline's own
code (which the library ordering described above lies outside), run after
this cause was identified, found **one further instance and it is on
the scoring path**: the single-sample GSEA rank table ordered tied genes by
the same architecture-dependent rule. Ties there are not incidental — log
expression has a floor, every gene resting on it shares an average rank, and
all 944 samples carry ties, 35,394,448 of 38,747,424 rank entries in total.
The two orderings produce different scores for 15,060 of 15,104
patient–signature cells across all 16 signatures. The largest score
difference between the two orderings is 0.0128, against a mean
per-signature score standard deviation of 0.0212 (macOS; 0.0130 and
0.0211 on Linux, where the default-sort scores differ too). The ordering is
confirmed architecture-dependent by hashing the rank table on both machines.
Under the default sort it hashes to `0bca973cfa26da89` on macOS and
`b98fb11e701e8e98` on Linux, and under the pinned sort to
`096172c056cc7a34` on both.
This does not touch any result above that is computed from the
mean-of-z-scores composite, which does not use the rank table. It does reach
the ssGSEA sensitivity arm, which has now been recomputed under the pinned
ordering: the uncorrected contrast moves from +0.180 to
+0.174, 15 of 16 signatures beat their null either way, and the
estimand stays undefined under ssGSEA for the reason given in Limitation 9.
A 0.61-standard-deviation movement in the per-cell scores (0.62 on Linux)
thus produces a
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
no unpinned order-dependent operation in the pipeline's own code. A separate
sweep of the sorting and grouping calls the pipeline makes into its
libraries found the library ordering described above at three call sites —
the random-patient split and the two site-classification fits — none of
which the index reads; that ordering was documented here rather than pinned
until 2026-09-24, and is pinned now, as described above.
