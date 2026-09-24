# Supplementary Note 2 — sensitivity analyses in full

*Supplementary material for “Proliferation, not immunity: H&E-inferred tumor
microenvironment signatures are immune-specific but carry no prognostic advantage
over random gene sets.”*

The full text of six analyses that the manuscript's Results summarize, in the
order the Results give them: two further checks on the scorer, the label-side
site control in NSCLC, the sensitivity of the index to the site-to-fold
partition and to other single analysis choices, the site control at smaller
cohort sizes, and the ancestry arm. No claim in the manuscript depends on
reading it.
Each section is the manuscript's own text, so “above”, “below”, “the first
Results section” and numbered Limitations refer to the manuscript.

*(Until 2026-09-24 these passages were part of the manuscript's Results. They
were lifted out word for word, with a summary left in each place, so that the
manuscript's counted body fits a 5,000-word limit.)*

## 1. A third scorer, and a second reliability estimator

**A third scorer places ssGSEA's reversal, not the index's sign, as the outlier.**
Two scorers cannot separate "not a mean composite" from "a rank-walk statistic",
so we added PLAGE — the leading singular vector of each set's standardized
genes, a data-weighted composite — and ran the same NSCLC arms under current
code. On the uncorrected contrast PLAGE gives +0.1353 [0.0889, 0.1834] with 13
of 16 signatures above their null, against +0.2963 [0.2442, 0.3514] and 16 of 16
for mean-z in the same run. The registered estimand is undefined under PLAGE as
it is under ssGSEA, for the same reason: the reliability estimator assumes a mean
composite. Under the scorer-agnostic reconstruction, though, PLAGE stays positive
(+0.0772, 11 of 16 above null) where ssGSEA reverses. So the reversal is specific
to ssGSEA rather than a property of every alternative to the mean, while the
attenuation of magnitude is common to both alternatives.

**The reliability estimator itself changes the magnitude more than the scorer
does.** Cronbach's α assumes every gene loads equally on one factor. Replacing
it with McDonald's ω from a one-factor model of the same residualized items
leaves the curated sets almost unchanged (median ω 0.977 against α 0.977) and
collapses the random sets (median ω 0.227 against α 0.836), so the
curated-versus-random reliability gap widens from 0.126 to 0.733. Reconstructed
with ω in place of α, the NSCLC index falls to 0.2147 [0.1739, 0.2554] from
0.2917 [0.2584, 0.3251] with α on the same estimator, and 15 of 16 signatures
still exceed their null. Read this as a bound rather than a better estimate: a
one-factor ω is only interpretable where one factor fits, and for random sets
after axis removal it does not, which is exactly why their ω is so low. Both
readings agree on the direction and on 15 or 16 of 16 signatures; they disagree
on magnitude by about a quarter.

## 2. The label-side site control in NSCLC

The same control on the NSCLC scope is a third the size: a median 14.5% of the
label across 31 usable sites and 796 patients, against 44% across 202 sites
pan-cancer (per signature in Supplementary Table S5). The gap is cancer type,
not the dummy coding. Adjusting each fit for its number of site terms lowers
the medians by similar small amounts, to an adjusted 42.4% pan-cancer and 11.2%
in NSCLC. Cancer type alone accounts for a median 38.5% of label variance
pan-cancer and 3.6% across the two NSCLC types, so pan-cancer the site term is
largely standing in for disease. Conditioned on type, site explains a median
9.9% in NSCLC and 5.5% pan-cancer: beyond disease, the NSCLC labels carry at
least as much site structure. The within-type permutation calibration agrees,
with a median of 7.0% across the 16 signatures in NSCLC (15 of 16 at
p ≤ 0.05). The two cohorts are therefore tabulated separately rather than
stacked.

## 3. Sensitivity to the choice of site-to-fold partition

Preserved-site cross-validation admits many valid assignments of sites to folds,
and the analysis uses one. We re-ran the primary endpoint under six assignments
in each cohort, varying only the partition with the null gene sets and bootstrap
draws held fixed, and later under 24 assignments in each cohort (per partition,
both cohorts and both sweeps, in Supplementary Table S6).

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

The two cohorts differ in how much of the total they attribute to partition
choice — 25.7% pan-cancer against 8.4% in NSCLC, in the cohort 7.6 times larger
— and that ordering runs against the intuition that a larger cohort is a safer
one. At six and five partitions we reported the contrast without building on it,
because so few partitions cannot establish it; the 24-partition sweep below
does. A chi-square interval for the true
between-partition standard deviation, at one fewer degree of freedom than the
number of partitions, spans [0.0056, 0.0218] in NSCLC and [0.0039, 0.0187]
pan-cancer; propagated through the variance share those become [3.4%, 35.5%]
and [11.0%, 74.0%]. Each standard deviation is pinned only to within a factor
of about four, which is what five draws buy, and both pairs of intervals overlap
across most of their range. The point estimates are ordered; these data do not
establish that the quantities behind them are.

What *is* precisely estimated is the bootstrap component, and it behaves as
sampling theory says it should: between the two cohorts it falls by a factor of
2.6490, against the 2.7556 predicted by the ratio of the square roots of the
sample sizes — agreement within 4%, and a useful check that the bootstrap is
doing what it claims. Over the same step the partition component falls by only
1.3630, roughly half the rate predicted either by patient count (2.7556) or by
tissue-source-site count (3.0171, for 619 sites against 68). We record that
asymmetry because it is present in the point estimates, and we offer no
explanation for it: neither candidate scaling accounts for it, and separating a
genuine design floor from five draws of a noisy variance estimate would take on
the order of 20 to 30 partitions per cohort. At six and five partitions, the
claim the analysis supported was therefore the narrow one it began with — the
reported bootstrap interval understates total uncertainty, by 4.5% in NSCLC and
16.0% pan-cancer — and not a general statement about how grouped-partition
variance scales with cohort size. We then ran 24 partitions per cohort.

**Twenty-four partitions per cohort.** These runs use the pinned tie ordering
described in Limitation 8, so they sit beside the figures above rather than
replacing them, and every one of the 48 was computed. Thirteen pan-cancer
partitions did not finish under the full configuration: eight exceeded the
computing cluster's four-hour limit, and five stopped in a secondary analysis
(the singular-value failure described below). They were completed with the two
secondary analyses switched off, since the index never reads them, and nine of
them also used a shortcut for the null refits that shares one factorization
among gene sets with the same missing values. Both variants were checked to
reproduce the index bit for bit, and the three partitions computed more than
once agree bit for bit, one of them under all three configurations.

In NSCLC the standard deviation across the 24 partitions is 0.0077, with a
chi-square 95% interval of [0.0060, 0.0107]. Against the same bootstrap standard
error, the NSCLC interval understates total uncertainty by 3.3%, and partition
choice makes up 6.4% [3.9%, 11.8%] of the NSCLC total variance. Pan-cancer, the
standard deviation across the 24 partitions is 0.0072 [0.0056, 0.0101]; the
pan-cancer interval understates total uncertainty by 19.4%, and partition choice
makes up 29.9% [20.5%, 45.6%] of the pan-cancer total variance. The
pinned-ordering values
replace neither the six-partition NSCLC share of 8.4% nor the five-partition
pan-cancer figures; each is reported with the sweep it came from, and the
24-partition values are the better estimates of both.

With 24 draws the scaling question has an answer. The bootstrap component again
falls with cohort size, by 2.652 against the 2.756 predicted by the square
roots of the sample sizes. The partition component does not fall: its ratio
between the cohorts is 1.060, with an F-based 95% interval of [0.697, 1.611].
That interval excludes both the patient-count prediction and the site-count
prediction (3.017), and it contains 1. The two standard deviations' intervals
overlap; the two shares' intervals do not. Over these two cohorts, the
spread across partitions behaves as a floor set by the partition design rather
than by cohort size. Because the bootstrap error shrinks as the cohort grows,
that floor is a larger share of the total uncertainty in the larger cohort,
which is why the ordering of the shares runs against intuition. Two cohorts
cannot show how the floor depends on the design, and we do not claim that
they do.

Two further observations. The sixth pan-cancer assignment was not computed,
and the cause was not the partition. The run stopped when a singular value
decomposition failed to converge inside the covariate-baseline regression, a
secondary analysis fitted under the random-patient split that the index does
not use. The design matrix there is well defined: it is rank-deficient because
some one-hot site columns are empty in a training fold, and LAPACK's alternative
driver decomposes it without difficulty. Repeating that regression for 24
random-patient splits on the same macOS machine reproduces the failure for
exactly that partition and for no other. On the Linux cluster, where the
random-patient split differs (Limitation 8), the same 24 splits fail that
regression for four other partitions, and the variance decomposition, another
secondary analysis, for a fifth: those are the five 24-partition runs above
that stopped, and the primary partition passes on both machines. The index
could therefore have been computed for all six assignments, and the reported
spread uses the five that completed. And this sensitivity is
specific to the index: for the split-robustness contrast Δr, which is a
difference *between* two split schemes rather than a contrast computed within
one, partition choice matters considerably more.

## 4. Sensitivity to other analysis choices

None of the analyses in this section was pre-registered. Each changes one choice
and recomputes only the index. In NSCLC the unchanged configuration returns the
corrected-ordering value, 0.2922 [0.2355, 0.3527], bit for bit, and every
variant is compared against it.

*Where the dominant expression axis is removed.* The registered index removes
it within each cancer type. Removing a single axis computed across the two
NSCLC types instead collapses the NSCLC index to −0.0096 [−0.0630, 0.0466], with
1 of 16 signatures above its null. Random gene sets become as predictable from
the image as the curated ones (median r 0.361 against 0.351), and 10 of the 16
nulls become degenerate. Pan-cancer, the same change moves the index the other
way, to 0.3512 [0.2692, 0.4353] against 0.2968 [0.2749, 0.3189], with 15 of 16
signatures above their null. There, curated and random sets both become far
more predictable (median r 0.725 and 0.452, against 0.391 and 0.100). The index
is therefore a statement about expression variation within a cancer type, which
is how it was registered. It does not survive an axis defined across types
unchanged, and the direction of the change depends on the cohort.

*Expression matching of the null.* Matching random sets on size alone, rather
than on size and mean expression, gives 0.3041 [0.2434, 0.3536]: expression
matching lowers the NSCLC index by 0.012. This decomposes the registered
estimand; it is not a test of its robustness.

*The null subsample inside the bootstrap.* Recomputing each null mean from 400
or from all 1,000 draws, rather than 200, leaves the estimate unchanged and the
interval nearly so: [0.2336, 0.3512] and [0.2333, 0.3511].

*The number of folds.* With 3 folds the NSCLC index is 0.2188 [0.1673, 0.2722],
with 15 of 16 signatures above their null; with 10 folds it is
0.3154 [0.2607, 0.3711]. The index rises with the training data each fold
provides, so its value is specific to the registered five folds. Its sign, and
the excess of curated over random sets, hold at every number of folds tried.

*Tumor-only expression.* Neither cohort restricted expression to the
primary-tumor sample, although every slide is a primary-tumor diagnostic
slide. In NSCLC, 51 of 944 patients were paired with a profile from other
tissue: 48 adjacent normal, 2 recurrence, and 1 patient with only a normal
sample. Pan-cancer, 578 of 7,168 profiles averaged several samples, and 40
patients had no tumor sample at all. Restricting expression to the tumor
sample raises the NSCLC index to 0.3448 [0.2731, 0.4079] (943 patients), so
adjacent-normal profiles were diluting the association. It moves the
pan-cancer index to 0.2937 [0.2701, 0.3149] (7,128 patients). All values
elsewhere in this paper use the original inputs.

## 5. The site control at smaller cohort sizes

An AUROC this close to 1 invites the question whether it reflects the embeddings
or the amount of data each site contributes, so we recomputed the control on
random patient subsamples. It does not depend on cohort size in any range we
can evaluate. Pan-cancer the median stays at 0.998 at the full 7,168 patients,
0.998 at half and 0.994 at a quarter (1,792 patients, 10 to 14 evaluable sites);
in NSCLC, 0.989, 0.980 and 0.894 to 0.972 over the same steps. This
stress test was run on Linux, where the control's folds are assigned
differently (Limitation 8), which is why its full-cohort NSCLC median reads
0.989 here and 0.992 above. Below a quarter
of either cohort no site retains the 20 patients the control requires, so the
limit we reach is site eligibility and not discriminability. Lowering the
threshold to 5 patients instead evaluates 47 NSCLC sites (median 0.987, 94%
above 0.9) and 320 pan-cancer sites (median 0.998, 97% above 0.9), so sites
contributing as few as five patients are still identifiable.

## 6. Ancestry

A pan-cancer ancestry arm (6,580 patients; EUR 5,312, AFR 594, ASIAN 502, AMR
172) is reported as Supplementary Table S8 with a negative headline: **the
contrast is not identifiable from tissue source site**. Only 0.6–1.0% of sites
carry at least 10 patients of both the reference and any comparison group
(ancestry × site Cramér's V = 0.537 against a permutation null with 95th
percentile 0.088). Ancestry calls are the published UCSF consensus set
(12). Any TCGA study reporting ancestry
subgroup differences in an
image model faces this constraint; most do not check it.
