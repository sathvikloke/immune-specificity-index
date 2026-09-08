# Immune-Specificity Index (ISI)

**Does H&E-inferred tumor microenvironment signal carry anything immune-specific
— or does it merely track a global expression axis that any random gene set of
matched size and expression would track equally well?**

This repository is the curated public snapshot of the code and frozen results
behind that audit. It is a *snapshot*, not a development repository: it is
published as a single commit so that no internal working material can survive in
its history.

---

## The answer, in two parts

**1. The immune-specific signal is real.** After residualizing every score on
the dominant within–cancer-type expression axis and disattenuating for the
measurement reliability of the comparator, curated TME signatures remain
substantially more image-predictable than size- and expression-matched random
gene sets.

| | pan-TCGA (n = 7,168, 31 types) | NSCLC (n = 944) |
|---|---|---|
| **ISI** | 0.2911 [0.2692, 0.3133] | 0.3182 [0.2614, 0.3763] |
| Signatures beating their null | 16 / 16 | 16 / 16 |
| Outcome (PFI, within type) | 10 / 32 | 0 / 32 (MDE 0.046) |
| Site AUROC (median) | 0.998 | 0.992 |

**2. That signal is not prognostic — and this is the part that matters
clinically.** Measured against the same matched-random comparator, the Hallmark
signatures that beat their null on progression-free interval are proliferative
and stromal: G2M checkpoint, E2F targets, angiogenesis, epithelial–mesenchymal
transition, and hypoxia. **No interferon, inflammatory, complement, or
allograft-rejection signature does.** Since those are precisely the programs
image-based TME inference is proposed to deliver, image-inferred immune scores
should not be used as prognostic biomarkers — not because the image reads them
poorly, but because the underlying transcriptomic scores are not prognostic once
compared fairly.

**A caveat the headline does not survive without.** The index is defined on a
mean-of-z composite. Under single-sample GSEA the uncorrected contrast survives
at about half its magnitude, but the *corrected* estimand does not exist at all:
disattenuation divides by the square root of the comparator's reliability, and a
rank-walk scorer drives random-set reliability from 0.799 to 0.258. Any index
built on a curated-versus-random comparison has to name its scorer as part of its
definition.

---

## What is here

```
05-PRE-REGISTRATION.md   the protocol
09-PAPER-DRAFT.md        the manuscript
14-SCIENCE-AUDIT.md      the internal audit, including everything that failed
CITATION.cff             citation metadata
LICENSE                  MIT, plus a note on data licensing
pipeline/
  src/aacr27/            the library
  scripts/               01-26, in dependency order
  tests/                 the test suite
  tools/                 detach helper
  results/               frozen result JSON/CSV + PROVENANCE.json
  README.md              pipeline documentation -- start here
  reproduce.sh           one-command reproduction
  requirements-lock.txt  exact pinned versions
  ENVIRONMENT.md         platform, BLAS and seeds of the reported run
```

**Deliberately not here:** raw input data (public, but large — `pipeline/README.md`
says where to get every source), internal handoff documents, operational and
personal material, and cluster job scripts.

## Reproducing

```
cd pipeline && bash reproduce.sh --check
```

`--check` verifies the environment, the frozen artifact hashes, the numeric
claims in the documents and the full test suite. It takes about four minutes and
exits 0. The full analysis run is a separate, much longer path — see
`pipeline/README.md`, and note that its total wall clock is an **estimate that
has never been measured end to end**.

Every reported number is hashed in `pipeline/results/PROVENANCE.json` and
mechanically compared against the documents that quote it, so a document and an
artifact cannot silently disagree.

## Two things to read honestly

**Platform.** The reported values were produced on macOS 15 / arm64 with
Python 3.13.9. The estimate is platform-dependent in its fourth digit; a
non-stable sort in fold assignment meant two machines built different
cross-validation partitions from byte-identical inputs. That is fixed, and both
the pre-fix and post-fix numbers are reported rather than one replacing the
other. See Limitations in the manuscript.

**Pre-registration.** The protocol in `05-PRE-REGISTRATION.md` was fixed before
the reported results were computed, and is timestamped by an annotated tag in the
private working repository. **That tag is not present here**, because this
snapshot is published as a fresh single-commit repository with no history — the
same property that keeps internal material out of it also means it carries no
commit predating the results. The protocol's *content* is public; its *timestamp*
is not independently verifiable from this repository alone. Stated plainly rather
than implied, because a pre-registration whose timestamp cannot be checked should
be read as self-attested.

## License

Code is MIT (`LICENSE`). The same file carries a **NOTE ON DATA**: the license
covers the code only, no dataset is redistributed here, and each source is used
under its own terms — Prov-GigaPath TCGA embeddings (CC-BY-4.0), MSigDB Hallmark
v2024.1 (CC-BY-4.0), Xena TOIL expression, PanCanAtlas ABSOLUTE purity, TCGA-CDR,
and published UCSF ancestry calls. BostonGene/MFP signature definitions are
deliberately **not** used, because that license restricts use to academic and
non-profit purposes and defines "Source Form" to include data.

## Citing

See `CITATION.cff`.
